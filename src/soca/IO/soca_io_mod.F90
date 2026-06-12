! (C) Copyright 2026 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! soca_io_mod
!
! Direct netcdf I/O for SOCA fields, replacing FMS register_restart_field /
! save_restart / restore_state. Enables LETKF parallel-ensemble I/O: FMS forces
! every read/write to be collective on the compute communicator, blocking
! "PE i alone reads/writes member i"; direct nf90_* lets a single PE do its
! own I/O without coordinating with the others.
!
! API: register-then-commit. Writer pattern:
!
!     type(soca_io_writer) :: w
!     call w%init(domain, "myfile.nc")
!     call w%enqueue("lonh", self%lonh)
!     call w%enqueue("lath", self%lath)
!     ...
!     call w%commit()   ! gathers from all PEs and writes the file
!
!
module soca_io_mod

use netcdf
use mpi
use kinds, only: kind_real
use mpp_mod, only: mpp_gather, mpp_scatter, mpp_broadcast, mpp_pe, mpp_root_pe, mpp_npes, &
                   mpp_get_current_pelist, mpp_error, FATAL
use mpp_domains_mod, only: domain2D, &
                           mpp_get_compute_domain, mpp_get_global_domain, &
                           mpp_get_data_domain
use fckit_configuration_module, only: fckit_configuration

implicit none
private

public :: soca_io_writer
public :: soca_io_reader
public :: soca_io_file_exists, soca_io_var_exists
public :: soca_io_config_from_yaml
public :: soca_io_ensemble_write_parallel
public :: soca_io_ensemble_batch_size
public :: soca_io_ensemble_root_pe
public :: soca_io_writers_commit_ensemble
public :: soca_io_readers_commit_ensemble

integer, parameter :: MAX_NAME       = 256

! Module-level I/O dispatch knobs, set once from geometry.io YAML at
! soca_geom_init via soca_io_config_from_yaml; persist for the run. Defaults are
! tuned for cluster-scale production (parallel write, scatter ensemble read,
! async MPI); single-state read stays strided to preserve direct-io behavior.
logical, save :: ensemble_write_parallel_   = .true.
logical, save :: ensemble_read_scatter_     = .true.
logical, save :: single_state_read_scatter_ = .false.
logical, save :: async_mpi_enabled_         = .true.
! Ensemble I/O batch size: how many members are processed (built, committed,
! freed) in one bulk read/write pass. 0 (the default) means a single batch over
! all members -- maximum I/O overlap, peak memory ~ all members resident at
! once. A positive N caps the working set to N members per pass, trading I/O
! concurrency for a memory ceiling (set N ~ npes for one wave per PE). Consumed
! by the soca_fields ensemble drivers.
integer, save :: ensemble_batch_size_       = 0
integer, parameter :: MAX_FILE_NDIMS = 8

! Reader is stateless across commits: every reader_commit nf90_opens the file,
! reads all enqueued vars, and nf90_closes. Within a single commit the var
! loop reuses the open ncid but nf90_inq_varid + inquire_variable +
! inquire_dimension run inline per var (microseconds; not worth caching).
! Holding NetCDF4 handles open across commits is what bloats LETKF memory:
! each open HDF5 file holds MB of metadata + chunk cache, and dozens of
! ensemble-member files per PE accumulates to GB-scale per-PE growth.

type :: var_entry
  character(len=MAX_NAME) :: name = ''
  integer :: ndims = 0          ! 1, 2 or 3
  integer :: nlevels = 0        ! 1 if ndims<=2, else nz
  ! Pointers to the caller's buffer (no copy). 1D: global on every PE
  ! (no gather, written directly by root). 2D/3D: the caller's halo-inclusive
  ! data-domain array; commit extracts the compute-domain slice into a per-var
  ! tile for mpp_gather. Caller must keep these alive and unmutated through
  ! commit (matches the FMS register_restart_field contract).
  real(kind=kind_real), pointer :: data1d(:)     => null()
  real(kind=kind_real), pointer :: data2d(:,:)   => null()
  real(kind=kind_real), pointer :: data3d(:,:,:) => null()
  character(len=MAX_NAME) :: long_name = ''
  character(len=MAX_NAME) :: units = 'none'
  ! Inter-stage gather buffers used by writer_stage_gather -> writer_stage_write.
  ! Allocated lazily on the writer's root_pe in stage_gather (size nx_g x ny_g
  ! for 2D, nx_g x ny_g x nlevels for 3D); freed in stage_close. Holding them
  ! per-var (rather than per-writer reused) lets writer_stage_gather populate
  ! every var before any disk write starts -- so in the ensemble path the
  ! writes in stage_write can run concurrently across writer PEs with no MPI.
  real(kind=kind_real), allocatable :: gbuf_2d(:,:)
  real(kind=kind_real), allocatable :: gbuf_3d(:,:,:)
end type var_entry

! Tracks one unique axis (X/Y/Z) by (size, domain_key) and remembers the netcdf
! dim/coord-var ids assigned during commit. Same dedup as FMS unique_axes(): a
! var's dim reuses an existing axis if its (size, domain) tuple matches, else a
! new axis is appended.
type :: axis_entry
  integer :: size       = 0
  integer :: domain_key = 0      ! 0 = no domain (used for 1D vars); 1 = main domain
  integer :: dimid      = -1     ! netcdf dim id (valid only on root after define)
  integer :: varid      = -1     ! 1D coord var id (valid only on root)
end type axis_entry

type :: soca_io_writer
  character(len=:), allocatable :: filename
  type(domain2D), pointer :: domain => null()
  integer :: isc, iec, jsc, jec        ! local compute domain (1-based as mpp returns)
  integer :: isd, ied, jsd, jed        ! data domain (for compute-slice offsets)
  integer :: isg, ieg, jsg, jeg        ! global domain
  integer :: nx_g, ny_g                ! global x/y sizes
  type(var_entry), allocatable :: vars(:)
  integer :: nvars = 0
  ! Writer's root PE: -1 sentinel means writer_init defaults it to mpp_root_pe()
  ! at use time (single-shot path). The ensemble orchestrator sets a rotated
  ! root_pe before stage dispatch so different members write from different PEs.
  integer :: root_pe = -1
  ! Inter-stage state populated by writer_stage_define and consumed by the later
  ! stages. ncid is held open between define and close (only meaningful on root).
  integer :: ncid = -1
  integer, allocatable :: varids(:)
  ! Per-direction unique-axis tables built in stage_define; stage_close frees.
  type(axis_entry), allocatable :: x_axes(:), y_axes(:), z_axes(:)
  integer, allocatable :: var_x_idx(:), var_y_idx(:), var_z_idx(:)
  integer :: nx_axes = 0, ny_axes = 0, nz_axes = 0
contains
  procedure :: init => writer_init
  procedure :: enqueue_1d => writer_enqueue_1d
  procedure :: enqueue_2d => writer_enqueue_2d
  procedure :: enqueue_3d => writer_enqueue_3d
  generic   :: enqueue => enqueue_1d, enqueue_2d, enqueue_3d
  procedure :: commit => writer_commit
end type soca_io_writer


! Read-side var_entry: holds POINTERS to caller buffers; commit fills them in
! place. Caller must keep buffers alive and unmutated between enqueue and commit.
type :: read_entry
  character(len=MAX_NAME) :: name = ''
  integer :: ndims = 0           ! 1, 2, 3, or 4
  real(kind=kind_real), pointer :: data1d(:)       => null()
  real(kind=kind_real), pointer :: data2d(:,:)     => null()
  real(kind=kind_real), pointer :: data3d(:,:,:)   => null()
  real(kind=kind_real), pointer :: data4d(:,:,:,:) => null()
  ! Inter-stage staged buffers for the scatter path: reader_pe pulls the
  ! WHOLE global field into these during reader_stage_read; reader_stage_distribute
  ! ships the per-PE compute tile out via MPI_Scatterv and the buffers are freed
  ! in reader_stage_close. Only allocated on reader_pe (1x1 dummies elsewhere).
  real(kind=kind_real), allocatable :: gbuf_2d(:,:)
  real(kind=kind_real), allocatable :: gbuf_3d(:,:,:)
  real(kind=kind_real), allocatable :: gbuf_4d(:,:,:,:)
end type read_entry

type :: soca_io_reader
  character(len=:), allocatable :: filename
  type(domain2D), pointer :: domain => null()
  integer :: isc, iec, jsc, jec        ! compute domain (1-based, as mpp returns)
  integer :: isd, ied, jsd, jed        ! data domain (1-based)
  integer :: isg, ieg, jsg, jeg        ! global
  integer :: nx_g, ny_g                ! global x/y sizes
  type(read_entry), allocatable :: vars(:)
  integer :: nvars = 0
  ! Reader's source PE for the scatter path: -1 sentinel means reader_init
  ! defaults it to mpp_root_pe(). The ensemble orchestrator sets a rotated
  ! reader_pe per member so different files are read concurrently across PEs.
  integer :: reader_pe = -1
contains
  procedure :: init => reader_init
  procedure :: enqueue_1d => reader_enqueue_1d
  procedure :: enqueue_2d => reader_enqueue_2d
  procedure :: enqueue_3d => reader_enqueue_3d
  procedure :: enqueue_4d => reader_enqueue_4d
  generic   :: enqueue => enqueue_1d, enqueue_2d, enqueue_3d, enqueue_4d
  procedure :: commit => reader_commit
end type soca_io_reader

contains

!==============================================================================
! init: prepare a writer for a specific file. The domain pointer is stored, so
! the caller must keep the domain alive until commit returns. Optional root_pe
! overrides the default mpp_root_pe() target for gather+write; the ensemble
! orchestrator uses this to stride writes across writer PEs.
!==============================================================================
subroutine writer_init(self, domain, filename, root_pe)
  class(soca_io_writer), intent(inout) :: self
  type(domain2D), target, intent(in)   :: domain
  character(len=*),       intent(in)   :: filename
  integer, optional,      intent(in)   :: root_pe

  self%filename = filename
  self%domain => domain

  call mpp_get_compute_domain(self%domain, self%isc, self%iec, self%jsc, self%jec)
  call mpp_get_data_domain   (self%domain, self%isd, self%ied, self%jsd, self%jed)
  call mpp_get_global_domain (self%domain, self%isg, self%ieg, self%jsg, self%jeg)
  self%nx_g = self%ieg - self%isg + 1
  self%ny_g = self%jeg - self%jsg + 1

  if (allocated(self%vars)) deallocate(self%vars)
  allocate(self%vars(64))
  self%nvars = 0

  if (present(root_pe)) then
    self%root_pe = root_pe
  else
    self%root_pe = mpp_root_pe()
  end if
  self%ncid = -1
  self%nx_axes = 0
  self%ny_axes = 0
  self%nz_axes = 0
end subroutine writer_init


!==============================================================================
! enqueue_1d / enqueue_2d / enqueue_3d: register a variable for writing. The
! writer holds a pointer to the caller's buffer (no copy); the caller must
! keep the source array alive and unmutated until commit() returns. The actual
! argument must satisfy the TARGET-association rules (declare allocatables
! with the TARGET attribute at the call site).
! 1D assumed global-on-every-PE (no gather; PE 0 writes directly).
! 2D/3D assumed compute-domain-decomposed (mpp_gather to PE 0; the
! compute-domain slice is extracted from the halo-inclusive caller buffer
! inside commit, one var/level at a time).
!==============================================================================
subroutine writer_enqueue_1d(self, name, src, long_name, units)
  class(soca_io_writer),         intent(inout) :: self
  character(len=*),              intent(in)    :: name
  real(kind=kind_real), target,  intent(in)    :: src(:)
  character(len=*), optional,    intent(in)    :: long_name, units

  call check_buf_1d('writer_enqueue_1d', name, size(src))
  call grow_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name    = name
  self%vars(self%nvars)%ndims   = 1
  self%vars(self%nvars)%nlevels = 1
  self%vars(self%nvars)%data1d => src
  self%vars(self%nvars)%long_name = name
  self%vars(self%nvars)%units     = 'none'
  if (present(long_name)) self%vars(self%nvars)%long_name = long_name
  if (present(units))     self%vars(self%nvars)%units     = units
end subroutine writer_enqueue_1d

! 2D/3D enqueue. Caller passes the whole halo-inclusive array (e.g. self%lon,
! shape (isd:ied, jsd:jed)). We hold a pointer; the compute-slice extraction
! happens lazily in commit so only one var's tile is in flight at a time.
subroutine writer_enqueue_2d(self, name, src, long_name, units)
  class(soca_io_writer),         intent(inout) :: self
  character(len=*),              intent(in)    :: name
  real(kind=kind_real), target,  intent(in)    :: src(:,:)
  character(len=*), optional,    intent(in)    :: long_name, units

  call check_buf_2d('writer_enqueue_2d', name, size(src, 1), size(src, 2), &
                    self%ied - self%isd + 1, self%jed - self%jsd + 1)
  call grow_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name    = name
  self%vars(self%nvars)%ndims   = 2
  self%vars(self%nvars)%nlevels = 1
  self%vars(self%nvars)%data2d => src
  self%vars(self%nvars)%long_name = name
  self%vars(self%nvars)%units     = 'none'
  if (present(long_name)) self%vars(self%nvars)%long_name = long_name
  if (present(units))     self%vars(self%nvars)%units     = units
end subroutine writer_enqueue_2d

subroutine writer_enqueue_3d(self, name, src, long_name, units)
  class(soca_io_writer),         intent(inout) :: self
  character(len=*),              intent(in)    :: name
  real(kind=kind_real), target,  intent(in)    :: src(:,:,:)
  character(len=*), optional,    intent(in)    :: long_name, units

  call check_buf_2d('writer_enqueue_3d', name, size(src, 1), size(src, 2), &
                    self%ied - self%isd + 1, self%jed - self%jsd + 1)
  call check_buf_1d('writer_enqueue_3d (z)', name, size(src, 3))
  call grow_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name    = name
  self%vars(self%nvars)%ndims   = 3
  self%vars(self%nvars)%nlevels = size(src, 3)
  self%vars(self%nvars)%data3d => src
  self%vars(self%nvars)%long_name = name
  self%vars(self%nvars)%units     = 'none'
  if (present(long_name)) self%vars(self%nvars)%long_name = long_name
  if (present(units))     self%vars(self%nvars)%units     = units
end subroutine writer_enqueue_3d

subroutine grow_if_needed(self)
  class(soca_io_writer), intent(inout) :: self
  type(var_entry), allocatable :: tmp(:)
  if (self%nvars < size(self%vars)) return
  allocate(tmp(2 * size(self%vars)))
  tmp(1:self%nvars) = self%vars(1:self%nvars)
  call move_alloc(tmp, self%vars)
end subroutine grow_if_needed


!==============================================================================
! commit: single-shot orchestrator. Sequences the four stages on one writer
! (define -> gather -> write -> close). The original monolithic writer_commit
! was split into stages so the ensemble path (soca_io_writers_commit_ensemble)
! can interleave them across multiple writers: define-all -> gather-all ->
! write-all -> close-all. Single-shot semantics are unchanged.
!==============================================================================
subroutine writer_commit(self)
  class(soca_io_writer), intent(inout) :: self

  call writer_stage_define(self)
  call writer_stage_gather(self)
  call writer_stage_write(self)
  call writer_stage_close(self)
end subroutine writer_commit


!==============================================================================
! Stage 1 (define): build axis tables, then on the writer's root_pe create the
! netCDF file, define every dim/coord-var/data-var, exit define mode, and write
! the Time + axis coord-var data. Across multiple writers with different
! root_pes (ensemble path), the file-create work runs concurrently because each
! writer's root is on a different rank and no MPI is involved.
!==============================================================================
subroutine writer_stage_define(self)
  type(soca_io_writer), intent(inout) :: self

  integer, parameter :: MAX_AXES_PER_DIR = 32
  integer :: v, dom_key, dimid_t, varid_t
  logical :: is_root

  is_root = (mpp_pe() == self%root_pe)

  ! Allocate per-direction unique-axis tables on every PE (cheap; need them on
  ! every PE because find_or_add_axis is invoked outside the is_root branch).
  if (allocated(self%x_axes)) deallocate(self%x_axes)
  if (allocated(self%y_axes)) deallocate(self%y_axes)
  if (allocated(self%z_axes)) deallocate(self%z_axes)
  if (allocated(self%var_x_idx)) deallocate(self%var_x_idx)
  if (allocated(self%var_y_idx)) deallocate(self%var_y_idx)
  if (allocated(self%var_z_idx)) deallocate(self%var_z_idx)
  allocate(self%x_axes(MAX_AXES_PER_DIR), self%y_axes(MAX_AXES_PER_DIR), self%z_axes(MAX_AXES_PER_DIR))
  allocate(self%var_x_idx(self%nvars), self%var_y_idx(self%nvars), self%var_z_idx(self%nvars))
  self%var_x_idx = 0; self%var_y_idx = 0; self%var_z_idx = 0
  self%nx_axes = 0; self%ny_axes = 0; self%nz_axes = 0

  do v = 1, self%nvars
    ! 1D vars: no domain (key=0; global on every PE). 2D/3D vars share self%domain (key=1).
    dom_key = 1
    if (self%vars(v)%ndims == 1) dom_key = 0
    select case (self%vars(v)%ndims)
    case (1)
      call find_or_add_axis(self%x_axes, self%nx_axes, &
                            size(self%vars(v)%data1d), dom_key, self%var_x_idx(v))
    case (2)
      call find_or_add_axis(self%x_axes, self%nx_axes, self%nx_g, dom_key, self%var_x_idx(v))
      call find_or_add_axis(self%y_axes, self%ny_axes, self%ny_g, dom_key, self%var_y_idx(v))
    case (3)
      call find_or_add_axis(self%x_axes, self%nx_axes, self%nx_g, dom_key, self%var_x_idx(v))
      call find_or_add_axis(self%y_axes, self%ny_axes, self%ny_g, dom_key, self%var_y_idx(v))
      call find_or_add_axis(self%z_axes, self%nz_axes, &
                            size(self%vars(v)%data3d, 3), dom_key, self%var_z_idx(v))
    end select
  end do

  if (.not. is_root) return

  if (allocated(self%varids)) deallocate(self%varids)
  allocate(self%varids(self%nvars))

  call ncc(nf90_create(self%filename, &
      ior(NF90_CLOBBER, ior(NF90_NETCDF4, NF90_CLASSIC_MODEL)), self%ncid), &
      'nf90_create '//trim(self%filename))

  call define_axis_dims_and_coords(self%ncid, self%x_axes, self%nx_axes, 'xaxis_', 'X')
  call define_axis_dims_and_coords(self%ncid, self%y_axes, self%ny_axes, 'yaxis_', 'Y')
  call define_axis_dims_and_coords(self%ncid, self%z_axes, self%nz_axes, 'zaxis_', 'Z')

  call ncc(nf90_def_dim(self%ncid, 'Time', NF90_UNLIMITED, dimid_t), 'def_dim Time')
  call ncc(nf90_def_var(self%ncid, 'Time', NF90_DOUBLE, [dimid_t], varid_t), 'def_var Time')
  call ncc(nf90_put_att(self%ncid, varid_t, 'long_name',     'Time'),       'att Time:long_name')
  call ncc(nf90_put_att(self%ncid, varid_t, 'units',         'time level'), 'att Time:units')
  call ncc(nf90_put_att(self%ncid, varid_t, 'cartesian_axis','T'),          'att Time:cartesian_axis')

  do v = 1, self%nvars
    select case (self%vars(v)%ndims)
    case (1)
      call ncc(nf90_def_var(self%ncid, trim(self%vars(v)%name), NF90_DOUBLE, &
          [self%x_axes(self%var_x_idx(v))%dimid, dimid_t], self%varids(v)), &
          'def_var '//trim(self%vars(v)%name))
    case (2)
      call ncc(nf90_def_var(self%ncid, trim(self%vars(v)%name), NF90_DOUBLE, &
          [self%x_axes(self%var_x_idx(v))%dimid, &
           self%y_axes(self%var_y_idx(v))%dimid, dimid_t], &
          self%varids(v)), 'def_var '//trim(self%vars(v)%name))
    case (3)
      call ncc(nf90_def_var(self%ncid, trim(self%vars(v)%name), NF90_DOUBLE, &
          [self%x_axes(self%var_x_idx(v))%dimid, &
           self%y_axes(self%var_y_idx(v))%dimid, &
           self%z_axes(self%var_z_idx(v))%dimid, dimid_t], &
          self%varids(v)), 'def_var '//trim(self%vars(v)%name))
    end select
    call ncc(nf90_put_att(self%ncid, self%varids(v), 'long_name', trim(self%vars(v)%long_name)), &
        'att '//trim(self%vars(v)%name)//':long_name')
    call ncc(nf90_put_att(self%ncid, self%varids(v), 'units',     trim(self%vars(v)%units)), &
        'att '//trim(self%vars(v)%name)//':units')
  end do

  call ncc(nf90_enddef(self%ncid), 'enddef')

  call ncc(nf90_put_var(self%ncid, varid_t, [1.0_kind_real], start=[1], count=[1]), 'put Time')
  call put_axis_coord_data(self%ncid, self%x_axes, self%nx_axes)
  call put_axis_coord_data(self%ncid, self%y_axes, self%ny_axes)
  call put_axis_coord_data(self%ncid, self%z_axes, self%nz_axes)
end subroutine writer_stage_define


!==============================================================================
! Stage 2 (gather): for each 2D/3D var, extract the compute-domain tile from
! the caller's halo-inclusive buffer and mpp_gather it to self%root_pe; the
! gathered global field lands in self%vars(v)%gbuf_2d / gbuf_3d (held until
! stage_write). 1D vars are global on every PE and stage_write reads them
! directly; this stage is a no-op for 1D.
!
! Non-root PEs allocate 1x1 dummies so the recv buffer argument to mpp_gather
! is allocated (Fortran assumed-shape contract).
!==============================================================================
subroutine writer_stage_gather(self)
  type(soca_io_writer), intent(inout) :: self

  integer, allocatable :: pelist(:)
  integer :: v, nx_c, ny_c, i_off, j_off, nlev
  real(kind=kind_real), allocatable :: tile2(:,:), tile3(:,:,:)
  logical :: is_root

  is_root = (mpp_pe() == self%root_pe)
  call mpi_pelist(pelist)

  nx_c  = self%iec - self%isc + 1
  ny_c  = self%jec - self%jsc + 1
  i_off = self%isc - self%isd + 1
  j_off = self%jsc - self%jsd + 1

  do v = 1, self%nvars
    select case (self%vars(v)%ndims)
    case (1)
      ! Nothing to do; data1d is global on every PE.
    case (2)
      if (.not. allocated(tile2)) allocate(tile2(nx_c, ny_c))
      if (is_root) then
        if (.not. allocated(self%vars(v)%gbuf_2d)) &
            allocate(self%vars(v)%gbuf_2d(self%nx_g, self%ny_g))
      else
        if (.not. allocated(self%vars(v)%gbuf_2d)) allocate(self%vars(v)%gbuf_2d(1, 1))
      end if
      tile2 = self%vars(v)%data2d(i_off : i_off + nx_c - 1, &
                                  j_off : j_off + ny_c - 1)
      call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
                      tile2, self%vars(v)%gbuf_2d, is_root)
    case (3)
      nlev = self%vars(v)%nlevels
      if (allocated(tile3)) then
        if (size(tile3, 3) /= nlev) deallocate(tile3)
      end if
      if (.not. allocated(tile3)) allocate(tile3(nx_c, ny_c, nlev))
      if (is_root) then
        if (.not. allocated(self%vars(v)%gbuf_3d)) &
            allocate(self%vars(v)%gbuf_3d(self%nx_g, self%ny_g, nlev))
      else
        if (.not. allocated(self%vars(v)%gbuf_3d)) allocate(self%vars(v)%gbuf_3d(1, 1, 1))
      end if
      tile3 = self%vars(v)%data3d(i_off : i_off + nx_c - 1, &
                                  j_off : j_off + ny_c - 1, :)
      ! Single 3D mpp_gather per 3D var: one collective replaces nlevels 2D
      ! gathers, and root receives the assembled global field directly.
      call mpp_gather(self%isc, self%iec, self%jsc, self%jec, nlev, pelist, &
                      tile3, self%vars(v)%gbuf_3d, is_root)
    end select
  end do

  if (allocated(tile2))  deallocate(tile2)
  if (allocated(tile3))  deallocate(tile3)
  if (allocated(pelist)) deallocate(pelist)
end subroutine writer_stage_gather


!==============================================================================
! Stage 3 (write): on self%root_pe, dump each var to disk via nf90_put_var.
! For 1D vars the source is data1d directly (global on every PE; root just
! reads its own copy). For 2D/3D the source is the gbuf_2d / gbuf_3d that
! stage_gather assembled. No MPI -> across multiple writers with different
! root_pes (ensemble path) this stage runs concurrently per PE.
!==============================================================================
subroutine writer_stage_write(self)
  type(soca_io_writer), intent(inout) :: self

  integer :: v, nlev
  logical :: is_root

  is_root = (mpp_pe() == self%root_pe)
  if (.not. is_root) return

  do v = 1, self%nvars
    select case (self%vars(v)%ndims)
    case (1)
      call ncc(nf90_put_var(self%ncid, self%varids(v), self%vars(v)%data1d, &
          start=[1, 1], count=[size(self%vars(v)%data1d), 1]), &
          'put '//trim(self%vars(v)%name))
    case (2)
      call ncc(nf90_put_var(self%ncid, self%varids(v), self%vars(v)%gbuf_2d, &
          start=[1, 1, 1], count=[self%nx_g, self%ny_g, 1]), &
          'put '//trim(self%vars(v)%name))
    case (3)
      nlev = self%vars(v)%nlevels
      call ncc(nf90_put_var(self%ncid, self%varids(v), self%vars(v)%gbuf_3d, &
          start=[1, 1, 1, 1], count=[self%nx_g, self%ny_g, nlev, 1]), &
          'put '//trim(self%vars(v)%name))
    end select
  end do
end subroutine writer_stage_write


!==============================================================================
! Stage 4 (close): on self%root_pe close the file; everyone deallocates the
! inter-stage buffers and the writer's working set. No MPI -> concurrent across
! writers in the ensemble path.
!==============================================================================
subroutine writer_stage_close(self)
  type(soca_io_writer), intent(inout) :: self

  integer :: v
  logical :: is_root

  is_root = (mpp_pe() == self%root_pe)

  if (is_root .and. self%ncid >= 0) then
    call ncc(nf90_close(self%ncid), 'nf90_close')
  end if
  self%ncid = -1

  if (allocated(self%varids))    deallocate(self%varids)
  if (allocated(self%x_axes))    deallocate(self%x_axes)
  if (allocated(self%y_axes))    deallocate(self%y_axes)
  if (allocated(self%z_axes))    deallocate(self%z_axes)
  if (allocated(self%var_x_idx)) deallocate(self%var_x_idx)
  if (allocated(self%var_y_idx)) deallocate(self%var_y_idx)
  if (allocated(self%var_z_idx)) deallocate(self%var_z_idx)
  self%nx_axes = 0; self%ny_axes = 0; self%nz_axes = 0

  if (allocated(self%vars)) then
    do v = 1, self%nvars
      if (allocated(self%vars(v)%gbuf_2d)) deallocate(self%vars(v)%gbuf_2d)
      if (allocated(self%vars(v)%gbuf_3d)) deallocate(self%vars(v)%gbuf_3d)
    end do
    deallocate(self%vars)
  end if
  self%nvars = 0
end subroutine writer_stage_close


!==============================================================================
! Bulk-write orchestrator: take an array of pre-configured writers (each one
! init'd against the same domain with its own filename and rotated root_pe)
! and run all four stages across the whole set.
!
! Phase ordering:
!   1. define-all  (each writer's root_pe creates its file; concurrent, no MPI)
!   2. gather-all  (sync: per-writer per-var mpp_gather -- world collectives
!                   sequenced; async: all (writer, var) MPI_Igatherv batched
!                   into one wave + MPI_Waitall so the runtime can overlap)
!   3. write-all   (each writer's root_pe nf90_put_var's its data; concurrent)
!   4. close-all   (concurrent nf90_close)
!
! No inter-phase barriers are needed: the gather phase is self-synchronizing
! (collective + MPI_Waitall), and define/write/close are independent per
! writer root_pe, so each rank runs them at its own pace.
!
! Memory: each writer's root_pe holds one member's worth of every var's
! gbuf_2d/gbuf_3d between phases 2 and 3 -- peak ~global state size per writer.
!==============================================================================
subroutine soca_io_writers_commit_ensemble(writers)
  type(soca_io_writer), intent(inout) :: writers(:)

  integer :: m

  if (size(writers) == 0) return

  ! All writers in a batch must share one FMS domain: the async gather
  ! allgathers writers(1)'s compute-domain bounds once and reuses them to unpack
  ! every member's tiles. A mismatched domain would silently corrupt the field.
  do m = 2, size(writers)
    if (.not. associated(writers(m)%domain, writers(1)%domain)) &
        call mpp_error(FATAL, &
            'soca_io_writers_commit_ensemble: all writers in a batch must share one domain')
  end do

  ! Phase 1: file create + dim/var defs. Concurrent across root_pes; no MPI.
  do m = 1, size(writers)
    call writer_stage_define(writers(m))
  end do

  ! Phase 2: gather. Sync path issues per-writer per-var mpp_gather (each is a
  ! world-comm collective, so they serialize across writers). Async path
  ! batches every (writer, var) MPI_Igatherv into one wave + MPI_Waitall.
  if (async_mpi_enabled_) then
    call writer_stage_gather_all_async(writers)
  else
    do m = 1, size(writers)
      call writer_stage_gather(writers(m))
    end do
  end if

  ! Phase 3: each root_pe writes its writer's data. No MPI -> concurrent disk.
  do m = 1, size(writers)
    call writer_stage_write(writers(m))
  end do

  ! Phase 4: each root_pe closes its file. No MPI -> concurrent.
  do m = 1, size(writers)
    call writer_stage_close(writers(m))
  end do
end subroutine soca_io_writers_commit_ensemble


!==============================================================================
! Async cross-writer gather: pack every PE's compute-domain tile for every
! (writer, 2D-or-3D var) into a per-(request) send buffer, post one
! MPI_Igatherv per (writer, var) targeted at writers(m)%root_pe, MPI_Waitall,
! then each root unpacks its received tiles into the var's gbuf. The runtime
! can overlap Igathers whose roots differ; sync mpp_gather would force one
! world-comm collective per gather, serializing the whole batch.
!
! Each 3D var ships as a single block (i fastest, then j, then k); 2D vars get
! one Igatherv each. Per-sender memory cost is the sum of (writers, vars) of
! compute-tile-size * nz_var bytes -- the whole batch in flight at once.
!==============================================================================
subroutine writer_stage_gather_all_async(writers)
  type(soca_io_writer), intent(inout) :: writers(:)

  integer :: nprocs, ierr, r, m, v, idx, rq_count, max_reqs, req_idx
  integer :: my_isc, my_iec, my_jsc, my_jec, my_count2d
  integer :: nz_v, full_count, isg, jsg, nx_g, ny_g
  integer :: i, j, k, ig, jg, kk, mp_comm
  integer, allocatable :: all_isc(:), all_iec(:), all_jsc(:), all_jec(:)
  integer, allocatable :: tile_counts2d(:)
  integer, allocatable :: reqs(:)
  integer, allocatable :: req_writer(:), req_var(:)
  ! Per-request count/displacement arrays. MPI_Igatherv is nonblocking, so the
  ! count/displ arguments must stay valid (unmutated) until MPI_Waitall. Each
  ! request has its own column rc(:,rq)/disp(:,rq) -- a single shared rc(:)
  ! recomputed each loop iteration would be aliased across in-flight gathers.
  integer, allocatable :: rc(:,:), disp(:,:)
  type :: gather_buf_t
    real(kind=kind_real), allocatable :: buf(:)
  end type
  type(gather_buf_t), allocatable :: sends(:), recvs(:)
  logical :: is_root_for_this

  if (size(writers) == 0) return
  call current_mpi_comm(mp_comm)
  nprocs = mpp_npes()

  ! All writers in the ensemble share the same FMS domain (the geometry's),
  ! so gather the compute-domain bounds once and reuse for every Igatherv.
  call mpp_get_compute_domain(writers(1)%domain, my_isc, my_iec, my_jsc, my_jec)
  call gather_compute_domains(writers(1)%domain, mp_comm, all_isc, all_iec, all_jsc, all_jec)

  allocate(tile_counts2d(nprocs))
  do r = 1, nprocs
    tile_counts2d(r) = (all_iec(r) - all_isc(r) + 1) * (all_jec(r) - all_jsc(r) + 1)
  end do
  my_count2d = (my_iec - my_isc + 1) * (my_jec - my_jsc + 1)

  ! Allocate the gbufs on each writer's root_pe up front (one per 2D/3D var).
  ! These buffers persist through stage_write and are freed in stage_close.
  do m = 1, size(writers)
    is_root_for_this = (mpp_pe() == writers(m)%root_pe)
    if (.not. is_root_for_this) cycle
    do v = 1, writers(m)%nvars
      select case (writers(m)%vars(v)%ndims)
      case (2)
        if (.not. allocated(writers(m)%vars(v)%gbuf_2d)) &
            allocate(writers(m)%vars(v)%gbuf_2d(writers(m)%nx_g, writers(m)%ny_g))
      case (3)
        if (.not. allocated(writers(m)%vars(v)%gbuf_3d)) &
            allocate(writers(m)%vars(v)%gbuf_3d( &
                writers(m)%nx_g, writers(m)%ny_g, writers(m)%vars(v)%nlevels))
      end select
    end do
  end do

  ! Count requests: one per (writer, gatherable var). 1D vars contribute no
  ! request (data1d is global on every PE).
  max_reqs = 0
  do m = 1, size(writers)
    do v = 1, writers(m)%nvars
      if (writers(m)%vars(v)%ndims >= 2) max_reqs = max_reqs + 1
    end do
  end do
  if (max_reqs == 0) goto 9999

  allocate(reqs(max_reqs))
  allocate(req_writer(max_reqs), req_var(max_reqs))
  allocate(sends(max_reqs), recvs(max_reqs))
  allocate(rc(nprocs, max_reqs), disp(nprocs, max_reqs))
  rq_count = 0

  ! Phase A/B: pack send + post Igatherv per (writer, var).
  do m = 1, size(writers)
    is_root_for_this = (mpp_pe() == writers(m)%root_pe)
    do v = 1, writers(m)%nvars
      if (writers(m)%vars(v)%ndims < 2) cycle
      rq_count = rq_count + 1
      req_writer(rq_count) = m
      req_var(rq_count)    = v

      if (writers(m)%vars(v)%ndims == 2) then
        nz_v = 1
      else
        nz_v = writers(m)%vars(v)%nlevels
      end if

      do r = 1, nprocs
        rc(r, rq_count) = tile_counts2d(r) * nz_v
      end do
      disp(1, rq_count) = 0
      do r = 2, nprocs
        disp(r, rq_count) = disp(r-1, rq_count) + rc(r-1, rq_count)
      end do

      allocate(sends(rq_count)%buf(my_count2d * nz_v))
      idx = 0
      if (writers(m)%vars(v)%ndims == 2) then
        do j = 1, my_jec - my_jsc + 1
          do i = 1, my_iec - my_isc + 1
            idx = idx + 1
            sends(rq_count)%buf(idx) = writers(m)%vars(v)%data2d( &
                writers(m)%isc - writers(m)%isd + i, &
                writers(m)%jsc - writers(m)%jsd + j)
          end do
        end do
      else
        do k = 1, nz_v
          do j = 1, my_jec - my_jsc + 1
            do i = 1, my_iec - my_isc + 1
              idx = idx + 1
              sends(rq_count)%buf(idx) = writers(m)%vars(v)%data3d( &
                  writers(m)%isc - writers(m)%isd + i, &
                  writers(m)%jsc - writers(m)%jsd + j, k)
            end do
          end do
        end do
      end if

      if (is_root_for_this) then
        full_count = 0
        do r = 1, nprocs
          full_count = full_count + rc(r, rq_count)
        end do
        allocate(recvs(rq_count)%buf(full_count))
      else
        ! Non-root recvbuf is not significant to MPI_Igatherv, but the actual
        ! argument must still be an allocated array (1-elem dummy mirrors the
        ! reader scatter path's sendbuf handling).
        allocate(recvs(rq_count)%buf(1))
      end if

      call MPI_Igatherv(sends(rq_count)%buf, my_count2d * nz_v, MPI_DOUBLE_PRECISION, &
          recvs(rq_count)%buf, rc(:, rq_count), disp(:, rq_count), MPI_DOUBLE_PRECISION, &
          writers(m)%root_pe, mp_comm, reqs(rq_count), ierr)
    end do
  end do

  ! Phase C: wait for everything.
  call MPI_Waitall(rq_count, reqs, MPI_STATUSES_IGNORE, ierr)

  ! Phase D: each root unpacks its tiles into the per-var gbuf. The pack loop
  ! above wrote each sender's tile contiguously; Igatherv lays them back-to-
  ! back at disp(r) in the recvbuf in rank order.
  do req_idx = 1, rq_count
    m = req_writer(req_idx)
    v = req_var(req_idx)
    if (mpp_pe() /= writers(m)%root_pe) cycle
    isg = writers(m)%isg
    jsg = writers(m)%jsg
    nx_g = writers(m)%nx_g
    ny_g = writers(m)%ny_g
    if (writers(m)%vars(v)%ndims == 2) then
      idx = 0
      do r = 1, nprocs
        do jg = all_jsc(r), all_jec(r)
          do ig = all_isc(r), all_iec(r)
            idx = idx + 1
            writers(m)%vars(v)%gbuf_2d(ig - isg + 1, jg - jsg + 1) = recvs(req_idx)%buf(idx)
          end do
        end do
      end do
    else
      idx = 0
      do r = 1, nprocs
        do kk = 1, writers(m)%vars(v)%nlevels
          do jg = all_jsc(r), all_jec(r)
            do ig = all_isc(r), all_iec(r)
              idx = idx + 1
              writers(m)%vars(v)%gbuf_3d(ig - isg + 1, jg - jsg + 1, kk) = recvs(req_idx)%buf(idx)
            end do
          end do
        end do
      end do
    end if
  end do

  do i = 1, rq_count
    if (allocated(sends(i)%buf)) deallocate(sends(i)%buf)
    if (allocated(recvs(i)%buf)) deallocate(recvs(i)%buf)
  end do
  deallocate(reqs, req_writer, req_var, sends, recvs, rc, disp)

9999 continue
  deallocate(all_isc, all_iec, all_jsc, all_jec, tile_counts2d)
end subroutine writer_stage_gather_all_async


!==============================================================================
! Reader: init / enqueue_* / commit. Caller buffer stays in place; enqueue
! records a pointer, commit fills the compute-domain interior. Halos are left
! untouched -- the caller refreshes them via mpp_update_domains (same as FMS).
!
! Optional reader_pe sets the source PE for the scatter path (defaults to
! mpp_root_pe() in single-state mode); the ensemble orchestrator pre-sets a
! rotated reader_pe per member so different files are read concurrently.
!==============================================================================
subroutine reader_init(self, domain, filename, reader_pe)
  class(soca_io_reader), intent(inout) :: self
  type(domain2D), target, intent(in)   :: domain
  character(len=*),       intent(in)   :: filename
  integer, optional,      intent(in)   :: reader_pe

  self%filename = filename
  self%domain => domain

  call mpp_get_compute_domain(self%domain, self%isc, self%iec, self%jsc, self%jec)
  call mpp_get_data_domain   (self%domain, self%isd, self%ied, self%jsd, self%jed)
  call mpp_get_global_domain (self%domain, self%isg, self%ieg, self%jsg, self%jeg)
  self%nx_g = self%ieg - self%isg + 1
  self%ny_g = self%jeg - self%jsg + 1

  if (allocated(self%vars)) deallocate(self%vars)
  allocate(self%vars(64))
  self%nvars = 0

  if (present(reader_pe)) then
    self%reader_pe = reader_pe
  else
    self%reader_pe = mpp_root_pe()
  end if
end subroutine reader_init


! 1D vars are global on every PE (no scatter): PE 0 reads, broadcasts.
subroutine reader_enqueue_1d(self, name, dst)
  class(soca_io_reader),                  intent(inout) :: self
  character(len=*),                       intent(in)    :: name
  real(kind=kind_real), target,           intent(inout) :: dst(:)

  call check_buf_1d('reader_enqueue_1d', name, size(dst))
  call grow_reader_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name  = name
  self%vars(self%nvars)%ndims = 1
  self%vars(self%nvars)%data1d => dst
end subroutine reader_enqueue_1d


! 2D vars are domain-decomposed: each PE pulls its own compute slice into the
! caller's data-domain buffer. Caller buffer must be sized exactly
! (data_xsize, data_ysize); halos are left undefined for mpp_update_domains.
subroutine reader_enqueue_2d(self, name, dst)
  class(soca_io_reader),                  intent(inout) :: self
  character(len=*),                       intent(in)    :: name
  real(kind=kind_real), target,           intent(inout) :: dst(:,:)

  call check_buf_2d('reader_enqueue_2d', name, size(dst, 1), size(dst, 2), &
                    self%ied - self%isd + 1, self%jed - self%jsd + 1)
  call grow_reader_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name  = name
  self%vars(self%nvars)%ndims = 2
  self%vars(self%nvars)%data2d => dst
end subroutine reader_enqueue_2d


! 3D buffer: (data_xsize, data_ysize, nlevels). Spatial dims are
! domain-decomposed; the third dim is held entire on every PE.
subroutine reader_enqueue_3d(self, name, dst)
  class(soca_io_reader),                  intent(inout) :: self
  character(len=*),                       intent(in)    :: name
  real(kind=kind_real), target,           intent(inout) :: dst(:,:,:)

  call check_buf_2d('reader_enqueue_3d', name, size(dst, 1), size(dst, 2), &
                    self%ied - self%isd + 1, self%jed - self%jsd + 1)
  call check_buf_1d('reader_enqueue_3d (z)', name, size(dst, 3))
  call grow_reader_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name  = name
  self%vars(self%nvars)%ndims = 3
  self%vars(self%nvars)%data3d => dst
end subroutine reader_enqueue_3d


! 4D buffer: (data_xsize, data_ysize, n3, n4). Spatial dims are
! domain-decomposed; the trailing two dims are held entire on every PE.
! Used for CICE category+level fields (Tsnz_h etc.).
subroutine reader_enqueue_4d(self, name, dst)
  class(soca_io_reader),                  intent(inout) :: self
  character(len=*),                       intent(in)    :: name
  real(kind=kind_real), target,           intent(inout) :: dst(:,:,:,:)

  call check_buf_2d('reader_enqueue_4d', name, size(dst, 1), size(dst, 2), &
                    self%ied - self%isd + 1, self%jed - self%jsd + 1)
  call check_buf_1d('reader_enqueue_4d (n3)', name, size(dst, 3))
  call check_buf_1d('reader_enqueue_4d (n4)', name, size(dst, 4))
  call grow_reader_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name  = name
  self%vars(self%nvars)%ndims = 4
  self%vars(self%nvars)%data4d => dst
end subroutine reader_enqueue_4d


subroutine grow_reader_if_needed(self)
  class(soca_io_reader), intent(inout) :: self
  type(read_entry), allocatable :: tmp(:)
  if (self%nvars < size(self%vars)) return
  allocate(tmp(2 * size(self%vars)))
  tmp(1:self%nvars) = self%vars(1:self%nvars)
  call move_alloc(tmp, self%vars)
end subroutine grow_reader_if_needed


!==============================================================================
! Read all enqueued vars. Single-state path dispatches on the YAML knob
! geometry.io.'single state read': 'strided' (default, every PE reads its own
! tile direct) or 'scatter' (reader_pe reads the whole field then MPI_Scatterv
! to tile owners). The strided path mirrors FMS MPP_READ_2DDECOMP -- no PE-0
! bottleneck, N parallel reads. The scatter path concentrates the disk read
! on one PE: useful when 200-500 PEs concurrently striding the same file
! stresses the filesystem.
!==============================================================================
subroutine reader_commit(self)
  class(soca_io_reader), intent(inout) :: self

  if (single_state_read_scatter_) then
    call reader_stage_read(self)
    call reader_stage_distribute(self)
    call reader_stage_close(self)
  else
    call commit_reader_strided(self)
    ! release pointers; caller's buffers untouched (they hold the read data)
    if (allocated(self%vars)) deallocate(self%vars)
    self%nvars = 0
  end if
end subroutine reader_commit


!==============================================================================
! Per-PE strided read implementation. Opens the file, pulls each var's
! compute-domain tile via nf90_get_var(start, count), closes. Stateless --
! no module-level handle cache.
!==============================================================================
subroutine commit_reader_strided(self)
  class(soca_io_reader), intent(inout) :: self

  integer :: ncid, v, n3, n4
  integer :: nx_c, ny_c, i_off, j_off, i_start, j_start
  real(kind=kind_real), allocatable :: tile2(:,:), tile3(:,:,:), tile4(:,:,:,:)

  call ncc(nf90_open(self%filename, NF90_NOWRITE, ncid), &
      'nf90_open '//trim(self%filename))

  nx_c    = self%iec - self%isc + 1
  ny_c    = self%jec - self%jsc + 1
  i_off   = self%isc - self%isd + 1     ! 1-based offset into data-domain buf
  j_off   = self%jsc - self%jsd + 1
  i_start = self%isc - self%isg + 1     ! 1-based start in the global file dim
  j_start = self%jsc - self%jsg + 1

  ! tile2 size is invariant across vars; tile3/tile4 are reallocated only when
  ! the trailing dims change between vars (rare in typical state I/O).
  do v = 1, self%nvars
    select case (self%vars(v)%ndims)
    case (1)
      call read_var_strided(ncid, self%vars(v)%name, &
          1, 1, size(self%vars(v)%data1d), 1, dst1=self%vars(v)%data1d)

    case (2)
      if (.not. allocated(tile2)) allocate(tile2(nx_c, ny_c))
      call read_var_strided(ncid, self%vars(v)%name, &
          i_start, j_start, nx_c, ny_c, dst2=tile2)
      self%vars(v)%data2d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1) = tile2

    case (3)
      n3 = size(self%vars(v)%data3d, 3)
      if (allocated(tile3)) then
        if (size(tile3, 3) /= n3) deallocate(tile3)
      end if
      if (.not. allocated(tile3)) allocate(tile3(nx_c, ny_c, n3))
      call read_var_strided(ncid, self%vars(v)%name, &
          i_start, j_start, nx_c, ny_c, dst3=tile3)
      self%vars(v)%data3d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1, :) = tile3

    case (4)
      n3 = size(self%vars(v)%data4d, 3)
      n4 = size(self%vars(v)%data4d, 4)
      if (allocated(tile4)) then
        if (size(tile4, 3) /= n3 .or. size(tile4, 4) /= n4) deallocate(tile4)
      end if
      if (.not. allocated(tile4)) allocate(tile4(nx_c, ny_c, n3, n4))
      call read_var_strided(ncid, self%vars(v)%name, &
          i_start, j_start, nx_c, ny_c, dst4=tile4)
      self%vars(v)%data4d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1, :, :) = tile4
    end select
  end do

  call ncc(nf90_close(ncid), 'nf90_close '//trim(self%filename))

  if (allocated(tile2)) deallocate(tile2)
  if (allocated(tile3)) deallocate(tile3)
  if (allocated(tile4)) deallocate(tile4)
end subroutine commit_reader_strided


!==============================================================================
! Stage 1 (read): on self%reader_pe, pull the whole global field for every
! enqueued var into a per-var staged buffer. No MPI -> across multiple readers
! with different reader_pes (ensemble path) this stage runs concurrently per
! PE (independent files, independent disk I/O).
!
! 1D vars: read directly into data1d on reader_pe (reader_stage_distribute
! broadcasts it). 2D/3D/4D: alloc gbuf_{2,3,4}d on reader_pe; non-reader PEs
! hold nothing (the scatter recvbuf is allocated on the fly in distribute).
!==============================================================================
subroutine reader_stage_read(self)
  type(soca_io_reader), intent(inout) :: self

  integer :: ncid, v, n, n3, n4
  logical :: is_root

  is_root = (mpp_pe() == self%reader_pe)
  if (.not. is_root) return

  call ncc(nf90_open(self%filename, NF90_NOWRITE, ncid), &
      'nf90_open '//trim(self%filename))

  do v = 1, self%nvars
    select case (self%vars(v)%ndims)
    case (1)
      n = size(self%vars(v)%data1d)
      call read_var_strided(ncid, self%vars(v)%name, &
          1, 1, n, 1, dst1=self%vars(v)%data1d)
    case (2)
      if (.not. allocated(self%vars(v)%gbuf_2d)) &
          allocate(self%vars(v)%gbuf_2d(self%nx_g, self%ny_g))
      call read_var_strided(ncid, self%vars(v)%name, &
          1, 1, self%nx_g, self%ny_g, dst2=self%vars(v)%gbuf_2d)
    case (3)
      n3 = size(self%vars(v)%data3d, 3)
      if (.not. allocated(self%vars(v)%gbuf_3d)) &
          allocate(self%vars(v)%gbuf_3d(self%nx_g, self%ny_g, n3))
      call read_var_strided(ncid, self%vars(v)%name, &
          1, 1, self%nx_g, self%ny_g, dst3=self%vars(v)%gbuf_3d)
    case (4)
      n3 = size(self%vars(v)%data4d, 3)
      n4 = size(self%vars(v)%data4d, 4)
      if (.not. allocated(self%vars(v)%gbuf_4d)) &
          allocate(self%vars(v)%gbuf_4d(self%nx_g, self%ny_g, n3, n4))
      call read_var_strided(ncid, self%vars(v)%name, &
          1, 1, self%nx_g, self%ny_g, dst4=self%vars(v)%gbuf_4d)
    end select
  end do

  call ncc(nf90_close(ncid), 'nf90_close '//trim(self%filename))
end subroutine reader_stage_read


!==============================================================================
! Stage 2 (distribute, sync): per-var send the relevant slice from reader_pe
! to every PE in the world communicator. 1D vars are MPI_Bcast. 2D/3D/4D vars
! are packed rank-ordered on reader_pe and shipped via MPI_Scatterv; each
! receiver unpacks its compute tile into data2d/data3d/data4d.
!
! Per-var Allgather of the compute-domain bounds: needed because the pack
! order on reader_pe must match the unpack order on every receiver. The
! Allgather happens once per (reader, var); for the async batched path the
! same bounds are reused across all vars (see reader_stage_distribute_all_async).
!==============================================================================
subroutine reader_stage_distribute(self)
  type(soca_io_reader), intent(inout) :: self

  integer :: nprocs, ierr, v, r, i, j, k, k4, n3, n4, idx, full_count
  integer :: my_isc, my_iec, my_jsc, my_jec, my_count2d
  integer :: nx_c, ny_c, i_off, j_off, mp_comm
  integer, allocatable :: all_isc(:), all_iec(:), all_jsc(:), all_jec(:)
  integer, allocatable :: tile_counts2d(:), rc(:), disp(:)
  real(kind=kind_real), allocatable :: sendbuf(:), recvbuf(:)
  logical :: is_root

  is_root = (mpp_pe() == self%reader_pe)
  nprocs  = mpp_npes()
  call current_mpi_comm(mp_comm)

  call mpp_get_compute_domain(self%domain, my_isc, my_iec, my_jsc, my_jec)
  call gather_compute_domains(self%domain, mp_comm, all_isc, all_iec, all_jsc, all_jec)

  allocate(tile_counts2d(nprocs))
  do r = 1, nprocs
    tile_counts2d(r) = (all_iec(r) - all_isc(r) + 1) * (all_jec(r) - all_jsc(r) + 1)
  end do
  my_count2d = (my_iec - my_isc + 1) * (my_jec - my_jsc + 1)

  nx_c  = self%iec - self%isc + 1
  ny_c  = self%jec - self%jsc + 1
  i_off = self%isc - self%isd + 1
  j_off = self%jsc - self%jsd + 1

  allocate(rc(nprocs), disp(nprocs))

  do v = 1, self%nvars
    select case (self%vars(v)%ndims)
    case (1)
      ! 1D: every PE wants the same global vector. Bcast directly into data1d.
      call MPI_Bcast(self%vars(v)%data1d, size(self%vars(v)%data1d), &
          MPI_DOUBLE_PRECISION, self%reader_pe, mp_comm, ierr)

    case (2, 3, 4)
      if (self%vars(v)%ndims == 2) then
        n3 = 1; n4 = 1
      else if (self%vars(v)%ndims == 3) then
        n3 = size(self%vars(v)%data3d, 3); n4 = 1
      else
        n3 = size(self%vars(v)%data4d, 3); n4 = size(self%vars(v)%data4d, 4)
      end if

      do r = 1, nprocs
        rc(r) = tile_counts2d(r) * n3 * n4
      end do
      disp(1) = 0
      do r = 2, nprocs
        disp(r) = disp(r-1) + rc(r-1)
      end do

      if (is_root) then
        full_count = 0
        do r = 1, nprocs
          full_count = full_count + rc(r)
        end do
        allocate(sendbuf(full_count))
        idx = 0
        ! Pack rank-by-rank, k4 fastest-outer, then k, then j, then i (matches unpack)
        do r = 1, nprocs
          do k4 = 1, n4
            do k = 1, n3
              do j = all_jsc(r), all_jec(r)
                do i = all_isc(r), all_iec(r)
                  idx = idx + 1
                  select case (self%vars(v)%ndims)
                  case (2)
                    sendbuf(idx) = self%vars(v)%gbuf_2d(i - self%isg + 1, j - self%jsg + 1)
                  case (3)
                    sendbuf(idx) = self%vars(v)%gbuf_3d(i - self%isg + 1, j - self%jsg + 1, k)
                  case (4)
                    sendbuf(idx) = self%vars(v)%gbuf_4d(i - self%isg + 1, j - self%jsg + 1, k, k4)
                  end select
                end do
              end do
            end do
          end do
        end do
      else
        allocate(sendbuf(1))
      end if

      allocate(recvbuf(my_count2d * n3 * n4))
      call MPI_Scatterv(sendbuf, rc, disp, MPI_DOUBLE_PRECISION, &
          recvbuf, my_count2d * n3 * n4, MPI_DOUBLE_PRECISION, &
          self%reader_pe, mp_comm, ierr)

      ! Unpack recvbuf into the caller's data buffer compute slice. Halos
      ! untouched (caller refreshes via mpp_update_domains).
      idx = 0
      do k4 = 1, n4
        do k = 1, n3
          do j = 1, ny_c
            do i = 1, nx_c
              idx = idx + 1
              select case (self%vars(v)%ndims)
              case (2)
                self%vars(v)%data2d(i_off + i - 1, j_off + j - 1) = recvbuf(idx)
              case (3)
                self%vars(v)%data3d(i_off + i - 1, j_off + j - 1, k) = recvbuf(idx)
              case (4)
                self%vars(v)%data4d(i_off + i - 1, j_off + j - 1, k, k4) = recvbuf(idx)
              end select
            end do
          end do
        end do
      end do

      deallocate(sendbuf, recvbuf)
    end select
  end do

  deallocate(all_isc, all_iec, all_jsc, all_jec, tile_counts2d, rc, disp)
end subroutine reader_stage_distribute


!==============================================================================
! Stage 3 (close): drop the staged global buffers and the reader's working set.
!==============================================================================
subroutine reader_stage_close(self)
  type(soca_io_reader), intent(inout) :: self
  integer :: v
  if (allocated(self%vars)) then
    do v = 1, self%nvars
      if (allocated(self%vars(v)%gbuf_2d)) deallocate(self%vars(v)%gbuf_2d)
      if (allocated(self%vars(v)%gbuf_3d)) deallocate(self%vars(v)%gbuf_3d)
      if (allocated(self%vars(v)%gbuf_4d)) deallocate(self%vars(v)%gbuf_4d)
    end do
    deallocate(self%vars)
  end if
  self%nvars = 0
end subroutine reader_stage_close


!==============================================================================
! Bulk-read orchestrator. Two modes governed by the geometry.io knob
! 'ensemble read':
!   - strided: loop readers, each one does commit_reader_strided (every PE
!     reads its compute-domain tile direct from that member's file). All
!     members done sequentially; per-member read is parallel across PEs.
!   - scatter: read-all -> distribute-all. Phase 1: every reader_pe reads its
!     member's whole global field into staged buffers, no MPI -> reader PEs
!     hit different files concurrently. Phase 2: per-reader MPI_Scatterv from
!     reader_pe (sync: per-reader collectives serialize; async: all
!     (reader, var) Iscatterv + Waitall batched).
!==============================================================================
subroutine soca_io_readers_commit_ensemble(readers)
  type(soca_io_reader), intent(inout) :: readers(:)

  integer :: m

  if (size(readers) == 0) return

  ! All readers in a batch must share one FMS domain (see the writer
  ! orchestrator note): the async scatter allgathers readers(1)'s
  ! compute-domain bounds once and reuses them for every member.
  do m = 2, size(readers)
    if (.not. associated(readers(m)%domain, readers(1)%domain)) &
        call mpp_error(FATAL, &
            'soca_io_readers_commit_ensemble: all readers in a batch must share one domain')
  end do

  if (.not. ensemble_read_scatter_) then
    ! Strided ensemble path: loop members, run the existing per-PE strided
    ! reader for each. Members are sequenced; within each member every PE
    ! reads its own tile.
    do m = 1, size(readers)
      call commit_reader_strided(readers(m))
      if (allocated(readers(m)%vars)) deallocate(readers(m)%vars)
      readers(m)%nvars = 0
    end do
    return
  end if

  ! Scatter ensemble path.
  ! Phase 1: each reader_pe reads its file. No MPI -> reads concurrent.
  do m = 1, size(readers)
    call reader_stage_read(readers(m))
  end do

  ! Phase 2: scatter. Sync per-reader serializes; async batches all
  ! (reader, var) Iscatterv/Ibcast + Waitall.
  if (async_mpi_enabled_) then
    call reader_stage_distribute_all_async(readers)
  else
    do m = 1, size(readers)
      call reader_stage_distribute(readers(m))
    end do
  end if

  ! Phase 3: free staged buffers + drop working set.
  do m = 1, size(readers)
    call reader_stage_close(readers(m))
  end do
end subroutine soca_io_readers_commit_ensemble


!==============================================================================
! Async cross-reader scatter: per (reader, var) pack the relevant rank-ordered
! slices on the reader_pe, post one MPI_Iscatterv per var (MPI_Ibcast for 1D),
! MPI_Waitall, then each receiver unpacks its recvbuf into the caller buffer
! compute slice. Sync mpp/MPI_Scatterv serializes across readers because each
! is a world-comm collective; async batches them so the runtime can overlap
! scatters whose roots differ.
!==============================================================================
subroutine reader_stage_distribute_all_async(readers)
  type(soca_io_reader), intent(inout) :: readers(:)

  integer :: nprocs, ierr, r, m, v, idx, full_count
  integer :: my_isc, my_iec, my_jsc, my_jec, my_count2d, mp_comm
  integer :: n3, n4, i, j, k, k4, nx_c, ny_c, i_off, j_off
  integer, allocatable :: all_isc(:), all_iec(:), all_jsc(:), all_jec(:)
  integer, allocatable :: tile_counts2d(:)
  ! Per-request count/displacement columns: MPI_Iscatterv is nonblocking, so
  ! these must stay valid until MPI_Waitall (one column per in-flight request).
  integer, allocatable :: rc(:,:), disp(:,:)
  integer :: rq_count, max_reqs, ux
  integer, allocatable :: reqs(:), req_reader(:), req_var(:)
  type :: scatter_buf_t
    real(kind=kind_real), allocatable :: buf(:)
  end type
  type(scatter_buf_t), allocatable :: sends(:), recvs(:)
  logical :: is_root_for_this

  if (size(readers) == 0) return
  call current_mpi_comm(mp_comm)
  nprocs = mpp_npes()

  call mpp_get_compute_domain(readers(1)%domain, my_isc, my_iec, my_jsc, my_jec)
  call gather_compute_domains(readers(1)%domain, mp_comm, all_isc, all_iec, all_jsc, all_jec)

  allocate(tile_counts2d(nprocs))
  do r = 1, nprocs
    tile_counts2d(r) = (all_iec(r) - all_isc(r) + 1) * (all_jec(r) - all_jsc(r) + 1)
  end do
  my_count2d = (my_iec - my_isc + 1) * (my_jec - my_jsc + 1)

  ! Count requests: one per (reader, var). 1D and 2D/3D/4D all count.
  max_reqs = 0
  do m = 1, size(readers)
    max_reqs = max_reqs + readers(m)%nvars
  end do
  if (max_reqs == 0) goto 9999

  allocate(reqs(max_reqs), req_reader(max_reqs), req_var(max_reqs))
  allocate(sends(max_reqs), recvs(max_reqs))
  allocate(rc(nprocs, max_reqs), disp(nprocs, max_reqs))
  rq_count = 0

  ! Phase A: pack send buffers + post Iscatterv / Ibcast.
  do m = 1, size(readers)
    is_root_for_this = (mpp_pe() == readers(m)%reader_pe)
    do v = 1, readers(m)%nvars
      rq_count = rq_count + 1
      req_reader(rq_count) = m
      req_var(rq_count)    = v

      if (readers(m)%vars(v)%ndims == 1) then
        ! 1D: post a single MPI_Ibcast of data1d (every PE wants the same).
        call MPI_Ibcast(readers(m)%vars(v)%data1d, size(readers(m)%vars(v)%data1d), &
            MPI_DOUBLE_PRECISION, readers(m)%reader_pe, mp_comm, reqs(rq_count), ierr)
        cycle
      end if

      if (readers(m)%vars(v)%ndims == 2) then
        n3 = 1; n4 = 1
      else if (readers(m)%vars(v)%ndims == 3) then
        n3 = size(readers(m)%vars(v)%data3d, 3); n4 = 1
      else
        n3 = size(readers(m)%vars(v)%data4d, 3); n4 = size(readers(m)%vars(v)%data4d, 4)
      end if

      do r = 1, nprocs
        rc(r, rq_count) = tile_counts2d(r) * n3 * n4
      end do
      disp(1, rq_count) = 0
      do r = 2, nprocs
        disp(r, rq_count) = disp(r-1, rq_count) + rc(r-1, rq_count)
      end do

      if (is_root_for_this) then
        full_count = 0
        do r = 1, nprocs
          full_count = full_count + rc(r, rq_count)
        end do
        allocate(sends(rq_count)%buf(full_count))
        idx = 0
        do r = 1, nprocs
          do k4 = 1, n4
            do k = 1, n3
              do j = all_jsc(r), all_jec(r)
                do i = all_isc(r), all_iec(r)
                  idx = idx + 1
                  select case (readers(m)%vars(v)%ndims)
                  case (2)
                    sends(rq_count)%buf(idx) = readers(m)%vars(v)%gbuf_2d( &
                        i - readers(m)%isg + 1, j - readers(m)%jsg + 1)
                  case (3)
                    sends(rq_count)%buf(idx) = readers(m)%vars(v)%gbuf_3d( &
                        i - readers(m)%isg + 1, j - readers(m)%jsg + 1, k)
                  case (4)
                    sends(rq_count)%buf(idx) = readers(m)%vars(v)%gbuf_4d( &
                        i - readers(m)%isg + 1, j - readers(m)%jsg + 1, k, k4)
                  end select
                end do
              end do
            end do
          end do
        end do
      else
        allocate(sends(rq_count)%buf(1))
      end if

      allocate(recvs(rq_count)%buf(my_count2d * n3 * n4))
      call MPI_Iscatterv(sends(rq_count)%buf, rc(:, rq_count), disp(:, rq_count), MPI_DOUBLE_PRECISION, &
          recvs(rq_count)%buf, my_count2d * n3 * n4, MPI_DOUBLE_PRECISION, &
          readers(m)%reader_pe, mp_comm, reqs(rq_count), ierr)
    end do
  end do

  ! Phase B: wait.
  call MPI_Waitall(rq_count, reqs, MPI_STATUSES_IGNORE, ierr)

  ! Phase C: unpack recvbufs into caller compute slices (1D was bcast in-place).
  do idx = 1, rq_count
    m = req_reader(idx)
    v = req_var(idx)
    if (readers(m)%vars(v)%ndims == 1) cycle

    if (readers(m)%vars(v)%ndims == 2) then
      n3 = 1; n4 = 1
    else if (readers(m)%vars(v)%ndims == 3) then
      n3 = size(readers(m)%vars(v)%data3d, 3); n4 = 1
    else
      n3 = size(readers(m)%vars(v)%data4d, 3); n4 = size(readers(m)%vars(v)%data4d, 4)
    end if

    nx_c  = readers(m)%iec - readers(m)%isc + 1
    ny_c  = readers(m)%jec - readers(m)%jsc + 1
    i_off = readers(m)%isc - readers(m)%isd + 1
    j_off = readers(m)%jsc - readers(m)%jsd + 1

    ux = 0
    do k4 = 1, n4
      do k = 1, n3
        do j = 1, ny_c
          do i = 1, nx_c
            ux = ux + 1
            select case (readers(m)%vars(v)%ndims)
            case (2)
              readers(m)%vars(v)%data2d(i_off + i - 1, j_off + j - 1) = recvs(idx)%buf(ux)
            case (3)
              readers(m)%vars(v)%data3d(i_off + i - 1, j_off + j - 1, k) = recvs(idx)%buf(ux)
            case (4)
              readers(m)%vars(v)%data4d(i_off + i - 1, j_off + j - 1, k, k4) = recvs(idx)%buf(ux)
            end select
          end do
        end do
      end do
    end do
  end do

  do idx = 1, size(sends)
    if (allocated(sends(idx)%buf)) deallocate(sends(idx)%buf)
    if (allocated(recvs(idx)%buf)) deallocate(recvs(idx)%buf)
  end do
  deallocate(reqs, req_reader, req_var, sends, recvs, rc, disp)

9999 continue
  deallocate(all_isc, all_iec, all_jsc, all_jec, tile_counts2d)
end subroutine reader_stage_distribute_all_async




!==============================================================================
! Find an existing axis matching (size, domain_key) or append a new one.
! Sets `idx` to the (1-based) position in the axes() array.
!==============================================================================
subroutine find_or_add_axis(axes, n_axes, size_, dom_key, idx)
  type(axis_entry), intent(inout) :: axes(:)
  integer,          intent(inout) :: n_axes
  integer,          intent(in)    :: size_, dom_key
  integer,          intent(out)   :: idx
  integer :: i

  do i = 1, n_axes
    if (axes(i)%size == size_ .and. axes(i)%domain_key == dom_key) then
      idx = i
      return
    end if
  end do

  if (n_axes >= size(axes)) then
    call mpp_error(FATAL, 'soca_io_mod: too many unique axes (raise MAX_AXES_PER_DIR)')
  end if
  n_axes = n_axes + 1
  axes(n_axes)%size       = size_
  axes(n_axes)%domain_key = dom_key
  axes(n_axes)%dimid      = -1
  axes(n_axes)%varid      = -1
  idx = n_axes
end subroutine find_or_add_axis


!==============================================================================
! Define dims and 1D coordinate variables for one axis direction (X/Y/Z).
! Names are <prefix><N> where N is the 1-based axis index, matching FMS.
!==============================================================================
subroutine define_axis_dims_and_coords(ncid, axes, n_axes, prefix, cart_axis)
  integer,             intent(in)    :: ncid
  type(axis_entry),    intent(inout) :: axes(:)
  integer,             intent(in)    :: n_axes
  character(len=*),    intent(in)    :: prefix    ! 'xaxis_' / 'yaxis_' / 'zaxis_'
  character(len=1),    intent(in)    :: cart_axis ! 'X' / 'Y' / 'Z'
  integer :: j
  character(len=16) :: name

  do j = 1, n_axes
    if (j < 10) then
      write(name, '(A,I1)') trim(prefix), j
    else
      write(name, '(A,I2)') trim(prefix), j
    end if
    call ncc(nf90_def_dim(ncid, trim(name), axes(j)%size, axes(j)%dimid), &
        'def_dim '//trim(name))
    call ncc(nf90_def_var(ncid, trim(name), NF90_DOUBLE, [axes(j)%dimid], axes(j)%varid), &
        'def_var '//trim(name))
    call ncc(nf90_put_att(ncid, axes(j)%varid, 'long_name',      trim(name)),    'att '//trim(name))
    call ncc(nf90_put_att(ncid, axes(j)%varid, 'units',          'none'),        'att '//trim(name))
    call ncc(nf90_put_att(ncid, axes(j)%varid, 'cartesian_axis', cart_axis),     'att '//trim(name))
  end do
end subroutine define_axis_dims_and_coords


!==============================================================================
! Write the index-sequence values 1.0, 2.0, ..., size(j) for each coord var.
!==============================================================================
subroutine put_axis_coord_data(ncid, axes, n_axes)
  integer,             intent(in)    :: ncid
  type(axis_entry),    intent(inout) :: axes(:)
  integer,             intent(in)    :: n_axes
  integer :: j, i
  real(kind=kind_real), allocatable :: idxbuf(:)

  do j = 1, n_axes
    if (allocated(idxbuf)) then
      if (size(idxbuf) /= axes(j)%size) deallocate(idxbuf)
    end if
    if (.not. allocated(idxbuf)) allocate(idxbuf(axes(j)%size))
    do i = 1, axes(j)%size
      idxbuf(i) = real(i, kind=kind_real)
    end do
    call ncc(nf90_put_var(ncid, axes(j)%varid, idxbuf), 'put coord var')
  end do
  if (allocated(idxbuf)) deallocate(idxbuf)
end subroutine put_axis_coord_data


!==============================================================================
! Replacements for FMS file_exist / field_exist (FMS I/O metadata helpers).
!==============================================================================
logical function soca_io_file_exists(filename)
  character(len=*), intent(in) :: filename
  inquire(file=trim(filename), exist=soca_io_file_exists)
end function soca_io_file_exists

logical function soca_io_var_exists(filename, varname)
  character(len=*), intent(in) :: filename, varname
  integer :: ncid, varid, status
  soca_io_var_exists = .false.
  if (.not. soca_io_file_exists(filename)) return
  status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
  if (status /= NF90_NOERR) return
  status = nf90_inq_varid(ncid, trim(varname), varid)
  soca_io_var_exists = (status == NF90_NOERR)
  status = nf90_close(ncid)
end function soca_io_var_exists


!==============================================================================
! Strided read of one variable into a caller-owned tile. Builds start/count
! from (i_start, j_start, nx, ny) plus the file's actual dim sizes, tolerating
! degenerate size-1 dims the buffer doesn't enumerate (same middle-dim trick
! FMS register_restart_field + restore_state did transparently). For 2D+:
! count(1)=nx, count(2)=ny, count(ndims)=1 (trailing time), middle entries from
! middle_dims; total element count is preserved so netcdf-fortran reads into
! the smaller-rank buffer correctly. 1D (dst1): count(1)=nx, j_start/ny ignored.
! Concrete cases:
!   - 3D state Salt(time, Layer, lath, lonh) -> 3D tile (nx, ny, Layer)
!     file_ndims=4, count=[nx, ny, Layer, 1].
!   - 5D CICE Tsnz_h(time, nc=5, nksnow=1, nj, ni) -> 3D tile (nx, ny, 5)
!     file_ndims=5, count=[nx, ny, 1, 5, 1] = [nx, ny, nksnow, nc, time].
!==============================================================================
subroutine read_var_strided(ncid, name, i_start, j_start, nx, ny, &
                            dst1, dst2, dst3, dst4)
  integer,                       intent(in)  :: ncid
  integer,                       intent(in)  :: i_start, j_start, nx, ny
  character(len=*),              intent(in)  :: name
  real(kind=kind_real), optional, intent(out) :: dst1(:)
  real(kind=kind_real), optional, intent(out) :: dst2(:,:)
  real(kind=kind_real), optional, intent(out) :: dst3(:,:,:)
  real(kind=kind_real), optional, intent(out) :: dst4(:,:,:,:)

  integer :: varid, file_ndims, dd, sz, dst_rank
  integer :: file_dimids(MAX_FILE_NDIMS)
  integer :: st(MAX_FILE_NDIMS), ct(MAX_FILE_NDIMS)
  integer :: total_ct, expected_total
  character(len=8) :: total_str

  call ncc(nf90_inq_varid(ncid, trim(name), varid), 'inq '//trim(name))
  call ncc(nf90_inquire_variable(ncid, varid, ndims=file_ndims, dimids=file_dimids), &
           'inquire '//trim(name))
  if (file_ndims > MAX_FILE_NDIMS) call mpp_error(FATAL, &
      'soca_io_mod read_var_strided: '//trim(name)//' exceeds MAX_FILE_NDIMS')

  if (present(dst1)) then
    dst_rank = 1
  else if (present(dst2)) then
    dst_rank = 2
  else if (present(dst3)) then
    dst_rank = 3
  else if (present(dst4)) then
    dst_rank = 4
  else
    call mpp_error(FATAL, 'soca_io_mod read_var_strided: no destination provided')
  end if
  if (file_ndims < dst_rank) call mpp_error(FATAL, &
      'soca_io_mod read_var_strided: '//trim(name)//' has fewer file dims than destination rank')

  st(1:file_ndims) = 1
  ct(1:file_ndims) = 1
  if (dst_rank == 1) then
    ct(1) = nx
  else
    st(1) = i_start
    st(2) = j_start
    ct(1) = nx
    ct(2) = ny
    if (file_ndims == dst_rank) then
      ! All file dims are spatial (no trailing time/record). Fill from file.
      do dd = 3, file_ndims
        call ncc(nf90_inquire_dimension(ncid, file_dimids(dd), len=sz), &
                 'inq dim '//trim(name))
        ct(dd) = sz
      end do
    else
      ! file_ndims > dst_rank: one trailing time/record dim plus any squeezable
      ! middle dims (e.g. CICE Tsnz_h's nksnow=1 between (ni,nj) and nc).
      ct(file_ndims) = 1
      do dd = 3, file_ndims - 1
        call ncc(nf90_inquire_dimension(ncid, file_dimids(dd), len=sz), &
                 'inq dim '//trim(name))
        ct(dd) = sz
      end do
    end if
    ! Total-count check: catches both per-dim size mismatches (e.g. file
    ! z=75 vs destination z=25) and the silent partial-fill that motivated
    ! this rewrite.
    total_ct = 1
    do dd = 1, file_ndims
      total_ct = total_ct * ct(dd)
    end do
    expected_total = nx * ny
    if (present(dst3)) expected_total = expected_total * size(dst3, 3)
    if (present(dst4)) expected_total = expected_total * size(dst4, 3) * size(dst4, 4)
    if (total_ct /= expected_total) then
      write(total_str, '(i0)') expected_total
      call mpp_error(FATAL, 'soca_io_mod read_var_strided: '//trim(name)// &
          ' file element count does not match destination ('//trim(total_str)//')')
    end if
  end if

  if (present(dst2)) then
    call ncc(nf90_get_var(ncid, varid, dst2, start=st(1:file_ndims), count=ct(1:file_ndims)), &
        'get '//trim(name))
  else if (present(dst3)) then
    call ncc(nf90_get_var(ncid, varid, dst3, start=st(1:file_ndims), count=ct(1:file_ndims)), &
        'get '//trim(name))
  else if (present(dst4)) then
    call ncc(nf90_get_var(ncid, varid, dst4, start=st(1:file_ndims), count=ct(1:file_ndims)), &
        'get '//trim(name))
  else if (present(dst1)) then
    call ncc(nf90_get_var(ncid, varid, dst1, start=st(1:file_ndims), count=ct(1:file_ndims)), &
        'get '//trim(name))
  end if
end subroutine read_var_strided


!==============================================================================
! Fetch mpp's current pelist (= geometry's f_comm pelist) for mpp_gather. Using
! MPI_COMM_WORLD's pelist is wrong in ensemble mode where each task has its
! own size-1 mpp world.
!==============================================================================
subroutine mpi_pelist(pelist)
  integer, allocatable, intent(out) :: pelist(:)
  allocate(pelist(mpp_npes()))
  call mpp_get_current_pelist(pelist)
end subroutine mpi_pelist


!==============================================================================
! Fetch the MPI communicator handle that the current mpp pelist sits on. Used
! by the async ensemble I/O paths so MPI_I{gatherv,scatterv,bcast} run on the
! same comm as mpp's collectives (== the geometry's f_comm). Using
! MPI_COMM_WORLD here would deadlock in the LETKF "nens per MPI task" split
! where each task has its own size-1 mpp world.
!==============================================================================
subroutine current_mpi_comm(mp_comm)
  integer, intent(out) :: mp_comm
  integer, allocatable :: pelist(:)
  character(len=128) :: name
  ! mpp_get_current_pelist requires a pelist out-arg sized mpp_npes(); the async
  ! I/O callers only need the communicator handle, so it stays local here.
  allocate(pelist(mpp_npes()))
  call mpp_get_current_pelist(pelist, name, mp_comm)
  deallocate(pelist)
end subroutine current_mpi_comm


!==============================================================================
! Allgather every PE's compute-domain bounds (isc/iec/jsc/jec) in mp_comm rank
! order. One packed collective replaces the four separate Allgathers the bulk
! read/write paths used to issue; results are indexed by rank, so all_isc(r) is
! rank (r-1)'s compute-domain start -- matching the rank-ordered Igatherv /
! Iscatterv buffer layout. The bounds are static, but the gather is cheap and
! avoids any cross-call caching hazard with multiple geometries.
!==============================================================================
subroutine gather_compute_domains(domain, mp_comm, all_isc, all_iec, all_jsc, all_jec)
  type(domain2D),       intent(in)  :: domain
  integer,              intent(in)  :: mp_comm
  integer, allocatable, intent(out) :: all_isc(:), all_iec(:), all_jsc(:), all_jec(:)
  integer :: nprocs, ierr, my(4)
  integer, allocatable :: gathered(:)

  nprocs = mpp_npes()
  call mpp_get_compute_domain(domain, my(1), my(2), my(3), my(4))
  allocate(gathered(4 * nprocs))
  call MPI_Allgather(my, 4, MPI_INTEGER, gathered, 4, MPI_INTEGER, mp_comm, ierr)
  allocate(all_isc(nprocs), all_iec(nprocs), all_jsc(nprocs), all_jec(nprocs))
  all_isc = gathered(1::4)
  all_iec = gathered(2::4)
  all_jsc = gathered(3::4)
  all_jec = gathered(4::4)
  deallocate(gathered)
end subroutine gather_compute_domains


!==============================================================================
! Netcdf error check. Aborts on error with mpp_error(FATAL, ...).
!==============================================================================
subroutine ncc(status, where)
  integer,           intent(in) :: status
  character(len=*),  intent(in) :: where
  if (status /= NF90_NOERR) then
    call mpp_error(FATAL, 'soca_io_mod ['//trim(where)//']: '//trim(nf90_strerror(status)))
  end if
end subroutine ncc


!==============================================================================
! Resolve the module-level I/O dispatch knobs from the geometry YAML. Looked-up
! keys (all optional, defaults retained on miss):
!     io:
!       ensemble write: parallel | sequential   (default parallel)
!       ensemble read:  scatter  | strided      (default scatter)
!       single state read: strided | scatter    (default strided)
!       async mpi: true | false                 (default true)
! Called once from soca_geom_init with the geometry's fckit_configuration; the
! resolved values persist module-level for the rest of the run.
!==============================================================================
subroutine soca_io_config_from_yaml(f_conf)
  type(fckit_configuration), intent(in) :: f_conf
  type(fckit_configuration) :: io_conf
  character(len=:), allocatable :: sval
  logical :: lval, ok
  integer :: ival

  if (.not. f_conf%has("io")) return
  ok = f_conf%get("io", io_conf)
  if (.not. ok) return

  if (io_conf%has("ensemble write")) then
    ok = io_conf%get("ensemble write", sval)
    if (ok) then
      select case (trim(sval))
      case ("parallel");   ensemble_write_parallel_ = .true.
      case ("sequential"); ensemble_write_parallel_ = .false.
      case default
        call mpp_error(FATAL, "soca_io_mod: geometry.io.'ensemble write' must be " // &
            "'parallel' or 'sequential', got '" // trim(sval) // "'")
      end select
    end if
  end if

  if (io_conf%has("ensemble read")) then
    ok = io_conf%get("ensemble read", sval)
    if (ok) then
      select case (trim(sval))
      case ("scatter"); ensemble_read_scatter_ = .true.
      case ("strided"); ensemble_read_scatter_ = .false.
      case default
        call mpp_error(FATAL, "soca_io_mod: geometry.io.'ensemble read' must be " // &
            "'scatter' or 'strided', got '" // trim(sval) // "'")
      end select
    end if
  end if

  if (io_conf%has("single state read")) then
    ok = io_conf%get("single state read", sval)
    if (ok) then
      select case (trim(sval))
      case ("scatter"); single_state_read_scatter_ = .true.
      case ("strided"); single_state_read_scatter_ = .false.
      case default
        call mpp_error(FATAL, "soca_io_mod: geometry.io.'single state read' must be " // &
            "'scatter' or 'strided', got '" // trim(sval) // "'")
      end select
    end if
  end if

  if (io_conf%has("async mpi")) then
    ok = io_conf%get("async mpi", lval)
    if (ok) async_mpi_enabled_ = lval
  end if

  if (io_conf%has("ensemble batch size")) then
    ok = io_conf%get("ensemble batch size", ival)
    if (ok) then
      if (ival < 0) then
        call mpp_error(FATAL, "soca_io_mod: geometry.io.'ensemble batch size' must be " // &
            ">= 0 (0 = single batch over all members)")
      end if
      ensemble_batch_size_ = ival
    end if
  end if
end subroutine soca_io_config_from_yaml


!==============================================================================
! Public getter for the ensemble-write dispatch knob. Consumed by fields_mod
! when it picks between the single-shot and parallel ensemble write paths.
! (ensemble_read_scatter_ / single_state_read_scatter_ / async_mpi_enabled_ are
! read directly within this module, so they need no public getters.)
!==============================================================================
logical function soca_io_ensemble_write_parallel()
  soca_io_ensemble_write_parallel = ensemble_write_parallel_
end function soca_io_ensemble_write_parallel

! Ensemble I/O batch size knob (0 = single batch over all members).
integer function soca_io_ensemble_batch_size()
  soca_io_ensemble_batch_size = ensemble_batch_size_
end function soca_io_ensemble_batch_size


!==============================================================================
! Strided writer/reader-PE rotation: map a 1-based member index m to a rank in
! [0, mpp_npes()-1] using floor((m-1) * npes / nmembers). Spreads ensemble
! members across the world communicator so concurrent disk I/O does not
! collide on PE 0; when nmembers > npes the surplus members wrap onto earlier
! PEs in stride. Reused by the writer and reader ensemble orchestrators.
!
! Examples (npes=8):
!   nmembers=4  -> roots [0, 2, 4, 6]
!   nmembers=8  -> roots [0..7]
!   nmembers=20 -> roots [0,0,0,1,1,2,2,2,3,3,4,4,4,5,5,6,6,6,7,7]
!==============================================================================
function soca_io_ensemble_root_pe(m, nmembers) result(pe)
  integer, intent(in) :: m
  integer, intent(in) :: nmembers
  integer :: pe
  pe = (m - 1) * mpp_npes() / nmembers
end function soca_io_ensemble_root_pe


!==============================================================================
! Defensive buffer-shape checks, called from each enqueue. Catches an
! unallocated actual (size 0) and a caller that mis-sized the data-domain
! buffer (e.g. used compute extents instead of data-domain extents). Pure
! integer comparisons -- no perf cost.
!==============================================================================
subroutine check_buf_1d(role, name, n)
  character(len=*), intent(in) :: role, name
  integer,          intent(in) :: n
  character(len=32) :: ns
  if (n <= 0) then
    write(ns, '(I0)') n
    call mpp_error(FATAL, 'soca_io_mod '//trim(role)//': buffer for "'// &
        trim(name)//'" has non-positive size '//trim(ns)// &
        ' (unallocated or empty?)')
  end if
end subroutine check_buf_1d

subroutine check_buf_2d(role, name, n1, n2, expected_n1, expected_n2)
  character(len=*), intent(in) :: role, name
  integer,          intent(in) :: n1, n2, expected_n1, expected_n2
  character(len=128) :: s
  if (n1 /= expected_n1 .or. n2 /= expected_n2) then
    write(s, '(A,I0,A,I0,A,I0,A,I0,A)') &
        ' got (', n1, ',', n2, '); expected data-domain (', &
        expected_n1, ',', expected_n2, ')'
    call mpp_error(FATAL, 'soca_io_mod '//trim(role)//': buffer for "'// &
        trim(name)//'" wrong shape -- '//trim(s))
  end if
end subroutine check_buf_2d

end module soca_io_mod
