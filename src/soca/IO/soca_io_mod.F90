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
use kinds, only: kind_real
use mpp_mod, only: mpp_gather, mpp_scatter, mpp_broadcast, mpp_pe, mpp_root_pe, mpp_npes, &
                   mpp_get_current_pelist, mpp_error, FATAL
use mpp_domains_mod, only: domain2D, &
                           mpp_get_compute_domain, mpp_get_global_domain, &
                           mpp_get_data_domain

implicit none
private

public :: soca_io_writer
public :: soca_io_reader
public :: soca_io_file_exists, soca_io_var_exists

integer, parameter :: MAX_NAME       = 256
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
! the caller must keep the domain alive until commit returns.
!==============================================================================
subroutine writer_init(self, domain, filename)
  class(soca_io_writer), intent(inout) :: self
  type(domain2D), target, intent(in)   :: domain
  character(len=*),       intent(in)   :: filename

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
! commit: PE 0 creates the file structure, then each var is mpp_gather'd and
! PE 0 writes via nf90_put_var. Equivalent to FMS mpp_io threading=MPP_SINGLE
! -- the goal is a clean, debuggable, FMS-free baseline, not a speedup.
!==============================================================================
subroutine writer_commit(self)
  class(soca_io_writer), intent(inout) :: self

  integer :: ncid, dimid_t, varid_t, v
  integer, allocatable :: varids(:)
  integer, allocatable :: var_x_idx(:), var_y_idx(:), var_z_idx(:)
  type(axis_entry), allocatable :: x_axes(:), y_axes(:), z_axes(:)
  integer :: nx_axes, ny_axes, nz_axes
  logical :: is_root
  real(kind=kind_real), allocatable :: gbuf2d(:,:), gbuf3d(:,:,:)
  integer, allocatable :: pelist(:)
  integer :: dom_key
  integer, parameter :: MAX_AXES_PER_DIR = 32
  integer :: nx_c, ny_c, i_off, j_off, nlev
  real(kind=kind_real), allocatable :: tile2(:,:), tile3(:,:,:)

  is_root = (mpp_pe() == mpp_root_pe())
  call mpi_pelist(pelist)

  nx_c  = self%iec - self%isc + 1
  ny_c  = self%jec - self%jsc + 1
  i_off = self%isc - self%isd + 1
  j_off = self%jsc - self%jsd + 1

  ! Build per-direction unique-axis tables (FMS algorithm): match each var's
  ! dim against existing axes by (size, domain_key); reuse on match, append on
  ! miss.
  allocate(x_axes(MAX_AXES_PER_DIR), y_axes(MAX_AXES_PER_DIR), z_axes(MAX_AXES_PER_DIR))
  allocate(var_x_idx(self%nvars), var_y_idx(self%nvars), var_z_idx(self%nvars))
  var_x_idx = 0; var_y_idx = 0; var_z_idx = 0
  nx_axes = 0; ny_axes = 0; nz_axes = 0

  do v = 1, self%nvars
    ! 1D vars: no domain (global on every PE, key=0). 2D/3D vars share
    ! self%domain (key=1).
    dom_key = 1
    if (self%vars(v)%ndims == 1) dom_key = 0

    ! Every var contributes an X axis (Fortran first dim). For 2D/3D the local
    ! buffer is only the compute slice, so use the writer's GLOBAL extents for
    ! the axis size; 1D buffers are already global.
    select case (self%vars(v)%ndims)
    case (1)
      call find_or_add_axis(x_axes, nx_axes, size(self%vars(v)%data1d), dom_key, var_x_idx(v))
    case (2)
      call find_or_add_axis(x_axes, nx_axes, self%nx_g, dom_key, var_x_idx(v))
      call find_or_add_axis(y_axes, ny_axes, self%ny_g, dom_key, var_y_idx(v))
    case (3)
      call find_or_add_axis(x_axes, nx_axes, self%nx_g, dom_key, var_x_idx(v))
      call find_or_add_axis(y_axes, ny_axes, self%ny_g, dom_key, var_y_idx(v))
      call find_or_add_axis(z_axes, nz_axes, size(self%vars(v)%data3d, 3), dom_key, var_z_idx(v))
    end select
  end do

  ! Phase 1: PE 0 defines the file structure -- dims, coord vars, data vars.
  if (is_root) then
    allocate(varids(self%nvars))

    call ncc(nf90_create(self%filename, &
        ior(NF90_CLOBBER, ior(NF90_NETCDF4, NF90_CLASSIC_MODEL)), ncid), &
        'nf90_create '//trim(self%filename))

    call define_axis_dims_and_coords(ncid, x_axes, nx_axes, 'xaxis_', 'X')
    call define_axis_dims_and_coords(ncid, y_axes, ny_axes, 'yaxis_', 'Y')
    call define_axis_dims_and_coords(ncid, z_axes, nz_axes, 'zaxis_', 'Z')

    call ncc(nf90_def_dim(ncid, 'Time', NF90_UNLIMITED, dimid_t), 'def_dim Time')
    call ncc(nf90_def_var(ncid, 'Time', NF90_DOUBLE, [dimid_t], varid_t), 'def_var Time')
    call ncc(nf90_put_att(ncid, varid_t, 'long_name',     'Time'),       'att Time:long_name')
    call ncc(nf90_put_att(ncid, varid_t, 'units',         'time level'), 'att Time:units')
    call ncc(nf90_put_att(ncid, varid_t, 'cartesian_axis','T'),          'att Time:cartesian_axis')

    do v = 1, self%nvars
      select case (self%vars(v)%ndims)
      case (1)
        ! Fortran dim list: [xaxis, Time] -> file order (Time, xaxis)
        call ncc(nf90_def_var(ncid, trim(self%vars(v)%name), NF90_DOUBLE, &
            [x_axes(var_x_idx(v))%dimid, dimid_t], varids(v)), &
            'def_var '//trim(self%vars(v)%name))
      case (2)
        ! Fortran [xaxis, yaxis, Time] -> file (Time, yaxis, xaxis)
        call ncc(nf90_def_var(ncid, trim(self%vars(v)%name), NF90_DOUBLE, &
            [x_axes(var_x_idx(v))%dimid, y_axes(var_y_idx(v))%dimid, dimid_t], &
            varids(v)), 'def_var '//trim(self%vars(v)%name))
      case (3)
        ! Fortran [xaxis, yaxis, zaxis, Time] -> file (Time, zaxis, yaxis, xaxis)
        call ncc(nf90_def_var(ncid, trim(self%vars(v)%name), NF90_DOUBLE, &
            [x_axes(var_x_idx(v))%dimid, y_axes(var_y_idx(v))%dimid, &
             z_axes(var_z_idx(v))%dimid, dimid_t], varids(v)), &
            'def_var '//trim(self%vars(v)%name))
      end select
      call ncc(nf90_put_att(ncid, varids(v), 'long_name', trim(self%vars(v)%long_name)), &
          'att '//trim(self%vars(v)%name)//':long_name')
      call ncc(nf90_put_att(ncid, varids(v), 'units',     trim(self%vars(v)%units)), &
          'att '//trim(self%vars(v)%name)//':units')
    end do

    call ncc(nf90_enddef(ncid), 'enddef')

    ! Coordinate-var data is just the index sequence 1..size, matching FMS.
    call ncc(nf90_put_var(ncid, varid_t, [1.0_kind_real], start=[1], count=[1]), 'put Time')
    call put_axis_coord_data(ncid, x_axes, nx_axes)
    call put_axis_coord_data(ncid, y_axes, ny_axes)
    call put_axis_coord_data(ncid, z_axes, nz_axes)
  end if

  ! Phase 2: gather and write each user variable. The compute-domain tile is
  ! extracted from the caller's halo-inclusive buffer one var (and, for 3D,
  ! one level) at a time, so peak local memory is a single (nx_c, ny_c) tile.
  ! Non-root allocates a 1x1 dummy for gbuf2d so the actual argument to
  ! mpp_gather is always allocated (assumed-shape dummy requires it).
  if (is_root) then
    allocate(gbuf2d(self%nx_g, self%ny_g))
  else
    allocate(gbuf2d(1, 1))
  end if
  allocate(tile2(nx_c, ny_c))

  do v = 1, self%nvars
    if (self%vars(v)%ndims == 1) then
      if (is_root) then
        call ncc(nf90_put_var(ncid, varids(v), self%vars(v)%data1d, &
            start=[1, 1], count=[size(self%vars(v)%data1d), 1]), &
            'put '//trim(self%vars(v)%name))
      end if
    else if (self%vars(v)%ndims == 2) then
      tile2 = self%vars(v)%data2d(i_off : i_off + nx_c - 1, &
                                  j_off : j_off + ny_c - 1)
      call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
                      tile2, gbuf2d, is_root)
      if (is_root) then
        call ncc(nf90_put_var(ncid, varids(v), gbuf2d, &
            start=[1, 1, 1], count=[self%nx_g, self%ny_g, 1]), &
            'put '//trim(self%vars(v)%name))
      end if
    else
      ! Single 3D mpp_gather per 3D var: one collective replaces nlevels 2D
      ! gathers, and root receives the assembled global field directly into
      ! gbuf3d (no per-level gbuf2d->gbuf3d memcpy). Reuse tile3/gbuf3d across
      ! 3D vars when nlevels matches (typical for ocean state).
      nlev = self%vars(v)%nlevels
      if (is_root) then
        if (allocated(gbuf3d)) then
          if (size(gbuf3d, 3) /= nlev) deallocate(gbuf3d)
        end if
        if (.not. allocated(gbuf3d)) allocate(gbuf3d(self%nx_g, self%ny_g, nlev))
      else
        ! Non-root dummy so the actual argument to mpp_gather is allocated.
        if (.not. allocated(gbuf3d)) allocate(gbuf3d(1, 1, 1))
      end if
      if (allocated(tile3)) then
        if (size(tile3, 3) /= nlev) deallocate(tile3)
      end if
      if (.not. allocated(tile3)) allocate(tile3(nx_c, ny_c, nlev))
      tile3 = self%vars(v)%data3d(i_off : i_off + nx_c - 1, &
                                  j_off : j_off + ny_c - 1, :)
      call mpp_gather(self%isc, self%iec, self%jsc, self%jec, nlev, pelist, &
                      tile3, gbuf3d, is_root)
      if (is_root) then
        call ncc(nf90_put_var(ncid, varids(v), gbuf3d, &
            start=[1, 1, 1, 1], count=[self%nx_g, self%ny_g, nlev, 1]), &
            'put '//trim(self%vars(v)%name))
      end if
    end if
  end do

  ! Phase 3: close.
  if (is_root) then
    call ncc(nf90_close(ncid), 'nf90_close')
    deallocate(varids)
  end if
  deallocate(gbuf2d, tile2, pelist, x_axes, y_axes, z_axes, var_x_idx, var_y_idx, var_z_idx)
  if (allocated(gbuf3d)) deallocate(gbuf3d)
  if (allocated(tile3)) deallocate(tile3)

  ! drop pointer entries; caller's data is unaffected
  if (allocated(self%vars)) deallocate(self%vars)
  self%nvars = 0
end subroutine writer_commit


!==============================================================================
! Reader: init / enqueue_* / commit. Caller buffer stays in place; enqueue
! records a pointer, commit fills the compute-domain interior. Halos are left
! untouched -- the caller refreshes them via mpp_update_domains (same as FMS).
!==============================================================================
subroutine reader_init(self, domain, filename)
  class(soca_io_reader), intent(inout) :: self
  type(domain2D), target, intent(in)   :: domain
  character(len=*),       intent(in)   :: filename

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
! Read all enqueued vars. Every PE opens the file NF90_NOWRITE and pulls only
! its compute-domain tile via nf90_get_var(start, count) -- mirrors FMS's
! MPP_READ_2DDECOMP: no PE-0 bottleneck, no mpp_broadcast, N parallel reads.
! Classic / 64-bit-offset netcdf allows concurrent read-only opens; library
! state is process-local. 1D vars also read independently on every PE.
!==============================================================================
subroutine reader_commit(self)
  class(soca_io_reader), intent(inout) :: self

  call commit_reader_strided(self)

  ! release pointers; caller's buffers untouched (they hold the read data)
  if (allocated(self%vars)) deallocate(self%vars)
  self%nvars = 0
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
! PE-0 read + per-PE scatter. PE 0 nf90_get_var's the global field; mpp_scatter
! sends each PE its compute-domain slice (1D vars are broadcast). Mirrors FMS
! 2024.02 fms_netcdf_domain_io.F90:domain_read_3d.
!
! TODO: currently unused -- single-state reads use commit_reader_strided. Will
! be exercised by parallel-ensemble I/O, where one reader PE per member
! scatters that member to its compute-PE group.
!==============================================================================
subroutine commit_reader_scatter(self)
  class(soca_io_reader), intent(inout) :: self

  integer :: ncid, v, n, n3, n4, k4
  integer :: nx_c, ny_c, i_off, j_off
  integer :: is_f, ie_f, js_f, je_f  ! PE tile in 1-based file-space indices
  integer, allocatable :: pelist(:)
  real(kind=kind_real), allocatable :: gbuf2(:,:), gbuf3(:,:,:), gbuf4(:,:,:,:)
  real(kind=kind_real), allocatable :: tile2(:,:), tile3(:,:,:)
  logical :: is_root

  is_root = (mpp_pe() == mpp_root_pe())
  call mpi_pelist(pelist)

  ncid = -1
  if (is_root) call ncc(nf90_open(self%filename, NF90_NOWRITE, ncid), &
      'nf90_open '//trim(self%filename))

  nx_c   = self%iec - self%isc + 1
  ny_c   = self%jec - self%jsc + 1
  i_off  = self%isc - self%isd + 1
  j_off  = self%jsc - self%jsd + 1
  ! Map PE compute indices (isg-based) to 1-based file/global-buffer indices.
  ! FMS 2025.02 removed the ishift/jshift optional args from mpp_scatter, so
  ! we pre-apply the shift in the indices we pass.
  is_f   = self%isc - self%isg + 1
  ie_f   = self%iec - self%isg + 1
  js_f   = self%jsc - self%jsg + 1
  je_f   = self%jec - self%jsg + 1

  do v = 1, self%nvars
    select case (self%vars(v)%ndims)
    case (1)
      ! Global on every PE: PE 0 reads, broadcasts.
      n = size(self%vars(v)%data1d)
      if (is_root) call read_var_strided(ncid, self%vars(v)%name, &
          1, 1, n, 1, dst1=self%vars(v)%data1d)
      call mpp_broadcast(self%vars(v)%data1d, n, mpp_root_pe())

    case (2)
      allocate(tile2(nx_c, ny_c))
      if (is_root) then
        allocate(gbuf2(self%nx_g, self%ny_g))
        call read_var_strided(ncid, self%vars(v)%name, &
            1, 1, self%nx_g, self%ny_g, dst2=gbuf2)
      else
        allocate(gbuf2(1, 1))  ! dummy: mpp_scatter only reads input_data on root
      end if
      call mpp_scatter(is_f, ie_f, js_f, je_f, &
                       pelist, tile2, gbuf2, is_root)
      deallocate(gbuf2)
      self%vars(v)%data2d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1) = tile2
      deallocate(tile2)

    case (3)
      n3 = size(self%vars(v)%data3d, 3)
      allocate(tile3(nx_c, ny_c, n3))
      if (is_root) then
        allocate(gbuf3(self%nx_g, self%ny_g, n3))
        call read_var_strided(ncid, self%vars(v)%name, &
            1, 1, self%nx_g, self%ny_g, dst3=gbuf3)
      else
        allocate(gbuf3(1, 1, 1))
      end if
      call mpp_scatter(is_f, ie_f, js_f, je_f, n3, &
                       pelist, tile3, gbuf3, is_root)
      deallocate(gbuf3)
      self%vars(v)%data3d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1, :) = tile3
      deallocate(tile3)

    case (4)
      ! mpp_scatter is 2D/3D only; loop the outer (4th) dim and call 3D scatter.
      n3 = size(self%vars(v)%data4d, 3)
      n4 = size(self%vars(v)%data4d, 4)
      allocate(tile3(nx_c, ny_c, n3))
      if (is_root) then
        allocate(gbuf4(self%nx_g, self%ny_g, n3, n4))
        call read_var_strided(ncid, self%vars(v)%name, &
            1, 1, self%nx_g, self%ny_g, dst4=gbuf4)
      else
        allocate(gbuf3(1, 1, 1))  ! dummy for non-root in the 3D scatter call
      end if
      do k4 = 1, n4
        if (is_root) then
          call mpp_scatter(is_f, ie_f, js_f, je_f, n3, &
                           pelist, tile3, gbuf4(:,:,:,k4), is_root)
        else
          call mpp_scatter(is_f, ie_f, js_f, je_f, n3, &
                           pelist, tile3, gbuf3, is_root)
        end if
        self%vars(v)%data4d(i_off : i_off + nx_c - 1, &
                            j_off : j_off + ny_c - 1, :, k4) = tile3
      end do
      if (is_root) deallocate(gbuf4)
      if (allocated(gbuf3)) deallocate(gbuf3)
      deallocate(tile3)
    end select
  end do

  if (is_root) call ncc(nf90_close(ncid), 'nf90_close '//trim(self%filename))

  if (allocated(pelist)) deallocate(pelist)
end subroutine commit_reader_scatter




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
