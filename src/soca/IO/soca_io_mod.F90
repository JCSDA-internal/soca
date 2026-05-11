! (C) Copyright 2026 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! soca_io_mod
!
! Direct netcdf I/O for SOCA fields, replacing the FMS register_restart_field /
! save_restart / restore_state pattern. Designed so that LETKF parallel-ensemble
! I/O (issue #1125) becomes possible: FMS forces every read/write to be
! collective on the model's compute communicator, which prevents
! "PE i alone reads/writes member i" patterns. With direct nf90_* calls,
! a single PE can do its own I/O without coordinating with other PEs.
!
! API: register-then-commit lifecycle. The SOCA-side caller pattern is:
!
!     type(soca_io_writer) :: w
!     call w%init(domain, "myfile.nc")
!     call w%enqueue("lonh", self%lonh)
!     call w%enqueue("lath", self%lath)
!     ...
!     call w%commit()   ! gathers from all PEs and writes the file
!
! enqueue copies data in so the caller can mutate or free its source
! arrays before commit. Failure mode is mpp_error(FATAL).
!
! Read path: two implementations selectable via SOCA_IO_READ_MODE env var
! (default broadcast; see read_mode docstring below for rationale and
! benchmarks):
!   broadcast - PE 0 opens, nf90_get_var the full global field, then
!               mpp_broadcast to all PEs which slice their compute tile.
!   strided   - every PE opens NF90_NOWRITE independently and reads only
!               its compute-domain tile via nf90_get_var(start, count).
! Both paths share the same module-level handle / var-metadata cache so
! every (PE, filename) pair pays the nf90_open and nf90_inq cost once.
!
! Write path: PE 0 creates the file structure, then per-variable mpp_gather
! pulls each tile to PE 0, which writes the global field via nf90_put_var.
! Output is file-compatible with the legacy FMS path (same variable names,
! same dtypes, Time leading axis) but uses cleaner dimension names
! (xh/yh/Time) in place of FMS's auto-numbered xaxis_N/yaxis_N. The
! "checksum" attribute FMS attached to each variable is omitted.

module soca_io_mod

use netcdf
use mpi
use kinds, only: kind_real
use mpp_mod, only: mpp_gather, mpp_broadcast, mpp_pe, mpp_root_pe, mpp_npes, &
                   mpp_get_current_pelist, mpp_error, FATAL
use mpp_domains_mod, only: domain2D, &
                           mpp_get_compute_domain, mpp_get_global_domain, &
                           mpp_get_data_domain

implicit none
private

public :: soca_io_writer
public :: soca_io_reader
public :: soca_io_file_exists, soca_io_var_exists
public :: soca_io_close_all

integer, parameter :: MAX_NAME = 64

! Module-level cache of nf90_open(NF90_NOWRITE) handles, keyed by filename.
! Mirrors FMS fms_io's get_file_unit + files_read(i)%var(j) tables: a file
! is opened once on first reader_commit and the ncid is reused across every
! subsequent reader_commit on the same path. Per-variable metadata (varid,
! file_ndims, middle-dim sizes) is also cached so each repeat read only
! does the nf90_get_var, skipping inq_varid + inquire_variable +
! inquire_dimension. Closes happen only on soca_io_close_all (called from
! soca_geom_end) or process exit.
!
! Read-only by design -- if the caller needs to read a file it has just
! written, it must call soca_io_close_all first to drop the stale handle.
integer, parameter :: MAX_FILE_NDIMS = 8
type :: cached_var
  character(len=MAX_NAME) :: name = ''
  integer :: varid = -1
  integer :: file_ndims = 0
  integer :: middle_dims(MAX_FILE_NDIMS) = 1   ! sizes for dims 3..ndims-1
end type cached_var
type :: cached_open
  character(len=:), allocatable :: filename
  integer :: ncid = -1
  type(cached_var), allocatable :: vars(:)
  integer :: nvars = 0
end type cached_open
integer, parameter :: MAX_CACHED_OPENS = 256
integer, parameter :: VAR_CACHE_INIT   = 16
type(cached_open), save :: read_cache(MAX_CACHED_OPENS)
integer, save :: n_read_cached = 0

! Read strategy: PE-0-only read + mpp_broadcast (default) vs per-PE strided
! reads (each PE opens + nf90_get_var its tile). Override with
! SOCA_IO_READ_MODE=strided|broadcast at process start. Latched on first
! reader_commit.
!
! Broadcast is the default because real DA cycling reads files that were
! just written (model output, prev-cycle analyses) -- the page cache is
! always cold. On cold cache, 8 concurrent strided readers thrash the
! kernel readahead and pay 8x the open-syscall cost, while broadcast's
! single PE-0 sequential read keeps readahead happy and moves data to
! peers over fast in-node MPI. Measured on rancor 1deg 20-mem LETKF
! (cold cache): broadcast 13.5 s vs strided 30.5 s. Warm cache reverses
! the result (strided 3.9 s vs broadcast 13.3 s) but cycling never sees
! warm cache. On a parallel filesystem (Lustre/GPFS) or multi-node setup
! where the broadcast would have to cross the interconnect, strided may
! win again -- toggle the env var to test.
integer, parameter :: READ_MODE_STRIDED   = 1
integer, parameter :: READ_MODE_BROADCAST = 2
integer, save :: read_mode = 0   ! 0 = not yet resolved

type :: var_entry
  character(len=MAX_NAME) :: name = ''
  integer :: ndims = 0          ! 1, 2 or 3
  integer :: nlevels = 0        ! 1 if ndims<=2, else nz
  ! Data is copied into these buffers at enqueue time so the caller is free
  ! to deallocate / reuse their source arrays before commit. For 1D vars the
  ! data is global on every PE (no gather). For 2D/3D vars the buffer holds
  ! the COMPUTE-DOMAIN slice (1-based) and is gathered to the root in commit.
  real(kind=kind_real), allocatable :: data1d(:)
  real(kind=kind_real), allocatable :: data2d(:,:)
  real(kind=kind_real), allocatable :: data3d(:,:,:)
  character(len=MAX_NAME) :: long_name = ''
  character(len=MAX_NAME) :: units = 'none'
  character(len=1)        :: cartesian_axis = ' '  ! ' ', 'X', 'Y', 'Z', 'T'
end type var_entry

! Tracks one unique axis (X, Y or Z) by (size, domain_key) and remembers
! the netcdf dim/coord-var ids assigned to it during commit. Mirrors
! FMS's unique_axes(): each direction has a counter; a var's dim reuses
! an existing axis if its (size, domain) tuple matches one already seen,
! else a new axis is created.
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


! Read-side var_entry: the reader holds POINTERS to the caller's buffers
! and fills them in commit. Caller must keep buffers alive and not mutate
! them between enqueue and commit.
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
! init: prepare a writer for a specific file. domain pointer is stored;
! caller must keep the domain alive until commit returns.
!==============================================================================
subroutine writer_init(self, domain, filename)
  class(soca_io_writer), intent(inout) :: self
  type(domain2D), target, intent(in)   :: domain
  character(len=*),       intent(in)   :: filename

  self%filename = filename
  self%domain => domain

  call mpp_get_compute_domain(self%domain, self%isc, self%iec, self%jsc, self%jec)
  call mpp_get_global_domain(self%domain,  self%isg, self%ieg, self%jsg, self%jeg)
  self%nx_g = self%ieg - self%isg + 1
  self%ny_g = self%jeg - self%jsg + 1

  if (allocated(self%vars)) deallocate(self%vars)
  allocate(self%vars(64))
  self%nvars = 0
end subroutine writer_init


!==============================================================================
! enqueue_1d / enqueue_2d / enqueue_3d: register a variable for writing.
! The data is held by reference; do not deallocate before calling commit.
! 1D data is assumed global-on-every-PE (no gather, PE 0 writes directly).
! 2D/3D data is assumed compute-domain-decomposed (mpp_gather to PE 0).
!==============================================================================
subroutine writer_enqueue_1d(self, name, src, long_name, units, cartesian_axis)
  class(soca_io_writer),       intent(inout) :: self
  character(len=*),            intent(in)    :: name
  real(kind=kind_real),        intent(in)    :: src(:)
  character(len=*), optional,  intent(in)    :: long_name, units
  character(len=1), optional,  intent(in)    :: cartesian_axis

  call grow_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name    = name
  self%vars(self%nvars)%ndims   = 1
  self%vars(self%nvars)%nlevels = 1
  allocate(self%vars(self%nvars)%data1d(size(src)))
  self%vars(self%nvars)%data1d = src
  self%vars(self%nvars)%long_name = name
  self%vars(self%nvars)%units     = 'none'
  if (present(long_name))      self%vars(self%nvars)%long_name      = long_name
  if (present(units))          self%vars(self%nvars)%units          = units
  if (present(cartesian_axis)) self%vars(self%nvars)%cartesian_axis = cartesian_axis
end subroutine writer_enqueue_1d

! 2D / 3D enqueue. Caller passes the whole halo-inclusive array as it was
! allocated (e.g. self%lon, shape (isd:ied, jsd:jed)). Inside the routine
! the assumed-shape array re-bases to (1:nx_data, 1:ny_data); we use the
! domain's data-domain bounds to compute the offset of the compute slice
! within that 1-based shape and copy it out into a 1-based buffer.
subroutine writer_enqueue_2d(self, name, src, long_name, units)
  class(soca_io_writer),       intent(inout) :: self
  character(len=*),            intent(in)    :: name
  real(kind=kind_real),        intent(in)    :: src(:,:)
  character(len=*), optional,  intent(in)    :: long_name, units

  integer :: nx, ny, i_off, j_off, isd, ied, jsd, jed

  call grow_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name    = name
  self%vars(self%nvars)%ndims   = 2
  self%vars(self%nvars)%nlevels = 1
  nx = self%iec - self%isc + 1
  ny = self%jec - self%jsc + 1
  allocate(self%vars(self%nvars)%data2d(nx, ny))
  call mpp_get_data_domain(self%domain, isd, ied, jsd, jed)
  i_off = self%isc - isd + 1
  j_off = self%jsc - jsd + 1
  self%vars(self%nvars)%data2d(:,:) = src(i_off : i_off + nx - 1, &
                                          j_off : j_off + ny - 1)
  self%vars(self%nvars)%long_name = name
  self%vars(self%nvars)%units     = 'none'
  if (present(long_name)) self%vars(self%nvars)%long_name = long_name
  if (present(units))     self%vars(self%nvars)%units     = units
end subroutine writer_enqueue_2d

subroutine writer_enqueue_3d(self, name, src, long_name, units)
  class(soca_io_writer),       intent(inout) :: self
  character(len=*),            intent(in)    :: name
  real(kind=kind_real),        intent(in)    :: src(:,:,:)
  character(len=*), optional,  intent(in)    :: long_name, units

  integer :: nx, ny, nz, i_off, j_off, isd, ied, jsd, jed

  call grow_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name    = name
  self%vars(self%nvars)%ndims   = 3
  self%vars(self%nvars)%nlevels = size(src, 3)
  nx = self%iec - self%isc + 1
  ny = self%jec - self%jsc + 1
  nz = size(src, 3)
  allocate(self%vars(self%nvars)%data3d(nx, ny, nz))
  call mpp_get_data_domain(self%domain, isd, ied, jsd, jed)
  i_off = self%isc - isd + 1
  j_off = self%jsc - jsd + 1
  self%vars(self%nvars)%data3d(:,:,:) = src(i_off : i_off + nx - 1, &
                                            j_off : j_off + ny - 1, :)
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
! commit: PE 0 creates the file structure, then we gather each variable in
! turn via mpp_gather and PE 0 writes via nf90_put_var. This is the
! FMS-equivalent path (mpp_io threading=MPP_SINGLE); speedup vs FMS is
! not the goal -- a clean, debuggable, FMS-free baseline is.
!==============================================================================
subroutine writer_commit(self)
  class(soca_io_writer), intent(inout) :: self

  integer :: ncid, dimid_t, varid_t, v, k
  integer, allocatable :: varids(:)
  integer, allocatable :: var_x_idx(:), var_y_idx(:), var_z_idx(:)
  type(axis_entry), allocatable :: x_axes(:), y_axes(:), z_axes(:)
  integer :: nx_axes, ny_axes, nz_axes
  logical :: is_root
  real(kind=kind_real), allocatable :: gbuf2d(:,:), gbuf3d(:,:,:)
  integer, allocatable :: pelist(:)
  integer :: nprocs, dom_key
  integer, parameter :: MAX_AXES_PER_DIR = 32

  is_root = (mpp_pe() == mpp_root_pe())
  call mpi_pelist(pelist, nprocs)

  ! ----------------------------------------------------------------------
  ! Build per-direction unique-axis tables matching FMS's algorithm:
  ! a var's dim at position k is matched against existing axes by
  ! (size, domain_key); same tuple -> reuse; new tuple -> append.
  ! ----------------------------------------------------------------------
  allocate(x_axes(MAX_AXES_PER_DIR), y_axes(MAX_AXES_PER_DIR), z_axes(MAX_AXES_PER_DIR))
  allocate(var_x_idx(self%nvars), var_y_idx(self%nvars), var_z_idx(self%nvars))
  var_x_idx = 0; var_y_idx = 0; var_z_idx = 0
  nx_axes = 0; ny_axes = 0; nz_axes = 0

  do v = 1, self%nvars
    ! 1D vars carry no domain (they are global-on-every-PE); 2D/3D vars
    ! all share self%domain in this writer (key = 1).
    dom_key = 1
    if (self%vars(v)%ndims == 1) dom_key = 0

    ! Every var contributes an X axis (Fortran first dim). 2D/3D vars hold
    ! only their compute slice locally, so we use the writer's GLOBAL extents
    ! for the axis sizes; 1D vars are already global on every PE so size()
    ! is fine.
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

  ! ----------------------------------------------------------------------
  ! Phase 1: PE 0 defines the file structure -- dims, coord vars, data vars.
  ! ----------------------------------------------------------------------
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

  ! ----------------------------------------------------------------------
  ! Phase 2: gather and write each user variable.
  ! ----------------------------------------------------------------------
  if (is_root) allocate(gbuf2d(self%nx_g, self%ny_g))

  do v = 1, self%nvars
    if (self%vars(v)%ndims == 1) then
      if (is_root) then
        call ncc(nf90_put_var(ncid, varids(v), self%vars(v)%data1d, &
            start=[1, 1], count=[size(self%vars(v)%data1d), 1]), &
            'put '//trim(self%vars(v)%name))
      end if
    else if (self%vars(v)%ndims == 2) then
      call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
                      self%vars(v)%data2d, gbuf2d, is_root)
      if (is_root) then
        call ncc(nf90_put_var(ncid, varids(v), gbuf2d, &
            start=[1, 1, 1], count=[self%nx_g, self%ny_g, 1]), &
            'put '//trim(self%vars(v)%name))
      end if
    else
      if (is_root) then
        if (allocated(gbuf3d)) deallocate(gbuf3d)
        allocate(gbuf3d(self%nx_g, self%ny_g, self%vars(v)%nlevels))
      end if
      do k = 1, self%vars(v)%nlevels
        call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
                        self%vars(v)%data3d(:,:,k), gbuf2d, is_root)
        if (is_root) gbuf3d(:,:,k) = gbuf2d
      end do
      if (is_root) then
        call ncc(nf90_put_var(ncid, varids(v), gbuf3d, &
            start=[1, 1, 1, 1], count=[self%nx_g, self%ny_g, self%vars(v)%nlevels, 1]), &
            'put '//trim(self%vars(v)%name))
      end if
    end if
  end do

  ! ----------------------------------------------------------------------
  ! Phase 3: close.
  ! ----------------------------------------------------------------------
  if (is_root) then
    call ncc(nf90_close(ncid), 'nf90_close')
    deallocate(varids, gbuf2d)
    if (allocated(gbuf3d)) deallocate(gbuf3d)
  end if
  deallocate(pelist, x_axes, y_axes, z_axes, var_x_idx, var_y_idx, var_z_idx)

  ! release references; caller's data unaffected
  if (allocated(self%vars)) deallocate(self%vars)
  self%nvars = 0
end subroutine writer_commit


!==============================================================================
! Reader: init / enqueue_*/ commit. The caller's buffer stays in place;
! enqueue records a pointer, commit fills the buffer in the compute-domain
! interior. Halos are left untouched -- the caller refreshes them via
! mpp_update_domains the same way it did under FMS.
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


! 1D vars are global on every PE (no scatter): PE 0 reads, broadcasts,
! every PE has the full array in its own buffer.
subroutine reader_enqueue_1d(self, name, dst)
  class(soca_io_reader),                  intent(inout) :: self
  character(len=*),                       intent(in)    :: name
  real(kind=kind_real), target,           intent(inout) :: dst(:)

  call grow_reader_if_needed(self)
  self%nvars = self%nvars + 1
  self%vars(self%nvars)%name  = name
  self%vars(self%nvars)%ndims = 1
  self%vars(self%nvars)%data1d => dst
end subroutine reader_enqueue_1d


! 2D vars are domain-decomposed: PE 0 reads the global field, broadcasts,
! every PE copies its compute slice into the caller's data-domain buffer.
! Caller buffer must be sized exactly (data_xsize, data_ysize); halos
! remain undefined and the caller refreshes them with mpp_update_domains.
subroutine reader_enqueue_2d(self, name, dst)
  class(soca_io_reader),                  intent(inout) :: self
  character(len=*),                       intent(in)    :: name
  real(kind=kind_real), target,           intent(inout) :: dst(:,:)

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
! Read all enqueued vars. Every PE opens the file in NF90_NOWRITE and pulls
! only its compute-domain tile via nf90_get_var(start, count) -- this mirrors
! FMS's MPP_READ_2DDECOMP path: no PE-0 bottleneck, no mpp_broadcast, N
! parallel reads of (nx_c x ny_c x nz) slabs. Classic / 64-bit-offset netcdf
! files permit concurrent read-only opens; the library state is process-local.
! For 1D vars (file-global and identical on every PE) we still let every PE
! read the whole array independently rather than broadcast from root.
!==============================================================================
subroutine reader_commit(self)
  class(soca_io_reader), intent(inout) :: self

  if (read_mode == 0) call resolve_read_mode()
  select case (read_mode)
  case (READ_MODE_STRIDED)
    call commit_reader_strided(self)
  case (READ_MODE_BROADCAST)
    call commit_reader_broadcast(self)
  case default
    call mpp_error(FATAL, 'soca_io_mod reader_commit: unknown read_mode')
  end select

  ! release pointers; caller's buffers untouched (they hold the read data)
  if (allocated(self%vars)) deallocate(self%vars)
  self%nvars = 0
end subroutine reader_commit


!==============================================================================
! Per-PE strided read implementation. Every PE opens the file (cached) and
! pulls only its compute-domain tile via nf90_get_var(start, count).
!==============================================================================
subroutine commit_reader_strided(self)
  class(soca_io_reader), intent(inout) :: self

  integer :: ncid, file_idx, v, n3, n4
  integer :: nx_c, ny_c, i_off, j_off, i_start, j_start
  real(kind=kind_real), allocatable :: tile2(:,:), tile3(:,:,:), tile4(:,:,:,:)

  call get_read_ncid(self%filename, ncid, file_idx)

  nx_c    = self%iec - self%isc + 1
  ny_c    = self%jec - self%jsc + 1
  i_off   = self%isc - self%isd + 1     ! 1-based offset into data-domain buf
  j_off   = self%jsc - self%jsd + 1
  i_start = self%isc - self%isg + 1     ! 1-based start in the global file dim
  j_start = self%jsc - self%jsg + 1

  do v = 1, self%nvars
    select case (self%vars(v)%ndims)
    case (1)
      call read_var_strided(ncid, file_idx, self%vars(v)%name, &
          1, 1, size(self%vars(v)%data1d), 1, dst1=self%vars(v)%data1d)

    case (2)
      allocate(tile2(nx_c, ny_c))
      call read_var_strided(ncid, file_idx, self%vars(v)%name, &
          i_start, j_start, nx_c, ny_c, dst2=tile2)
      self%vars(v)%data2d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1) = tile2
      deallocate(tile2)

    case (3)
      n3 = size(self%vars(v)%data3d, 3)
      allocate(tile3(nx_c, ny_c, n3))
      tile3 = 0.0_kind_real
      call read_var_strided(ncid, file_idx, self%vars(v)%name, &
          i_start, j_start, nx_c, ny_c, dst3=tile3)
      self%vars(v)%data3d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1, :) = tile3
      deallocate(tile3)

    case (4)
      n3 = size(self%vars(v)%data4d, 3)
      n4 = size(self%vars(v)%data4d, 4)
      allocate(tile4(nx_c, ny_c, n3, n4))
      tile4 = 0.0_kind_real
      call read_var_strided(ncid, file_idx, self%vars(v)%name, &
          i_start, j_start, nx_c, ny_c, dst4=tile4)
      self%vars(v)%data4d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1, :, :) = tile4
      deallocate(tile4)
    end select
  end do
end subroutine commit_reader_strided


!==============================================================================
! PE-0 read + broadcast implementation, mirroring FMS mpp_io's MPP_SINGLE
! threading path: PE 0 nf90_get_var the entire global field, mpp_broadcast
! to all peers, each PE then slices its compute tile out of the global buf.
! Only PE 0 opens the file (and only PE 0's read_cache fills up).
!==============================================================================
subroutine commit_reader_broadcast(self)
  class(soca_io_reader), intent(inout) :: self

  integer :: ncid, file_idx, v, n, n3, n4
  integer :: nx_c, ny_c, i_off, j_off, ig0, jg0
  real(kind=kind_real), allocatable :: gbuf2(:,:), gbuf3(:,:,:), gbuf4(:,:,:,:)
  logical :: is_root

  is_root = (mpp_pe() == mpp_root_pe())

  ncid = -1
  file_idx = -1
  if (is_root) call get_read_ncid(self%filename, ncid, file_idx)

  nx_c = self%iec - self%isc + 1
  ny_c = self%jec - self%jsc + 1
  i_off = self%isc - self%isd + 1
  j_off = self%jsc - self%jsd + 1
  ig0   = self%isc - self%isg + 1
  jg0   = self%jsc - self%jsg + 1

  do v = 1, self%nvars
    select case (self%vars(v)%ndims)
    case (1)
      n = size(self%vars(v)%data1d)
      if (is_root) call read_var_strided(ncid, file_idx, self%vars(v)%name, &
          1, 1, n, 1, dst1=self%vars(v)%data1d)
      call mpp_broadcast(self%vars(v)%data1d, n, mpp_root_pe())

    case (2)
      allocate(gbuf2(self%nx_g, self%ny_g))
      if (is_root) call read_var_strided(ncid, file_idx, self%vars(v)%name, &
          1, 1, self%nx_g, self%ny_g, dst2=gbuf2)
      call mpp_broadcast(gbuf2, self%nx_g * self%ny_g, mpp_root_pe())
      self%vars(v)%data2d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1) &
        = gbuf2(ig0 : ig0 + nx_c - 1, jg0 : jg0 + ny_c - 1)
      deallocate(gbuf2)

    case (3)
      n3 = size(self%vars(v)%data3d, 3)
      allocate(gbuf3(self%nx_g, self%ny_g, n3))
      gbuf3 = 0.0_kind_real
      if (is_root) call read_var_strided(ncid, file_idx, self%vars(v)%name, &
          1, 1, self%nx_g, self%ny_g, dst3=gbuf3)
      call mpp_broadcast(gbuf3, self%nx_g * self%ny_g * n3, mpp_root_pe())
      self%vars(v)%data3d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1, :) &
        = gbuf3(ig0 : ig0 + nx_c - 1, jg0 : jg0 + ny_c - 1, :)
      deallocate(gbuf3)

    case (4)
      n3 = size(self%vars(v)%data4d, 3)
      n4 = size(self%vars(v)%data4d, 4)
      allocate(gbuf4(self%nx_g, self%ny_g, n3, n4))
      gbuf4 = 0.0_kind_real
      if (is_root) call read_var_strided(ncid, file_idx, self%vars(v)%name, &
          1, 1, self%nx_g, self%ny_g, dst4=gbuf4)
      call mpp_broadcast(gbuf4, self%nx_g * self%ny_g * n3 * n4, mpp_root_pe())
      self%vars(v)%data4d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1, :, :) &
        = gbuf4(ig0 : ig0 + nx_c - 1, jg0 : jg0 + ny_c - 1, :, :)
      deallocate(gbuf4)
    end select
  end do
end subroutine commit_reader_broadcast


!==============================================================================
! Latch the read mode from the SOCA_IO_READ_MODE env var on first commit.
! Default = strided; broadcast available for A/B benchmarking. Unrecognized
! values abort so typos surface immediately.
!==============================================================================
subroutine resolve_read_mode()
  character(len=32) :: mode_str
  integer :: status

  call get_environment_variable("SOCA_IO_READ_MODE", mode_str, status=status)
  if (status /= 0) then
    read_mode = READ_MODE_BROADCAST
    return
  end if
  select case (trim(mode_str))
  case ("broadcast", "")
    read_mode = READ_MODE_BROADCAST
  case ("strided")
    read_mode = READ_MODE_STRIDED
  case default
    call mpp_error(FATAL, 'soca_io_mod: SOCA_IO_READ_MODE="'//trim(mode_str)//&
        '" unrecognized (expected broadcast or strided)')
  end select
end subroutine resolve_read_mode


!==============================================================================
! Look up (or open and cache) the read handle for filename. Returns both the
! ncid and the file_idx into read_cache(); read_var_strided uses file_idx to
! look up per-var cached metadata. Per-PE cache.
!==============================================================================
subroutine get_read_ncid(filename, ncid, file_idx)
  character(len=*), intent(in)  :: filename
  integer,          intent(out) :: ncid, file_idx
  integer :: i

  do i = 1, n_read_cached
    if (allocated(read_cache(i)%filename)) then
      if (read_cache(i)%filename == filename) then
        ncid     = read_cache(i)%ncid
        file_idx = i
        return
      end if
    end if
  end do

  if (n_read_cached >= MAX_CACHED_OPENS) call mpp_error(FATAL, &
      'soca_io_mod get_read_ncid: cache full (increase MAX_CACHED_OPENS)')
  call ncc(nf90_open(filename, NF90_NOWRITE, ncid), 'nf90_open '//trim(filename))
  n_read_cached = n_read_cached + 1
  read_cache(n_read_cached)%filename = filename
  read_cache(n_read_cached)%ncid     = ncid
  allocate(read_cache(n_read_cached)%vars(VAR_CACHE_INIT))
  read_cache(n_read_cached)%nvars    = 0
  file_idx = n_read_cached
end subroutine get_read_ncid


!==============================================================================
! Fetch cached metadata for (file_idx, varname), or do the inq calls and
! cache them on first use. middle_dims(1:max(0, ndims-3)) gives the file's
! dim sizes for axes 3..ndims-1 (the layer / category / "extra" dims).
!==============================================================================
subroutine get_var_meta(file_idx, name, varid, file_ndims, middle_dims)
  integer,          intent(in)  :: file_idx
  character(len=*), intent(in)  :: name
  integer,          intent(out) :: varid, file_ndims
  integer,          intent(out) :: middle_dims(MAX_FILE_NDIMS)

  integer :: i, dd, sz, ncid, nv
  integer :: file_dimids(MAX_FILE_NDIMS)
  type(cached_var), allocatable :: tmp(:)

  nv = read_cache(file_idx)%nvars
  do i = 1, nv
    if (read_cache(file_idx)%vars(i)%name == name) then
      varid       = read_cache(file_idx)%vars(i)%varid
      file_ndims  = read_cache(file_idx)%vars(i)%file_ndims
      middle_dims = read_cache(file_idx)%vars(i)%middle_dims
      return
    end if
  end do

  ! First-touch -- do the inq calls and stash the result.
  ncid = read_cache(file_idx)%ncid
  call ncc(nf90_inq_varid(ncid, trim(name), varid), 'inq '//trim(name))
  call ncc(nf90_inquire_variable(ncid, varid, ndims=file_ndims, dimids=file_dimids), &
           'inquire '//trim(name))
  if (file_ndims > MAX_FILE_NDIMS) call mpp_error(FATAL, &
      'soca_io_mod get_var_meta: '//trim(name)//' exceeds MAX_FILE_NDIMS')
  middle_dims = 1
  do dd = 3, file_ndims - 1
    call ncc(nf90_inquire_dimension(ncid, file_dimids(dd), len=sz), 'dim '//trim(name))
    middle_dims(dd) = sz
  end do

  if (nv >= size(read_cache(file_idx)%vars)) then
    allocate(tmp(2 * size(read_cache(file_idx)%vars)))
    tmp(1:nv) = read_cache(file_idx)%vars(1:nv)
    call move_alloc(tmp, read_cache(file_idx)%vars)
  end if
  nv = nv + 1
  read_cache(file_idx)%nvars               = nv
  read_cache(file_idx)%vars(nv)%name        = name
  read_cache(file_idx)%vars(nv)%varid       = varid
  read_cache(file_idx)%vars(nv)%file_ndims  = file_ndims
  read_cache(file_idx)%vars(nv)%middle_dims = middle_dims
end subroutine get_var_meta


!==============================================================================
! Close every cached read handle and drop the table. Call from geometry
! shutdown (soca_geom%end) so the next geometry can open files fresh, and
! so a writer that re-creates a previously-read file won't see a stale
! cached ncid pointing at a deleted inode.
!==============================================================================
subroutine soca_io_close_all()
  integer :: i, status
  do i = 1, n_read_cached
    if (read_cache(i)%ncid >= 0) then
      status = nf90_close(read_cache(i)%ncid)
      read_cache(i)%ncid = -1
    end if
    if (allocated(read_cache(i)%filename)) deallocate(read_cache(i)%filename)
    if (allocated(read_cache(i)%vars))     deallocate(read_cache(i)%vars)
    read_cache(i)%nvars = 0
  end do
  n_read_cached = 0
end subroutine soca_io_close_all


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
    allocate(idxbuf(axes(j)%size))
    do i = 1, axes(j)%size
      idxbuf(i) = real(i, kind=kind_real)
    end do
    call ncc(nf90_put_var(ncid, axes(j)%varid, idxbuf), 'put coord var')
    deallocate(idxbuf)
  end do
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
! Strided read of a single variable into a caller-owned contiguous tile.
! Builds start/count from (i_start, j_start, nx, ny) and the variable's file
! shape, tolerating "extra" degenerate file dims (size 1) that the destination
! buffer doesn't enumerate. Same middle-dim trick FMS register_restart_field+
! restore_state did transparently. Concrete cases:
!   - 3D state Salt(time, Layer, lath, lonh) -> 3D tile (nx, ny, Layer)
!     file_ndims=4, count=[nx, ny, Layer, 1].
!   - 5D CICE Tsnz_h(time, nc=5, nksnow=1, nj, ni) -> 3D tile (nx, ny, 5)
!     file_ndims=5, count=[nx, ny, 1, 5, 1] = [nx, ny, nksnow=1, nc=5, time=1].
! In all 2D+ cases count(1)=nx, count(2)=ny, count(ndims)=1 and the middle
! count entries come from the file's actual dim sizes; total element count is
! preserved so netcdf-fortran reads into the smaller-rank buffer correctly.
! For 1D file vars (handled via dst1) count(1)=nx and we ignore j_start/ny.
!==============================================================================
subroutine read_var_strided(ncid, file_idx, name, i_start, j_start, nx, ny, &
                            dst1, dst2, dst3, dst4)
  integer,                       intent(in)  :: ncid, file_idx
  integer,                       intent(in)  :: i_start, j_start, nx, ny
  character(len=*),              intent(in)  :: name
  real(kind=kind_real), optional, intent(out) :: dst1(:)
  real(kind=kind_real), optional, intent(out) :: dst2(:,:)
  real(kind=kind_real), optional, intent(out) :: dst3(:,:,:)
  real(kind=kind_real), optional, intent(out) :: dst4(:,:,:,:)

  integer :: varid, file_ndims, dd
  integer :: middle_dims(MAX_FILE_NDIMS)
  integer, allocatable :: ct(:), st(:)

  call get_var_meta(file_idx, name, varid, file_ndims, middle_dims)
  allocate(st(file_ndims), ct(file_ndims))
  st = 1
  ct = 1
  if (present(dst1)) then
    ! 1D var: count(1)=nx, rest stays 1
    ct(1) = nx
    if (file_ndims >= 2) ct(file_ndims) = 1
  else
    if (file_ndims < 2) call mpp_error(FATAL, &
        'soca_io_mod read_var_strided: '//trim(name)//' has fewer than 2 file dims')
    st(1) = i_start
    st(2) = j_start
    ct(1) = nx
    ct(2) = ny
    ! File may have an optional trailing time dim and any number of middle dims.
    ! Trailing-time case: ct(file_ndims)=1 + middle dims from file (e.g. Salt's
    ! Layer dim, Tsnz_h's nc/nksnow). Pure 2D file (no time, no middle): just
    ! the leading x/y stride is enough.
    if (file_ndims >= 3) then
      ct(file_ndims) = 1
      do dd = 3, file_ndims - 1
        ct(dd) = middle_dims(dd)
      end do
    end if
  end if

  if (present(dst2)) then
    call ncc(nf90_get_var(ncid, varid, dst2, start=st, count=ct), 'get '//trim(name))
  else if (present(dst3)) then
    call ncc(nf90_get_var(ncid, varid, dst3, start=st, count=ct), 'get '//trim(name))
  else if (present(dst4)) then
    call ncc(nf90_get_var(ncid, varid, dst4, start=st, count=ct), 'get '//trim(name))
  else if (present(dst1)) then
    call ncc(nf90_get_var(ncid, varid, dst1, start=st, count=ct), 'get '//trim(name))
  end if
  deallocate(ct, st)
end subroutine read_var_strided


!==============================================================================
! Helper: fetch mpp's current pelist (= geometry's f_comm pelist) for use
! with mpp_gather. Using MPI_COMM_WORLD's pelist is wrong in ensemble mode
! where each MPI task has its own size-1 mpp world.
!==============================================================================
subroutine mpi_pelist(pelist, nprocs)
  integer, allocatable, intent(out) :: pelist(:)
  integer,              intent(out) :: nprocs
  nprocs = mpp_npes()
  allocate(pelist(nprocs))
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

end module soca_io_mod
