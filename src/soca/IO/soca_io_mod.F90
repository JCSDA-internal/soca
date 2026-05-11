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
!     call w%init(domain, "myfile.nc", method=soca_io_method_from_config(f_conf))
!     call w%enqueue("lonh", self%lonh)
!     call w%enqueue("lath", self%lath)
!     ...
!     call w%commit()   ! gathers from all PEs and writes the file
!
! The caller code is identical for SOCA_IO_SERIAL and SOCA_IO_PARALLEL.
! All gather/write strategy lives inside writer_commit's select case;
! enqueue copies data in so the caller can mutate or free its source
! arrays before commit, and commit's failure mode (mpp_error(FATAL))
! is the same for both methods. Don't add method-conditional logic to
! init or enqueue -- if a strategy needs more state, store it in the
! writer at init time and read it inside commit_<method>.
!
! SOCA_IO_SERIAL: all gathers and nf90_put_var calls happen on PE 0.
! Output is file-compatible with the legacy FMS path (same variable
! names, same dtypes, Time leading axis) but uses cleaner dimension
! names (xh/yh/Time) in place of FMS's auto-numbered xaxis_N/yaxis_N.
! The "checksum" attribute FMS attached to each variable is omitted.
!
! SOCA_IO_PARALLEL (TODO) replaces the per-variable mpp_gather loop with
! MPI_Igather/Igatherv posted across all variables before any wait, with
! per-var destination buffers. Optionally distributes writer-root duty
! round-robin across PEs (each writer-PE serializes its assigned vars'
! nf90_put_var via token-passing, since HDF5 forbids concurrent
! same-file writes -- verified in ~/work/jedi/io_spike/).

module soca_io_mod

use netcdf
use mpi
use kinds, only: kind_real
use mpp_mod, only: mpp_gather, mpp_broadcast, mpp_pe, mpp_root_pe, mpp_npes, &
                   mpp_get_current_pelist, mpp_error, FATAL
use mpp_domains_mod, only: domain2D, &
                           mpp_get_compute_domain, mpp_get_global_domain, &
                           mpp_get_data_domain
use fckit_configuration_module, only: fckit_configuration

implicit none
private

public :: soca_io_writer
public :: soca_io_reader
public :: soca_io_method_from_config
public :: soca_io_file_exists, soca_io_var_exists

! YAML schema (under the geometry block):
!     io:
!       method: serial | parallel
! `serial`   = PE 0 gathers and writes all variables (default if io: absent)
! `parallel` = parallel non-blocking gathers and a token-rotated serial
!              nf90_put_var on multiple writer PEs (NOT YET IMPLEMENTED)
integer, parameter, public :: SOCA_IO_SERIAL   = 1
integer, parameter, public :: SOCA_IO_PARALLEL = 2

integer, parameter :: MAX_NAME = 64

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
  integer :: method = SOCA_IO_SERIAL
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
  integer :: method = SOCA_IO_SERIAL
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
! Read the I/O method from a YAML config block:
!     io:
!       method: serial | parallel
! Defaults to SOCA_IO_SERIAL if the io: block is absent. An unrecognized
! method aborts with mpp_error(FATAL) so a typo surfaces immediately rather
! than silently using the wrong path.
!==============================================================================
function soca_io_method_from_config(f_conf) result(method)
  type(fckit_configuration), intent(in) :: f_conf
  integer :: method
  type(fckit_configuration) :: io_conf
  character(len=:), allocatable :: str
  logical :: ok

  method = SOCA_IO_SERIAL
  if (.not. f_conf%has("io")) return
  ok = f_conf%get("io", io_conf)
  if (.not. ok) return
  if (.not. io_conf%has("method")) return
  ok = io_conf%get("method", str)
  if (.not. ok) return
  select case (trim(str))
  case ("serial")
    method = SOCA_IO_SERIAL
  case ("parallel")
    method = SOCA_IO_PARALLEL
  case default
    call mpp_error(FATAL, &
        'soca_io_mod: unknown io.method "'//trim(str)//&
        '" (expected serial or parallel)')
  end select
end function soca_io_method_from_config


!==============================================================================
! init: prepare a writer for a specific file. domain pointer is stored;
! caller must keep the domain alive until commit returns.
!==============================================================================
subroutine writer_init(self, domain, filename, method)
  class(soca_io_writer), intent(inout) :: self
  type(domain2D), target, intent(in)   :: domain
  character(len=*),       intent(in)   :: filename
  integer, optional,      intent(in)   :: method

  self%filename = filename
  self%method = SOCA_IO_SERIAL
  if (present(method)) self%method = method
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
! commit: dispatches to the actual write implementation.
!==============================================================================
subroutine writer_commit(self)
  class(soca_io_writer), intent(inout) :: self

  select case (self%method)
  case (SOCA_IO_SERIAL)
    call commit_serial(self)
  case (SOCA_IO_PARALLEL)
    call mpp_error(FATAL, &
        'soca_io_mod: io.method=parallel not yet implemented')
  case default
    call mpp_error(FATAL, 'soca_io_mod: unknown io method')
  end select

  ! release references; caller's data unaffected
  if (allocated(self%vars)) deallocate(self%vars)
  self%nvars = 0
end subroutine writer_commit


!==============================================================================
! Serial implementation: PE 0 creates the file structure, then we gather each
! variable in turn via mpp_gather and PE 0 writes via nf90_put_var. This is
! the FMS-equivalent path; speedup vs FMS is not the goal -- a clean,
! debuggable, FMS-free baseline is.
!==============================================================================
subroutine commit_serial(self)
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
end subroutine commit_serial


!==============================================================================
! Reader: init / enqueue_*/ commit. The caller's buffer stays in place;
! enqueue records a pointer, commit fills the buffer in the compute-domain
! interior. Halos are left untouched -- the caller refreshes them via
! mpp_update_domains the same way it did under FMS.
!==============================================================================
subroutine reader_init(self, domain, filename, method)
  class(soca_io_reader), intent(inout) :: self
  type(domain2D), target, intent(in)   :: domain
  character(len=*),       intent(in)   :: filename
  integer, optional,      intent(in)   :: method

  self%filename = filename
  self%method = SOCA_IO_SERIAL
  if (present(method)) self%method = method
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


subroutine reader_commit(self)
  class(soca_io_reader), intent(inout) :: self

  select case (self%method)
  case (SOCA_IO_SERIAL)
    call commit_serial_read(self)
  case (SOCA_IO_PARALLEL)
    call mpp_error(FATAL, &
        'soca_io_mod: io.method=parallel not yet implemented')
  case default
    call mpp_error(FATAL, 'soca_io_mod: unknown io method')
  end select

  ! release pointers; caller's buffers untouched (they hold the read data)
  if (allocated(self%vars)) deallocate(self%vars)
  self%nvars = 0
end subroutine reader_commit


!==============================================================================
! Serial read: PE 0 nf90_open + nf90_get_var the global field for each var,
! MPI_Bcast to all PEs, each PE copies its compute slice into the caller's
! data-domain buffer. Matches the FMS register_restart_field+restore_state
! pattern bandwidth-wise (FMS does the same: read global on PE 0, broadcast,
! each PE indexes its slice).
!==============================================================================
subroutine commit_serial_read(self)
  class(soca_io_reader), intent(inout) :: self

  integer :: ncid, v, n, n3, n4
  integer :: varid
  integer :: i_off, j_off, nx_c, ny_c, ig0, jg0
  real(kind=kind_real), allocatable :: gbuf2d(:,:)
  real(kind=kind_real), allocatable :: gbuf3d(:,:,:)
  real(kind=kind_real), allocatable :: gbuf4d(:,:,:,:)
  logical :: is_root

  is_root = (mpp_pe() == mpp_root_pe())

  if (is_root) then
    call ncc(nf90_open(self%filename, NF90_NOWRITE, ncid), &
        'nf90_open '//trim(self%filename))
  end if

  ! gbuf2d holds the global 2D field during the read+broadcast for 2D vars;
  ! every PE allocates because it's the destination of mpp_broadcast.
  allocate(gbuf2d(self%nx_g, self%ny_g))

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
      if (is_root) then
        call ncc(nf90_inq_varid(ncid, trim(self%vars(v)%name), varid), &
            'inq '//trim(self%vars(v)%name))
        call ncc(nf90_get_var(ncid, varid, self%vars(v)%data1d, &
            start=[1, 1], count=[n, 1]), &
            'get '//trim(self%vars(v)%name))
      end if
      call mpp_broadcast(self%vars(v)%data1d, n, mpp_root_pe())

    case (2)
      if (is_root) then
        call ncc(nf90_inq_varid(ncid, trim(self%vars(v)%name), varid), &
            'inq '//trim(self%vars(v)%name))
        call ncc(nf90_get_var(ncid, varid, gbuf2d, &
            start=[1, 1, 1], count=[self%nx_g, self%ny_g, 1]), &
            'get '//trim(self%vars(v)%name))
      end if
      call mpp_broadcast(gbuf2d, self%nx_g * self%ny_g, mpp_root_pe())
      self%vars(v)%data2d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1) &
        = gbuf2d(ig0 : ig0 + nx_c - 1, jg0 : jg0 + ny_c - 1)

    case (3)
      n3 = size(self%vars(v)%data3d, 3)
      if (allocated(gbuf3d)) deallocate(gbuf3d)
      allocate(gbuf3d(self%nx_g, self%ny_g, n3))
      gbuf3d = 0.0_kind_real
      if (is_root) call read_global_var(ncid, self%vars(v)%name, &
          self%nx_g, self%ny_g, gbuf3d=gbuf3d)
      call mpp_broadcast(gbuf3d, self%nx_g * self%ny_g * n3, mpp_root_pe())
      self%vars(v)%data3d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1, :) &
        = gbuf3d(ig0 : ig0 + nx_c - 1, jg0 : jg0 + ny_c - 1, :)

    case (4)
      n3 = size(self%vars(v)%data4d, 3)
      n4 = size(self%vars(v)%data4d, 4)
      if (allocated(gbuf4d)) deallocate(gbuf4d)
      allocate(gbuf4d(self%nx_g, self%ny_g, n3, n4))
      gbuf4d = 0.0_kind_real
      if (is_root) call read_global_var(ncid, self%vars(v)%name, &
          self%nx_g, self%ny_g, gbuf4d=gbuf4d)
      call mpp_broadcast(gbuf4d, self%nx_g * self%ny_g * n3 * n4, mpp_root_pe())
      self%vars(v)%data4d(i_off : i_off + nx_c - 1, &
                          j_off : j_off + ny_c - 1, :, :) &
        = gbuf4d(ig0 : ig0 + nx_c - 1, jg0 : jg0 + ny_c - 1, :, :)
    end select
  end do

  if (is_root) then
    call ncc(nf90_close(ncid), 'nf90_close '//trim(self%filename))
  end if
  deallocate(gbuf2d)
  if (allocated(gbuf3d)) deallocate(gbuf3d)
  if (allocated(gbuf4d)) deallocate(gbuf4d)
end subroutine commit_serial_read


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
! Read a global var on PE 0 into either a 3D or 4D destination buffer,
! tolerating "extra" degenerate file dims (size 1) that the buffer doesn't
! enumerate. This is what FMS register_restart_field+restore_state did
! transparently. Concrete cases:
!   - 3D state field Salt(time, Layer, lath, lonh) -> 3D buffer (lonh, lath, Layer)
!     file_ndims=4, count=[nx_g, ny_g, Layer, 1].
!   - 5D CICE field Tsnz_h(time, nc=5, nksnow=1, nj, ni) -> 3D buffer (ni, nj, nc=5)
!     file_ndims=5, count=[nx_g, ny_g, nksnow=1, nc=5, time=1] = [72, 35, 1, 5, 1].
! In all cases count[1]=nx_g, count[2]=ny_g, count[ndims]=1, and the middle
! count entries come straight from the file's actual dim sizes; total
! element count is preserved so netcdf-fortran reads into the smaller-rank
! buffer correctly. If the file's middle-dim product exceeds the buffer's
! rest-extent the read is partial and the trailing portion of the buffer
! (which the caller initialized to zero before commit) stays zero --
! this is the original FMS behavior for the GDAS Tsnz_h case.
!==============================================================================
subroutine read_global_var(ncid, name, nx_g, ny_g, gbuf3d, gbuf4d)
  integer,                       intent(in)    :: ncid, nx_g, ny_g
  character(len=*),              intent(in)    :: name
  real(kind=kind_real), optional, intent(inout) :: gbuf3d(:,:,:)
  real(kind=kind_real), optional, intent(inout) :: gbuf4d(:,:,:,:)

  integer :: varid, file_ndims, dd, sz
  integer :: file_dimids(8)
  integer, allocatable :: ct(:), st(:)

  call ncc(nf90_inq_varid(ncid, trim(name), varid), 'inq '//trim(name))
  call ncc(nf90_inquire_variable(ncid, varid, ndims=file_ndims, dimids=file_dimids), &
           'inquire '//trim(name))
  if (file_ndims < 3) call mpp_error(FATAL, &
      'soca_io_mod read_global_var: '//trim(name)//' has fewer than 3 file dims')
  allocate(ct(file_ndims), st(file_ndims))
  st = 1
  ct(1) = nx_g
  ct(2) = ny_g
  ct(file_ndims) = 1
  do dd = 3, file_ndims - 1
    call ncc(nf90_inquire_dimension(ncid, file_dimids(dd), len=sz), 'dim '//trim(name))
    ct(dd) = sz
  end do

  if (present(gbuf3d)) then
    call ncc(nf90_get_var(ncid, varid, gbuf3d, start=st, count=ct), 'get '//trim(name))
  else if (present(gbuf4d)) then
    call ncc(nf90_get_var(ncid, varid, gbuf4d, start=st, count=ct), 'get '//trim(name))
  end if
  deallocate(ct, st)
end subroutine read_global_var


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
