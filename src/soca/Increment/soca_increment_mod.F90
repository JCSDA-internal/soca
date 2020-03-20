! (C) Copyright 2020-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_increment_mod

use soca_fields_mod
use soca_geom_iter_mod, only : soca_geom_iter
use kinds, only: kind_real
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use random_mod, only: normal_distribution

implicit none
private

type, public, extends(soca_fields) :: soca_increment

contains
  procedure :: dirac     => soca_increment_dirac
  procedure :: random    => soca_increment_random

  procedure :: getpoint  => soca_increment_getpoint
  procedure :: setpoint  => soca_increment_setpoint
end type


contains

! ------------------------------------------------------------------------------
!> initialize fields with random normal distribution
subroutine soca_increment_random(self)
  class(soca_increment), intent(inout) :: self
  integer, parameter :: rseed = 1 ! constant for reproducability of tests
    ! NOTE: random seeds are not quite working the way expected,
    !  it is only set the first time normal_distribution() is called with a seed
  integer :: z, i

  type(soca_field), pointer :: field

  ! set random values
  do i = 1, size(self%fields)
    field => self%fields(i)
    ! TODO remove this once increment / state are fully separated
    ! NOTE: can't randomize "hocn", testIncrementInterpAD fails
    if (field%name == 'hocn') cycle
    call normal_distribution(field%val,  0.0_kind_real, 1.0_kind_real, rseed)
  end do

  ! mask out land, set to zero
  do i=1,size(self%fields)
    field => self%fields(i)
    if (.not. associated(field%mask) ) cycle
    do z=1,field%nz
      field%val(:,:,z) = field%val(:,:,z) * field%mask(:,:)
    end do
  end do

  ! update domains
  call self%update_halos()
end subroutine soca_increment_random

! ------------------------------------------------------------------------------

subroutine soca_increment_getpoint(self, geoiter, values)
  class(soca_increment), intent(   in) :: self
  type(soca_geom_iter),  intent(   in) :: geoiter
  real(kind=kind_real),  intent(inout) :: values(:)

  integer :: ff, ii, nz
  type(soca_field), pointer :: field

  ! get values
  ! TODO generalize field names
  ii = 0
  do ff = 1, size(self%fields)
    field => self%fields(ff)
    select case(field%name)
    case("tocn", "socn", "ssh", "hocn", "cicen", "hicen","hsnon")
      nz = field%nz
      values(ii+1:ii+nz) = field%val(geoiter%iind, geoiter%jind,:)
      ii = ii + nz
    end select
  end do
end subroutine soca_increment_getpoint

! ------------------------------------------------------------------------------

subroutine soca_increment_setpoint(self, geoiter, values)
  class(soca_increment), intent(inout) :: self
  type(soca_geom_iter),  intent(   in) :: geoiter
  real(kind=kind_real),  intent(   in) :: values(:)

  integer :: ff, ii, nz
  type(soca_field), pointer :: field

  ! Set values
  ! TODO generalize field names
  ii = 0
  do ff = 1, size(self%fields)
    field => self%fields(ff)
    select case(field%name)
    case("tocn", "socn", "ssh", "hocn", "cicen", "hicen","hsnon")
      nz = field%nz
      field%val(geoiter%iind, geoiter%jind,:) = values(ii+1:ii+nz)
      ii = ii + nz
    end select
  end do
end subroutine soca_increment_setpoint


  ! ------------------------------------------------------------------------------
! TODO, generalize by removing the hardcoded int=>field_name
subroutine soca_increment_dirac(self, f_conf)
  class(soca_increment),        intent(inout) :: self
  type(fckit_configuration), value, intent(in):: f_conf   !< Configuration

  integer :: isc, iec, jsc, jec
  integer :: ndir,n, z
  integer,allocatable :: ixdir(:),iydir(:),izdir(:),ifdir(:)
  type(fckit_mpi_comm) :: f_comm

  type(soca_field), pointer :: field

  ! Get MPI communicator
  f_comm = fckit_mpi_comm()

  ! Get Diracs size
  ndir = f_conf%get_size("ixdir")
  if (( f_conf%get_size("iydir") /= ndir ) .or. &
      ( f_conf%get_size("izdir") /= ndir ) .or. &
      ( f_conf%get_size("ifdir") /= ndir )) &
      call abor1_ftn('soca_fields_dirac: inconsistent sizes for ixdir, iydir, izdir, and ifdir')

  ! Allocation
  allocate(ixdir(ndir))
  allocate(iydir(ndir))
  allocate(izdir(ndir))
  allocate(ifdir(ndir))

  ! Get Diracs positions
  call f_conf%get_or_die("ixdir", ixdir)
  call f_conf%get_or_die("iydir", iydir)
  call f_conf%get_or_die("izdir", izdir)
  call f_conf%get_or_die("ifdir", ifdir)

  ! get PE domain bounds
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  ! Setup Diracs
  call self%zeros()
  do n=1,ndir
      ! skip this index if not in the bounds of this PE
      if (ixdir(n) > iec .or. ixdir(n) < isc) cycle
      if (iydir(n) > jec .or. iydir(n) < jsc) cycle

    field => null()
    select case(ifdir(n))
    case (1)
      call self%get("tocn", field)
    case (2)
      call self%get("socn", field)
    case (3)
      call self%get("ssh", field)
    case (4)
      call self%get("cicen", field)
    case (5)
      call self%get("hicen", field)
    case default
      ! TODO print error that out of range
    end select
    if (associated(field)) then
      z = 1
      if (field%nz > 1) z = izdir(n)
      field%val(ixdir(n),iydir(n),izdir(n)) = 1.0
    end if
  end do
end subroutine soca_increment_dirac

end module soca_increment_mod