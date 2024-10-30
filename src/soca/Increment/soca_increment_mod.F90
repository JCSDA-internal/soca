! (C) Copyright 2020-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Increment fields
module soca_increment_mod

use atlas_module, only: atlas_field
use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables
use random_mod, only: normal_distribution

! soca modules
use soca_fields_mod, only : soca_field, soca_fields
use soca_geom_mod, only : soca_geom


implicit none
private

!-------------------------------------------------------------------------------
!> Increment fields.
!!
!! Any procedures that are shared with soca_state are implemented
!! in the soca_fields base class
type, public, extends(soca_fields) :: soca_increment

contains

  !> \name math operators
  !! \{

  !> \copybrief soca_increment_dirac \see soca_increment_dirac
  procedure :: dirac       => soca_increment_dirac

  !> \copybrief soca_increment_random \see soca_increment_random
  procedure :: random      => soca_increment_random

  !> \}

  !> \name background error decorrelation length scales
  !! \{

  !> \copybrief soca_horiz_scales \see soca_horiz_scales
  procedure :: horiz_scales       => soca_horiz_scales

  !> \copybrief soca_vert_scales \see soca_vert_scales
  procedure :: vert_scales       => soca_vert_scales

  !> \}

end type


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> initialize fields with random normal distribution
!!
!! \note "hocn" field, if present, is NOT randomized, because doing so
!!   causes problems
!! \relates soca_increment_mod::soca_increment
subroutine soca_increment_random(self)
  class(soca_increment), target, intent(inout) :: self

  integer, parameter :: rseed = 1 ! constant for reproducability of tests
    ! NOTE: random seeds are not quite working the way expected,
    !  it is only set the first time normal_distribution() is called with a seed
  integer :: jz, i

  type(soca_field), pointer :: field

  ! set random values
  do i = 1, size(self%fields)
    field => self%fields(i)
    ! TODO remove this once increment / state are fully separated
    ! NOTE: can't randomize "hocn", testIncrementInterpAD fails
    if (field%name == "sea_water_cell_thickness") cycle
    call normal_distribution(field%val,  0.0_kind_real, 1.0_kind_real, rseed)
  end do

  ! mask out land, set to zero
  do i=1,size(self%fields)
    field => self%fields(i)
    if (.not. associated(field%mask) ) cycle
    do jz=1,field%nz
      field%val(:,:,jz) = field%val(:,:,jz) * field%mask(:,:)
    end do
  end do

  ! update domains
  call self%update_halos()
end subroutine soca_increment_random


! ------------------------------------------------------------------------------
!> Apply a dirac increment
!!
!! \raises abor1_ftn aborts if there is an error in the input configuration
!! \todo generalize by removing the hardcoded int=>field_name
!! \relates soca_increment_mod::soca_increment
subroutine soca_increment_dirac(self, f_conf)
  class(soca_increment),        intent(inout) :: self
  type(fckit_configuration), value, intent(in):: f_conf   !< Configuration

  integer :: isc, iec, jsc, jec
  integer :: ndir,n, jz
  integer,allocatable :: ixdir(:),iydir(:),izdir(:),ifdir(:)

  type(soca_field), pointer :: field

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

    ! TODO this list is getting long, change it so that the field name
    ! is directly used in the yaml?
    field => null()
    select case(ifdir(n))
    case (1)
      call self%get("sea_water_potential_temperature", field)
    case (2)
      call self%get("sea_water_salinity", field)
    case (3)
      call self%get("sea_surface_height_above_geoid", field)
    case (4)
      call self%get("sea_ice_area_fraction", field)
    case (5)
      call self%get("sea_ice_thickness", field)
    case (6)
      call self%get("mass_concentration_of_chlorophyll_in_sea_water", field)
    case (7)
      call self%get("molar_concentration_of_biomass_in_sea_water_in_p_units", field)
    case (8)
      call self%get("eastward_sea_water_velocity", field)
    case (9)
      call self%get("northward_sea_water_velocity", field)
    case default
      ! TODO print error that out of range
    end select
    if (associated(field)) then
      jz = 1
      if (field%nz > 1) jz = izdir(n)
      field%val(ixdir(n),iydir(n),izdir(n)) = 1.0
    end if
  end do
end subroutine soca_increment_dirac


! ------------------------------------------------------------------------------
!> compute the horizontal decorelation length scales
!! NOTE: this function should be moved somehwere else, it does not belong in Increment!
!! \relates soca_increment_mod::soca_increment
subroutine soca_horiz_scales(self, f_conf)
  class(soca_increment),        intent(inout) :: self
  type(fckit_configuration), value, intent(in):: f_conf   !< Configuration

  integer :: i, j, jz
  type(fckit_configuration) :: subconf
  real(kind=kind_real) :: r_base, r_mult, r_min_grid, r_min, r_max

  type(atlas_field) :: aField
  real(kind=kind_real), pointer :: aFieldPtr(:,:)
  real(kind=kind_real), allocatable :: rossby(:,:)

  ! get a copy of the rossby radius from atlas
  allocate(rossby(self%geom%isd:self%geom%ied,self%geom%jsd:self%geom%jed))
  rossby = 0.0
  aField = self%geom%fieldset%field("rossby_radius")
  call aField%data(aFieldPtr)
  do j=self%geom%jsc,self%geom%jec
    do i=self%geom%isc,self%geom%iec
      rossby(i,j) = aFieldPtr(1, self%geom%atlas_ij2idx(i,j))
    end do
  end do
  call aField%final()


  ! NOTE, this is duplicated code also present in soca_covariance_mod and possibly elsewhere.
  ! This does not belong in soca_increment_mod and should be moved out

  ! rh is calculated as follows :
  ! 1) rh = "base value" + rossby_radius * "rossby mult"
  ! 2) minimum value of "min grid mult" * grid_size is imposed
  ! 3) min/max are imposed based on "min value" and "max value"
  ! 4) converted from a gaussian sigma to Gaspari-Cohn cutoff distance
  do i=1,size(self%fields)
    ! get parameters for correlation lengths
    call f_conf%get_or_die(trim(self%fields(i)%name), subconf)
    if (.not. subconf%get("base value", r_base)) r_base = 0.0
    if (.not. subconf%get("rossby mult", r_mult)) r_mult = 0.0
    if (.not. subconf%get("min grid mult", r_min_grid)) r_min_grid = 1.0
    if (.not. subconf%get("min value", r_min)) r_min = 0.0
    if (.not. subconf%get("max value", r_max)) r_max = huge(r_max)

    self%fields(i)%val(:,:,1) = r_base + r_mult*rossby(:,:)
    if (r_min_grid .gt. 0.0) then
      self%fields(i)%val(:,:,1) = max(self%fields(i)%val(:,:,1), sqrt(self%geom%cell_area)*r_min_grid)
    end if
    self%fields(i)%val(:,:,1) = min(r_max, self%fields(i)%val(:,:,1))
    self%fields(i)%val(:,:,1) = max(r_min, self%fields(i)%val(:,:,1))
    self%fields(i)%val(:,:,1) = 3.57_kind_real * self%fields(i)%val(:,:,1) ! convert from gaussian sigma to
                                                                           ! Gaspari-Cohn half width

    do jz=2,self%fields(i)%nz
      self%fields(i)%val(:,:,jz) = self%fields(i)%val(:,:,1)
    end do

  end do
end subroutine soca_horiz_scales


! ------------------------------------------------------------------------------
!> compute the vertical decorelation length scales
!!
!! \relates soca_increment_mod::soca_increment
subroutine soca_vert_scales(self, vert)
  class(soca_increment), intent(inout) :: self
  real(kind=kind_real),  intent(in)    :: vert

  integer :: i, jz

  ! compute scales
  do i=1,size(self%fields)
    do jz=1,self%fields(i)%nz
      self%fields(i)%val(:,:,jz) = 3.57_kind_real*self%geom%mask2d(:,:)*vert
    end do
  end do
end subroutine soca_vert_scales
! ------------------------------------------------------------------------------

end module soca_increment_mod
