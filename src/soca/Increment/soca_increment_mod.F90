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
  integer :: i, j, k, n, idx

  type(soca_field), pointer :: field
  type(atlas_field) :: afield
  real(kind=kind_real), pointer :: fdata(:,:)
  real(kind=kind_real), allocatable :: tmp3d(:,:,:)


  do n = 1, self%afieldset%size()
    afield = self%afieldset%field(n)
    field => self%fields(n) ! TODO remove this dependency
    call afield%data(fdata)

    ! NOTE: can't randomize "hocn", testIncrementInterpAD fails
    if (afield%name() == "sea_water_cell_thickness") then
      cycle
    end if

    ! allocate and fill with random values
    ! NOTE this is done in a weird way to keep answers from changing when it was refactored
    allocate(tmp3d(self%geom%isd:self%geom%ied, self%geom%jsd:self%geom%jed, afield%shape(1)))
    call normal_distribution(tmp3d,  0.0_kind_real, 1.0_kind_real, rseed)
    do j=self%geom%jsc,self%geom%jec
      do i=self%geom%isc,self%geom%iec
        idx = self%geom%atlas_ij2idx(i,j)
        do k=1,afield%shape(1)
          fdata(k,idx) = tmp3d(i,j,k)
        end do
      end do
    end do
    deallocate(tmp3d)

    ! mask land
    if (associated(field%mask)) then
      do j=self%geom%jsc,self%geom%jec
        do i=self%geom%isc,self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          do k=1,afield%shape(1)
            fdata(k,idx) = fdata(k,idx) * field%mask(i,j)
          end do
        end do
      end do
    end if

    call afield%set_dirty()
  end do
  call afield%final()
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

  type(atlas_field) :: field
  real(kind=kind_real), pointer :: fdata(:,:)

  ! Define field name mapping
  character(len=:), allocatable :: field_names(:)
  allocate(character(len=50)::field_names(9))
  field_names(1) = "sea_water_potential_temperature"
  field_names(2) = "sea_water_salinity"
  field_names(3) = "sea_surface_height_above_geoid"
  field_names(4) = "sea_ice_area_fraction"
  field_names(5) = "sea_ice_thickness"
  field_names(6) = "mass_concentration_of_chlorophyll_in_sea_water"
  field_names(7) = "molar_concentration_of_biomass_in_sea_water_in_p_units"
  field_names(8) = "eastward_sea_water_velocity"
  field_names(9) = "northward_sea_water_velocity"

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

  ! set all fields to zero
  do n=1,self%afieldset%size()
    field = self%afieldset%field(n)
    call field%data(fdata)
    fdata = 0.0
  end do

  ! Setup Diracs
  do n=1,ndir
    ! skip this index if not in the bounds of this PE
     if (ixdir(n) > iec .or. ixdir(n) < isc) cycle
     if (iydir(n) > jec .or. iydir(n) < jsc) cycle

    ! get field
    if (ifdir(n) <= 0 .or. ifdir(n) > 9) cycle
    field = self%afieldset%field(field_names(ifdir(n)))
    call field%data(fdata)

    ! set dirac
    fdata(izdir(n), self%geom%atlas_ij2idx(ixdir(n),iydir(n))) = 1.0

  end do
  call field%final()
end subroutine soca_increment_dirac


! ------------------------------------------------------------------------------
!> compute the horizontal decorelation length scales
!! NOTE: this function should be moved somehwere else, it does not belong in Increment!
!! \relates soca_increment_mod::soca_increment
subroutine soca_horiz_scales(self, f_conf)
  class(soca_increment),        intent(inout) :: self
  type(fckit_configuration), value, intent(in):: f_conf   !< Configuration

  integer :: n, i, j
  type(fckit_configuration) :: subconf
  real(kind=kind_real) :: r_base, r_mult, r_min_grid, r_min, r_max, val

  type(atlas_field) :: afield, area, rossby
  real(kind=kind_real), pointer :: data_field(:,:), data_area(:,:), data_rossby(:,:)

  ! get a copy of the input atlas fields needed
  rossby = self%geom%fieldset%field("rossby_radius")
  area = self%geom%fieldset%field("area")
  call rossby%data(data_rossby)
  call area%data(data_area)

  ! NOTE, this is duplicated code also present in soca_covariance_mod and possibly elsewhere.
  ! This does not belong in soca_increment_mod and should be moved out

  ! rh is calculated as follows :
  ! 1) rh = "base value" + rossby_radius * "rossby mult"
  ! 2) minimum value of "min grid mult" * grid_size is imposed
  ! 3) min/max are imposed based on "min value" and "max value"
  ! 4) converted from a gaussian sigma to Gaspari-Cohn cutoff distance
  do n=1, self%afieldset%size()
    afield = self%afieldset%field(n)
    call afield%data(data_field)

    ! get parameters for correlation lengths
    call f_conf%get_or_die(trim(afield%name()), subconf)
    if (.not. subconf%get("base value", r_base)) r_base = 0.0
    if (.not. subconf%get("rossby mult", r_mult)) r_mult = 0.0
    if (.not. subconf%get("min grid mult", r_min_grid)) r_min_grid = 1.0
    if (.not. subconf%get("min value", r_min)) r_min = 0.0
    if (.not. subconf%get("max value", r_max)) r_max = huge(r_max)

    do i=1, afield%shape(2)
      val = r_base + r_mult*data_rossby(1, i)
      if (r_min_grid > 0.0) val = max(val, sqrt(data_area(1, i))*r_min_grid)
      val = min(r_max, val)
      val = max(r_min, val)
      val = 3.57_kind_real * val ! convert from gaussian sigma to Gaspari-Cohn half width
      do j=1, afield%shape(1)
        data_field(j, i) = val
      end do
    end do

    call afield%set_dirty(rossby%dirty() .or. area%dirty())
  end do
  call afield%final()
  call rossby%final()
  call area%final()
end subroutine soca_horiz_scales


! ------------------------------------------------------------------------------
!> compute the vertical decorelation length scales
!!
!! \relates soca_increment_mod::soca_increment
subroutine soca_vert_scales(self, vert)
  class(soca_increment), intent(inout) :: self
  real(kind=kind_real),  intent(in)    :: vert

  type(atlas_field) :: field, mask
  real(kind=kind_real), pointer :: data_field(:,:), data_mask(:,:)

  integer :: n, i, k

  ! get a copy of the input atlas fields needed
  mask = self%geom%fieldset%field("mask_h")
  call mask%data(data_mask)

  ! compute scales
  do n=1,self%afieldset%size()
    field=self%afieldset%field(n)
    call field%data(data_field)
    do i=1,field%shape(2)
      do k=1,field%shape(1)
        data_field(k,i) = 3.57_kind_real * data_mask(1,i) * vert
      end do
    end do

  end do
  call field%final()
  call mask%final()
end subroutine soca_vert_scales
! ------------------------------------------------------------------------------

end module soca_increment_mod
