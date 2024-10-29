! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> bakground error bounds
module soca_bkgerrutil_mod

use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use soca_fields_mod, only: soca_fields, soca_field
use soca_utils, only: soca_adjust

implicit none
private

!> bounds for background error
!!
!! Used by soca_bkgerrgodas_mod::soca_bkgerrgodas_config and soca_bkgerr_mod::soca_bkgerr_config
!! \todo cleanup and remove hardcoded variable names
type, public :: soca_bkgerr_bounds_type
  real(kind=kind_real) :: t_min, t_max
  real(kind=kind_real) :: s_min, s_max
  real(kind=kind_real) :: ssh_min, ssh_max
  real(kind=kind_real) :: cicen_min, cicen_max
  real(kind=kind_real) :: hicen_min, hicen_max
  real(kind=kind_real) :: swh_min, swh_max
  real(kind=kind_real) :: chl_min, chl_max
  real(kind=kind_real) :: biop_min, biop_max
contains

  !> \copybrief soca_bkgerr_readbounds \see soca_bkgerr_readbounds
  procedure :: read => soca_bkgerr_readbounds

  !> \copybrief soca_bkgerr_applybounds \see soca_bkgerr_applybounds
  procedure :: apply => soca_bkgerr_applybounds

end type soca_bkgerr_bounds_type

contains


! ------------------------------------------------------------------------------
!> Read bounds from config
!!
!! \relates soca_bkgerrutil_mod::soca_bkgerr_bounds_type
subroutine soca_bkgerr_readbounds(self, f_conf)
  class(soca_bkgerr_bounds_type), intent(inout) :: self
  type(fckit_configuration),      intent(in)    :: f_conf

  ! Get bounds from configuration
  if(.not. f_conf%get("t_min", self%t_min)) self%t_min = 0.0
  if(.not. f_conf%get("t_max", self%t_max)) self%t_min = huge(0.0)
  if(.not. f_conf%get("s_min", self%s_min)) self%s_min = 0.0
  if(.not. f_conf%get("s_max", self%s_max)) self%s_max = huge(0.0)
  if(.not. f_conf%get("ssh_min", self%ssh_min)) self%ssh_min = 0.0
  if(.not. f_conf%get("ssh_max", self%ssh_max)) self%ssh_max = huge(0.0)
  if(.not. f_conf%get("cicen_min", self%cicen_min)) self%cicen_min = 0.0
  if(.not. f_conf%get("cicen_max", self%cicen_max)) self%cicen_max = huge(0.0)
  if(.not. f_conf%get("hicen_min", self%hicen_min)) self%hicen_min = 0.0
  if(.not. f_conf%get("hicen_max", self%hicen_max)) self%hicen_max = huge(0.0)
  if(.not. f_conf%get("chl_min", self%chl_min)) self%chl_min = 0.0
  if(.not. f_conf%get("chl_max", self%chl_max)) self%chl_max =  huge(0.0)
  if(.not. f_conf%get("biop_min", self%biop_min)) self%biop_min = 0.0
  if(.not. f_conf%get("biop_max", self%biop_max)) self%biop_max = huge(0.0)
  if(.not. f_conf%get("swh_min", self%swh_min)) self%swh_min = 0.0
  if(.not. f_conf%get("swh_max", self%swh_max)) self%swh_max = huge(0.0)
end subroutine soca_bkgerr_readbounds


! ------------------------------------------------------------------------------
!> Setup the static background error
!!
!! \relates soca_bkgerrutil_mod::soca_bkgerr_bounds_type
subroutine soca_bkgerr_applybounds(self, fld)
  class(soca_bkgerr_bounds_type), intent(inout) :: self
  type(soca_fields), target,      intent(inout) :: fld !< fields to apply bound to

  type(soca_field), pointer :: field

  integer :: isc, iec, jsc, jec, i, j, n
  real(kind=kind_real) :: vmin, vmax

  ! Apply config bounds to background error
  isc = fld%geom%isc ; iec = fld%geom%iec
  jsc = fld%geom%jsc ; jec = fld%geom%jec

  do n=1,size(fld%fields)
    field => fld%fields(n)
    select case(field%name)
    case ("sea_water_potential_temperature")
      vmin = self%t_min!! \relates soca_bkgerrutil_mod::soca_bkgerr_bounds_type
      vmax = self%t_max
    case ("sea_water_salinity")
      vmin = self%s_min
      vmax = self%s_max
    case ("sea_surface_height_above_geoid")
      vmin = self%ssh_min
      vmax = self%ssh_max
    case ("sea_ice_area_fraction")
      vmin = self%cicen_min
      vmax = self%cicen_max
    case ("sea_ice_thickness")
      vmin = self%hicen_min
      vmax = self%hicen_max
    case ("mass_concentration_of_chlorophyll_in_sea_water")
      vmin = self%chl_min
      vmax = self%chl_max
    case ("molar_concentration_of_biomass_in_sea_water_in_p_units")
      vmin = self%biop_min
      vmax = self%biop_max
    case ("sea_surface_wave_significant_height")
      vmin = self%swh_min
      vmax = self%swh_max
    case default
      cycle
    end select

    do i = isc, iec
      do j = jsc, jec
        field%val(i,j,:) = soca_adjust(field%val(i,j,:), vmin, vmax)
      end do
    end do
  end do

end subroutine soca_bkgerr_applybounds

end module soca_bkgerrutil_mod
