! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerrutil_mod

use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use soca_fields_mod, only: soca_fields, soca_field
use soca_utils, only: soca_adjust

implicit none

private
public :: soca_bkgerr_bounds_type, &
          soca_bkgerr_readbounds, soca_bkgerr_applybounds

type :: soca_bkgerr_bounds_type
   real(kind=kind_real) :: t_min, t_max
   real(kind=kind_real) :: s_min, s_max
   real(kind=kind_real) :: ssh_min, ssh_max
   real(kind=kind_real) :: cicen_min, cicen_max
   real(kind=kind_real) :: hicen_min, hicen_max
   real(kind=kind_real) :: chl_min, chl_max
   real(kind=kind_real) :: biop_min, biop_max
 contains
   procedure :: read => soca_bkgerr_readbounds
   procedure :: apply => soca_bkgerr_applybounds
end type soca_bkgerr_bounds_type

contains

! ------------------------------------------------------------------------------
!> Read bounds from config
subroutine soca_bkgerr_readbounds(self, f_conf)
  class(soca_bkgerr_bounds_type), intent(inout) :: self
  type(fckit_configuration),      intent(in)    :: f_conf

  ! Get bounds from configuration
  call f_conf%get_or_die("t_min", self%t_min)
  call f_conf%get_or_die("t_max", self%t_max)
  call f_conf%get_or_die("s_min", self%s_min)
  call f_conf%get_or_die("s_max", self%s_max)
  call f_conf%get_or_die("ssh_min", self%ssh_min)
  call f_conf%get_or_die("ssh_max", self%ssh_max)
  call f_conf%get_or_die("cicen_min", self%cicen_min)
  call f_conf%get_or_die("cicen_max", self%cicen_max)
  call f_conf%get_or_die("hicen_min", self%hicen_min)
  call f_conf%get_or_die("hicen_max", self%hicen_max)
  call f_conf%get_or_die("chl_min", self%chl_min)
  call f_conf%get_or_die("chl_max", self%chl_max)
  call f_conf%get_or_die("biop_min", self%biop_min)
  call f_conf%get_or_die("biop_max", self%biop_max)
end subroutine soca_bkgerr_readbounds

! ------------------------------------------------------------------------------
!> Setup the static background error
subroutine soca_bkgerr_applybounds(self, fld)
  class(soca_bkgerr_bounds_type), intent(inout) :: self
  type(soca_fields),              intent(inout) :: fld

  type(soca_field), pointer :: field

  integer :: isc, iec, jsc, jec, i, j, n
  real(kind=kind_real) :: vmin, vmax

  ! Apply config bounds to background error
  isc = fld%geom%isc ; iec = fld%geom%iec
  jsc = fld%geom%jsc ; jec = fld%geom%jec

  do n=1,size(fld%fields)
    field => fld%fields(n)
    select case(field%name)
    case ("tocn")
      vmin = self%t_min
      vmax = self%t_max
    case ("socn")
      vmin = self%s_min
      vmax = self%s_max
    case ("ssh")
      vmin = self%ssh_min
      vmax = self%ssh_max
    case ("cicen")
      vmin = self%cicen_min
      vmax = self%cicen_max
    case ("hicen")
      vmin = self%hicen_min
      vmax = self%hicen_max
    case ("chl")
      vmin = self%chl_min
      vmax = self%chl_max
    case ("biop")
      vmin = self%biop_min
      vmax = self%biop_max
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
