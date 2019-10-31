! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_ocnsfc_mod

use fckit_configuration_module, only: fckit_configuration
use fms_io_mod, only: fms_io_init, fms_io_exit, &
                      register_restart_field, restart_file_type, &
                      restore_state, free_restart_type
use fms_mod, only: read_data
use kinds, only: kind_real
use MOM_forcing_type, only: forcing
use random_mod, only: normal_distribution
use soca_geom_mod, only: soca_geom

implicit none

private
public :: soca_ocnsfc_type

type :: soca_ocnsfc_type
   type(soca_geom),      pointer     :: geom

   real(kind=kind_real), allocatable :: sw_rad(:,:)
   real(kind=kind_real), allocatable :: lw_rad(:,:)
   real(kind=kind_real), allocatable :: latent_heat(:,:)
   real(kind=kind_real), allocatable :: sens_heat(:,:)
   real(kind=kind_real), allocatable :: fric_vel(:,:)
 contains
   procedure :: create => soca_ocnsfc_create
   procedure :: delete => soca_ocnsfc_delete
   procedure :: zeros => soca_ocnsfc_zeros
   procedure :: abs => soca_ocnsfc_abs
   procedure :: random => soca_ocnsfc_random
   procedure :: copy => soca_ocnsfc_copy
   procedure :: add => soca_ocnsfc_add
   procedure :: schur => soca_ocnsfc_schur
   procedure :: sub => soca_ocnsfc_sub
   procedure :: mul => soca_ocnsfc_mul
   procedure :: axpy => soca_ocnsfc_axpy
   procedure :: diff_incr => soca_ocnsfc_diff_incr
   procedure :: read_restart => soca_ocnsfc_read_rst
   procedure :: read_diag => soca_ocnsfc_read_diag
   procedure :: getforcing => soca_ocnsfc_getforcing
   procedure :: pushforcing => soca_ocnsfc_pushforcing
end type soca_ocnsfc_type

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_create(self, geom)
  class(soca_ocnsfc_type), intent(inout) :: self
  type(soca_geom), pointer,   intent(in) :: geom

  integer :: isd, ied, jsd, jed

  ! Associate geometry
  self%geom => geom

  ! Indices for data domain (with halo)
  isd = geom%isd ; ied = geom%ied
  jsd = geom%jsd ; jed = geom%jed

  ! Allocate ocean state
  if (.not.allocated(self%sw_rad)) allocate(self%sw_rad(isd:ied,jsd:jed))
  if (.not.allocated(self%lw_rad)) allocate(self%lw_rad(isd:ied,jsd:jed))
  if (.not.allocated(self%latent_heat)) allocate(self%latent_heat(isd:ied,jsd:jed))
  if (.not.allocated(self%sens_heat)) allocate(self%sens_heat(isd:ied,jsd:jed))
  if (.not.allocated(self%fric_vel)) allocate(self%fric_vel(isd:ied,jsd:jed))

end subroutine soca_ocnsfc_create

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_delete(self)
  class(soca_ocnsfc_type), intent(inout) :: self

  ! Deallocate all ocean surface state
  deallocate(self%sw_rad)
  deallocate(self%lw_rad)
  deallocate(self%latent_heat)
  deallocate(self%sens_heat)
  deallocate(self%fric_vel)

end subroutine soca_ocnsfc_delete

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_zeros(self)
  class(soca_ocnsfc_type), intent(inout) :: self

  self%sw_rad      = 0.0_kind_real
  self%lw_rad      = 0.0_kind_real
  self%latent_heat = 0.0_kind_real
  self%sens_heat   = 0.0_kind_real
  self%fric_vel    = 0.0_kind_real

end subroutine soca_ocnsfc_zeros

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_abs(self)
  class(soca_ocnsfc_type), intent(inout) :: self

  self%sw_rad      = abs(self%sw_rad)
  self%lw_rad      = abs(self%lw_rad)
  self%latent_heat = abs(self%latent_heat)
  self%sens_heat   = abs(self%sens_heat)
  self%fric_vel    = abs(self%fric_vel)

end subroutine soca_ocnsfc_abs

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_random(self)
  class(soca_ocnsfc_type), intent(inout) :: self

  integer :: rseed = 1

  ! set random values
  call normal_distribution(self%sw_rad, 0.0_kind_real, 1.0_kind_real, rseed)
  call normal_distribution(self%lw_rad, 0.0_kind_real, 1.0_kind_real, rseed)
  call normal_distribution(self%latent_heat, 0.0_kind_real, 1.0_kind_real, rseed)
  call normal_distribution(self%sens_heat, 0.0_kind_real, 1.0_kind_real, rseed)
  call normal_distribution(self%fric_vel, 0.0_kind_real, 1.0_kind_real, rseed)

  ! mask out land, set to zero
  self%sw_rad = self%sw_rad * self%geom%mask2d
  self%lw_rad = self%lw_rad * self%geom%mask2d
  self%latent_heat = self%latent_heat * self%geom%mask2d
  self%sens_heat = self%sens_heat * self%geom%mask2d
  self%fric_vel = self%fric_vel * self%geom%mask2d

end subroutine soca_ocnsfc_random

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_copy(self, rhs)
  class(soca_ocnsfc_type), intent(inout) :: self
  class(soca_ocnsfc_type),    intent(in) :: rhs

  ! associate geometry
  self%geom => rhs%geom

  ! copy fields
  self%sw_rad      = rhs%sw_rad
  self%lw_rad      = rhs%lw_rad
  self%latent_heat = rhs%latent_heat
  self%sens_heat   = rhs%sens_heat
  self%fric_vel    = rhs%fric_vel

end subroutine soca_ocnsfc_copy

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_add(self, other)
  class(soca_ocnsfc_type), intent(inout) :: self
  class(soca_ocnsfc_type),    intent(in) :: other

  self%sw_rad      = self%sw_rad      + other%sw_rad
  self%lw_rad      = self%lw_rad      + other%lw_rad
  self%latent_heat = self%latent_heat + other%latent_heat
  self%sens_heat   = self%sens_heat   + other%sens_heat
  self%fric_vel    = self%fric_vel    + other%fric_vel

end subroutine soca_ocnsfc_add

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_schur(self, other)
  class(soca_ocnsfc_type), intent(inout) :: self
  class(soca_ocnsfc_type),    intent(in) :: other

  self%sw_rad      = self%sw_rad      * other%sw_rad
  self%lw_rad      = self%lw_rad      * other%lw_rad
  self%latent_heat = self%latent_heat * other%latent_heat
  self%sens_heat   = self%sens_heat   * other%sens_heat
  self%fric_vel    = self%fric_vel    * other%fric_vel

end subroutine soca_ocnsfc_schur

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_sub(self, other)
  class(soca_ocnsfc_type), intent(inout) :: self
  class(soca_ocnsfc_type),    intent(in) :: other

  self%sw_rad      = self%sw_rad      - other%sw_rad
  self%lw_rad      = self%lw_rad      - other%lw_rad
  self%latent_heat = self%latent_heat - other%latent_heat
  self%sens_heat   = self%sens_heat   - other%sens_heat
  self%fric_vel    = self%fric_vel    - other%fric_vel

end subroutine soca_ocnsfc_sub

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_mul(self, zz)
  class(soca_ocnsfc_type), intent(inout) :: self
  real(kind=kind_real),       intent(in) :: zz

  self%sw_rad      = zz * self%sw_rad
  self%lw_rad      = zz * self%lw_rad
  self%latent_heat = zz * self%latent_heat
  self%sens_heat   = zz * self%sens_heat
  self%fric_vel    = zz * self%fric_vel

end subroutine soca_ocnsfc_mul

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_axpy(self, zz, other)
  class(soca_ocnsfc_type), intent(inout) :: self
  real(kind=kind_real),       intent(in) :: zz
  class(soca_ocnsfc_type),    intent(in) :: other

  self%sw_rad      = self%sw_rad      + zz * other%sw_rad
  self%lw_rad      = self%lw_rad      + zz * other%lw_rad
  self%latent_heat = self%latent_heat + zz * other%latent_heat
  self%sens_heat   = self%sens_heat   + zz * other%sens_heat
  self%fric_vel    = self%fric_vel    + zz * other%fric_vel

end subroutine soca_ocnsfc_axpy

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_diff_incr(self, x1, x2)
  class(soca_ocnsfc_type), intent(inout) :: self
  class(soca_ocnsfc_type),    intent(in) :: x1
  class(soca_ocnsfc_type),    intent(in) :: x2

  self%sw_rad      = x1%sw_rad      - x2%sw_rad
  self%lw_rad      = x1%lw_rad      - x2%lw_rad
  self%latent_heat = x1%latent_heat - x2%latent_heat
  self%sens_heat   = x1%sens_heat   - x2%sens_heat
  self%fric_vel    = x1%fric_vel    - x2%fric_vel

end subroutine soca_ocnsfc_diff_incr

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_getforcing(self, fluxes)
  class(soca_ocnsfc_type), intent(inout) :: self
  type(forcing),              intent(in) :: fluxes !< Thermodynamic forcing

  ! Get ocnsfc from mom6 forcing
  self%sw_rad      = - real(fluxes%sw, kind=kind_real)
  self%lw_rad      = - real(fluxes%lw, kind=kind_real)
  self%latent_heat = - real(fluxes%latent, kind=kind_real)
  self%sens_heat   = - real(fluxes%sens, kind=kind_real)
  self%fric_vel    = real(fluxes%ustar, kind=kind_real)

end subroutine soca_ocnsfc_getforcing

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_pushforcing(self, fluxes)
  class(soca_ocnsfc_type), intent(in) :: self
  type(forcing),        intent(inout) :: fluxes !< Thermodynamic forcing

  ! Push ocnsfc into mom6 forcing
  fluxes%sw     = - real(self%sw_rad, kind=8)
  fluxes%lw     = - real(self%lw_rad, kind=8)
  fluxes%latent = - real(self%latent_heat, kind=8)
  fluxes%sens   = - real(self%sens_heat, kind=8)
  fluxes%ustar  = real(self%fric_vel, kind=8)

end subroutine soca_ocnsfc_pushforcing

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_read_rst(self, f_conf, geom, fldnames)
  class(soca_ocnsfc_type),   intent(inout) :: self
  type(fckit_configuration), intent(in)    :: f_conf
  type(soca_geom),           intent(in)    :: geom
  character(len=5),          intent(in)    :: fldnames(:)

  integer, parameter :: max_string_length=800
  integer :: idr, i
  character(len=max_string_length) :: filename, basename
  type(restart_file_type) :: restart
  character(len=:), allocatable :: str

  if ( f_conf%has("sfc_filename") ) then
      call f_conf%get_or_die("basename", str)
      basename = str
      call f_conf%get_or_die("sfc_filename", str)
     filename = trim(basename)//trim(str)
     deallocate(str)
  else
     call self%zeros()
     return
  end if

  call fms_io_init()
  do i = 1, size(fldnames)
     select case(fldnames(i))
     case('sw')
        idr = register_restart_field(restart, filename, 'sw_rad', &
                                     self%sw_rad(:,:), &
                                     domain=geom%Domain%mpp_domain)
     case('lw')
        idr = register_restart_field(restart, filename, 'lw_rad', &
                                     self%lw_rad(:,:), &
                                     domain=geom%Domain%mpp_domain)
     case('lhf')
        idr = register_restart_field(restart, filename, 'latent_heat', &
                                     self%latent_heat(:,:), &
                                     domain=geom%Domain%mpp_domain)
     case('shf')
        idr = register_restart_field(restart, filename, 'sens_heat', &
                                     self%sens_heat(:,:), &
                                     domain=geom%Domain%mpp_domain)
     case('us')
        idr = register_restart_field(restart, filename, 'fric_vel', &
                                     self%fric_vel(:,:), &
                                     domain=geom%Domain%mpp_domain)
     end select
  end do
  call restore_state(restart, directory='')
  call free_restart_type(restart)
  call fms_io_exit()

end subroutine soca_ocnsfc_read_rst

! ------------------------------------------------------------------------------
subroutine soca_ocnsfc_read_diag(self, f_conf, geom, fldnames)
  class(soca_ocnsfc_type),   intent(inout) :: self
  type(fckit_configuration), intent(in)    :: f_conf
  type(soca_geom),           intent(in)    :: geom
  character(len=5),          intent(in)    :: fldnames(:)

  integer, parameter :: max_string_length=800
  integer :: i
  character(len=max_string_length) :: filename
  character(len=:), allocatable :: str

  if ( f_conf%has("filename") ) then
      call f_conf%get_or_die("filename", str)
      filename = str
      deallocate(str)
  else
     call self%zeros()
     return
  end if

  call fms_io_init()
  do i = 1, size(fldnames)
     select case(fldnames(i))
     case('sw')
        call read_data(filename,"sw_rad", &
                       self%sw_rad(:,:), &
                       domain=geom%Domain%mpp_domain)
     case('lw')
        call read_data(filename,"lw_rad", &
                       self%lw_rad(:,:), &
                       domain=geom%Domain%mpp_domain)
     case('lhf')
        call read_data(filename,"latent_heat", &
                       self%latent_heat(:,:), &
                       domain=geom%Domain%mpp_domain)
     case('shf')
        call read_data(filename,"sens_heat", &
                       self%sens_heat(:,:), &
                       domain=geom%Domain%mpp_domain)
     case('us')
        call read_data(filename,"fric_vel", &
                       self%fric_vel(:,:), &
                       domain=geom%Domain%mpp_domain)
     end select
  end do
  call fms_io_exit()

end subroutine soca_ocnsfc_read_diag

end module soca_ocnsfc_mod
