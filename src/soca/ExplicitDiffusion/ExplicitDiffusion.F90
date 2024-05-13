! (C) Copyright 2023-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_explicitdiffusion_c

use iso_c_binding

use atlas_module, only : atlas_fieldset
use fckit_configuration_module, only: fckit_configuration

use soca_geom_mod
use soca_geom_mod_c
use soca_increment_mod
use soca_increment_reg, only : soca_increment_registry
use soca_diffusion_mod, only : soca_diffusion

implicit none
private

#define LISTED_TYPE soca_diffusion
#include "oops/util/linkedList_i.f"
type(registry_t), public :: soca_diffusion_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine soca_explicitdiffusion_setup_c(c_key_self, c_geom, c_conf) bind(c, name='soca_explicitdiffusion_setup_f90')
  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in) :: c_geom
  type(c_ptr),    intent(in) :: c_conf

  type(soca_diffusion), pointer :: self
  type(soca_geom), pointer :: geom

  call soca_diffusion_registry%init()
  call soca_diffusion_registry%add(c_key_self)
  call soca_diffusion_registry%get(c_key_self, self)

  call soca_geom_registry%get(c_geom, geom)

  call self%init(geom, fckit_configuration(c_conf))
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_explicitdiffusion_calibrate_c(c_key_self, c_conf) bind(c, name='soca_explicitdiffusion_calibrate_f90')
  integer(c_int), intent(inout) :: c_key_self
  type(c_ptr),       intent(in) :: c_conf

  type(soca_diffusion), pointer :: self

  call soca_diffusion_registry%get(c_key_self, self)

  call self%calibrate(fckit_configuration(c_conf))
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_explicitdiffusion_multiply_c(c_key_self, c_key_dx, c_sqrt) bind(c, name='soca_explicitdiffusion_multiply_f90')
  integer(c_int), intent(inout)  :: c_key_self
  integer(c_int),    intent(in)  :: c_key_dx
  logical(c_bool),   intent(in)  :: c_sqrt

  type(soca_diffusion), pointer :: self
  type(soca_increment), pointer :: dx
  logical :: sqrt

  sqrt = c_sqrt
  call soca_diffusion_registry%get(c_key_self, self)
  call soca_increment_registry%get(c_key_dx, dx)

  call dx%sync_from_atlas()
  call self%multiply(dx, sqrt)
  call dx%sync_to_atlas()
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_explicitdiffusion_writeparams_c(c_key_self, c_conf) &
    bind (c, name='soca_explicitdiffusion_writeparams_f90')
  integer(c_int),  intent(inout) :: c_key_self
  type(c_ptr),        intent(in) :: c_conf

  type(soca_diffusion), pointer :: self
  type(fckit_configuration) :: f_conf

  f_conf = fckit_configuration(c_conf)
  call soca_diffusion_registry%get(c_key_self, self)
  call self%write_params(f_conf)
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_explicitdiffusion_readparams_c(c_key_self, c_conf) &
    bind (c, name='soca_explicitdiffusion_readparams_f90')
    integer(c_int),  intent(inout) :: c_key_self
    type(c_ptr),        intent(in) :: c_conf

    type(soca_diffusion), pointer :: self
    type(fckit_configuration) :: f_conf

    f_conf = fckit_configuration(c_conf)
    call soca_diffusion_registry%get(c_key_self, self)
    call self%read_params(f_conf)
end subroutine

! ------------------------------------------------------------------------------

end module