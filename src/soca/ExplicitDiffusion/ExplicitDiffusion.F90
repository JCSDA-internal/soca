! (C) Copyright 2023-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_explicitdiffusion_c

use iso_c_binding

use soca_geom_mod
use soca_geom_mod_c
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

subroutine soca_explicitdiffusion_setup_c(c_key_self, c_geom) bind(c, name='soca_explicitdiffusion_setup_f90')
  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in) :: c_geom

  type(soca_diffusion), pointer :: self
  type(soca_geom), pointer :: geom

  call soca_diffusion_registry%init()
  call soca_diffusion_registry%add(c_key_self)
  call soca_diffusion_registry%get(c_key_self, self)

  call soca_geom_registry%get(c_geom, geom)

  call self%init(geom)
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_explicitdiffusion_calibrate_c(c_key_self) bind(c, name='soca_explicitdiffusion_calibrate_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(soca_diffusion), pointer :: self

  call soca_diffusion_registry%get(c_key_self, self)

  call self%calibrate() 
end subroutine

! ------------------------------------------------------------------------------
end module