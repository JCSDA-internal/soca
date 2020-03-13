! (C) Copyright 2017-2020 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_localization_mod_c

use iso_c_binding
use soca_geom_mod, only: soca_geom
use soca_geom_mod_c, only: soca_geom_registry
use soca_fields_mod, only: soca_field
use soca_fields_mod_c, only: soca_fields_registry
use soca_covariance_mod, only: soca_cov
use soca_covariance_mod_c, only: soca_cov_registry

implicit none

private

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
subroutine soca_localization_mult(c_key_conf, c_key_xincr) &
  bind(c,name='soca_localization_mult_f90')

  integer(c_int), intent(in) :: c_key_conf
  integer(c_int), intent(in) :: c_key_xincr

  type(soca_cov),   pointer :: conf   !< Config structure
  type(soca_field), pointer :: xincr

  call soca_cov_registry%get(c_key_conf,conf)
  call soca_field_registry%get(c_key_xincr,xincr)
  call abor1_ftn("localization: not implemented")

end subroutine soca_localization_mult

! ------------------------------------------------------------------------------
subroutine soca_localization_setup(c_key_conf, c_model, c_key_geom) &
  bind(c,name='soca_localization_setup_f90')

  integer(c_int), intent(inout) :: c_key_conf
  type(c_ptr),    intent(in)    :: c_model    !< The configuration
  integer(c_int), intent(in)    :: c_key_geom !< Geometry

  type(soca_cov),  pointer :: conf !< covar structure
  type(soca_geom), pointer :: geom !< Geometry

  call abor1_ftn("localization: not implemented")

  call soca_cov_registry%init()
  call soca_cov_registry%add(c_key_conf)
  call soca_cov_registry%get(c_key_conf, conf)
  call soca_geom_registry%get(c_key_geom, geom)

  return
end subroutine soca_localization_setup

! ------------------------------------------------------------------------------
subroutine soca_localization_delete(c_key_self) &
  bind(c,name='soca_localization_delete_f90')

  integer(c_int), intent(inout) :: c_key_self

  type(soca_cov), pointer :: self

  call abor1_ftn("localization: not implemented")
  call soca_cov_registry%get(c_key_self, self)

end subroutine soca_localization_delete

end module soca_localization_mod_c
