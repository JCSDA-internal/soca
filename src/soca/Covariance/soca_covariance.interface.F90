! (C) Copyright 2017-2020 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_covariance_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use oops_variables_mod
use soca_geom_mod, only : soca_geom
use soca_geom_mod_c, only : soca_geom_registry
use soca_increment_mod
use soca_increment_reg
use soca_state_mod
use soca_state_reg
use soca_covariance_mod, only: soca_cov, soca_cov_setup, soca_cov_delete, &
                               soca_cov_C_mult, soca_cov_sqrt_C_mult

implicit none

private
public :: soca_cov_registry

#define LISTED_TYPE soca_cov

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: soca_cov_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup for the SOCA model's background error covariance matrix

subroutine c_soca_b_setup(c_key_self, c_conf, c_key_geom, c_key_bkg, c_vars) &
     & bind (c,name='soca_b_setup_f90')
  integer(c_int), intent(inout) :: c_key_self   !< The background covariance structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int),    intent(in) :: c_key_geom   !< Geometry
  integer(c_int),    intent(in) :: c_key_bkg    !< Background
  type(c_ptr),value, intent(in) :: c_vars       !< List of variables

  type(soca_cov),   pointer :: self
  type(soca_geom),  pointer :: geom
  type(soca_state), pointer :: bkg
  type(oops_variables)      :: vars

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_cov_registry%init()
  call soca_cov_registry%add(c_key_self)
  call soca_cov_registry%get(c_key_self, self)
  call soca_state_registry%get(c_key_bkg,bkg)
  vars = oops_variables(c_vars)
  call soca_cov_setup(self, fckit_configuration(c_conf), geom, bkg, vars)

end subroutine c_soca_b_setup

! ------------------------------------------------------------------------------
!> Delete for the SOCA model's background error covariance matrix

subroutine c_soca_b_delete(c_key_self) bind (c,name='soca_b_delete_f90')
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure

  type(soca_cov),       pointer :: self

  call soca_cov_registry%get(c_key_self,self)
  call soca_cov_delete(self)
  call soca_cov_registry%remove(c_key_self)

end subroutine c_soca_b_delete

! ------------------------------------------------------------------------------

!> Multiply by covariance

subroutine c_soca_b_mult(c_key_self, c_key_in, c_key_out) bind(c,name='soca_b_mult_f90')
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure
  integer(c_int), intent(in)    :: c_key_in    !<    "   to Increment in
  integer(c_int), intent(in)    :: c_key_out   !<    "   to Increment out

  type(soca_cov),       pointer :: self
  type(soca_increment), pointer :: xin
  type(soca_increment), pointer :: xout

  call soca_cov_registry%get(c_key_self, self)
  call soca_increment_registry%get(c_key_in, xin)
  call soca_increment_registry%get(c_key_out, xout)

  call xout%copy(xin)              !< xout = xin
  call soca_cov_C_mult(self, xout) !< xout = C.xout

end subroutine c_soca_b_mult


! ------------------------------------------------------------------------------

!> Generate randomized C^1/2 x increment

subroutine c_soca_b_randomize(c_key_self, c_key_out) bind(c,name='soca_b_randomize_f90')
  integer(c_int), intent(in) :: c_key_self  !< covar config structure
  integer(c_int), intent(in) :: c_key_out   !< Randomized increment

  type(soca_cov),       pointer :: self
  type(soca_increment), pointer :: xout

  call soca_cov_registry%get(c_key_self, self)
  call soca_increment_registry%get(c_key_out, xout)

  ! Randomize increment
  call soca_cov_sqrt_C_mult(self, xout) !< xout = C^1/2.xout

end subroutine c_soca_b_randomize

end module soca_covariance_mod_c
