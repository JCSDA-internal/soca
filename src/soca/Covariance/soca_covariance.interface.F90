! (C) Copyright 2017-2021 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for soca_covariance_mod::soca_cov
module soca_covariance_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use oops_variables_mod, only: oops_variables

! soca modules
use soca_covariance_mod, only: soca_cov
use soca_geom_mod_c, only : soca_geom_registry
use soca_geom_mod, only : soca_geom
use soca_increment_mod, only : soca_increment
use soca_increment_reg, only : soca_increment_registry
use soca_state_mod, only : soca_state
use soca_state_reg, only : soca_state_registry

implicit none
private

#define LISTED_TYPE soca_cov

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for soca_cov
type(registry_t), public :: soca_cov_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> C++ interface for soca_covariance_mod::soca_cov::setup()
!!
!! Setup for the SOCA model's background error covariance matrix
subroutine soca_b_setup_c(c_key_self, c_conf, c_key_geom, c_key_bkg, c_vars) &
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
  call self%setup(fckit_configuration(c_conf), geom, bkg, vars)

end subroutine soca_b_setup_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_covariance_mod::soca_cov::delete()
!!
!! Delete for the SOCA model's background error covariance matrix
subroutine soca_b_delete_c(c_key_self) bind (c,name='soca_b_delete_f90')
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure

  type(soca_cov),       pointer :: self

  call soca_cov_registry%get(c_key_self,self)
  call self%delete()
  call soca_cov_registry%remove(c_key_self)

end subroutine soca_b_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_covariance_mod::soca_cov::mult()
subroutine soca_b_mult_c(c_key_self, c_key_in, c_key_out) bind(c,name='soca_b_mult_f90')
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
  call self%mult(xout) !< xout = C.xout

end subroutine soca_b_mult_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_covariance_mod::soca_cov::sqrt_c_mult()
!!
!! Generate randomized C^1/2 x increment
subroutine soca_b_randomize_c(c_key_self, c_key_out) bind(c,name='soca_b_randomize_f90')
  integer(c_int), intent(in) :: c_key_self  !< covar config structure
  integer(c_int), intent(in) :: c_key_out   !< Randomized increment

  type(soca_cov),       pointer :: self
  type(soca_increment), pointer :: xout

  call soca_cov_registry%get(c_key_self, self)
  call soca_increment_registry%get(c_key_out, xout)

  ! Randomize increment
  call self%sqrt_C_mult(xout) !< xout = C^1/2.xout

end subroutine soca_b_randomize_c

end module soca_covariance_mod_c
