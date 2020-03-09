! (C) Copyright 2017-2020 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_covariance_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use oops_variables_mod
use soca_geom_mod, only : soca_geom
use soca_fields_mod, only: soca_fields
use soca_covariance_mod, only: soca_cov, soca_cov_setup, soca_cov_delete, &
                               soca_cov_C_mult, soca_cov_sqrt_C_mult

implicit none
private

! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------
!> Setup for the SOCA model's background error covariance matrix
subroutine c_soca_b_setup(c_self, c_conf, c_geom, c_bkg, c_vars) &
     & bind (c,name='soca_b_setup_f90')
  type(c_ptr),    intent(inout) :: c_self   !< The background covariance structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  type(c_ptr),       intent(in) :: c_geom       !< Geometry
  type(c_ptr),       intent(in) :: c_bkg    !< Background
  type(c_ptr),value, intent(in) :: c_vars       !< List of variables

  type(soca_cov),   pointer :: self
  type(soca_geom),  pointer :: geom
  type(soca_fields),pointer :: bkg
  type(oops_variables)      :: vars

  allocate(self)
  c_self = c_loc(self)
  call c_f_pointer(c_geom, geom)
  call c_f_pointer(c_bkg, bkg)
  vars = oops_variables(c_vars)

  call soca_cov_setup(self, fckit_configuration(c_conf), geom, bkg, vars)

end subroutine c_soca_b_setup

! ------------------------------------------------------------------------------
!> Delete for the SOCA model's background error covariance matrix

subroutine c_soca_b_delete(c_self) bind (c,name='soca_b_delete_f90')
  type(c_ptr), intent(inout) :: c_self  !< The background covariance structure

  type(soca_cov),       pointer :: self

  call c_f_pointer(c_self, self)

  call soca_cov_delete(self)
  deallocate(self)

end subroutine c_soca_b_delete

! ------------------------------------------------------------------------------

!> Multiply by covariance

subroutine c_soca_b_mult(c_self, c_in, c_out) bind(c,name='soca_b_mult_f90')
  type(c_ptr), intent(inout) :: c_self  !< The background covariance structure
  type(c_ptr), intent(in)    :: c_in    !<    "   to Increment in
  type(c_ptr), intent(in)    :: c_out   !<    "   to Increment out

  type(soca_cov),    pointer :: self
  type(soca_fields), pointer :: xin
  type(soca_fields), pointer :: xout

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_in, xin)
  call c_f_pointer(c_out, xout)

  call xout%copy(xin)              !< xout = xin
  call soca_cov_C_mult(self, xout) !< xout = C.xout

end subroutine c_soca_b_mult


! ------------------------------------------------------------------------------

!> Generate randomized C^1/2 x increment

subroutine c_soca_b_randomize(c_self, c_out) bind(c,name='soca_b_randomize_f90')
  type(c_ptr), intent(in) :: c_self  !< covar config structure
  type(c_ptr), intent(in) :: c_out   !< Randomized increment

  type(soca_cov),   pointer :: self
  type(soca_fields),pointer :: xout

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_out, xout)

  ! Randomize increment
  call soca_cov_sqrt_C_mult(self, xout) !< xout = C^1/2.xout

end subroutine c_soca_b_randomize

end module soca_covariance_mod_c
