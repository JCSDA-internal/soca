
! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! ------------------------------------------------------------------------------

!> Setup for the SOCA model's background error covariance matrix

subroutine c_soca_b_setup(c_key_self, c_conf, c_key_geom, c_key_bkg) &
     & bind (c,name='soca_b_setup_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_geom_mod
  use soca_fields

  implicit none
  integer(c_int), intent(inout) :: c_key_self   !< The background covariance structure
  type(c_ptr),       intent(in) :: c_conf        !< The configuration
  integer(c_int),    intent(in) :: c_key_geom    !< Geometry
  integer(c_int),    intent(in) :: c_key_bkg     !< Background  

  type(soca_cov),   pointer :: self
  type(soca_geom),  pointer :: geom
  type(soca_field), pointer :: bkg

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_cov_registry%init()
  call soca_cov_registry%add(c_key_self)
  call soca_cov_registry%get(c_key_self, self)
  call soca_field_registry%get(c_key_bkg,bkg)
  call soca_cov_setup(self, c_conf, geom, bkg)

end subroutine c_soca_b_setup

! ------------------------------------------------------------------------------
!> Delete for the SOCA model's background error covariance matrix

subroutine c_soca_b_delete(c_key_self) bind (c,name='soca_b_delete_f90')

  use iso_c_binding
  use soca_covariance_mod

  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure
  type(soca_cov),       pointer :: self

  call soca_cov_registry%get(c_key_self,self)
  call soca_cov_delete(self)
  call soca_cov_registry%remove(c_key_self)
  
end subroutine c_soca_b_delete

! ------------------------------------------------------------------------------

!> Multiply by covariance

subroutine c_soca_b_mult(c_key_self, c_key_in, c_key_out) bind(c,name='soca_b_mult_f90')  
  !> xout = K D C^1/2 C^1/2^T D K^T xin 
  use iso_c_binding
  use soca_covariance_mod
  use soca_fields
  use kinds
  use mpi

  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure
  integer(c_int), intent(in)    :: c_key_in    !<    "   to Increment in
  integer(c_int), intent(in)    :: c_key_out   !<    "   to Increment out 

  type(soca_cov),   pointer :: self
  type(soca_field), pointer :: xin
  type(soca_field), pointer :: xout

  call soca_cov_registry%get(c_key_self, self)
  call soca_field_registry%get(c_key_in, xin)
  call soca_field_registry%get(c_key_out, xout)

  call copy(xout,xin)              !< xout = xin
  call soca_cov_C_mult(self, xout) !< xout = C.xout

end subroutine c_soca_b_mult

! ------------------------------------------------------------------------------

!> Setup linearization parameters (traj, ...)

subroutine c_soca_b_linearize(c_key_self, c_key_geom) bind(c,name='soca_b_linearize_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_geom_mod
  use soca_fields

  implicit none

  integer(c_int), intent(inout) :: c_key_self   !< The trajectory covariance structure
  integer(c_int), intent(in)    :: c_key_geom   !< Geometry

  ! Do nothing, current cov setup does not depend on the
  ! trajectory

end subroutine c_soca_b_linearize

! ------------------------------------------------------------------------------

!> Generate randomized C^1/2 x increment

subroutine c_soca_b_randomize(c_key_self, c_key_out) bind(c,name='soca_b_randomize_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_fields
  use random_vectors_mod
  use kinds

  implicit none

  integer(c_int), intent(in) :: c_key_self  !< covar config structure
  integer(c_int), intent(in) :: c_key_out   !< Randomized increment

  type(soca_cov),   pointer :: self
  type(soca_field), pointer :: xout
  type(soca_field)          :: xtmp  

  call soca_cov_registry%get(c_key_self, self)
  call soca_field_registry%get(c_key_out, xout)

  !call random(xout)                !< xtmp = random

  ! Randomize increment
  call create(xtmp,xout)
  call random(xtmp)                !< xtmp = random

  ! Re-scale increment
  xtmp%tocn = self%pert_scale%T*xtmp%tocn
  xtmp%socn = self%pert_scale%S*xtmp%socn
  xtmp%ssh = self%pert_scale%SSH*xtmp%ssh
  xtmp%cicen = self%pert_scale%AICE*xtmp%cicen  
  xtmp%hicen = self%pert_scale%AICE*xtmp%hicen
  
  ! Apply sqrt convolution to increment
  call soca_cov_sqrt_C_mult(self, xtmp) !< xtmp = C.xtmp
  call copy(xout,xtmp)             !< xout = xtmp
  call delete(xtmp)

end subroutine c_soca_b_randomize

! ------------------------------------------------------------------------------
