! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! ------------------------------------------------------------------------------

!> Setup for the SOCA model's background error covariance matrix

subroutine c_soca_b_setup(c_key_self, c_conf, c_key_geom) &
     & bind (c,name='soca_b_setup_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_geom_mod

  implicit none
  integer(c_int), intent(inout) :: c_key_self   !< The background covariance structure
  type(c_ptr), intent(in)    :: c_conf     !< The configuration
  integer(c_int), intent(in) :: c_key_geom !< Geometry
  type(soca_3d_covar_config), pointer :: self
  type(soca_geom),  pointer :: geom

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_3d_cov_registry%init()
  call soca_3d_cov_registry%add(c_key_self)
  call soca_3d_cov_registry%get(c_key_self, self)

  call soca_3d_covar_setup(c_conf, geom, self)

end subroutine c_soca_b_setup

! ------------------------------------------------------------------------------
!> Delete for the SOCA model's background error covariance matrix

subroutine c_soca_b_delete(c_key_self) bind (c,name='soca_b_delete_f90')

  use iso_c_binding
  use soca_covariance_mod

  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure
  type(soca_3d_covar_config), pointer :: self

  call soca_3d_cov_registry%get(c_key_self,self)
  call soca_3d_covar_delete(c_key_self)
  call soca_3d_cov_registry%remove(c_key_self)

end subroutine c_soca_b_delete

! ------------------------------------------------------------------------------

!> Multiply by inverse of covariance

subroutine c_soca_b_inv_mult(c_key_conf, c_key_in, c_key_out) bind(c,name='soca_b_invmult_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_fields
  use kinds

  implicit none
  integer(c_int), intent(in) :: c_key_conf  !< covar config structure
  integer(c_int), intent(in) :: c_key_in    !< Streamfunction: psi
  integer(c_int), intent(in) :: c_key_out   !< Streamfunction: psi
  type(soca_3d_covar_config), pointer :: conf
  type(soca_field), pointer :: xin
  type(soca_field), pointer :: xout

  call soca_3d_cov_registry%get(c_key_conf,conf)
  call soca_field_registry%get(c_key_in,xin)
  call soca_field_registry%get(c_key_out,xout)

  call zeros(xout)
  call copy(xout,xin) ! HACK, B^-1=Id

end subroutine c_soca_b_inv_mult

! ------------------------------------------------------------------------------

!> Multiply by covariance

subroutine c_soca_b_mult(c_key_conf, c_key_in, c_key_out, c_key_traj) bind(c,name='soca_b_mult_f90')
  ! xout = K D C^1/2 C^1/2^T D K xin 
  use iso_c_binding
  use soca_covariance_mod
  use soca_fields
  use kinds
  use soca_Butils
  implicit none
  integer(c_int), intent(in) :: c_key_conf  !< 
  integer(c_int), intent(in) :: c_key_in    !< 
  integer(c_int), intent(in) :: c_key_out   !< 
  integer(c_int), intent(in) :: c_key_traj  !< 
  type(soca_3d_covar_config), pointer :: conf
  type(soca_field), pointer :: xin
  type(soca_field), pointer :: xout
  type(soca_field), pointer :: traj  
  type(soca_field)          :: xtmp
  type(soca_field)          :: xtmp2  
  integer :: ncat, k, iter

  call soca_3d_cov_registry%get(c_key_conf,conf)
  call soca_field_registry%get(c_key_in,xin)
  call soca_field_registry%get(c_key_out,xout)
  print *,c_key_traj
  !call soca_field_registry%get(c_key_traj,traj)  

  print *,"============ IN B MULT ============="

  call create(xtmp,xin)  
  call copy(xtmp,xin)
  call copy(xout,xin)

  !call soca_3d_covar_K_mult_ad(xtmp, traj)  !xtmp=K^T.xtmp
  !call soca_3d_covar_D_mult(xtmp, conf)     ! xtmp = D.xtmp
  !call soca_3d_covar_mult(xtmp,xout,conf)   ! xout = C.xtmp  
  !call soca_3d_covar_D_mult(xtmp, conf)     ! xtmp = D.xtmp  
  !call soca_3d_covar_K_mult(xtmp, traj)  !xtmp=K^T.xtmp

  !call copy(xout,xtmp)
  
  call delete(xtmp)

end subroutine c_soca_b_mult

! ------------------------------------------------------------------------------

!> Multiply by covariance

subroutine c_soca_b_linearize(c_key_self, c_key_geom) bind(c,name='soca_b_linearize_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_geom_mod
  use soca_fields
  
  implicit none

  integer(c_int), intent(inout) :: c_key_self   !< The trajectory covariance structure
  integer(c_int), intent(in) :: c_key_geom !< Geometry
  type(soca_field), pointer  :: self  !< Trajectory
  type(soca_geom),  pointer  :: geom

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_field_registry%get(c_key_self, self)

  print *,self%ssh
  print *,'=traj'
  print *,'key=',c_key_self
  read(*,*)

end subroutine c_soca_b_linearize

! ------------------------------------------------------------------------------

!> Generate randomized increment

subroutine c_soca_b_randomize(c_key_conf, c_key_out) bind(c,name='soca_b_randomize_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_fields
  use random_vectors_mod
  use kinds

  implicit none
  integer(c_int), intent(in) :: c_key_conf  !< covar config structure
  integer(c_int), intent(in) :: c_key_out   !< Randomized increment
  type(soca_3d_covar_config), pointer :: conf
  type(soca_field), pointer :: xout
  real(kind=kind_real) :: prms ! Control vector
  !real(kind=kind_real), allocatable :: xctl(:,:,:) ! Control vector

  call soca_3d_cov_registry%get(c_key_conf,conf)
  call soca_field_registry%get(c_key_out,xout)

  !allocate(xctl(conf%nx, conf%ny, 2))

  !call random_vector(xctl(:,:,:))
  !call zeros(xout)
  call ones(xout)
  call random(xout)
  call fldrms(xout, prms)

  !call soca_3d_covar_sqrt_mult(conf%nx,conf%ny,xout,xctl,conf)

  !deallocate(xctl)

end subroutine c_soca_b_randomize

! ------------------------------------------------------------------------------
