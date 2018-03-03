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
!real(kind=kind_real), allocatable :: xctl(:,:,:) ! Control vector

call soca_3d_cov_registry%get(c_key_conf,conf)
call soca_field_registry%get(c_key_in,xin)
call soca_field_registry%get(c_key_out,xout)

!allocate(xctl(conf%nx, conf%ny, 2))
!xctl(:,:,:)=0.0_kind_real

print *,"[[[[[[[[[[[[[[[[[[[[[[[[ IN B INV MULT ]]]]]]]]]]]]]]]]]]]]]]]]"

!call soca_3d_covar_sqrt_inv_mult(conf%nx,conf%ny,xctl,xin,conf)
!call zeros(xout)
call ones(xout)
call self_schur(xout, xin)
!call soca_3d_covar_sqrt_inv_mult_ad(conf%nx,conf%ny,xctl,xout,conf)

!deallocate(xctl)

end subroutine c_soca_b_inv_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by covariance

subroutine c_soca_b_mult(c_key_conf, c_key_in, c_key_out) bind(c,name='soca_b_mult_f90')

use iso_c_binding
use soca_covariance_mod
use soca_fields
use kinds
use soca_Butils
implicit none
integer(c_int), intent(in) :: c_key_conf  !< covar config structure
integer(c_int), intent(in) :: c_key_in    !< Streamfunction: psi
integer(c_int), intent(in) :: c_key_out   !< Streamfunction: psi
type(soca_3d_covar_config), pointer :: conf
type(soca_field), pointer :: xin
type(soca_field), pointer :: xout
!real(kind=kind_real), allocatable :: xctl(:,:,:) ! Control vector
real(kind=kind_real), allocatable :: dy(:,:,:), Bdy(:,:,:)
integer :: nx, ny, ncat, nk, k

real(kind=kind_real) :: Lx=5.0, Ly=1.0

call soca_3d_cov_registry%get(c_key_conf,conf)
call soca_field_registry%get(c_key_in,xin)
call soca_field_registry%get(c_key_out,xout)

call zeros(xout)
print *,"[[[[[[[[[[[[[[[[[[[[[[[[ IN B MULT ]]]]]]]]]]]]]]]]]]]]]]]]"

nx = xin%geom%ocean%nx
ny = xin%geom%ocean%ny
ncat = xin%geom%ocean%ncat

!cicen
do k=2, 6
   print *,'category:',k
   call gauss(xin%cicen(:,:,k), xout%cicen(:,:,k), xin%geom%ocean%lon, xin%geom%ocean%lat, lx, ly)
end do

!ssh
print *,'ssh'
call gauss(xin%ssh, xout%ssh, xin%geom%ocean%lon, xin%geom%ocean%lat, lx, ly)

end subroutine c_soca_b_mult

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
