! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Structure holding configuration variables for the 3d error
!! covariance matrices of the SOCA analysis.

module soca_covariance_mod

use kinds
implicit none

!> Fortran derived type to hold configuration data for the SOCA background/model covariance
type :: soca_3d_covar_config
  integer :: nx !< Zonal grid dimension
  integer :: ny !< Meridional grid dimension
  real(kind=kind_real)    :: sigma        !< Standard deviation
  real(kind=kind_real)    :: vert_corr    !< Vertical correlation between levels
  real(kind=kind_real), allocatable :: sqrt_merid(:,:) !< sqrt(meridional correlation matrix)
  real(kind=kind_real), allocatable :: sqrt_inv_merid(:,:) !< Its inverse
  real(kind=kind_real), allocatable :: sqrt_zonal(:) !< Spectral weights for sqrt(zonal corr)
end type soca_3d_covar_config

#define LISTED_TYPE soca_3d_covar_config

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: soca_3d_cov_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "util/linkedList_c.f"
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

!> Setup for the SOCA model's 3d error covariance matrices (B and Q_i)

!> This routine queries the configuration for the parameters that define the
!! covariance matrix, and stores the relevant values in the
!! error covariance structure.

subroutine soca_3d_covar_setup(c_model, geom, config)

use soca_constants
use soca_geom_mod
use iso_c_binding
use config_mod
use fft_mod
use kinds
use fckit_log_module, only : fckit_log

implicit none
type(c_ptr), intent(in)   :: c_model  !< The configuration
type(soca_geom), intent(in) :: geom     !< Geometry
type(soca_3d_covar_config), intent(inout) :: config !< The covariance structure
real(kind=kind_real) :: corr_length_scale

!config%nx         = geom%nx
!config%ny         = geom%ny
config%sigma      = config_get_real(c_model,"standard_deviation")
config%vert_corr  = config_get_real(c_model,"vertical_correlation")
corr_length_scale = config_get_real(c_model,"horizontal_length_scale")

return
end subroutine soca_3d_covar_setup

! ------------------------------------------------------------------------------

!> Delete for the SOCA model's 3d error covariance matrices
!!$
subroutine soca_3d_covar_delete(c_key_conf)

use iso_c_binding

implicit none
integer(c_int), intent(inout) :: c_key_conf !< The model covariance structure

type(soca_3d_covar_config), pointer :: conf !< covar structure

call soca_3d_cov_registry%get(c_key_conf, conf)

!deallocate(conf%sqrt_zonal)
!deallocate(conf%sqrt_merid)
!deallocate(conf%sqrt_inv_merid)
call soca_3d_cov_registry%remove(c_key_conf)

end subroutine soca_3d_covar_delete

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse(sqrt(C)), where C is 3d covariance matrix

!!$subroutine soca_3d_covar_sqrt_inv_mult(kx,ky,xctl,xincr,config)
!!$use iso_c_binding
!!$use fft_mod
!!$use kinds
!!$use soca_fields
!!$
!!$implicit none
!!$integer(c_int), intent(in)    :: kx            !< Zonal grid dimension
!!$integer(c_int), intent(in)    :: ky            !< Meridional grid dimension
!!$real(c_double), intent(inout) :: xctl(kx,ky,2) !< inv(sqrt(C)) times psi
!!$type(soca_field), intent(in)    :: xincr         !< Streamfunction: psi
!!$type(soca_3d_covar_config), intent(in) :: config !< covar config structure
!!$
!!$real(kind=kind_real) :: zfour(kx+2), work(ky)
!!$integer :: i, j, k, iri, m
!!$real(kind=kind_real) :: zc, zero, one
!!$
!!$!--- multiply by standard deviation
!!$
!!$zc = 1.0_kind_real/config%sigma
!!$do k=1,2
!!$  do j=1,ky
!!$    do i=1,kx
!!$      xctl(i,j,k) = zc * xincr%x(i,j,k)
!!$    enddo
!!$  enddo
!!$enddo
!!$
!!$!--- multiply by inverse square-root of zonal correlation matrix
!!$
!!$do k=1,2
!!$  do j=1,ky
!!$    call fft_fwd(kx,xctl(:,j,k),zfour)
!!$    do m=0,kx/2
!!$      do iri=1,2
!!$        zfour(2*m+iri) = zfour(2*m+iri) / config%sqrt_zonal(m)
!!$      enddo
!!$    enddo
!!$    call fft_inv(kx,zfour,xctl(:,j,k))
!!$  enddo
!!$enddo
!!$
!!$!--- multiply by inverse square-root of meridional correlation matrix
!!$
!!$zero = 0.0_kind_real
!!$one = 1.0_kind_real
!!$do k=1,2
!!$  do i=1,kx
!!$    call DSYMV('L',ky,one,config%sqrt_inv_merid,ky,xctl(i,1,k),kx,zero,work,1)
!!$    do j=1,ky
!!$      xctl(i,j,k) = work(j)
!!$    enddo
!!$  enddo
!!$enddo
!!$
!!$!--- multiply by inverse symmetric square-root of vertical correlation matrix
!!$
!!$zc = 1.0_kind_real / sqrt(1.0_kind_real-config%vert_corr*config%vert_corr)
!!$do j=1,ky
!!$  do i=1,kx
!!$    xctl(i,j,2) = zc * (xctl(i,j,2) - config%vert_corr * xctl(i,j,1))
!!$  enddo
!!$enddo
!!$
!!$end subroutine soca_3d_covar_sqrt_inv_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse(sqrt(C)) - Adjoint

!!$subroutine soca_3d_covar_sqrt_inv_mult_ad(kx,ky,xctl,xincr,config)
!!$use iso_c_binding
!!$use fft_mod
!!$use kinds
!!$use soca_fields
!!$
!!$implicit none
!!$integer(c_int), intent(in) :: kx               !< Zonal grid dimension
!!$integer(c_int), intent(in) :: ky               !< Meridional grid dimension
!!$type(soca_field), intent(inout) :: xincr         !< sqrt(C) times streamfunction
!!$real(c_double), intent(in) :: xctl(kx,ky,2)    !< Streamfunction
!!$type(soca_3d_covar_config), intent(in) :: config !< covariance config structure
!!$
!!$real(kind=kind_real), allocatable :: xout(:,:,:)
!!$real(kind=kind_real) :: zfour(kx+2), work(ky)
!!$integer :: i, j, k, iri, m
!!$real(kind=kind_real) :: zc, zero, one
!!$
!!$!--- adoint of multiplication by inverse symmetric square-root of vertical
!!$!--- correlation matrix
!!$
!!$allocate(xout(kx,ky,2))
!!$
!!$zc = sqrt(1.0_kind_real-config%vert_corr*config%vert_corr)
!!$do j=1,ky
!!$  do i=1,kx
!!$    xout(i,j,1) = xctl(i,j,1) -config%vert_corr * xctl(i,j,2) &
!!$                    & *(1.0_kind_real/zc)
!!$    xout(i,j,2) = xctl(i,j,2)*(1.0_kind_real/zc)
!!$  enddo
!!$enddo
!!$
!!$!--- adjoint multiplication by inverse sqrt of meridional correlation matrix
!!$
!!$zero = 0.0_kind_real
!!$one = 1.0_kind_real
!!$do k=1,2
!!$  do i=1,kx
!!$    call DSYMV('L',ky,one,config%sqrt_inv_merid,ky,xout(i,1,k), &
!!$        &      kx,zero,work,1)
!!$    do j=1,ky
!!$      xout(i,j,k) = work(j)
!!$    enddo
!!$  enddo
!!$enddo
!!$
!!$!--- multiply by inverse square-root of zonal correlation matrix
!!$
!!$do k=1,2
!!$  do j=1,ky
!!$    call fft_fwd(kx,xout(:,j,k),zfour)
!!$    do m=0,kx/2
!!$      do iri=1,2
!!$        zfour(2*m+iri) = zfour(2*m+iri) / config%sqrt_zonal(m)
!!$      enddo
!!$    enddo
!!$    call fft_inv(kx,zfour,xout(:,j,k))
!!$  enddo
!!$enddo
!!$
!!$!--- adjoint of multiplication by standard deviation
!!$
!!$do k=1,2
!!$  do j=1,ky
!!$    do i=1,kx
!!$      xincr%x(i,j,k) = xincr%x(i,j,k) + xout(i,j,k) &
!!$                    & *(1.0_kind_real/config%sigma)
!!$    enddo
!!$  enddo
!!$enddo
!!$
!!$deallocate(xout)
!!$
!!$end subroutine soca_3d_covar_sqrt_inv_mult_ad

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C), where C is a 3d covariance matrix

!!$subroutine soca_3d_covar_sqrt_mult(kx,ky,xincr,xctl,config)
!!$use iso_c_binding
!!$use fft_mod
!!$use kinds
!!$use soca_fields
!!$
!!$implicit none
!!$integer(c_int), intent(in) :: kx               !< Zonal grid dimension
!!$integer(c_int), intent(in) :: ky               !< Meridional grid dimension
!!$type(soca_field), intent(inout) :: xincr         !< sqrt(C) times streamfunction
!!$real(c_double), intent(in) :: xctl(kx,ky,2)    !< Streamfunction
!!$type(soca_3d_covar_config), intent(in) :: config !< covariance config structure
!!$
!!$integer :: i, j, k, iri, m
!!$real(kind=kind_real) :: zc, zero, one
!!$real(kind=kind_real) :: zfour(kx+2), work(ky)
!!$
!!$!--- multiply by symmetric square-root of vertical correlation matrix
!!$
!!$zc = sqrt(1.0_kind_real-config%vert_corr*config%vert_corr)
!!$do j=1,ky
!!$  do i=1,kx
!!$    xincr%x(i,j,1) = xctl(i,j,1)
!!$    xincr%x(i,j,2) = config%vert_corr * xctl(i,j,1) + zc * xctl(i,j,2)
!!$  enddo
!!$enddo
!!$
!!$!--- multiply by square-root of meridional correlation matrix
!!$
!!$zero = 0.0_kind_real
!!$one = 1.0_kind_real
!!$do k=1,2
!!$  do i=1,kx
!!$    call DSYMV('L',ky,one,config%sqrt_merid,ky,xincr%x(i,1,k),kx,zero,work,1)
!!$    do j=1,ky
!!$      xincr%x(i,j,k) = work(j)
!!$    enddo
!!$  enddo
!!$enddo
!!$
!!$!--- multiply by square-root of zonal correlation matrix
!!$
!!$do k=1,2
!!$  do j=1,ky
!!$    call fft_fwd(kx,xincr%x(:,j,k),zfour)
!!$    do m=0,kx/2
!!$      do iri=1,2
!!$        zfour(2*m+iri) = zfour(2*m+iri) * config%sqrt_zonal(m)
!!$      enddo
!!$    enddo
!!$    call fft_inv(kx,zfour,xincr%x(:,j,k))
!!$  enddo
!!$enddo
!!$
!!$!--- multiply by standard deviation
!!$
!!$do k=1,2
!!$  do j=1,ky
!!$    do i=1,kx
!!$      xincr%x(i,j,k) = xincr%x(i,j,k) * config%sigma
!!$    enddo
!!$  enddo
!!$enddo
!!$
!!$end subroutine soca_3d_covar_sqrt_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C) - Adjoint

!!$subroutine soca_3d_covar_sqrt_mult_ad(kx,ky,xincr,xctl,config)
!!$use iso_c_binding
!!$use kinds
!!$use soca_fields
!!$
!!$implicit none
!!$integer(c_int), intent(in)    :: kx            !< Zonal grid spacing
!!$integer(c_int), intent(in)    :: ky            !< Meridional grid spacing
!!$real(c_double), intent(inout) :: xctl(kx,ky,2) !< Result
!!$type(soca_field), intent(in)    :: xincr         !< Streamfunction: psi
!!$type(soca_3d_covar_config), intent(in) :: config !< covar config structure
!!$
!!$real(kind=kind_real), allocatable :: xout(:,:,:)
!!$integer :: i, j, k, iri, m
!!$real(kind=kind_real) :: zc, zero, one
!!$real(kind=kind_real) :: zgrid(kx), zfour(kx+2), work(ky)
!!$
!!$!--- adjoint of multiplication by standard deviation
!!$
!!$allocate(xout(kx,ky,2))
!!$
!!$do k=1,2
!!$  do j=1,ky
!!$    zgrid(:) = xincr%x(:,j,k) * config%sigma
!!$    call fft_fwd(kx,zgrid,zfour)
!!$    do m=0,kx/2
!!$      do iri=1,2
!!$        zfour(2*m+iri) = zfour(2*m+iri) * config%sqrt_zonal(m)
!!$      enddo
!!$    enddo
!!$    call fft_inv(kx,zfour,xout(:,j,k))
!!$  enddo
!!$enddo
!!$
!!$!--- adjoint of multiplication by square-root of meridional correlation matrix
!!$
!!$zero = 0.0_kind_real
!!$one = 1.0_kind_real
!!$do k=1,2
!!$  do i=1,kx
!!$    call DSYMV('L',ky,ONE,CONFIg%sqrt_merid,ky,xout(i,1,k),kx,zero,work,1)
!!$    do j=1,ky
!!$      xout(i,j,k) = work(j)
!!$    enddo
!!$  enddo
!!$enddo
!!$
!!$!--- adjoint multiplication by symmetric square-root of vert correlation matrix
!!$
!!$zc = sqrt(1.0_kind_real-config%vert_corr*config%vert_corr)
!!$do j=1,ky
!!$  do i=1,kx
!!$    xctl(i,j,1) = xctl(i,j,1) + xout(i,j,1) &
!!$              & +config%vert_corr * xout(i,j,2)
!!$    xctl(i,j,2) = xctl(i,j,2) + xout(i,j,2) * zc
!!$  enddo
!!$enddo
!!$
!!$deallocate(xout)
!!$
!!$end subroutine soca_3d_covar_sqrt_mult_ad

! ------------------------------------------------------------------------------

end module soca_covariance_mod
