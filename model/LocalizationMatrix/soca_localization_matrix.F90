! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module SOCALOCALIZATION
use iso_c_binding
use soca_fields
use kinds
use soca_constants
use soca_geom_mod
use config_mod
use soca_covariance_mod

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine soca_localization_mult(c_key_conf, c_key_xincr) bind(c,name='soca_localization_mult_f90')
implicit none
integer(c_int), intent(in) :: c_key_conf
integer(c_int), intent(in) :: c_key_xincr

type(soca_3d_covar_config), pointer :: conf   !< Config structure
type(soca_field), pointer :: xincr
real(kind=kind_real), allocatable :: xctl(:,:,:) ! Control vector

call soca_3d_cov_registry%get(c_key_conf,conf)
call soca_field_registry%get(c_key_xincr,xincr)

!allocate(xctl(conf%nx, conf%ny, 2))

!xctl(:,:,:)=0.0_kind_real
!call soca_3d_covar_sqrt_mult_ad(conf%nx,conf%ny,xincr,xctl,conf)
!call zeros(xincr1)
!call soca_3d_covar_sqrt_mult(conf%nx,conf%ny,xincr,xctl,conf)

!deallocate(xctl)
end subroutine soca_localization_mult

! ------------------------------------------------------------------------------

subroutine soca_localization_setup(c_key_conf, c_model, c_key_geom) bind(c,name='soca_localization_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_conf
type(c_ptr), intent(in)    :: c_model  !< The configuration
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(soca_3d_covar_config), pointer :: conf !< covar structure
type(soca_geom), pointer :: geom     !< Geometry

call soca_3d_cov_registry%init()
call soca_3d_cov_registry%add(c_key_conf)
call soca_3d_cov_registry%get(c_key_conf, conf)
call soca_geom_registry%get(c_key_geom, geom)
call soca_3d_covar_setup(c_model, geom, conf)

return
end subroutine soca_localization_setup

! ------------------------------------------------------------------------------

subroutine soca_localization_delete(c_key_self) bind(c,name='soca_localization_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(soca_3d_covar_config), pointer :: self

call soca_3d_cov_registry%get(c_key_self, self)
!call soca_3d_covar_delete(self)

end subroutine soca_localization_delete

! ------------------------------------------------------------------------------

end module SOCALOCALIZATION
