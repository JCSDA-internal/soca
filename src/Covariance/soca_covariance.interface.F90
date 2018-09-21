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
  
  type(soca_3d_covar_config), pointer :: self
  type(soca_geom),            pointer :: geom
  type(soca_field),           pointer :: bkg

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_3d_cov_registry%init()
  call soca_3d_cov_registry%add(c_key_self)
  call soca_3d_cov_registry%get(c_key_self, self)
  call soca_field_registry%get(c_key_bkg,bkg)

  call soca_3d_covar_setup(c_conf, geom, self, bkg)

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
  
end subroutine c_soca_b_delete

