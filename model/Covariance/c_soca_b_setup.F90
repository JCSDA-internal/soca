! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Setup for the model's background error covariance matrix

subroutine c_soca_b_setup(c_key_conf, c_model, c_key_geom) &
          & bind (c,name='soca_b_setup_f90')

use iso_c_binding
use soca_covariance_mod
use soca_geom_mod

implicit none
integer(c_int), intent(inout) :: c_key_conf   !< The background covariance structure
type(c_ptr), intent(in)    :: c_model  !< The configuration
integer(c_int), intent(in) :: c_key_geom !< Geometry

type(soca_3d_covar_config), pointer :: conf
type(soca_geom),  pointer :: geom

! ------------------------------------------------------------------------------

call soca_geom_registry%get(c_key_geom, geom)
call soca_3d_cov_registry%init()
call soca_3d_cov_registry%add(c_key_conf)
call soca_3d_cov_registry%get(c_key_conf, conf)

call soca_3d_covar_setup(c_model, geom, conf)

end subroutine c_soca_b_setup
