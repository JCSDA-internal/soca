
!> Setup for the model's background error covariance matrix

subroutine c_mom5cice5_b_setup(c_key_conf, c_model, c_key_geom) &
          & bind (c,name='mom5cice5_b_setup_f90')

use iso_c_binding
use mom5cice5_covariance_mod
use mom5cice5_geom_mod

implicit none
integer(c_int), intent(inout) :: c_key_conf   !< The background covariance structure
type(c_ptr), intent(in)    :: c_model  !< The configuration
integer(c_int), intent(in) :: c_key_geom !< Geometry

type(mom5cice5_3d_covar_config), pointer :: conf
type(mom5cice5_geom),  pointer :: geom

! ------------------------------------------------------------------------------

call mom5cice5_geom_registry%get(c_key_geom, geom)
call mom5cice5_3d_cov_registry%init()
call mom5cice5_3d_cov_registry%add(c_key_conf)
call mom5cice5_3d_cov_registry%get(c_key_conf, conf)

call mom5cice5_3d_covar_setup(c_model, geom, conf)

end subroutine c_mom5cice5_b_setup
