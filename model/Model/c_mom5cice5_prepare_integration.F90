
subroutine c_mom5cice5_prepare_integration(c_key_conf, c_key_state) &
         & bind(c,name='mom5cice5_prepare_integration_f90')

use iso_c_binding
use mom5cice5_fields
use mom5cice5_configs
use mom5cice5_constants

implicit none
integer(c_int), intent(in) :: c_key_conf  !< Configuration structure
integer(c_int), intent(in) :: c_key_state !< Model fields

type(mom5cice5_config), pointer :: conf
type(mom5cice5_field), pointer  :: flds

! ------------------------------------------------------------------------------

call mom5cice5_field_registry%get(c_key_state,flds)
call mom5cice5_config_registry%get(c_key_conf, conf)

! -- Do stuff below

end subroutine c_mom5cice5_prepare_integration
