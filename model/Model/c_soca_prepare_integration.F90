
subroutine c_soca_prepare_integration(c_key_conf, c_key_state) &
         & bind(c,name='soca_prepare_integration_f90')

use iso_c_binding
use soca_fields
use soca_configs
use soca_constants

implicit none
integer(c_int), intent(in) :: c_key_conf  !< Configuration structure
integer(c_int), intent(in) :: c_key_state !< Model fields

type(soca_config), pointer :: conf
type(soca_field), pointer  :: flds

! ------------------------------------------------------------------------------

call soca_field_registry%get(c_key_state,flds)
call soca_config_registry%get(c_key_conf, conf)

! -- Do stuff below

end subroutine c_soca_prepare_integration
