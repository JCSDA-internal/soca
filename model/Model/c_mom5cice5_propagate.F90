
!> Perform a timestep of the model

subroutine c_mom5cice5_propagate(c_key_conf, c_key_state) bind(c,name='mom5cice5_propagate_f90')

use iso_c_binding
use mom5cice5_fields
use mom5cice5_configs

implicit none
integer(c_int), intent(in) :: c_key_conf  !< Config structure
integer(c_int), intent(in) :: c_key_state !< Model fields

type(mom5cice5_config), pointer :: conf
type(mom5cice5_field),  pointer :: flds

! ------------------------------------------------------------------------------

call mom5cice5_config_registry%get(c_key_conf, conf)
call mom5cice5_field_registry%get(c_key_state,flds)

call propagate(flds, conf)

! ------------------------------------------------------------------------------
return
end subroutine c_mom5cice5_propagate
