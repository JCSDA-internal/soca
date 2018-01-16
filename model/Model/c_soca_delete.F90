!> Delete the model

subroutine c_soca_delete(c_key_conf) bind (c,name='soca_delete_f90')

  use soca_configs
  use iso_c_binding

  implicit none
  integer(c_int), intent(inout) :: c_key_conf !< Key to configuration structure
  type(soca_config), pointer :: conf

  ! ------------------------------------------------------------------------------
  print *,'=============== model delete'
  call soca_config_registry%get(c_key_conf, conf)  
  call soca_config_registry%remove(c_key_conf)



  ! ------------------------------------------------------------------------------
  return
end subroutine c_soca_delete
