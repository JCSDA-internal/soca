!> Delete the model

subroutine c_mom5cice5_delete(c_key_conf) bind (c,name='mom5cice5_delete_f90')

  use mom5cice5_configs
  use iso_c_binding

  implicit none
  integer(c_int), intent(inout) :: c_key_conf !< Key to configuration structure
  type(mom5cice5_config), pointer :: conf

  ! ------------------------------------------------------------------------------

  call mom5cice5_config_registry%get(c_key_conf, conf)
  
  call mom5cice5_config_registry%remove(c_key_conf)

  ! ------------------------------------------------------------------------------
  return
end subroutine c_mom5cice5_delete
