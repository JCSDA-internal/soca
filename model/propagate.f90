!> Perform a timestep of the model

subroutine propagate(flds,config)

  use mom5cice5_fields
  use mom5cice5_configs
  use mom5cice5_constants
  use kinds

  implicit none
  type(mom5cice5_field),  intent(inout) :: flds
  type(mom5cice5_config), intent(in)    :: config

  print *,'Not implemented: Propagate state with identity ...'

  ! ------------------------------------------------------------------------------
  return
end subroutine propagate
