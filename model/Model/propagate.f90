!> Perform a timestep of the model

subroutine propagate(flds,config)

  use soca_fields
  use soca_configs
  use soca_constants
  use kinds

  implicit none
  type(soca_field),  intent(inout) :: flds
  type(soca_config), intent(in)    :: config

  flds%cicen=0.9_kind_real*flds%cicen
  print *,'Not implemented: Propagate state with identity ...'

  ! ------------------------------------------------------------------------------
  return
end subroutine propagate
