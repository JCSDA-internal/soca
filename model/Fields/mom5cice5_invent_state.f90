

!> Invent an initial state for the model.

!! Two slightly different initial states may be created (according to whether
!! or not ctype is set to 'f').

subroutine invent_state(flds,config)

  use mom5cice5_fields
  use iso_c_binding
  use config_mod
  use fckit_log_module, only : log
  use mom5cice5_constants
  use kinds
  
  implicit none

  type(mom5cice5_field), intent(inout) :: flds    !< Model fields
  type(c_ptr), intent(in)       :: config  !< Configuration structure

  ! ------------------------------------------------------------------------------
  call zeros(flds)
  call random_number(flds%cicen)
  call random_number(flds%vicen)
  call random_number(flds%vsnon)
  ! ------------------------------------------------------------------------------
  
  return
end subroutine invent_state
