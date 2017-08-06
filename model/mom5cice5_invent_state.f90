

!> Invent an initial state for the QG model.

!> This routine invent an initial state for the QG model. It is used to
!! initialise the "truth run". The initial state consists of a horizontally
!! uniform wind in each layer, with a vertical shear sufficient to produce
!! baroclinic instability. Povided the orography is non-zero and is not
!! symmetrically place in the domain, this is sufficient to generate a
!! non-trivial flow after a few days of integration.
!!
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
  ! ------------------------------------------------------------------------------
  
  return
end subroutine invent_state
