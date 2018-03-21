! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Invent an initial state for the model.

!! Two slightly different initial states may be created (according to whether
!! or not ctype is set to 'f').

subroutine invent_state(flds,config)

  use soca_fields
  use iso_c_binding
  use config_mod
  use fckit_log_module, only : log
  use soca_constants
  use kinds
  
  implicit none

  type(soca_field), intent(inout) :: flds    !< Model fields
  type(c_ptr), intent(in)       :: config  !< Configuration structure

  ! ------------------------------------------------------------------------------

  call zeros(flds)
  call random_number(flds%cicen)
  

  ! ------------------------------------------------------------------------------
  
  return
end subroutine invent_state
