! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Invent an initial state for the model.

!! Two slightly different initial states may be created (according to whether
!! or not ctype is set to 'f').

subroutine invent_state(flds)

  use soca_fields
  
  implicit none

  type(soca_field), intent(inout) :: flds    !< Model fields

  ! ------------------------------------------------------------------------------

  call zeros(flds)
  call random_number(flds%cicen)
  

  ! ------------------------------------------------------------------------------
  
  return
end subroutine invent_state
