! (C) Copyright 2020-2024 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> State fields
module soca_state_mod

! soca modules
use soca_fields_mod

implicit none
private

!-------------------------------------------------------------------------------
!> State fields.
!!
!! Any procedures that are shared with soca_increment are implemented
!! in the soca_fields base class
type, public, extends(soca_fields) :: soca_state

end type

end module
