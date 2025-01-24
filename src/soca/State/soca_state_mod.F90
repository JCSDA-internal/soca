! (C) Copyright 2020-2024 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> State fields
module soca_state_mod

use logger_mod
use kinds, only: kind_real
use oops_variables_mod
use atlas_module, only: atlas_field

! soca modules
use soca_geom_mod
use soca_fields_mod
use soca_increment_mod

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
