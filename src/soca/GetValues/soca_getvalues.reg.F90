! (C) Copyright 2020-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------
module soca_getvalues_reg

use soca_getvalues_mod

implicit none
private

!> Linked list interface - defines registry_t type
#define LISTED_TYPE soca_getvalues
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t), public:: soca_getvalues_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

end module soca_getvalues_reg
