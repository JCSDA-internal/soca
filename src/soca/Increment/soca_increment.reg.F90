! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> registry for soca_increment_mod::soca_increment instances for use in
!! Fortran/C++ interface of soca_increment_mod_c
module soca_increment_reg

use soca_increment_mod

implicit none
private

#define LISTED_TYPE soca_increment

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for soca_increment instances
type(registry_t), public :: soca_increment_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

end module soca_increment_reg