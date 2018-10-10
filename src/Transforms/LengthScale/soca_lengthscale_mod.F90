!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bkgerr_mod

  use kinds
  implicit none

  !> Fortran derived type to hold configuration D
  type :: soca_bkgerr_config

  end type soca_bkgerr_config

#define LISTED_TYPE soca_bkgerr_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_bkgerr_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

  subroutine soca_bkgerr_setup(c_conf, config)
    use iso_c_binding
    use config_mod
    use kinds

    implicit none

    type(c_ptr),              intent(in) :: c_conf   !< The configuration
    type(soca_bkgerr_config), intent(inout) :: config   !< Config parameters for D

  end subroutine soca_bkgerr_setup
end module soca_bkgerr_mod
