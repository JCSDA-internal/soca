!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
!> Structure holding configuration variables for the  model

module soca_configs

  use kinds
  implicit none
  private
  public :: soca_config
  public :: soca_config_registry

  !> Fortran derived type to hold configuration data for the  model
  type :: soca_config
     integer :: nx     !< Zonal grid dimension
     integer :: ny     !< Meridional grid dimension
     ! dimensional parameters
     real(kind=kind_real) :: dt0       !< dimensional time (seconds)
  end type soca_config

#define LISTED_TYPE soca_config

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_config_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

end module soca_configs
