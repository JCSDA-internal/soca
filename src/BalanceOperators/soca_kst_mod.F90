!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_kst_mod

  use kinds
  implicit none

  !> Fortran derived type to hold configuration Kst
  type :: soca_kst_config
    real(kind=kind_real) :: dsdtmax !> 1.0 [psu/K]
    real(kind=kind_real) :: dsdzmin !> 3.0e-3 [psu/m] 
    real(kind=kind_real) :: dtdzmin !> 1.0e-3 [K/m]
 end type soca_kst_config

#define LISTED_TYPE soca_kst_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
 type(registry_t) :: soca_kst_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

  subroutine soca_kst_setup(c_conf, config)
    use iso_c_binding
    use config_mod
    use kinds
    
    implicit none
    
    type(c_ptr),              intent(in) :: c_conf   !< The configuration
    type(soca_kst_config), intent(inout) :: config   !< Config parameters for Kst

    config%dsdtmax      = config_get_real(c_conf,"dsdtmax")
    config%dsdzmin      = config_get_real(c_conf,"dsdzmin")
    config%dtdzmin      = config_get_real(c_conf,"dtdzmin")
    print *,'dsdtmax=',config%dsdtmax
    print *,'dsdzmin=',config%dsdzmin
    print *,'dtdzmin=',config%dtdzmin
    read(*,*)
  end subroutine soca_kst_setup

end module soca_kst_mod
