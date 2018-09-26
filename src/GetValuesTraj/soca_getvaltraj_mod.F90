!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!> Fortran module handling interpolation trajectory

module soca_getvaltraj_mod

  !General JEDI uses
  use kinds
  use iso_c_binding
  use soca_bumpinterp2d_mod
  
  implicit none
  private

  public soca_getvaltraj
  public soca_getvaltraj_registry
  public c_soca_getvaltraj_setup, c_soca_getvaltraj_delete

  type :: soca_getvaltraj
     integer                 :: nobs
     type(soca_bumpinterp2d) :: horiz_interp
     logical                 :: interph_initialized = .false.
     integer                 :: obstype_index
  end type soca_getvaltraj

#define LISTED_TYPE soca_getvaltraj

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_getvaltraj_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine c_soca_getvaltraj_setup(c_key_self) bind(c,name='soca_getvaltraj_setup_f90')

    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(soca_getvaltraj), pointer :: self

    ! Init, add and get key
    ! ---------------------
    call soca_getvaltraj_registry%init()
    call soca_getvaltraj_registry%add(c_key_self)
    call soca_getvaltraj_registry%get(c_key_self,self)

    self%interph_initialized = .false.
    self%nobs = 0
    self%obstype_index = c_key_self

  end subroutine c_soca_getvaltraj_setup

  ! ------------------------------------------------------------------------------

  subroutine c_soca_getvaltraj_delete(c_key_self) bind(c,name='soca_getvaltraj_delete_f90')

    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(soca_getvaltraj), pointer :: self

    ! Get key
    call soca_getvaltraj_registry%get(c_key_self, self)

    if (self%interph_initialized) then
       self%nobs = 0
       self%interph_initialized = .false.      
    endif

    ! Remove key
    call soca_getvaltraj_registry%remove(c_key_self)
    
  end subroutine c_soca_getvaltraj_delete

  ! ------------------------------------------------------------------------------

end module soca_getvaltraj_mod
