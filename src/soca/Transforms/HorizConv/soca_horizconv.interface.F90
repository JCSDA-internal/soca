!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module c_soca_horizconv_mod
  use iso_c_binding
  use soca_horizconv_mod
  use soca_fields_mod_c
  use soca_fields

  implicit none

  private
  public :: soca_horizconv_registry

#define LISTED_TYPE soca_horizconv_type

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_horizconv_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------
  
  ! ------------------------------------------------------------------------------
  !> Constructor for D (standard deviation of background error)
  subroutine c_soca_horizconv_setup(c_key_self, c_conf, c_key_bkg) &
       &bind(c,name='soca_horizconv_setup_f90')
    integer(c_int), intent(inout) :: c_key_self   !< The D structure
    type(c_ptr),       intent(in) :: c_conf       !< The configuration
    integer(c_int), intent(in)    :: c_key_bkg    !< Background field

    type(soca_field), pointer :: bkg
    type(soca_horizconv_type), pointer :: self

    call soca_horizconv_registry%init()
    call soca_horizconv_registry%add(c_key_self)
    call soca_horizconv_registry%get(c_key_self, self)
    call soca_field_registry%get(c_key_bkg, bkg)

    call soca_horizconv_setup(c_conf, self, bkg)

  end subroutine c_soca_horizconv_setup

  ! ------------------------------------------------------------------------------
  !> Destructor for D
  subroutine c_soca_horizconv_delete(c_key_self) bind(c,name='soca_horizconv_delete_f90')
    integer(c_int), intent(inout) :: c_key_self
    type(soca_horizconv_type), pointer :: self

    call soca_horizconv_registry%get(c_key_self, self)
    call soca_horizconv_registry%remove(c_key_self)

  end subroutine c_soca_horizconv_delete

  ! ------------------------------------------------------------------------------
  !> Multiplication forward and adjoint
  subroutine c_soca_horizconv_mult_f90(c_key_self, c_key_a, c_key_m)&
       &bind(c,name='soca_horizconv_mult_f90')
    integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
    integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out 
    integer(c_int), intent(in) :: c_key_self 

    type(soca_field), pointer :: dxa
    type(soca_field), pointer :: dxm
    type(soca_horizconv_type), pointer :: self

    call soca_field_registry%get(c_key_a,dxa)
    call soca_field_registry%get(c_key_m,dxm)
    call soca_horizconv_registry%get(c_key_self,self)  

    !< Computes dxm = D dxa
    call copy(dxm, dxa)
    call soca_horizconv_mult(self, dxa, dxm)

  end subroutine c_soca_horizconv_mult_f90

end module c_soca_horizconv_mod
