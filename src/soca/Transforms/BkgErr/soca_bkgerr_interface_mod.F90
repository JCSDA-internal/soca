! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerr_interface_mod

  use iso_c_binding
  use soca_fields_mod, only: soca_field, copy, delete
  use soca_fields_interface_mod, only: soca_field_registry
  use soca_bkgerr_mod, only: soca_bkgerr_config, &
                             soca_bkgerr_setup, soca_bkgerr_mult

implicit none

private
public :: soca_bkgerr_registry

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
!> Constructor for D (standard deviation of background error)
subroutine c_soca_bkgerr_setup(c_key_self, c_conf, c_key_bkg) &
  bind(c,name='soca_bkgerr_setup_f90')

  use fckit_configuration_module, only: fckit_configuration

  integer(c_int), intent(inout) :: c_key_self   !< The D structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int), intent(in)    :: c_key_bkg    !< Background field

  type(soca_field), pointer :: bkg
  type(soca_bkgerr_config), pointer :: self

  call soca_bkgerr_registry%init()
  call soca_bkgerr_registry%add(c_key_self)
  call soca_bkgerr_registry%get(c_key_self, self)
  call soca_field_registry%get(c_key_bkg, bkg)

  call soca_bkgerr_setup(fckit_configuration(c_conf), self, bkg)

end subroutine c_soca_bkgerr_setup

! ------------------------------------------------------------------------------
!> Destructor for D
subroutine c_soca_bkgerr_delete(c_key_self) bind(c,name='soca_bkgerr_delete_f90')

  integer(c_int), intent(inout) :: c_key_self
  type(soca_bkgerr_config), pointer :: self

  call soca_bkgerr_registry%get(c_key_self, self)
  if (associated(self%bkg)) nullify(self%bkg)
  call delete(self%std_bkgerr)

  call soca_bkgerr_registry%remove(c_key_self)

end subroutine c_soca_bkgerr_delete

! ------------------------------------------------------------------------------
!> Multiplication forward and adjoint
subroutine c_soca_bkgerr_mult_f90(c_key_self, c_key_a, c_key_m)&
  bind(c,name='soca_bkgerr_mult_f90')

  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out
  integer(c_int), intent(in) :: c_key_self

  type(soca_field), pointer :: dxa
  type(soca_field), pointer :: dxm
  type(soca_bkgerr_config), pointer :: self

  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_bkgerr_registry%get(c_key_self,self)

  !< Computes dxm = D dxa
  call copy(dxm, dxa)
  call soca_bkgerr_mult(self, dxa, dxm)

end subroutine c_soca_bkgerr_mult_f90

end module soca_bkgerr_interface_mod
