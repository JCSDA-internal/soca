! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerrgodas_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use soca_state_mod
use soca_state_reg
use soca_increment_mod
use soca_increment_reg
use soca_bkgerrgodas_mod, only: soca_bkgerrgodas_config, &
                                soca_bkgerrgodas_setup, soca_bkgerrgodas_mult

implicit none

private
public :: soca_bkgerrgodas_registry

#define LISTED_TYPE soca_bkgerrgodas_config

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: soca_bkgerrgodas_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
!> Constructor for D (standard deviation of background error)
subroutine c_soca_bkgerrgodas_setup(c_key_self, c_conf, c_key_bkg) &
  bind(c,name='soca_bkgerrgodas_setup_f90')

  integer(c_int), intent(inout) :: c_key_self   !< The D structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int), intent(in)    :: c_key_bkg    !< Background field

  type(soca_state), pointer :: bkg
  type(soca_bkgerrgodas_config), pointer :: self

  call soca_bkgerrgodas_registry%init()
  call soca_bkgerrgodas_registry%add(c_key_self)
  call soca_bkgerrgodas_registry%get(c_key_self, self)
  call soca_state_registry%get(c_key_bkg, bkg)

  call soca_bkgerrgodas_setup(fckit_configuration(c_conf), self, bkg)

end subroutine c_soca_bkgerrgodas_setup

! ------------------------------------------------------------------------------
!> Destructor for D
subroutine c_soca_bkgerrgodas_delete(c_key_self) bind(c,name='soca_bkgerrgodas_delete_f90')

  integer(c_int), intent(inout) :: c_key_self
  type(soca_bkgerrgodas_config), pointer :: self

  call soca_bkgerrgodas_registry%get(c_key_self, self)
  if (associated(self%bkg)) nullify(self%bkg)
  call self%std_bkgerr%delete()

  call soca_bkgerrgodas_registry%remove(c_key_self)

end subroutine c_soca_bkgerrgodas_delete

! ------------------------------------------------------------------------------
!> Multiplication forward and adjoint
subroutine c_soca_bkgerrgodas_mult_f90(c_key_self, c_key_a, c_key_m)&
  bind(c,name='soca_bkgerrgodas_mult_f90')

  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out
  integer(c_int), intent(in) :: c_key_self

  type(soca_increment), pointer :: dxa
  type(soca_increment), pointer :: dxm
  type(soca_bkgerrgodas_config), pointer :: self

  call soca_increment_registry%get(c_key_a,dxa)
  call soca_increment_registry%get(c_key_m,dxm)
  call soca_bkgerrgodas_registry%get(c_key_self,self)

  !< Computes dxm = D dxa
  call dxm%copy(dxa)
  call soca_bkgerrgodas_mult(self, dxa, dxm)

end subroutine c_soca_bkgerrgodas_mult_f90

end module soca_bkgerrgodas_mod_c
