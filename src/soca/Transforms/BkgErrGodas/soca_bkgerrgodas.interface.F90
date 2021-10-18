! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interface for soca_bkgerrgodas_mod::soca_bkgerrgodas
module soca_bkgerrgodas_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration

! soca modules
use soca_bkgerrgodas_mod, only: soca_bkgerrgodas
use soca_geom_mod, only: soca_geom
use soca_geom_mod_c, only: soca_geom_registry
use soca_increment_mod, only: soca_increment
use soca_increment_reg, only: soca_increment_registry
use soca_state_mod, only: soca_state
use soca_state_reg, only: soca_state_registry

implicit none
private

#define LISTED_TYPE soca_bkgerrgodas

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for soca_bkgerrgodas
type(registry_t), public :: soca_bkgerrgodas_registry


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
!> C++ interface for soca_bkgerrgodas_mod::soca_bkgerrgodas::setup()
!!
!! Constructor for D (standard deviation of background error)
subroutine soca_bkgerrgodas_setup_c(c_key_self, c_conf, c_key_bkg, c_key_geom) &
  bind(c,name='soca_bkgerrgodas_setup_f90')

  integer(c_int), intent(inout) :: c_key_self   !< The D structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int), intent(in)    :: c_key_bkg    !< Background field
  integer(c_int), intent(in)    :: c_key_geom    !< Geometry

  type(soca_state), pointer :: bkg
  type(soca_geom), pointer :: geom
  type(soca_bkgerrgodas), pointer :: self

  call soca_bkgerrgodas_registry%init()
  call soca_bkgerrgodas_registry%add(c_key_self)
  call soca_bkgerrgodas_registry%get(c_key_self, self)
  call soca_state_registry%get(c_key_bkg, bkg)
  call soca_geom_registry%get(c_key_geom, geom)

  call self%setup(fckit_configuration(c_conf), bkg, geom)

end subroutine soca_bkgerrgodas_setup_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_bkgerrgodas_mod::soca_bkgerrgodas destructor
subroutine soca_bkgerrgodas_delete_c(c_key_self) &
  bind(c,name='soca_bkgerrgodas_delete_f90')

  integer(c_int), intent(inout) :: c_key_self
  type(soca_bkgerrgodas), pointer :: self

  call soca_bkgerrgodas_registry%get(c_key_self, self)
  if (associated(self%bkg)) nullify(self%bkg)
  call self%std_bkgerr%delete()

  call soca_bkgerrgodas_registry%remove(c_key_self)

end subroutine soca_bkgerrgodas_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_bkgerrgodas_mod::soca_bkgerrgodas::mult()
!!
!! Multiplication forward and adjoint
subroutine soca_bkgerrgodas_mult_c(c_key_self, c_key_a, c_key_m)&
  bind(c,name='soca_bkgerrgodas_mult_f90')

  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out
  integer(c_int), intent(in) :: c_key_self

  type(soca_increment), pointer :: dxa
  type(soca_increment), pointer :: dxm
  type(soca_bkgerrgodas), pointer :: self

  call soca_increment_registry%get(c_key_a,dxa)
  call soca_increment_registry%get(c_key_m,dxm)
  call soca_bkgerrgodas_registry%get(c_key_self,self)

  !< Computes dxm = D dxa
  call dxm%copy(dxa)
  call self%mult(dxa, dxm)

end subroutine soca_bkgerrgodas_mult_c

end module soca_bkgerrgodas_mod_c
