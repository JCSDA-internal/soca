! (C) Copyright 2017-2024 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interface for soca_balance_mod::soca_balance
module soca_balance_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding

! soca modules
use soca_balance_mod, only: soca_balance
use soca_geom_mod_c, only: soca_geom_registry
use soca_geom_mod, only: soca_geom
use soca_increment_mod, only: soca_increment
use soca_increment_reg, only: soca_increment_registry
use soca_state_mod, only: soca_state
use soca_state_reg, only: soca_state_registry

implicit none
private


#define LISTED_TYPE soca_balance

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for soca_balance_mod::soca_balance
type(registry_t), public :: soca_balance_registry


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
!> C++ interface for soca_balance_mod::soca_balance::setup()
subroutine soca_balance_setup_c(c_key_self, c_conf, c_key_traj, c_key_geom) &
  bind(c,name='soca_balance_setup_f90')

  integer(c_int), intent(inout) :: c_key_self   !< The D structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int), intent(in)    :: c_key_traj   !< Background field
  integer(c_int), intent(in)    :: c_key_geom   !< Geometry

  type(soca_state), pointer :: traj
  type(soca_balance), pointer :: self
  type(soca_geom), pointer :: geom

  call soca_balance_registry%init()
  call soca_balance_registry%add(c_key_self)
  call soca_balance_registry%get(c_key_self, self)
  call soca_state_registry%get(c_key_traj, traj)
  call soca_geom_registry%get(c_key_geom, geom)

  call self%setup(fckit_configuration(c_conf), traj, geom)
end subroutine soca_balance_setup_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_balance_mod::soca_balance::delete()
subroutine soca_balance_delete_c(c_key_self) &
  bind(c,name='soca_balance_delete_f90')

  integer(c_int), intent(inout) :: c_key_self

  type(soca_balance), pointer :: self

  call soca_balance_registry%get(c_key_self,self)
  call self%delete()
  call soca_balance_registry%remove(c_key_self)
end subroutine soca_balance_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_balance_mod::soca_balance::mult()
subroutine soca_balance_mult_c(c_key_self, c_key_a, c_key_m)&
  bind(c,name='soca_balance_mult_f90')

  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out
  integer(c_int), intent(in) :: c_key_self

  type(soca_increment), pointer :: dxa
  type(soca_increment), pointer :: dxm
  type(soca_balance), pointer :: self

  call soca_increment_registry%get(c_key_a,dxa)
  call soca_increment_registry%get(c_key_m,dxm)
  call soca_balance_registry%get(c_key_self,self)

  !< Computes dxm = K dxa
  call self%mult(dxa, dxm)
end subroutine soca_balance_mult_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_balance_mod::soca_balance::multinv()
!!
!! Multiplication inverse
subroutine soca_balance_multinv_c(c_key_self, c_key_m, c_key_a)&
  bind(c,name='soca_balance_multinv_f90')

  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out
  integer(c_int), intent(in) :: c_key_self

  type(soca_increment), pointer :: dxa
  type(soca_increment), pointer :: dxm
  type(soca_balance), pointer :: self

  call soca_increment_registry%get(c_key_a,dxa)
  call soca_increment_registry%get(c_key_m,dxm)
  call soca_balance_registry%get(c_key_self,self)

  !< Computes dxa = K^-1 dxm
  call self%multinv(dxa, dxm)
end subroutine soca_balance_multinv_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_balance_mod::soca_balance::multad()
!!
!! Multiplication adjoint
subroutine soca_balance_multad_c(c_key_self, c_key_m, c_key_a)&
  bind(c,name='soca_balance_multad_f90')

  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out
  integer(c_int), intent(in) :: c_key_self

  type(soca_increment), pointer :: dxa
  type(soca_increment), pointer :: dxm
  type(soca_balance), pointer :: self

  call soca_increment_registry%get(c_key_a,dxa)
  call soca_increment_registry%get(c_key_m,dxm)
  call soca_balance_registry%get(c_key_self,self)

  !< Computes dxa = K^T dxm
  call self%multad(dxa, dxm)
end subroutine soca_balance_multad_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_balance_mod::soca_balance::multinvad()
!!
!! Multiplication inverse adjoint
subroutine soca_balance_multinvad_c(c_key_self, c_key_a, c_key_m)&
  bind(c,name='soca_balance_multinvad_f90')

  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out
  integer(c_int), intent(in) :: c_key_self

  type(soca_increment), pointer :: dxa
  type(soca_increment), pointer :: dxm
  type(soca_balance), pointer :: self

  call soca_increment_registry%get(c_key_a,dxa)
  call soca_increment_registry%get(c_key_m,dxm)
  call soca_balance_registry%get(c_key_self,self)

  !< Computes dxm = (K^-1)^T dxa
  call self%multinvad(dxa, dxm)
end subroutine soca_balance_multinvad_c

end module soca_balance_mod_c
