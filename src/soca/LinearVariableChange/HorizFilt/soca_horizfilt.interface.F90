! (C) Copyright 2017-2021 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interface for soca_horizfilt_mod::soca_horizfilt
module soca_horizfilt_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use oops_variables_mod, only: oops_variables

! soca modules
use soca_geom_mod_c, only: soca_geom_registry
use soca_geom_mod, only : soca_geom
use soca_horizfilt_mod, only: soca_horizfilt
use soca_increment_mod, only: soca_increment
use soca_increment_reg, only: soca_increment_registry
use soca_state_mod, only: soca_state
use soca_state_reg, only: soca_state_registry

implicit none
private


#define LISTED_TYPE soca_horizfilt

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for soca_horizfilt_mod::soca_horizfilt
type(registry_t), public :: soca_horizfilt_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> C++ interface for soca_horizfilt_mod::soca_horizfilt::setup()
!!
!! Setup for the filtering operator
subroutine soca_horizfilt_setup_c(c_key_self, &
                                  c_conf, &
                                  c_key_geom, &
                                  c_key_traj, &
                                  c_vars) &
        & bind (c,name='soca_horizfilt_setup_f90')
  integer(c_int), intent(inout) :: c_key_self   !< The filtering structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int),    intent(in) :: c_key_geom   !< Geometry
  integer(c_int),    intent(in) :: c_key_traj   !< Trajectory
  type(c_ptr),value, intent(in) :: c_vars       !< List of variables

  type(soca_horizfilt), pointer :: self
  type(soca_geom),           pointer :: geom
  type(soca_state),          pointer :: traj
  type(oops_variables)               :: vars

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_state_registry%get(c_key_traj, traj)
  call soca_horizfilt_registry%init()
  call soca_horizfilt_registry%add(c_key_self)
  call soca_horizfilt_registry%get(c_key_self, self)
  vars = oops_variables(c_vars)
  call self%setup(fckit_configuration(c_conf), geom, traj, vars)

end subroutine soca_horizfilt_setup_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_horizfilt_mod::soca_horizfilt::delete()
!!
!! Delete filtering operator
subroutine soca_horizfilt_delete_c(c_key_self) bind (c,name='soca_horizfilt_delete_f90')
  integer(c_int), intent(inout) :: c_key_self  !< The filtering structure

  type(soca_horizfilt),       pointer :: self

  call soca_horizfilt_registry%get(c_key_self,self)
  call self%delete()
  call soca_horizfilt_registry%remove(c_key_self)

end subroutine soca_horizfilt_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_horizfilt_mod::soca_horizfilt::mult()
!!
!! Multiply
subroutine soca_horizfilt_mult_c(c_key_self, c_key_in, c_key_out, c_key_geom) bind(c,name='soca_horizfilt_mult_f90')
  integer(c_int), intent(inout) :: c_key_self  !< The filtering structure
  integer(c_int), intent(in)    :: c_key_in    !<    "   to Increment in
  integer(c_int), intent(in)    :: c_key_out   !<    "   to Increment out
  integer(c_int), intent(in)    :: c_key_geom  !< Geometry

  type(soca_horizfilt), pointer :: self
  type(soca_increment),      pointer :: xin
  type(soca_increment),      pointer :: xout
  type(soca_geom),           pointer :: geom

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_horizfilt_registry%get(c_key_self, self)
  call soca_increment_registry%get(c_key_in, xin)
  call soca_increment_registry%get(c_key_out, xout)

  call self%mult(xin, xout, geom) !< xout = C.xout

end subroutine soca_horizfilt_mult_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_horizfilt_mod::soca_horizfilt::multad()
!!
!! Multiply adjoint
subroutine soca_horizfilt_mult_ad_c(c_key_self, c_key_in, c_key_out, c_key_geom) &
      bind(c,name='soca_horizfilt_multad_f90')
  integer(c_int), intent(inout) :: c_key_self  !< The filtering structure
  integer(c_int), intent(in)    :: c_key_in    !<    "   to Increment in
  integer(c_int), intent(in)    :: c_key_out   !<    "   to Increment out
  integer(c_int), intent(in)    :: c_key_geom  !< Geometry

  type(soca_horizfilt), pointer :: self
  type(soca_increment),      pointer :: xin
  type(soca_increment),      pointer :: xout
  type(soca_geom),           pointer :: geom

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_horizfilt_registry%get(c_key_self, self)
  call soca_increment_registry%get(c_key_in, xin)
  call soca_increment_registry%get(c_key_out, xout)

  call self%multad(xin, xout, geom) !< xout = C^T.xout

end subroutine soca_horizfilt_mult_ad_c

end module soca_horizfilt_mod_c
