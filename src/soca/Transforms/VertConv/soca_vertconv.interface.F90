! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_vertconv_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use soca_fields_mod, only: soca_field, copy
use soca_fields_mod_c, only: soca_field_registry
use soca_vertconv_mod, only: soca_vertconv, soca_conv_setup, &
                             soca_conv, soca_conv_ad

implicit none

private
public :: soca_vertconv_registry

#define LISTED_TYPE soca_vertconv

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: soca_vertconv_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
!> Constructor for Vertconv
subroutine c_soca_vertconv_setup(c_key_self, c_conf, c_key_traj, c_key_bkg) &
  bind(c,name='soca_vertconv_setup_f90')

  integer(c_int), intent(inout) :: c_key_self   !< The Vertconv structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int),    intent(in) :: c_key_traj   !< trajectory
  integer(c_int),    intent(in) :: c_key_bkg    !< background

  type(soca_vertconv), pointer :: self
  type(soca_field), pointer :: traj
  type(soca_field), pointer :: bkg
  real(kind=kind_real), allocatable :: jac(:)

  integer :: isc, iec, jsc, jec, i, j, k, nl

  call soca_vertconv_registry%init()
  call soca_vertconv_registry%add(c_key_self)
  call soca_vertconv_registry%get(c_key_self, self)
  call soca_field_registry%get(c_key_traj, traj)
  call soca_field_registry%get(c_key_bkg, bkg)

  call soca_conv_setup (self, bkg, traj, fckit_configuration(c_conf))

end subroutine c_soca_vertconv_setup

! ------------------------------------------------------------------------------
!> Destructor for Vertconv
subroutine c_soca_vertconv_delete(c_key_self) bind(c,name='soca_vertconv_delete_f90')

  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure

  type(soca_vertconv), pointer :: self

  ! Deallocate trajectory and backgroun
  ! TODO
  ! Deallocate ocean depth array
  ! TODO

  call soca_vertconv_registry%get(c_key_self, self)

  if (associated(self%traj)) nullify(self%traj)
  if (associated(self%bkg)) nullify(self%bkg)

  call soca_vertconv_registry%remove(c_key_self)

end subroutine c_soca_vertconv_delete

! ------------------------------------------------------------------------------
!> Multiplication
subroutine c_soca_vertconv_mult_f90(c_key_a, c_key_m, c_key_traj, c_key_self)&
  bind(c,name='soca_vertconv_mult_f90')

  integer(c_int), intent(in) :: c_key_a     !< Increment in
  integer(c_int), intent(in) :: c_key_m     !< Increment out
  integer(c_int), intent(in) :: c_key_traj  !< trajectory
  integer(c_int), intent(in) :: c_key_self  !< config

  type(soca_field), pointer :: dxa  ! in
  type(soca_field), pointer :: dxm  ! out
  type(soca_field), pointer :: traj
  type(soca_vertconv),   pointer :: self

  call soca_field_registry%get(c_key_a, dxa)
  call soca_field_registry%get(c_key_m, dxm)
  call soca_field_registry%get(c_key_traj, traj)
  call soca_vertconv_registry%get(c_key_self, self)

  !< Computes dxm = Vertconv dxa

  ! dxm = dxa
  call copy(dxm, dxa)

  ! Apply forward convolution operator to T & S
  call soca_conv(self, dxm, dxa)

end subroutine c_soca_vertconv_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_vertconv_multad_f90(c_key_m, c_key_a, c_key_traj, c_key_self)&
  bind(c,name='soca_vertconv_multad_f90')

  integer(c_int), intent(in) :: c_key_a     !< Increment out
  integer(c_int), intent(in) :: c_key_m     !< Increment in
  integer(c_int), intent(in) :: c_key_traj  !< Trajectory
  integer(c_int), intent(in) :: c_key_self  !< config

  type(soca_field),      pointer :: dxa
  type(soca_field),      pointer :: dxm
  type(soca_field),      pointer :: traj
  type(soca_vertconv),   pointer :: self

  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_field_registry%get(c_key_traj,traj)
  call soca_vertconv_registry%get(c_key_self, self)

  ! dxa = dxm
  call copy(dxa,dxm)

  ! Apply adjoint of convolution operator
  call soca_conv_ad(self, dxm, dxa)

end subroutine c_soca_vertconv_multad_f90

end module soca_vertconv_mod_c
