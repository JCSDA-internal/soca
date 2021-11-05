! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for soca_vertconv_mod::soca_vertconv
module soca_vertconv_mod_c

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds, only: kind_real

! soca modules
use soca_geom_mod_c, only: soca_geom_registry
use soca_geom_mod, only: soca_geom
use soca_increment_mod, only: soca_increment
use soca_increment_reg, only: soca_increment_registry
use soca_state_mod, only: soca_state
use soca_state_reg, only: soca_state_registry
use soca_vertconv_mod, only: soca_vertconv

implicit none
private

#define LISTED_TYPE soca_vertconv

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for soca_vertconv
type(registry_t), public :: soca_vertconv_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"


! ------------------------------------------------------------------------------
!> C++ interface for soca_vertconv_mod::soca_vertconv::setup()
!!
!! Constructor for Vertconv
subroutine soca_vertconv_setup_c(c_key_self, c_conf, c_key_bkg, c_key_geom) &
  bind(c,name='soca_vertconv_setup_f90')

  integer(c_int), intent(inout) :: c_key_self   !< The Vertconv structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int),    intent(in) :: c_key_bkg    !< background
  integer(c_int),    intent(in) :: c_key_geom   !< geometry

  type(soca_vertconv), pointer :: self
  type(soca_state), pointer :: bkg
  type(soca_geom), pointer :: geom

  call soca_vertconv_registry%init()
  call soca_vertconv_registry%add(c_key_self)
  call soca_vertconv_registry%get(c_key_self, self)
  call soca_state_registry%get(c_key_bkg, bkg)
  call soca_geom_registry%get(c_key_geom, geom)

  call self%setup(bkg, geom, fckit_configuration(c_conf))

end subroutine soca_vertconv_setup_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_vertconv_mod::soca_vertconv destructor
!!
!! Destructor for Vertconv
subroutine soca_vertconv_delete_c(c_key_self) bind(c,name='soca_vertconv_delete_f90')

  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure

  type(soca_vertconv), pointer :: self

  ! Deallocate background
  ! TODO
  ! Deallocate ocean depth array
  ! TODO

  call soca_vertconv_registry%get(c_key_self, self)

  if (associated(self%bkg)) nullify(self%bkg)

  call soca_vertconv_registry%remove(c_key_self)

end subroutine soca_vertconv_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_vertconv_mod::soca_vertconv::mult()
!!
!! Multiplication
subroutine soca_vertconv_mult_c(c_key_a, c_key_m, c_key_self)&
  bind(c,name='soca_vertconv_mult_f90')

  integer(c_int), intent(in) :: c_key_a     !< Increment in
  integer(c_int), intent(in) :: c_key_m     !< Increment out
  integer(c_int), intent(in) :: c_key_self  !< config

  type(soca_increment), pointer :: dxa  ! in
  type(soca_increment), pointer :: dxm  ! out
  type(soca_vertconv),  pointer :: self

  call soca_increment_registry%get(c_key_a, dxa)
  call soca_increment_registry%get(c_key_m, dxm)
  call soca_vertconv_registry%get(c_key_self, self)

  !< Computes dxm = Vertconv dxa

  ! dxm = dxa
  call dxm%copy( dxa)

  ! Apply forward convolution operator to T & S
  call self%mult(dxm, dxa)

end subroutine soca_vertconv_mult_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_vertconv_mod::soca_vertconv::mult_ad()
!!
!! Multiplication adjoint
subroutine soca_vertconv_multad_c(c_key_m, c_key_a, c_key_self)&
  bind(c,name='soca_vertconv_multad_f90')

  integer(c_int), intent(in) :: c_key_a     !< Increment out
  integer(c_int), intent(in) :: c_key_m     !< Increment in
  integer(c_int), intent(in) :: c_key_self  !< config

  type(soca_increment),   pointer :: dxa
  type(soca_increment),   pointer :: dxm
  type(soca_vertconv),    pointer :: self

  call soca_increment_registry%get(c_key_a,dxa)
  call soca_increment_registry%get(c_key_m,dxm)
  call soca_vertconv_registry%get(c_key_self, self)

  ! dxa = dxm
  call dxa%copy(dxm)

  ! Apply adjoint of convolution operator
  call self%mult_ad(dxm, dxa)

end subroutine soca_vertconv_multad_c

end module soca_vertconv_mod_c
