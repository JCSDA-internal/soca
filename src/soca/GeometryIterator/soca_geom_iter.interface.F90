! (C) Copyright 2019-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for soca_geom_iter_mod::soca_geom_iter
module soca_geom_iter_mod_c

use iso_c_binding
use kinds
use soca_geom_iter_mod
use soca_geom_mod_c,  only : soca_geom_registry
use soca_geom_mod, only: soca_geom

implicit none
private

#define LISTED_TYPE soca_geom_iter

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for soca_geom_iter instances
type(registry_t), public :: soca_geom_iter_registry


contains


!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
!> C++ interface for soca_geom_iter_mod::soca_geom_iter::setup()
subroutine soca_geom_iter_setup_c(c_key_self, c_key_geom, c_iindex, c_jindex) bind(c, name='soca_geom_iter_setup_f90')
  integer(c_int), intent(inout) :: c_key_self !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_geom !< Geometry
  integer(c_int), intent(   in) :: c_iindex    !< Index
  integer(c_int), intent(   in) :: c_jindex    !< Index

  ! Local variables
  type(soca_geom_iter),     pointer :: self
  type(soca_geom),          pointer :: geom

  ! Interface
  call soca_geom_iter_registry%init()
  call soca_geom_iter_registry%add(c_key_self)
  call soca_geom_iter_registry%get(c_key_self, self)
  call soca_geom_registry%get(c_key_geom, geom)

  ! Call Fortran
  call self%setup(geom, c_iindex, c_jindex)

end subroutine soca_geom_iter_setup_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_geom_iter_mod::soca_geom_iter::clone()
subroutine soca_geom_iter_clone_c(c_key_self, c_key_other) bind(c, name='soca_geom_iter_clone_f90')
  integer(c_int), intent(inout) :: c_key_self  !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_other !< Other geometry iterator

  ! Local variables
  type(soca_geom_iter), pointer :: self, other

  ! Interface
  call soca_geom_iter_registry%get(c_key_other, other)
  call soca_geom_iter_registry%init()
  call soca_geom_iter_registry%add(c_key_self)
  call soca_geom_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%clone(other)

end subroutine soca_geom_iter_clone_c


! ------------------------------------------------------------------------------
!> !> C++ interface for deleting soca_geom_iter_mod::soca_geom_iter
subroutine soca_geom_iter_delete_c(c_key_self) bind(c, name='soca_geom_iter_delete_f90')
    integer(c_int), intent(inout) :: c_key_self !< Geometry iterator

    ! Clear interface
    call soca_geom_iter_registry%remove(c_key_self)

end subroutine soca_geom_iter_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_geom_iter_mod::soca_geom_iter::equals()
subroutine soca_geom_iter_equals_c(c_key_self, c_key_other, c_equals) bind(c, name='soca_geom_iter_equals_f90')
  integer(c_int), intent(inout) :: c_key_self  !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_other !< Other geometry iterator
  integer(c_int), intent(inout) :: c_equals    !< Equality flag

  ! Local variables
  type(soca_geom_iter),pointer :: self,other

  ! Interface
  call soca_geom_iter_registry%get(c_key_self, self)
  call soca_geom_iter_registry%get(c_key_other, other)

  ! Call Fortran
  call self%equals(other, c_equals)

end subroutine soca_geom_iter_equals_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_geom_iter_mod::soca_geom_iter::current()
subroutine soca_geom_iter_current_c(c_key_self, c_lon, c_lat) bind(c, name='soca_geom_iter_current_f90')
  integer(c_int), intent(   in) :: c_key_self !< Geometry iterator
  real(c_double), intent(inout) :: c_lat      !< Latitude
  real(c_double), intent(inout) :: c_lon      !< Longitude

  ! Local variables
  type(soca_geom_iter), pointer :: self

  ! Interface
  call soca_geom_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%current(c_lon, c_lat)

end subroutine soca_geom_iter_current_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_geom_iter_mod::soca_geom_iter::next()
subroutine soca_geom_iter_next_c(c_key_self) bind(c, name='soca_geom_iter_next_f90')
  integer(c_int), intent(in) :: c_key_self !< Geometry iterator

  ! Local variables
  type(soca_geom_iter), pointer :: self

  ! Interface
  call soca_geom_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%next()
end subroutine soca_geom_iter_next_c


! ------------------------------------------------------------------------------
!> C++ interface to get cell area from soca_geom_iter_mod::soca_geom_iter
subroutine soca_geom_iter_get_area_c(c_key_self, c_val) bind(c, name='soca_geom_iter_get_area_f90')
  integer(c_int), intent(   in) :: c_key_self !< Geometry iterator
  real(c_double), intent(inout) :: c_val

  type(soca_geom_iter), pointer :: self
  call soca_geom_iter_registry%get(c_key_self, self)

  c_val = self%geom%cell_area(self%iind,self%jind)
end subroutine


! ------------------------------------------------------------------------------
!> C++ interface to get rossby radius from soca_geom_iter_mod::soca_geom_iter
subroutine soca_geom_iter_get_rossby_c(c_key_self, c_val) bind(c, name='soca_geom_iter_get_rossby_f90')
  integer(c_int), intent(   in) :: c_key_self !< Geometry iterator
  real(c_double), intent(inout) :: c_val

  type(soca_geom_iter), pointer :: self
  call soca_geom_iter_registry%get(c_key_self, self)

  c_val = self%geom%rossby_radius(self%iind,self%jind)
end subroutine


end module soca_geom_iter_mod_c
