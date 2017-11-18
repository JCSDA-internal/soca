
!> Fortran module for observation operators for the model
module soca_obsoper_mod

  use iso_c_binding
  use config_mod
  use soca_vars_mod
  use kinds

  implicit none
  private
  public :: soca_obsoper, soca_oper_setup, soca_oper_delete
  public :: soca_obsoper_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type for observations for the model
  type :: soca_obsoper
     character(len=30)    :: request
     type(soca_vars) :: varin
     integer              :: ncol
  end type soca_obsoper

#define LISTED_TYPE soca_obsoper

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_obsoper_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine soca_oper_setup(self, c_conf, svars, ncol)
    implicit none
    type(soca_obsoper), intent(inout) :: self
    type(c_ptr), intent(in)    :: c_conf
    character(len=*), intent(in) :: svars(:)
    integer :: ncol

    self%request = config_get_string(c_conf, len(self%request), "ObsType")
    call soca_vars_setup(self%varin, svars)
    self%ncol = ncol

  end subroutine soca_oper_setup

  ! ------------------------------------------------------------------------------

  subroutine soca_oper_delete(self)
    implicit none
    type(soca_obsoper), intent(inout) :: self

    deallocate(self%varin%fldnames)

  end subroutine soca_oper_delete

  ! ------------------------------------------------------------------------------

  subroutine c_soca_obsoper_inputs(c_key_self, c_key_vars) bind(c,name='soca_obsoper_inputs_f90')
    implicit none
    integer(c_int), intent(in)    :: c_key_self
    integer(c_int), intent(inout) :: c_key_vars

    type(soca_obsoper), pointer :: self
    type(soca_vars), pointer :: vars

    call soca_obsoper_registry%get(c_key_self, self)
    call soca_vars_registry%init()
    call soca_vars_registry%add(c_key_vars)
    call soca_vars_registry%get(c_key_vars, vars)
    call soca_vars_clone(self%varin, vars)

  end subroutine c_soca_obsoper_inputs

  ! ------------------------------------------------------------------------------

end module soca_obsoper_mod
