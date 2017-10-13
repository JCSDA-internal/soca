
!> Fortran module for fake observations for the MOM5CICE5 model
module mom5cice5_obsoper_mod

  use iso_c_binding
  use config_mod
  use mom5cice5_vars_mod
  use kinds

  implicit none
  private
  public :: mom5cice5_obsoper, mom5cice5_oper_setup, mom5cice5_oper_delete
  public :: mom5cice5_obsoper_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type for fake observations for the MOM5CICE5 model
  type :: mom5cice5_obsoper
     character(len=30) :: request
     type(mom5cice5_vars) :: varin
     integer :: ncol
  end type mom5cice5_obsoper

#define LISTED_TYPE mom5cice5_obsoper

  !> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: mom5cice5_obsoper_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "util/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_oper_setup(self, c_conf, svars, ncol)
    implicit none
    type(mom5cice5_obsoper), intent(inout) :: self
    type(c_ptr), intent(in)    :: c_conf
    character(len=*), intent(in) :: svars(:)
    integer :: ncol

    self%request = config_get_string(c_conf, len(self%request), "ObsType")
    call mom5cice5_vars_setup(self%varin, svars)
    self%ncol = ncol

  end subroutine mom5cice5_oper_setup

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_oper_delete(self)
    implicit none
    type(mom5cice5_obsoper), intent(inout) :: self

    deallocate(self%varin%fldnames)

  end subroutine mom5cice5_oper_delete

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_obsoper_inputs(c_key_self, c_key_vars) bind(c,name='mom5cice5_obsoper_inputs_f90')
    implicit none
    integer(c_int), intent(in)    :: c_key_self
    integer(c_int), intent(inout) :: c_key_vars

    type(mom5cice5_obsoper), pointer :: self
    type(mom5cice5_vars), pointer :: vars

    call mom5cice5_obsoper_registry%get(c_key_self, self)
    call mom5cice5_vars_registry%init()
    call mom5cice5_vars_registry%add(c_key_vars)
    call mom5cice5_vars_registry%get(c_key_vars, vars)
    call mom5cice5_vars_clone(self%varin, vars)

  end subroutine c_mom5cice5_obsoper_inputs

  ! ------------------------------------------------------------------------------

end module mom5cice5_obsoper_mod
