
!> Fortran module to handle variables for MOM5 & CICE5 model
module mom5cice5_vars_mod

  use iso_c_binding
  use config_mod

  implicit none
  private
  public :: mom5cice5_vars, mom5cice5_vars_setup, mom5cice5_vars_clone
  public :: mom5cice5_vars_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to represent MOM5 & CICE5 model variables
  type :: mom5cice5_vars
     integer :: nv
     character(len=5), allocatable :: fldnames(:) !< Variable identifiers
  end type mom5cice5_vars

#define LISTED_TYPE mom5cice5_vars

  !> Linked list interface - defines registry_t type
#include "linkedList_i.f"

  !> Global registry
  type(registry_t) :: mom5cice5_vars_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_vars_setup(self, cvars)
    implicit none
    type(mom5cice5_vars), intent(inout) :: self
    character(len=5), intent(in) :: cvars(:)
    integer :: jj

    self%nv = size(cvars)

    do jj=1,self%nv
       if (cvars(jj)/="cicen" .and. cvars(jj)/="hicen" .and. cvars(jj)/="vicen" &
            .and. cvars(jj)/="hsnon" .and. cvars(jj)/="vsnon".and. cvars(jj)/="tsfcn" &
            .and. cvars(jj)/="qsnon" .and. cvars(jj)/="sicnk".and. cvars(jj)/="sssoc" &
            .and. cvars(jj)/="qicnk" .and. cvars(jj)/="tlioc".and. cvars(jj)/="sstoc") then        

          call abor1_ftn ("mom5cice5_vars_setup: unknown field")
       end if
    enddo
    allocate(self%fldnames(self%nv))
    self%fldnames(:)=cvars(:)

  end subroutine mom5cice5_vars_setup

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_vars_create(c_key_self, c_conf) bind(c,name='mom5cice5_var_create_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr), intent(in)    :: c_conf

    type(mom5cice5_vars), pointer :: self
    character(len=2) :: svar

    call mom5cice5_vars_registry%init()
    call mom5cice5_vars_registry%add(c_key_self)
    call mom5cice5_vars_registry%get(c_key_self, self)

    svar = config_get_string(c_conf,len(svar),"variables")
    print *,'svar=',svar
    select case (svar)
    case ("ci","cv")
       self%nv = 12
       allocate(self%fldnames(self%nv))
       self%fldnames(1) = "cicen"
       self%fldnames(2) = "hicen"
       self%fldnames(3) = "vicen"
       self%fldnames(4) = "hsnon"
       self%fldnames(5) = "vsnon"
       self%fldnames(6) = "tsfcn"
       self%fldnames(7) = "qsnon"
       self%fldnames(8) = "sicnk"
       self%fldnames(9) = "sssoc"
       self%fldnames(10) = "qicnk"
       self%fldnames(11) = "tlioc"
       self%fldnames(12) = "sstoc"
    case default
       call abor1_ftn("c_mom5cice5_vars_create: undefined variables")
    end select

    return
  end subroutine c_mom5cice5_vars_create

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_vars_clone(c_key_self, c_key_other) bind(c,name='mom5cice5_var_clone_f90')
    implicit none
    integer(c_int), intent(in)    :: c_key_self
    integer(c_int), intent(inout) :: c_key_other

    type(mom5cice5_vars), pointer :: self, other

    call mom5cice5_vars_registry%get(c_key_self, self)
    call mom5cice5_vars_registry%add(c_key_other)
    call mom5cice5_vars_registry%get(c_key_other, other)

    call mom5cice5_vars_clone(self, other)

  end subroutine c_mom5cice5_vars_clone

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_vars_clone(self, other)
    implicit none
    type(mom5cice5_vars), intent(in)    :: self
    type(mom5cice5_vars), intent(inout) :: other

    other%nv = self%nv
    allocate(other%fldnames(other%nv))
    other%fldnames(:)=self%fldnames(:)

  end subroutine mom5cice5_vars_clone

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_vars_delete(c_key_self) bind(c,name='mom5cice5_var_delete_f90')

    implicit none
    integer(c_int), intent(inout) :: c_key_self

    type(mom5cice5_vars), pointer :: self
    call mom5cice5_vars_registry%get(c_key_self, self)
    deallocate(self%fldnames)
    call mom5cice5_vars_registry%remove(c_key_self)

    return
  end subroutine c_mom5cice5_vars_delete

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_vars_info(c_key_self, c_nv) bind(c,name='mom5cice5_var_info_f90')
    implicit none
    integer(c_int), intent(in)    :: c_key_self
    integer(c_int), intent(inout) :: c_nv
    type(mom5cice5_vars), pointer :: self

    call mom5cice5_vars_registry%get(c_key_self, self)

    c_nv = self%nv

    return
  end subroutine c_mom5cice5_vars_info

  ! ------------------------------------------------------------------------------

end module mom5cice5_vars_mod
