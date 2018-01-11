! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module to handle variables for MOM5 & CICE5 model
module soca_vars_mod

  use iso_c_binding
  use config_mod

  implicit none
  private
  public :: soca_vars, soca_vars_create

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to represent MOM5 & CICE5 model variables
  type :: soca_vars
     integer                       :: nv          !< Number of variable type
     character(len=5), allocatable :: fldnames(:) !< Variable identifiers
  end type soca_vars

#define LISTED_TYPE soca_vars

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_vars_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

  ! ------------------------------------------------------------------------------

!!$  subroutine soca_vars_setup(self, cvars)
!!$    implicit none
!!$    type(soca_vars), intent(inout) :: self
!!$    character(len=5), intent(in) :: cvars(:)
!!$    !character(*), intent(in) :: cvars(:)
!!$    integer :: jj
!!$
!!$    self%nv = size(cvars)
!!$    do jj=1,self%nv
!!$       if (cvars(jj)/="cicen" .and. cvars(jj)/="hicen" .and. cvars(jj)/="vicen" &
!!$            .and. cvars(jj)/="hsnon" .and. cvars(jj)/="vsnon".and. cvars(jj)/="tsfcn" &
!!$            .and. cvars(jj)/="qsnon" .and. cvars(jj)/="sicnk".and. cvars(jj)/="socn" &
!!$            .and. cvars(jj)/="qicnk" .and. cvars(jj)/="tocn") then            
!!$          
!!$          call abor1_ftn ("soca_vars_setup: unknown field")
!!$       end if
!!$    enddo
!!$    allocate(self%fldnames(self%nv))
!!$    self%fldnames(:)=cvars(:)
!!$
!!$  end subroutine soca_vars_setup

  ! ------------------------------------------------------------------------------

!!$  subroutine c_soca_vars_create(c_key_self, c_conf) bind(c,name='soca_var_create_f90')
!!$    implicit none
!!$    integer(c_int), intent(inout) :: c_key_self
!!$    type(c_ptr), intent(in)    :: c_conf
!!$
!!$    type(soca_vars), pointer :: self
!!$    character(len=2) :: svar
!!$    integer :: nzi, nzo, ncat
!!$    
!!$    call soca_vars_registry%init()
!!$    call soca_vars_registry%add(c_key_self)
!!$    call soca_vars_registry%get(c_key_self, self)
!!$
!!$    svar = config_get_string(c_conf,len(svar),"variables")
!!$
!!$    select case (svar)
!!$    case ("nl","tl","cv","ci","x")
!!$       self%nv = 11
!!$       allocate(self%fldnames(self%nv))
!!$       self%fldnames(1) = "cicen"
!!$       self%fldnames(2) = "hicen"
!!$       self%fldnames(3) = "hsnon"
!!$       self%fldnames(4) = "tsfcn"
!!$       self%fldnames(5) = "qsnon"
!!$       self%fldnames(6) = "sicnk"
!!$       self%fldnames(7) = "qicnk"
!!$       self%fldnames(8) = "socn"
!!$       self%fldnames(9) = "tocn"
!!$       self%fldnames(10) = "ssh"
!!$    case default
!!$       call abor1_ftn("c_soca_vars_create: undefined variables")
!!$    end select
!!$
!!$    return
!!$  end subroutine c_soca_vars_create

  subroutine soca_vars_create(self, kvars)
    implicit none
    type(soca_vars), intent(inout) :: self
    integer(c_int), dimension(*), intent(in) :: kvars
    integer :: ii, jj

    if (kvars(1)<1 .or. kvars(1)>11) call abor1_ftn ("soca_vars_create: error variables")
    if (kvars(kvars(1)+2)/=999) call abor1_ftn ("soca_vars_create: error check")

    self%nv = 0

    do jj=1,kvars(1)
       ii=jj+1
       if (kvars(ii)<1 .or. kvars(ii)>11) call abor1_ftn ("soca_vars_create: unknown index")
       self%nv=self%nv+1
    enddo

    allocate(self%fldnames(self%nv))

    ii = 0
    do jj=1,kvars(1)
       if (kvars(jj+1)/=11) ii=ii+1
       if (kvars(jj+1)==1) self%fldnames(ii) = "cicen"
       if (kvars(jj+1)==2) self%fldnames(ii) = "hicen"
       if (kvars(jj+1)==3) self%fldnames(ii) = "hsnon"
       if (kvars(jj+1)==4) self%fldnames(ii) = "tsfcn"
       if (kvars(jj+1)==5) self%fldnames(ii) = "qsnon"
       if (kvars(jj+1)==6) self%fldnames(ii) = "sicnk"
       if (kvars(jj+1)==7) self%fldnames(ii) = "qicnk"
       if (kvars(jj+1)==8) self%fldnames(ii) = "socn"
       if (kvars(jj+1)==9) self%fldnames(ii) = "tocn"
       if (kvars(jj+1)==10) self%fldnames(ii) = "ssh"
    enddo
    
  end subroutine soca_vars_create
  
  ! ------------------------------------------------------------------------------

!!$  subroutine c_soca_vars_clone(c_key_self, c_key_other) bind(c,name='soca_var_clone_f90')
!!$    implicit none
!!$    integer(c_int), intent(in)    :: c_key_self
!!$    integer(c_int), intent(inout) :: c_key_other
!!$
!!$    type(soca_vars), pointer :: self, other
!!$
!!$    call soca_vars_registry%get(c_key_self, self)
!!$    call soca_vars_registry%add(c_key_other)
!!$    call soca_vars_registry%get(c_key_other, other)
!!$
!!$    call soca_vars_clone(self, other)
!!$
!!$  end subroutine c_soca_vars_clone

  ! ------------------------------------------------------------------------------

!!$  subroutine soca_vars_clone(self, other)
!!$    implicit none
!!$    type(soca_vars), intent(in)    :: self
!!$    type(soca_vars), intent(inout) :: other
!!$
!!$    other%nv = self%nv
!!$    allocate(other%fldnames(other%nv))
!!$    other%fldnames(:)=self%fldnames(:)
!!$
!!$  end subroutine soca_vars_clone
!!$
!!$  ! ------------------------------------------------------------------------------
!!$
!!$  subroutine c_soca_vars_delete(c_key_self) bind(c,name='soca_var_delete_f90')
!!$
!!$    implicit none
!!$    integer(c_int), intent(inout) :: c_key_self
!!$
!!$    type(soca_vars), pointer :: self
!!$    call soca_vars_registry%get(c_key_self, self)
!!$    deallocate(self%fldnames)
!!$    !deallocate(self%ns)
!!$    call soca_vars_registry%remove(c_key_self)
!!$
!!$    return
!!$  end subroutine c_soca_vars_delete

  ! ------------------------------------------------------------------------------

!!$  subroutine c_soca_vars_info(c_key_self, c_nv) bind(c,name='soca_var_info_f90')
!!$    implicit none
!!$    integer(c_int), intent(in)    :: c_key_self
!!$    integer(c_int), intent(inout) :: c_nv
!!$    type(soca_vars), pointer :: self
!!$
!!$    call soca_vars_registry%get(c_key_self, self)
!!$
!!$    c_nv = self%nv
!!$
!!$    return
!!$  end subroutine c_soca_vars_info

  ! ------------------------------------------------------------------------------

end module soca_vars_mod
