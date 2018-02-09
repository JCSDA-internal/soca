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
       if (kvars(jj+1)/=10) ii=ii+1
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

end module soca_vars_mod
