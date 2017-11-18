
!> Fortran module handling observation locations

module soca_locs_mod

  use iso_c_binding
  use soca_obs_vectors
  use kinds

  implicit none
  private
  public :: soca_locs, soca_loc_setup
  public :: soca_locs_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to hold observation locations
  type :: soca_locs
     integer                           :: nloc     !< Number of obs loc in ]t,t+dt] (see ObsSpace) CONFUSING ... CHECK
     real(kind=kind_real), allocatable :: xyz(:,:) !< Need to be allocated as (3, nloc) 3: lon, lat, lev
  end type soca_locs

#define LISTED_TYPE soca_locs

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_locs_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine soca_loc_setup(self, lvec)
    implicit none
    type(soca_locs), intent(inout) :: self
    type(obs_vect), intent(in) :: lvec
    integer :: jc, jo

    self%nloc=lvec%nobs
    allocate(self%xyz(3,self%nloc))
    do jo=1,self%nloc
       do jc=1,3
          self%xyz(jc,jo)=lvec%values(jc,jo)
       enddo
    enddo

  end subroutine soca_loc_setup

  ! ------------------------------------------------------------------------------

  subroutine c_soca_loc_delete(key) bind(c,name='soca_loc_delete_f90')

    implicit none
    integer(c_int), intent(inout) :: key
    type(soca_locs), pointer :: self

    call soca_locs_registry%get(key,self)
    deallocate(self%xyz)
    call soca_locs_registry%remove(key)

  end subroutine c_soca_loc_delete

  ! ------------------------------------------------------------------------------

  subroutine c_soca_loc_nobs(key, kobs) bind(c,name='soca_loc_nobs_f90')

    implicit none
    integer(c_int), intent(in) :: key
    integer(c_int), intent(inout) :: kobs
    type(soca_locs), pointer :: self

    call soca_locs_registry%get(key,self)
    kobs = self%nloc

  end subroutine c_soca_loc_nobs

  ! ------------------------------------------------------------------------------

end module soca_locs_mod
