
!> Fortran module handling observation locations

module mom5cice5_locs_mod

  use iso_c_binding
  use mom5cice5_obs_vectors
  use kinds

  implicit none
  private
  public :: mom5cice5_locs, mom5cice5_loc_setup
  public :: mom5cice5_locs_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to hold observation locations
  type :: mom5cice5_locs
     integer                           :: nloc     !< Number of obs loc in ]t,t+dt] (see ObsSpace) CONFUSING ... CHECK
     real(kind=kind_real), allocatable :: xyz(:,:) !< Need to be allocated as (3, nloc) 3: lon, lat, lev
  end type mom5cice5_locs

#define LISTED_TYPE mom5cice5_locs

  !> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: mom5cice5_locs_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "util/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_loc_setup(self, lvec)
    implicit none
    type(mom5cice5_locs), intent(inout) :: self
    type(obs_vect), intent(in) :: lvec
    integer :: jc, jo

    self%nloc=lvec%nobs
    allocate(self%xyz(3,self%nloc))
    do jo=1,self%nloc
       do jc=1,3
          self%xyz(jc,jo)=lvec%values(jc,jo)
       enddo
    enddo

  end subroutine mom5cice5_loc_setup

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_loc_delete(key) bind(c,name='mom5cice5_loc_delete_f90')

    implicit none
    integer(c_int), intent(inout) :: key
    type(mom5cice5_locs), pointer :: self

    call mom5cice5_locs_registry%get(key,self)
    deallocate(self%xyz)
    call mom5cice5_locs_registry%remove(key)

  end subroutine c_mom5cice5_loc_delete

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_loc_nobs(key, kobs) bind(c,name='mom5cice5_loc_nobs_f90')

    implicit none
    integer(c_int), intent(in) :: key
    integer(c_int), intent(inout) :: kobs
    type(mom5cice5_locs), pointer :: self

    call mom5cice5_locs_registry%get(key,self)
    kobs = self%nloc

  end subroutine c_mom5cice5_loc_nobs

  ! ------------------------------------------------------------------------------

end module mom5cice5_locs_mod
