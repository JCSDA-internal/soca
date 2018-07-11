
!> Fortran module handling interpolation trajectory

module soca_getvaltraj_mod

  !General JEDI uses
  use kinds
  use iso_c_binding
  use type_bump, only: bump_type

  implicit none
  private

  public soca_getvaltraj
  public soca_getvaltraj_registry
  public c_soca_getvaltraj_setup, c_soca_getvaltraj_delete

  type :: soca_getvaltraj
     !integer :: nobs
     !type(bump_type) :: bump
     !logical :: lalloc = .false.
  end type soca_getvaltraj

#define LISTED_TYPE soca_getvaltraj

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_getvaltraj_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine c_soca_getvaltraj_setup(c_key_self) bind(c,name='soca_getvaltraj_setup_f90')

    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(soca_getvaltraj), pointer :: self

    ! Init, add and get key
    ! ---------------------
    call soca_getvaltraj_registry%init()
    call soca_getvaltraj_registry%add(c_key_self)
    call soca_getvaltraj_registry%get(c_key_self,self)

    print *,'======================================='
    print*, 'dh: getvaltraj_setup', c_key_self
    print *,'======================================='

    !self%lalloc = .false.
    !self%nobs = 0
    !self%ngrid = 0

  end subroutine c_soca_getvaltraj_setup

  ! ------------------------------------------------------------------------------

  subroutine c_soca_getvaltraj_delete(c_key_self) bind(c,name='soca_getvaltraj_delete_f90')

    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(soca_getvaltraj), pointer :: self

    ! Get key
    call soca_getvaltraj_registry%get(c_key_self, self)

    !if (self%lalloc) then
    !  self%nobs = 0
    !  self%ngrid = 0
    !  if (allocated(self%pt)) deallocate(self%pt)
    !  if (allocated(self%q)) deallocate(self%q)
    !  call self%bump%dealloc
    !  self%lalloc = .false.
    !endif

    ! Remove key
    call soca_getvaltraj_registry%remove(c_key_self)

  end subroutine c_soca_getvaltraj_delete

  ! ------------------------------------------------------------------------------

end module soca_getvaltraj_mod
