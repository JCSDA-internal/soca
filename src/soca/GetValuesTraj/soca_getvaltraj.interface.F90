! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module handling interpolation trajectory

module soca_getvaltraj_mod_c

  use iso_c_binding
  use soca_getvaltraj_mod, only: soca_getvaltraj

  implicit none
  private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
!> Setup trajectory for interpolation
subroutine c_soca_getvaltraj_setup(c_self) &
  & bind(c,name='soca_getvaltraj_setup_f90')
  type(c_ptr), intent(inout) :: c_self

  type(soca_getvaltraj), pointer :: self

  allocate(self)
  c_self = c_loc(self)

  self%interph_initialized = .false.
  self%nobs = 0

end subroutine c_soca_getvaltraj_setup

! ------------------------------------------------------------------------------
!> Release memory
subroutine c_soca_getvaltraj_delete(c_self) &
  & bind(c,name='soca_getvaltraj_delete_f90')
  type(c_ptr), intent(inout) :: c_self

  type(soca_getvaltraj), pointer :: self

  call c_f_pointer(c_self, self)

  ! Clean up
  call self%horiz_interp%delete()
  self%nobs = 0
  self%interph_initialized = .false.
  deallocate(self)

end subroutine c_soca_getvaltraj_delete

end module soca_getvaltraj_mod_c
