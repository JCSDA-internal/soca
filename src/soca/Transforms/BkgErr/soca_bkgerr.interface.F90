! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerr_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use soca_fields_mod, only: soca_fields
use soca_bkgerr_mod, only: soca_bkgerr_config, &
                           soca_bkgerr_setup, soca_bkgerr_mult

implicit none

private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Constructor for D (standard deviation of background error)
subroutine c_soca_bkgerr_setup(c_self, c_conf, c_bkg) &
  & bind(c,name='soca_bkgerr_setup_f90')
  type(c_ptr), intent(inout) :: c_self   !< The D structure
  type(c_ptr),    intent(in) :: c_conf       !< The configuration
  type(c_ptr),    intent(in) :: c_bkg    !< Background field

  type(soca_fields), pointer :: bkg
  type(soca_bkgerr_config), pointer :: self

  allocate(self)
  c_self = c_loc(self)
  call c_f_pointer(c_bkg, bkg)

  call soca_bkgerr_setup(fckit_configuration(c_conf), self, bkg)
end subroutine c_soca_bkgerr_setup

! ------------------------------------------------------------------------------
!> Destructor for D
subroutine c_soca_bkgerr_delete(c_self) &
  & bind(c,name='soca_bkgerr_delete_f90')
  type(c_ptr), intent(inout) :: c_self

  type(soca_bkgerr_config), pointer :: self

  call c_f_pointer(c_self, self)

  if (associated(self%bkg)) nullify(self%bkg)
  call self%std_bkgerr%delete()
  deallocate(self)

end subroutine c_soca_bkgerr_delete

! ------------------------------------------------------------------------------
!> Multiplication forward and adjoint
subroutine c_soca_bkgerr_mult_f90(c_self, c_a, c_m) &
  & bind(c,name='soca_bkgerr_mult_f90')
  type(c_ptr),    intent(in) :: c_self
  type(c_ptr),    intent(in) :: c_a     !<    "   to Increment in
  type(c_ptr), intent(inout) :: c_m     !<    "   to Increment out

  type(soca_fields),        pointer :: dxa, dxm
  type(soca_bkgerr_config), pointer :: self

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_a, dxa)
  call c_f_pointer(c_m, dxm)

  !< Computes dxm = D dxa
  call dxm%copy(dxa)
  call soca_bkgerr_mult(self, dxa, dxm)

end subroutine c_soca_bkgerr_mult_f90

end module soca_bkgerr_mod_c
