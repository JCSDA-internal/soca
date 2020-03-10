! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerrgodas_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use soca_fields_mod, only: soca_fields
use soca_bkgerrgodas_mod, only: soca_bkgerrgodas_config, &
                                soca_bkgerrgodas_setup, soca_bkgerrgodas_mult

implicit none
private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Constructor for D (standard deviation of background error)
subroutine c_soca_bkgerrgodas_setup(c_self, c_conf, c_bkg) &
  & bind(c,name='soca_bkgerrgodas_setup_f90')
  type(c_ptr),  intent(inout) :: c_self   !< The D structure
  type(c_ptr),     intent(in) :: c_conf       !< The configuration
  type(c_ptr),     intent(in) :: c_bkg    !< Background field

  type(soca_bkgerrgodas_config), pointer :: self
  type(soca_fields),             pointer :: bkg

  allocate(self)
  c_self = c_loc(self)
  call c_f_pointer(c_bkg, bkg)

  call soca_bkgerrgodas_setup(fckit_configuration(c_conf), self, bkg)

end subroutine c_soca_bkgerrgodas_setup

! ------------------------------------------------------------------------------
!> Destructor for D
subroutine c_soca_bkgerrgodas_delete(c_self) &
  & bind(c,name='soca_bkgerrgodas_delete_f90')
  type(c_ptr), intent(inout) :: c_self

  type(soca_bkgerrgodas_config), pointer :: self

  call c_f_pointer(c_self, self)

  if (associated(self%bkg)) nullify(self%bkg)
  call self%std_bkgerr%delete()
  deallocate(self)

end subroutine c_soca_bkgerrgodas_delete

! ------------------------------------------------------------------------------
!> Multiplication forward and adjoint
subroutine c_soca_bkgerrgodas_mult_f90(c_self, c_a, c_m) &
  & bind(c,name='soca_bkgerrgodas_mult_f90')
  type(c_ptr),    intent(in) :: c_self
  type(c_ptr),    intent(in) :: c_a     !<    "   to Increment in
  type(c_ptr), intent(inout) :: c_m     !<    "   to Increment out

  type(soca_bkgerrgodas_config), pointer :: self
  type(soca_fields),             pointer :: dxa, dxm

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_a, dxa)
  call c_f_pointer(c_m, dxm)

  !< Computes dxm = D dxa
  call dxm%copy(dxa)
  call soca_bkgerrgodas_mult(self, dxa, dxm)

end subroutine c_soca_bkgerrgodas_mult_f90

end module soca_bkgerrgodas_mod_c
