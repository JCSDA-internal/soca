! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_balance_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use soca_fields_mod, only: soca_fields
use soca_balance_mod, only: soca_balance_config, &
                            soca_balance_setup, soca_balance_delete, &
                            soca_balance_mult, soca_balance_multad, &
                            soca_balance_multinv, soca_balance_multinvad

implicit none
private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Constructor for D (standard deviation of background error)
subroutine c_soca_balance_setup(c_self, c_conf, c_traj) &
  & bind(c,name='soca_balance_setup_f90')
  type(c_ptr), intent(inout) :: c_self
  type(c_ptr),    intent(in) :: c_conf       !< The configuration
  type(c_ptr),    intent(in) :: c_traj   !< Background field

  type(soca_fields), pointer :: traj
  type(soca_balance_config), pointer :: self

  ! create the fortran version of ourself
  allocate(self)
  c_self = c_loc(self)
  call c_f_pointer(c_traj, traj)

  call soca_balance_setup(fckit_configuration(c_conf), self, traj)

end subroutine c_soca_balance_setup

! ------------------------------------------------------------------------------
!> Destructor for D
subroutine c_soca_balance_delete(c_self) &
  & bind(c,name='soca_balance_delete_f90')
  type(c_ptr), intent(inout) :: c_self

  type(soca_balance_config), pointer :: self

  call c_f_pointer(c_self, self)
  call soca_balance_delete(self)
  deallocate(self)

end subroutine c_soca_balance_delete

! ------------------------------------------------------------------------------
!> Multiplication forward
subroutine c_soca_balance_mult_f90(c_self, c_a, c_m)&
  & bind(c,name='soca_balance_mult_f90')
  type(c_ptr),    intent(in) :: c_self
  type(c_ptr),    intent(in) :: c_a     !<    "   to Increment in
  type(c_ptr), intent(inout) :: c_m     !<    "   to Increment out

  type(soca_fields),         pointer :: dxa, dxm
  type(soca_balance_config), pointer :: self

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_a, dxa)
  call c_f_pointer(c_m, dxm)

  !< Computes dxm = K dxa
  call soca_balance_mult(self, dxa, dxm)

end subroutine c_soca_balance_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication inverse
subroutine c_soca_balance_multinv_f90(c_self, c_m, c_a)&
  & bind(c,name='soca_balance_multinv_f90')
  type(c_ptr),    intent(in) :: c_self
  type(c_ptr),    intent(in) :: c_m     !<    "   to Increment out
  type(c_ptr), intent(inout) :: c_a     !<    "   to Increment in

  type(soca_fields),         pointer :: dxa, dxm
  type(soca_balance_config), pointer :: self

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_a, dxa)
  call c_f_pointer(c_m, dxm)

  !< Computes dxa = K^-1 dxm
  call soca_balance_multinv(self, dxa, dxm)

end subroutine c_soca_balance_multinv_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_balance_multad_f90(c_self, c_m, c_a)&
  & bind(c,name='soca_balance_multad_f90')
  type(c_ptr),    intent(in) :: c_self
  type(c_ptr),    intent(in) :: c_m     !<    "   to Increment out
  type(c_ptr), intent(inout) :: c_a     !<    "   to Increment in

  type(soca_fields),         pointer :: dxa, dxm
  type(soca_balance_config), pointer :: self

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_a, dxa)
  call c_f_pointer(c_m, dxm)

  !< Computes dxa = K^T dxm
  call soca_balance_multad(self, dxa, dxm)

end subroutine c_soca_balance_multad_f90

! ------------------------------------------------------------------------------
!> Multiplication inverse adjoint
subroutine c_soca_balance_multinvad_f90(c_self, c_a, c_m)&
  & bind(c,name='soca_balance_multinvad_f90')
  type(c_ptr),    intent(in) :: c_self
  type(c_ptr),    intent(in) :: c_a     !<    "   to Increment in
  type(c_ptr), intent(inout) :: c_m     !<    "   to Increment out

  type(soca_fields),         pointer :: dxa, dxm
  type(soca_balance_config), pointer :: self

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_a, dxa)
  call c_f_pointer(c_m, dxm)

  !< Computes dxm = (K^-1)^T dxa
  call soca_balance_multinvad(self, dxa, dxm)

end subroutine c_soca_balance_multinvad_f90

end module soca_balance_mod_c
