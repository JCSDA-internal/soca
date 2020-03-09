! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_balance_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use soca_fields_mod, only: soca_fields
use soca_fields_mod_c, only: soca_field_registry
use soca_balance_mod, only: soca_balance_config, &
                            soca_balance_setup, soca_balance_delete, &
                            soca_balance_mult, soca_balance_multad, &
                            soca_balance_multinv, soca_balance_multinvad

implicit none

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Constructor for D (standard deviation of background error)
subroutine c_soca_balance_setup(c_self, c_conf, c_key_traj) &
  bind(c,name='soca_balance_setup_f90')

  type(c_ptr),   intent(inout) :: c_self
  type(c_ptr),      intent(in) :: c_conf       !< The configuration
  integer(c_int), intent(in)   :: c_key_traj   !< Background field

  type(soca_fields), pointer :: traj
  type(soca_balance_config), pointer :: self

  ! create the fortran version of ourself
  allocate(self)
  c_self = c_loc(self)

  call soca_field_registry%get(c_key_traj, traj)

  call soca_balance_setup(fckit_configuration(c_conf), self, traj)

end subroutine c_soca_balance_setup

! ------------------------------------------------------------------------------
!> Destructor for D
subroutine c_soca_balance_delete(c_self) bind(c,name='soca_balance_delete_f90')

  type(c_ptr), intent(in) :: c_self

  type(soca_balance_config), pointer :: self

  call c_f_pointer(c_self, self)
  call soca_balance_delete(self)
  deallocate(self)

end subroutine c_soca_balance_delete

! ------------------------------------------------------------------------------
!> Multiplication forward
subroutine c_soca_balance_mult_f90(c_self, c_key_a, c_key_m)&
  bind(c,name='soca_balance_mult_f90')

  type(c_ptr),    intent(in) :: c_self
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out

  type(soca_fields), pointer :: dxa
  type(soca_fields), pointer :: dxm
  type(soca_balance_config), pointer :: self

  call c_f_pointer(c_self, self)
  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)

  !< Computes dxm = K dxa
  call soca_balance_mult(self, dxa, dxm)

end subroutine c_soca_balance_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication inverse
subroutine c_soca_balance_multinv_f90(c_self, c_key_m, c_key_a)&
  bind(c,name='soca_balance_multinv_f90')

  type(c_ptr),    intent(in) :: c_self
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out


  type(soca_fields), pointer :: dxa
  type(soca_fields), pointer :: dxm
  type(soca_balance_config), pointer :: self

  call c_f_pointer(c_self, self)
  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)

  !< Computes dxa = K^-1 dxm
  call soca_balance_multinv(self, dxa, dxm)

end subroutine c_soca_balance_multinv_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_balance_multad_f90(c_self, c_key_m, c_key_a)&
  bind(c,name='soca_balance_multad_f90')

  type(c_ptr),    intent(in) :: c_self
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out

  type(soca_fields), pointer :: dxa
  type(soca_fields), pointer :: dxm
  type(soca_balance_config), pointer :: self

  call c_f_pointer(c_self, self)
  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)

  !< Computes dxa = K^T dxm
  call soca_balance_multad(self, dxa, dxm)

end subroutine c_soca_balance_multad_f90

! ------------------------------------------------------------------------------
!> Multiplication inverse adjoint
subroutine c_soca_balance_multinvad_f90(c_self, c_key_a, c_key_m)&
  bind(c,name='soca_balance_multinvad_f90')

  type(c_ptr),    intent(in) :: c_self
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out

  type(soca_fields), pointer :: dxa
  type(soca_fields), pointer :: dxm
  type(soca_balance_config), pointer :: self

  call c_f_pointer(c_self, self)
  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)

  !< Computes dxm = (K^-1)^T dxa
  call soca_balance_multinvad(self, dxa, dxm)

end subroutine c_soca_balance_multinvad_f90

end module soca_balance_mod_c
