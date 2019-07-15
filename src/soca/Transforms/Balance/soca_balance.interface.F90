
!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

! ------------------------------------------------------------------------------
!> Constructor for D (standard deviation of background error)
subroutine c_soca_balance_setup(c_key_self, c_conf, c_key_traj) &
     &bind(c,name='soca_balance_setup_f90')
  use iso_c_binding
  use soca_balance_mod
  use soca_fields_mod_c

  integer(c_int), intent(inout) :: c_key_self   !< The D structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int), intent(in)    :: c_key_traj   !< Background field

  type(soca_field), pointer :: traj
  type(soca_balance_config), pointer :: self

  call soca_balance_registry%init()
  call soca_balance_registry%add(c_key_self)
  call soca_balance_registry%get(c_key_self, self)
  call soca_field_registry%get(c_key_traj, traj)

  call soca_balance_setup(c_conf, self, traj)
  
end subroutine c_soca_balance_setup

! ------------------------------------------------------------------------------
!> Destructor for D
subroutine c_soca_balance_delete(c_key_self) bind(c,name='soca_balance_delete_f90')
  use iso_c_binding
  use soca_balance_mod

  implicit none
  integer(c_int), intent(inout) :: c_key_self

  type(soca_balance_config), pointer :: self
  
  call soca_balance_registry%get(c_key_self,self)
  call soca_balance_delete(self)
  call soca_balance_registry%remove(c_key_self)

end subroutine c_soca_balance_delete

! ------------------------------------------------------------------------------
!> Multiplication forward
subroutine c_soca_balance_mult_f90(c_key_self, c_key_a, c_key_m)&
     &bind(c,name='soca_balance_mult_f90')
  use iso_c_binding
  use soca_balance_mod
  use soca_fields_mod_c
  use kinds
  use soca_kst_mod
  
  implicit none
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out 
  integer(c_int), intent(in) :: c_key_self 

  type(soca_field), pointer :: dxa
  type(soca_field), pointer :: dxm
  type(soca_balance_config), pointer :: self
  
  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_balance_registry%get(c_key_self,self)  

  !< Computes dxm = K dxa
  call soca_balance_mult(self, dxa, dxm)
  
end subroutine c_soca_balance_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication inverse
subroutine c_soca_balance_multinv_f90(c_key_self, c_key_m, c_key_a)&
     &bind(c,name='soca_balance_multinv_f90')
  use iso_c_binding
  use soca_balance_mod
  use soca_fields_mod_c
  use kinds
  use soca_kst_mod
  
  implicit none
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out 
  integer(c_int), intent(in) :: c_key_self 

  type(soca_field), pointer :: dxa
  type(soca_field), pointer :: dxm
  type(soca_balance_config), pointer :: self
  
  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_balance_registry%get(c_key_self,self)  

  !< Computes dxa = K^-1 dxm
  call soca_balance_multinv(self, dxa, dxm)
  
end subroutine c_soca_balance_multinv_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_balance_multad_f90(c_key_self, c_key_m, c_key_a)&
     &bind(c,name='soca_balance_multad_f90')
  use iso_c_binding
  use soca_balance_mod
  use soca_fields_mod_c
  use kinds
  use soca_kst_mod
  
  implicit none
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out 
  integer(c_int), intent(in) :: c_key_self 

  type(soca_field), pointer :: dxa
  type(soca_field), pointer :: dxm
  type(soca_balance_config), pointer :: self
  
  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_balance_registry%get(c_key_self,self)  

  !< Computes dxa = K^T dxm
  call soca_balance_multad(self, dxa, dxm)
  
end subroutine c_soca_balance_multad_f90

! ------------------------------------------------------------------------------
!> Multiplication inverse adjoint
subroutine c_soca_balance_multinvad_f90(c_key_self, c_key_a, c_key_m)&
     &bind(c,name='soca_balance_multinvad_f90')
  use iso_c_binding
  use soca_balance_mod
  use soca_fields_mod_c
  use kinds
  use soca_kst_mod
  
  implicit none
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out 
  integer(c_int), intent(in) :: c_key_self 

  type(soca_field), pointer :: dxa
  type(soca_field), pointer :: dxm
  type(soca_balance_config), pointer :: self
  
  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_balance_registry%get(c_key_self,self)  

  !< Computes dxm = (K^-1)^T dxa
  call soca_balance_multinvad(self, dxa, dxm)
  
end subroutine c_soca_balance_multinvad_f90
