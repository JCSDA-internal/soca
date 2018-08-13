!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

! ------------------------------------------------------------------------------
!> Constructor for Kst
subroutine c_soca_kst_setup(c_key_self, c_conf) bind(c,name='soca_kst_setup_f90')
  use iso_c_binding
  use soca_kst_mod

  integer(c_int), intent(inout) :: c_key_self   !< The Kst structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  
  type(soca_kst_config), pointer :: self

  call soca_kst_registry%init()
  call soca_kst_registry%add(c_key_self)
  call soca_kst_registry%get(c_key_self, self)

  call soca_kst_setup(c_conf, self)
  
end subroutine c_soca_kst_setup

! ------------------------------------------------------------------------------
!> Destructor for Kst
subroutine c_soca_kst_delete(c_key_self) bind(c,name='soca_kst_delete_f90')
  use iso_c_binding
  use soca_kst_mod

  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure
  type(soca_kst_config), pointer :: self

  print *,'key:',c_key_self
  call soca_kst_registry%get(c_key_self,self)
  call soca_kst_registry%remove(c_key_self)
  
end subroutine c_soca_kst_delete

! ------------------------------------------------------------------------------
!> Linearization setup for Kst
subroutine c_soca_kst_linearize_f90() bind(c,name='soca_kst_linearize_f90')
  use iso_c_binding
end subroutine c_soca_kst_linearize_f90

! ------------------------------------------------------------------------------
!> Multiplication
subroutine c_soca_kst_mult_f90() bind(c,name='soca_kst_mult_f90')
  use iso_c_binding  
end subroutine c_soca_kst_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_kst_multad_f90() bind(c,name='soca_kst_multad_f90')
  use iso_c_binding  
end subroutine c_soca_kst_multad_f90
