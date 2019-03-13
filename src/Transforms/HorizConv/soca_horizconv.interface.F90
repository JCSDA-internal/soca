!
! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

! ------------------------------------------------------------------------------
!> Constructor for Horizconv
subroutine c_soca_horizconv_setup(c_key_self, c_conf, c_key_bkg) &
     & bind(c,name='soca_horizconv_setup_f90')
  use iso_c_binding
  use soca_horizconv_mod
  use config_mod
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices
  
  implicit none
  
  integer(c_int), intent(inout) :: c_key_self   !< The Horizconv structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  !integer(c_int),    intent(in) :: c_key_traj   !< trajectory
  integer(c_int),    intent(in) :: c_key_bkg    !< background  

  type(soca_horizconv), pointer :: self
  !type(soca_field), pointer :: traj
  type(soca_field), pointer :: bkg  
  real(kind=kind_real), allocatable :: jac(:)

  integer :: isc, iec, jsc, jec, i, j, k, nl
  
  call soca_horizconv_registry%init()
  call soca_horizconv_registry%add(c_key_self)
  call soca_horizconv_registry%get(c_key_self, self)
  !call soca_field_registry%get(c_key_traj, traj)
  call soca_field_registry%get(c_key_bkg, bkg)  

  call soca_conv_setup (self, bkg, c_conf)

end subroutine c_soca_horizconv_setup

! ------------------------------------------------------------------------------
!> Destructor for Horizconv
subroutine c_soca_horizconv_delete(c_key_self) bind(c,name='soca_horizconv_delete_f90')
  use iso_c_binding
  use soca_horizconv_mod

  implicit none

  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure

  type(soca_horizconv), pointer :: self

  ! Deallocate trajectory and backgroun
  ! TODO
  ! Deallocate ocean depth array
  ! TODO
  
  call soca_horizconv_registry%get(c_key_self, self)  
  
  !if (associated(self%traj)) nullify(self%traj)
  if (associated(self%bkg)) nullify(self%bkg)  
  
  call soca_horizconv_registry%remove(c_key_self)
  
end subroutine c_soca_horizconv_delete

! ------------------------------------------------------------------------------
!> Multiplication
subroutine c_soca_horizconv_mult_f90(c_key_a, c_key_m, c_key_self)&
     & bind(c,name='soca_horizconv_mult_f90')
  use iso_c_binding
  use soca_horizconv_mod
  use soca_fields
  use soca_utils
  use kinds
  use config_mod
  
  implicit none

  integer(c_int), intent(in) :: c_key_a     !< Increment in
  integer(c_int), intent(in) :: c_key_m     !< Increment out 
  !integer(c_int), intent(in) :: c_key_traj  !< trajectory
  integer(c_int), intent(in) :: c_key_self  !< config
  
  type(soca_field), pointer :: dxa  ! in
  type(soca_field), pointer :: dxm  ! out
  !type(soca_field), pointer :: traj  
  type(soca_horizconv),   pointer :: self
  
  call soca_field_registry%get(c_key_a, dxa)
  call soca_field_registry%get(c_key_m, dxm)
  !call soca_field_registry%get(c_key_traj, traj)
  call soca_horizconv_registry%get(c_key_self, self)    

  !< Computes dxm = Horizconv dxa

  ! dxm = dxa
  call copy(dxm, dxa) 

  ! Apply forward convolution operator to T & S
  call soca_conv(self, dxm, dxa)

end subroutine c_soca_horizconv_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_horizconv_multad_f90(c_key_m, c_key_a, c_key_self)&
     & bind(c,name='soca_horizconv_multad_f90')
  use iso_c_binding
  use soca_horizconv_mod
  use soca_fields
  use kinds
  
  implicit none

  integer(c_int), intent(in) :: c_key_a     !< Increment out
  integer(c_int), intent(in) :: c_key_m     !< Increment in 
  !integer(c_int), intent(in) :: c_key_traj  !< Trajectory
  integer(c_int), intent(in) :: c_key_self  !< config

  type(soca_field),      pointer :: dxa
  type(soca_field),      pointer :: dxm
  !type(soca_field),      pointer :: traj
  type(soca_horizconv),   pointer :: self

  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  !call soca_field_registry%get(c_key_traj,traj)
  call soca_horizconv_registry%get(c_key_self, self)    

  ! dxa = dxm
  call copy(dxa,dxm)
  
  ! Apply adjoint of convolution operator
  call soca_conv_ad(self, dxm, dxa)

end subroutine c_soca_horizconv_multad_f90
