!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

! ------------------------------------------------------------------------------
!> Constructor for Vertconv
subroutine c_soca_vertconv_setup(c_key_self, c_conf, c_key_traj, c_key_bkg) &
     & bind(c,name='soca_vertconv_setup_f90')
  use iso_c_binding
  use soca_vertconv_mod
  use config_mod
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices
  
  implicit none
  
  integer(c_int), intent(inout) :: c_key_self   !< The Vertconv structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int),    intent(in) :: c_key_traj   !< trajectory
  integer(c_int),    intent(in) :: c_key_bkg    !< background  

  type(soca_vertconv), pointer :: self
  type(soca_field), pointer :: traj
  type(soca_field), pointer :: bkg  
  real(kind=kind_real), allocatable :: jac(:)

  integer :: isc, iec, jsc, jec, i, j, k, nl
  
  call soca_vertconv_registry%init()
  call soca_vertconv_registry%add(c_key_self)
  call soca_vertconv_registry%get(c_key_self, self)
  call soca_field_registry%get(c_key_traj, traj)
  call soca_field_registry%get(c_key_traj, bkg)  

  nl = size(bkg%hocn,3)
  
  ! Get configuration for vertical convolution
  self%lz   = config_get_real(c_conf, "Lz")

  ! Store trajectory and background
  call create_copy(self%traj, traj)
  call create_copy(self%bkg, bkg)
  
  ! Indices for compute domain (no halo)
  call geom_get_domain_indices(bkg%geom%ocean, "compute", isc, iec, jsc, jec)
  self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec
  
  ! Initialize local ocean depth from layer thickness
  allocate(self%z(isc:iec, jsc:jec, nl))
  do i = isc, iec
     do j = jsc, jec
        if (bkg%geom%ocean%mask2d(i,j).eq.1) then
           do k = 1, nl
              self%z(i,j,k) = sum(bkg%hocn(i,j,k:))
           end do
        end if
     end do
  end do

end subroutine c_soca_vertconv_setup

! ------------------------------------------------------------------------------
!> Destructor for Vertconv
subroutine c_soca_vertconv_delete(c_key_self) bind(c,name='soca_vertconv_delete_f90')
  use iso_c_binding
  use soca_vertconv_mod

  implicit none

  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure

  type(soca_vertconv), pointer :: self

  ! Deallocate trajectory and backgroun
  ! TODO
  ! Deallocate ocean depth array
  ! TODO
  
  call soca_vertconv_registry%get(c_key_self, self)  
  call soca_vertconv_registry%remove(c_key_self)

end subroutine c_soca_vertconv_delete

! ------------------------------------------------------------------------------
!> Multiplication
subroutine c_soca_vertconv_mult_f90(c_key_a, c_key_m, c_key_traj, c_key_self)&
     & bind(c,name='soca_vertconv_mult_f90')
  use iso_c_binding
  use soca_vertconv_mod
  use soca_fields
  use soca_utils
  use kinds
  use config_mod
  
  implicit none

  integer(c_int), intent(in) :: c_key_a     !< Increment in
  integer(c_int), intent(in) :: c_key_m     !< Increment out 
  integer(c_int), intent(in) :: c_key_traj  !< trajectory
  integer(c_int), intent(in) :: c_key_self  !< config
  
  type(soca_field), pointer :: dxa  ! in
  type(soca_field), pointer :: dxm  ! out
  type(soca_field), pointer :: traj  
  type(soca_vertconv),   pointer :: self
  
  call soca_field_registry%get(c_key_a, dxa)
  call soca_field_registry%get(c_key_m, dxm)
  call soca_field_registry%get(c_key_traj, traj)
  call soca_vertconv_registry%get(c_key_self, self)    

  !< Computes dxm = Vertconv dxa

  ! dxm = dxa
  !call copy(dxm, dxa) 

  ! Apply forward convolution operator to T & S
  call soca_conv(self, dxm, dxa)

end subroutine c_soca_vertconv_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_vertconv_multad_f90(c_key_m, c_key_a, c_key_traj, c_key_self)&
     & bind(c,name='soca_vertconv_multad_f90')
  use iso_c_binding
  use soca_vertconv_mod
  use soca_fields
  use kinds
  
  implicit none

  integer(c_int), intent(in) :: c_key_a     !< Increment out
  integer(c_int), intent(in) :: c_key_m     !< Increment in 
  integer(c_int), intent(in) :: c_key_traj  !< Trajectory
  integer(c_int), intent(in) :: c_key_self  !< config

  type(soca_field),      pointer :: dxa
  type(soca_field),      pointer :: dxm
  type(soca_field),      pointer :: traj
  type(soca_vertconv),   pointer :: self

  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_field_registry%get(c_key_traj,traj)
  call soca_vertconv_registry%get(c_key_self, self)    

  ! dxa = dxm
  !call copy(dxa,dxm)
  
  ! Apply adjoint of convolution operator
  call soca_conv_ad(self, dxm, dxa)

end subroutine c_soca_vertconv_multad_f90
