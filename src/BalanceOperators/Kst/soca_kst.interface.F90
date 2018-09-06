!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

! ------------------------------------------------------------------------------
!> Constructor for Kst
subroutine c_soca_kst_setup(c_key_self, c_conf, c_key_traj) bind(c,name='soca_kst_setup_f90')
  use iso_c_binding
  use soca_kst_mod
  use config_mod
  use soca_fields

  implicit none
  
  integer(c_int), intent(inout) :: c_key_self   !< The Kst structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int),    intent(in) :: c_key_traj  !< trajectory

  type(soca_kst), pointer :: self
  type(soca_field), pointer :: traj  
  real(kind=kind_real), allocatable :: jac(:)

  integer :: isc, iec, jsc, jec, i, j, k
  
  call soca_kst_registry%init()
  call soca_kst_registry%add(c_key_self)
  call soca_kst_registry%get(c_key_self, self)
  call soca_field_registry%get(c_key_traj, traj)

  ! Get configuration for S(T)
  self%dsdtmax      = config_get_real(c_conf,"dsdtmax")
  self%dsdzmin      = config_get_real(c_conf,"dsdzmin")
  self%dtdzmin      = config_get_real(c_conf,"dtdzmin")

  ! Allocate and compute the jacobian
  isc = traj%geom%ocean%G%isc
  iec = traj%geom%ocean%G%iec
  jsc = traj%geom%ocean%G%jsc
  jec = traj%geom%ocean%G%jec
  
  allocate(self%jacobian(isc:iec,jsc:jec,traj%geom%ocean%nzo))
  allocate(jac(traj%geom%ocean%nzo))
  self%jacobian=0.0
  do i = isc, iec
     do j = jsc, jec
        jac=0.0
        call soca_soft_jacobian(jac,&
             &traj%tocn(i,j,:),&
             &traj%socn(i,j,:),&
             &traj%hocn(i,j,:),&
             &self%dsdtmax, self%dsdzmin, self%dtdzmin)
        self%jacobian(i,j,:) = jac(:)
     end do
  end do
  deallocate(jac)

end subroutine c_soca_kst_setup

! ------------------------------------------------------------------------------
!> Destructor for Kst
subroutine c_soca_kst_delete(c_key_self) bind(c,name='soca_kst_delete_f90')
  use iso_c_binding
  use soca_kst_mod

  implicit none

  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure

  type(soca_kst), pointer :: self

  call soca_kst_registry%get(c_key_self, self)  
  if (allocated(self%jacobian)) deallocate(self%jacobian)
  call soca_kst_registry%remove(c_key_self)

end subroutine c_soca_kst_delete

! ------------------------------------------------------------------------------
!> Multiplication
subroutine c_soca_kst_mult_f90(c_key_a, c_key_m, c_key_traj, c_key_self)&
     & bind(c,name='soca_kst_mult_f90')
  use iso_c_binding
  use soca_kst_mod
  use soca_fields
  use soca_utils
  use kinds
  use config_mod
    
  !use soca_balanceop

  implicit none

  integer(c_int), intent(in) :: c_key_a     !< Increment in
  integer(c_int), intent(in) :: c_key_m     !< Increment out 
  integer(c_int), intent(in) :: c_key_traj  !< trajectory
  integer(c_int), intent(in) :: c_key_self  !< config
  
  type(soca_field), pointer :: dxa
  type(soca_field), pointer :: dxm
  type(soca_field), pointer :: traj  
  type(soca_kst),   pointer :: self
  
  real(kind=kind_real), allocatable :: dtv(:), dsv(:)
  integer :: isc, iec, jsc, jec, i, j, k
  
  call soca_field_registry%get(c_key_a, dxa)
  call soca_field_registry%get(c_key_m, dxm)
  call soca_field_registry%get(c_key_traj, traj)
  call soca_kst_registry%get(c_key_self, self)    

  !< Computes dxm = Kst dxa

  ! Indices for compute domain (no halo)
  isc = traj%geom%ocean%G%isc
  iec = traj%geom%ocean%G%iec
  jsc = traj%geom%ocean%G%jsc
  jec = traj%geom%ocean%G%jec

  call copy(dxm,dxa)  

  ! T-S balance    !< dxm = K dxa

  allocate(dsv(traj%geom%ocean%nzo),dtv(traj%geom%ocean%nzo))  
  dsv=0.0
  dtv=0.0    
  do i = isc, iec
     do j = jsc, jec
        dtv = dxa%tocn(i,j,:)
        call soca_soft_tl (dsv, dtv, self%jacobian(i,j,:))
        dxm%socn(i,j,:) = dxm%socn(i,j,:) + dsv
     end do
  end do
  deallocate(dsv,dtv)

end subroutine c_soca_kst_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_kst_multad_f90(c_key_m, c_key_a, c_key_traj, c_key_self)&
     & bind(c,name='soca_kst_multad_f90')
  use iso_c_binding
  use soca_kst_mod
  use soca_fields
  use kinds
  !use soca_balanceop

  implicit none

  integer(c_int), intent(in) :: c_key_a     !< Increment out
  integer(c_int), intent(in) :: c_key_m     !< Increment in 
  integer(c_int), intent(in) :: c_key_traj  !< Trajectory
  integer(c_int), intent(in) :: c_key_self  !< config

  type(soca_field),      pointer :: dxa
  type(soca_field),      pointer :: dxm
  type(soca_field),      pointer :: traj
  type(soca_kst),   pointer :: self

  real(kind=kind_real), allocatable :: dtv(:), dsv(:)
  integer :: isc, iec, jsc, jec, i, j, k

  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_field_registry%get(c_key_traj,traj)
  call soca_kst_registry%get(c_key_self, self)    

  !< Computes dxa = K^T dxm
  
  ! Indices for compute domain (no halo)
  isc = traj%geom%ocean%G%isc
  iec = traj%geom%ocean%G%iec
  jsc = traj%geom%ocean%G%jsc
  jec = traj%geom%ocean%G%jec

  call copy(dxa,dxm)
  
  ! T-S balance
  allocate(dsv(traj%geom%ocean%nzo),dtv(traj%geom%ocean%nzo))
  dsv=0.0
  dtv=0.0
  do i = isc, iec
     do j = jsc, jec
        dsv = dxm%socn(i,j,:)          
        call soca_soft_ad (dsv, dtv, self%jacobian(i,j,:))
        dxa%tocn(i,j,:) = dxa%tocn(i,j,:) + dtv
     end do
  end do
  deallocate(dtv,dsv)

end subroutine c_soca_kst_multad_f90
