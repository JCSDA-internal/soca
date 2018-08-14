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

  call soca_kst_registry%remove(c_key_self)

end subroutine c_soca_kst_delete

! ------------------------------------------------------------------------------
!> Multiplication
subroutine c_soca_kst_mult_f90(c_key_a, c_key_m, c_key_traj)&
     &bind(c,name='soca_kst_mult_f90')
  use iso_c_binding
  use soca_kst_mod
  use soca_fields
  use soca_utils
  use kinds
  !use soca_balanceop

  implicit none
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out 
  integer(c_int), intent(in) :: c_key_traj  !<    "   to trajectory

  type(soca_field), pointer :: dxa
  type(soca_field), pointer :: dxm
  type(soca_field), pointer :: traj  

  real(kind=kind_real), allocatable :: dtv(:), dsv(:)
  integer :: isc, iec, jsc, jec, i, j, k
  
  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_field_registry%get(c_key_traj,traj)  

  !< Computes dxm = Kst dxa
  
  ! Indices for compute domain (no halo)
  call soca_compute_domain(isc, iec, jsc, jec, traj%geom)  
  
  ! T-S balance    !< dxm = K dxa
  allocate(dsv(traj%geom%ocean%nzo),dtv(traj%geom%ocean%nzo))
  dsv=0.0
  dtv=0.0    
  do i = isc, iec
     do j = jsc, jec
        dtv = dxa%tocn(i,j,:)
        call soca_soft_tl (dsv,dtv,&
             &traj%tocn(i,j,:),&
             &traj%socn(i,j,:),&
             &traj%hocn(i,j,:))
        dxm%socn(i,j,:) = dxm%socn(i,j,:) + dsv
     end do
  end do
  deallocate(dsv,dtv)

end subroutine c_soca_kst_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_kst_multad_f90(c_key_m, c_key_a, c_key_traj) bind(c,name='soca_kst_multad_f90')
  use iso_c_binding
  use soca_kst_mod
  use soca_fields
  use kinds
  !use soca_balanceop

  implicit none
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment out
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment in 
  integer(c_int), intent(in) :: c_key_traj  !<    "   to trajectory

  type(soca_field), pointer :: dxa
  type(soca_field), pointer :: dxm
  type(soca_field), pointer :: traj  

  real(kind=kind_real), allocatable :: dtv(:), dsv(:)
  integer :: isc, iec, jsc, jec, i, j, k

  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_field_registry%get(c_key_traj,traj)  

  !< Computes dxa = K^T dxm
  
  ! Indices for compute domain (no halo)
  call soca_compute_domain(isc, iec, jsc, jec, traj%geom)  

  ! T-S balance
  allocate(dsv(traj%geom%ocean%nzo),dtv(traj%geom%ocean%nzo))
  dsv=0.0
  dtv=0.0
  do i = isc, iec
     do j = jsc, jec
        dsv = dxm%socn(i,j,:)          
        call soca_soft_ad (dsv,dtv,&
             &traj%tocn(i,j,:),&
             &traj%socn(i,j,:),&
             &traj%hocn(i,j,:))
        dxa%tocn(i,j,:) = dxa%tocn(i,j,:) + dtv
     end do
  end do
  deallocate(dtv,dsv)

end subroutine c_soca_kst_multad_f90
