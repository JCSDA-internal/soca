!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

! ------------------------------------------------------------------------------
!> Constructor for Ktc
subroutine c_soca_ktc_setup(c_key_self, c_conf) bind(c,name='soca_ktc_setup_f90')
  use iso_c_binding
  use soca_ktc_mod

  integer(c_int), intent(inout) :: c_key_self   !< The Ktc structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration

  type(soca_ktc_config), pointer :: self

  call soca_ktc_registry%init()
  call soca_ktc_registry%add(c_key_self)
  call soca_ktc_registry%get(c_key_self, self)

  call soca_ktc_setup(c_conf, self)

end subroutine c_soca_ktc_setup

! ------------------------------------------------------------------------------
!> Destructor for Ktc
subroutine c_soca_ktc_delete(c_key_self) bind(c,name='soca_ktc_delete_f90')
  use iso_c_binding
  use soca_ktc_mod

  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure

  call soca_ktc_registry%remove(c_key_self)

end subroutine c_soca_ktc_delete

! ------------------------------------------------------------------------------
!> Multiplication
subroutine c_soca_ktc_mult_f90(c_key_a, c_key_m, c_key_traj)&
     &bind(c,name='soca_ktc_mult_f90')
  use iso_c_binding
  use soca_ktc_mod
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

  integer :: isc, iec, jsc, jec, i, j, k
  real(kind=kind_real), allocatable :: dcn(:), cnb(:)
  real(kind=kind_real) :: dt  
  real(kind=kind_real) :: tb   !< Background potential temperature [C]
  real(kind=kind_real) :: sb   !< Background practical salinity [psu]

  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_field_registry%get(c_key_traj,traj)  

  !< Computes dxm = Ktc dxa

  ! Indices for compute domain (no halo)
  isc = traj%geom%ocean%G%isc
  iec = traj%geom%ocean%G%iec
  jsc = traj%geom%ocean%G%jsc
  jec = traj%geom%ocean%G%jec

  call copy(dxm,dxa)  

  ! T-C balance    !< dxm = K dxa
  allocate(dcn(traj%geom%ocean%ncat),cnb(traj%geom%ocean%ncat))
  do i = isc, iec
     do j = jsc, jec
        if (sum(traj%cicen(i,j,2:)).gt.0.0) then
           dxm%tocn(i,j,1) = dxm%tocn(i,j,1) - sum(dxa%cicen(i,j,2:))
        end if
     end do
  end do

  deallocate(dcn)
  
end subroutine c_soca_ktc_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_ktc_multad_f90(c_key_m, c_key_a, c_key_traj) bind(c,name='soca_ktc_multad_f90')
  use iso_c_binding
  use soca_ktc_mod
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

  integer :: isc, iec, jsc, jec, i, j, k
  real(kind=kind_real), allocatable :: dcn(:), cnb(:)
  real(kind=kind_real) :: dt    
  real(kind=kind_real) :: tb   !< Background potential temperature [C]
  real(kind=kind_real) :: sb   !< Background practical salinity [psu]

  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_field_registry%get(c_key_traj,traj)  

  !< Computes dxa = K^T dxm

  ! Indices for compute domain (no halo)
  isc = traj%geom%ocean%G%isc
  iec = traj%geom%ocean%G%iec
  jsc = traj%geom%ocean%G%jsc
  jec = traj%geom%ocean%G%jec

  call copy(dxa,dxm)

  allocate(dcn(traj%geom%ocean%ncat),cnb(traj%geom%ocean%ncat))
  do i = isc, iec
     do j = jsc, jec
        if (sum(traj%cicen(i,j,2:)).gt.0.0) then
           !dxa%cicen(i,j,1) = 0.0
           do k = 1, traj%geom%ocean%ncat+1
              dxa%cicen(i,j,k) = dxa%cicen(i,j,k)-dxm%tocn(i,j,1)
           end do
        end if
     end do
  end do

  deallocate(dcn, cnb)


end subroutine c_soca_ktc_multad_f90
