!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

! ------------------------------------------------------------------------------
!> Constructor for Ksshts
subroutine c_soca_ksshts_setup(c_key_self, c_conf) bind(c,name='soca_ksshts_setup_f90')
  use iso_c_binding
  use soca_ksshts_mod

  integer(c_int), intent(inout) :: c_key_self   !< The Ksshts structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration

  type(soca_ksshts_config), pointer :: self

  call soca_ksshts_registry%init()
  call soca_ksshts_registry%add(c_key_self)
  call soca_ksshts_registry%get(c_key_self, self)

  call soca_ksshts_setup(c_conf, self)

end subroutine c_soca_ksshts_setup

! ------------------------------------------------------------------------------
!> Destructor for Ksshts
subroutine c_soca_ksshts_delete(c_key_self) bind(c,name='soca_ksshts_delete_f90')
  use iso_c_binding
  use soca_ksshts_mod

  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure

  call soca_ksshts_registry%remove(c_key_self)

end subroutine c_soca_ksshts_delete

! ------------------------------------------------------------------------------
!> Multiplication
subroutine c_soca_ksshts_mult_f90(c_key_a, c_key_m, c_key_traj)&
     &bind(c,name='soca_ksshts_mult_f90')
  use iso_c_binding
  use soca_ksshts_mod
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
  real(kind=kind_real) :: tb, sb, dt, ds,deta, z, p, h

  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_field_registry%get(c_key_traj,traj)  

  !< Computes dxm = Ksshts dxa

  ! Indices for compute domain (no halo)
  isc = traj%geom%ocean%G%isc
  iec = traj%geom%ocean%G%iec
  jsc = traj%geom%ocean%G%jsc
  jec = traj%geom%ocean%G%jec

  call copy(dxm,dxa)
  ! Steric height/density balance
  do i = isc, iec
     do j = jsc, jec
        do k = 1, traj%geom%ocean%nzo
           tb=traj%tocn(i,j,k)
           sb=traj%socn(i,j,k)
           dt=dxa%tocn(i,j,k)
           ds=dxa%socn(i,j,k)
           if (k.eq.1) then
              z=traj%hocn(i,j,k)
           else
              z=sum(traj%hocn(i,j,1:k-1))+0.5_kind_real*traj%hocn(i,j,k)
           end if
           h=traj%hocn(i,j,k)
           p=z
           call soca_steric_tl(deta, dt, ds, tb, sb, p, h,&
                &traj%geom%ocean%lon(i,j), traj%geom%ocean%lat(i,j))
           dxm%ssh(i,j)=dxm%ssh(i,j)+deta    
        end do
     end do
  end do


end subroutine c_soca_ksshts_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_ksshts_multad_f90(c_key_m, c_key_a, c_key_traj) bind(c,name='soca_ksshts_multad_f90')
  use iso_c_binding
  use soca_ksshts_mod
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
  real(kind=kind_real) :: tb, sb, dt, ds,deta, z, p, h

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
  ! Steric height/density balance
  do i = isc, iec
     do j = jsc, jec
        do k = traj%geom%ocean%nzo, 1, -1
           tb=traj%tocn(i,j,k)
           sb=traj%socn(i,j,k)
           if (k.eq.1) then
              z=traj%hocn(i,j,k)
           else
              z=sum(traj%hocn(i,j,1:k-1))+0.5_kind_real*traj%hocn(i,j,k)
           end if
           h=traj%hocn(i,j,k)
           p=z
           deta=dxm%ssh(i,j)
           call soca_steric_ad(deta, dt, ds, tb, sb, p, h,&
                &traj%geom%ocean%lon(i,j), traj%geom%ocean%lat(i,j))
           dxa%tocn(i,j,k) = dxa%tocn(i,j,k)+dt
           dxa%socn(i,j,k) = dxa%socn(i,j,k)+ds
        end do
     end do
  end do
  
end subroutine c_soca_ksshts_multad_f90
