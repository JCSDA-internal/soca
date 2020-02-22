! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_balance_mod

use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use soca_fields_mod, only: soca_fields
use soca_kst_mod, only: soca_kst, soca_soft_jacobian
use soca_ksshts_mod, only: soca_ksshts, soca_steric_jacobian

implicit none

private
public :: soca_balance_config, &
          soca_balance_setup, soca_balance_delete, &
          soca_balance_mult, soca_balance_multad, &
          soca_balance_multinv, soca_balance_multinvad

!> Fortran derived type to hold configuration D
type :: soca_balance_config
   type(soca_fields), pointer :: traj                !> Trajectory
   integer                    :: isc, iec, jsc, jec  !> Compute domain
   type(soca_kst)             :: kst                 !> T/S balance
   type(soca_ksshts)          :: ksshts              !> SSH/T/S balance
   real(kind=kind_real), allocatable :: kct(:,:)     !> C/T Jacobian
end type soca_balance_config

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Initialization of the balance operator and its trajectory
subroutine soca_balance_setup(f_conf, self, traj)
  type(fckit_configuration),   intent(in)  :: f_conf
  type(soca_balance_config), intent(inout) :: self
  type(soca_fields),   target, intent(in)  :: traj

  integer :: isc, iec, jsc, jec, i, j, k, nl
  real(kind=kind_real), allocatable :: jac(:)

  ! Number of ocean layer
  nl = size(traj%hocn,3)

  ! Store trajectory
  self%traj => traj

  ! Indices for compute domain
  isc=traj%geom%isc; iec=traj%geom%iec
  jsc=traj%geom%jsc; jec=traj%geom%jec

  self%isc=isc; self%iec=iec
  self%jsc=jsc; self%jec=jec

  ! Get configuration for Kst

  call f_conf%get_or_die("dsdtmax", self%kst%dsdtmax)
  call f_conf%get_or_die("dsdzmin", self%kst%dsdzmin)
  call f_conf%get_or_die("dtdzmin", self%kst%dtdzmin)
  call f_conf%get_or_die("nlayers", self%kst%nlayers) ! Set jac to 0 in the
                                                      ! nlayers top layers

  ! Compute and store Jacobian of Kst
  allocate(self%kst%jacobian(isc:iec,jsc:jec,traj%geom%nzo))
  allocate(jac(nl))
  self%kst%jacobian=0.0
  do i = isc, iec
     do j = jsc, jec
        jac=0.0
        call soca_soft_jacobian(jac,&
             &traj%tocn(i,j,:),&
             &traj%socn(i,j,:),&
             &traj%hocn(i,j,:),&
             &self%kst%dsdtmax, self%kst%dsdzmin, self%kst%dtdzmin)
        ! Set Jacobian to 0 above mixed layer
        do k=1,nl
           if (self%traj%layer_depth(i,j,k)<self%traj%mld(i,j)) then
              jac(k) = 0.0_kind_Real
           end if
        end do
        self%kst%jacobian(i,j,:) = jac(:)
        self%kst%jacobian(i,j,1:self%kst%nlayers) =  0.0_kind_real
     end do
  end do
  deallocate(jac)

  ! Compute Jacobian of Ksshts
  allocate(self%ksshts%kssht(isc:iec,jsc:jec,traj%geom%nzo))
  allocate(self%ksshts%ksshs(isc:iec,jsc:jec,traj%geom%nzo))
  allocate(jac(2))
  self%ksshts%kssht=0.0
  self%ksshts%ksshs=0.0
  do i = isc, iec
     do j = jsc, jec
        do k = 1, nl
           call soca_steric_jacobian (jac, &
                &traj%tocn(i,j,k),&
                &traj%socn(i,j,k),&
                &self%traj%layer_depth(i,j,k),&
                &traj%hocn(i,j,k),&
                &traj%geom%lon(i,j),&
                &traj%geom%lat(i,j))
           self%ksshts%kssht(i,j,k) = jac(1)
           self%ksshts%ksshs(i,j,k) = jac(2)
        end do
     end do
  end do
  deallocate(jac)

  ! Compute Kst
  allocate(self%kct(isc:iec,jsc:jec))
  self%kct = 0.0_kind_real
  do i = isc, iec
     do j = jsc, jec
        if (sum(traj%seaice%cicen(i,j,2:)) > 1.0e-3_kind_real) then
           self%kct = -0.01d0 ! TODO: Insert regression coef
        end if
     end do
  end do

end subroutine soca_balance_setup

! ------------------------------------------------------------------------------
!> Destructor for the balance oprator
subroutine soca_balance_delete(self)
  type(soca_balance_config), intent(inout) :: self

  nullify(self%traj)
  deallocate(self%kst%jacobian)
  deallocate(self%ksshts%kssht)
  deallocate(self%ksshts%ksshs)
  deallocate(self%kct)

end subroutine soca_balance_delete

! ------------------------------------------------------------------------------
! Apply forward balance operator
subroutine soca_balance_mult(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_fields),         intent(in) :: dxa
  type(soca_fields),      intent(inout) :: dxm

  integer :: i, j, k
  real(kind=kind_real) :: deta, dxc

  !>    [ I       0   0  0 ]
  !>    [ Kst     I   0  0 ]
  !> K= [ Ketat Ketas I  0 ]
  !>    [ Kct     0   0  I ]

  do i = self%isc, self%iec
     do j = self%jsc, self%jec
        ! Temperature
        dxc = sum(dxa%seaice%cicen(i,j,2:))
        dxm%tocn(i,j,1) = dxa%tocn(i,j,1) + self%kct(i,j) * dxc
        dxm%tocn(i,j,:) = dxa%tocn(i,j,:)

        ! Salinity
        dxm%socn(i,j,:) = dxa%socn(i,j,:) +&
             &self%kst%jacobian(i,j,:) * dxa%tocn(i,j,:)

        ! SSH
        deta = 0.0_kind_real
        do k = 1, size(self%traj%hocn,3)
           deta = deta + self%ksshts%kssht(i,j,k) * dxa%tocn(i,j,k) +&
                &self%ksshts%ksshs(i,j,k) * dxa%socn(i,j,k)
        end do
        dxm%ssh(i,j) = dxa%ssh(i,j) + deta

        ! Ice fraction
        dxm%seaice%cicen(i,j,:) =  dxa%seaice%cicen(i,j,:)
        do k = 1, size(self%traj%seaice%hicen,3)
           dxm%seaice%cicen(i,j,k+1) =  dxm%seaice%cicen(i,j,k+1) +&
                & self%kct(i,j) * dxa%tocn(i,j,1)
        end do

        ! Ice thickness
        dxm%seaice%hicen(i,j,:) =  dxa%seaice%hicen(i,j,:)
     end do
  end do
  ! Surface fields
  call dxm%ocnsfc%copy(dxa%ocnsfc)

end subroutine soca_balance_mult

! ------------------------------------------------------------------------------
! Apply backward balance operator
subroutine soca_balance_multad(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_fields),         intent(in) :: dxm
  type(soca_fields),      intent(inout) :: dxa

  integer :: i, j

  do i = self%isc, self%iec
     do j = self%jsc, self%jec
        ! Temperature
        dxa%tocn(i,j,1) = dxm%tocn(i,j,1) + &
             &self%kst%jacobian(i,j,1) * dxm%socn(i,j,1) + &
             &self%ksshts%kssht(i,j,1) * dxm%ssh(i,j) +&
             &self%kct(i,j) * sum(dxm%seaice%cicen(i,j,2:))
        dxa%tocn(i,j,2:) = dxm%tocn(i,j,2:) + &
             &self%kst%jacobian(i,j,2:) * dxm%socn(i,j,2:) + &
             &self%ksshts%kssht(i,j,2:) * dxm%ssh(i,j)
        ! Salinity
        dxa%socn(i,j,:) = dxm%socn(i,j,:) + &
             &self%ksshts%ksshs(i,j,:) * dxm%ssh(i,j)
        ! SSH
        dxa%ssh(i,j)    = dxm%ssh(i,j)
        ! Ice fraction
        dxa%seaice%cicen(i,j,:) =  dxm%seaice%cicen(i,j,:)
        ! Ice thickness
        dxa%seaice%hicen(i,j,:) =  dxm%seaice%hicen(i,j,:)
     end do
  end do
  ! Surface fields
  call dxa%ocnsfc%copy(dxm%ocnsfc)

end subroutine soca_balance_multad

! ------------------------------------------------------------------------------
! Apply inverse of the forward balance operator
subroutine soca_balance_multinv(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_fields),         intent(in) :: dxm
  type(soca_fields),      intent(inout) :: dxa

  real(kind=kind_real) :: deta
  integer :: i, j, k

  do i = self%isc, self%iec
     do j = self%jsc, self%jec
        ! Temperature
        dxa%tocn(i,j,:) = dxm%tocn(i,j,:)
        ! Salinity
        dxa%socn(i,j,:) = dxm%socn(i,j,:) -&
             &self%kst%jacobian(i,j,:) * dxm%tocn(i,j,:)
        ! SSH
        deta = 0.0d0
        do k = 1, size(self%traj%hocn,3)
           deta = deta + ( self%ksshts%ksshs(i,j,k) * self%kst%jacobian(i,j,k) - &
                &self%ksshts%kssht(i,j,k) ) * dxm%tocn(i,j,k) - &
                &self%ksshts%ksshs(i,j,k) * dxm%socn(i,j,k)
        end do
        dxa%ssh(i,j)    = dxm%ssh(i,j) + deta
        ! Ice fraction
        dxa%seaice%cicen(i,j,:) =  dxm%seaice%cicen(i,j,:)
        do k = 1, size(self%traj%seaice%hicen,3)
           dxa%seaice%cicen(i,j,k+1) =  dxa%seaice%cicen(i,j,k+1) -&
                & self%kct(i,j) * dxm%tocn(i,j,1)
        end do
        ! Ice thickness
        dxa%seaice%hicen(i,j,:) =  dxm%seaice%hicen(i,j,:)
     end do
  end do
  ! Surface fields
  call dxa%ocnsfc%copy(dxm%ocnsfc)
end subroutine soca_balance_multinv

! ------------------------------------------------------------------------------
! Apply inverse of the backward balance operator
subroutine soca_balance_multinvad(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_fields),      intent(inout) :: dxm
  type(soca_fields),         intent(in) :: dxa

  integer :: i, j

  do i = self%isc, self%iec
     do j = self%jsc, self%jec
        ! Ice thickness
        dxm%seaice%hicen(i,j,:) =  dxa%seaice%hicen(i,j,:)
        ! Ice fraction
        dxm%seaice%cicen(i,j,:) =  dxa%seaice%cicen(i,j,:)
        ! Temperature
        dxm%tocn(i,j,1) = dxa%tocn(i,j,1) &
             & - self%kst%jacobian(i,j,1) * dxa%socn(i,j,1) &
             & + ( self%ksshts%ksshs(i,j,1) * self%kst%jacobian(i,j,1) &
             &     - self%ksshts%kssht(i,j,1) ) * dxa%ssh(i,j) &
             & - self%kct(i,j) * sum(dxa%seaice%cicen(i,j,2:))
        dxm%tocn(i,j,2:) = dxa%tocn(i,j,2:) &
             & - self%kst%jacobian(i,j,2:) * dxa%socn(i,j,2:) &
             & + ( self%ksshts%ksshs(i,j,2:) * self%kst%jacobian(i,j,2:) &
             &     - self%ksshts%kssht(i,j,2:) ) * dxa%ssh(i,j)
        ! Salinity
        dxm%socn(i,j,:) = dxa%socn(i,j,:) - &
             &self%ksshts%ksshs(i,j,:) * dxa%ssh(i,j)
        ! SSH
        dxm%ssh(i,j)    = dxa%ssh(i,j)
     end do
  end do
  ! Surface fields
  call dxm%ocnsfc%copy(dxa%ocnsfc)

end subroutine soca_balance_multinvad

end module soca_balance_mod
