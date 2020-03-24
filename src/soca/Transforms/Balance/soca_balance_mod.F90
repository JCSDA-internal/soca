! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_balance_mod

use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use soca_fields_mod
use soca_increment_mod
use soca_state_mod
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
   type(soca_state ), pointer :: traj                !> Trajectory
   integer                    :: isc, iec, jsc, jec  !> Compute domain
   type(soca_kst)             :: kst                 !> T/S balance
   type(soca_ksshts)          :: ksshts              !> SSH/T/S balance
   real(kind=kind_real), allocatable :: kct(:,:)     !> C/T Jacobian
end type soca_balance_config

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Initialization of the balance operator and its trajectory.
!> balances always used: T,S,SSH
!> optional balances depending on input fields: cicen
subroutine soca_balance_setup(f_conf, self, traj)
  type(fckit_configuration),   intent(in)  :: f_conf
  type(soca_balance_config), intent(inout) :: self
  type(soca_state),   target, intent(in)  :: traj

  integer :: isc, iec, jsc, jec, i, j, k, nl
  real(kind=kind_real), allocatable :: jac(:)
  type(soca_field), pointer :: tocn, socn, hocn, cicen, mld, layer_depth

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

  ! Get required fields
  call traj%get("tocn", tocn)
  call traj%get("socn", socn)
  call traj%get("hocn", hocn)
  call traj%get("mld", mld)
  call traj%get("layer_depth", layer_depth)
  if (traj%has("cicen"))  call traj%get("cicen", cicen)

  ! allocate space
  nl = hocn%nz
  allocate(self%kst%jacobian(isc:iec,jsc:jec,traj%geom%nzo))
  allocate(jac(nl))
  self%kst%jacobian=0.0

  ! Compute and store Jacobian of Kst
  do i = isc, iec
     do j = jsc, jec
        jac=0.0
        call soca_soft_jacobian(jac,&
             &tocn%val(i,j,:),&
             &socn%val(i,j,:),&
             &hocn%val(i,j,:),&
             &self%kst%dsdtmax, self%kst%dsdzmin, self%kst%dtdzmin)
        ! Set Jacobian to 0 above mixed layer
        do k=1,nl
           if (layer_depth%val(i,j,k) < mld%val(i,j,1)) then
              jac(k) = 0.0_kind_Real
           end if
        end do
        self%kst%jacobian(i,j,:) = jac(:)
        self%kst%jacobian(i,j,1:self%kst%nlayers) =  0.0_kind_real
     end do
  end do
  deallocate(jac)

  ! Compute Jacobian of Ksshts
  allocate(self%ksshts%kssht, mold=self%kst%jacobian)
  allocate(self%ksshts%ksshs, mold=self%kst%jacobian)
  allocate(jac(2))
  self%ksshts%kssht=0.0
  self%ksshts%ksshs=0.0
  do i = isc, iec
     do j = jsc, jec
        do k = 1, nl
           call soca_steric_jacobian (jac, &
                tocn%val(i,j,k), &
                socn%val(i,j,k), &
                &layer_depth%val(i,j,k),&
                &hocn%val(i,j,k),&
                &traj%geom%lon(i,j),&
                &traj%geom%lat(i,j))
           self%ksshts%kssht(i,j,k) = jac(1)
           self%ksshts%ksshs(i,j,k) = jac(2)
        end do
     end do
  end do
  deallocate(jac)

  ! Compute Kct
  if (traj%has("cicen")) then
    allocate(self%kct(isc:iec,jsc:jec))
    self%kct = 0.0_kind_real
    do i = isc, iec
      do j = jsc, jec
          if (sum(cicen%val(i,j,:)) > 1.0e-3_kind_real) then
            self%kct = -0.01d0 ! TODO: Insert regression coef
          end if
      end do
    end do
  end if

end subroutine soca_balance_setup

! ------------------------------------------------------------------------------
!> Destructor for the balance oprator
subroutine soca_balance_delete(self)
  type(soca_balance_config), intent(inout) :: self

  ! the following always exist
  nullify(self%traj)
  deallocate(self%kst%jacobian)
  deallocate(self%ksshts%kssht)
  deallocate(self%ksshts%ksshs)

  ! only exists if cicen was given
  if (allocated(self%kct)) deallocate(self%kct)

end subroutine soca_balance_delete

! ------------------------------------------------------------------------------
! Apply forward balance operator
subroutine soca_balance_mult(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_increment),      intent(in) :: dxa
  type(soca_increment),   intent(inout) :: dxm

  type(soca_field), pointer :: fld_m, fld_a
  type(soca_field), pointer :: tocn_a, socn_a

  integer :: i, j, k, n

  !>    [ I       0   0  0 ]
  !>    [ Kst     I   0  0 ]
  !> K= [ Ketat Ketas I  0 ]
  !>    [ Kct     0   0  I ]

  call dxa%get("tocn",tocn_a)
  call dxa%get("socn",socn_a)

  do n=1, size(dxm%fields)
    fld_m => dxm%fields(n)
    fld_a => dxa%fields(n)

    do i = self%isc, self%iec
      do j = self%jsc, self%jec
        select case(fld_m%name)
        case default
          fld_m%val(i,j,:) = fld_a%val(i,j,:)

        case("socn") ! Salinity
          fld_m%val(i,j,:) = fld_a%val(i,j,:) + &
            & self%kst%jacobian(i,j,:) * tocn_a%val(i,j,:)

        case ("ssh") ! SSH
          fld_m%val(i,j,:) = fld_a%val(i,j,:)
          do k = 1, tocn_a%nz
            fld_m%val(i,j,:) = fld_m%val(i,j,:) + &
              & self%ksshts%kssht(i,j,k) * tocn_a%val(i,j,k) + &
              & self%ksshts%ksshs(i,j,k) * socn_a%val(i,j,k)
          end do

        case ("cicen") ! Ice fraction
          do k = 1, fld_m%nz
            fld_m%val(i,j,k) = fld_a%val(i,j,k) + &
              & self%kct(i,j) * tocn_a%val(i,j,1)
          end do

        end select
      end do
    end do
  end do
end subroutine soca_balance_mult

! ------------------------------------------------------------------------------
! Apply backward balance operator
subroutine soca_balance_multad(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_increment),      intent(in) :: dxm
  type(soca_increment),   intent(inout) :: dxa

  type(soca_field), pointer :: fld_a, fld_m
  type(soca_field), pointer :: socn_m, ssh_m, cicen_m
  integer :: i, j, n

  cicen_m => null()

  call dxm%get("socn", socn_m)
  call dxm%get("ssh",  ssh_m)
  if (dxm%has("cicen")) call dxm%get("cicen",cicen_m)

  do n = 1, size(dxa%fields)
    fld_a => dxa%fields(n)
    fld_m => dxm%fields(n)

    do i = self%isc, self%iec
      do j = self%jsc, self%jec
        select case(fld_a%name)
        case default
          fld_a%val(i,j,:) = fld_m%val(i,j,:)

        case ("tocn") ! Temperature
          fld_a%val(i,j,:) = fld_m%val(i,j,:) + &
            & self%kst%jacobian(i,j,:) * socn_m%val(i,j,:) + &
            & self%ksshts%kssht(i,j,:) * ssh_m%val(i,j,1)

          if (associated(cicen_m)) then ! use cicen only if present
            fld_a%val(i,j,1) = fld_a%val(i,j,1) + &
              & self%kct(i,j) * sum(cicen_m%val(i,j,:))
          end if

        case ("socn") ! Salinity
          fld_a%val(i,j,:) = fld_m%val(i,j,:) + &
            & self%ksshts%ksshs(i,j,:) * ssh_m%val(i,j, 1)

        end select
      end do
    end do
  end do
end subroutine soca_balance_multad

! ------------------------------------------------------------------------------
! Apply inverse of the forward balance operator
subroutine soca_balance_multinv(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_increment),      intent(in) :: dxm
  type(soca_increment),   intent(inout) :: dxa

  integer :: i, j, k, n
  type(soca_field), pointer :: fld_m, fld_a
  type(soca_field), pointer :: tocn_m, socn_m

  call dxm%get("tocn", tocn_m)
  call dxm%get("socn", socn_m)

  do n = 1, size(dxa%fields)
    fld_a => dxa%fields(n)
    fld_m => dxm%fields(n)

    do i = self%isc, self%iec
      do j = self%jsc, self%jec
        select case(fld_a%name)
        case default
          fld_a%val(i,j,:) = fld_m%val(i,j,:)

        case ('socn') ! Salinity
          fld_a%val(i,j,:) = fld_m%val(i,j,:) - &
            & self%kst%jacobian(i,j,:) * tocn_m%val(i,j,:)

        case ('ssh') ! SSH
          fld_a%val(i,j, :) = fld_m%val(i,j, :)
          do k = 1, tocn_m%nz
            fld_a%val(i,j,:) = fld_a%val(i,j,:) + &
              & ( self%ksshts%ksshs(i,j,k) * self%kst%jacobian(i,j,k) - &
              & self%ksshts%kssht(i,j,k) ) *  tocn_m%val(i,j,k) - &
              & self%ksshts%ksshs(i,j,k) * socn_m%val(i,j,k)
          end do

        case ('cicen') ! Ice fraction
          fld_a%val(i,j,:) =  fld_m%val(i,j,:)
          do k = 1, fld_m%nz
            fld_a%val(i,j,k) = fld_a%val(i,j,k) - &
              & self%kct(i,j) * tocn_m%val(i,j,1)
          end do

        end select
      end do
    end do
  end do
end subroutine soca_balance_multinv

! ------------------------------------------------------------------------------
! Apply inverse of the backward balance operator
subroutine soca_balance_multinvad(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_increment),   intent(inout) :: dxm
  type(soca_increment),      intent(in) :: dxa

  integer :: i, j, n
  type(soca_field), pointer :: fld_a, fld_m
  type(soca_field), pointer :: socn_a, ssh_a, cicen_a

  cicen_a => null()

  call dxa%get("socn", socn_a)
  call dxa%get("ssh",  ssh_a)
  if (dxa%has("cicen")) call dxa%get("cicen",cicen_a)

  do n = 1, size(dxm%fields)
    fld_m => dxm%fields(n)
    fld_a => dxa%fields(n)

    do i = self%isc, self%iec
      do j = self%jsc, self%jec
        select case (fld_m%name)
        case default
          fld_m%val(i,j,:) = fld_a%val(i,j,:)

        case ('tocn') ! Temperature
          fld_m%val(i,j,:) = fld_a%val(i,j,:) &
            & - self%kst%jacobian(i,j,:) * socn_a%val(i,j,:) &
            & + ( self%ksshts%ksshs(i,j,:) * self%kst%jacobian(i,j,:) &
            &     - self%ksshts%kssht(i,j,:) ) * ssh_a%val(i,j,1)

          if (associated(cicen_a)) then ! use cicen only if present
            fld_m%val(i,j,1) = fld_m%val(i,j,1) &
              & - self%kct(i,j) * sum(cicen_a%val(i,j,:))
          end if

        case ('socn') ! Salinity
          fld_m%val(i,j,:) = fld_a%val(i,j,:) - &
            & self%ksshts%ksshs(i,j,:) * ssh_a%val(i,j,1)

        end select
      end do
    end do
  end do
end subroutine soca_balance_multinvad

end module soca_balance_mod
