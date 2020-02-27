! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_balance_mod

use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use soca_fields_mod, only: soca_fields, soca_field
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
  type(soca_field), pointer :: tocn, socn, hocn, cicen

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
  call traj%get("tocn", tocn)
  call traj%get("socn", socn)
  call traj%get("hocn", hocn)
  call traj%get("cicen", cicen)
  nl = hocn%nz

  allocate(self%kst%jacobian(isc:iec,jsc:jec,traj%geom%nzo))
  allocate(jac(nl))
  self%kst%jacobian=0.0
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
                tocn%val(i,j,k), &
                socn%val(i,j,k), &
                &self%traj%layer_depth(i,j,k),&
                &hocn%val(i,j,k),&
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
        if (sum(cicen%val(i,j,2:)) > 1.0e-3_kind_real) then
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

  type(soca_field), pointer :: fld_m, fld_a
  type(soca_field), pointer :: tocn_m, tocn_a
  type(soca_field), pointer :: socn_m, socn_a
  type(soca_field), pointer :: ssh_m, ssh_a
  type(soca_field), pointer :: cicen_m, cicen_a

  integer :: i, j, k
  real(kind=kind_real) :: deta, dxc

  !>    [ I       0   0  0 ]
  !>    [ Kst     I   0  0 ]
  !> K= [ Ketat Ketas I  0 ]
  !>    [ Kct     0   0  I ]

  ! TODO skip parts of balance for variables that aren't present (i.e u,v)
  call dxm%get("tocn",tocn_m)
  call dxa%get("tocn",tocn_a)
  call dxm%get("socn",socn_m)
  call dxa%get("socn",socn_a)
  call dxm%get("ssh",ssh_m)
  call dxa%get("ssh",ssh_a)
  call dxm%get("cicen",cicen_m)
  call dxa%get("cicen",cicen_a)  

  do i = self%isc, self%iec
     do j = self%jsc, self%jec
        ! Temperature
        dxc = sum(cicen_a%val(i,j,2:))
        tocn_m%val(i,j,1) = tocn_a%val(i,j,1) + self%kct(i,j) * dxc
        tocn_m%val(i,j,:) = tocn_a%val(i,j,:)

        ! Salinity
        socn_m%val(i,j,:) = socn_a%val(i,j,:) +&
             &self%kst%jacobian(i,j,:) * tocn_a%val(i,j,:)

        ! SSH
        deta = 0.0_kind_real
        do k = 1, tocn_a%nz
           deta = deta + self%ksshts%kssht(i,j,k) * tocn_a%val(i,j,k) +&
                &self%ksshts%ksshs(i,j,k) * socn_a%val(i,j,k)
        end do
        ssh_m%val(i,j,:) = ssh_a%val(i,j,:) + deta

        ! Ice fraction
        cicen_m%val(i,j,:) = cicen_a%val(i,j,:)
        do k = 1, cicen_m%nz-1
           cicen_m%val(i,j,k+1) = cicen_m%val(i,j,k+1) +&
                & self%kct(i,j) * tocn_a%val(i,j,1)
        end do
     end do
  end do

  ! copy surface fields 
  ! TODO do this to all fields not explicitly handled above?
  do i=1, size(dxa%fields)
    fld_a => dxa%fields(i)
    call dxm%get(fld_a%name, fld_m)
    select case(fld_a%name)
    case ('sw','lw','lhf','shf','us','hicen')
      fld_m%val = fld_a%val
    end select
  end do
  
end subroutine soca_balance_mult

! ------------------------------------------------------------------------------
! Apply backward balance operator
subroutine soca_balance_multad(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_fields),         intent(in) :: dxm
  type(soca_fields),      intent(inout) :: dxa

  type(soca_field), pointer :: fld_a, fld_m
  type(soca_field), pointer :: tocn_a, tocn_m
  type(soca_field), pointer :: socn_a, socn_m
  type(soca_field), pointer :: ssh_a, ssh_m
  type(soca_field), pointer :: cicen_m  
  integer :: i, j

  call dxa%get("tocn", tocn_a)
  call dxm%get("tocn", tocn_m)
  call dxa%get("socn", socn_a)
  call dxm%get("socn", socn_m)
  call dxa%get("ssh",  ssh_a)
  call dxm%get("ssh",  ssh_m)
  call dxm%get("cicen",cicen_m)

  do i = self%isc, self%iec
     do j = self%jsc, self%jec
        ! Temperature
        tocn_a%val(i,j,1) = tocn_m%val(i,j,1) + &
             &self%kst%jacobian(i,j,1) * socn_m%val(i,j,1) + &
             &self%ksshts%kssht(i,j,1) * ssh_m%val(i,j,1) +&
             &self%kct(i,j) * sum(cicen_m%val(i,j,2:))
        tocn_a%val(i,j,2:) = tocn_m%val(i,j,2:) + &
             &self%kst%jacobian(i,j,2:) * socn_m%val(i,j,2:) + &
             &self%ksshts%kssht(i,j,2:) * ssh_m%val(i,j,1)
        ! Salinity
        socn_a%val(i,j,:) = socn_m%val(i,j,:) + &
             &self%ksshts%ksshs(i,j,:) * ssh_m%val(i,j, 1)
        ! SSH
        ssh_a%val(i,j,:) = ssh_m%val(i,j,:)
     end do
  end do

  ! copy surface fields 
  ! TODO do this to all fields not explicitly handled above?
  do i=1, size(dxm%fields)
    fld_m => dxm%fields(i)
    call dxa%get(fld_m%name, fld_a)
    select case(fld_m%name)
    case ('sw','lw','lhf','shf','us','cicen','hicen')
      call fld_a%copy(fld_m)
    end select
  end do

end subroutine soca_balance_multad

! ------------------------------------------------------------------------------
! Apply inverse of the forward balance operator
subroutine soca_balance_multinv(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_fields),         intent(in) :: dxm
  type(soca_fields),      intent(inout) :: dxa

  real(kind=kind_real) :: deta
  integer :: i, j, k
  type(soca_field), pointer :: fld_m, fld_a
  type(soca_field), pointer :: tocn_m, tocn_a
  type(soca_field), pointer :: socn_m, socn_a
  type(soca_field), pointer :: ssh_m,  ssh_a
  type(soca_field), pointer :: cicen_m, cicen_a

  call dxm%get("tocn", tocn_m)
  call dxa%get("tocn", tocn_a)
  call dxm%get("socn", socn_m)
  call dxa%get("socn", socn_a)
  call dxm%get("ssh",  ssh_m)
  call dxa%get("ssh",  ssh_a)
  call dxm%get("cicen",cicen_m)
  call dxa%get("cicen",cicen_a)

  do i = self%isc, self%iec
     do j = self%jsc, self%jec
        ! Temperature
        tocn_a%val(i,j,:) = tocn_m%val(i,j,:)
        
        ! Salinity
        socn_a%val(i,j,:) = socn_m%val(i,j,:) -&
             &self%kst%jacobian(i,j,:) * tocn_m%val(i,j,:)
        ! SSH
        deta = 0.0d0
        do k = 1, tocn_a%nz
           deta = deta + ( self%ksshts%ksshs(i,j,k) * self%kst%jacobian(i,j,k) - &
                &self%ksshts%kssht(i,j,k) ) * tocn_m%val(i,j,k) - &
                &self%ksshts%ksshs(i,j,k) * socn_m%val(i,j,k)
        end do
        ssh_a%val(i,j, :) = ssh_m%val(i,j, :) + deta
        ! Ice fraction
        cicen_a%val(i,j,:) =  cicen_m%val(i,j,:)
        do k = 1, cicen_a%nz-1
           cicen_a%val(i,j,k+1) = cicen_a%val(i,j,k+1) -&
                & self%kct(i,j) * tocn_m%val(i,j,1)
        end do
     end do
  end do

  ! copy surface fields 
  ! TODO do this to all fields not explicitly handled above?
  do i=1, size(dxm%fields)
    fld_m => dxm%fields(i)
    call dxa%get(fld_m%name, fld_a)
    select case(fld_m%name)
    case ('sw','lw','lhf','shf','us', 'hicen')
      call fld_a%copy(fld_m)
    end select
  end do
end subroutine soca_balance_multinv

! ------------------------------------------------------------------------------
! Apply inverse of the backward balance operator
subroutine soca_balance_multinvad(self, dxa, dxm)
  type(soca_balance_config), intent(in) :: self
  type(soca_fields),      intent(inout) :: dxm
  type(soca_fields),         intent(in) :: dxa

  integer :: i, j
  type(soca_field), pointer :: fld_a, fld_m
  type(soca_field), pointer :: tocn_a, tocn_m
  type(soca_field), pointer :: socn_a, socn_m
  type(soca_field), pointer :: ssh_a,  ssh_m
  type(soca_field), pointer :: cicen_a  

  call dxm%get("tocn", tocn_m)  
  call dxa%get("tocn", tocn_a)
  call dxm%get("socn", socn_m)  
  call dxa%get("socn", socn_a)
  call dxm%get("ssh",  ssh_m)  
  call dxa%get("ssh",  ssh_a)
  call dxa%get("cicen",cicen_a)

  do i = self%isc, self%iec
     do j = self%jsc, self%jec
        ! Temperature
        tocn_m%val(i,j,1) = tocn_a%val(i,j,1) &
             & - self%kst%jacobian(i,j,1) * socn_a%val(i,j,1) &
             & + ( self%ksshts%ksshs(i,j,1) * self%kst%jacobian(i,j,1) &
             &     - self%ksshts%kssht(i,j,1) ) * ssh_a%val(i,j,1) &
             & - self%kct(i,j) * sum(cicen_a%val(i,j,2:))
        tocn_m%val(i,j,2:) = tocn_a%val(i,j,2:) &
             & - self%kst%jacobian(i,j,2:) * socn_a%val(i,j,2:) &
             & + ( self%ksshts%ksshs(i,j,2:) * self%kst%jacobian(i,j,2:) &
             &     - self%ksshts%kssht(i,j,2:) ) * ssh_a%val(i,j,1)
        ! Salinity
        socn_m%val(i,j,:) = socn_a%val(i,j,:) - &
             &self%ksshts%ksshs(i,j,:) * ssh_a%val(i,j,1)
        ! SSH
        ssh_m%val(i,j,:) = ssh_a%val(i,j,:)
     end do
  end do

  ! copy surface fields 
  ! TODO do this to all fields not explicitly handled above?
  do i=1, size(dxa%fields)
    fld_a => dxa%fields(i)
    call dxm%get(fld_a%name, fld_m)
    select case(fld_a%name)
    case ('sw','lw','lhf','shf','us','hicen','cicen')
      fld_m%val = fld_a%val
    end select
  end do

end subroutine soca_balance_multinvad

end module soca_balance_mod
