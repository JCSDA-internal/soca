!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_balance_mod

  use kinds
  use soca_fields
  use soca_kst_mod
  use soca_ksshts_mod
  
  implicit none

  !> Fortran derived type to hold configuration D
  type :: soca_balance_config
     type(soca_field),  pointer :: traj                !> Trajectory
     integer                    :: isc, iec, jsc, jec  !> Compute domain
     type(soca_kst)             :: kst                 !> T/S balance
     type(soca_ksshts)          :: ksshts              !> SSH/T/S balance
     real(kind=kind_real), allocatable :: kct(:,:)     !> C/T Jacobian     
  end type soca_balance_config

#define LISTED_TYPE soca_balance_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_balance_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------
  subroutine soca_balance_setup(c_conf, self, traj)

    use kinds
    use iso_c_binding
    use config_mod
    use soca_fields
    use soca_model_geom_type, only : geom_get_domain_indices
    use soca_kst_mod
    use datetime_mod
    use mpi
    
    implicit none

    type(soca_balance_config), intent(inout) :: self
    type(soca_field),    target, intent(in)  :: traj
    type(c_ptr),                 intent(in)  :: c_conf

    integer :: isc, iec, jsc, jec, i, j, k, nl
    real(kind=kind_real), allocatable :: jac(:)
    
    ! Number of ocean layer
    nl = size(traj%hocn,3)

    ! Store trajectory
    self%traj => traj

    ! Indices for compute domain
    call geom_get_domain_indices(traj%geom%ocean, "compute", isc, iec, jsc, jec)
    self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec

    ! Get configuration for Kst
    self%kst%dsdtmax      = config_get_real(c_conf,"dsdtmax")
    self%kst%dsdzmin      = config_get_real(c_conf,"dsdzmin")
    self%kst%dtdzmin      = config_get_real(c_conf,"dtdzmin")

    ! Compute and store Jacobian of Kst
    allocate(self%kst%jacobian(isc:iec,jsc:jec,traj%geom%ocean%nzo))
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
       end do
    end do
    deallocate(jac)
    
    ! Compute Jacobian of Ksshts
    allocate(self%ksshts%kssht(isc:iec,jsc:jec,traj%geom%ocean%nzo))
    allocate(self%ksshts%ksshs(isc:iec,jsc:jec,traj%geom%ocean%nzo))
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
                  &traj%geom%ocean%lon(i,j),&
                  &traj%geom%ocean%lat(i,j))
             ! Set Jacobian to 0 in mixed layer
             !if (self%traj%layer_depth(i,j,k)<self%traj%mld(i,j)) then
             !   jac = 0.0_kind_Real
             !end if             
             self%ksshts%kssht(i,j,k) = jac(1)
             self%ksshts%ksshs(i,j,k) = jac(2)
          end do
       end do
    end do
    deallocate(jac)

    ! Compute Kst
    allocate(self%kct(isc:iec,jsc:jec))    
    self%kct = -0.01d0 ! HaHa    
    
  end subroutine soca_balance_setup

  ! ------------------------------------------------------------------------------
  subroutine soca_balance_delete(self)

    use kinds
    use iso_c_binding

    implicit none
    
    type(soca_balance_config), intent(inout) :: self

    nullify(self%traj)
    deallocate(self%kst%jacobian)
    deallocate(self%ksshts%kssht)
    deallocate(self%ksshts%ksshs)    
    deallocate(self%kct)
    
  end subroutine soca_balance_delete

  ! ------------------------------------------------------------------------------
  subroutine soca_balance_mult(self, dxa, dxm)

    !>    [ I       0   0  0 ]
    !>    [ Kst     I   0  0 ]
    !> K= [ Ketat Ketas I  0 ]
    !>    [ Kct     0   0  I ]
    
    use kinds
    use soca_model_geom_type, only : geom_get_domain_indices

    implicit none

    type(soca_balance_config), intent(in) :: self    
    type(soca_field),          intent(in) :: dxa
    type(soca_field),       intent(inout) :: dxm

    integer :: i, j, k
    real(kind=kind_real) :: deta, dxc

    do i = self%isc, self%iec
       do j = self%jsc, self%jec
          ! Temperature
          dxc = sum(dxa%cicen(i,j,2:))
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
          dxm%cicen(i,j,:) =  dxa%cicen(i,j,:)
          do k = 1, size(self%traj%hicen,3)
             dxm%cicen(i,j,k+1) =  dxm%cicen(i,j,k+1) +&
                  & self%kct(i,j) * dxa%tocn(i,j,1)
          end do
          ! Ice thickness
          dxm%hicen(i,j,:) =  dxa%hicen(i,j,:)
       end do
    end do

  end subroutine soca_balance_mult
  
  ! ------------------------------------------------------------------------------
  subroutine soca_balance_multad(self, dxa, dxm)

    use kinds
    use soca_model_geom_type, only : geom_get_domain_indices

    implicit none

    type(soca_balance_config), intent(in) :: self
    type(soca_field),          intent(in) :: dxm    
    type(soca_field),       intent(inout) :: dxa

    real(kind=kind_real) :: dxc
    integer :: i, j, k

    do i = self%isc, self%iec
       do j = self%jsc, self%jec
          ! Temperature          
          dxa%tocn(i,j,1) = dxm%tocn(i,j,1) + &
               &self%kst%jacobian(i,j,1) * dxm%socn(i,j,1) + &
               &self%ksshts%kssht(i,j,1) * dxm%ssh(i,j) +&
               &self%kct(i,j) * sum(dxm%cicen(i,j,2:))
          dxa%tocn(i,j,2:) = dxm%tocn(i,j,2:) + &
               &self%kst%jacobian(i,j,2:) * dxm%socn(i,j,2:) + &
               &self%ksshts%kssht(i,j,2:) * dxm%ssh(i,j)
          ! Salinity
          dxa%socn(i,j,:) = dxm%socn(i,j,:) + &
               &self%ksshts%ksshs(i,j,:) * dxm%ssh(i,j)
          ! SSH
          dxa%ssh(i,j)    = dxm%ssh(i,j)
          ! Ice fraction
          dxa%cicen(i,j,:) =  dxm%cicen(i,j,:)
          ! Ice thickness
          dxa%hicen(i,j,:) =  dxm%hicen(i,j,:)
       end do
    end do

  end subroutine soca_balance_multad

  ! ------------------------------------------------------------------------------  
  subroutine soca_balance_multinv(self, dxa, dxm)

    use kinds
    use soca_model_geom_type, only : geom_get_domain_indices

    implicit none

    type(soca_balance_config), intent(in) :: self
    type(soca_field),          intent(in) :: dxm    
    type(soca_field),       intent(inout) :: dxa

    real(kind=kind_real) :: dxc, deta
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
          dxa%cicen(i,j,:) =  dxm%cicen(i,j,:)
          do k = 1, size(self%traj%hicen,3)
             dxa%cicen(i,j,k+1) =  dxa%cicen(i,j,k+1) -&
                  & self%kct(i,j) * dxm%tocn(i,j,1)
          end do          
          ! Ice thickness
          dxa%hicen(i,j,:) =  dxm%hicen(i,j,:)
       end do
    end do

  end subroutine soca_balance_multinv

  ! ------------------------------------------------------------------------------
  subroutine soca_balance_multinvad(self, dxa, dxm)

    use kinds
    use soca_model_geom_type, only : geom_get_domain_indices

    implicit none

    type(soca_balance_config), intent(in) :: self
    type(soca_field),       intent(inout) :: dxm    
    type(soca_field),          intent(in) :: dxa

    real(kind=kind_real) :: dxc
    integer :: i, j, k

    do i = self%isc, self%iec
       do j = self%jsc, self%jec
          ! Ice thickness
          dxm%hicen(i,j,:) =  dxa%hicen(i,j,:)
          ! Ice fraction
          dxm%cicen(i,j,:) =  dxa%cicen(i,j,:)          
          ! Temperature          
          dxm%tocn(i,j,1) = dxa%tocn(i,j,1) &
               & - self%kst%jacobian(i,j,1) * dxa%socn(i,j,1) &
               & + ( self%ksshts%ksshs(i,j,1) * self%kst%jacobian(i,j,1) &
               &     - self%ksshts%kssht(i,j,1) ) * dxa%ssh(i,j) &
               & - self%kct(i,j) * sum(dxa%cicen(i,j,2:))
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

  end subroutine soca_balance_multinvad
  
end module soca_balance_mod


