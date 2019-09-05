! (C) Copyright 2017- UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Structure holding configuration variables for the 3d error
!! horizfilt matrices of the SOCA analysis.

module soca_horizfilt_mod
  use config_mod
  use fckit_geometry_module, only: sphere_distance
  use iso_c_binding
  use kinds
  use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad
  use oobump_mod, only: bump_read_conf
  use soca_fields_mod
  use soca_geom_mod_c
  use soca_geom_mod, only : soca_geom
  use soca_utils
  use type_bump
  use type_nam
  use random_mod
  use variables_mod

  implicit none

  private
  public :: soca_horizfilt_setup, soca_horizfilt_delete
  public :: soca_horizfilt_mult, soca_horizfilt_multad

  !> Fortran derived type to hold configuration data for horizfilt
  type, public :: soca_horizfilt_type
     type(soca_field),         pointer :: bkg            !< Background field (or first guess)
     type(oops_vars)                   :: vars           !< Apply filtering to vars
     real(kind=kind_real), allocatable :: wgh(:,:,:,:)     !< Filtering weight
     integer :: isc, iec, jsc, jec  !< indices of compute domain
     integer :: isd, ied, jsd, jed  !< indices of data domain
   contains
     procedure :: setup => soca_horizfilt_setup
     procedure :: delete => soca_horizfilt_delete
     procedure :: mult => soca_horizfilt_mult
     procedure :: multad => soca_horizfilt_multad
  end type soca_horizfilt_type

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------

  !> Setup for the horizfilt operator

  subroutine soca_horizfilt_setup(self, c_conf, geom, vars)
    class(soca_horizfilt_type), intent(inout) :: self   !< The horizfilt structure
    type(c_ptr),                   intent(in) :: c_conf !< The configuration
    type(soca_geom),               intent(in) :: geom   !< Geometry
    type(oops_vars),               intent(in) :: vars   !< List of variables

    character(len=3)  :: domain
    integer :: isc, iec, jsc, jec, i, j, ivar, ii, jj
    logical :: init_seaice, init_ocean
    real(kind=kind_real) :: sum_dist, dist(-1:1,-1:1), sum_w

    ! Setup list of variables to apply filtering on
    self%vars = vars

    ! Indices for compute/data domain
    self%isd = geom%isd ;  self%ied = geom%ied ; self%jsd = geom%jsd; self%jed = geom%jed
    self%isc = geom%isc ;  self%iec = geom%iec ; self%jsc = geom%jsc; self%jec = geom%jec

    ! Allocate and compute filtering weights
    allocate(self%wgh(self%isd:self%ied,self%jsd:self%jed,-1:1,-1:1))

    ! Compute distance based weights
    self%wgh = 0.0_kind_real
    do i = self%isc, self%iec
       do j = self%jsc, self%jec
          ! Great circle distance
          do ii = -1,1
             do jj = -1,1
                dist(ii,jj) = sphere_distance(geom%lon(i,j), geom%lat(i,j), geom%lon(i+ii,j+jj), geom%lat(i+ii,j+jj) )
             end do
          end do

          ! Weighted to distance
          sum_dist=sum(dist)
          do ii = -1,1
             do jj = -1,1
                self%wgh(i,j,ii,jj) = sum_dist - dist(ii,jj)
             end do
          end do

          ! Normalize
          sum_w = sum(self%wgh(i,j,:,:))
          do ii = -1,1
             do jj = -1,1
                self%wgh(i,j,ii,jj) = self%wgh(i,j,ii,jj) / sum_w
             end do
          end do

       end do
    end do

  end subroutine soca_horizfilt_setup

  ! ------------------------------------------------------------------------------

  !> Delete horizfilt

  subroutine soca_horizfilt_delete(self)
    class(soca_horizfilt_type), intent(inout) :: self       !< The horizfilt structure

    deallocate(self%wgh)
    nullify(self%bkg)

  end subroutine soca_horizfilt_delete

  ! ------------------------------------------------------------------------------

  subroutine soca_horizfilt_mult(self, dxin, dxout, geom)
    class(soca_horizfilt_type), intent(inout) :: self  !< The horizfilt structure
    type(soca_field),              intent(in) :: dxin  !< Input: Increment
    type(soca_field),           intent(inout) :: dxout !< Output: filtered Increment
    type(soca_geom),               intent(in) :: geom

    integer :: k, ivar
    real(kind=kind_real), allocatable, dimension(:,:) :: dxi, dxo

    allocate(dxi(self%isd:self%ied,self%jsd:self%jed))
    allocate(dxo(self%isd:self%ied,self%jsd:self%jed))

    do ivar = 1, self%vars%nv
       select case (trim(self%vars%fldnames(ivar)))

       case ("ssh")
          call soca_filt2d(self, dxin%ssh, dxout%ssh, geom)

       case ("tocn")
          do k = 1, geom%nzo
             dxi = dxin%tocn(:,:,k)
             call soca_filt2d(self, dxi, dxo, geom)
             dxout%tocn(:,:,k) = dxo
          end do

       case ("socn")
          do k = 1, geom%nzo
             dxi = dxin%socn(:,:,k)
             call soca_filt2d(self, dxi, dxo, geom)
             dxout%socn(:,:,k) = dxo
          end do

       case ("cicen")
          do k = 1, geom%ncat
             dxi = dxin%seaice%cicen(:,:,k)
             call soca_filt2d(self, dxi, dxo, geom)
             dxout%seaice%cicen(:,:,k) = dxo
          end do

       case ("hicen")
          do k = 1, geom%ncat
             dxi = dxin%seaice%hicen(:,:,k)
             call soca_filt2d(self, dxi, dxo, geom)
             dxout%seaice%hicen(:,:,k) = dxo
          end do

       end select
    end do
    deallocate(dxi, dxo)

  end subroutine soca_horizfilt_mult

  ! ------------------------------------------------------------------------------

  subroutine soca_horizfilt_multad(self, dxin, dxout, geom)
    class(soca_horizfilt_type), intent(inout) :: self  !< The horizfilt structure
    type(soca_field),              intent(in) :: dxin  !< Input:
    type(soca_field),           intent(inout) :: dxout !< Output:
    type(soca_geom),               intent(in) :: geom

    integer :: k, ivar
    real(kind=kind_real), allocatable, dimension(:,:) :: dxi, dxo

    allocate(dxi(self%isd:self%ied,self%jsd:self%jed))
    allocate(dxo(self%isd:self%ied,self%jsd:self%jed))

    do ivar = 1, self%vars%nv
       select case (trim(self%vars%fldnames(ivar)))

       case ("ssh")
          call soca_filt2d(self, dxin%ssh, dxout%ssh, geom)

       case ("tocn")
          do k = 1, geom%nzo
             dxi = dxin%tocn(:,:,k)
             call soca_filt2d_ad(self, dxi, dxo, geom)
             dxout%tocn(:,:,k) = dxo
          end do

       case ("socn")
          do k = 1, geom%nzo
             dxi = dxin%socn(:,:,k)
             call soca_filt2d_ad(self, dxi, dxo, geom)
             dxout%socn(:,:,k) = dxo
          end do

       case ("cicen")
          do k = 1, geom%ncat
             dxi = dxin%seaice%cicen(:,:,k)
             call soca_filt2d_ad(self, dxi, dxo, geom)
             dxout%seaice%cicen(:,:,k) = dxo
          end do

       case ("hicen")
          do k = 1, geom%ncat
             dxi = dxin%seaice%hicen(:,:,k)
             call soca_filt2d_ad(self, dxi, dxo, geom)
             dxout%seaice%hicen(:,:,k) = dxo
          end do

       end select
    end do
    deallocate(dxi, dxo)

  end subroutine soca_horizfilt_multad

  ! ------------------------------------------------------------------------------

  subroutine soca_filt2d(self, dxin, dxout, geom)
    class(soca_horizfilt_type),           intent(in) :: self
    real(kind=kind_real),    allocatable, intent(in) :: dxin(:,:)
    real(kind=kind_real), allocatable, intent(inout) :: dxout(:,:)
    type(soca_geom),                      intent(in) :: geom

    integer :: i, j
    real(kind=kind_real), allocatable :: dxtmp(:,:)

    ! Make a temporary copy of dxin
    allocate(dxtmp(self%isd:self%ied,self%jsd:self%jed))
    dxtmp = dxin

    ! Update halo points of input array
    dxtmp=0.0
    dxtmp(self%isc:self%iec,self%jsc:self%jec)=1.0
    print *,'before dxtmp:',sum(dxtmp)
    call mpp_update_domains(dxtmp, geom%Domain%mpp_domain, complete=.true.)
    print *,'after dxtmp:',sum(dxtmp)
    ! 9-point weighted average
    dxout = 0.0_kind_real
    do i = self%isc, self%iec
       do j = self%jsc, self%jec
          if (geom%mask2d(i,j)==1) then
             dxout(i,j) = &
               self%wgh(i,j,-1,1)*geom%mask2d(i-1,j+1)*dxtmp(i-1,j+1) + &
               self%wgh(i,j,0,1)*geom%mask2d(i,j+1)*dxtmp(i,j+1) + &
               self%wgh(i,j,1,1)*geom%mask2d(i+1,j+1)*dxtmp(i+1,j+1) + &
               self%wgh(i,j,-1,0)*geom%mask2d(i-1,j)*dxtmp(i-1,j) + &
               self%wgh(i,j,0,0)*geom%mask2d(i,j)*dxtmp(i,j) + &
               self%wgh(i,j,1,0)*geom%mask2d(i+1,j)*dxtmp(i+1,j) + &
               self%wgh(i,j,-1,-1)*geom%mask2d(i-1,j-1)*dxtmp(i-1,j-1) + &
               self%wgh(i,j,0,-1)*geom%mask2d(i,j-1)*dxtmp(i,j-1) + &
               self%wgh(i,j,1,-1)*geom%mask2d(i+1,j-1)*dxtmp(i+1,j-1)
          end if
       end do
    end do

    !dxout = 0.0
    !dxout(self%isc:self%iec,self%jsc:self%jec)=1.0
    print *,'before:',sum(dxout), sum(dxin)
    ! Update halo
    call mpp_update_domains(dxout, geom%Domain%mpp_domain, complete=.true.)
    print *,'after:',sum(dxout), sum(dxin)
    !print *,'dxout=',dxout(:1,1:10)
    !print *,'dxin=',dxin(:1,1:10)
    !print *,'dxtmp=',dxtmp(:1,1:10)

  end subroutine soca_filt2d

  ! ------------------------------------------------------------------------------

  subroutine soca_filt2d_ad(self, dxin, dxout, geom)
    class(soca_horizfilt_type),           intent(in) :: self
    real(kind=kind_real),    allocatable, intent(in) :: dxin(:,:)
    real(kind=kind_real), allocatable, intent(inout) :: dxout(:,:)
    type(soca_geom),                      intent(in) :: geom

    integer :: i, j

    dxout = 0.0_kind_real
    ! Adjoint of 9-point weighted average
    do i =  self%iec, self%isc, -1
       do j = self%jec, self%jsc, -1
          if (geom%mask2d(i,j)==1) then
             dxout(i-1,j+1) = dxout(i-1,j+1) + self%wgh(i,j,-1, 1)*geom%mask2d(i-1,j+1)*dxin(i,j)
             dxout(i,j+1)   = dxout(i,j+1)   + self%wgh(i,j, 0, 1)*geom%mask2d(i,j+1)*dxin(i,j)
             dxout(i+1,j+1) = dxout(i+1,j+1) + self%wgh(i,j, 1, 1)*geom%mask2d(i+1,j+1)*dxin(i,j)
             dxout(i-1,j)   = dxout(i-1,j)   + self%wgh(i,j,-1, 0)*geom%mask2d(i-1,j)*dxin(i,j)
             dxout(i,j)     = dxout(i,j)     + self%wgh(i,j, 0, 0)*geom%mask2d(i,j)*dxin(i,j)
             dxout(i+1,j)   = dxout(i+1,j)   + self%wgh(i,j, 1, 0)*geom%mask2d(i+1,j)*dxin(i,j)
             dxout(i-1,j-1) = dxout(i-1,j-1) + self%wgh(i,j,-1,-1)*geom%mask2d(i-1,j-1)*dxin(i,j)
             dxout(i,j-1)   = dxout(i,j-1)   + self%wgh(i,j, 0,-1)*geom%mask2d(i,j-1)*dxin(i,j)
             dxout(i+1,j-1) = dxout(i+1,j-1) + self%wgh(i,j, 1,-1)*geom%mask2d(i+1,j-1)*dxin(i,j)
          end if
       end do
    end do

    ! Update inner grid from halo
    !call mpp_update_domains(dxout, geom%Domain%mpp_domain, complete=.true.)

    !call mpp_update_domains_ad(dxout, geom%Domain%mpp_domain)

  end subroutine soca_filt2d_ad
end module soca_horizfilt_mod
