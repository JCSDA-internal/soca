! (C) Copyright 2017-2021 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> horizontal filtering
module soca_horizfilt_mod

use atlas_module, only: atlas_geometry
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad
use oops_variables_mod, only: oops_variables

! soca modules
use soca_fields_mod, only: soca_field
use soca_geom_mod, only : soca_geom
use soca_increment_mod, only: soca_increment
use soca_state_mod, only: soca_state

implicit none
private

!> Variable transform: horizontal filtering
type, public :: soca_horizfilt
  type(oops_variables)              :: vars           !< Apply filtering to vars
  real(kind=kind_real), allocatable :: wgh(:,:,:,:)   !< Filtering weight
  real(kind=kind_real) :: scale_flow  !< Used with "flow" filter, sea surface height decorrelation scale
  real(kind=kind_real) :: scale_dist
  real(kind=kind_real) :: niter !< number of iterations of filter to perform
  !> \name indices of compute domain
  !! \{
  integer :: isc, iec, jsc, jec
  !> \}

  !> \name indices of data domain
  !! \{
  integer :: isd, ied, jsd, jed
  !> \}

contains

  !> \copybrief soca_horizfilt_setup \see soca_horizfilt_setup
  procedure :: setup => soca_horizfilt_setup

  !> \copybrief soca_horizfilt_delete \see soca_horizfilt_delete
  procedure :: delete => soca_horizfilt_delete

  !> \copybrief soca_horizfilt_mult \see soca_horizfilt_mult
  procedure :: mult => soca_horizfilt_mult

  !> \copybrief soca_horizfilt_multad \see soca_horizfilt_multad
  procedure :: multad => soca_horizfilt_multad
end type soca_horizfilt


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Setup for the horizfilt operator
!!
!! \relates soca_horizfilt_mod::soca_horizfilt
subroutine soca_horizfilt_setup(self, f_conf, geom, traj, vars)
  class(soca_horizfilt), intent(inout) :: self   !< The horizfilt structure
  type(fckit_configuration),     intent(in) :: f_conf !< The configuration
  type(soca_geom),               intent(in) :: geom   !< Geometry
  type(soca_state),              intent(in) :: traj   !< Trajectory
  type(oops_variables),          intent(in) :: vars   !< List of variables

  type(soca_field), pointer :: ssh

  integer :: i, j, ii, jj
  real(kind=kind_real) :: dist(-1:1,-1:1), sum_w, r_dist, r_flow
  type(atlas_geometry) :: ageometry

  ! Setup list of variables to apply filtering on
  self%vars = vars

  ! Get filter length scales from configuration
  if (.not. f_conf%get("scale_dist", self%scale_dist)) self%scale_dist = -1
  if (.not. f_conf%get("scale_flow", self%scale_flow)) self%scale_flow = -1
  call f_conf%get_or_die("niter", self%niter)

  ! scale the gaussian length scales depending on the number of iterations
  ! NOTE: numerical instability creeps in once niter is large and scale_dist
  !  is much smaller than the size of a grid box
  self%scale_dist = self%scale_dist / sqrt(self%niter)
  self%scale_flow = self%scale_flow / sqrt(self%niter)

  ! Indices for compute/data domain
  self%isd = geom%isd ;  self%ied = geom%ied ; self%jsd = geom%jsd; self%jed = geom%jed
  self%isc = geom%isc ;  self%iec = geom%iec ; self%jsc = geom%jsc; self%jec = geom%jec

  ! Allocate and compute filtering weights
  allocate(self%wgh(self%isd:self%ied,self%jsd:self%jed,-1:1,-1:1))

  ! Create UnitSphere geometry
  ageometry = atlas_geometry("Earth")

  ! Compute distance based weights
  self%wgh = 0.0_kind_real
  r_dist = 1.0
  r_flow = 1.0
  do j = self%jsc, self%jec
      do i = self%isc, self%iec
        do ii = -1,1
            do jj = -1,1
              ! Great circle distance
              if(self%scale_dist > 0) then
                r_dist = ageometry%distance(geom%lon(i,j), geom%lat(i,j), geom%lon(i+ii,j+jj), geom%lat(i+ii,j+jj) )
                r_dist = exp(-0.5 * (r_dist/self%scale_dist) ** 2)
              end if

              ! flow based distance (ssh difference)
              if(self%scale_flow > 0) then
                call traj%get("ssh", ssh)
                r_flow = abs(ssh%val(i,j,1) - ssh%val(i+ii,j+jj,1))
                r_flow = exp(-0.5 * ((r_flow / self%scale_flow) ** 2))
              end if

              ! multiply together and apply the land mask
              dist(ii,jj) = geom%mask2d(i+ii,j+jj) * r_dist * r_flow
            end do
        end do

        ! Normalize
        sum_w = sum(dist)
        if (sum_w>0.0_kind_real) then
          self%wgh(i,j,:,:) = dist / sum_w
        endif

      end do
  end do

end subroutine soca_horizfilt_setup


! ------------------------------------------------------------------------------
!> Delete horizfilt
!!
!! \relates soca_horizfilt_mod::soca_horizfilt
subroutine soca_horizfilt_delete(self)
  class(soca_horizfilt), intent(inout) :: self       !< The horizfilt structure

  deallocate(self%wgh)

end subroutine soca_horizfilt_delete


! ------------------------------------------------------------------------------
!> Forward filtering
!!
!! \relates soca_horizfilt_mod::soca_horizfilt
subroutine soca_horizfilt_mult(self, dxin, dxout, geom)
  class(soca_horizfilt), intent(inout) :: self  !< The horizfilt structure
  type(soca_increment),          intent(in) :: dxin  !< Input: Increment
  type(soca_increment),       intent(inout) :: dxout !< Output: filtered Increment
  type(soca_geom),               intent(in) :: geom

  type(soca_field), pointer :: field_i, field_o

  integer :: k, ivar
  real(kind=kind_real), allocatable, dimension(:,:) :: dxi, dxo

  allocate(dxi(self%isd:self%ied,self%jsd:self%jed))
  allocate(dxo(self%isd:self%ied,self%jsd:self%jed))

  do ivar = 1, self%vars%nvars()
    call dxin%get(trim(self%vars%variable(ivar)),  field_i)
    call dxout%get(trim(self%vars%variable(ivar)), field_o)
    do k = 1, field_i%nz
      dxi = field_i%val(:,:,k)
      call soca_filt2d(self, dxi, dxo, geom)
      field_o%val(:,:,k) = dxo
    end do
  end do
  deallocate(dxi, dxo)

end subroutine soca_horizfilt_mult


! ------------------------------------------------------------------------------
!> Backward filtering
!!
!! \relates soca_horizfilt_mod::soca_horizfilt
subroutine soca_horizfilt_multad(self, dxin, dxout, geom)
  class(soca_horizfilt), intent(inout) :: self  !< The horizfilt structure
  type(soca_increment),          intent(in) :: dxin  !< Input:
  type(soca_increment),       intent(inout) :: dxout !< Output:
  type(soca_geom),               intent(in) :: geom

  type(soca_field), pointer :: field_i, field_o
  integer :: k, ivar
  real(kind=kind_real), allocatable, dimension(:,:) :: dxi, dxo

  allocate(dxi(self%isd:self%ied,self%jsd:self%jed))
  allocate(dxo(self%isd:self%ied,self%jsd:self%jed))

  do ivar = 1, self%vars%nvars()
    call dxin%get(trim(self%vars%variable(ivar)),  field_i)
    call dxout%get(trim(self%vars%variable(ivar)), field_o)
      do k = 1, field_i%nz
        dxi = field_i%val(:,:,k)
        call soca_filt2d_ad(self, dxi, dxo, geom)
        field_o%val(:,:,k) = dxo
      end do
  end do
  deallocate(dxi, dxo)

end subroutine soca_horizfilt_multad


! ------------------------------------------------------------------------------
!> Forward filtering for 2D array
!!
!! used by soca_horizfilt_mod::soca_horizfilt::mult()
!! \relates soca_horizfilt_mod::soca_horizfilt
subroutine soca_filt2d(self, dxin, dxout, geom)
  class(soca_horizfilt),           intent(in) :: self
  real(kind=kind_real),    allocatable, intent(in) :: dxin(:,:)
  real(kind=kind_real), allocatable, intent(inout) :: dxout(:,:)
  type(soca_geom),                      intent(in) :: geom

  integer :: i, j
  real(kind=kind_real), allocatable :: dxtmp(:,:)

  ! Make a temporary copy of dxin
  allocate(dxtmp(self%isd:self%ied,self%jsd:self%jed))
  dxtmp = 0.0_kind_real
  dxtmp(self%isc:self%iec,self%jsc:self%jec) = dxin(self%isc:self%iec,self%jsc:self%jec)

  ! Update halo points of input array
  call mpp_update_domains(dxtmp, geom%Domain%mpp_domain, complete=.true.)

  ! 9-point distance weighted average
  dxout = 0.0_kind_real
  do j = self%jsc, self%jec
      do i = self%isc, self%iec
        if (geom%mask2d(i,j)==1) then
            dxout(i,j) = &
              self%wgh(i,j,-1,1)*dxtmp(i-1,j+1) + &
              self%wgh(i,j,0,1)*dxtmp(i,j+1) + &
              self%wgh(i,j,1,1)*dxtmp(i+1,j+1) + &
              self%wgh(i,j,-1,0)*dxtmp(i-1,j) + &
              self%wgh(i,j,0,0)*dxtmp(i,j) + &
              self%wgh(i,j,1,0)*dxtmp(i+1,j) + &
              self%wgh(i,j,-1,-1)*dxtmp(i-1,j-1) + &
              self%wgh(i,j,0,-1)*dxtmp(i,j-1) + &
              self%wgh(i,j,1,-1)*dxtmp(i+1,j-1)
        end if
      end do
  end do

  ! Update halo
  call mpp_update_domains(dxout, geom%Domain%mpp_domain, complete=.true.)

  deallocate(dxtmp)

end subroutine soca_filt2d


! ------------------------------------------------------------------------------
!> Backward filtering for 2D array
!!
!! used by soca_horizfilt_mod::soca_horizfilt::multad()
!! \relates soca_horizfilt_mod::soca_horizfilt
subroutine soca_filt2d_ad(self, dxin, dxout, geom)
  class(soca_horizfilt),           intent(in) :: self
  real(kind=kind_real),    allocatable, intent(in) :: dxin(:,:)
  real(kind=kind_real), allocatable, intent(inout) :: dxout(:,:)
  type(soca_geom),                      intent(in) :: geom

  integer :: i, j
  real(kind=kind_real), allocatable :: dxtmp(:,:)

  ! Make a temporary copy of dxin
  allocate(dxtmp(self%isd:self%ied,self%jsd:self%jed))
  dxtmp = 0.0_kind_real
  dxtmp(self%isc:self%iec,self%jsc:self%jec) = dxin(self%isc:self%iec,self%jsc:self%jec)

  dxout = 0.0_kind_real
  ! Adjoint of 9-point weighted average
  do j = self%jec, self%jsc, -1
      do i =  self%iec, self%isc, -1
        if (geom%mask2d(i,j)==1) then
            dxout(i-1,j+1) = dxout(i-1,j+1) + self%wgh(i,j,-1, 1)*dxtmp(i,j)
            dxout(i,j+1)   = dxout(i,j+1)   + self%wgh(i,j, 0, 1)*dxtmp(i,j)
            dxout(i+1,j+1) = dxout(i+1,j+1) + self%wgh(i,j, 1, 1)*dxtmp(i,j)
            dxout(i-1,j)   = dxout(i-1,j)   + self%wgh(i,j,-1, 0)*dxtmp(i,j)
            dxout(i,j)     = dxout(i,j)     + self%wgh(i,j, 0, 0)*dxtmp(i,j)
            dxout(i+1,j)   = dxout(i+1,j)   + self%wgh(i,j, 1, 0)*dxtmp(i,j)
            dxout(i-1,j-1) = dxout(i-1,j-1) + self%wgh(i,j,-1,-1)*dxtmp(i,j)
            dxout(i,j-1)   = dxout(i,j-1)   + self%wgh(i,j, 0,-1)*dxtmp(i,j)
            dxout(i+1,j-1) = dxout(i+1,j-1) + self%wgh(i,j, 1,-1)*dxtmp(i,j)
        end if
      end do
  end do

  ! Adjoint of halo update
  call mpp_update_domains_ad(dxout, geom%Domain%mpp_domain, complete=.true.)

  deallocate(dxtmp)

end subroutine soca_filt2d_ad

end module soca_horizfilt_mod
