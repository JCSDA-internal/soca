! (C) Copyright 2017- UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Structure holding configuration variables for the 3d error
!! horizfilt matrices of the SOCA analysis.

module soca_horizfilt_mod
  use config_mod
  use iso_c_binding
  use kinds
  use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad
  use oobump_mod, only: bump_read_conf
  use soca_fields
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
     real(kind=kind_real), allocatable :: wgh(:,:)       !< Filtering weight
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
    integer :: isc, iec, jsc, jec, i, j, ivar
    logical :: init_seaice, init_ocean
    real(kind=kind_real) :: sum_wgh

    ! Setup list of variables to apply filtering on
    self%vars = vars

    ! Indices for compute/data domain
    self%isd = geom%isd ;  self%ied = geom%ied ; self%jsd = geom%jsd; self%jed = geom%jed
    self%isc = geom%isc ;  self%iec = geom%iec ; self%jsc = geom%jsc; self%jec = geom%jec

    ! Allocate and compute filtering weights
    allocate(self%wgh(self%isd:self%ied,self%jsd:self%jed))

    ! TODO (Guillaume): Compute distance based weights
    self%wgh = 1.0_kind_real

    ! Normalize
    do i = self%isc, self%iec
       do j = self%jsc, self%jec
          sum_wgh = &
               &self%wgh(i-1,j+1) + self%wgh(i,j+1) + self%wgh(i+1,j+1) + &
               &self%wgh(i-1,j)   + self%wgh(i,j)   + self%wgh(i+1,j)   + &
               &self%wgh(i-1,j-1) + self%wgh(i,j-1) + self%wgh(i+1,j-1)
          self%wgh(i,j) = self%wgh(i,j)/sum_wgh
       end do
    end do

    self%wgh = 1.0_kind_real
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

    !real(kind=kind_real), allocatable :: dx2d(:,:)

    !allocate(dx2d(self%isd:self%ied,self%jsd:self%jed))
    !dx2d = dxin%ssh
    call soca_filt2d(self, dxin%ssh, dxout%ssh, geom)
    !dxout%ssh = dx2d
    !deallocate(dx2d)

  end subroutine soca_horizfilt_mult

  ! ------------------------------------------------------------------------------

  subroutine soca_horizfilt_multad(self, dxin, dxout, geom)
    class(soca_horizfilt_type), intent(inout) :: self  !< The horizfilt structure
    type(soca_field),              intent(in) :: dxin  !< Input: Increment
    type(soca_field),           intent(inout) :: dxout !< Output: filtered Increment
    type(soca_geom),               intent(in) :: geom

    real(kind=kind_real), allocatable :: dx2d(:,:)

    !allocate(dx2d(self%isd:self%ied,self%jsd:self%jed))
    !dx2d = dxin%ssh
    call soca_filt2d_ad(self, dxin%ssh, dxout%ssh, geom)
    !call soca_filt2d(self, dxin%ssh, dxout%ssh, geom)
    !dxout%ssh = dx2d
    !deallocate(dx2d)

  end subroutine soca_horizfilt_multad

  ! ------------------------------------------------------------------------------

  subroutine soca_filt2d(self, dxin, dxout, geom)
    class(soca_horizfilt_type), intent(in) :: self
    real(kind=kind_real),       intent(in) :: dxin(:,:)
    real(kind=kind_real),    intent(inout) :: dxout(:,:)
    type(soca_geom),            intent(in) :: geom

    integer :: i, j
    real(kind=kind_real), allocatable :: dxs(:,:)

    !allocate(dxs(self%isd:self%ied,self%jsd:self%jed))
    !dxs = dx

    ! Fill halo
    !call mpp_update_domains(dxout, geom%G%Domain%mpp_domain)

    ! 9-point weighted average
    do i = self%isc, self%iec
       do j = self%jsc, self%jec
          dxout(i,j) = &
           dxin(i-1,j+1) + dxin(i,j+1) + dxin(i+1,j+1) + &
           dxin(i-1,j)   + dxin(i,j)   + dxin(i+1,j)     + &
           dxin(i-1,j-1) + dxin(i,j-1) + dxin(i+1,j-1)
       end do
    end do

    ! Update halo
    !call mpp_update_domains(dxout, geom%G%Domain%mpp_domain)

  end subroutine soca_filt2d

  ! ------------------------------------------------------------------------------

  subroutine soca_filt2d_ad(self, dxin, dxout, geom)
    class(soca_horizfilt_type), intent(in) :: self
    real(kind=kind_real),       intent(in) :: dxin(:,:)
    real(kind=kind_real),    intent(inout) :: dxout(:,:)
    type(soca_geom),            intent(in) :: geom

    integer :: i, j
    real(kind=kind_real), allocatable :: dxs(:,:)

    allocate(dxs(self%isd:self%ied,self%jsd:self%jed))
    dxs = dxin

    ! Update inner grid from halo
    !call mpp_update_domains_ad(dxout, geom%G%Domain%mpp_domain)

    dxout = 0.0_kind_real
    ! Adjoint of 9-point weighted average
    do i =  self%iec-1, self%isc+1, -1
       do j = self%jec-1, self%jsc+1, -1
          dxout(i-i,j+1) = dxout(i-i,j+1) + dxs(i,j)
          dxout(i,j+1)   = dxout(i,j+1)   + dxs(i,j)
          dxout(i+1,j+1) = dxout(i+1,j+1) + dxs(i,j)
          dxout(i-1,j)   = dxout(i-1,j)   + dxs(i,j)
          dxout(i,j)     = dxout(i,j)     + dxs(i,j)
          dxout(i+1,j)   = dxout(i+1,j)   + dxs(i,j)
          dxout(i-1,j-1) = dxout(i-1,j-1) + dxs(i,j)
          dxout(i,j-1)   = dxout(i,j-1)   + dxs(i,j)
          dxout(i+1,j-1) = dxout(i+1,j-1) + dxs(i,j)
          dxs(i,j) = 0.0_kind_real
       end do
    end do


  end subroutine soca_filt2d_ad
end module soca_horizfilt_mod
