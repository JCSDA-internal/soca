!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bkgerr_mod
  use config_mod
  use datetime_mod
  use iso_c_binding
  use kinds
  use soca_bkgerrutil_mod
  use soca_fields
  use soca_utils
  use soca_omb_stats_mod
  use fckit_mpi_module

  implicit none

  !> Fortran derived type to hold configuration D
  type :: soca_bkgerr_config
     type(soca_field),         pointer :: bkg
     type(soca_field)                  :: std_bkgerr
     type(soca_bkgerr_bounds_type)     :: bounds         ! Bounds for bkgerr
     real(kind=kind_real)              :: std_sst
     real(kind=kind_real)              :: std_sss
     integer                           :: isc, iec, jsc, jec
  end type soca_bkgerr_config

#define LISTED_TYPE soca_bkgerr_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_bkgerr_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------
  !> Setup the static background error
  subroutine soca_bkgerr_setup(c_conf, self, bkg)
    type(soca_bkgerr_config), intent(inout) :: self
    type(soca_field),    target, intent(in) :: bkg
    type(c_ptr),                 intent(in) :: c_conf

    integer :: isc, iec, jsc, jec, i, j, k, nl
    type(datetime) :: vdate
    character(len=800) :: fname = 'soca_bkgerrsoca.nc'
    logical :: read_from_file = .false.

    ! Get number of ocean levels
    nl = size(bkg%hocn,3)

    ! Allocate memory for bkgerror
    call create_copy(self%std_bkgerr, bkg)

    ! Read variance
    ! Precomputed from an ensemble of (K^-1 dx)
    call read_file(self%std_bkgerr, c_conf, vdate)

    ! Convert to standard deviation
    self%std_bkgerr%tocn = sqrt(self%std_bkgerr%tocn)
    self%std_bkgerr%socn = sqrt(self%std_bkgerr%socn)
    self%std_bkgerr%ssh = sqrt(self%std_bkgerr%ssh)

    ! Get bounds from configuration
    call self%bounds%read(c_conf)

    ! Get constand background error for sst and sss
    if (config_element_exists(c_conf,"fixed_std_sst")) then
       self%std_sst = config_get_real(c_conf,"fixed_std_sst")
       self%std_bkgerr%tocn(:,:,1) = self%std_sst
    end if
    if (config_element_exists(c_conf,"fixed_std_sss")) then
       self%std_sss = config_get_real(c_conf,"fixed_std_sss")
       self%std_bkgerr%socn(:,:,1) = self%std_sss
    end if

    ! Invent background error for ocnsfc fields: set it
    ! to 10% of the background for now ...
    ! TODO: Read background error for ocnsfc from file
    call self%std_bkgerr%ocnsfc%copy(bkg%ocnsfc)
    call self%std_bkgerr%ocnsfc%abs()
    call self%std_bkgerr%ocnsfc%mul(0.1_kind_real)

    ! Associate background
    self%bkg => bkg

    ! Indices for compute domain (no halo)
    isc=bkg%geom%isc; iec=bkg%geom%iec
    jsc=bkg%geom%jsc; jec=bkg%geom%jec

    self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec

    ! Apply config bounds to background error
    call self%bounds%apply(self%std_bkgerr)

    ! Save filtered background error
    call soca_fld2file(self%std_bkgerr, fname)

  end subroutine soca_bkgerr_setup

  ! ------------------------------------------------------------------------------
  !> Apply background error: dxm = D dxa
  subroutine soca_bkgerr_mult(self, dxa, dxm)
    type(soca_bkgerr_config),    intent(in) :: self
    type(soca_field),            intent(in) :: dxa
    type(soca_field),         intent(inout) :: dxm

    integer :: isc, iec, jsc, jec, i, j, k

    ! Indices for compute domain (no halo)
    isc=self%bkg%geom%isc; iec=self%bkg%geom%iec
    jsc=self%bkg%geom%jsc; jec=self%bkg%geom%jec

    do i = isc, iec
       do j = jsc, jec
          dxm%ssh(i,j) = self%std_bkgerr%ssh(i,j) * dxa%ssh(i,j)
          dxm%tocn(i,j,:) = self%std_bkgerr%tocn(i,j,:) * dxa%tocn(i,j,:)
          dxm%socn(i,j,:) = self%std_bkgerr%socn(i,j,:) * dxa%socn(i,j,:)
       end do
    end do
    ! Surface fields
    call dxm%ocnsfc%copy(dxa%ocnsfc)
    call dxm%ocnsfc%schur(self%std_bkgerr%ocnsfc)

    ! Sea-ice
    call dxm%seaice%copy(dxa%seaice)
    call dxm%seaice%schur(self%std_bkgerr%seaice)

  end subroutine soca_bkgerr_mult

  ! ------------------------------------------------------------------------------

end module soca_bkgerr_mod


