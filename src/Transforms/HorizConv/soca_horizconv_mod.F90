!
! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_horizconv_mod
  use config_mod
  use iso_c_binding  
  use kinds
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices
  use soca_bumpcorr2d_mod
  use tools_func
  use type_mpl
  
  implicit none

  !> Fortran derived type to hold the setup for Horizconv
  type :: soca_horizconv
     type(soca_bumpcorr2d_type) :: conv   !> Bump convolution data structure 
     character(len=1024)        :: prefix !> File prefix for bump    
     type(soca_field),  pointer :: bkg    !> Background
  end type soca_horizconv

#define LISTED_TYPE soca_horizconv

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_horizconv_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------
  ! Setup for the horizical convolution
  subroutine soca_conv_setup (self, bkg, c_conf)
    type(soca_horizconv),   intent(inout) :: self
    type(soca_field), target, intent(in) :: bkg
    type(c_ptr),              intent(in) :: c_conf

    integer :: isc, iec, jsc, jec, i, j, k, nl
    character(len=1024)        :: test
    
    nl = size(bkg%hocn,3)
  
    ! Get configuration for horizontal smoothing decorrelation length scale
    self%prefix = config_get_string(c_conf, len(self%prefix), "filter_prefix")

    ! Associate background
    self%bkg => bkg

    ! Initialize bump
    call self%conv%initialize(bkg%geom, c_conf, self%prefix)
    
  end subroutine soca_conv_setup

  ! ------------------------------------------------------------------------------
  !> Apply forward convolution  
  subroutine soca_conv (self, convdx, dx)
    type(soca_horizconv), intent(inout) :: self
    type(soca_field),        intent(in) :: dx
    type(soca_field),     intent(inout) :: convdx

    real(kind=kind_real), allocatable :: incr2d(:,:)
    integer :: isc, iec, jsc, jec
    
    ! Make copy of dx
    call copy(convdx, dx)

    ! Hardcode var to smooth    
    call self%conv%apply(convdx%socn(:,:,1), self%bkg%geom)
    
  end subroutine soca_conv

  ! ------------------------------------------------------------------------------
  !> Apply backward convolution
  subroutine soca_conv_ad (self, convdx, dx)
    type(soca_horizconv), intent(inout) :: self
    type(soca_field),     intent(inout) :: dx     ! OUT
    type(soca_field),        intent(in) :: convdx ! IN

    ! Make copy of convdx
    call copy(dx, convdx)
    !call self%conv%apply(convdx, self%bkg%geom)
    ! Hardcode var to smooth    
    call self%conv%applyad(dx%socn(:,:,1), self%bkg%geom)
    
  end subroutine soca_conv_ad

end module soca_horizconv_mod
