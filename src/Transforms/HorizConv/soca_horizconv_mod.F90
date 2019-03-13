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
  use tools_func
  use type_mpl
  
  implicit none

  !> Fortran derived type to hold the setup for Horizconv
  type :: soca_horizconv
     real(kind=kind_real)      :: lz                 !> Horizical decorrelation [m]
     !type(soca_field),pointer  :: traj               !> Trajectory
     type(soca_field), pointer :: bkg                !> Background     
     integer                   :: isc, iec, jsc, jec !> Compute domain 
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
    !type(soca_field), target, intent(in) :: traj
    type(c_ptr),              intent(in) :: c_conf

    integer :: isc, iec, jsc, jec, i, j, k, nl
    
    nl = size(bkg%hocn,3)
  
    ! Get configuration for horizical convolution
    self%lz      = config_get_real(c_conf, "Lz")

    ! Store trajectory and background
    !self%traj => traj
    self%bkg => bkg
    
    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(bkg%geom%ocean, "compute", &
         &isc, iec, jsc, jec)
    self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec
  
  end subroutine soca_conv_setup

  ! ------------------------------------------------------------------------------
  !> Apply forward convolution  
  subroutine soca_conv (self, convdx, dx)
    type(soca_horizconv), intent(in) :: self
    type(soca_field),    intent(in) :: dx
    type(soca_field),   intent(inout) :: convdx


  end subroutine soca_conv

  ! ------------------------------------------------------------------------------
  !> Apply backward convolution
  subroutine soca_conv_ad (self, convdx, dx)
    type(soca_horizconv), intent(in) :: self
    type(soca_field), intent(inout) :: dx     ! OUT
    type(soca_field),    intent(in) :: convdx ! IN

  end subroutine soca_conv_ad

end module soca_horizconv_mod
