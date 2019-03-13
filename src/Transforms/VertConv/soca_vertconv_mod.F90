!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_vertconv_mod
  use config_mod
  use iso_c_binding  
  use kinds
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices  
  use tools_func
  use type_mpl
  
  implicit none

  !> Fortran derived type to hold the setup for Vertconv
  type :: soca_vertconv
     real(kind=kind_real)      :: lz                 !> Vertical decorrelation [m]
     type(soca_field),pointer  :: traj               !> Trajectory
     type(soca_field), pointer :: bkg                !> Background     
     integer                   :: isc, iec, jsc, jec !> Compute domain 
  end type soca_vertconv

#define LISTED_TYPE soca_vertconv

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_vertconv_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------
  ! Setup for the vertical convolution
  ! TODO: Investigate computing and storing weights in vertconc data structure
  subroutine soca_conv_setup (self, bkg, traj, c_conf)
    type(soca_vertconv),   intent(inout) :: self
    type(soca_field), target, intent(in) :: bkg
    type(soca_field), target, intent(in) :: traj
    type(c_ptr),              intent(in) :: c_conf

    integer :: isc, iec, jsc, jec, i, j, k, nl
    
    nl = size(bkg%hocn,3)
  
    ! Get configuration for vertical convolution
    self%lz      = config_get_real(c_conf, "Lz")

    ! Store trajectory and background
    self%traj => traj
    self%bkg => bkg
    
    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(bkg%geom%ocean, "compute", &
         &isc, iec, jsc, jec)
    self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec
  
  end subroutine soca_conv_setup

  ! ------------------------------------------------------------------------------
  !> Apply forward convolution  
  subroutine soca_conv (self, convdx, dx)
    type(soca_vertconv), intent(in) :: self
    type(soca_field),    intent(in) :: dx
    type(soca_field),   intent(inout) :: convdx

    real(kind=kind_real), allocatable :: z(:), zp(:), lzd(:)
    real(kind=kind_real) :: dist2, coef, lz
    integer :: nl, j, k, id, jd
    type(mpl_type) :: mpl
    
    lz = self%lz
    nl = size(self%bkg%layer_depth,3)

    allocate(z(nl), zp(nl), lzd(nl))    
    do id = self%isc, self%iec
       do jd = self%jsc, self%jec
          if (self%bkg%geom%ocean%mask2d(id,jd).eq.1) then
             z(:) = self%bkg%layer_depth(id,jd,:)
             zp = z
             lzd = max(lz,4.0_kind_real*abs(self%bkg%hocn(id,jd,:)))
             do j = 1, nl
                convdx%tocn(id,jd,j) = 0.0d0
                convdx%socn(id,jd,j) = 0.0d0             
                do k = 1,nl
                   dist2 = abs(z(j)-zp(k))
                   coef = gc99(mpl, dist2/lzd(k))
                   convdx%tocn(id,jd,j) = convdx%tocn(id,jd,j) &
                        &+ dx%tocn(id,jd,k)*coef
                   convdx%socn(id,jd,j) = convdx%socn(id,jd,j) &
                        &+ dx%socn(id,jd,k)*coef
                end do
             end do
          end if
       end do
    end do
    deallocate(z, zp, lzd)

  end subroutine soca_conv

  ! ------------------------------------------------------------------------------
  !> Apply backward convolution
  subroutine soca_conv_ad (self, convdx, dx)
    type(soca_vertconv), intent(in) :: self
    type(soca_field), intent(inout) :: dx     ! OUT
    type(soca_field),    intent(in) :: convdx ! IN

    real(kind=kind_real), allocatable :: z(:), zp(:), lzd(:)
    real(kind=kind_real) :: dist2, coef, lz
    integer :: nl, j, k, id, jd
    type(mpl_type) :: mpl

    lz = self%lz
    nl = size(self%bkg%layer_depth,3)
    allocate(z(nl), zp(nl), lzd(nl))
    do id = self%isc, self%iec
       do jd = self%jsc, self%jec
          z(:) = self%bkg%layer_depth(id,jd,:)
          zp = z
          if (self%bkg%geom%ocean%mask2d(id,jd).eq.1) then
             dx%tocn(id,jd,:) = 0.0d0
             dx%socn(id,jd,:) = 0.0d0
             lzd = max(lz,4.0_kind_real*abs(self%bkg%hocn(id,jd,:)))
             do j = nl, 1, -1
                do k = nl, 1, -1
                   dist2 = abs(z(j)-zp(k))
                   coef = gc99(mpl, dist2/lzd(k))
                   dx%tocn(id,jd,k) = dx%tocn(id,jd,k) + coef*convdx%tocn(id,jd,j)
                   dx%socn(id,jd,k) = dx%socn(id,jd,k) + coef*convdx%socn(id,jd,j)
                end do
             end do
          end if
       end do
    end do
    deallocate(z, zp, lzd)

  end subroutine soca_conv_ad

end module soca_vertconv_mod
