! (C) Copyright 2017- UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_vertconv_mod
  use config_mod
  use iso_c_binding
  use kinds
  use soca_fields
  use tools_func
  use type_mpl

  implicit none

  !> Fortran derived type to hold the setup for Vertconv
  type :: soca_vertconv
     real(kind=kind_real)      :: lz_min             !> Vertical decorrelation minimum [m]
     real(kind=kind_real)      :: lz_mld             !> if /= 0, Use MLD to calculate Lz
     real(kind=kind_real)      :: lz_mld_max         !> if calculating Lz from MLD, max value to use
     real(kind=kind_real)      :: scale_layer_thick  !> Set the minimum decorrelation scale
                                                     !> as a multiple of the layer thickness
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

    ! Get configuration for vertical convolution
    self%lz_min       = config_get_real(c_conf, "Lz_min")
    self%lz_mld       = config_get_int(c_conf, "Lz_mld")
    if ( self%lz_mld /= 0) then
      self%lz_mld_max   = config_get_real(c_conf, "Lz_mld_max")
    end if
    self%scale_layer_thick = config_get_real(c_conf, "scale_layer_thick")

    ! Store trajectory and background
    self%traj => traj
    self%bkg => bkg

    ! Indices for compute domain (no halo)
    self%isc=bkg%geom%isc; self%iec=bkg%geom%iec
    self%jsc=bkg%geom%jsc; self%jec=bkg%geom%jec

  end subroutine soca_conv_setup

  ! ------------------------------------------------------------------------------
  !> Calculate vertical correlation lengths for a given column
  subroutine soca_calc_lz(self, i, j, lz)
    type(soca_vertconv), intent(in) :: self
    integer, intent(in) :: i, j
    real(kind=kind_real), intent(inout) :: lz(:)
    real(kind=kind_real) :: mld, z
    integer :: k

    ! minium scale is based on layer thickness
    lz = self%lz_min
    lz = max(lz, self%scale_layer_thick*abs(self%bkg%hocn(i,j,:)))

    ! if the upper Lz should be calculated from the MLD
    ! interpolate values from top to bottom of ML
    if ( self%lz_mld /= 0 ) then
       mld = self%bkg%mld(i,j)
       mld = min( mld, self%lz_mld_max)
       mld = max( mld, self%lz_min)
       do k=1, size(lz)
          z = self%bkg%layer_depth(i,j, k)
          if (z >= mld) exit  ! end of ML, exit loop
          lz(k) = max(lz(k), lz(k) + (mld - lz(k)) * (1.0 - z/mld))
       end do
    end if

  end subroutine soca_calc_lz

  ! ------------------------------------------------------------------------------
  !> Apply forward convolution
  subroutine soca_conv (self, convdx, dx)
    type(soca_vertconv), intent(in) :: self
    type(soca_field),    intent(in) :: dx
    type(soca_field),   intent(inout) :: convdx

    real(kind=kind_real), allocatable :: z(:), lz(:)
    real(kind=kind_real) :: dist2, coef, mld
    integer :: nl, j, k, id, jd
    type(mpl_type) :: mpl

    nl = size(self%bkg%layer_depth,3)

    allocate(z(nl), lz(nl))
    do id = self%isc, self%iec
       do jd = self%jsc, self%jec
          ! skip land
          if (self%bkg%geom%mask2d(id,jd) /= 1) cycle

          ! get correlation lengths
          call soca_calc_lz(self, id, jd, lz)

          ! perform convolution
          z(:) = self%bkg%layer_depth(id,jd,:)
          do j = 1, nl
            convdx%tocn(id,jd,j) = 0.0d0
            convdx%socn(id,jd,j) = 0.0d0
            do k = 1,nl
              dist2 = abs(z(j)-z(k))
              coef = gc99(mpl, dist2/lz(k))
              convdx%tocn(id,jd,j) = convdx%tocn(id,jd,j) &
                &+ dx%tocn(id,jd,k)*coef
              convdx%socn(id,jd,j) = convdx%socn(id,jd,j) &
                &+ dx%socn(id,jd,k)*coef
              end do
          end do
       end do
    end do
    deallocate(z, lz)
  end subroutine soca_conv


  ! ------------------------------------------------------------------------------
  !> Apply backward convolution
  subroutine soca_conv_ad (self, convdx, dx)
    type(soca_vertconv), intent(in) :: self
    type(soca_field), intent(inout) :: dx     ! OUT
    type(soca_field),    intent(in) :: convdx ! IN

    real(kind=kind_real), allocatable :: z(:), lz(:)
    real(kind=kind_real) :: dist2, coef
    integer :: nl, j, k, id, jd
    type(mpl_type) :: mpl

    nl = size(self%bkg%layer_depth,3)
    allocate(z(nl), lz(nl))

    do id = self%isc, self%iec
       do jd = self%jsc, self%jec
          ! skip land
          if (self%bkg%geom%mask2d(id,jd) /= 1) cycle

          ! get correlation lengths
          call soca_calc_lz(self, id, jd, lz)

          ! perform convolution
          z(:) = self%bkg%layer_depth(id,jd,:)
          dx%tocn(id,jd,:) = 0.0d0
          dx%socn(id,jd,:) = 0.0d0
          do j = nl, 1, -1
             do k = nl, 1, -1
                dist2 = abs(z(j)-z(k))
                coef = gc99(mpl, dist2/lz(k))
                dx%tocn(id,jd,k) = dx%tocn(id,jd,k) + coef*convdx%tocn(id,jd,j)
                dx%socn(id,jd,k) = dx%socn(id,jd,k) + coef*convdx%socn(id,jd,j)
             end do
          end do
       end do
    end do
    deallocate(z, lz)

  end subroutine soca_conv_ad


end module soca_vertconv_mod
