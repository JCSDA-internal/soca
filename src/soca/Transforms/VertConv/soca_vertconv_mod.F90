! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_vertconv_mod
   
use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use type_mpl, only: mpl_type
use tools_func, only: fit_func
use soca_fields_mod, only: soca_fields, soca_field

implicit none

private
public :: soca_vertconv, &
          soca_conv_setup, soca_conv, soca_conv_ad, soca_calc_lz

!> Fortran derived type to hold the setup for Vertconv
type :: soca_vertconv
   real(kind=kind_real)      :: lz_min             !> Vertical decorrelation minimum [m]
   real(kind=kind_real)      :: lz_mld             !> if /= 0, Use MLD to calculate Lz
   real(kind=kind_real)      :: lz_mld_max         !> if calculating Lz from MLD, max value to use
   real(kind=kind_real)      :: scale_layer_thick  !> Set the minimum decorrelation scale
                                                   !> as a multiple of the layer thickness
   type(soca_fields),pointer :: traj               !> Trajectory
   type(soca_fields),pointer :: bkg                !> Background
   integer                   :: isc, iec, jsc, jec !> Compute domain
end type soca_vertconv

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! Setup for the vertical convolution
! TODO: Investigate computing and storing weights in vertconc data structure
subroutine soca_conv_setup (self, bkg, traj, f_conf)
  type(fckit_configuration), intent(in) :: f_conf
  type(soca_vertconv),    intent(inout) :: self
  type(soca_fields), target, intent(in) :: bkg
  type(soca_fields), target, intent(in) :: traj

  ! Get configuration for vertical convolution
  call f_conf%get_or_die("Lz_min", self%lz_min )
  call f_conf%get_or_die("Lz_mld", self%lz_mld )
  if ( self%lz_mld /= 0) &
      call f_conf%get_or_die("Lz_mld_max", self%lz_mld_max )
  call f_conf%get_or_die("scale_layer_thick", self%scale_layer_thick )

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
  type(soca_fields),   intent(in) :: dx
  type(soca_fields),intent(inout) :: convdx

  real(kind=kind_real), allocatable :: z(:), lz(:)
  real(kind=kind_real) :: dist2, coef
  integer :: nl, j, k, id, jd, n
  type(mpl_type) :: mpl

  type(soca_field), pointer :: field_dx, field_convdx

  nl = size(self%bkg%layer_depth,3)

  allocate(z(nl), lz(nl))

  do n=1,size(dx%fields)
    select case(dx%fields(n)%name)
    case ("tocn","socn")
      call dx%get(dx%fields(n)%name, field_dx)
      call convdx%get(dx%fields(n)%name, field_convdx)
      do id = self%isc, self%iec
        do jd = self%jsc, self%jec
          ! skip land
          if (self%bkg%geom%mask2d(id,jd) /= 1) cycle

          ! get correlation lengths
          call soca_calc_lz(self, id, jd, lz)

          ! perform convolution
          z(:) = self%bkg%layer_depth(id,jd,:)
          do j = 1, nl
            field_convdx%val(id,jd,j) = 0.0d0
            do k = 1,nl
              dist2 = abs(z(j)-z(k))
              coef = fit_func(mpl, 'gc99', dist2/lz(k))
              field_convdx%val(id,jd,j) = field_convdx%val(id,jd,j) &
                &+ field_dx%val(id,jd,k)*coef
            end do
          end do
        end do
      end do
    end select
  end do
  deallocate(z, lz)
end subroutine soca_conv

! ------------------------------------------------------------------------------
!> Apply backward convolution
subroutine soca_conv_ad (self, convdx, dx)
  type(soca_vertconv), intent(in) :: self
  type(soca_fields),intent(inout) :: dx     ! OUT
  type(soca_fields),   intent(in) :: convdx ! IN

  real(kind=kind_real), allocatable :: z(:), lz(:)
  real(kind=kind_real) :: dist2, coef
  integer :: nl, j, k, id, jd, n
  type(mpl_type) :: mpl
  type(soca_field), pointer :: field_dx, field_convdx

  nl = size(self%bkg%layer_depth,3)
  allocate(z(nl), lz(nl))

  do n=1,size(dx%fields)
    select case(dx%fields(n)%name)
    case ("tocn","socn")
      call dx%get(dx%fields(n)%name, field_dx)
      call convdx%get(dx%fields(n)%name, field_convdx)
      do id = self%isc, self%iec
        do jd = self%jsc, self%jec
          ! skip land
          if (self%bkg%geom%mask2d(id,jd) /= 1) cycle

          ! get correlation lengths
          call soca_calc_lz(self, id, jd, lz)

          ! perform convolution
          z(:) = self%bkg%layer_depth(id,jd,:)
          field_dx%val(id,jd,:) = 0.0d0
          do j = nl, 1, -1
            do k = nl, 1, -1
              dist2 = abs(z(j)-z(k))
              coef = fit_func(mpl, 'gc99', dist2/lz(k))
              field_dx%val(id,jd,k) = field_dx%val(id,jd,k) + coef*field_convdx%val(id,jd,j)
            end do
          end do
        end do
      end do
    end select
  end do
  deallocate(z, lz)

end subroutine soca_conv_ad

end module soca_vertconv_mod
