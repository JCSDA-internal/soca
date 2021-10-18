! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> variable transform: vertical convolution
module soca_vertconv_mod

use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use tools_func, only: fit_func
use type_mpl, only: mpl_type
use type_probe, only: probe

! soca modules
use soca_fields_mod, only: soca_field
use soca_geom_mod, only: soca_geom
use soca_increment_mod, only:soca_increment
use soca_state_mod, only: soca_state

implicit none
private

!> Variable transform for vertical convolution
!!
!! \note this only operates on tocn and socn
type, public :: soca_vertconv
  real(kind=kind_real)      :: lz_min             !> Vertical decorrelation minimum [m]
  real(kind=kind_real)      :: lz_mld             !> if /= 0, Use MLD to calculate Lz
  real(kind=kind_real)      :: lz_mld_max         !> if calculating Lz from MLD, max value to use
  real(kind=kind_real)      :: scale_layer_thick  !> Set the minimum decorrelation scale
                                                  !> as a multiple of the layer thickness
  type(soca_state), pointer :: bkg                !> Background
  type(soca_geom),  pointer :: geom               !> Geometry
  integer                   :: isc, iec, jsc, jec !> Compute domain

contains
  !> \copybrief soca_conv_setup \see soca_conv_setup
  procedure :: setup => soca_conv_setup

  !> \copybrief soca_conv_mult \see soca_conv_mult
  procedure :: mult => soca_conv_mult

  !> \copybrief soca_conv_mult_ad \see soca_conv_mult_ad
  procedure :: mult_ad => soca_conv_mult_ad

end type soca_vertconv


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup for the vertical convolution
!!
!! \todo: Investigate computing and storing weights in vertconc data structure
!! \relates soca_vertconv_mod::soca_vertconv
subroutine soca_conv_setup (self, bkg, geom, f_conf)
  class(soca_vertconv),    intent(inout) :: self
  type(soca_state),  target, intent(in) :: bkg !< T/S background
  type(fckit_configuration), intent(in) :: f_conf !< yaml configuration
  type(soca_geom),   target, intent(in) :: geom !< input geometry

  ! Get configuration for vertical convolution
  call f_conf%get_or_die("Lz_min", self%lz_min )
  call f_conf%get_or_die("Lz_mld", self%lz_mld )
  if ( self%lz_mld /= 0) &
      call f_conf%get_or_die("Lz_mld_max", self%lz_mld_max )
  call f_conf%get_or_die("scale_layer_thick", self%scale_layer_thick )

  ! Associate background and geometry
  self%bkg => bkg
  self%geom => geom

  ! Indices for compute domain (no halo)
  self%isc=geom%isc; self%iec=geom%iec
  self%jsc=geom%jsc; self%jec=geom%jec

end subroutine soca_conv_setup


! ------------------------------------------------------------------------------
!> Calculate vertical correlation lengths for a given column
!!
!! \relates soca_vertconv_mod::soca_vertconv
subroutine soca_calc_lz(self, i, j, lz)
  class(soca_vertconv), intent(in) :: self
  integer, intent(in) :: i !< i index of grid point
  integer, intent(in) :: j !< j index of grid point
  real(kind=kind_real), intent(inout) :: lz(:) !< output correlation legnths

  real(kind=kind_real) :: mld, z
  integer :: k
  type(soca_field), pointer :: hocn, mld_fld, layer_depth

  ! minium scale is based on layer thickness
  call self%bkg%get("hocn", hocn)
  call self%bkg%get("mld", mld_fld)
  call self%bkg%get("layer_depth", layer_depth)
  lz = self%lz_min
  lz = max(lz, self%scale_layer_thick*abs(hocn%val(i,j,:)))

  ! if the upper Lz should be calculated from the MLD
  ! interpolate values from top to bottom of ML
  if ( self%lz_mld /= 0 ) then
     mld = mld_fld%val(i,j, 1)
     mld = min( mld, self%lz_mld_max)
     mld = max( mld, self%lz_min)
     do k=1, size(lz)
        z = layer_depth%val(i,j, k)
        if (z >= mld) exit  ! end of ML, exit loop
        lz(k) = max(lz(k), lz(k) + (mld - lz(k)) * (1.0 - z/mld))
     end do
  end if

end subroutine soca_calc_lz


! ------------------------------------------------------------------------------
!> Apply forward convolution
!!
!! \relates soca_vertconv_mod::soca_vertconv
subroutine soca_conv_mult (self, convdx, dx)
  class(soca_vertconv), intent(inout) :: self
  type(soca_increment),   intent(in) :: dx !< input increment to convolve
  type(soca_increment),intent(inout) :: convdx !< output increment

  real(kind=kind_real), allocatable :: z(:), lz(:)
  real(kind=kind_real) :: dist2, coef
  integer :: nl, j, k, id, jd, n
  type(mpl_type) :: mpl

  type(soca_field), pointer :: field_dx, field_convdx, layer_depth

  call probe%get_instance('soca')

  call self%bkg%get("layer_depth", layer_depth)
  nl = layer_depth%nz

  allocate(z(nl), lz(nl))

  do n=1,size(dx%fields)
    ! TODO remove these hardcoded values, use the yaml file
    select case(dx%fields(n)%name)
    case ("tocn", "socn")
      call dx%get(dx%fields(n)%name, field_dx)
      call convdx%get(dx%fields(n)%name, field_convdx)
      do id = self%isc, self%iec
        do jd = self%jsc, self%jec
          ! skip land
          if (self%geom%mask2d(id,jd) /= 1) cycle

          ! get correlation lengths
          call soca_calc_lz(self, id, jd, lz)

          ! perform convolution
          z(:) = layer_depth%val(id,jd,:)
          do j = 1, nl
            field_convdx%val(id,jd,j) = 0.0d0
            do k = 1,nl
              dist2 = abs(z(j)-z(k))
              coef = fit_func(mpl, dist2/lz(k))
              field_convdx%val(id,jd,j) = field_convdx%val(id,jd,j) &
                &+ field_dx%val(id,jd,k)*coef
            end do
          end do
        end do
      end do
    end select
  end do
  deallocate(z, lz)
end subroutine soca_conv_mult


! ------------------------------------------------------------------------------
!> Apply backward convolution
!!
!! \relates soca_vertconv_mod::soca_vertconv
subroutine soca_conv_mult_ad (self, convdx, dx)
  class(soca_vertconv), intent(inout) :: self
  type(soca_increment),intent(inout) :: dx     !< output increment
  type(soca_increment),   intent(in) :: convdx !< input increment

  real(kind=kind_real), allocatable :: z(:), lz(:)
  real(kind=kind_real) :: dist2, coef
  integer :: nl, j, k, id, jd, n
  type(mpl_type) :: mpl
  type(soca_field), pointer :: field_dx, field_convdx, layer_depth

  call probe%get_instance('soca')

  call self%bkg%get("layer_depth", layer_depth)
  nl = layer_depth%nz
  allocate(z(nl), lz(nl))

  do n=1,size(dx%fields)
    select case(dx%fields(n)%name)
   ! TODO remove these hardcoded values, use the yaml file
    case ("tocn", "socn")
      call dx%get(dx%fields(n)%name, field_dx)
      call convdx%get(dx%fields(n)%name, field_convdx)
      do id = self%isc, self%iec
        do jd = self%jsc, self%jec
          ! skip land
          if (self%geom%mask2d(id,jd) /= 1) cycle

          ! get correlation lengths
          call soca_calc_lz(self, id, jd, lz)

          ! perform convolution
          z(:) = layer_depth%val(id,jd,:)
          field_dx%val(id,jd,:) = 0.0d0
          do j = nl, 1, -1
            do k = nl, 1, -1
              dist2 = abs(z(j)-z(k))
              coef = fit_func(mpl, dist2/lz(k))
              field_dx%val(id,jd,k) = field_dx%val(id,jd,k) + coef*field_convdx%val(id,jd,j)
            end do
          end do
        end do
      end do
    end select
  end do
  deallocate(z, lz)

end subroutine soca_conv_mult_ad

end module soca_vertconv_mod
