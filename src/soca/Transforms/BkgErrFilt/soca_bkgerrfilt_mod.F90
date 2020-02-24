! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerrfilt_mod

use fckit_configuration_module, only: fckit_configuration
use datetime_mod, only: datetime
use kinds, only: kind_real
use soca_fields_mod, only: soca_fields, soca_field, soca_fld2file

implicit none

private
public :: soca_bkgerrfilt_config, &
          soca_bkgerrfilt_setup, soca_bkgerrfilt_mult

!> Fortran derived type to hold configuration
type :: soca_bkgerrfilt_config
   type(soca_fields),   pointer :: bkg
   type(soca_fields)            :: filt
   real(kind=kind_real)         :: efold_z           ! E-folding scale
   real(kind=kind_real)         :: scale             ! Rescaling factor
   real(kind=kind_real)         :: ocn_depth_min     ! Minimum depth
   integer                      :: isc, iec, jsc, jec
end type soca_bkgerrfilt_config

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup the static background error
subroutine soca_bkgerrfilt_setup(f_conf, self, bkg)
  type(fckit_configuration),    intent(in)    :: f_conf
  type(soca_bkgerrfilt_config), intent(inout) :: self
  type(soca_fields), target,    intent(in)    :: bkg

  integer :: isc, iec, jsc, jec, i, j, k, nl
  real(kind=kind_real) :: efold
  character(len=800) :: fname = 'soca_bkgerrfilt.nc'
  type(soca_field), pointer :: tocn

  ! Get number of ocean levels
  nl = size(bkg%hocn,3)

  ! Allocate memory for bkgerrfiltor and set to zero
  !call create_copy(self%filt, bkg)
  call self%filt%copy(bkg)
  call self%filt%zeros()

  ! Read parameters from config
  call f_conf%get_or_die("ocean_depth_min", self%ocn_depth_min)
  call f_conf%get_or_die("rescale_bkgerr", self%scale)
  call f_conf%get_or_die("efold_z", self%efold_z)

  ! Associate background
  self%bkg => bkg

  call self%filt%get("tocn", tocn)

  ! Setup rescaling and masks
  isc=bkg%geom%isc ; self%isc=isc ; iec=bkg%geom%iec ; self%iec=iec
  jsc=bkg%geom%jsc ; self%jsc=jsc ; jec=bkg%geom%jec ; self%jec=jec
  do i = isc, iec
     do j = jsc, jec
        if (sum(bkg%hocn(i,j,:)).gt.self%ocn_depth_min) then
           self%filt%ssh(i,j) = self%scale
           do k = 1, nl
              if (bkg%hocn(i,j,k).gt.1e-3_kind_real) then
                 ! Only apply if layer is thick enough
                 efold = self%scale*exp(-self%bkg%layer_depth(i,j,k)/self%efold_z)
              else
                 ! Set to zero if layer is too thin
                 efold = 0.0_kind_real
              end if
              tocn%val(i,j,k) = efold
              self%filt%socn(i,j,k) = efold
           end do
        else
           ! Set to zero if ocean is too shallow
           self%filt%ssh(i,j)    = 0.0_kind_real
           tocn%val(i,j,:) = 0.0_kind_real
           self%filt%socn(i,j,:) = 0.0_kind_real
        end if

        ! Do nothing for sea-ice
        self%filt%seaice%cicen(i,j,:) =  1.0_kind_real
        self%filt%seaice%hicen(i,j,:) =  1.0_kind_real
     end do
  end do

  ! Save filtered background error
  call soca_fld2file(self%filt, fname)

end subroutine soca_bkgerrfilt_setup

! ------------------------------------------------------------------------------
!> Apply background error: dxm = D dxa
subroutine soca_bkgerrfilt_mult(self, dxa, dxm)
  type(soca_bkgerrfilt_config), intent(in) :: self
  type(soca_fields),            intent(in) :: dxa
  type(soca_fields),         intent(inout) :: dxm

  integer :: i, j
  type(soca_field), pointer :: field, field_a, field_m


  call self%filt%get("tocn", field)
  call dxa%get("tocn", field_a)
  call dxm%get("tocn", field_m)

  do i = self%isc, self%iec
     do j = self%jsc, self%jec
        if (self%bkg%geom%mask2d(i,j).eq.1) then
           dxm%ssh(i,j) = self%filt%ssh(i,j) * dxa%ssh(i,j)
           field_m%val(i,j,:) = field%val(i,j,:) * field_a%val(i,j,:)
           dxm%socn(i,j,:) = self%filt%socn(i,j,:)  * dxa%socn(i,j,:)

           dxm%seaice%cicen(i,j,:) =  self%filt%seaice%cicen(i,j,:) * dxa%seaice%cicen(i,j,:)
           dxm%seaice%hicen(i,j,:) =  self%filt%seaice%hicen(i,j,:) * dxa%seaice%hicen(i,j,:)
        else
           dxm%ssh(i,j) = 0.0_kind_real
           field_m%val(i,j,:) = 0.0_kind_real
           dxm%socn(i,j,:) = 0.0_kind_real

           dxm%seaice%cicen(i,j,:) = 0.0_kind_real
           dxm%seaice%hicen(i,j,:) = 0.0_kind_real
        end if
     end do
  end do
  ! Surface fields
  call dxm%ocnsfc%copy(dxa%ocnsfc)

end subroutine soca_bkgerrfilt_mult

! ------------------------------------------------------------------------------

end module soca_bkgerrfilt_mod


