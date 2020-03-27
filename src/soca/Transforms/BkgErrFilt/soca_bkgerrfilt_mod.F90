! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerrfilt_mod

use fckit_configuration_module, only: fckit_configuration
use datetime_mod, only: datetime
use kinds, only: kind_real
use soca_fields_mod
use soca_increment_mod
use soca_state_mod

implicit none

private
public :: soca_bkgerrfilt_config, &
          soca_bkgerrfilt_setup, soca_bkgerrfilt_mult

!> Fortran derived type to hold configuration
type :: soca_bkgerrfilt_config
   type(soca_state),    pointer :: bkg
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
  type(soca_state),  target,    intent(in)    :: bkg

  integer :: isc, iec, jsc, jec, i, j, k
  real(kind=kind_real) :: efold
  character(len=800) :: fname = 'soca_bkgerrfilt.nc'
  type(soca_field), pointer :: tocn, socn, ssh, hocn, layer_depth

  ! Allocate memory for bkgerrfiltor and set to zero
  call self%filt%copy(bkg)
  call self%filt%zeros()

  ! Read parameters from config
  call f_conf%get_or_die("ocean_depth_min", self%ocn_depth_min)
  call f_conf%get_or_die("rescale_bkgerr", self%scale)
  call f_conf%get_or_die("efold_z", self%efold_z)

  ! Associate background
  self%bkg => bkg

  call self%filt%get("tocn", tocn)
  call self%filt%get("socn", socn)
  call self%filt%get("ssh", ssh)
  call bkg%get("hocn", hocn)
  call bkg%get("layer_depth", layer_depth)

  ! Setup rescaling and masks
  isc=bkg%geom%isc ; self%isc=isc ; iec=bkg%geom%iec ; self%iec=iec
  jsc=bkg%geom%jsc ; self%jsc=jsc ; jec=bkg%geom%jec ; self%jec=jec
  do i = isc, iec
     do j = jsc, jec
        if (sum(hocn%val(i,j,:)).gt.self%ocn_depth_min) then
           ssh%val(i,j,:) = self%scale
           do k = 1, hocn%nz
              if (hocn%val(i,j,k).gt.1e-3_kind_real) then
                 ! Only apply if layer is thick enough
                 efold = self%scale*exp(-layer_depth%val(i,j,k)/self%efold_z)
              else
                 ! Set to zero if layer is too thin
                 efold = 0.0_kind_real
              end if
              tocn%val(i,j,k) = efold
              socn%val(i,j,k) = efold
           end do
        else
           ! Set to zero if ocean is too shallow
           ssh%val(i,j,:)  = 0.0_kind_real
           tocn%val(i,j,:) = 0.0_kind_real
           socn%val(i,j,:) = 0.0_kind_real
        end if
     end do
  end do

  ! set other things to 1
  do i=1,size(self%filt%fields)
    select case(self%filt%fields(i)%name)
    case ('ssh','tocn','socn')
      continue
    case default
      self%filt%fields(i)%val = 1.0_kind_real
    end select
  end do

  ! Save filtered background error
  call self%filt%write_file(fname)

end subroutine soca_bkgerrfilt_setup

! ------------------------------------------------------------------------------
!> Apply background error: dxm = D dxa
subroutine soca_bkgerrfilt_mult(self, dxa, dxm)
  type(soca_bkgerrfilt_config), intent(in) :: self
  type(soca_increment),         intent(in) :: dxa
  type(soca_increment),      intent(inout) :: dxm

  integer :: i, j, n
  type(soca_field), pointer :: field_f, field_a, field_m

  ! make sure fields are the right shape
  call dxa%check_congruent(dxm)
  call dxa%check_subset(self%filt)

  ! multiply
  do n=1,size(dxa%fields)
    field_a => dxa%fields(n)
    call self%filt%get(field_a%name, field_f)
    call dxm%get(field_a%name, field_m)
    do i = self%isc, self%iec
      do j = self%jsc, self%jec
        if (self%bkg%geom%mask2d(i,j).eq.1) then
          field_m%val(i,j,:) = field_f%val(i,j,:) * field_a%val(i,j,:)
        else
          field_m%val(i,j,:) = 0.0_kind_real
        end if
      end do
    end do
  end do
end subroutine soca_bkgerrfilt_mult

! ------------------------------------------------------------------------------

end module soca_bkgerrfilt_mod


