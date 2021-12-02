! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> variable transform: background error
module soca_bkgerr_mod

use datetime_mod, only: datetime
use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real

! soca modules
use soca_bkgerrutil_mod, only: soca_bkgerr_bounds_type
use soca_fields_mod, only: soca_fields, soca_field
use soca_geom_mod, only: soca_geom
use soca_increment_mod, only: soca_increment
use soca_state_mod, only: soca_state

implicit none
private


!> Variable transform for background error
type, public :: soca_bkgerr
  type(soca_fields) :: std_bkgerr

  ! private members
  type(soca_bkgerr_bounds_type) , private  :: bounds !< Bounds for bkgerr
  type(soca_geom),  pointer, private       :: geom !< geometry

contains

  !> \copybrief soca_bkgerr_setup \see soca_bkgerr_setup
  procedure :: setup => soca_bkgerr_setup

  !> \copybrief soca_bkgerr_mult \see soca_bkgerr_mult
  procedure :: mult => soca_bkgerr_mult
end type soca_bkgerr


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Setup the static background error
!!
!! \note the precomputed standard devations in std_bkgerr are only used
!! for tocn, socn, and ssh.
!! \relates soca_bkgerr_mod::soca_bkgerr
subroutine soca_bkgerr_setup(self, f_conf, bkg, geom)
  class(soca_bkgerr), target,  intent(inout) :: self
  type(fckit_configuration),   intent(in)    :: f_conf !< configuration
  type(soca_state),    target, intent(in)    :: bkg !< background
  type(soca_geom),     target, intent(in)    :: geom !< geometry

  type(soca_field), pointer :: field, field_bkg
  real(kind=kind_real) :: std
  integer :: i
  type(datetime) :: vdate
  character(len=800) :: fname = 'soca_bkgerrsoca.nc'

  self%geom => geom

  ! Allocate memory for bkgerror
  call self%std_bkgerr%copy(bkg)
  !call create_copy(self%std_bkgerr, bkg)

  ! Read variance
  ! Precomputed from an ensemble of (K^-1 dx)
  call self%std_bkgerr%read(f_conf, vdate)

  ! Convert to standard deviation
  do i=1,size(self%std_bkgerr%fields)
    field => self%std_bkgerr%fields(i)
    select case(field%name)
    case ("tocn", "socn", "ssh")
      field%val = sqrt(field%val)
    end select
  end do

  ! Get bounds from configuration
  call self%bounds%read(f_conf)

  ! Get constand background error for sst and sss
  if ( f_conf%has("fixed_std_sst") ) then
    call f_conf%get_or_die("fixed_std_sst", std)
    call self%std_bkgerr%get("tocn", field)
    field%val(:,:,1) = std
  end if
  if ( f_conf%has("fixed_std_sss") ) then
      call f_conf%get_or_die("fixed_std_sss", std)
      call self%std_bkgerr%get("socn", field)
      field%val(:,:,1) = std
  end if

  ! Invent background error for ocnsfc and ocn_bgc fields:
  ! set it to 10% or 20% of the background for now ...
  ! TODO: Read background error for ocnsfc and ocn_bgc from
  ! files
  do i=1,size(self%std_bkgerr%fields)
    field => self%std_bkgerr%fields(i)
    select case(field%name)
    case ('sw','lw','lhf','shf','us')
      call bkg%get(field%name, field_bkg)
      field%val = abs(field_bkg%val) * 0.1_kind_real
    case ('chl','biop')
      call bkg%get(field%name, field_bkg)
      field%val = abs(field_bkg%val) * 0.2_kind_real
    end select
  end do

  ! Apply config bounds to background error
  call self%bounds%apply(self%std_bkgerr)

  ! Save filtered background error
  call self%std_bkgerr%write_file(fname)

end subroutine soca_bkgerr_setup


! ------------------------------------------------------------------------------
!> Apply background error: dxm = D dxa
!!
!! \relates soca_bkgerr_mod::soca_bkgerr
subroutine soca_bkgerr_mult(self, dxa, dxm)
  class(soca_bkgerr),           intent(in)    :: self
  type(soca_increment), target, intent(in)    :: dxa !< input increment
  type(soca_increment),         intent(inout) :: dxm !< output increment

  type(soca_field), pointer :: field_m, field_a, field_e

  integer :: isc, iec, jsc, jec, i, j, n

  ! make sure fields are correct shape
  call dxa%check_congruent(dxm)
  call dxa%check_subset(self%std_bkgerr)

  ! Indices for compute domain (no halo)
  isc=self%geom%isc; iec=self%geom%iec
  jsc=self%geom%jsc; jec=self%geom%jec

  ! multiply
  do n=1,size(dxa%fields)
    field_a => dxa%fields(n)
    call self%std_bkgerr%get(field_a%name, field_e)
    call dxm%get(field_a%name, field_m)
    do i = isc, iec
      do j = jsc, jec
        field_m%val(i,j,:) = field_e%val(i,j,:) * field_a%val(i,j,:)
      end do
    end do
  end do
end subroutine soca_bkgerr_mult

! ------------------------------------------------------------------------------

end module soca_bkgerr_mod
