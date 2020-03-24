! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerr_mod

use fckit_configuration_module, only: fckit_configuration
use datetime_mod, only: datetime
use kinds, only: kind_real
use soca_fields_mod
use soca_state_mod
use soca_increment_mod
use soca_bkgerrutil_mod, only: soca_bkgerr_bounds_type

implicit none

private
public :: soca_bkgerr_config, &
          soca_bkgerr_setup, soca_bkgerr_mult

!> Fortran derived type to hold configuration D
type :: soca_bkgerr_config
   type(soca_state),         pointer :: bkg
   type(soca_fields)                 :: std_bkgerr
   type(soca_bkgerr_bounds_type)     :: bounds         ! Bounds for bkgerr
   real(kind=kind_real)              :: std_sst
   real(kind=kind_real)              :: std_sss
   integer                           :: isc, iec, jsc, jec
end type soca_bkgerr_config

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup the static background error
subroutine soca_bkgerr_setup(f_conf, self, bkg)
  type(fckit_configuration),   intent(in) :: f_conf
  type(soca_bkgerr_config), intent(inout) :: self
  type(soca_state),    target, intent(in) :: bkg

  type(soca_field), pointer :: field, field_bkg

  integer :: isc, iec, jsc, jec, i
  type(datetime) :: vdate
  character(len=800) :: fname = 'soca_bkgerrsoca.nc'

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
    call f_conf%get_or_die("fixed_std_sst", self%std_sst)
    call self%std_bkgerr%get("tocn", field)
    field%val(:,:,1) = self%std_sst
  end if
  if ( f_conf%has("fixed_std_sss") ) then
      call f_conf%get_or_die("fixed_std_sss", self%std_sss)
      call self%std_bkgerr%get("socn", field)
      field%val(:,:,1) = self%std_sss
  end if

  ! Invent background error for ocnsfc fields: set it
  ! to 10% of the background for now ...
  ! TODO: Read background error for ocnsfc from file
  do i=1,size(self%std_bkgerr%fields)
    field => self%std_bkgerr%fields(i)
    select case(field%name)
    case ('sw','lw','lhf','shf','us')
      call bkg%get(field%name, field_bkg)
      field%val = abs(field_bkg%val) * 0.1_kind_real
    end select
  end do

  ! Associate background
  self%bkg => bkg

  ! Indices for compute domain (no halo)
  isc=bkg%geom%isc; iec=bkg%geom%iec
  jsc=bkg%geom%jsc; jec=bkg%geom%jec
  self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec

  ! Apply config bounds to background error
  call self%bounds%apply(self%std_bkgerr)

  ! Save filtered background error
  call self%std_bkgerr%write_file(fname)

end subroutine soca_bkgerr_setup

! ------------------------------------------------------------------------------
!> Apply background error: dxm = D dxa
subroutine soca_bkgerr_mult(self, dxa, dxm)
  type(soca_bkgerr_config),    intent(in) :: self
  type(soca_increment),        intent(in) :: dxa
  type(soca_increment),     intent(inout) :: dxm

  type(soca_field), pointer :: field_m, field_a, field_e

  integer :: isc, iec, jsc, jec, i, j, n

  ! make sure fields are correct shape
  call dxa%check_congruent(dxm)
  call dxa%check_subset(self%std_bkgerr)

  ! Indices for compute domain (no halo)
  isc=self%bkg%geom%isc; iec=self%bkg%geom%iec
  jsc=self%bkg%geom%jsc; jec=self%bkg%geom%jec

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
