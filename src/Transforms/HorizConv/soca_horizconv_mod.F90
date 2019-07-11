!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_horizconv_mod
  use config_mod
  use datetime_mod
  use iso_c_binding
  use soca_fields

  implicit none

  private
  public :: soca_horizconv_setup, soca_horizconv_mult, soca_horizconv_type
    
  !> Fortran derived type to hold configuration D
  type :: soca_horizconv_type

  end type soca_horizconv_type

contains

  ! ------------------------------------------------------------------------------
  !> Setup the static background error
  subroutine soca_horizconv_setup(c_conf, self, bkg)
    type(soca_horizconv_type), intent(inout) :: self
    type(soca_field),    target, intent(in) :: bkg
    type(c_ptr),                 intent(in) :: c_conf

    integer :: isc, iec, jsc, jec, i, j, k, nl

  end subroutine soca_horizconv_setup

  ! ------------------------------------------------------------------------------
  !> Apply background error: dxm = D dxa
  subroutine soca_horizconv_mult(self, dxa, dxm)
    type(soca_horizconv_type),    intent(in) :: self
    type(soca_field),            intent(in) :: dxa
    type(soca_field),         intent(inout) :: dxm

    integer :: isc, iec, jsc, jec, i, j, k

  end subroutine soca_horizconv_mult

  ! ------------------------------------------------------------------------------

end module soca_horizconv_mod


