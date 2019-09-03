! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module handling interpolation trajectory

module soca_getvaltraj_mod

use soca_bumpinterp2d_mod, only: soca_bumpinterp2d

implicit none

private
public :: soca_getvaltraj

type :: soca_getvaltraj
 integer                 :: nobs
 type(soca_bumpinterp2d), allocatable :: horiz_interp(:)
 integer                 :: bumpid
 logical                 :: interph_initialized = .false.
 integer                 :: obstype_index
end type soca_getvaltraj

end module soca_getvaltraj_mod
