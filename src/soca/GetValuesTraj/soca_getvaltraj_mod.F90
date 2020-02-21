! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module handling interpolation trajectory

module soca_getvaltraj_mod

use unstructured_interpolation_mod

implicit none

private
public :: soca_getvaltraj

type :: soca_getvaltraj
 integer                 :: nobs
 type(unstrc_interp)     :: horiz_interp
 type(unstrc_interp)     :: horiz_interp_masked
 logical                 :: interph_initialized = .false.
end type soca_getvaltraj

end module soca_getvaltraj_mod
