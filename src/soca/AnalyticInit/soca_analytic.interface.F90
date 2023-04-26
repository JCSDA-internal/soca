! (C) Copyright 2021-2021 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


!> C++ interfaces for soca_analytic_mod
module soca_analytic_mod_c

use iso_c_binding

USE ufo_geovals_mod,            ONLY : ufo_geovals
USE ufo_geovals_mod_c,          ONLY : ufo_geovals_registry
USE ufo_sampled_locations_mod,  ONLY : ufo_sampled_locations

use soca_analytic_mod

implicit none
private

contains

!> C++ interface for soca_analytic_mod::soca_analytic_geovals()
subroutine soca_analytic_geovals_c(c_key_geovals, c_locs) &
        bind (c,name='soca_analytic_geovals_f90')

    integer(c_int),  intent(inout) :: c_key_geovals
    type(c_ptr), value, intent(in) :: c_locs

    type(ufo_geovals), pointer :: geovals
    type(ufo_sampled_locations) :: locs

    call ufo_geovals_registry%get(c_key_geovals, geovals)
    locs = ufo_sampled_locations(c_locs)

    call soca_analytic_geovals(geovals, locs)

end subroutine

end module
