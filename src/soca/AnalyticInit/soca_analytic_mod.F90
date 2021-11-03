! (C) Copyright 2021-2021 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_analytic_mod

use kinds, only : kind_real
use ufo_geovals_mod, only : ufo_geovals
use ufo_locations_mod, only : ufo_locations
use soca_state_mod, only : soca_state


implicit none
private

public :: soca_analytic_geovals
public :: soca_analytic_state

contains

subroutine soca_analytic_geovals(geovals, locs)
    type(ufo_geovals), intent(inout) :: geovals
    type(ufo_locations),  intent(in) :: locs

    real(kind=kind_real), allocatable :: lons(:), lats(:)
    integer :: ivar, iloc, ival
    real(kind=kind_real) :: val
    character(len=:), allocatable :: name

    allocate(lons(locs%nlocs()))
    allocate(lats(locs%nlocs()))
    call locs%get_lons(lons)
    call locs%get_lats(lats)

    do ivar = 1, geovals%nvar
        name = geovals%variables(ivar)
        do iloc = 1, geovals%geovals(ivar)%nlocs
            do ival = 1, geovals%geovals(ivar)%nval
                val = soca_analytic_val(&
                    name, lats(iloc), lons(iloc), ival*1.0_kind_real)
                geovals%geovals(ivar)%vals(ival, iloc) = val
            end do
        end do
    end do

end subroutine

subroutine soca_analytic_state(state)
    class(soca_state), intent(inout) :: state

    integer :: i, j, f, k
    real(kind=kind_real) :: val

    do f = 1, size(state%fields)
        do i = state%geom%isc, state%geom%iec
            do j = state%geom%jsc, state%geom%jec
                do k = 1, state%fields(f)%nz
                ! TODO use the depth instead of the depth level
                    val = soca_analytic_val(&
                        state%fields(f)%name, state%fields(f)%lat(i,j), &
                        state%fields(f)%lon(i,j), k*1.0_kind_real)
                    state%fields(f)%val(i,j,k) = val
                end do
            end do
        end do
    end do
end subroutine

function soca_analytic_val(var, lat, lon, depth) result(val)
    character(len=:), allocatable, intent(inout) :: var
    real(kind=kind_real), intent(in) :: lat, lon, depth
    real(kind=kind_real) :: val

    integer :: i, j
    real(kind=kind_real) :: rvar
    integer :: ivar

    ! create hash from string
    ivar = 0
    do i=1,len(var)
        ivar = ieor(ivar, ishft(ichar(var(i:i)), mod(i-1,4)*8))
    end do
    rvar = mod(abs(ieor(ivar, z'5D7A9F43')), 1000) / 1000.0

    val = (sind(lon) + lat/90.0) * (0.5/depth) + rvar

end function


end module