! (C) Copyright 2020-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_getvalues_mod

use soca_geom_mod, only: soca_geom
use soca_state_mod, only: soca_state
use soca_fields_mod, only: soca_fields
use datetime_mod, only: datetime
use kinds, only: kind_real
use ufo_geovals_mod, only: ufo_geovals
use ufo_locs_mod, only: ufo_locs, ufo_locs_time_mask
use unstructured_interpolation_mod, only: unstrc_interp
use fckit_mpi_module, only: fckit_mpi_comm
use fckit_log_module, only : fckit_log

implicit none
private

type, public :: soca_getvalues
  type(unstrc_interp) :: horiz_interp
  type(unstrc_interp) :: horiz_interp_masked
contains

  ! constructors / destructors
  procedure :: create => soca_getvalues_create
  procedure :: delete => soca_getvalues_delete

  ! apply interpolation
  procedure :: fill_geovals=> soca_getvalues_fillgeovals

end type

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine soca_getvalues_create(self, geom, locs)
  class(soca_getvalues), intent(inout) :: self
  type(soca_geom),          intent(in) :: geom
  type(ufo_locs),           intent(in) :: locs

  integer :: isc, iec, jsc, jec, ni, nj
  type(fckit_mpi_comm) :: f_comm
  integer :: nn, ngrid_in, ngrid_out
  character(len=8) :: wtype = 'barycent'
  real(kind=kind_real), allocatable :: lats_in(:), lons_in(:)

  ! Indices for compute domain (no halo)
  isc = geom%isc ; iec = geom%iec
  jsc = geom%jsc ; jec = geom%jec

  f_comm = fckit_mpi_comm()
  ngrid_out = locs%nlocs
  nn = 3

  ! create interpolation weights for fields that do NOT use the land mask
  ni = iec - isc + 1
  nj = jec - jsc + 1
  ngrid_in = ni * nj
  allocate(lats_in(ngrid_in), lons_in(ngrid_in))

  lons_in = reshape(geom%lon(isc:iec,jsc:jec), (/ngrid_in/))
  lats_in = reshape(geom%lat(isc:iec,jsc:jec), (/ngrid_in/))
  call self%horiz_interp%create(f_comm, nn, wtype, &
                           ngrid_in, lats_in, lons_in, &
                           ngrid_out, locs%lat, locs%lon)

  ! create interpolation weights for fields that DO use the land mask
  deallocate(lats_in, lons_in)
  ngrid_in = count(geom%mask2d(isc:iec,jsc:jec) > 0)
  lons_in = pack(geom%lon(isc:iec,jsc:jec), mask=geom%mask2d(isc:iec,jsc:jec) > 0)
  lats_in = pack(geom%lat(isc:iec,jsc:jec), mask=geom%mask2d(isc:iec,jsc:jec) > 0)
  call self%horiz_interp_masked%create(f_comm, nn, wtype, &
                           ngrid_in, lats_in, lons_in, &
                           ngrid_out, locs%lat, locs%lon)

end subroutine soca_getvalues_create

!------------------------------------------------------------------------------
subroutine soca_getvalues_delete(self)
  class(soca_getvalues), intent(inout) :: self

end subroutine soca_getvalues_delete

!------------------------------------------------------------------------------
subroutine soca_getvalues_fillgeovals(self, geom, state, t1, t2, locs, geovals)
  class(soca_getvalues), intent(inout) :: self
  type(soca_geom),          intent(in) :: geom
  type(soca_state),         intent(in) :: state
  type(datetime),           intent(in) :: t1
  type(datetime),           intent(in) :: t2
  type(ufo_locs),           intent(in) :: locs
  type(ufo_geovals),     intent(inout) :: geovals

  logical, allocatable :: time_mask(:)
  real(kind=kind_real), allocatable :: geovals_all(:,:), geovals_tmp(:)
  integer :: gvar, nval

  ! Get mask for locations in this time window
  ! ------------------------------------------
  call ufo_locs_time_mask(locs, t1, t2, time_mask)

  ! Allocate geovals
  ! ----------------
  if (.not. geovals%linit) then
    do gvar = 1, geovals%nvar
      ! Set number of levels/categories (nval)
      nval = nlev_from_ufovar(state, geovals%variables(gvar))
      if (nval==0) cycle ! should we abort?

      ! Allocate GeoVaLs
      geovals%geovals(gvar)%nval = nval
      allocate(geovals%geovals(gvar)%vals(geovals%geovals(gvar)%nval, geovals%geovals(gvar)%nlocs))
      geovals%geovals(gvar)%vals = 0.0_kind_real
      geovals%linit = .true.
    end do
   endif
   geovals%linit = .true.



end subroutine soca_getvalues_fillgeovals

! ------------------------------------------------------------------------------
!> Get 3rd dimension of fld
function nlev_from_ufovar(fld, var) result(nval)
  type(soca_state), intent(in) :: fld
  character(len=*),  intent(in) :: var
  integer                       :: nval

  integer :: i

  ! fields that are not derived (i.e. the "cf_name" is set for a given field)
  do i=1,size(fld%fields)
    if (fld%fields(i)%cf_name == var) then
      nval = fld%fields(i)%nz
      return
    end if
  end do

  ! otherwise, the ufovar name is a derived type
  select case (var)
  case ("sea_surface_temperature", &
        "surface_temperature_where_sea", &
        "sea_surface_salinity", &
        "sea_floor_depth_below_sea_surface", &
        "sea_area_fraction")
     nval = 1

  case ("sea_water_salinity")
     nval = fld%geom%nzo

  case default
     nval = 0
     call fckit_log%debug("soca_interpfields_mod:nlef_from_ufovar geoval does not exist: "//var)

  end select

end function nlev_from_ufovar


end module soca_getvalues_mod
