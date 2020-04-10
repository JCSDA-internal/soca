! (C) Copyright 2020-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_getvalues_mod

use soca_geom_mod, only: soca_geom
use soca_state_mod, only: soca_state
use soca_fields_mod, only: soca_fields, soca_field
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

  call self%horiz_interp%delete()
  call self%horiz_interp_masked%delete()

end subroutine soca_getvalues_delete

!------------------------------------------------------------------------------
subroutine soca_getvalues_fillgeovals(self, geom, fld, t1, t2, locs, geovals)
  class(soca_getvalues), intent(inout) :: self
  type(soca_geom),          intent(in) :: geom
  type(soca_state),         intent(in) :: fld
  type(datetime),           intent(in) :: t1
  type(datetime),           intent(in) :: t2
  type(ufo_locs),           intent(in) :: locs
  type(ufo_geovals),     intent(inout) :: geovals

  logical, allocatable :: time_mask(:)
  logical :: masked
  integer :: ivar, nlocs, n
  integer :: ival, nval, indx
  integer :: isc, iec, jsc, jec
  integer :: ns
  real(kind=kind_real), allocatable :: gom_window(:)
  real(kind=kind_real), allocatable :: fld3d(:,:,:), fld3d_un(:)
  type(soca_field), pointer :: fldptr

  ! Indices for compute domain (no halo)
  isc = geom%isc ; iec = geom%iec
  jsc = geom%jsc ; jec = geom%jec

  ! Get mask for locations in this time window
  call ufo_locs_time_mask(locs, t1, t2, time_mask)

  ! Allocate temporary geoval and 3d field for the current time window
  do ivar = 1, geovals%nvar

    ! Set number of levels/categories (nval)
    nval = nlev_from_ufovar(fld, geovals%variables(ivar))
    print *,ivar, trim(geovals%variables(ivar)),nval,geovals%geovals(ivar)%nlocs
    read(*,*)
    if (nval==0) cycle ! should we abort?

    ! Allocate geovals
    geovals%geovals(ivar)%nval = nval
    if (.not. geovals%linit) then
      allocate(geovals%geovals(ivar)%vals(nval, geovals%geovals(ivar)%nlocs))
      geovals%geovals(ivar)%vals = 0.0_kind_real
    end if

    ! Return if no observations
    if ( geovals%geovals(ivar)%nlocs == 0 ) return

    allocate(gom_window(locs%nlocs))
    allocate(fld3d(isc:iec,jsc:jec,1:nval))
    nullify(fldptr)

    ! Extract fld3d from field
    masked = .true. ! by default fields are assumed to need a land mask applied,
    ! (Currently only sea area fraction is unmasked)

    ! if we are lucky and this variable is a non-derived type, check the fields structure
    do n=1,size(fld%fields)
      if (fld%fields(n)%cf_name == geovals%variables(ivar)) then
        fld3d = fld%fields(n)%val(isc:iec,jsc:jec,1:nval)
      end if
    end do

    ! otherwise, we are dealing with a derived type, prepare for a long "select case" statement
    select case (trim(geovals%variables(ivar)))
    case ("sea_water_salinity")
      call fld%get("socn", fldptr)
      fld3d = fldptr%val(isc:iec,jsc:jec,1:nval)

    case ("sea_surface_temperature")
      call fld%get("tocn", fldptr)
      fld3d(isc:iec,jsc:jec,1) = fldptr%val(isc:iec,jsc:jec,1)

      ! TODO: Move unit change elsewhere, check if is COARDS.
    case ("surface_temperature_where_sea")
      call fld%get("tocn", fldptr)
      fld3d(isc:iec,jsc:jec,1) = fldptr%val(isc:iec,jsc:jec,1) + 273.15_kind_real

    case ("sea_surface_salinity")
      call fld%get("socn", fldptr)
      fld3d(isc:iec,jsc:jec,1) = fldptr%val(isc:iec,jsc:jec,1)

    case ("sea_floor_depth_below_sea_surface")
      call fld%get("hocn", fldptr)
      fld3d(isc:iec,jsc:jec,1) = sum(fldptr%val(isc:iec,jsc:jec,:),dim=3)

    case ("sea_area_fraction")
      fld3d(isc:iec,jsc:jec,1) = real(fld%geom%mask2d(isc:iec,jsc:jec),kind=kind_real)
      masked = .false.

    case default
      call fckit_log%debug("soca_interpfields_mod:interp geoval does not exist")
    end select

    ! Apply forward interpolation: Model ---> Obs
    do ival = 1, nval
      if (masked) then
        ns = count(fld%geom%mask2d(isc:iec,jsc:jec) > 0 )
        if (.not. allocated(fld3d_un)) allocate(fld3d_un(ns))
        fld3d_un = pack(fld3d(isc:iec,jsc:jec,ival), mask=fld%geom%mask2d(isc:iec,jsc:jec) > 0)
        call self%horiz_interp_masked%apply(fld3d_un, gom_window)
      else
        ns = (iec - isc + 1) * (jec - jsc + 1)
        if (.not. allocated(fld3d_un)) allocate(fld3d_un(ns))
        fld3d_un = reshape(fld3d(isc:iec,jsc:jec,ival), (/ns/))
        call self%horiz_interp%apply(fld3d_un(1:ns), gom_window)
      end if

      ! Fill proper geoval according to time window
      do indx = 1, locs%nlocs
        if (time_mask(indx)) then
          geovals%geovals(ivar)%vals(ival, indx) = gom_window(indx)
        end if
      end do
    end do
print *,'*********************',geovals%geovals(ivar)%vals, trim(geovals%variables(ivar))
    ! Deallocate temporary arrays
    deallocate(fld3d_un)
    deallocate(fld3d)
    deallocate(gom_window)
  end do

  ! If we reach this point, geovals has been initialized
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
