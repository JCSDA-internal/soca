! (C) Copyright 2020-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_getvalues_mod

use soca_geom_mod, only: soca_geom
use soca_fields_mod, only: soca_fields, soca_field
use datetime_mod, only: datetime
use kinds, only: kind_real
use ufo_geovals_mod, only: ufo_geovals
use ufo_locations_mod
use unstructured_interpolation_mod, only: unstrc_interp
use fckit_log_module, only : fckit_log
use iso_c_binding

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
  procedure :: fill_geovals_ad=> soca_getvalues_fillgeovals_ad

end type

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine soca_getvalues_create(self, geom, locs)
  class(soca_getvalues), intent(inout) :: self
  type(soca_geom),          intent(in) :: geom
  type(ufo_locations),      intent(in) :: locs

  integer :: isc, iec, jsc, jec, ni, nj
  integer :: nn, ngrid_in, ngrid_out
  character(len=8) :: wtype = 'barycent'
  real(kind=kind_real), allocatable :: lats_in(:), lons_in(:)
  real(8), allocatable, dimension(:) :: locs_lons, locs_lats

  ! Indices for compute domain (no halo)
  isc = geom%isc ; iec = geom%iec
  jsc = geom%jsc ; jec = geom%jec

  ! get location lat/lons
  ngrid_out = locs%nlocs()
  allocate(locs_lons(ngrid_out), locs_lats(ngrid_out))
  call locs%get_lons(locs_lons)
  call locs%get_lats(locs_lats)

  ! create interpolation weights for fields that do NOT use the land mask
  nn = 3
  ni = iec - isc + 1
  nj = jec - jsc + 1
  ngrid_in = ni * nj
  allocate(lats_in(ngrid_in), lons_in(ngrid_in))

  lons_in = reshape(geom%lon(isc:iec,jsc:jec), (/ngrid_in/))
  lats_in = reshape(geom%lat(isc:iec,jsc:jec), (/ngrid_in/))
  call self%horiz_interp%create(geom%f_comm, nn, wtype, &
                           ngrid_in, lats_in, lons_in, &
                           ngrid_out, locs_lats, locs_lons)

  ! create interpolation weights for fields that DO use the land mask
  deallocate(lats_in, lons_in)
  ngrid_in = count(geom%mask2d(isc:iec,jsc:jec) > 0)
  lons_in = pack(geom%lon(isc:iec,jsc:jec), mask=geom%mask2d(isc:iec,jsc:jec) > 0)
  lats_in = pack(geom%lat(isc:iec,jsc:jec), mask=geom%mask2d(isc:iec,jsc:jec) > 0)
  call self%horiz_interp_masked%create(geom%f_comm, nn, wtype, &
                           ngrid_in, lats_in, lons_in, &
                           ngrid_out, locs_lats, locs_lons)
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
  class(soca_fields),       intent(in) :: fld
  type(datetime),           intent(in) :: t1
  type(datetime),           intent(in) :: t2
  type(ufo_locations),      intent(in) :: locs
  type(ufo_geovals),     intent(inout) :: geovals

  logical(c_bool), allocatable :: time_mask(:)
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
  allocate(time_mask(locs%nlocs()))
  call locs%get_timemask(t1,t2,time_mask)

  ! Allocate temporary geoval and 3d field for the current time window
  do ivar = 1, geovals%nvar

    call fld%get(geovals%variables(ivar), fldptr)
    if (fldptr%metadata%dummy_atm) cycle ! TODO remove this hack
    nval = fldptr%nz

    ! Allocate geovals
    if (.not. geovals%linit) then
      geovals%geovals(ivar)%nval = nval
      allocate(geovals%geovals(ivar)%vals(nval, geovals%geovals(ivar)%nlocs))
      geovals%geovals(ivar)%vals = 0.0_kind_real
    end if

    ! Return if no observations
    if ( geovals%geovals(ivar)%nlocs == 0 ) return

    allocate(gom_window(locs%nlocs()))
    allocate(fld3d(isc:iec,jsc:jec,1:nval))

    masked = fldptr%metadata%masked
    fld3d = fldptr%val(isc:iec,jsc:jec,1:nval)

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
      do indx = 1, locs%nlocs()
        if (time_mask(indx)) then
          geovals%geovals(ivar)%vals(ival, indx) = gom_window(indx)
        end if
      end do
    end do

    ! Deallocate temporary arrays
    deallocate(fld3d_un)
    deallocate(fld3d)
    deallocate(gom_window)
  end do

  ! If we reach this point, geovals has been initialized
  geovals%linit = .true.
end subroutine soca_getvalues_fillgeovals

!------------------------------------------------------------------------------
subroutine soca_getvalues_fillgeovals_ad(self, geom, incr, t1, t2, locs, geovals)
  class(soca_getvalues), intent(inout) :: self
  type(soca_geom),          intent(in) :: geom
  class(soca_fields),    intent(inout) :: incr
  type(datetime),           intent(in) :: t1
  type(datetime),           intent(in) :: t2
  type(ufo_locations),      intent(in) :: locs
  type(ufo_geovals),     intent(in) :: geovals

  logical(c_bool), allocatable :: time_mask(:)
  logical :: masked
  integer :: ivar, nlocs, n
  integer :: ival, nval, indx
    integer :: ni, nj
  integer :: isc, iec, jsc, jec
  integer :: ns
  real(kind=kind_real), allocatable :: gom_window(:,:)
    real(kind=kind_real), allocatable :: gom_window_ival(:)
  real(kind=kind_real), allocatable :: incr3d(:,:,:), incr3d_un(:)
  type(soca_field), pointer :: field
    logical :: found

  ! Indices for compute domain (no halo)
  isc = geom%isc ; iec = geom%iec
  jsc = geom%jsc ; jec = geom%jec

  ! Get mask for locations in this time window
  allocate(time_mask(locs%nlocs()))
  call locs%get_timemask(t1,t2,time_mask)
  allocate(gom_window_ival(locs%nlocs()))

  do ivar = 1, geovals%nvar
    call incr%get(geovals%variables(ivar), field)
    nval = field%nz

    ! Allocate temporary geoval and 3d field for the current time window
    allocate(gom_window(nval,locs%nlocs()))
    allocate(incr3d(isc:iec,jsc:jec,1:nval))
    incr3d = 0.0_kind_real
    gom_window = 0.0_kind_real

    ! determine if this variable should use the masked grid
    masked = field%metadata%masked

    ! Apply backward interpolation: Obs ---> Model
    if (masked) then
      ns = count(incr%geom%mask2d(isc:iec,jsc:jec) > 0)
    else
      ni = iec - isc + 1
      nj = jec - jsc + 1
      ns = ni * nj
    end if
    if (.not.allocated(incr3d_un)) allocate(incr3d_un(ns))

    do ival = 1, nval
      ! Fill proper geoval according to time window
      do indx = 1, locs%nlocs()
        if (time_mask(indx)) then
          gom_window(ival, indx) = geovals%geovals(ivar)%vals(ival, indx)
        end if
      end do
      gom_window_ival = gom_window(ival,1:locs%nlocs())

      if (masked) then
        incr3d_un = pack(incr3d(isc:iec,jsc:jec,ival), mask = incr%geom%mask2d(isc:iec,jsc:jec) >0)
        call self%horiz_interp_masked%apply_ad(incr3d_un, gom_window_ival)
        incr3d(isc:iec,jsc:jec,ival) = unpack(incr3d_un, &
              mask = incr%geom%mask2d(isc:iec,jsc:jec) >0, &
              field = incr3d(isc:iec,jsc:jec,ival))
      else
        incr3d_un = reshape(incr3d(isc:iec,jsc:jec,ival), (/ns/))
        call self%horiz_interp%apply_ad(incr3d_un(1:ns), gom_window_ival)
        incr3d(isc:iec,jsc:jec,ival) = reshape(incr3d_un(1:ns),(/ni,nj/))
      end if
    end do

    field%val(isc:iec, jsc:jec, 1:nval) = field%val(isc:iec, jsc:jec, 1:nval) + &
      incr3d(isc:iec, jsc:jec, 1:nval)

    ! Deallocate temporary arrays
    deallocate(incr3d)
    deallocate(gom_window)

  end do

  deallocate(gom_window_ival)

end subroutine soca_getvalues_fillgeovals_ad

end module soca_getvalues_mod
