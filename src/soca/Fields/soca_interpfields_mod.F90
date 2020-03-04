! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------
!> Interpolation interface for Fields

module soca_interpfields_mod

use kinds, only: kind_real
use fckit_log_module, only : fckit_log
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
use oops_variables_mod
use ufo_geovals_mod, only: ufo_geovals
use ufo_locs_mod, only: ufo_locs
use soca_getvaltraj_mod, only: soca_getvaltraj
use soca_fields_mod, only: soca_fields, soca_field
use soca_geom_mod, only : soca_geom
use unstructured_interpolation_mod

implicit none
private
public :: getvalues, getvalues_ad

interface getvalues
   procedure getvalues_traj, getvalues_notraj
end interface getvalues

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Compute interpolation weights
subroutine initialize_interph(fld, locs, horiz_interp, horiz_interp_masked)
  type(soca_fields),   intent(in) :: fld
  type(ufo_locs),      intent(in) :: locs
  type(unstrc_interp), intent(out) :: horiz_interp
  type(unstrc_interp), intent(out) :: horiz_interp_masked

  integer :: isc, iec, jsc, jec, ni, nj
  type(fckit_mpi_comm) :: f_comm
  integer :: nn, ngrid_in, ngrid_out
  character(len=8) :: wtype = 'barycent'
  real(kind=kind_real), allocatable :: lats_in(:), lons_in(:)

  ! Indices for compute domain (no halo)
  isc = fld%geom%isc ; iec = fld%geom%iec
  jsc = fld%geom%jsc ; jec = fld%geom%jec

  f_comm = fckit_mpi_comm()
  ngrid_out = locs%nlocs
  nn = 3

  ! create interpolation weights for fields that do NOT use the land mask
  ni = iec - isc + 1
  nj = jec - jsc + 1
  ngrid_in = ni * nj
  allocate(lats_in(ngrid_in), lons_in(ngrid_in))
  lons_in = reshape(fld%geom%lon(isc:iec,jsc:jec), (/ngrid_in/))
  lats_in = reshape(fld%geom%lat(isc:iec,jsc:jec), (/ngrid_in/))
  call horiz_interp%create(f_comm, nn, wtype, &
                           ngrid_in, lats_in, lons_in, &
                           ngrid_out, locs%lat, locs%lon)

  ! create interpolation weights for fields that DO use the land mask
  deallocate(lats_in, lons_in)
  ngrid_in = count(fld%geom%mask2d(isc:iec,jsc:jec) > 0)
  lons_in = pack(fld%geom%lon(isc:iec,jsc:jec), mask=fld%geom%mask2d(isc:iec,jsc:jec) > 0)
  lats_in = pack(fld%geom%lat(isc:iec,jsc:jec), mask=fld%geom%mask2d(isc:iec,jsc:jec) > 0)
  call horiz_interp_masked%create(f_comm, nn, wtype, &
                           ngrid_in, lats_in, lons_in, &
                           ngrid_out, locs%lat, locs%lon)
end subroutine initialize_interph

! ------------------------------------------------------------------------------
!> Apply forward interpolation (tl or nl)
subroutine getvalues_traj(fld, locs, vars, geoval, traj, interp_type)
  type(soca_fields),     intent(inout) :: fld
  type(ufo_locs),           intent(in) :: locs
  type(oops_variables),     intent(in) :: vars
  type(ufo_geovals),     intent(inout) :: geoval
  type(soca_getvaltraj), intent(inout) :: traj
  character(2),   optional, intent(in) :: interp_type

  type(fckit_mpi_comm) :: f_comm
  integer :: allpes_nlocs, nlocs

  ! Get local obs in [t, t+dt[
  nlocs = locs%nlocs

  ! Get global nlocs in [t, t+dt[
  f_comm = fckit_mpi_comm()
  call f_comm%allreduce(nlocs, allpes_nlocs, fckit_mpi_sum())

  ! Initialize traj and interp
  if (.not.(traj%interph_initialized)) then
     traj%nobs = locs%nlocs
     call initialize_interph(fld, locs, traj%horiz_interp, traj%horiz_interp_masked)
     traj%interph_initialized = .true.
  end if

  select case (interp_type)
  case('tl')
     ! Apply interpolation with TL transform
     ! TODO: pass in "traj" and do something with it, at some point?
     call interp(fld, locs, vars, geoval, traj%horiz_interp, traj%horiz_interp_masked)
  case('nl')
     ! Apply interpolation with NL transform
     call interp(fld, locs, vars, geoval, traj%horiz_interp, traj%horiz_interp_masked)
  end select

end subroutine getvalues_traj

! ------------------------------------------------------------------------------
!> Apply forward interpolation
subroutine getvalues_notraj(fld, locs, vars, geoval)
  type(soca_fields),  intent(inout) :: fld
  type(ufo_locs),        intent(in) :: locs
  type(oops_variables),  intent(in) :: vars
  type(ufo_geovals),  intent(inout) :: geoval

  type(unstrc_interp) :: horiz_interp
  type(unstrc_interp) :: horiz_interp_masked

  call initialize_interph(fld, locs, horiz_interp, horiz_interp_masked)

  ! Apply interpolation with NL transform
  call interp(fld, locs, vars, geoval, horiz_interp, horiz_interp_masked)

end subroutine getvalues_notraj

! ------------------------------------------------------------------------------
!> Interace to forward interpolation (NL and TL)
subroutine interp(fld, locs, vars, geoval, horiz_interp, horiz_interp_masked)
  type(soca_fields),       intent(inout) :: fld
  type(ufo_locs),             intent(in) :: locs
  type(oops_variables),       intent(in) :: vars
  type(ufo_geovals),       intent(inout) :: geoval
  type(unstrc_interp),     intent(inout) :: horiz_interp
  type(unstrc_interp),     intent(inout) :: horiz_interp_masked

  logical :: masked
  integer :: ivar, nlocs, n
  integer :: ival, nval, indx
  integer :: isc, iec, jsc, jec
  integer :: ns
  real(kind=kind_real), allocatable :: gom_window(:,:)
  real(kind=kind_real), allocatable :: fld3d(:,:,:), fld3d_un(:)
  type(soca_field), pointer :: fldptr

  ! Indices for compute domain (no halo)
  isc = fld%geom%isc ; iec = fld%geom%iec
  jsc = fld%geom%jsc ; jec = fld%geom%jec

  ! Loop through ufo vars
  do ivar = 1, vars%nvars()
    ! Set number of levels/categories (nval)
    nval = nlev_from_ufovar(fld, vars%variable(ivar))
     if (nval==0) cycle

    ! Allocate GeoVaLs (fields at locations)
    geoval%geovals(ivar)%nval = nval
    if (.not.(allocated(geoval%geovals(ivar)%vals))) then
      ! Number of obs in pe
      nlocs = geoval%geovals(ivar)%nlocs
      allocate(geoval%geovals(ivar)%vals(nval,nlocs))
      geoval%linit = .true.
    end if

    ! Allocate temporary geoval and 3d field for the current time window
    allocate(gom_window(nval,locs%nlocs))
    allocate(fld3d(isc:iec,jsc:jec,1:nval))
    nullify(fldptr)

    ! Extract fld3d from field
    masked = .true. ! by default fields are assumed to need a land mask applied,
                    ! (Currently only sea area fraction is unmasked)

    ! if we are lucky and this variable is a non-derived type, check the fields structure
    do n=1,size(fld%fields)
      if (fld%fields(n)%cf_name == vars%variable(ivar)) then
        fld3d = fld%fields(n)%val(isc:iec,jsc:jec,1:nval)
      end if
    end do

    ! otherwise, we are dealing with a derived type, prepare for a long "select case" statement
    select case (trim(vars%variable(ivar)))
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
          call horiz_interp_masked%apply(fld3d_un, gom_window(ival,:))
        else
          ns = (iec - isc + 1) * (jec - jsc + 1)
          if (.not. allocated(fld3d_un)) allocate(fld3d_un(ns))
          fld3d_un = reshape(fld3d(isc:iec,jsc:jec,ival), (/ns/))
          call horiz_interp%apply(fld3d_un(1:ns), gom_window(ival,:))
        end if

        ! Fill proper geoval according to time window
        do indx = 1, locs%nlocs
           geoval%geovals(ivar)%vals(ival, locs%indx(indx)) = gom_window(ival, indx)
        end do
     end do

     ! Deallocate temporary arrays
     deallocate(fld3d_un)
     deallocate(fld3d)
     deallocate(gom_window)

  end do

end subroutine interp

! ------------------------------------------------------------------------------
!> Apply backward interpolation
subroutine getvalues_ad(incr, locs, vars, geoval, traj)
  type(soca_fields),             intent(inout) :: incr
  type(ufo_locs),                   intent(in) :: locs
  type(oops_variables),             intent(in) :: vars
  type(ufo_geovals),             intent(inout) :: geoval
  type(soca_getvaltraj), target, intent(inout) :: traj

  logical :: masked
  integer :: ivar, nval, ival, indx
  integer :: isc, iec, jsc, jec
  integer :: ni, nj, ns, n
  real(kind=kind_real), allocatable :: gom_window(:,:)
  real(kind=kind_real), allocatable :: incr3d(:,:,:), incr3d_un(:)
  real(kind=kind_real), allocatable :: gom_window_ival(:)
  type(soca_field), pointer :: field
  logical :: found


  ! Indices for compute domain (no halo)
  isc = incr%geom%isc ; iec = incr%geom%iec
  jsc = incr%geom%jsc ; jec = incr%geom%jec

  allocate(gom_window_ival(locs%nlocs))

  do ivar = 1, vars%nvars()
    ! Set number of levels/categories (nval)
    nval = nlev_from_ufovar(incr, vars%variable(ivar))

    ! Allocate temporary geoval and 3d field for the current time window
    allocate(gom_window(nval,locs%nlocs))
    allocate(incr3d(isc:iec,jsc:jec,1:nval))
    incr3d = 0.0_kind_real

    ! determine if this variable should use the masked grid
    ! (currently all of them, perhaps have atm vars use unmasked interp at some point??)
    masked = .true.

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
      do indx = 1, locs%nlocs
        gom_window(ival, indx) = geoval%geovals(ivar)%vals(ival, locs%indx(indx))
      end do
      gom_window_ival = gom_window(ival,1:locs%nlocs)

      if (masked) then
        incr3d_un = pack(incr3d(isc:iec,jsc:jec,ival), mask = incr%geom%mask2d(isc:iec,jsc:jec) >0)
        call traj%horiz_interp_masked%apply_ad(incr3d_un, gom_window_ival)
        incr3d(isc:iec,jsc:jec,ival) = unpack(incr3d_un, &
              mask = incr%geom%mask2d(isc:iec,jsc:jec) >0, &
              field = incr3d(isc:iec,jsc:jec,ival))
      else
        incr3d_un = reshape(incr3d(isc:iec,jsc:jec,ival), (/ns/))
        call traj%horiz_interp%apply_ad(incr3d_un(1:ns), gom_window_ival)
        incr3d(isc:iec,jsc:jec,ival) = reshape(incr3d_un(1:ns),(/ni,nj/))
      end if
    end do

    ! if we are lucky and this variable is a non-derived type, check the fields structure
    found = .false.
    do n=1,size(incr%fields)
      if (incr%fields(n)%cf_name == vars%variable(ivar)) then
        incr%fields(n)%val(isc:iec,jsc:jec,1:nval) = &
          incr%fields(n)%val(isc:iec,jsc:jec,1:nval) + incr3d(isc:iec,jsc:jec,1:nval)
        found = .true.
        exit
      end if
    end do

    ! otherwise, we are dealing with a derived type, prepare for a long "select case" statement
    if (.not. found) then
      select case (trim(vars%variable(ivar)))
      case ("sea_water_salinity")
        call incr%get("socn", field)
        field%val(isc:iec,jsc:jec,1:nval) = field%val(isc:iec,jsc:jec,1:nval) + incr3d

      case ("sea_surface_temperature")
        call incr%get("tocn", field)
        field%val(isc:iec,jsc:jec,1:nval) = field%val(isc:iec,jsc:jec,1:nval) + incr3d

      case ("sea_surface_salinity")
        call incr%get("socn", field)
        field%val(isc:iec,jsc:jec,1) = field%val(isc:iec,jsc:jec,1) + incr3d(isc:iec,jsc:jec,1)

      case default
        call abor1_ftn("soca_interpfields_mod:getvalues_ad geoval does not exist")

      end select
    end if

    ! Deallocate temporary arrays
    deallocate(incr3d)
    deallocate(gom_window)

  end do

  deallocate(gom_window_ival)

end subroutine getvalues_ad

! ------------------------------------------------------------------------------
!> Get 3rd dimension of fld
function nlev_from_ufovar(fld, var) result(nval)
  type(soca_fields), intent(in) :: fld
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

end module soca_interpfields_mod
