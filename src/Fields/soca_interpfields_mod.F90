!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! ------------------------------------------------------------------------------
!> Interpolation interface for Fields

module soca_interpfields_mod

  use variables_mod, only: oops_vars
  use ufo_geovals_mod
  use ufo_locs_mod
  use soca_getvaltraj_mod
  use soca_bumpinterp2d_mod
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices
  use kinds
  use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum

  implicit none
  private
  public :: getvalues, getvalues_ad

  interface getvalues
     procedure getvalues_traj, getvalues_notraj
  end interface getvalues

contains
  ! ------------------------------------------------------------------------------
  !> Compute interpolation weights
  subroutine initialize_interph(fld, locs, horiz_interp, bumpid)
    type(soca_field),         intent(in) :: fld
    type(ufo_locs),           intent(in) :: locs
    type(soca_bumpinterp2d), intent(out) :: horiz_interp
    integer,                  intent(in) :: bumpid

    integer :: nobs
    integer :: isc, iec, jsc, jec

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(fld%geom%ocean, "compute", isc, iec, jsc, jec)

    ! Compute interpolation weights
    call horiz_interp%initialize(&
            &      fld%geom%ocean%lon(isc:iec,jsc:jec),&
            &      fld%geom%ocean%lat(isc:iec,jsc:jec),&
            &      fld%geom%ocean%mask2d(isc:iec,jsc:jec),&
            &      locs%lon, locs%lat, bumpid)

  end subroutine initialize_interph

  ! ------------------------------------------------------------------------------
  !> Apply forward interpolation (tl or nl)
  subroutine getvalues_traj(fld, locs, vars, geovals, traj, interp_type)
    type(soca_field),      intent(inout) :: fld
    type(ufo_locs),           intent(in) :: locs
    type(oops_vars),          intent(in) :: vars
    type(ufo_geovals),     intent(inout) :: geovals
    type(soca_getvaltraj), intent(inout) :: traj
    character(2),   optional, intent(in) :: interp_type

    integer, save :: bumpid = 1000
    type(fckit_mpi_comm) :: f_comm
    integer :: allpes_nlocs, nlocs
    integer :: isc, iec, jsc, jec

    ! Sanity check for fields
    call check(fld)

    ! Get local obs in [t, t+dt[
    nlocs = locs%nlocs

    ! Get global nlocs in [t, t+dt[
    f_comm = fckit_mpi_comm()
    call f_comm%allreduce(nlocs, allpes_nlocs, fckit_mpi_sum())

    ! Initialize traj and interp
    if (.not.(traj%interph_initialized)) then
       traj%bumpid = bumpid
       traj%nobs = locs%nlocs
       if (traj%nobs>0) traj%noobs = .false.
       call initialize_interph(fld, locs, traj%horiz_interp(1), traj%bumpid)
       !call traj%horiz_interp(1)%info()

       ! Allocate T & S trajectory for the current time slot
       call geom_get_domain_indices(fld%geom%ocean, "compute", isc, iec, jsc, jec)
       allocate(traj%temp(isc:iec,jsc:jec,fld%geom%ocean%nzo))
       allocate(traj%salt(isc:iec,jsc:jec,fld%geom%ocean%nzo))       

       traj%interph_initialized = .true.       
       bumpid = bumpid + 1
    end if

    ! Apply interpolation
    call interp(fld, locs, vars, geovals, traj%horiz_interp(1), interp_type)

  end subroutine getvalues_traj

  ! ------------------------------------------------------------------------------
  !> Apply forward interpolation
  subroutine getvalues_notraj(fld, locs, vars, geovals)
    type(soca_field),   intent(inout) :: fld
    type(ufo_locs),        intent(in) :: locs
    type(oops_vars),       intent(in) :: vars
    type(ufo_geovals),  intent(inout) :: geovals

    type(soca_bumpinterp2d) :: horiz_interp
    integer, save :: bumpid = 2000

    call check(fld)

    call initialize_interph(fld, locs, horiz_interp, bumpid)
    call interp(fld, locs, vars, geovals, horiz_interp, interp_type='nl')
    bumpid = bumpid + 1

  end subroutine getvalues_notraj

  ! ------------------------------------------------------------------------------
  !> Apply backward interpolation
  subroutine getvalues_ad(fld, locs, vars, geovals, traj)
    type(soca_field),              intent(inout) :: fld
    type(ufo_locs),                   intent(in) :: locs
    type(oops_vars),                  intent(in) :: vars
    type(ufo_geovals),             intent(inout) :: geovals
    type(soca_getvaltraj), target, intent(inout) :: traj

    type(soca_bumpinterp2d), pointer :: horiz_interp_p
    integer :: icat, ilev, ivar, nobs, nval, ival, indx
    character(len=160) :: record
    integer :: isc, iec, jsc, jec
    real(kind=kind_real), allocatable :: gom_window(:,:)
    real(kind=kind_real), allocatable :: fld3d(:,:,:)

    horiz_interp_p => traj%horiz_interp(1)

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(fld%geom%ocean, "compute", isc, iec, jsc, jec)

    do ivar = 1, vars%nv
       ! Set number of levels/categories (nval)
       call nlev_from_ufovar(fld, vars, ivar, nval)

       ! Allocate temporary geoval and 3d field for the current time window
       allocate(gom_window(nval,locs%nlocs))
       allocate(fld3d(isc:iec,jsc:jec,1:nval))
       fld3d = 0.0_kind_real

       ! Apply backward interpolation: Obs ---> Model
       do ival = 1, nval
          ! Fill proper geoval according to time window
          do indx = 1, locs%nlocs
             gom_window(ival, indx) = geovals%geovals(ivar)%vals(ival, locs%indx(indx))
          end do
          
          ! Apply adjoint of transform <--- TODO: I think it goes here ...
          select case (trim(vars%fldnames(ivar)))
          case ("sea_water_potential_temperature")
             !TODO: gom = transform_ad(gom)

          case ("sea_water_practical_salinity", "sea_water_salinity")          
             !TODO: gom = transform_ad(gom)
          end select
          
          call horiz_interp_p%applyad(fld3d(:,:,ival), gom_window(ival,1:locs%nlocs))
       end do

       ! Copy fld3d into field
       select case (trim(vars%fldnames(ivar)))
       case ("sea_ice_category_area_fraction")
          fld%cicen(isc:iec,jsc:jec,2:nval+1) = fld%cicen(isc:iec,jsc:jec,2:nval+1) +&
               &fld3d

       case ("sea_ice_category_thickness")
          fld%hicen(isc:iec,jsc:jec,1:nval) = fld%hicen(isc:iec,jsc:jec,1:nval) +&
               &fld3d

       case ("sea_surface_height_above_geoid")
          fld%ssh(isc:iec,jsc:jec) = fld%ssh(isc:iec,jsc:jec) +&
               &fld3d(isc:iec,jsc:jec,1)

       case ("sea_water_potential_temperature")
          fld%tocn(isc:iec,jsc:jec,1:nval) = fld%tocn(isc:iec,jsc:jec,1:nval) +&
               &fld3d

       case ("sea_water_practical_salinity", "sea_water_salinity")
          fld%socn(isc:iec,jsc:jec,1:nval) = fld%socn(isc:iec,jsc:jec,1:nval) +&
               &fld3d

       case ("sea_water_cell_thickness")
          fld%hocn(isc:iec,jsc:jec,1:nval) = fld%hocn(isc:iec,jsc:jec,1:nval) +&
               &fld3d

       case ("sea_surface_temperature")
          fld%tocn(isc:iec,jsc:jec,1) = fld%tocn(isc:iec,jsc:jec,1) +&
               &fld3d(isc:iec,jsc:jec,1)

       case ("sea_surface_salinity")
          fld%socn(isc:iec,jsc:jec,1) = fld%socn(isc:iec,jsc:jec,1) +&
               &fld3d(isc:iec,jsc:jec,1)

       ! Cool skin
       case ("net_downwelling_shortwave_radiation")
          fld%ocnsfc%sw_rad(isc:iec,jsc:jec) = fld%ocnsfc%sw_rad(isc:iec,jsc:jec) +&
               &fld3d(isc:iec,jsc:jec,1)

       case ("net_downwelling_longwave_radiation")
          fld%ocnsfc%lw_rad(isc:iec,jsc:jec) = fld%ocnsfc%lw_rad(isc:iec,jsc:jec) +&
               &fld3d(isc:iec,jsc:jec,1)

       case ("upward_latent_heat_flux_in_air")
          fld%ocnsfc%latent_heat(isc:iec,jsc:jec) = fld%ocnsfc%latent_heat(isc:iec,jsc:jec) +&
               &fld3d(isc:iec,jsc:jec,1)

       case ("upward_sensible_heat_flux_in_air")
          fld%ocnsfc%sens_heat(isc:iec,jsc:jec) = fld%ocnsfc%sens_heat(isc:iec,jsc:jec) +&
               &fld3d(isc:iec,jsc:jec,1)

       case ("friction_velocity_over_water")
          fld%ocnsfc%fric_vel(isc:iec,jsc:jec) = fld%ocnsfc%fric_vel(isc:iec,jsc:jec) +&
               &fld3d(isc:iec,jsc:jec,1)
       end select

       ! Deallocate temporary arrays
       deallocate(fld3d)
       deallocate(gom_window)

    end do

  end subroutine getvalues_ad

  ! ------------------------------------------------------------------------------
  !> Interace to bump forward interpolation
  subroutine interp(fld, locs, vars, geovals, horiz_interp, interp_type)
    type(soca_field),         intent(inout) :: fld
    type(ufo_locs),              intent(in) :: locs
    type(oops_vars),             intent(in) :: vars
    type(ufo_geovals),        intent(inout) :: geovals
    type(soca_bumpinterp2d),  intent(inout) :: horiz_interp
    character(2),      optional, intent(in) :: interp_type

    integer :: icat, ilev, ivar, nlocs, nlocs_window
    integer :: ival, nval, indx
    character(len=160) :: record
    integer :: isc, iec, jsc, jec
    real(kind=kind_real), allocatable :: gom_window(:,:)
    real(kind=kind_real), allocatable :: fld3d(:,:,:)
    integer :: iii(3), bumpid=1111

    ! Sanity check
    call check(fld)

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(fld%geom%ocean, "compute", isc, iec, jsc, jec)

    ! Loop through ufo vars
    do ivar = 1, vars%nv
       ! Set number of levels/categories (nval)
       call nlev_from_ufovar(fld, vars, ivar, nval)

       ! Allocate GeoVaLs (fields at locations)
       geovals%geovals(ivar)%nval = nval
       if (.not.(allocated(geovals%geovals(ivar)%vals))) then
          ! Number of obs in pe
          nlocs = geovals%geovals(ivar)%nlocs

          allocate(geovals%geovals(ivar)%vals(nval,nlocs))
          geovals%linit = .true.
       end if

       ! Allocate temporary geoval and 3d field for the current time window
       allocate(gom_window(nval,locs%nlocs))
       allocate(fld3d(isc:iec,jsc:jec,1:nval))

       ! Extract fld3d from field
       select case (trim(vars%fldnames(ivar)))

       case ("sea_ice_category_area_fraction")
          fld3d = fld%cicen(isc:iec,jsc:jec,2:nval+1)

       case ("sea_ice_category_thickness")
          fld3d = fld%hicen(isc:iec,jsc:jec,1:nval)

       case ("sea_surface_height_above_geoid")
          fld3d(isc:iec,jsc:jec,1) = fld%ssh(isc:iec,jsc:jec)

       case ("sea_water_potential_temperature")
          fld3d = fld%tocn(isc:iec,jsc:jec,1:nval)

       case ("sea_water_practical_salinity", "sea_water_salinity")
          fld3d = fld%socn(isc:iec,jsc:jec,1:nval)

       case ("sea_water_cell_thickness")
          fld3d = fld%hocn(isc:iec,jsc:jec,1:nval)

       case ("sea_surface_temperature")
          fld3d(isc:iec,jsc:jec,1) = fld%tocn(isc:iec,jsc:jec,1)

       case ("sea_surface_salinity")
          fld3d(isc:iec,jsc:jec,1) = fld%socn(isc:iec,jsc:jec,1)

       case ("sea_area_fraction")
          fld3d(isc:iec,jsc:jec,1) = real(fld%geom%ocean%mask2d(isc:iec,jsc:jec),kind=kind_real)

       case ("net_downwelling_shortwave_radiation")
          fld3d(isc:iec,jsc:jec,1) = fld%ocnsfc%sw_rad(isc:iec,jsc:jec)

       case ("upward_latent_heat_flux_in_air")
          fld3d(isc:iec,jsc:jec,1) = fld%ocnsfc%latent_heat(isc:iec,jsc:jec)
          
       case ("upward_sensible_heat_flux_in_air")
          fld3d(isc:iec,jsc:jec,1) = fld%ocnsfc%sens_heat(isc:iec,jsc:jec)

       case ("net_downwelling_longwave_radiation")
          fld3d(isc:iec,jsc:jec,1) = fld%ocnsfc%lw_rad(isc:iec,jsc:jec)

       case ("friction_velocity_over_water")
          fld3d(isc:iec,jsc:jec,1) = fld%ocnsfc%fric_vel(isc:iec,jsc:jec)

       case default
          call abor1_ftn("soca_interpfields_mod: geoval does not exist")
       end select

       ! Apply transform
       select case (trim(vars%fldnames(ivar)))
       case ("sea_water_potential_temperature")
          select case (interp_type)
          case ("nl")
             !TODO: fld3d = transform_nl(fld3d)
          case ("tl")
             !TODO: fld3d = transform_tl(fld3d)
          end select
          
       case ("sea_water_practical_salinity", "sea_water_salinity")          
          select case (interp_type)
          case ("nl")
             !TODO: fld3d = transform_nl(fld3d)
          case ("tl")
             !TODO: fld3d = transform_tl(fld3d)
          end select       
          
       end select
       
       ! Apply forward interpolation: Model ---> Obs
       do ival = 1, nval
          call horiz_interp%apply(fld3d(isc:iec,jsc:jec,ival), gom_window(ival,:))
          ! Fill proper geoval according to time window
          do indx = 1, locs%nlocs
             geovals%geovals(ivar)%vals(ival, locs%indx(indx)) = gom_window(ival, indx)
          end do
       end do

       ! Deallocate temporary arrays
       deallocate(fld3d)
       deallocate(gom_window)

    end do

  end subroutine interp

  ! ------------------------------------------------------------------------------
  !> Get 3rd dimension of fld
  subroutine nlev_from_ufovar(fld, vars, index_vars, nval)
    type(soca_field), intent(in) :: fld
    type(oops_vars),  intent(in) :: vars
    integer,          intent(in) :: index_vars
    integer,         intent(out) :: nval

    ! Get number of levels or categories (nval)
    select case (trim(vars%fldnames(index_vars)))
    case ("sea_ice_category_area_fraction", &
          "sea_ice_category_thickness")
       nval = fld%geom%ocean%ncat

    case ("sea_surface_height_above_geoid", &
          "sea_surface_temperature", &
          "sea_surface_salinity", &
          "sea_area_fraction", &
          "net_downwelling_shortwave_radiation", &
          "upward_latent_heat_flux_in_air", &
          "upward_sensible_heat_flux_in_air", &
          "net_downwelling_longwave_radiation", &
          "friction_velocity_over_water")
       nval = 1

    case ("sea_water_potential_temperature", &
          "sea_water_practical_salinity", &
          "sea_water_salinity", &
          "sea_water_cell_thickness")
       nval = fld%geom%ocean%nzo

    case default
       call abor1_ftn("soca_interpfields_mod: Could not set nval")

    end select

  end subroutine nlev_from_ufovar

end module soca_interpfields_mod
