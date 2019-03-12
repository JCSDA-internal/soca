!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! ------------------------------------------------------------------------------
!> Interpolation interface for Fields

module soca_interpfields_mod

  use ufo_geovals_mod
  use ufo_vars_mod
  use ufo_locs_mod
  use fckit_log_module, only : fckit_log, log
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

    !if (locs%nlocs==0) return

    ! Compute interpolation weights
    call horiz_interp%initialize(&
            &      fld%geom%ocean%lon(isc:iec,jsc:jec),&
            &      fld%geom%ocean%lat(isc:iec,jsc:jec),&
            &      fld%geom%ocean%mask2d(isc:iec,jsc:jec),&
            &      locs%lon, locs%lat, bumpid)

  end subroutine initialize_interph
  
  ! ------------------------------------------------------------------------------
  !> Apply forward linearized interpolation (linearized about traj)
  subroutine getvalues_traj(fld, locs, vars, geovals, traj)
    type(soca_field),      intent(inout) :: fld
    type(ufo_locs),           intent(in) :: locs
    type(ufo_vars),           intent(in) :: vars    
    type(ufo_geovals),     intent(inout) :: geovals
    type(soca_getvaltraj), intent(inout) :: traj    

    integer, save :: bumpid = 1000
    type(fckit_mpi_comm) :: f_comm
    integer :: allpes_nlocs, nlocs

    ! Sanity check for fields
    call check(fld)

    ! Get local obs in [t, t+dt[
    nlocs = locs%nlocs
    
    ! Get global nlocs in [t, t+dt[
    f_comm = fckit_mpi_comm()
    call f_comm%allreduce(nlocs, allpes_nlocs, fckit_mpi_sum())

    ! Return if allpes_nlocs == 0 ???

    ! Initialize traj and interp
    if (.not.(traj%interph_initialized)) then
       traj%bumpid = bumpid
       traj%nobs = locs%nlocs
       if (traj%nobs>0) traj%noobs = .false.
       !if (.not.traj%noobs) return        ! Exit if no obs
       call initialize_interph(fld, locs, traj%horiz_interp(1), traj%bumpid)
       call traj%horiz_interp(1)%info()
       traj%interph_initialized = .true.
       bumpid = bumpid + 1
    end if

    ! Apply interpolation
    call interp_tl(fld, locs, vars, geovals, traj%horiz_interp(1))    

  end subroutine getvalues_traj

  ! ------------------------------------------------------------------------------
  !> Apply (possibly nonlinear) interpolation  
  subroutine getvalues_notraj(fld, locs, vars, geovals)
    type(soca_field),   intent(inout) :: fld
    type(ufo_locs),        intent(in) :: locs
    type(ufo_vars),        intent(in) :: vars    
    type(ufo_geovals),  intent(inout) :: geovals

    type(soca_bumpinterp2d) :: horiz_interp    
    integer, save :: bumpid = 2000

    if (locs%nlocs==0) return
    call check(fld)

    call initialize_interph(fld, locs, horiz_interp, bumpid)
    call interp_tl(fld, locs, vars, geovals, horiz_interp)
    bumpid = bumpid + 1

  end subroutine getvalues_notraj
  
  ! ------------------------------------------------------------------------------
  !> Apply backward interpolation
  subroutine getvalues_ad(fld, locs, vars, geovals, traj)
    type(soca_field),              intent(inout) :: fld
    type(ufo_locs),                   intent(in) :: locs
    type(ufo_vars),                   intent(in) :: vars        
    type(ufo_geovals),             intent(inout) :: geovals
    type(soca_getvaltraj), target, intent(inout) :: traj    

    type(soca_bumpinterp2d), pointer :: horiz_interp_p
    integer :: icat, ilev, ivar, nobs, nval, ival, indx
    character(len=160) :: record
    integer :: isc, iec, jsc, jec
    real(kind=kind_real), allocatable :: gom_window(:,:)
    real(kind=kind_real), allocatable :: fld3d(:,:,:)        

    !if (traj%noobs) return
    
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
          call horiz_interp_p%applyad(fld3d(:,:,ival), gom_window(ival,1:locs%nlocs))
       end do

       ! Copy fld3d into field
       select case (trim(vars%fldnames(ivar)))
       case ("sea_ice_area_fraction")
          fld%cicen(isc:iec,jsc:jec,2:nval+1) = fld%cicen(isc:iec,jsc:jec,2:nval+1) +&
               &fld3d

       case ("sea_ice_thickness")
          fld%hicen(isc:iec,jsc:jec,1:nval) = fld%hicen(isc:iec,jsc:jec,1:nval) +&
               &fld3d

       case ("sea_surface_height_above_geoid")
          fld%ssh(isc:iec,jsc:jec) = fld%ssh(isc:iec,jsc:jec) +&
               &fld3d(isc:iec,jsc:jec,1)

       case ("sea_water_potential_temperature")
          fld%tocn(isc:iec,jsc:jec,1:nval) = fld%tocn(isc:iec,jsc:jec,1:nval) +&
               &fld3d

       case ("sea_water_practical_salinity")
          fld%socn(isc:iec,jsc:jec,1:nval) = fld%socn(isc:iec,jsc:jec,1:nval) +&
               &fld3d

       case ("sea_water_cell_thickness")
          fld%hocn(isc:iec,jsc:jec,1:nval) = fld%hocn(isc:iec,jsc:jec,1:nval) +&
               &fld3d  

       case ("sea_surface_temperature")
          fld%tocn(isc:iec,jsc:jec,1) = fld%tocn(isc:iec,jsc:jec,1) +&
               &fld3d(isc:iec,jsc:jec,1)

       end select

       ! Deallocate temporary arrays
       deallocate(fld3d)
       deallocate(gom_window)

    end do

  end subroutine getvalues_ad

  ! ------------------------------------------------------------------------------
  !> Interace to bump forward interpolation
  subroutine interp_tl(fld, locs, vars, geovals, horiz_interp)
    type(soca_field),         intent(inout) :: fld
    type(ufo_locs),              intent(in) :: locs
    type(ufo_vars),              intent(in) :: vars    
    type(ufo_geovals),        intent(inout) :: geovals
    type(soca_bumpinterp2d),  intent(inout) :: horiz_interp

    integer :: icat, ilev, ivar, nobs, nobs_window
    integer :: ival, nval, indx    
    character(len=160) :: record
    integer :: isc, iec, jsc, jec
    real(kind=kind_real), allocatable :: gom_window(:,:)
    real(kind=kind_real), allocatable :: fld3d(:,:,:)        
    integer :: iii(3), bumpid=1111
    ! Exit if no obs
    !if (locs%nlocs==0) return

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
          nobs = geovals%geovals(ivar)%nobs
          
          allocate(geovals%geovals(ivar)%vals(nval,nobs))
          geovals%lalloc = .true.       
          geovals%linit = .true.
       end if
       
       ! Allocate temporary geoval and 3d field for the current time window
       allocate(gom_window(nval,locs%nlocs))
       allocate(fld3d(isc:iec,jsc:jec,1:nval))

       call log%info(trim(vars%fldnames(ivar)),newl=.true.)

       ! Extract fld3d from field 
       select case (trim(vars%fldnames(ivar)))

       case ("sea_ice_area_fraction")
          fld3d = fld%cicen(isc:iec,jsc:jec,2:nval+1)

       case ("sea_ice_thickness")
          fld3d = fld%hicen(isc:iec,jsc:jec,1:nval)

       case ("sea_surface_height_above_geoid")
          fld3d(isc:iec,jsc:jec,1) = fld%ssh(isc:iec,jsc:jec) 

       case ("sea_water_potential_temperature")
          fld3d = fld%tocn(isc:iec,jsc:jec,1:nval)

       case ("sea_water_practical_salinity")
          fld3d = fld%socn(isc:iec,jsc:jec,1:nval)

       case ("sea_water_cell_thickness")
          fld3d = fld%hocn(isc:iec,jsc:jec,1:nval)

       case ("sea_surface_temperature")
          fld3d(isc:iec,jsc:jec,1) = fld%tocn(isc:iec,jsc:jec,1)

       case ("sea_area_fraction")
          fld3d(isc:iec,jsc:jec,1) = real(fld%geom%ocean%mask2d(isc:iec,jsc:jec),kind=kind_real)

       case default
          call abor1_ftn("soca_interpfields_mod: geoval does not exist")
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

  end subroutine interp_tl

  ! ------------------------------------------------------------------------------
  !> Get 3rd dimension of fld
  subroutine nlev_from_ufovar(fld, vars, index_vars, nval)
    type(soca_field), intent(in) :: fld    
    type(ufo_vars),   intent(in) :: vars
    integer,          intent(in) :: index_vars    
    integer,         intent(out) :: nval
    
    ! Get number of levels or categories (nval)
    select case (trim(vars%fldnames(index_vars)))
    case ("sea_ice_area_fraction", &
          "sea_ice_thickness")
       nval = fld%geom%ocean%ncat

    case ("sea_surface_height_above_geoid", &
          "sea_surface_temperature", &
          "sea_area_fraction")
       nval = 1

    case ("sea_water_potential_temperature", &
          "sea_water_practical_salinity", &
          "sea_water_cell_thickness")
       nval = fld%geom%ocean%nzo
       
    case default
       call abor1_ftn("soca_interpfields_mod: Could not set nval")

    end select

  end subroutine nlev_from_ufovar
  
end module soca_interpfields_mod
