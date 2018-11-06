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
  use ioda_locs_mod
  use fckit_log_module, only : fckit_log  
  use soca_getvaltraj_mod
  use soca_bumpinterp2d_mod
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices
  use kinds

  implicit none
  private
  public :: getvalues, getvalues_ad

  interface getvalues
     procedure getvalues_traj, getvalues_notraj
  end interface getvalues
  
contains
  ! ------------------------------------------------------------------------------
  subroutine initialize_interph(fld, locs, horiz_interp, bumpid)    
  
    implicit none

    type(soca_field),         intent(in) :: fld
    type(ioda_locs),          intent(in) :: locs
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
  subroutine getvalues_traj(fld, locs, vars, geovals, traj)

    implicit none

    type(soca_field),              intent(inout) :: fld
    type(ioda_locs),                  intent(in) :: locs
    type(ufo_vars),                   intent(in) :: vars    
    type(ufo_geovals),             intent(inout) :: geovals
    type(soca_getvaltraj), target, intent(inout) :: traj    

    integer, save :: bumpid = 1000

    call check(fld)    
    if (.not.(traj%interph_initialized)) then
       traj%bumpid = bumpid
       call initialize_interph(fld, locs, traj%horiz_interp, traj%bumpid)
       call traj%horiz_interp%info()
       traj%interph_initialized = .true.
       bumpid = bumpid + 1
    end if
    call interp_tl(fld, locs, vars, geovals, traj%horiz_interp)    

  end subroutine getvalues_traj

  ! ------------------------------------------------------------------------------
  subroutine getvalues_notraj(fld, locs, vars, geovals)

    implicit none

    type(soca_field),   intent(inout) :: fld
    type(ioda_locs),       intent(in) :: locs
    type(ufo_vars),        intent(in) :: vars    
    type(ufo_geovals),  intent(inout) :: geovals

    type(soca_bumpinterp2d) :: horiz_interp    
    integer, save :: bumpid = 2000
    
    call check(fld)    
    call initialize_interph(fld, locs, horiz_interp, bumpid)
    call interp_tl(fld, locs, vars, geovals, horiz_interp)    
    bumpid = bumpid + 1
    
  end subroutine getvalues_notraj
  
  ! ------------------------------------------------------------------------------

  subroutine getvalues_ad(fld, locs, vars, geovals, traj)
    
    implicit none

    type(soca_field),              intent(inout) :: fld
    type(ioda_locs),                  intent(in) :: locs
    type(ufo_vars),                   intent(in) :: vars        
    type(ufo_geovals),             intent(inout) :: geovals
    type(soca_getvaltraj), target, intent(inout) :: traj    

    type(soca_bumpinterp2d), pointer :: horiz_interp_p
    integer :: icat, ilev, ivar, nobs, nval    
    character(len=160) :: record
    integer :: isc, iec, jsc, jec


    horiz_interp_p => traj%horiz_interp

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(fld%geom%ocean, "compute", isc, iec, jsc, jec)

    do ivar = 1, vars%nv
       write(record,*) "getvalues_ad: ",trim(vars%fldnames(ivar))
       call fckit_log%info(record)

       select case (trim(vars%fldnames(ivar)))
       case ("ice_concentration")
          do icat = 1,fld%geom%ocean%ncat
             call horiz_interp_p%applyad(fld%cicen(isc:iec,jsc:jec,icat+1), geovals%geovals(ivar)%vals(icat,:))
          enddo

       case ("ice_thickness")
          do icat = 1,fld%geom%ocean%ncat
             call horiz_interp_p%applyad(fld%hicen(isc:iec,jsc:jec,icat), geovals%geovals(ivar)%vals(icat,:))
          enddo

       case ("sea_surface_height_above_geoid","steric_height") !!!! steric height sould be  different case
          call horiz_interp_p%applyad(fld%ssh(isc:iec,jsc:jec), geovals%geovals(ivar)%vals(1,:))

       case ("ocean_potential_temperature")
          do ilev = 1, fld%geom%ocean%nzo
             call horiz_interp_p%applyad(fld%tocn(isc:iec,jsc:jec,ilev), geovals%geovals(ivar)%vals(ilev,:))
          end do

       case ("ocean_salinity")
          do ilev = 1, fld%geom%ocean%nzo
             call horiz_interp_p%applyad(fld%socn(isc:iec,jsc:jec,ilev), geovals%geovals(ivar)%vals(ilev,:))
          end do

       case ("ocean_upper_level_temperature")
          call horiz_interp_p%applyad(fld%tocn(isc:iec,jsc:jec,1), geovals%geovals(ivar)%vals(1,:))

       end select
    end do

  end subroutine getvalues_ad

  ! ------------------------------------------------------------------------------
  subroutine interp_tl(fld, locs, vars, geovals, horiz_interp)

    implicit none

    type(soca_field),      intent(inout) :: fld
    type(ioda_locs),          intent(in) :: locs
    type(ufo_vars),           intent(in) :: vars    
    type(ufo_geovals),     intent(inout) :: geovals
    type(soca_bumpinterp2d), intent(inout)  :: horiz_interp

    integer :: icat, ilev, ivar, nobs, nval    
    character(len=160) :: record
    integer :: isc, iec, jsc, jec

    call check(fld)    

    nobs = locs%nlocs

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(fld%geom%ocean, "compute", isc, iec, jsc, jec)

    do ivar = 1, vars%nv
       write(record,*) "interp_tl: ",trim(vars%fldnames(ivar))
       call fckit_log%info(record)

       select case (trim(vars%fldnames(ivar)))
       case ("ice_concentration","ice_thickness")
          nval = fld%geom%ocean%ncat

       case ("steric_height",&
            &"sea_surface_height_above_geoid",&
            &"ocean_upper_level_temperature")
          nval = 1

       case ("ocean_potential_temperature","ocean_salinity")
          nval = fld%geom%ocean%nzo

       case default
          write(record,*) "interp_tl: Doing nothing "
          call fckit_log%info(record)          

       end select

       ! Allocate GeoVaLs (fields at locations)
       if (nval.eq.0) call abor1_ftn("Wrong nval: nval = 0")
       geovals%geovals(ivar)%nval = nval
       if (allocated(geovals%geovals(ivar)%vals)) then
          deallocate(geovals%geovals(ivar)%vals)
       end if
       allocate(geovals%geovals(ivar)%vals(nval,nobs))
       geovals%geovals(ivar)%vals(:,:)=0.0_kind_real    
       geovals%lalloc = .true.       
       geovals%linit = .true.

       select case (trim(vars%fldnames(ivar)))

       case ("ice_concentration")
          do icat = 1,fld%geom%ocean%ncat
             call horiz_interp%apply(fld%cicen(isc:iec,jsc:jec,icat+1)*fld%geom%ocean%mask2d(isc:iec,jsc:jec),&
                  &geovals%geovals(ivar)%vals(icat,:))
          end do
       case ("ice_thickness")
          do icat = 1,fld%geom%ocean%ncat
             call horiz_interp%apply(fld%hicen(isc:iec,jsc:jec,icat)*fld%geom%ocean%mask2d(isc:iec,jsc:jec),&
                  &geovals%geovals(ivar)%vals(icat,:))
          end do

       case ("sea_surface_height_above_geoid","steric_height")
          call horiz_interp%apply(fld%ssh(isc:iec,jsc:jec)*fld%geom%ocean%mask2d(isc:iec,jsc:jec),&
               &geovals%geovals(ivar)%vals(1,:))

       case ("ocean_potential_temperature")
          do ilev = 1, fld%geom%ocean%nzo
             call horiz_interp%apply(fld%tocn(isc:iec,jsc:jec,ilev), geovals%geovals(ivar)%vals(ilev,:))
          end do

       case ("ocean_salinity")
          do ilev = 1, fld%geom%ocean%nzo
             call horiz_interp%apply(fld%socn(isc:iec,jsc:jec,ilev), geovals%geovals(ivar)%vals(ilev,:))
          end do

       case ("ocean_layer_thickness")
          do ilev = 1, fld%geom%ocean%nzo
             call horiz_interp%apply(fld%hocn(isc:iec,jsc:jec,ilev), geovals%geovals(ivar)%vals(ilev,:))
          end do

       case ("ocean_upper_level_temperature")          
          call horiz_interp%apply(fld%tocn(isc:iec,jsc:jec,1), geovals%geovals(ivar)%vals(1,:))

       end select

    end do

  end subroutine interp_tl

  
end module soca_interpfields_mod
