! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Handle fields for the  model

module soca_fields_io_mod

  use config_mod
  use soca_geom_mod
  use soca_fields
  use ufo_vars_mod
  use soca_bumpinterp2d_mod
  use soca_getvaltraj_mod
  use kinds
  use iso_c_binding
  use fms_mod,          only: read_data, write_data, set_domain
  use fms_io_mod,       only : fms_io_init, fms_io_exit,&
       &register_restart_field, restart_file_type,&
       &restore_state, query_initialized,&
       &free_restart_type, save_restart
  use datetime_mod
  use duration_mod  
  use fckit_log_module, only : log

  implicit none

  private

  public :: soca_fld2file, soca_write_restart, soca_genfilename

contains

  ! ------------------------------------------------------------------------------
  !> Save soca fields to file using fms write_data
  subroutine soca_fld2file(fld, filename)

    implicit none
    type(soca_field), intent(in) :: fld    !< Fields
    character(len=800), intent(in) :: filename
    integer :: ii
    character(len=1024):: buf

    call check(fld)

    call fms_io_init()
    call set_domain( fld%geom%ocean%G%Domain%mpp_domain )
    do ii = 1, fld%nf
       select case(fld%fldnames(ii))

       case ('ssh')
          call write_data( filename, "ssh", fld%ssh, fld%geom%ocean%G%Domain%mpp_domain)
          call write_data( filename, "rossby_radius", fld%geom%ocean%rossby_radius, fld%geom%ocean%G%Domain%mpp_domain)
       case ('tocn')
          call write_data( filename, "temp", fld%tocn, fld%geom%ocean%G%Domain%mpp_domain)
       case ('socn')
          call write_data( filename, "salt", fld%socn, fld%geom%ocean%G%Domain%mpp_domain)
       case ('hocn')
          call write_data( filename, "h", fld%hocn, fld%geom%ocean%G%Domain%mpp_domain)
       case ('hicen')
          call write_data( filename, "hicen", fld%hicen, fld%geom%ocean%G%Domain%mpp_domain)
       case ('hsnon')
          call write_data(filename, "hsnon", fld%hsnon, fld%geom%ocean%G%Domain%mpp_domain)
       case ('cicen')
          call write_data(filename, "cicen", fld%cicen, fld%geom%ocean%G%Domain%mpp_domain)
       case ('qicnk')
          call write_data(filename, "qicnk", fld%qicnk, fld%geom%ocean%G%Domain%mpp_domain)
       case ('sicnk')
          call write_data(filename, "sicnk", fld%sicnk, fld%geom%ocean%G%Domain%mpp_domain)
       case ('qsnon')
          call write_data(filename, "qsnon", fld%qsnon, fld%geom%ocean%G%Domain%mpp_domain)
       case ('tsfcn')
          call write_data(filename, "tsfcn", fld%tsfcn, fld%geom%ocean%G%Domain%mpp_domain)
       case default

       end select

    end do
    call fms_io_exit()

  end subroutine soca_fld2file

  ! ------------------------------------------------------------------------------
  !> Save soca fields in a restart format  
  subroutine soca_write_restart(fld, c_conf, vdate)

    implicit none

    type(soca_field), intent(inout) :: fld      !< Fields
    type(c_ptr),         intent(in) :: c_conf   !< Configuration
    type(datetime),   intent(inout) :: vdate    !< DateTime

    integer, parameter :: max_string_length=800
    character(len=max_string_length) :: ocn_filename
    character(len=max_string_length) :: ice_filename, basename, incr_filename    
    character(len=20) :: sdate
    character(len=1024)  :: buf
    integer :: iread, ii

    type(restart_file_type) :: ice_restart
    type(restart_file_type) :: ocean_restart
    integer :: idr, idr_ocean

    integer            :: nobs, nval, pe, ierror


    ! Generate file names
    ocn_filename = genfilename(c_conf,max_string_length,vdate,"ocn")
    ice_filename = genfilename(c_conf,max_string_length,vdate,"ice")

    call fms_io_init()
    ! Ocean State
    idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'ave_ssh', fld%ssh(:,:), &
         domain=fld%geom%ocean%G%Domain%mpp_domain)
    idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'Temp', fld%tocn(:,:,:), &
         domain=fld%geom%ocean%G%Domain%mpp_domain)             
    idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'Salt', fld%socn(:,:,:), &
         domain=fld%geom%ocean%G%Domain%mpp_domain)             
    idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'h', fld%hocn(:,:,:), &
         domain=fld%geom%ocean%G%Domain%mpp_domain)             

    ! Sea-Ice
    idr = register_restart_field(ice_restart, ice_filename, 'part_size', fld%AOGCM%Ice%part_size, &
         domain=fld%geom%ocean%G%Domain%mpp_domain)
    idr = register_restart_field(ice_restart, ice_filename, 'h_ice', fld%AOGCM%Ice%h_ice, &
         domain=fld%geom%ocean%G%Domain%mpp_domain)
    idr = register_restart_field(ice_restart, ice_filename, 'h_snow', fld%AOGCM%Ice%h_snow, &
         domain=fld%geom%ocean%G%Domain%mpp_domain)
    idr = register_restart_field(ice_restart, ice_filename, 'enth_ice', fld%AOGCM%Ice%enth_ice, &
         domain=fld%geom%ocean%G%Domain%mpp_domain)
    idr = register_restart_field(ice_restart, ice_filename, 'enth_snow', fld%AOGCM%Ice%enth_snow, &
         domain=fld%geom%ocean%G%Domain%mpp_domain)
    idr = register_restart_field(ice_restart, ice_filename, 'T_skin', fld%AOGCM%Ice%T_skin, &
         domain=fld%geom%ocean%G%Domain%mpp_domain)
    idr = register_restart_field(ice_restart, ice_filename, 'sal_ice', fld%AOGCM%Ice%sal_ice, &
         domain=fld%geom%ocean%G%Domain%mpp_domain)

    call save_restart(ocean_restart, directory='')
    call save_restart(ice_restart, directory='')
    call free_restart_type(ice_restart)
    call free_restart_type(ocean_restart)
    call fms_io_exit()

    return

  end subroutine soca_write_restart

  ! ------------------------------------------------------------------------------
  !> Generate filename (based on oops/qg)
  function soca_genfilename (c_conf,length,vdate, domain_type)
    use iso_c_binding
    use datetime_mod
    use duration_mod
    type(c_ptr),                intent(in) :: c_conf
    integer,                    intent(in) :: length
    type(datetime),             intent(in) :: vdate
    character(len=3), optional, intent(in) :: domain_type    
    character(len=length)                  :: genfilename

    character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
         & prefix, mmb
    type(datetime) :: rdate
    type(duration) :: step
    integer lenfn

    ! here we should query the length and then allocate "string".
    ! But Fortran 90 does not allow variable-length allocatable strings.
    ! config_get_string checks the string length and aborts if too short.
    fdbdir = config_get_string(c_conf,len(fdbdir),"datadir")
    expver = config_get_string(c_conf,len(expver),"exp")
    typ    = config_get_string(c_conf,len(typ)   ,"type")

    if (present(domain_type)) then
       expver = trim(domain_type)//"."//expver
    else
       expver = "ocn.ice."//expver
    end if
    if (typ=="ens") then
       mmb = config_get_string(c_conf, len(mmb), "member")
       lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
       prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ) // "." // TRIM(mmb)
    else
       lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
       prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ)
    endif

    if (typ=="fc" .or. typ=="ens") then
       referencedate = config_get_string(c_conf,len(referencedate),"date")
       call datetime_to_string(vdate, validitydate)
       call datetime_create(TRIM(referencedate), rdate)
       call datetime_diff(vdate, rdate, step)
       call duration_to_string(step, sstep)
       lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
       genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
    endif

    if (typ=="an") then
       call datetime_to_string(vdate, validitydate)
       lenfn = lenfn + 1 + LEN_TRIM(validitydate)
       genfilename = TRIM(prefix) // "." // TRIM(validitydate)
    endif

    if (typ=="incr") then
       genfilename = 'test-incr.nc'
    endif

    if (lenfn>length) &
         & call abor1_ftn("fields:genfilename: filename too long")

  end function soca_genfilename

end module soca_fields_io_mod
