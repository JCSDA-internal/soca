! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_fieldsutils_mod
  use fckit_configuration_module, only: fckit_configuration
  use datetime_mod
  use duration_mod
  use iso_c_binding
  use kinds

  implicit none

  private
  public :: fldinfo, soca_genfilename

  interface fldinfo
     module procedure fldinfo3d, fldinfo2d
  end interface fldinfo

contains

  ! ------------------------------------------------------------------------------

  subroutine fldinfo3d(fld, info)
    real(kind=kind_real),  intent(in) :: fld(:,:,:)
    real(kind=kind_real), intent(out) :: info(3)

    info(1) = minval(fld)
    info(2) = maxval(fld)
    info(3) = sum(fld)/size(fld,3)

  end subroutine fldinfo3d

  ! ------------------------------------------------------------------------------

  subroutine fldinfo2d(fld, info)
    real(kind=kind_real),  intent(in) :: fld(:,:)
    real(kind=kind_real), intent(out) :: info(3)

    info(1) = minval(fld)
    info(2) = maxval(fld)
    info(3) = sum(fld)

  end subroutine fldinfo2d

  ! ------------------------------------------------------------------------------
  !> Generate filename (based on oops/qg)
  function soca_genfilename (c_conf,length,vdate, domain_type)
    type(c_ptr),                intent(in) :: c_conf
    integer,                    intent(in) :: length
    type(datetime),             intent(in) :: vdate
    character(len=3), optional, intent(in) :: domain_type

    character(len=length)                  :: soca_genfilename
    character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
         & prefix, mmb
    type(datetime) :: rdate
    type(duration) :: step
    integer lenfn
    type(fckit_configuration) :: f_conf

    f_conf = fckit_configuration(c_conf)

    if ( f_conf%has("datadir") ) &
        call f_conf%get_or_die("datadir", fdbdir)
    if ( f_conf%has("exp") ) &
        call f_conf%get_or_die("exp", expver)
    if ( f_conf%has("type") ) &
        call f_conf%get_or_die("type", typ)

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
       soca_genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
    endif

    if (typ=="an") then
       call datetime_to_string(vdate, validitydate)
       lenfn = lenfn + 1 + LEN_TRIM(validitydate)
       soca_genfilename = TRIM(prefix) // "." // TRIM(validitydate)
    endif

    if (typ=="incr") then
       soca_genfilename = 'test-incr.nc'
    endif

    if (lenfn>length) &
         & call abor1_ftn("fields:genfilename: filename too long")

  end function soca_genfilename

end module soca_fieldsutils_mod
