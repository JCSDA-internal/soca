! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_fieldsutils_mod

use fckit_configuration_module, only: fckit_configuration
use datetime_mod, only: datetime, &
                        datetime_to_string, datetime_create, datetime_diff
use duration_mod, only: duration, duration_to_string
use kinds, only: kind_real

implicit none

private
public :: fldinfo, soca_genfilename

interface fldinfo
   module procedure fldinfo3d, fldinfo2d
end interface fldinfo

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

subroutine fldinfo3d(fld, mask, info)
  real(kind=kind_real),  intent(in) :: fld(:,:,:)
  logical,               intent(in) :: mask(:,:)
  real(kind=kind_real), intent(out) :: info(3)

  integer :: z
  real(kind=kind_real) :: tmp(3,size(fld, dim=3))

  ! calculate the min/max/sum separately for each masked level
  do z = 1, size(tmp, dim=2)
     tmp(1,z) = minval(fld(:,:,z), mask=mask)
     tmp(2,z) = maxval(fld(:,:,z), mask=mask)
     tmp(3,z) = sum(   fld(:,:,z), mask=mask) / size(fld, dim=3)
  end do

  ! then combine the min/max/sum over all levels
  info(1) = minval(tmp(1,:))
  info(2) = maxval(tmp(2,:))
  info(3) = sum(   tmp(3,:))
end subroutine fldinfo3d

! ------------------------------------------------------------------------------

subroutine fldinfo2d(fld, mask, info)
  real(kind=kind_real),  intent(in) :: fld(:,:)
  logical,               intent(in) :: mask(:,:)
  real(kind=kind_real), intent(out) :: info(3)

  info(1) = minval(fld, mask=mask)
  info(2) = maxval(fld, mask=mask)
  info(3) = sum(   fld, mask=mask)

end subroutine fldinfo2d

! ------------------------------------------------------------------------------
!> Generate filename (based on oops/qg)
function soca_genfilename (f_conf,length,vdate,domain_type)
  type(fckit_configuration),  intent(in) :: f_conf
  integer,                    intent(in) :: length
  type(datetime),             intent(in) :: vdate
  character(len=3), optional, intent(in) :: domain_type

  character(len=length)                  :: soca_genfilename
  character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
       & prefix, mmb
  type(datetime) :: rdate
  type(duration) :: step
  integer lenfn
  character(len=:), allocatable :: str

  call f_conf%get_or_die("datadir", str)
  fdbdir = str
  call f_conf%get_or_die("exp", str)
  expver = str
  call f_conf%get_or_die("type", str)
  typ = str

  if (present(domain_type)) then
     expver = trim(domain_type)//"."//expver
  else
     expver = "ocn.ice."//expver
  end if
  if (typ=="ens") then
     call f_conf%get_or_die("member", str)
     mmb = str
     lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
     prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ) // "." // TRIM(mmb)
  else
     lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
     prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ)
  endif

  if (typ=="fc" .or. typ=="ens") then
     call f_conf%get_or_die("date", str)
     referencedate = str
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

   if ( allocated(str) ) deallocate(str)

end function soca_genfilename

end module soca_fieldsutils_mod
