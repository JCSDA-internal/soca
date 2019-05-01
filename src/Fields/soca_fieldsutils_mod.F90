! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_fieldsutils_mod
  use kinds
  
  implicit none

  private
  public :: fldinfo

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
    info(3) =    sum(fld)/size(fld,3)

  end subroutine fldinfo3d

  ! ------------------------------------------------------------------------------

  subroutine fldinfo2d(fld, info)
    real(kind=kind_real),  intent(in) :: fld(:,:)
    real(kind=kind_real), intent(out) :: info(3)

    info(1) = minval(fld)
    info(2) = maxval(fld)
    info(3) =    sum(fld)

  end subroutine fldinfo2d

end module soca_fieldsutils_mod
