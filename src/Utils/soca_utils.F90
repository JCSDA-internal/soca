! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_utils

  implicit none
  private
  public :: inside_polygon, write2pe
contains

  ! ------------------------------------------------------------------------------  
  !
  !> PRIVATE FUNCTION COPIED FROM FMS-GFDL/horiz_interp/horiz_interp_bilinear.F90
  !>
  !> The function will return true if the point x,y is inside a polygon, or
  !> NO if it is not.  If the point is exactly on the edge of a polygon,
  !> the function will return .true.
  !>
  !> real polyx(:) : longitude coordinates of corners 
  !> real polyx(:) : latitude  coordinates of corners
  !> real x,y      : point to be tested 
  !> ??? How to deal with truncation error.
  !>
  function inside_polygon(polyx, polyy, x, y)
    real, dimension(:), intent(in) :: polyx, polyy
    real,               intent(in) :: x, y
    logical                        :: inside_polygon
    integer                        :: i, j, nedges
    real                           :: xx

    inside_polygon = .false.
    nedges = size(polyx(:))
    j = nedges
    do i = 1, nedges
       if( (polyy(i) < y .AND. polyy(j) >= y) .OR. (polyy(j) < y .AND. polyy(i) >= y) ) then
          xx = polyx(i)+(y-polyy(i))/(polyy(j)-polyy(i))*(polyx(j)-polyx(i))
          if( xx == x ) then
             inside_polygon = .true.
             return
          else if( xx < x ) then
             inside_polygon = .not. inside_polygon
          endif
       endif
       j = i
    enddo

    return

  end function inside_polygon
  ! ------------------------------------------------------------------------------

  subroutine write2pe(vec,varname,filename,append)

    use netcdf
    use kinds
    
    implicit none

    real(kind=kind_real), intent(in) :: vec(:)
    character(len=256),   intent(in) :: varname
    character(len=256),   intent(in) :: filename
    logical,              intent(in) :: append

    !netcdf stuff    
    integer(kind=4) :: iNcid
    integer(kind=4) :: iDim_ID
    integer(kind=4) :: iVar_ID
    integer         :: ndims=1, ns

    ns=size(vec)
    if (append) then  ! If file exists, append to it
       call nc_check( nf90_open(filename, NF90_WRITE, iNcid) )
       call nc_check( nf90_inquire(iNcid, nDimensions = ndims) )
       call nc_check( nf90_inq_dimid(iNcid, "ns", iDim_ID) )
       call nc_check( nf90_redef(iNcid) )
    else    
       call nc_check(nf90_create(filename, NF90_CLOBBER, iNcid))
       call nc_check(nf90_def_dim(iNcid, "ns", ns, iDim_ID))
    end if

    ! Define of variables.
    call nc_check( nf90_def_var(iNcid, trim(varname), NF90_DOUBLE,  (/ iDim_ID /), iVar_ID) )

    ! End define mode.
    call nc_check(nf90_enddef(iNcid))

    ! Writing
    call nc_check(nf90_put_var(iNcid, iVar_ID , vec))     

    ! Close file.
    call nc_check(nf90_close(iNcid))
    
  end subroutine write2pe

  ! ------------------------------------------------------------------------------  
  
  subroutine nc_check(status)

    use netcdf
    IMPLICIT NONE
    integer(4), intent ( in) :: status
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine nc_check
  
end module soca_utils
