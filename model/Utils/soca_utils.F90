! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_utils

  implicit none

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

end module soca_utils
