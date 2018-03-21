!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_Butils
  implicit none

contains

  subroutine gauss(xincr, xctrl, lon, lat, lx, ly)
    use kinds  
    implicit none
    real(kind=kind_real), intent(IN) :: xincr(:, :)
    real(kind=kind_real), intent(IN) :: lon(:, :), lat(:, :)
    real(kind=kind_real), intent(IN) :: lx
    real(kind=kind_real), intent(IN) :: ly
    real(kind=kind_real), intent(OUT) :: xctrl(:, :)
    integer :: i, j, k, l
    integer :: nx, ny
    integer :: dim(2)
    real(kind=kind_real) :: wght
    real(kind=kind_real) :: dlon, dlat
    intrinsic SHAPE
    intrinsic ABS
    intrinsic EXP
    real(kind=kind_real) :: val
    dim = shape(xincr)
    nx = dim(1)
    ny = dim(2)
    do i=1,nx
       do j=1,ny
          val = 0.0
          do k=1,nx
             do l=1,ny
                if (xincr(k, l) .ne. 0.0) then
                   dlat = lat(k, l) - lat(i, j)
                   dlon = lon(i, j) - lon(k, l)
                   if  (abs(dlon) .gt. 180.0) dlon = dlon - 360.0
                   wght = exp(-((dlon/(2*lx))**2+(dlat/(2*ly))**2))
                   val = val + xincr(k, l)*wght
                end if
             end do
          end do
          xctrl(i, j) = val
       end do
    end do
  end subroutine gauss

  subroutine simple_Bdy(Bdy, dy, lon, lat)
    use kinds  
    implicit none

    real(kind=kind_real), intent(in) :: dy(:,:,:)
    real(kind=kind_real), intent(in) :: lon(:,:), lat(:,:)    
    real(kind=kind_real), allocatable, intent(inout) :: Bdy(:,:,:)


    integer :: ii, jj, k, i, j, nx, ny, nk, r
    real(kind=kind_real) :: Lx=5.0, Ly=1.0

    !real, allocatable :: xctrl(:,:), xincr(:,:), xincrb(:,:)

    nx=size(dy,1)
    ny=size(dy,2)
    nk=size(dy,2)

    do k=2, 6
       print *,'category:',k, maxval(dy(:,:,k))
       call gauss(dy(:,:,k), Bdy(:,:,k), lon, lat, lx, ly)
       print *,'Done with category:',k       
    end do

  end subroutine Simple_Bdy
end module soca_Butils
