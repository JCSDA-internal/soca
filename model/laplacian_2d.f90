
!> Horizontal Laplacian operator

!> Nothing fancy here.
!! It's just the standard 5-point finite-difference Laplacian.

subroutine laplacian_2d (x,del2x,nx,ny,deltax,deltay)

!--- The 2d Laplacian

use kinds

implicit none
integer, intent(in) :: nx         !< Zonal grid dimension
integer, intent(in) :: ny         !< Meridional grid dimension
real(kind=kind_real), intent(in)  :: x(nx,ny)     !< Streamfunction
real(kind=kind_real), intent(out) :: del2x(nx,ny) !< Result of applying Laplacian to x
real(kind=kind_real), intent(in) :: deltax        !< Zonal grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: deltay        !< Meridional grid spacing (non-dimensional)

!--- del-squared of the streamfunction (5-point laplacian)

del2x(:,:) = -2.0_kind_real*( 1.0_kind_real/(deltax*deltax) &
                             +1.0_kind_real/(deltay*deltay))*x(:,:)

del2x(1:nx-1,:) = del2x(1:nx-1,:) + (1.0_kind_real/(deltax*deltax))*x(2:nx  ,:)
del2x(nx    ,:) = del2x(nx    ,:) + (1.0_kind_real/(deltax*deltax))*x(1     ,:)
del2x(2:nx  ,:) = del2x(2:nx  ,:) + (1.0_kind_real/(deltax*deltax))*x(1:nx-1,:)
del2x(1     ,:) = del2x(1     ,:) + (1.0_kind_real/(deltax*deltax))*x(nx    ,:)

del2x(:,1:ny-1) = del2x(:,1:ny-1) + (1.0_kind_real/(deltay*deltay))*x(:,2:ny  )
del2x(:,2:ny  ) = del2x(:,2:ny  ) + (1.0_kind_real/(deltay*deltay))*x(:,1:ny-1)

end subroutine laplacian_2d
