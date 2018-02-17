module soca_Butils
  implicit none

contains

  function to_radian(degree) result(rad)
    ! degrees to radians
    use kinds      
    real(kind=kind_real),intent(in) :: degree
    real(kind=kind_real), parameter :: deg_to_rad = atan(1.0)/45 ! exploit intrinsic atan to generate pi/180 runtime constant
    real(kind=kind_real) :: rad
    
    rad = degree*deg_to_rad
  end function to_radian
 
  function haversine(deglat1,deglon1,deglat2,deglon2) result (dist)
    ! great circle distance -- adapted from Matlab
    use tools_const, only: deg2rad
    use kinds      
    implicit none

    real(kind=kind_real),intent(in) :: deglat1,deglon1,deglat2,deglon2
    real(kind=kind_real) :: a,c,dist,dlat,dlon,lat1,lat2
    real(kind=kind_real),parameter :: radius = 6372.8e3
    
    dlat = deg2rad*(deglat2-deglat1)
    dlon = deg2rad*(deglon2-deglon1)
    lat1 = deg2rad*(deglat1)
    lat2 = deg2rad*(deglat2)
    a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
    c = 2*asin(sqrt(a))
    dist = radius*c
  end function haversine


  subroutine naiveGauss (source, filtered, r, lon, lat)
    use kinds      
    implicit none
    integer, intent(in)                          :: r
    real(kind=kind_real), intent(in)                 :: source(:,:)
    real(kind=kind_real), intent(in)                 :: lon(:,:), lat(:,:)
    real(kind=kind_real), intent(out)                :: filtered(:,:)

    real(kind=kind_real), parameter                  :: PI = 4.*atan(1.)

    integer                                      :: i, j, k, l
    integer                                      :: ii, jj      ! inside the array
    integer                                      :: nx, ny
    integer                                      :: dim(2)
    real(kind=kind_real)                             :: val, wsum, lon1, lon2, lat1, lat2, L0=300.0e3
    real(kind=kind_real)                             :: dsq         ! distance squared
    real(kind=kind_real)                             :: wght        ! weight

    dim = shape(source)
    nx = dim(1)
    ny = dim(2)

    do i = 1, nx
       do j = 1, ny
          val = 0
          wsum = 0
          do k = i-r, i+r     ! inner loop over kernel
             do l = j-r, j+r  ! inner loop over kernel

                ii = min(nx, max(1,k))   ! make sure i is always inside the grid (this implies values are extendet (stretched at the boundaries))
                jj = min(ny, max(1,l))   ! make sure j is always inside the grid (this implies values are extendet (stretched at the boundaries))

                if (source(ii,jj).ne.0.0) then
                   lat1 = lat(i,j)
                   lon1 = lon(i,j)
                   lat2 = lat(ii,jj)
                   lon2 = lon(ii,jj)                                
                   dsq = haversine(lat1, lon1, lat2, lon2)
                   wght = exp(-(dsq / (2.d0*L0))**2) !/ (2.d0*PI*r**2)
                   val = val + source(ii,jj) * wght                
                   wsum = wsum + wght          
                end if
             end do
          end do
          filtered(i,j) = val / wsum
          
       end do
    end do
    
  end subroutine naiveGauss
  
  subroutine simple_Bdy(Bdy, dy, lon, lat)
    use kinds  
    implicit none

    real(kind=kind_real), intent(in) :: dy(:,:,:)
    real(kind=kind_real), intent(in) :: lon(:,:), lat(:,:)    
    real(kind=kind_real), allocatable, intent(inout) :: Bdy(:,:,:)
    integer :: ii, jj, k, i, j, nx, ny, nk, r
    real(kind=kind_real) :: dij, lon1, lon2, lat1, lat2, L=500.0e3, wgt, sumwgt
       
    nx=size(dy,1)
    ny=size(dy,2)
    nk=size(dy,2)
    
    Bdy=0.0

    r=200
    do k=1, 5
       call naiveGauss (dy(:,:,k), Bdy(:,:,k), r, lon, lat)
    end do
    
  end subroutine Simple_Bdy
end module soca_Butils
