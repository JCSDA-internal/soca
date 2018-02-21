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

  subroutine naiveGauss (source, filtered, lon, lat)
    use kinds      
    implicit none

    real(kind=kind_real), intent(in)                 :: source(:,:)
    real(kind=kind_real), intent(in)                 :: lon(:,:), lat(:,:)
    real(kind=kind_real), intent(out)                :: filtered(:,:)

    integer                                      :: i, j, k, l
    integer                                      :: nx, ny
    integer                                      :: dim(2)
    real(kind=kind_real)                             :: val, Lx=5.0, Ly=1
    real(kind=kind_real)                             :: dsq         ! distance squared
    real(kind=kind_real)                             :: wght        ! weight
    real(kind=kind_real)                             :: dlon, dlat    
    integer                             :: r=20    

    dim = shape(source)
    nx = dim(1)
    ny = dim(2)

    print *,'in simple B'
    do i = 1, nx
       do j = 1, ny
          val = 0.0
          do k = 1, nx    
             do l = 1, ny
                if (source(k,l).ne.0.0) then
                   dlat = lat(k,l)-lat(i,j)
                   dlon = abs(lon(i,j)-lon(k,l))
                   if (dlon>180.0) dlon=dlon-360.0                   
                   dsq = (dlon/(2*Lx))**2+(dlat/(2*Ly))**2
                   wght = exp(-dsq)
                   val = val + source(k,l) * wght
                end if
             end do
          end do
          filtered(i,j) = val
       end do
    end do

  end subroutine naiveGauss


  subroutine simple_Bdy(Bdy, dy, lon, lat)
    use kinds  
    implicit none

    real(kind=kind_real), intent(in) :: dy(:,:,:)
    real(kind=kind_real), intent(in) :: lon(:,:), lat(:,:)    
    real(kind=kind_real), allocatable, intent(inout) :: Bdy(:,:,:)
    real(kind=kind_real), allocatable :: dy_b(:,:), Bdy_b(:,:)
    
    integer :: ii, jj, k, i, j, nx, ny, nk, r
    real(kind=kind_real) :: dij, lon1, lon2, lat1, lat2, L=100.0e3, wgt, sumwgt
       
    nx=size(dy,1)
    ny=size(dy,2)
    nk=size(dy,2)

    allocate(dy_b(nx, ny), Bdy_b(nx, ny))

    Bdy=0.0
    do k=2, 6
       print *,k
       dy_b = 0.0
       Bdy_b = 0.0       
       call naiveGauss (dy(:,:,k), Bdy(:,:,k), lon, lat)
    end do

  end subroutine Simple_Bdy
end module soca_Butils
