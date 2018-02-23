module soca_Butils
  implicit none

contains

  !  Differentiation of gauss in reverse (adjoint) mode:
  !   gradient     of useful results: xctrl
  !   with respect to varying inputs: xctrl xincr
  !   RW status of diff variables: xctrl:in-out xincr:out
  subroutine gauss_b(xincrb, xctrlb, lon, lat, lx, ly)
    use kinds
    implicit none
    real :: xincrb(:, :)
    real, intent(IN) :: lon(:, :), lat(:, :)
    real, intent(IN) :: lx
    real, intent(IN) :: ly
    real :: xctrlb(:, :)
    integer :: i, j, k, l
    integer :: nx, ny
    integer :: dim(2)
    real :: wght
    real :: dlon, dlat
    real :: valb
    integer :: branch
    real :: val
    dim = shape(xctrlb)
    nx = dim(1)
    ny = dim(2)

    wght=0.0
    xincrb = 0.0
    do i=nx,1,-1
       do j=ny,1,-1
          valb = xctrlb(i, j)
          xctrlb(i, j) = 0.0
          do k=nx,1,-1
             do l=ny,1,-1
                xincrb(k, l) = xincrb(k, l) + wght*valb
                dlat = lat(k, l) - lat(i, j)
                dlon = lon(i, j) - lon(k, l)
                if (abs(dlon) .gt. 180.0) dlon = dlon - 360.0
                wght = exp(-((dlon/(2*lx))**2+(dlat/(2*ly))**2))                
             end do
          end do
       end do
    end do
  end subroutine gauss_b

  subroutine gauss(xincr, xctrl, lon, lat, lx, ly)

    implicit none
    real, intent(IN) :: xincr(:, :)
    real, intent(IN) :: lon(:, :), lat(:, :)
    real, intent(IN) :: lx
    real, intent(IN) :: ly
    real, intent(OUT) :: xctrl(:, :)
    integer :: i, j, k, l
    integer :: nx, ny
    integer :: dim(2)
    real :: wght
    real :: dlon, dlat
    intrinsic SHAPE
    intrinsic ABS
    intrinsic EXP
    real :: val
    dim = shape(xincr)
    nx = dim(1)
    ny = dim(2)
    do i=1,nx
       do j=1,ny
          val = 0.0
          do k=1,nx
             do l=1,ny
                !if (xincr(k, l) .ne. 0.0) then
                dlat = lat(k, l) - lat(i, j)
                dlon = lon(i, j) - lon(k, l)
                if (abs(dlon) .gt. 180.0) dlon = dlon - 360.0
                wght = exp(-((dlon/(2*lx))**2+(dlat/(2*ly))**2))
                val = val + xincr(k, l)*wght
                !end if
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
    real(kind=kind_real), allocatable :: dy_b(:,:), Bdy_b(:,:)

    integer :: ii, jj, k, i, j, nx, ny, nk, r
    real(kind=kind_real) :: Lx=5.0, Ly=1.0

    real, allocatable :: xctrl(:,:), xincr(:,:), xincrb(:,:)

    nx=size(dy,1)
    ny=size(dy,2)
    nk=size(dy,2)

    allocate(dy_b(nx, ny), Bdy_b(nx, ny), xctrl(nx, ny), xincr(nx, ny), xincrb(nx, ny))

    Bdy=0.0
    do k=2, 6
       print *,k
       dy_b = 0.0
       Bdy_b = 0.0

       xincrb = 0.0       
       xctrl = real(dy(:,:,k)) 
       call gauss_b(xincrb, xctrl, real(lon), real(lat), real(lx), real(ly))
       xctrl = 0.0
       call gauss(xincrb, xctrl, real(lon), real(lat), real(lx), real(ly))
       Bdy(:,:,k)=xctrl
       print *,xctrl
       !call gauss(dy(:,:,k), Bdy(:,:,k), lon, lat, lx, ly)
    end do

  end subroutine Simple_Bdy
end module soca_Butils
