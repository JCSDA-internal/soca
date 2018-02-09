module soca_interph_mod

  use kinds
  
  implicit none
  private

  type, public :: soca_hinterp
     integer, allocatable, dimension(:,:)  :: index
     integer :: nobs
     logical :: initialized
     logical :: alloc
   contains
     procedure :: interp_init
     procedure :: interp_compute_weight
     procedure :: interp_apply     
     procedure :: interp_exit
  end type soca_hinterp

contains

  !--------------------------------------------
  subroutine interp_init(self, nobs)

    implicit none

    integer, intent(in) :: nobs
    class(soca_hinterp), intent(out) :: self
    real(kind=kind_real), allocatable, dimension(:) :: mlon, mlat    

    self%nobs = nobs
    allocate(self%index(nobs,2))
    self%initialized = .false.
    self%alloc = .true.
    
  end subroutine interp_init

  !--------------------------------------------  
  subroutine interp_compute_weight(self, lon, lat, lono, lato)

    use kinds

    implicit none
    
    class(soca_hinterp), intent(inout) :: self    
    real(kind=kind_real), dimension(:,:), intent(in) :: lon, lat
    real(kind=kind_real), dimension(:), intent(in) :: lono, lato
    real(kind=kind_real), allocatable, dimension(:,:) :: wrk    
    integer :: nobs, ni, nj, k, ij(2), cnt

    ni = size(lon,1)
    nj = size(lon,2)

    allocate(wrk(ni,nj))

    ! A very dumb nearest interp
    !$OMP PARALLEL DO
    cnt = 0
    do k = 1, self%nobs
       if (lono(k).gt.-360.0) then
          wrk = abs(lon - lono(k)) + abs(lat - lato(k))
          self%index(k,:) = minloc(wrk)
          cnt = cnt + 1
       else
          self%index(k,:) = 0
       end if
    end do
    !$OMP END PARALLEL DO

    self%initialized = .false.
    
  end subroutine interp_compute_weight

  !--------------------------------------------  
  subroutine interp_apply(self, fld, obs)

    use kinds
    class(soca_hinterp), intent(in) :: self    
    real(kind=kind_real), dimension(:,:), intent(in) :: fld
    real(kind=kind_real), dimension(:), intent(out) :: obs    
    integer :: k

    do k = 1, self%nobs
       obs(k) = fld(self%index(k,1),self%index(k,2))
    end do
    
  end subroutine interp_apply
  
  !--------------------------------------------  
  subroutine interp_exit(self)

    implicit none
    
    class(soca_hinterp), intent(out) :: self

    deallocate(self%index)
    self%initialized = .false.
    self%alloc = .true.
    self%nobs = 0
    
  end subroutine interp_exit

end module soca_interph_mod

