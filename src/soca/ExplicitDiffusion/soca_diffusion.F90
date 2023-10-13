! (C) Copyright 2023-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_diffusion_mod

use atlas_module, only: atlas_fieldset, atlas_field, atlas_real
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_min, fckit_mpi_max, fckit_mpi_sum
use kinds, only: kind_real
use logger_mod
use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad
use random_mod
use fms_io_mod

use soca_increment_mod
use soca_geom_mod, only : soca_geom

implicit none
private

! ------------------------------------------------------------------------------
! because I'm TIRED of the unreadability of endless isc, jsc .....
#define DOMAIN                  self%geom%isc:self%geom%iec,self%geom%jsc:self%geom%jec
#define DOMAIN_WITH_HALO        self%geom%isd:self%geom%ied,self%geom%jsd:self%geom%jed
#define LOOP_DOMAIN_I           self%geom%isc, self%geom%iec
#define LOOP_DOMAIN_J           self%geom%jsc, self%geom%jec
#define LOOP_DOMAIN_WITH_HALO_I self%geom%isd, self%geom%ied
#define LOOP_DOMAIN_WITH_HALO_J self%geom%jsd, self%geom%jed

! ------------------------------------------------------------------------------

type, public :: soca_diffusion
 private
  ! grid metrics
  real(kind_real), allocatable :: inv_sqrt_area(:,:) !< 1/sqrt(area)
  real(kind_real), allocatable :: dx(:,:)
  real(kind_real), allocatable :: dy(:,:)
  real(kind_real), allocatable :: pmon_u(:,:) ! pm/pn at u points (pm = 1/dx, pn = 1/dy)
  real(kind_real), allocatable :: pnom_v(:,:) ! pn/pm at v points
  real(kind_real), allocatable :: mask(:,:)   ! 1.0 where water, 0.0 where land
  
  ! parameters calculated during calibration() / read()
  real(kind_real), allocatable :: KhDt(:,:)
  real(kind_real), allocatable :: normalization(:,:)
  integer :: n_iter = -1 !< set to -1 to indicate has not been initialized

  class(soca_geom), pointer :: geom    

contains
  procedure :: init => soca_diffusion_init
  procedure :: calibrate => soca_diffusion_calibrate
  procedure :: multiply => soca_diffusion_multiply
  procedure :: write_params => soca_diffusion_write_params
  procedure :: read_params => soca_diffusion_read_params
  
  procedure, private :: multiply_2D => soca_diffusion_multiply_2D
  procedure, private :: multiply_2D_tl => soca_diffusion_multiply_2D_tl
  procedure, private :: multiply_2D_ad => soca_diffusion_multiply_2D_ad
  procedure, private :: diffusion_steps => soca_diffusion_diffusion_steps
  procedure, private :: diffusion_steps_ad => soca_diffusion_diffusion_steps_ad
  procedure, private :: calc_stats => soca_diffusion_calc_stats
  
  procedure, private :: calc_norm_bruteforce => soca_diffusion_calc_norm_bruteforce
  procedure, private :: calc_norm_randomization => soca_diffusion_calc_norm_randomization
end type soca_diffusion

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! calculate the masked global min/max/mean values for a given field
! TODO: is there any chance we need to also operate on unmasked fields?!?
subroutine soca_diffusion_calc_stats(self, field, stats)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable, intent(in) :: field(:,:)
  real(kind=kind_real),             intent(out) :: stats(3)

  real(kind=kind_real) :: l_min, l_max, l_sum, l_count, g_count

  l_min =  minval(field(DOMAIN), &
                  mask=self%mask(DOMAIN)==1.0)
  l_max =  maxval(field(DOMAIN), &
                  mask=self%mask(DOMAIN)==1.0)
  l_sum =     sum(field(DOMAIN), &
                  mask=self%mask(DOMAIN)==1.0)
  l_count = count(self%mask(DOMAIN)==1.0)

  call self%geom%f_comm%allreduce(l_min, stats(1), fckit_mpi_min())
  call self%geom%f_comm%allreduce(l_max, stats(2), fckit_mpi_max())
  call self%geom%f_comm%allreduce(l_sum, stats(3), fckit_mpi_sum())
  call self%geom%f_comm%allreduce(l_count, g_count, fckit_mpi_sum())
  stats(3) = stats(3) / g_count
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_init(self, geom)
  class(soca_diffusion), intent(inout) :: self
  class(soca_geom), target, intent(in) :: geom

  real(kind=kind_real) :: stats(3) ! min, max, mean
  character(len=1024) :: str  
  integer :: i, j

  call oops_log%trace("soca_diffusion::init() starting", flush=.true.)
  self%geom => geom

  ! grid and derived grid parameters
  !---------------------------------------------------------------------------
  call oops_log%info("ExplicitDiffusion: Initializing grid...")
  
  allocate(self%mask(DOMAIN_WITH_HALO))
  allocate(self%dx(DOMAIN_WITH_HALO))
  allocate(self%dy(DOMAIN_WITH_HALO))
  ! NOTE MOM6 likes to leave unused halo regions undefined, so we explicitly set default values here
  self%mask = 0.0
  self%mask(DOMAIN) = self%geom%mask2d(DOMAIN)
  self%dx = 1.0e-5   
  self%dx(DOMAIN) = geom%dx(DOMAIN)
  self%dy = 1.0e-5
  self%dy(DOMAIN) = geom%dy(DOMAIN)
  call mpp_update_domains(self%mask, self%geom%Domain%mpp_domain, complete=.true.)
  call mpp_update_domains(self%dx, self%geom%Domain%mpp_domain, complete=.true.)
  call mpp_update_domains(self%dy, self%geom%Domain%mpp_domain, complete=.true.)

  call self%calc_stats(self%dx,stats)
    write (str, *) "  dx:   min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)
  call self%calc_stats(self%dy, stats)
    write (str, *) "  dy:   min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)
  call self%calc_stats(self%geom%cell_area, stats)
    write (str, *) "  area: min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)

  ! calculate derived parameters
  allocate(self%pmon_u(DOMAIN_WITH_HALO))
  allocate(self%pnom_v(DOMAIN_WITH_HALO))
  allocate(self%inv_sqrt_area(DOMAIN_WITH_HALO))
  self%pmon_u = 0.0
  self%pnom_v = 0.0
  self%inv_sqrt_area = 0.0

  do j = LOOP_DOMAIN_J
    do i = LOOP_DOMAIN_I
      self%inv_sqrt_area(i,j) = 1.0 / sqrt(self%dx(i,j)*self%dy(i,j))
      self%pmon_u(i,j) = (1.0/self%dx(i-1,j) + 1.0/self%dx(i,j)) / (1.0/self%dy(i-1,j) + 1.0/self%dy(i,j))
      self%pnom_v(i,j) = (1.0/self%dy(i,j-1) + 1.0/self%dy(i,j)) / (1.0/self%dx(i,j-1) + 1.0/self%dx(i,j))
    end do
  end do

  call mpp_update_domains(self%inv_sqrt_area, self%geom%Domain%mpp_domain, complete=.true.)
  call mpp_update_domains(self%pmon_u, self%geom%Domain%mpp_domain, complete=.true.)
  call mpp_update_domains(self%pnom_v, self%geom%Domain%mpp_domain, complete=.true.)
   
  call oops_log%trace("soca_diffusion::init() done", flush=.true.)
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_calibrate(self, f_conf)
  class(soca_diffusion),  intent(inout) :: self
  type(fckit_configuration), intent(in) :: f_conf

  real(kind=kind_real) :: stability_factor = 2.0
  real(kind=kind_real) :: stats(3) ! min, max, mean  
  character(len=1024) :: str  
  character(len=:), allocatable :: str2, str3
  integer :: i, j, idr
  type(restart_file_type) :: restart_file
  logical :: b

  real(kind=kind_real) :: fixed_scale
  real(kind=kind_real), allocatable :: hz_scales(:,:), r_tmp(:,:)

  call oops_log%trace("soca_diffusion::calibrate() starting", flush=.true.)
  call oops_log%info("ExplicitDiffusion: running calibration")

  ! Get input lengthscales. Either from:
  !  1) a fixed length scale used globally
  !  2) read in from a file
  ! the result is hz_scales containing the length scales (defined as 1 sigma of a guassian)
  allocate(hz_scales(DOMAIN_WITH_HALO))
  hz_scales = 1e10  
  if (.not. f_conf%has("scales.fixed value") .neqv. f_conf%has("scales.from file")) then
    ! that was an XOR opperation above, if you were curious
    call abor1_ftn("calibration.scales must define 1 of 'fixed value' or 'from file'")
  end if
  if ( f_conf%has("scales.fixed value")) then
    call oops_log%info("  Using fixed length scales")
    call f_conf%get_or_die("scales.fixed value", fixed_scale)
    hz_scales = fixed_scale
  else
    ! yeah, this is messy. Do it better when things are atlas-ified
    call oops_log%info("  Reading length scales from file")
    call f_conf%get_or_die("scales.from file.filename", str2)
    call f_conf%get_or_die("scales.from file.variable name", str3)
    call fms_io_init()
    idr = register_restart_field(restart_file, str2, str3, &
      hz_scales, domain=self%geom%Domain%mpp_domain)
    call restore_state(restart_file, directory='')
    call free_restart_type(restart_file)
    call fms_io_exit()
  end if
  call f_conf%get_or_die("scales.as gaussian", b)
  if (.not. b) then
    ! by default, a gaspari cohn half width is expected in the config.
    ! (but the rest of this code asssumes gaussian 1 sigma)
    ! Do the conversion if needed
    hz_scales = hz_scales / 3.57_kind_real
  end if
  call mpp_update_domains(hz_scales, self%geom%Domain%mpp_domain, complete=.true.)

  call self%calc_stats(hz_scales, stats)
  write (str, *) "  L_hz: min=", stats(1), "max=", stats(2), "mean=", stats(3)
  call oops_log%info(str)

  ! calculate the minimum number of iterations needed, rounding up to the
  ! nearest even number.
  !  M >= (L/grid_size)^2    
  allocate(r_tmp(DOMAIN_WITH_HALO))
  r_tmp = 0.0
  do j = LOOP_DOMAIN_J
    do i = LOOP_DOMAIN_I
      if (self%mask(i,j) == 0.0 ) cycle
      ! am I off by a factor of two here? eh, seems to be working
      r_tmp(i,j) = stability_factor * hz_scales(i,j)**2 * (1.0/(self%dx(i,j)**2) + 1.0/(self%dy(i,j)**2))
    end do
  end do
  call self%calc_stats(r_tmp, stats)
  self%n_iter = ceiling(stats(2))
  if (mod(self%n_iter,2) == 1) self%n_iter = self%n_iter + 1
  write (str, *) "  minimum iterations: ", self%n_iter
  call oops_log%info(str)

  ! calculate KhDt
  allocate(self%KhDt(DOMAIN_WITH_HALO))
  self%KhDt = hz_scales**2 / (2.0 * self%n_iter)

  ! calculate normalization
  call oops_log%info("Calculating normalization...")
  allocate(self%normalization(DOMAIN_WITH_HALO))
  self%normalization = 1.0
  call f_conf%get_or_die("normalization.method", str2)
  if (str2 == "brute force") then
    call self%calc_norm_bruteforce()
  else if (str2 == "randomization") then
    call f_conf%get_or_die("normalization.iterations", i)
    call self%calc_norm_randomization(i)
  else
    call abor1_ftn("ERROR: normalization.method must be 'brute force' or 'randomization'")
  end if

  call oops_log%trace("soca_diffusion::calibrate() done", flush=.true.)
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_multiply(self, dx)
  class(soca_diffusion), intent(inout) :: self
  type(soca_increment),  intent(inout) :: dx

  real(kind=kind_real), allocatable :: tmp2d(:,:)
  character(len=1024) :: str  
  integer :: f, z

  call oops_log%trace("soca_diffusion::multiply() starting", flush=.true.)

  if (self%n_iter <= 0) then
    ! uninitialized, calibrate or read should have been called before now
    call abor1_ftn("ERROR: soca_diffusion has not be initialized.")
  end if

  allocate(tmp2d(DOMAIN_WITH_HALO))

  do f=1, size(dx%fields)
    write (str, *) " multiplying field: ", dx%fields(f)%name
    do z = 1, dx%fields(f)%nz
      tmp2d = dx%fields(f)%val(:,:,z)
      call self%multiply_2D(tmp2d)
      ! NOTE: this is here because the SABERadjoint test doesn't take into account the halo correctly
      ! so we have to leave the halo alone
      dx%fields(f)%val(DOMAIN,z) = tmp2d(DOMAIN)
    end do    
    call oops_log%debug(str)
  end do

  call oops_log%trace("soca_diffusion::multiply() done", flush=.true.)
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_multiply_2D(self, field)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable :: field(:,:)

  call self%multiply_2D_ad(field)
  call self%multiply_2D_tl(field)

end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_multiply_2D_tl(self, field)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable :: field(:,:)
  
  !call mpp_update_domains(field, self%geom%Domain%mpp_domain, complete=.true.)

  ! apply grid metric
  field = field * self%inv_sqrt_area

  ! apply M/2 iterations of diffusion
  call self%diffusion_steps(field, self%n_iter / 2)

  ! apply normalization
  field = field * self%normalization
  
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_multiply_2D_ad(self, field)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable :: field(:,:)  

  ! more halo updates than we really need? but here just to be safe
  !call mpp_update_domains_ad(field, self%geom%Domain%mpp_domain, complete=.true.)

  ! apply normalization
  field = field * self%normalization

  ! TODO code the actual adjoint operator
  ! apply M/2 iterations of diffusion
  call self%diffusion_steps_ad(field, self%n_iter / 2)

  ! apply grid metric
  field = field * self%inv_sqrt_area
  
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_diffusion_steps(self, field, niter)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable, intent(inout) :: field(:,:)
  integer, intent(in) :: niter

  real(kind=kind_real), allocatable :: flux_x(:,:), flux_y(:,:), hfac(:,:)
  integer :: i, j, iter
 
  ! NOTE: flux_x(i,j) is the flux through the western edge of the grid cell.
  !  this is opposite of MOM6 conventions where u(i,j) points are east of t(i,j) points.
  !  (dosen't really matter, just something to note)
  allocate(flux_x(DOMAIN_WITH_HALO))
  allocate(flux_y(DOMAIN_WITH_HALO))
  flux_x = 0.0
  flux_y = 0.0

  ! calculate some needed constants
  allocate(hfac(DOMAIN_WITH_HALO))
  hfac = 0.0
  do j=LOOP_DOMAIN_J
    do i=LOOP_DOMAIN_I
      hfac(i,j) = (1.0/self%dy(i,j)/self%dx(i,j))
    end do
  end do

  call mpp_update_domains(field, self%geom%Domain%mpp_domain, complete=.true.)

  do iter=1,niter
      ! calculate diffusive flux on each edge of a grid box. masking out where there is land
    do j=LOOP_DOMAIN_J 
      do i=LOOP_DOMAIN_I+1 ! assume halo size is >= 1, and skip doing a halo update
        flux_x(i,j) = self%pmon_u(i,j) * 0.5 * (self%KhDt(i,j) + self%KhDt(i-1,j)) * (field(i,j) - field(i-1,j))
        flux_x(i,j) = flux_x(i,j) * self%mask(i,j) * self%mask(i-1,j)
      end do
    end do
    do j=LOOP_DOMAIN_J+1 ! assume halo size is >= 1, and skip doing a halo update
      do i=LOOP_DOMAIN_I
        flux_y(i,j) = self%pnom_v(i,j) * 0.5 * (self%KhDt(i,j) + self%KhDt(i,j-1)) * (field(i,j) - field(i,j-1))
        flux_y(i,j) = flux_y(i,j) * self%mask(i,j) * self%mask(i,j-1)
      end do
    end do

    ! time-step hz diffusion terms
    do j=LOOP_DOMAIN_J
      do i=LOOP_DOMAIN_I
        field(i,j) = field(i,j) + hfac(i,j) * &
          (flux_x(i+1, j) - flux_x(i,j) + flux_y(i, j+1) - flux_y(i,j))
      end do
    end do

    call mpp_update_domains(field, self%geom%Domain%mpp_domain, complete=.true.)
  end do
end subroutine


! ------------------------------------------------------------------------------
subroutine soca_diffusion_diffusion_steps_ad(self, field, niter)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable, intent(inout) :: field(:,:) 
  integer, intent(in) :: niter

  real(kind=kind_real), allocatable :: wrk_old(:,:), wrk_new(:,:), tmp(:,:)
  real(kind=kind_real), allocatable :: flux_x(:,:), flux_y(:,:), hfac(:,:)
  real(kind=kind_real) :: adfac
  integer :: i, j, iter
  
  ! NOTE: flux_x(i,j) is the flux through the western edge of the grid cell.
  !  this is opposite of MOM6 conventions where u(i,j) points are east of t(i,j) points.
  !  (dosen't really matter, just something to note)
  allocate(wrk_new(DOMAIN_WITH_HALO))
  allocate(wrk_old(DOMAIN_WITH_HALO))
  allocate(tmp(DOMAIN_WITH_HALO))
  allocate(flux_x(DOMAIN_WITH_HALO))
  allocate(flux_y(DOMAIN_WITH_HALO))
  wrk_old = 0.0
  wrk_new = 0.0
  flux_x = 0.0
  flux_y = 0.0

  ! calculate some needed constants
  allocate(hfac(DOMAIN_WITH_HALO))
  hfac = 0.0
  do j=LOOP_DOMAIN_J
    do i=LOOP_DOMAIN_I
      hfac(i,j) = (1.0/self%dy(i,j)/self%dx(i,j))
    end do
  end do  

  ! adjoint of convoled solution  
  ! TODO this should be called, why is it breaking?
  !call mpp_update_domains_ad(field, self%geom%Domain%mpp_domain, complete=.true.)
  do j=LOOP_DOMAIN_J
    do i=LOOP_DOMAIN_I
      wrk_old(i,j) = wrk_old(i,j) + field(i,j)
      field(i,j) = 0.0
    end do
  end do


  ! integrate adjoint hz diffusion terms
  do iter=1,niter
    tmp = wrk_new
    wrk_new = wrk_old
    wrk_old = tmp    
    !call mpp_update_domains_ad(wrk_new, self%geom%Domain%mpp_domain, complete=.true.)
    flux_x = 0.0
    flux_y = 0.0

    ! time-step adjoint hz diffusion terms
    do j=LOOP_DOMAIN_J
      do i=LOOP_DOMAIN_I
        adfac=hfac(i,j)*wrk_new(i,j)
        flux_y(i,j  ) = flux_y(i,j  ) - adfac
        flux_y(i,j+1) = flux_y(i,j+1) + adfac
        flux_x(i  ,j) = flux_x(i  ,j) - adfac
        flux_x(i+1,j) = flux_x(i+1,j) + adfac
        wrk_old(i,j) = wrk_new(i,j)
        !wrk_new(i,j) = 0.0
      end do
    end do
    wrk_new = 0.0

    ! compute adjoint diffusive flux
    do j=LOOP_DOMAIN_J+1
      do i=LOOP_DOMAIN_I
        flux_y(i,j) = flux_y(i,j) * self%mask(i,j) * self%mask(i,j-1)
        adfac = self%pnom_v(i,j) * 0.5*(self%KhDt(i,j-1)+self%KhDt(i,j)) * flux_y(i,j)
        wrk_old(i,j-1) = wrk_old(i,j-1) - adfac
        wrk_old(i,j  ) = wrk_old(i,j  ) + adfac
        flux_y(i,j) = 0.0
      end do
    end do
    do j=LOOP_DOMAIN_J
      do i=LOOP_DOMAIN_I+1
        flux_x(i,j) = flux_x(i,j) * self%mask(i,j) * self%mask(i-1,j)
        adfac = self%pmon_u(i,j) * 0.5*(self%KhDt(i-1,j)+self%KhDt(i,j)) * flux_x(i,j)
        wrk_old(i-1,j) = wrk_old(i-1,j) - adfac
        wrk_old(i,  j) = wrk_old(i,  j) + adfac
        flux_x(i,j) = 0.0
      end do
    end do
    call mpp_update_domains_ad(wrk_old, self%geom%Domain%mpp_domain, complete=.true.)
  end do

  ! set adjoint initial conditions
  field = field + wrk_old
  wrk_old = 0.0
  !call mpp_update_domains_ad(field, self%geom%Domain%mpp_domain, complete=.true.)

end subroutine

! ------------------------------------------------------------------------------
! Calculate the exact normalization weights using the brute force method
! (creating a dirac at every SINGLE point).
! You probably don't want to use this, except for testing. Use randomization.
! ------------------------------------------------------------------------------
subroutine soca_diffusion_calc_norm_bruteforce(self)
  class(soca_diffusion), intent(inout) :: self

  integer :: i, j
  logical :: local
  real(kind=kind_real), allocatable :: r_tmp(:,:), norm(:,:)

  call oops_log%info("Calculating normalization with BRUTEFORCE (eeeek!)")

  allocate(r_tmp(DOMAIN_WITH_HALO))
  allocate(norm(DOMAIN_WITH_HALO))

  self%normalization = 1.0
  do j=self%geom%jsg, self%geom%jeg
    do i=self%geom%isg, self%geom%ieg
      r_tmp = 0.0
      local = i >= self%geom%isc .and. i <= self%geom%iec .and. &
              j >= self%geom%jsc .and. j <= self%geom%jec
      if (local) r_tmp(i,j) = 1.0
      call self%multiply_2D(r_tmp)
      if (local) then
        if(self%mask(i,j) == 0.0) cycle
        norm(i,j) = 1.0 / sqrt(r_tmp(i,j))
      end if
    end do
  end do
  call mpp_update_domains(norm, self%geom%Domain%mpp_domain, complete=.true.)
  self%normalization = norm
end subroutine


! ------------------------------------------------------------------------------
! Estimate the normalization weights by creating random vectors (normally distributed)
! applying the diffusion TL, and keeping a running statistic of the variance of
! those results.
!
! Typically a good number of iterations is around 1000
! ------------------------------------------------------------------------------
subroutine soca_diffusion_calc_norm_randomization(self, iter)
  class(soca_diffusion), intent(inout) :: self
  integer, intent(in) :: iter

  real(kind=kind_real), allocatable :: field(:,:)
  real(kind=kind_real), allocatable :: s(:,:)
  real(kind=kind_real), allocatable :: m(:,:), new_m(:,:)

  integer :: n, n10pct
  character(len=1024) :: str  
  
  allocate(field(DOMAIN_WITH_HALO))
  allocate(s(DOMAIN_WITH_HALO))
  allocate(m(DOMAIN_WITH_HALO))
  allocate(new_m(DOMAIN_WITH_HALO))

  s = 0.0
  m = 0.0

  n10pct = iter/10

  do n=1,iter
    if (mod(n, n10pct) == 0) then
      write (str, *) "normalization: ", 10*n/n10pct, "% "
      call oops_log%info(str, flush=.true.)
    end if
  
    ! create a random vector
    call normal_distribution(field, 0.0_kind_real, 1.0_kind_real, n, .true.) 

    ! apply the diffusion TL
    call self%multiply_2D_tl(field)

    ! keep track of the stats needed for a running variance calculation
    ! (Welford 1962 algorithm)
    new_m = m + (field-m)/n
    s = s + (field - m)*(field - new_m)
    m = new_m
  end do
  
  ! calculate final variance
  field = (s/(iter-1)) 

  ! normalization (where ocean) is 1/sqrt(variance)
  where (self%mask == 1.0)  self%normalization = 1.0 / sqrt(field)
  
  call mpp_update_domains(self%normalization, self%geom%Domain%mpp_domain, complete=.true.)

end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_write_params(self, filename)
  class(soca_diffusion), intent(inout) :: self
  character(len=*),      intent(in)    :: filename

  type(restart_file_type) :: restart_file
  integer :: idr

  ! read from file
  call fms_io_init()  
  idr = register_restart_field(restart_file, filename, "iterations", &
    self%n_iter, domain=self%geom%Domain%mpp_domain)
  idr = register_restart_field(restart_file, filename, "khdt", &
    self%KhDt, domain=self%geom%Domain%mpp_domain)
  idr = register_restart_field(restart_file, filename, "normalization", &
    self%normalization, domain=self%geom%Domain%mpp_domain)  
  call save_restart(restart_file, directory='')
  call free_restart_type(restart_file)
  call fms_io_exit()  
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_read_params(self, filename)
  class(soca_diffusion), intent(inout) :: self
  character(len=*),      intent(in)    :: filename
  
  type(restart_file_type) :: restart_file
  integer :: idr

  allocate(self%KhDt(DOMAIN_WITH_HALO))
  allocate(self%normalization(DOMAIN_WITH_HALO))
 
  ! initialize with safe values, because not all halo points 
  ! are updated in a MOM6 halo update
  self%KhDt = 0.0
  self%normalization = 1.0

  ! read from file
  call fms_io_init()
  idr = register_restart_field(restart_file, filename, "iterations", &
    self%n_iter, domain=self%geom%Domain%mpp_domain)
  idr = register_restart_field(restart_file, filename, "khdt", &
    self%KhDt, domain=self%geom%Domain%mpp_domain)
  idr = register_restart_field(restart_file, filename, "normalization", &
    self%normalization, domain=self%geom%Domain%mpp_domain)  
  call restore_state(restart_file, directory='')
  call free_restart_type(restart_file)
  call fms_io_exit()

  ! update halos
  call mpp_update_domains(self%normalization, self%geom%Domain%mpp_domain, complete=.true.)
  call mpp_update_domains(self%KhDt, self%geom%Domain%mpp_domain, complete=.true.)
end subroutine

! ------------------------------------------------------------------------------

end module soca_diffusion_mod
