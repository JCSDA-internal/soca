! (C) Copyright 2023-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_diffusion_mod

! use atlas_module, only: atlas_geometry
use kinds, only: kind_real
use fckit_mpi_module, only: fckit_mpi_min, fckit_mpi_max, fckit_mpi_sum
use logger_mod
use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad

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

  !real(kind_real), allocatable :: w (:,:,:)
  
!   real(kind_real), allocatable :: weights_x (:,:,:)
!   real(kind_real), allocatable :: weights_y (:,:,:)  
!   real(kind_real), allocatable :: normalization (:,:)
  logical, allocatable :: mask(:,:)  
  class(soca_geom), pointer :: geom
  
  integer :: n_iter
!   logical :: normalize

contains
  procedure :: init => soca_diffusion_init
  procedure :: calibrate => soca_diffusion_calibrate

  procedure, private :: calc_stats => soca_diffusion_calc_stats
!   procedure :: mult => soca_diffusion_filter_mult
  
!   procedure, private :: calc_norm => soca_diffusion_filter_calcnorm_bruteforce
!   !procedure, private :: calc_norm => soca_diffusion_filter_calcnorm_randomization

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
                  mask=self%geom%mask2d(DOMAIN)==1.0)
  l_max =  maxval(field(DOMAIN), &
                  mask=self%geom%mask2d(DOMAIN)==1.0)
  l_sum =     sum(field(DOMAIN), &
                  mask=self%geom%mask2d(DOMAIN)==1.0)
  l_count = count(self%geom%mask2d(DOMAIN)==1.0)

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
  allocate(self%dx(DOMAIN_WITH_HALO))
  allocate(self%dy(DOMAIN_WITH_HALO))
  self%dx = geom%dx
  self%dy = geom%dy
  call self%calc_stats(geom%dx,stats)
    write (str, *) "  dx:   min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)
  call self%calc_stats(geom%dy, stats)
    write (str, *) "  dy:   min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)
  call self%calc_stats(geom%cell_area, stats)
    write (str, *) "  area: min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)

  ! calculate derive parameters
  allocate(self%mask(DOMAIN_WITH_HALO))
  self%mask = self%geom%mask2d == 1.0

  allocate(self%inv_sqrt_area(DOMAIN_WITH_HALO))
  do j = LOOP_DOMAIN_J
    do i = LOOP_DOMAIN_I
      self%inv_sqrt_area(i,j) = 1.0 / sqrt(self%geom%dx(i,j)*self%geom%dy(i,j))
    end do
  end do

  ! update halos
  call mpp_update_domains(self%inv_sqrt_area, self%geom%Domain%mpp_domain, complete=.true.)

  ! call calc_stats(self%inv_sqrt_area, self%geom, stats)
  ! write (str, *) "  debug: min=", stats(1), "max=", stats(2), "mean=", stats(3)
  ! call oops_log%info(str)

    
  call oops_log%trace("soca_diffusion::init() done", flush=.true.)
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_calibrate(self)
  class(soca_diffusion), intent(inout) :: self

  real(kind=kind_real) :: stats(3) ! min, max, mean
  character(len=1024) :: str  
  integer :: i, j

  real(kind=kind_real), allocatable :: hz_scales(:,:), r_tmp(:,:)

  call oops_log%trace("soca_diffusion::calibrate() starting", flush=.true.)
 
  ! TODO get input lengthscales
  allocate(hz_scales(DOMAIN_WITH_HALO))
  hz_scales = 500.0e3
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
      if (.not. self%mask(i,j) ) cycle
      r_tmp(i,j) = (hz_scales(i,j) / min(self%dx(i,j), self%dy(i,j))) ** 2
    end do
  end do
  call self%calc_stats(r_tmp, stats)
  self%n_iter = ceiling(stats(2))
  if (mod(self%n_iter,2) == 1) self%n_iter = self%n_iter + 1
  write (str, *) "  minimum iterations: ", self%n_iter
  call oops_log%info(str)

  ! TODO calculate Kh

  ! TODO calculate normalization

  call oops_log%trace("soca_diffusion::calibrate() done", flush=.true.)
end subroutine

! ------------------------------------------------------------------------------

end module soca_diffusion_mod
