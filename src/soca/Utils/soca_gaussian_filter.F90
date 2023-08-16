! (C) Copyright 2023-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_gaussian_filter_mod

use kinds, only: kind_real
use fckit_mpi_module, only: fckit_mpi_min, fckit_mpi_max, fckit_mpi_sum
use logger_mod
use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad

use soca_geom_mod, only : soca_geom

implicit none
private

! ------------------------------------------------------------------------------

type, public :: soca_gaussian_filter
private
  real(kind_real), allocatable :: weights_x (:,:,:)
  real(kind_real), allocatable :: weights_y (:,:,:)
  class(soca_geom), pointer :: geom
  integer :: n_iter

contains
  procedure :: init => soca_gaussian_filter_init
  procedure :: mult => soca_gaussian_filter_mult
  ! procedure :: mult_ad => soca_gaussian_filter_mult_ad

end type soca_gaussian_filter

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------
! calculate the masked global min/max/mean values for a given field
! TODO: is there any chance we need to also operate on unmasked fields?!?
subroutine calc_stats(field, geom, mask, stats)
  real(kind=kind_real),  intent(in) :: field(:,:)
  type(soca_geom),       intent(in) :: geom
  logical,               intent(in) :: mask(:,:)
  real(kind=kind_real), intent(out) :: stats(3)
  
  real(kind=kind_real) :: l_min, l_max, l_sum, l_count, g_count

  l_min =  minval(    field(geom%isc:geom%iec, geom%jsc:geom%jec), &
                  mask=mask(geom%isc:geom%iec, geom%jsc:geom%jec))
  l_max =  maxval(    field(geom%isc:geom%iec, geom%jsc:geom%jec), &
                  mask=mask(geom%isc:geom%iec, geom%jsc:geom%jec))
  l_sum =     sum(    field(geom%isc:geom%iec, geom%jsc:geom%jec), &
                  mask=mask(geom%isc:geom%iec, geom%jsc:geom%jec))
  l_count = count(     mask(geom%isc:geom%iec, geom%jsc:geom%jec))

  call geom%f_comm%allreduce(l_min, stats(1), fckit_mpi_min())
  call geom%f_comm%allreduce(l_max, stats(2), fckit_mpi_max())
  call geom%f_comm%allreduce(l_sum, stats(3), fckit_mpi_sum())
  call geom%f_comm%allreduce(l_count, g_count, fckit_mpi_sum())
  stats(3) = stats(3) / g_count 
end subroutine

! ------------------------------------------------------------------------------
! Initialize the gaussian filter with the given horizontal scales
! 1. determine the minimum number of iterations required
! 2. calculate the weights based on the number of iteraetions, and the hz scales
subroutine soca_gaussian_filter_init(self, geom, hz_scales)
  class(soca_gaussian_filter), intent(inout) :: self
  class(soca_geom),       target, intent(in) :: geom
  real(kind=kind_real), allocatable, intent(in) :: hz_scales(:,:)

  logical :: mask(geom%isd:geom%ied,geom%jsd:geom%jed)
  real(kind=kind_real) :: dx(geom%isd:geom%ied,geom%jsd:geom%jed)
  real(kind=kind_real) :: stats(3) ! min,max,mean  
  real(kind=kind_real) :: r_tmp(geom%isd:geom%ied,geom%jsd:geom%jed)
  integer :: n_min, i, j
  character(len=1024) :: str

  call oops_log%debug("soca_gaussian_filter::init ", flush=.true.)
  self%geom => geom

  ! get the grid spacing, and desired filter standard deviation
  ! TODO, get the correct dx/dy values!    
  call calc_stats(hz_scales, geom, mask, stats )
  write (str, *) "  Hz Scales:  min=", stats(1), "max=", stats(2), "mean=", stats(3)
  call oops_log%debug(str)

  dx = sqrt(geom%cell_area)
  mask = geom%mask2d == 1.0
  call calc_stats(dx, geom, mask, stats )
  write (str, *) "  dx lengths: min=", stats(1), "max=", stats(2), "mean=", stats(3)  
  call oops_log%debug(str)

  ! calculate the minimum number of iterations needed
  ! n_min = ceil( 3/2 * sigma^2 / dx^2 )
  r_tmp = 0.0
  where (dx > 0) r_tmp = 3.0/2.0 * hz_scales**2 / dx**2
  call calc_stats(r_tmp, geom, mask, stats )
  n_min = ceiling(stats(2)) ! use global max
  write (str, *) "  minimum iterations: ", n_min  
  call oops_log%debug(str)
  self%n_iter = n_min + 2 ! add 1, just for the heck of it

  
  ! calculate the weights for the kernel
  ! (For now, assume hz isotropic)
  ! w = sigma^2 / (2 * n * dx^2)
  allocate(self%weights_x(-1:1,geom%isd:geom%ied,geom%jsc:geom%jed))
  self%weights_x = 0.0
  do i = self%geom%isd, self%geom%ied
    do j = self%geom%jsd, self%geom%jed
      if (geom%mask2d(i,j) == 0.0) cycle

      if (dx(i,j) == 0.0) cycle ! why am i geting dx of 0.0 in some places??
      self%weights_x(-1,i,j) = hz_scales(i,j)**2 / (2.0 * self%n_iter * dx(i,j)**2)          
      self%weights_x(1,i,j) = self%weights_x(-1,i,j)
      self%weights_x(0,i,j) = 1.0 - 2.0 * self%weights_x(-1,i,j)
      
      ! redistribute weights at boundaries
      if(geom%mask2d(i-1, j) == 0.0) then
        self%weights_x(0,i,j) = self%weights_x(0,i,j) + self%weights_x(-1,i,j)
        self%weights_x(-1,i,j) = 0.0
      end if
      if(geom%mask2d(i+1, j) == 0.0) then
        self%weights_x(0,i,j) = self%weights_x(0,i,j) + self%weights_x(1,i,j)
        self%weights_x(1,i,j) = 0.0
      end if      
    end do
  end do  
  
  allocate(self%weights_y(-1:1,geom%isd:geom%ied,geom%jsc:geom%jed))
  self%weights_y = 0.0
  do i = self%geom%isd, self%geom%ied
    do j = self%geom%jsd, self%geom%jed
      if (geom%mask2d(i,j) == 0.0) cycle

      if (dx(i,j) == 0.0) cycle ! why am i geting dx of 0.0 in some places??
      self%weights_y(-1,i,j) = hz_scales(i,j)**2 / (2.0 * self%n_iter * dx(i,j)**2)    
      self%weights_y(1,i,j) = self%weights_y(-1,i,j)
      self%weights_y(0,i,j) = 1.0 - 2.0 * self%weights_y(-1,i,j)

      ! redistribute weights at boundaries
      if(geom%mask2d(i, j-1) == 0.0) then
        self%weights_y(0,i,j) = self%weights_y(0,i,j) + self%weights_y(-1,i,j)
        self%weights_y(-1,i,j) = 0.0
      end if
      if(geom%mask2d(i, j+1) == 0.0) then
        self%weights_y(0,i,j) = self%weights_y(0,i,j) + self%weights_y(1,i,j)
        self%weights_y(1,i,j) = 0.0
      end if   

    end do
  end do  


end subroutine

! ------------------------------------------------------------------------------

subroutine soca_gaussian_filter_mult(self, xin, xout )
  class(soca_gaussian_filter),       intent(inout) :: self
  real(kind=kind_real), allocatable, intent(in)    :: xin (:,:,:)
  real(kind=kind_real), allocatable, intent(inout) :: xout (:,:,:)

  character(len=1024) :: str
  integer :: iter, i, j
  real(kind=kind_real), allocatable :: xtmp(:,:,:)

  write (str, *) "soca_gaussian_filter:mult Running for ", self%n_iter, " iterations"
  call oops_log%debug(str, flush=.true.)
  allocate(xtmp(self%geom%isd:self%geom%ied, self%geom%jsd:self%geom%jed, size(xin,3)))

  xout = xin

  do iter=0, self%n_iter    
    
    ! x direction
    xtmp = xout
    do i = self%geom%isc, self%geom%iec
      do j = self%geom%jsc, self%geom%jec
        if (self%geom%mask2d(i,j) == 0.0) cycle
        xout(i,j,:) = &
            self%weights_x( 0, i, j) * xtmp(i,   j, :) &
          + self%weights_x(-1, i, j) * xtmp(i-1, j, :) &
          + self%weights_x(+1, i, j) * xtmp(i+1, j, :)        
      end do
    end do
    call mpp_update_domains(xout, self%geom%Domain%mpp_domain, complete=.true.)

    ! y direction
    xtmp = xout
    do i = self%geom%isc, self%geom%iec
      do j = self%geom%jsc, self%geom%jec
        if (self%geom%mask2d(i,j) == 0.0) cycle              
        xout(i,j,:) = &
            self%weights_y( 0, i, j) * xtmp(i,   j, :) &
          + self%weights_y(-1, i, j) * xtmp(i, j-1, :) &
          + self%weights_y(+1, i, j) * xtmp(i, j+1, :)
      end do
    end do
    call mpp_update_domains(xout, self%geom%Domain%mpp_domain, complete=.true.)
  end do
end subroutine

end module soca_gaussian_filter_mod
