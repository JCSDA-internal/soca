! (C) Copyright 2023-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_gaussian_filter_mod

use kinds, only: kind_real
use fckit_mpi_module, only: fckit_mpi_min, fckit_mpi_max, fckit_mpi_sum
use logger_mod

use soca_geom_mod, only : soca_geom

implicit none
private

! ------------------------------------------------------------------------------

type, public :: soca_gaussian_filter
private
  real(kind_real), allocatable :: weights (:,:,:)

contains
  procedure :: init => soca_gaussian_filter_init

end type soca_gaussian_filter

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

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

subroutine soca_gaussian_filter_init(self, geom, hz_scales)
  class(soca_gaussian_filter), intent(inout) :: self
  class(soca_geom),            intent(in)    :: geom
  real(kind=kind_real),        intent(in)    :: hz_scales(:,:)

  logical :: mask(geom%isd:geom%ied,geom%jsd:geom%jed)
  real(kind=kind_real) :: dx(geom%isd:geom%ied,geom%jsd:geom%jed)
  real(kind=kind_real) :: stats(3) ! min,max,mean  
  real(kind=kind_real) :: r_tmp(geom%isd:geom%ied,geom%jsd:geom%jed)
  integer :: n_min
  character(len=1024) :: str

  call oops_log%debug("soca_gaussian_filter::init ", flush=.true.)

  ! get the grid spacing, and desired filter standard deviation
  ! TODO, get the correct dx/dy values!    
  call calc_stats(hz_scales, geom, mask, stats )
  write (str, *), "  Hz Scales:  min=", stats(1), "max=", stats(2), "mean=", stats(3)
  call oops_log%debug(str)

  dx = sqrt(geom%cell_area)
  mask = geom%mask2d == 1.0
  call calc_stats(dx, geom, mask, stats )
  write (str, *), "  dx lengths: min=", stats(1), "max=", stats(2), "mean=", stats(3)  
  call oops_log%debug(str)

  ! calculate the minimum number of iterations needed
  ! n_min = ceil( 3/2 * sigma^2 / dx^2 )
  r_tmp = 0.0
  where (dx > 0) r_tmp = 3.0/2.0 * hz_scales**2 / dx**2
  call calc_stats(r_tmp, geom, mask, stats )
  n_min = ceiling(stats(2)) ! use global max
  write (str, *), "  minimum iterations: ", n_min
  call oops_log%debug(str)
  
  ! calculate the weights for the kernel
  ! w = sigma^2 / (2 * n * dx^2)
  allocate(self%weights(3,geom%isd:geom%ied,geom%jsd:geom%jed))
  self%weights = 0.0
  where (dx > 0) self%weights(1,:,:) = hz_scales**2 / (2.0 * n_min * dx**2)
  self%weights(3,:,:) = self%weights(1,:,:)
  self%weights(2,:,:) = 1.0 - 2.0 * self%weights(1,:,:)
  call calc_stats(self%weights(2,:,:), geom, mask, stats )
  write (str, *), "  cntr weights: min=", stats(1), "max=", stats(2), "mean=", stats(3)  
  call oops_log%debug(str)

end subroutine

! ------------------------------------------------------------------------------

end module soca_gaussian_filter_mod
