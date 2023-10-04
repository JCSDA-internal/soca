! (C) Copyright 2023-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_diffusion_mod

! use atlas_module, only: atlas_geometry
use kinds, only: kind_real
use fckit_mpi_module, only: fckit_mpi_min, fckit_mpi_max, fckit_mpi_sum
use logger_mod
! use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad

use soca_geom_mod, only : soca_geom

implicit none
private

! ------------------------------------------------------------------------------

type, public :: soca_diffusion
 private
  real(kind_real), allocatable :: w (:,:,:)
!   real(kind_real), allocatable :: weights_x (:,:,:)
!   real(kind_real), allocatable :: weights_y (:,:,:)  
!   real(kind_real), allocatable :: normalization (:,:)
  class(soca_geom), pointer :: geom
!   integer :: n_iter
!   logical :: normalize

contains
  procedure :: init => soca_diffusion_init
  procedure :: calibrate => soca_diffusion_calibrate
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
subroutine calc_stats(field, geom, stats)
  real(kind=kind_real), allocatable, intent(in) :: field(:,:)
  type(soca_geom),                   intent(in) :: geom
  real(kind=kind_real),             intent(out) :: stats(3)

  real(kind=kind_real) :: l_min, l_max, l_sum, l_count, g_count

  l_min =  minval(field(geom%isc:geom%iec, geom%jsc:geom%jec), &
                  mask=geom%mask2d(geom%isc:geom%iec, geom%jsc:geom%jec)==1.0)
  l_max =  maxval(field(geom%isc:geom%iec, geom%jsc:geom%jec), &
                  mask=geom%mask2d(geom%isc:geom%iec, geom%jsc:geom%jec)==1.0)
  l_sum =     sum(field(geom%isc:geom%iec, geom%jsc:geom%jec), &
                  mask=geom%mask2d(geom%isc:geom%iec, geom%jsc:geom%jec)==1.0)
  l_count = count(geom%mask2d(geom%isc:geom%iec, geom%jsc:geom%jec)==1.0)

  call geom%f_comm%allreduce(l_min, stats(1), fckit_mpi_min())
  call geom%f_comm%allreduce(l_max, stats(2), fckit_mpi_max())
  call geom%f_comm%allreduce(l_sum, stats(3), fckit_mpi_sum())
  call geom%f_comm%allreduce(l_count, g_count, fckit_mpi_sum())
  stats(3) = stats(3) / g_count
end subroutine

! ------------------------------------------------------------------------------
subroutine soca_diffusion_init(self, geom)
  class(soca_diffusion), intent(inout) :: self
  class(soca_geom), target, intent(in) :: geom

  real(kind=kind_real) :: stats(3) ! min, max, mean
  character(len=1024) :: str  

  call oops_log%trace("soca_diffusion::init() starting", flush=.true.)
  self%geom => geom

  ! TODO: move these to a separate "calibrate" subroutine
  ! stats
  call oops_log%info("ExplicitDiffusion: Initializing grid...")
  call calc_stats(geom%dx, self%geom, stats)
    write (str, *) "  dx:   min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)
  call calc_stats(geom%dy, self%geom, stats)
    write (str, *) "  dy:   min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)
  call calc_stats(geom%cell_area, self%geom, stats)
    write (str, *) "  area: min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)

  call oops_log%trace("soca_diffusion::init() done", flush=.true.)
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_diffusion_calibrate(self)
  class(soca_diffusion), intent(inout) :: self

  call oops_log%trace("soca_diffusion::calibrate() starting", flush=.true.)

  ! TODO calculate metrics
  ! TODO check the horizontal scales
  ! TODO calculate the minimum number of iterations needed
  ! TODO calculate Kh

  ! TODO calculate normalization

  call oops_log%trace("soca_diffusion::calibrate() done", flush=.true.)
end subroutine

! ------------------------------------------------------------------------------

end module soca_diffusion_mod
