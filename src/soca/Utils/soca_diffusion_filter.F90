! (C) Copyright 2023-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_diffusion_filter_mod

use atlas_module, only: atlas_geometry
use kinds, only: kind_real
use fckit_mpi_module, only: fckit_mpi_min, fckit_mpi_max, fckit_mpi_sum
use logger_mod
use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad

use soca_geom_mod, only : soca_geom

implicit none
private

! ------------------------------------------------------------------------------

type, public :: soca_diffusion_filter
private
  real(kind_real), allocatable :: weights_x (:,:,:)
  real(kind_real), allocatable :: weights_y (:,:,:)  
  real(kind_real), allocatable :: normalization (:,:)
  class(soca_geom), pointer :: geom
  integer :: n_iter
  logical :: normalize

contains
  procedure :: init => soca_diffusion_filter_init
  procedure :: mult => soca_diffusion_filter_mult
  
  procedure, private :: calc_norm => soca_diffusion_filter_calcnorm_bruteforce
  ! procedure :: mult_ad => soca_diffusion_filter_mult_ad

end type soca_diffusion_filter

! ------------------------------------------------------------------------------

contains


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
! calculate the grid spacing (dx,dy) based on the given lat lon values
subroutine calc_dx_dy(geom, dx, dy)
  class(soca_geom), intent(in) :: geom
  real(kind=kind_real), allocatable, intent(inout) :: dx(:,:), dy(:,:)

  integer :: i,j
  type(atlas_geometry) :: atlas_geom
  real(kind=kind_real) :: lat_m1, lat_p1, lon_m1, lon_p1

  dx = 1e10
  dy = 1e10
  atlas_geom = atlas_geometry("Earth")
  do j = geom%jsc, geom%jec
    do i = geom%isc, geom%iec
      if (geom%mask2d(i,j) == 0.0) then
        ! todo remove this?
        dx(i,j) = 1e10
        dy(i,j) = 1e10
        cycle
      end if
      
      lat_m1 = geom%lat(i,j-1)
      lat_p1 = geom%lat(i,j+1)
      lon_m1 = geom%lon(i,j-1)
      lon_p1 = geom%lon(i,j+1)    
      if (lon_m1 < -900) lon_m1 = geom%lon(i,j)
      if (lon_p1 < -900) lon_m1 = geom%lon(i,j)
      if (lat_m1 < -900) then
        lat_m1 = 2.0*geom%lat(i,j)-geom%lat(i+1,j)
      end if
      if (lat_p1 < -900) then
         lat_p1 = 2.0*geom%lat(i,j)-geom%lat(i-1,j)
         if (lat_p1 > 90) then
          lat_p1 = 2.0* 90 - lat_p1
          lon_p1 = lon_p1 + 180.0
         end if
      end if
      
      dx(i,j) = 0.5 * atlas_geom%distance(geom%lon(i-1, j), geom%lat(i-1,j), geom%lon(i+1,j), geom%lat(i+1,j))
      dy(i,j) = 0.5 * atlas_geom%distance(lon_m1, lat_m1, lon_p1, lat_p1)
    end do
  end do

  call mpp_update_domains(dx, geom%Domain%mpp_domain, complete=.true.)
  call mpp_update_domains(dy, geom%Domain%mpp_domain, complete=.true.)
end subroutine


! ------------------------------------------------------------------------------
! Initialize the diffusion filter with the given horizontal scales
! 1. determine the minimum number of iterations required
! 2. calculate the weights based on the number of iteraetions, and the hz scales
subroutine soca_diffusion_filter_init(self, geom, hz_scales)
  class(soca_diffusion_filter),       intent(inout) :: self
  class(soca_geom),          target, intent(in)    :: geom
  real(kind=kind_real), allocatable, intent(inout)    :: hz_scales(:,:)

  logical, allocatable :: mask(:,:)
  real(kind=kind_real), allocatable :: dx(:,:)
  real(kind=kind_real), allocatable :: dy(:,:)
  real(kind=kind_real) :: stats(3) ! min,max,mean
  real(kind=kind_real), allocatable :: r_tmp(:,:)
  integer :: n_min, i, j
  character(len=1024) :: str

  call oops_log%debug("soca_diffusion_filter::init ", flush=.true.)
  self%geom => geom
  mask = geom%mask2d == 1.0
  self%normalize = .true.

  ! calculate the grid spacing
  ! (note, I *could* just get dx/dy from mom6, but I'm calculating
  !  it here anyway in the hopes of being more generic someday)
  allocate(dx(geom%isd:geom%ied,geom%jsd:geom%jed))
  allocate(dy(geom%isd:geom%ied,geom%jsd:geom%jed))
  call calc_dx_dy(geom, dx, dy)
  call calc_stats(dx, geom, stats )
  write (str, *) "  dx lengths: min=", stats(1), "max=", stats(2), "mean=", stats(3)
  call oops_log%debug(str)
  call calc_stats(dy, geom, stats )
  write (str, *) "  dy lengths: min=", stats(1), "max=", stats(2), "mean=", stats(3)
  call oops_log%debug(str)

  ! check the horizontal scales
  call calc_stats(hz_scales, geom, mask, stats )
  write (str, *) "  Hz Scales:  min=", stats(1), "max=", stats(2), "mean=", stats(3)
  call oops_log%debug(str)

  ! calculate the minimum number of iterations needed
  ! M >= (L/grid_size)^2
  allocate(r_tmp(geom%isd:geom%ied,geom%jsd:geom%jed))
  r_tmp = 0.0
  do j=geom%jsd,geom%jed
    do i=geom%isd,geom%ied
      if (.not. mask(i,j)) cycle
      r_tmp(i,j) = (hz_scales(i,j) / min(dx(i,j),dy(i,j)))**2
    end do
  end do
  call calc_stats(r_tmp, geom, stats )
  n_min = ceiling(stats(2)) ! use global max
  self%n_iter = n_min + 2  ! add 1 or 2, just to be safe
  write (str, *) "  minimum iterations: ", self%n_iter
  call oops_log%debug(str)
  self%n_iter = n_min +5  ! add 1 or 2, just to be safe
  
  ! calculate the weights
  ! w = L^2 / 2M
  allocate(self%weights_x(-1:1,geom%isd:geom%ied,geom%jsc:geom%jed))
  allocate(self%weights_y(-1:1,geom%isd:geom%ied,geom%jsc:geom%jed))
  self%weights_x = 0.0
  self%weights_y = 0.0
  do j = self%geom%jsc, self%geom%jec
  do i = self%geom%isc, self%geom%iec
      if (geom%mask2d(i,j) == 0.0) cycle

      ! BC = zero derivative
      self%weights_x(-1,i,j) = merge(1.0, 0.0, geom%mask2d(i-1,j)==1.0)
      self%weights_x( 1,i,j) = merge(1.0, 0.0, geom%mask2d(i+1,j)==1.0)
      self%weights_x( 0,i,j) = -sum(self%weights_x(:,i,j))
      self%weights_y(-1,i,j) = merge(1.0, 0.0, geom%mask2d(i,j-1)==1.0)
      self%weights_y( 1,i,j) = merge(1.0, 0.0, geom%mask2d(i,j+1)==1.0)
      self%weights_y( 0,i,j) = -sum(self%weights_y(:,i,j))
      
      ! ! BC = zero value
      ! self%weights_x( 0,i,j) = -2.0
      ! self%weights_y( 0,i,j) = -2.0

      self%weights_x( :,i,j) = self%weights_x(:,i,j) * (hz_scales(i,j)/dx(i,j))**2 / (2.0*self%n_iter)
      self%weights_y( :,i,j) = self%weights_y(:,i,j) * (hz_scales(i,j)/dy(i,j))**2 / (2.0*self%n_iter) 
      
    end do
  end do

  if (self%normalize) then
    allocate(self%normalization(geom%isd:geom%ied, geom%jsd:geom%jed))    
    call self%calc_norm
    call calc_stats(self%normalization, geom, stats )
    write (str, *) "  normalization: min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%debug(str)
  end if

end subroutine

! ------------------------------------------------------------------------------
subroutine soca_diffusion_filter_mult(self, xin, xout )
  class(soca_diffusion_filter),       intent(inout) :: self
  real(kind=kind_real), allocatable, intent(in)    :: xin (:,:,:)
  real(kind=kind_real), allocatable, intent(inout) :: xout (:,:,:)

  character(len=1024) :: str
  integer :: iter, i, j, k
  real(kind=kind_real), allocatable :: xtmp(:,:,:)
  real(kind=kind_real) :: x0, xm1, xp1

  write (str, *) "soca_diffusion_filter:mult Running for ", self%n_iter, " iterations"
  call oops_log%debug(str, flush=.true.)
  allocate(xtmp(self%geom%isd:self%geom%ied, self%geom%jsd:self%geom%jed, size(xin,3)))

  xout = xin

  if (self%normalize) then
    do k=1, size(xout, dim=3) 
      xout(:,:,k) = xout(:,:,k) * self%normalization          
    end do
  end if

  do iter=1, self%n_iter

    ! x direction
    xtmp = xout
    do i = self%geom%isc, self%geom%iec
      do j = self%geom%jsc, self%geom%jec
        if (self%geom%mask2d(i,j) == 0.0) cycle        
        xout(i,j,:) = xtmp(i,j,:) &
           + self%weights_x(-1,i,j) * xtmp(i-1,j,:) &
           + self%weights_x( 0,i,j) * xtmp(i,j,:) & 
           + self%weights_x( 1,i,j) * xtmp(i+1,j,:)
      end do
    end do
    ! TODO, in the above code, go into the halo in the y direction by 1 
    ! so we can skip this first halo exchange    
    call mpp_update_domains(xout, self%geom%Domain%mpp_domain, complete=.true.)

    ! y direction
    xtmp = xout
    do i = self%geom%isc, self%geom%iec
      do j = self%geom%jsc, self%geom%jec
        if (self%geom%mask2d(i,j) == 0.0) cycle        
        xout(i,j,:) = xtmp(i,j,:) &
           + self%weights_y(-1,i,j) * xtmp(i,j-1,:) &
           + self%weights_y( 0,i,j) * xtmp(i,j,:) &
           + self%weights_y( 1,i,j) * xtmp(i,j+1,:)
      end do
    end do
    call mpp_update_domains(xout, self%geom%Domain%mpp_domain, complete=.true.)
  end do

  if (self%normalize) then
    do k=1, size(xout, dim=3)
      xout(:,:,k) = xout(:,:,k) * self%normalization
      !xout(:,:,k) =  self%normalization
    end do
  end if
end subroutine

! ------------------------------------------------------------------------------
subroutine soca_diffusion_filter_calcnorm_bruteforce(self)
  class(soca_diffusion_filter),       intent(inout) :: self

  integer :: i, j
  logical :: local
  real(kind=kind_real), allocatable :: r_tmp(:,:,:), r_tmp2(:,:,:)

  allocate(r_tmp(self%geom%isd:self%geom%ied, self%geom%jsd:self%geom%jed, 1))
  allocate(r_tmp2(self%geom%isd:self%geom%ied, self%geom%jsd:self%geom%jed, 1))

  self%normalize = .false.
  self%normalization = 1.0
  do j=self%geom%jsg, self%geom%jeg
    do i=self%geom%isg, self%geom%ieg
      r_tmp = 0.0
      local = i >= self%geom%isc .and. i <= self%geom%iec .and. &
              j >= self%geom%jsc .and. j <= self%geom%jec
      if (local) r_tmp(i,j,:) = 1.0
      call self%mult(r_tmp, r_tmp2)
      if (local) then
        if(self%geom%mask2d(i,j) == 0.0) cycle
        self%normalization(i,j) = 1.0 / sqrt(r_tmp2(i,j,1))
      end if
    end do
  end do

  self%normalize = .true.
end subroutine

end module soca_diffusion_filter_mod
