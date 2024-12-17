! (C) Copyright 2022-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_soca2cice_mod

use atlas_module, only: atlas_geometry, atlas_indexkdtree
use fckit_configuration_module, only: fckit_configuration
use fckit_exception_module, only: fckit_exception
use fckit_mpi_module, only: fckit_mpi_comm
use kinds, only: kind_real

use icepack_itd
use icepack_mushy_physics, only: liquidus_temperature_mush
use icepack_therm_shared, only: l_brine, icepack_ice_temperature
use icepack_parameters, only: icepack_init_parameters, icepack_recompute_constants
use icepack_parameters, only: ktherm, heat_capacity
use icepack_parameters, only: rhos, Lfresh, cp_ice

use soca_geom_mod, only: soca_geom
use soca_state_mod, only: soca_state
use soca_fields_mod, only: soca_field
use soca_ciceutils_mod, only: cice_state

implicit none
private

integer :: root=0

!> analysis to cice
!!
!! - forward: deaggregates a 2D analysis of sea-ice and inserts
!!            analysis in CICE restarts
!! - inverse: TODO(G), aggregates seaice variables along CICE sea-ice
!!            categories, save the aggregated variables in a file
!!            readable by soca

type, public :: soca_soca2cice_params
   real(kind=kind_real) :: seaice_edge
   logical :: shuffle
   logical :: rescale_prior
   real(kind=kind_real) :: rescale_min_hice
   real(kind=kind_real) :: rescale_min_hsno
end type soca_soca2cice_params

type, public :: soca_soca2cice
   type(fckit_mpi_comm) :: f_comm
   integer :: myrank
   integer :: ncat, ni, nj, ice_lev, sno_lev
   character(len=:), allocatable :: rst_filename
   character(len=:), allocatable :: rst_out_filename
   type(cice_state) :: cice
   type(atlas_indexkdtree) :: kdtree
   type(soca_soca2cice_params) :: arctic, antarctic
contains
  procedure :: setup => soca_soca2cice_setup
  procedure :: changevar => soca_soca2cice_changevar
  procedure, private :: shuffle_ice
  procedure, private :: check_ice_bounds
  procedure, private :: prior_dist_rescale
  procedure, private :: cleanup_ice
end type soca_soca2cice


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Initialization of the nonlinear state to cice change of variable.
!!
subroutine soca_soca2cice_setup(self, geom)
  class(soca_soca2cice), intent(inout) :: self
  type(soca_geom), target, intent(in)  :: geom !< geometry

  integer(kind=4) :: ncid
  integer(kind=4) :: dimid
  integer(kind=4) :: varid

  integer :: myrank, root=0, count
  real(kind=kind_real), allocatable :: buffer(:)

  real(kind=kind_real) :: aice, aice0
  integer, allocatable :: ij(:,:)
  integer :: l, n_src, i, j, n
  type(atlas_geometry) :: ageometry

  ! Communicator
  self%f_comm = geom%f_comm
  self%myrank = geom%f_comm%rank()

  ! Initialize icepack's global variables ...
  call icepack_init_parameters()
  call icepack_recompute_constants()
  l_brine = .true.
  ktherm = 2
  heat_capacity = .true.

  ! initialize cice
  call self%cice%init(geom, self%rst_filename, self%rst_out_filename, self%ice_lev, self%sno_lev)

  ! read cice fields from restart
  call self%cice%read(geom)

  ! Initialize kd-tree
  ageometry = atlas_geometry("UnitSphere")
  self%kdtree = atlas_indexkdtree(ageometry)
  call self%kdtree%reserve(self%cice%agg%n_src)
  call self%kdtree%build(self%cice%agg%n_src, self%cice%agg%lon, self%cice%agg%lat)

end subroutine soca_soca2cice_setup

! ------------------------------------------------------------------------------
!> soca state to model
!!
subroutine soca_soca2cice_changevar(self, geom, xa, xm)
  class(soca_soca2cice), intent(inout) :: self
  type(soca_geom), target, intent(in)  :: geom
  type(soca_state),         intent(in) :: xa
  type(soca_state),      intent(inout) :: xm

  ! fix bounds
  call self%check_ice_bounds(geom, xm)

  ! add ice in the background where needed
  if (self%arctic%shuffle .or. self%antarctic%shuffle) call self%shuffle_ice(geom, xm)

  ! de-aggregate using the prior distribution
  if (self%arctic%rescale_prior .or. self%antarctic%rescale_prior) call self%prior_dist_rescale(geom, xm)

  ! cleanup seaice state
  call self%cleanup_ice(geom, xm)

  ! write cice restart
  call self%cice%write(geom)

end subroutine soca_soca2cice_changevar

! ------------------------------------------------------------------------------
!> model to soca state
!!
subroutine soca_soca2cice_changevarinv(self, xa, xm)
  class(soca_soca2cice), intent(inout) :: self
  type(soca_state),    intent(in) :: xm
  type(soca_state), intent(inout) :: xa

  ! TODO (G): generate soca readable file from the cice restart

end subroutine soca_soca2cice_changevarinv

! ------------------------------------------------------------------------------
!> fix out of bounds values
!!
subroutine check_ice_bounds(self, geom, xm)
  class(soca_soca2cice), intent(inout) :: self
  type(soca_geom), target, intent(in)  :: geom
  type(soca_state),      intent(inout) :: xm

  type(soca_field), pointer :: aice_ana, hice_ana, hsno_ana
  integer :: i, j
  real(kind=kind_real) :: hice
  real(kind=kind_real) :: hice_max = 8.0

  ! pointers to soca fields (most likely an analysis)
  call xm%get("sea_ice_area_fraction",aice_ana)
  call xm%get("sea_ice_thickness",hice_ana)
  call xm%get("sea_ice_snow_thickness",hsno_ana)

  ! check seaice fraction bounds
  where (aice_ana%val<0_kind_real)
     aice_ana%val = 0_kind_real
  end where
  where (aice_ana%val>1_kind_real)
     aice_ana%val = 1_kind_real
  end where

  ! check seaice thickness bounds
  where (hice_ana%val<0_kind_real)
     hice_ana%val = 0_kind_real
  end where

  ! check snow thickness bounds
  where (hsno_ana%val<0_kind_real)
     hsno_ana%val = 0_kind_real
  end where

end subroutine check_ice_bounds

! ------------------------------------------------------------------------------
!> add seaice to the background
!!
subroutine shuffle_ice(self, geom, xm)
  class(soca_soca2cice), intent(inout) :: self
  type(soca_geom), target, intent(in)  :: geom
  type(soca_state),      intent(inout) :: xm

  real(kind=kind_real) :: aice, seaice_edge
  integer :: i, j, k, n, ii, jj
  type(soca_field), pointer :: t_ana, s_ana, aice_ana
  integer :: minidx(1), nn_max
  integer, allocatable :: idx(:)
  real(kind=kind_real), allocatable :: testmin(:)
  type(cice_state) :: cice_in

  ! Make sure the search tree is smaller than the data size
  nn_max = min(self%cice%agg%n_src, 9)
  allocate(idx(nn_max), testmin(nn_max))

  ! pointers to soca fields (most likely an analysis)
  call xm%get("sea_water_potential_temperature",t_ana)
  call xm%get("sea_water_salinity",s_ana)
  call xm%get("sea_ice_area_fraction",aice_ana)

  call cice_in%copydata(self%cice)
  do i = geom%isc, geom%iec
     do j = geom%jsc, geom%jec

        aice = aice_ana%val(i,j,1)    ! ice fraction analysis

        ! Skip if outside of domain
        if (geom%lat(i,j)>0.0_kind_real) then
          if (.not. self%arctic%shuffle) cycle
          seaice_edge = self%arctic%seaice_edge
        else
          if (.not. self%antarctic%shuffle) cycle
          seaice_edge = self%antarctic%seaice_edge
        endif
        if (self%cice%aice(i,j).gt.seaice_edge) cycle     ! skip if the background has more ice than the threshold
        if (aice.le.0.0_kind_real) cycle                  ! 0 ice analysis is treated elsewhere
        if (self%cice%agg%n_src == 0) cycle               ! skip if there are no points on this task with ice in the background
        ! find neighbors. TODO (G): add constraint for thickness and snow depth as well
        call self%kdtree%closestPoints(geom%lon(i,j), geom%lat(i,j), nn_max, idx)
        do k = 1, nn_max
           testmin(k) = abs(cice_in%aice(self%cice%agg%ij(1, idx(k)), self%cice%agg%ij(2, idx(k))) - aice)
        end do
        minidx = minloc(testmin) ! I know, I rock.
        ii = self%cice%agg%ij(1, idx(minidx(1)))
        jj = self%cice%agg%ij(2, idx(minidx(1)))

        ! update local no ice state with closest non-0 ice state
        self%cice%aice(i, j) = cice_in%aice(ii, jj)
        self%cice%aicen(i, j,:) = cice_in%aicen(ii, jj, :)
        self%cice%vicen(i, j,:) = cice_in%vicen(ii, jj, :)
        self%cice%vsnon(i, j,:) = cice_in%vsnon(ii, jj, :)
        self%cice%apnd(i, j,:) = cice_in%apnd(ii, jj, :)
        self%cice%hpnd(i, j,:) = cice_in%hpnd(ii, jj, :)
        self%cice%ipnd(i, j,:) = cice_in%ipnd(ii, jj, :)
        self%cice%tsfcn(i, j,:) = cice_in%tsfcn(ii, jj, :)

        do k = 1, self%ice_lev
           self%cice%qice(i, j,: , k) = cice_in%qice(ii, jj, :, k)
           self%cice%sice(i, j,: , k) = cice_in%sice(ii, jj, :, k)
        end do
        do k = 1, self%sno_lev
           self%cice%qsno(i, j,: , k) = cice_in%qsno(ii, jj, :, k)
        end do
     end do
  end do

end subroutine shuffle_ice

! ------------------------------------------------------------------------------
!> clean-up the CICE state
!!
subroutine cleanup_ice(self, geom, xm)
  class(soca_soca2cice), intent(inout) :: self
  type(soca_geom), target, intent(in)  :: geom
  type(soca_state),      intent(inout) :: xm

  integer :: i, j, k, n, n_src
  type(soca_field), pointer :: t_ana, s_ana, aice_ana, hice_ana, hsno_ana

  real(kind=kind_real), allocatable :: h_bounds(:)
  real(kind=kind_real), allocatable :: zTin(:), zTsn(:), temp_sno_test

  ! pointers to soca fields (most likely an analysis)
  call xm%get("sea_water_potential_temperature",t_ana)
  call xm%get("sea_water_salinity",s_ana)
  call xm%get("sea_ice_area_fraction",aice_ana)
  call xm%get("sea_ice_thickness", hice_ana)
  call xm%get("sea_ice_snow_thickness", hsno_ana)

  ! get thickness category bounds
  allocate(h_bounds(0:self%ncat))
  call icepack_init_itd(self%ncat, h_bounds) ! TODO (G): move that in setup

  ! TODO (G): re-bin sea-ice that is out of the thickness category

  ! reset sea-ice where ice fraction of the category is 0
  ! and check/fix snow temperature
  allocate(zTin(self%ice_lev), zTsn(self%sno_lev))
  do i = geom%isc, geom%iec
     do j = geom%jsc, geom%jec
        if (aice_ana%val(i,j,1).eq.0.0_kind_real) cycle

        do k = 1, self%ncat
           ! zero out enthalpy if no ice
           if (self%cice%aicen(i,j,k).eq.0.0_kind_real) then
              self%cice%vicen(i,j,k) = 0_kind_real
              self%cice%vsnon(i,j,k) = 0_kind_real

              self%cice%apnd(i,j,k) = 0_kind_real
              self%cice%hpnd(i,j,k) = 0_kind_real
              self%cice%ipnd(i,j,k) = 0_kind_real

              self%cice%qice(i,j,k,:) = 0_kind_real
              self%cice%sice(i,j,k,:) = 0_kind_real
              self%cice%qsno(i,j,k,:) = 0_kind_real

              self%cice%tsfcn(i,j,k) = liquidus_temperature_mush(s_ana%val(i,j,1))
           else
              ! Check snow temperature and adjust if out of wack
              ! for future ref, vertical geometry from top to bottom:
              ! T surface                (tsfcn)
              ! T snow level 1           (from qsno)
              ! ...
              ! T snow level sno_lev
              ! T ice level 1            (from qice)
              ! ...
              ! T ice level ice_lev
              ! Ocean T level 1          (from MOM)
              do n = 1, self%ice_lev
                 zTin(n) = icepack_ice_temperature(self%cice%qice(i,j,k,n), self%cice%sice(i,j,k,n))
              end do
              zTsn(1) = (Lfresh + self%cice%qsno(i,j,k,1)/rhos)/cp_ice
              temp_sno_test = 0.5_kind_real*(self%cice%tsfcn(i,j,k) + zTin(1))

              ! TODO (G): the 2 deg departure check is pulled out of thin ice ... or somethin' move
              !           to config?
              if (abs(zTsn(1)-temp_sno_test).lt.2.0_kind_real) cycle

              ! if we're here, there's a problem with snow, "fix" it!
              self%cice%qsno(i,j,k,1) = rhos*(temp_sno_test*cp_ice - Lfresh)

           end if
        end do
     end do
  end do

  ! re-compute aggregates = analysis that is effectively inserted in the restart
  do i = geom%isc, geom%iec
     do j = geom%jsc, geom%jec
        aice_ana%val(i,j,1) = sum(self%cice%aicen(i,j,:))
        hice_ana%val(i,j,1) = sum(self%cice%vicen(i,j,:))
        hsno_ana%val(i,j,1) = sum(self%cice%vsnon(i,j,:))
     end do
  end do

end subroutine cleanup_ice

! ------------------------------------------------------------------------------
!> add seaice to the background
subroutine prior_dist_rescale(self, geom, xm)
  class(soca_soca2cice), intent(inout) :: self
  type(soca_geom), target, intent(in)  :: geom
  type(soca_state),      intent(inout) :: xm

  real(kind=kind_real) :: alpha, hice, hsno, seaice_edge, rescale_min_hice, rescale_min_hsno
  type(soca_field), pointer :: s_ana, aice_ana, hice_ana, hsno_ana
  integer :: c, i, j

  call xm%get("sea_ice_area_fraction", aice_ana)
  call xm%get("sea_ice_thickness", hice_ana)
  call xm%get("sea_ice_snow_thickness", hsno_ana)
  call xm%get("sea_water_salinity", s_ana)

  do i = geom%isc, geom%iec
     do j = geom%jsc, geom%jec

        if (geom%lat(i,j)>0.0_kind_real) then
          if (.not. self%arctic%rescale_prior) cycle
          seaice_edge = self%arctic%seaice_edge
          rescale_min_hice = self%arctic%rescale_min_hice
          rescale_min_hsno = self%arctic%rescale_min_hsno
        else
          if (.not. self%antarctic%rescale_prior) cycle
          seaice_edge = self%antarctic%seaice_edge
          rescale_min_hice = self%antarctic%rescale_min_hice
          rescale_min_hsno = self%antarctic%rescale_min_hsno
        endif
        if (self%cice%aice(i,j).lt.seaice_edge) cycle ! Only rescale within the icepack

        ! rescale background to match aggregate ice concentration analysis
        alpha = aice_ana%val(i,j,1)/self%cice%aice(i,j)
        self%cice%aice(i,j) = alpha * self%cice%aice(i,j)
        do c = 1, self%ncat
           self%cice%aicen(i,j,c) = alpha*self%cice%aicen(i,j,c)
           self%cice%vicen(i,j,c) = alpha*self%cice%vicen(i,j,c)
           self%cice%vsnon(i,j,c) = alpha*self%cice%vsnon(i,j,c)
        end do

        ! adjust ice volume to match mean cell thickness
        hice = sum(self%cice%vicen(i,j,:))
        if (hice.gt.rescale_min_hice) then
           alpha = hice_ana%val(i,j,1)/hice
           self%cice%vicen(i,j,:) = alpha*self%cice%vicen(i,j,:)
        end if

        ! adjust snow volume to match mean cell thickness
        hsno = sum(self%cice%vsnon(i,j,:))
        if (hsno.gt.rescale_min_hsno) then
           alpha = hsno_ana%val(i,j,1)/hsno
           self%cice%vsnon(i,j,:) = alpha*self%cice%vsnon(i,j,:)
        end if

     end do
  end do

end subroutine prior_dist_rescale

! ------------------------------------------------------------------------------

end module soca_soca2cice_mod
