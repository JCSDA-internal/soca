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

use icepack_itd, only: icepack_init_itd, cleanup_itd
use icepack_warnings, only: icepack_warnings_flush, icepack_warnings_aborted
use icepack_tracers, only: icepack_init_tracer_indices, nt_tsfc, nt_qice, nt_qsno, nt_sice
use icepack_parameters, only: icepack_init_parameters, icepack_recompute_constants
use icepack_parameters, only: ktherm, heat_capacity
use icepack_mushy_physics, only: liquidus_temperature_mush
use icepack_therm_shared, only: l_brine

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
   integer :: ncat, ni, nj, ice_lev, sno_lev, shuffle_n
   real(kind=kind_real) :: dt
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
  type(soca_field), pointer :: s_ana, aice_ana
  integer :: minidx(1), nn_max
  integer, allocatable :: idx(:)
  real(kind=kind_real), allocatable :: testmin(:)
  type(cice_state) :: cice_in

  ! Make sure the search tree is smaller than the data size
  nn_max = min(self%cice%agg%n_src, self%shuffle_n)
  allocate(idx(nn_max), testmin(nn_max))

  ! pointers to soca fields (most likely an analysis)
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
        if (aice.le.0.0_kind_real) then
           self%cice%aicen(i,j,:) = 0_kind_real
           self%cice%vicen(i,j,:) = 0_kind_real
           self%cice%vsnon(i,j,:) = 0_kind_real
           self%cice%apnd(i,j,:) = 0_kind_real
           self%cice%hpnd(i,j,:) = 0_kind_real
           self%cice%ipnd(i,j,:) = 0_kind_real
           self%cice%qice(i,j,:,:) = 0_kind_real
           self%cice%sice(i,j,:,:) = 0_kind_real
           self%cice%qsno(i,j,:,:) = 0_kind_real
           self%cice%tsfcn(i,j,:) = liquidus_temperature_mush(s_ana%val(i,j,1))
        endif
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

  integer :: i, j, k, ntracers
  integer :: nt_tsfc_in, nt_qice_in, nt_qsno_in, nt_sice_in
  real(kind=kind_real) :: aice, aice0
  type(soca_field), pointer :: t_ana, s_ana, aice_ana, hice_ana, hsno_ana
  real(kind=kind_real), allocatable :: h_bounds(:)
  real(kind=kind_real), allocatable :: tracers(:,:)   ! (ntracers, ncat)
  logical, allocatable :: first_ice(:)                ! (ncat) ! For bgc and S tracers. set to true if zapping ice.
  integer, allocatable :: trcr_depend(:)              ! (ntracers), = 0 for aicen tracers, 1 for vicen, 2 for vsnon
  real(kind=kind_real), allocatable :: trcr_base(:,:) ! (ntracers, 3);  = 0 or 1 depending on tracer dependency
                                                      ! argument 2:  (1) aice, (2) vice, (3) vsno
  integer, allocatable :: n_trcr_strata(:)            ! number of underlying tracer layers
  integer, allocatable :: nt_strata(:,:)              ! indices of underlying tracer layers
  real(kind=kind_real), allocatable :: fiso_ocn(:)    ! isotope flux to ocean, not used as long as icepack's
                                                      ! tr_iso is false (default), so not initialized here.

  ! pointers to soca fields (most likely an analysis)
  call xm%get("sea_water_potential_temperature",t_ana)
  call xm%get("sea_water_salinity",s_ana)
  call xm%get("sea_ice_area_fraction",aice_ana)
  call xm%get("sea_ice_thickness", hice_ana)
  call xm%get("sea_ice_snow_thickness", hsno_ana)

  ! get thickness category bounds
  allocate(h_bounds(0:self%ncat))
  call icepack_init_itd(self%ncat, h_bounds) ! TODO (G): move that in setup
  ! initialize tracers (ice/snow temperature, ice and snow enthalpies, ice salinity)
  ntracers = 1+2*self%ice_lev+self%sno_lev
  allocate(tracers(ntracers, self%ncat))
  allocate(trcr_depend(ntracers), trcr_base(ntracers, 3))
  allocate(n_trcr_strata(ntracers), nt_strata(ntracers, 2))
  n_trcr_strata(:) = 0
  nt_strata(:,:) = 0
  trcr_base(:, :) = 0.0

  ! ice/snow surface temperature: ice area tracer
  nt_tsfc_in = 1
  ntracers = 1
  trcr_depend(nt_tsfc_in) = 0
  trcr_base(nt_tsfc_in, 1) = 1.0

  ! ice enthalpy: ice volume tracer
  nt_qice_in = ntracers + 1
  ntracers = ntracers + self%ice_lev
  do k = 1, self%ice_lev
   trcr_depend(nt_qice_in + k - 1) = 1
   trcr_base(nt_qice_in + k - 1, 2) = 1.0
  enddo

  ! snow enthalpy: snow volume tracer
  nt_qsno_in = ntracers + 1
  ntracers = ntracers + self%sno_lev
  do k = 1, self%sno_lev
    trcr_depend(nt_qsno_in + k - 1) = 2
    trcr_base(nt_qsno_in + k - 1, 3) = 1.0
  enddo

  ! ice salinity: ice volume tracer
  nt_sice_in = ntracers + 1
  ntracers = ntracers + self%ice_lev
  do k = 1, self%ice_lev
   trcr_depend(nt_sice_in + k - 1) = 1
   trcr_base(nt_sice_in + k - 1, 2) = 1.0
  enddo

  call icepack_init_tracer_indices(nt_tsfc_in=nt_tsfc_in, nt_qice_in=nt_qice_in, &
                                   nt_qsno_in=nt_qsno_in, nt_sice_in=nt_sice_in)
  allocate(first_ice(self%ncat))
  first_ice(:) = .true.

  do i = geom%isc, geom%iec
     do j = geom%jsc, geom%jec
        ! setup tracers at this gridpoint
        tracers(nt_tsfc,:) = self%cice%tsfcn(i,j,:)
        do k = 1, self%ice_lev
          tracers(nt_qice+k-1, :) = self%cice%qice(i,j,:,k)
          tracers(nt_sice+k-1, :) = self%cice%sice(i,j,:,k)
        enddo
        do k = 1, self%sno_lev
          tracers(nt_qsno+k-1, :) = self%cice%qsno(i,j,:,k)
        enddo
        ! call icepack_cleanup_itd: rebins thickness categories if necessary,
        ! eliminates very small ice areas while conserving mass and energy
        call cleanup_itd(self%dt, ntracers, self%ice_lev, self%sno_lev, self%ncat, &
                         h_bounds, self%cice%aicen(i,j,:), tracers, &
                         self%cice%vicen(i,j,:), self%cice%vsnon(i,j,:), &
                         ! ice and total water concentration are computed in the call using aicen
                         aice, aice0, &
                         ! number of aerosol tracers, bio tracers, bio layers
                         0, 0, 0, &
                         ! aerosol flag, topo pond flag, heat capacity flag, flag for zapping ice for bgc and s tracers
                         .false., .false., .true., first_ice, &
                         ! tracer indices and sizes used in rebinning
                         trcr_depend, trcr_base, n_trcr_strata, nt_strata, &
                         fiso_ocn=fiso_ocn)
        ! put tracers back
        self%cice%tsfcn(i,j,:) = tracers(nt_tsfc,:)
        do k = 1, self%ice_lev
          self%cice%qice(i,j,:,k) = tracers(nt_qice+k-1, :)
          self%cice%sice(i,j,:,k) = tracers(nt_sice+k-1, :)
        enddo
        do k = 1, self%sno_lev
          self%cice%qsno(i,j,:,k) = tracers(nt_qsno+k-1, :)
        enddo
        call icepack_warnings_flush(6)
        if (icepack_warnings_aborted()) then
           call abor1_ftn("Soca2Cice: icepack aborted during cleanup_itd")
        endif
        ! re-compute aggregates = analysis that is effectively inserted in the restart
        aice_ana%val(i,j,1) = sum(self%cice%aicen(i,j,:))
        hice_ana%val(i,j,1) = sum(self%cice%vicen(i,j,:))
        hsno_ana%val(i,j,1) = sum(self%cice%vsnon(i,j,:))
     end do
  end do

  deallocate(h_bounds, tracers, trcr_depend, trcr_base, n_trcr_strata, nt_strata, first_ice)

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
        if (self%cice%aice(i,j).le.seaice_edge) cycle ! Only rescale within the icepack

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
