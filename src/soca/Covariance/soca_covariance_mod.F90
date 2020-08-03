! (C) Copyright 2017-2020 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Structure holding configuration variables for the 3d error
!! covariance matrices of the SOCA analysis.

module soca_covariance_mod

use atlas_module, only: atlas_fieldset, atlas_field, atlas_real, atlas_integer, atlas_functionspace
use, intrinsic :: iso_c_binding, only : c_char
use fckit_configuration_module, only: fckit_configuration
use random_mod, only: normal_distribution
use oops_variables_mod
use type_bump, only: bump_type
use kinds, only: kind_real
use soca_fields_mod
use soca_increment_mod
use soca_state_mod

use soca_geom_mod, only : soca_geom

implicit none

private
public :: soca_cov, soca_cov_setup, soca_cov_delete, &
          soca_cov_C_mult, soca_cov_sqrt_C_mult

!> Fortran derived type to hold configuration data for the SOCA background/model covariance
type :: soca_pert
  real(kind=kind_real) :: T, S, SSH, AICE, HICE
end type soca_pert

type :: soca_cov
   type(bump_type),     pointer :: ocean_conv(:)  !< Ocean convolution op from bump
   type(bump_type),     pointer :: seaice_conv(:) !< Seaice convolution op from bump
   type(soca_state),    pointer :: bkg            !< Background field (or first guess)
   logical                      :: initialized = .false.
   type(soca_pert)              :: pert_scale
   real(kind=kind_real)         :: ocn_l0
   real(kind=kind_real)         :: ice_l0
   type(oops_variables)         :: vars           !< Apply B to vars
 contains
   procedure :: setup => soca_cov_setup
   procedure :: delete => soca_cov_delete
   procedure :: mult => soca_cov_C_mult
   procedure :: sqrt_C_mult => soca_cov_sqrt_C_mult
end type soca_cov

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Setup for the SOCA model's 3d error covariance matrices (B and Q_i)

!> This routine queries the configuration for the parameters that define the
!! covariance matrix, and stores the relevant values in the
!! error covariance structure.

subroutine soca_cov_setup(self, f_conf, geom, bkg, vars)
  class(soca_cov),        intent(inout) :: self   !< The covariance structure
  type(fckit_configuration), intent(in) :: f_conf !< The configuration
  type(soca_geom),           intent(in) :: geom   !< Geometry
  type(soca_state),  target, intent(in) :: bkg    !< Background
  type(oops_variables),      intent(in) :: vars   !< List of variables

  character(len=3)  :: domain
  integer :: isc, iec, jsc, jec, ivar
  logical :: init_seaice, init_ocean

  ! Setup list of variables to apply B on
  self%vars = vars

  ! Set default ensemble perturbation scales to 1.0, overwrite scales if they exist
  if (.not. f_conf%get("pert_T", self%pert_scale%T))       self%pert_scale%T = 1.0
  if (.not. f_conf%get("pert_S", self%pert_scale%S))       self%pert_scale%S = 1.0
  if (.not. f_conf%get("pert_SSH", self%pert_scale%SSH))   self%pert_scale%SSH = 1.0
  if (.not. f_conf%get("pert_AICE", self%pert_scale%AICE)) self%pert_scale%AICE = 1.0
  if (.not. f_conf%get("pert_HICE", self%pert_scale%HICE)) self%pert_scale%HICE = 1.0

  ! Setup ocean and ice decorrelation length scales
  self%ocn_l0 = 500.0d3
  self%ice_l0 = 500.0d3
  if ( f_conf%has("ocean_corr_scale") ) &
      call f_conf%get_or_die("ocean_corr_scale", self%ocn_l0)
  if ( f_conf%has("ice_corr_scale") ) &
      call f_conf%get_or_die("ice_corr_scale", self%ice_l0)

  ! Associate background
  self%bkg => bkg

  ! Indices for compute domain (no halo)
  isc = bkg%geom%isc ; iec = bkg%geom%iec
  jsc = bkg%geom%jsc ; jec = bkg%geom%jec

  ! Determine what convolution op to initialize
  init_seaice = .false.
  init_ocean = .false.
  do ivar = 1, self%vars%nvars()
     select case(trim(self%vars%variable(ivar)))
     case('cicen','hicen')
        init_seaice = .true.
     case('tocn', 'socn', 'ssh')
        init_ocean = .true.
     end select
  end do

  ! Initialize ocean bump if tocn or socn or ssh are in self%vars
  domain = 'ocn'
  allocate(self%ocean_conv(1))
  if (init_ocean) then
     call soca_bump_correlation(self, self%ocean_conv(1), geom, f_conf, domain)
  end if

  ! Initialize seaice bump if cicen or hicen are in self%vars
  domain = 'ice'
  allocate(self%seaice_conv(1))
  if (init_seaice) then
     call soca_bump_correlation(self, self%seaice_conv(1), geom, f_conf, domain)
  end if

  self%initialized = .true.

end subroutine soca_cov_setup

! ------------------------------------------------------------------------------

!> Delete for the SOCA model's 3d error covariance matrices

subroutine soca_cov_delete(self)
  class(soca_cov), intent(inout) :: self       !< The covariance structure

  call self%ocean_conv(1)%dealloc()
  call self%seaice_conv(1)%dealloc()
  deallocate(self%ocean_conv)
  deallocate(self%seaice_conv)
  nullify(self%bkg)
  self%initialized = .false.

end subroutine soca_cov_delete

! ------------------------------------------------------------------------------

subroutine soca_cov_C_mult(self, dx)
  class(soca_cov),      intent(inout) :: self !< The covariance structure
  type(soca_increment), intent(inout) :: dx   !< Input: Increment
                                          !< Output: C dx
  integer :: i, z
  type(soca_field), pointer :: field
  type(bump_type), pointer :: conv

  do i = 1, self%vars%nvars()
    if (.not. dx%has(self%vars%variable(i))) cycle ! why is this sometimes getting an "empty" list with "none" in it?
    call dx%get(trim(self%vars%variable(i)), field)

    ! TODO remove the hardcoded variables
    ! ice or ocean convolution ?
    select case(field%name)
    case ('tocn', 'socn', 'ssh', 'sw', 'lw', 'lhf', 'shf', 'us')
      conv => self%ocean_conv(1)
    case ('hicen','cicen')
      conv => self%seaice_conv(1)
    case default
      cycle
    end select

    ! apply convolution on each level
    do z = 1, field%nz
      call soca_2d_convol(field%val(:,:,z), conv, dx%geom)
    end do
  end do
end subroutine soca_cov_C_mult

! ------------------------------------------------------------------------------

subroutine soca_cov_sqrt_C_mult(self, dx)
  class(soca_cov),      intent(inout) :: self !< The covariance structure
  type(soca_increment), intent(inout) :: dx   !< Input: Increment
                                          !< Output: C^1/2 dx
  integer :: i, z
  type(soca_field), pointer :: field
  real(kind=kind_real) :: scale
  type(bump_type), pointer :: conv

  do i = 1, self%vars%nvars()
    conv => null()
    call dx%get(trim(self%vars%variable(i)), field)

    select case(field%name)
    case('tocn')
      scale = self%pert_scale%T
      conv => self%ocean_conv(1)
    case ('socn')
      scale = self%pert_scale%S
      conv => self%ocean_conv(1)
    case ('ssh', 'sw', 'lw', 'lhf', 'shf', 'us')
      scale = self%pert_scale%SSH
      conv => self%ocean_conv(1)
    case('cicen')
      scale = self%pert_scale%AICE
      conv => self%seaice_conv(1)
    case ('hicen')
      scale = self%pert_scale%HICE
      conv => self%seaice_conv(1)
    end select

    if (associated(conv)) then
      do z = 1,field%nz
        call soca_2d_sqrt_convol(field%val(:,:,z), conv, dx%geom, scale)
      end do
    end if

  end do
end subroutine soca_cov_sqrt_C_mult

! ------------------------------------------------------------------------------

subroutine soca_bump_correlation(self, horiz_convol, geom, f_conf, domain)
  class(soca_cov),        intent(inout) :: self   !< The covariance structure
  type(bump_type),        intent(inout) :: horiz_convol
  type(soca_geom),           intent(in) :: geom
  type(fckit_configuration), intent(in) :: f_conf !< Handle to configuration
  character(len=3),          intent(in) :: domain

  integer, pointer :: int_ptr_2(:,:)
  real(kind=kind_real), pointer :: real_ptr_1(:), real_ptr_2(:,:)
  type(atlas_functionspace) :: afunctionspace
  type(atlas_fieldset) :: afieldset, rh, rv
  type(atlas_field) :: afield
  type(fckit_configuration) :: f_grid

  ! Grid setup
  f_grid = fckit_configuration()
  call f_grid%set('verbosity', 'none')
  call f_grid%set('prefix', domain)
  call f_grid%set('nl', 1)
  call f_grid%set('nv', 1)
  call f_grid%set('variables', ['var'])
  call f_grid%set('nts', 1)
  call f_grid%set('timeslots', ['00'])

  ! Wrap functionspace
  afunctionspace = atlas_functionspace(geom%afunctionspace%c_ptr())

  ! Geometry fieldset setup
  afieldset = atlas_fieldset()

  ! Add area
  afield = geom%afunctionspace%create_field(name='area', kind=atlas_real(kind_real), levels=0)
  call afield%data(real_ptr_1)
  real_ptr_1 = pack(geom%cell_area(geom%isc:geom%iec,geom%jsc:geom%jec),.true.)
  call afieldset%add(afield)
  call afield%final()

  ! Add vertical unit
  afield = geom%afunctionspace%create_field(name='vunit', kind=atlas_real(kind_real), levels=1)
  call afield%data(real_ptr_2)
  real_ptr_2(1,:) = 1.0
  call afieldset%add(afield)
  call afield%final()

  ! Add geographical mask
  afield = geom%afunctionspace%create_field(name='gmask', kind=atlas_integer(kind(0)), levels=1)
  call afield%data(int_ptr_2)
  int_ptr_2(1,:) = int(pack(geom%mask2d(geom%isc:geom%iec,geom%jsc:geom%jec),.true.))
  call afieldset%add(afield)
  call afield%final()

  ! Create BUMP object
  call horiz_convol%create(geom%f_comm,afunctionspace,afieldset,f_conf,f_grid)

  if (horiz_convol%nam%new_nicas) then
     ! rh
     rh = atlas_fieldset()
     afield = geom%afunctionspace%create_field('var_00',kind=atlas_real(kind_real),levels=0)
     call rh%add(afield)
     call afield%data(real_ptr_1)
     if (domain.eq.'ocn') real_ptr_1 = self%ocn_l0+pack(geom%rossby_radius(geom%isc:geom%iec,geom%jsc:geom%jec), .true.)
     if (domain.eq.'ice') real_ptr_1 = self%ice_l0
     call afield%final()

     ! rv
     rv = atlas_fieldset()
     afield = geom%afunctionspace%create_field('var_00',kind=atlas_real(kind_real),levels=0)
     call rv%add(afield)
     call afield%data(real_ptr_1)
     real_ptr_1 = 1.0
     call afield%final()

     ! Copy length-scales into BUMP
     call horiz_convol%set_parameter('cor_rh', rh)
     call horiz_convol%set_parameter('cor_rv', rv)

     ! Clean up
     call rh%final()
     call rv%final()
  end if

  ! Run BUMP drivers
  call horiz_convol%run_drivers()

end subroutine soca_bump_correlation

! ------------------------------------------------------------------------------

subroutine soca_2d_convol(dx, horiz_convol, geom)
  real(kind=kind_real), intent(inout) :: dx(:,:)
  type(bump_type),      intent(inout) :: horiz_convol
  type(soca_geom),         intent(in) :: geom

  type(atlas_fieldset) :: tmp_incr

  ! Allocate ATLAS tmp_increment and make copy of dx
  call geom%struct2atlas(dx(:,:), tmp_incr)

  ! Apply 2D convolution
  call horiz_convol%apply_nicas(tmp_incr)

  ! Copy ATLAS tmp_incr to structured dx
  call geom%atlas2struct(dx(:,:), tmp_incr)

  ! Clean up
  call tmp_incr%final()

end subroutine soca_2d_convol

! ------------------------------------------------------------------------------

subroutine soca_2d_sqrt_convol(dx, horiz_convol, geom, pert_scale)
  real(kind=kind_real), intent(inout) :: dx(:,:)
  type(bump_type),      intent(inout) :: horiz_convol
  type(soca_geom),         intent(in) :: geom
  real(kind=kind_real),    intent(in) :: pert_scale

  type(atlas_fieldset) :: tmp_incr
  real(kind=kind_real), allocatable :: pcv(:)
  integer, parameter :: rseed = 1 ! constant for reproducability of tests
                                  ! TODO: pass seed through config
  integer :: nn

  ! Allocate ATLAS tmp_increment and make copy of dx
  call geom%struct2atlas(dx(:,:), tmp_incr)

  ! Get control variable size
  call horiz_convol%get_cv_size(nn)
  allocate(pcv(nn))
  pcv = 0.0_kind_real
  call normal_distribution(pcv, 0.0_kind_real, 1.0_kind_real, rseed)
  pcv = pert_scale * pcv

  ! Apply C^1/2
  call horiz_convol%apply_nicas_sqrt(pcv, tmp_incr)

  ! Back to structured grid
  call geom%atlas2struct(dx(:,:), tmp_incr)

  ! Clean up
  deallocate(pcv)
  call tmp_incr%final()

end subroutine soca_2d_sqrt_convol

end module soca_covariance_mod
