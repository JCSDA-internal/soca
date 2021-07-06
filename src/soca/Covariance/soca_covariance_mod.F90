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
use fckit_log_module
use random_mod, only: normal_distribution
use oops_variables_mod
use type_bump, only: bump_type
use tools_func, only: gau2gc
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
type :: soca_cov
   type(bump_type),     pointer :: conv(:)        !< convolution op from bump
   type(soca_state),    pointer :: bkg            !< Background field (or first guess)
   type(oops_variables)         :: vars           !< Apply B to vars

   real(kind=kind_real), allocatable :: pert_scale(:) !< index matches "vars"
   type(oops_variables), allocatable :: conv_vars(:)  !< index mathces "conv"

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

  type(fckit_configuration) :: f_conf2
  type(fckit_configuration), allocatable :: f_conf_list(:)
  character(len=:), allocatable :: domain_vars(:)
  character(len=:), allocatable :: domain
  integer :: i, isc, iec, jsc, jec, ivar

  ! Setup list of variables to apply B on
  self%vars = vars

  ! get perturbation scales (or set to 1.0)
  allocate(self%pert_scale(self%vars%nvars()))
  self%pert_scale = 1.0
  if (f_conf%get("perturbation scales", f_conf2)) then
    do ivar=1,self%vars%nvars()
      if ( .not. f_conf2%get(self%vars%variable(ivar), self%pert_scale(ivar))) then
        if (geom%f_comm%rank() == 0) call fckit_log%warning( &
          "WARNING: no pertubation scale given for '"  //trim(self%vars%variable(ivar)) &
           // "' using default of 1.0")
      end if
    end do
  end if

  ! Associate background
  self%bkg => bkg

  ! Indices for compute domain (no halo)
  isc = bkg%geom%isc ; iec = bkg%geom%iec
  jsc = bkg%geom%jsc ; jec = bkg%geom%jec

  ! Initialize bump
  call f_conf%get_or_die("bump", f_conf2)
  call f_conf%get_or_die("correlation", f_conf_list)
  allocate(self%conv(size(f_conf_list)))
  allocate(self%conv_vars(size(self%conv)))
  do i=1,size(f_conf_list)
    call f_conf_list(i)%get_or_die("name", domain)
    call f_conf_list(i)%get_or_die("vars", domain_vars)
    self%conv_vars(i) = oops_variables()
    call self%conv_vars(i)%push_back(domain_vars)

    call soca_bump_correlation(self, self%conv(i), geom, f_conf2, f_conf_list(i), domain)
  end do

end subroutine soca_cov_setup

! ------------------------------------------------------------------------------

!> Delete for the SOCA model's 3d error covariance matrices

subroutine soca_cov_delete(self)
  class(soca_cov), intent(inout) :: self       !< The covariance structure

  deallocate(self%conv)
  deallocate(self%conv_vars)
  deallocate(self%pert_scale)
  nullify(self%bkg)

end subroutine soca_cov_delete

! ------------------------------------------------------------------------------

subroutine soca_cov_C_mult(self, dx)
  class(soca_cov),      intent(inout) :: self !< The covariance structure
  type(soca_increment), intent(inout) :: dx   !< Input: Increment
                                          !< Output: C dx
  integer :: i, z, j, k
  type(soca_field), pointer :: field
  type(bump_type), pointer :: conv

  do i = 1, self%vars%nvars()
    if (.not. dx%has(self%vars%variable(i))) cycle ! why is this sometimes getting an "empty" list with "none" in it?
    call dx%get(trim(self%vars%variable(i)), field)

    ! determine which horizontal convolution to use
    nullify(conv)
    outer: do j=1,size(self%conv_vars)
      do k=1,self%conv_vars(j)%nvars()
        if (self%conv_vars(j)%variable(k) == field%name) then
          conv => self%conv(j)
          exit outer
        end if
      end do
    end do outer
    if ( .not. associated(conv)) then
      call abor1_ftn("No valid bump operator found for field '"//field%name//"'")
    end if

    ! safety check to make sure field is on the correct grid
    ! that the convolution was built with

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
  integer :: i, z, j, k
  type(soca_field), pointer :: field
  real(kind=kind_real) :: scale
  type(bump_type), pointer :: conv

  do i = 1, self%vars%nvars()
    conv => null()
    call dx%get(trim(self%vars%variable(i)), field)

    ! find matching index in self%vars and get the perturbation scale
    if (.not. allocated(self%pert_scale)) then
      call abor1_ftn("ERROR: cannot use sqrt_C_mult if no perturbation scales given")
    endif
    do j=1,self%vars%nvars()
      if (self%vars%variable(j) == field%name) exit
    end do
    scale = self%pert_scale(j)

    ! determine which horizontal convolution to use
    nullify(conv)
    outer: do j=1,size(self%conv_vars)
      do k=1,self%conv_vars(j)%nvars()
        if (self%conv_vars(j)%variable(k) == field%name) then
          conv => self%conv(j)
          exit outer
        end if
      end do
    end do outer
    if ( .not. associated(conv)) then
      call abor1_ftn("No valid bump operator found for "//field%name)
    end if

    ! apply convolution
    do z = 1,field%nz
      call soca_2d_sqrt_convol(field%val(:,:,z), conv, dx%geom, scale)
    end do

  end do
end subroutine soca_cov_sqrt_C_mult

! ------------------------------------------------------------------------------

subroutine soca_bump_correlation(self, horiz_convol, geom, f_conf_bump, f_conf_domain, domain)
  class(soca_cov),        intent(inout) :: self   !< The covariance structure
  type(bump_type),        intent(inout) :: horiz_convol
  type(soca_geom),           intent(in) :: geom
  type(fckit_configuration), intent(in) :: f_conf_bump, f_conf_domain
  character(len=3),          intent(in) :: domain

  integer :: i
  integer, pointer :: int_ptr_2(:,:)
  real(kind=kind_real), pointer :: real_ptr_1(:), real_ptr_2(:,:)
  real(kind=kind_real), allocatable :: lats(:), area(:)
  type(atlas_functionspace) :: afunctionspace
  type(atlas_fieldset) :: afieldset, rh, rv
  type(atlas_field) :: afield
  type(fckit_configuration) :: f_grid
  real(kind=kind_real) :: r_base, r_mult, r_min, r_max, r_min_grid


  ! Grid setup
  f_grid = fckit_configuration()
  call f_grid%set('prefix', domain)
  call f_grid%set('nl', 1)
  call f_grid%set('nv', 1)
  call f_grid%set('variables', ['var'])

  ! Wrap functionspace
  afunctionspace = atlas_functionspace(geom%afunctionspace%c_ptr())

  ! Geometry fieldset setup
  afieldset = atlas_fieldset()

  lats = pack(geom%lat(geom%isc:geom%iec,geom%jsc:geom%jec),.true.)

  ! Add area
  afield = geom%afunctionspace%create_field(name='area', kind=atlas_real(kind_real), levels=0)
  call afield%data(real_ptr_1)
  area = pack(geom%cell_area(geom%isc:geom%iec,geom%jsc:geom%jec),.true.)
  real_ptr_1 = area
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
  call horiz_convol%create(geom%f_comm,afunctionspace,afieldset,f_conf_bump,f_grid)

  if (horiz_convol%nam%new_nicas) then
    ! get parameters for correlation lengths
    if (.not. f_conf_domain%get('base value', r_base))    r_base = 0.0
    if (.not. f_conf_domain%get('rossby mult', r_mult))   r_mult = 0.0
    if (.not. f_conf_domain%get('min grid mult', r_min_grid))  r_min_grid = 1.0
    if (.not. f_conf_domain%get('min value', r_min))      r_min  = 0.0
    if (.not. f_conf_domain%get('max value', r_max))      r_max  = huge(r_max)

    ! rh is calculated as follows :
    ! 1) rh = "base value" + rossby_radius * "rossby mult"
    ! 2) minimum value of "min grid mult" * grid_size is imposed
    ! 3) min/max are imposed based on "min value" and "max value"
    ! 4) converted from a gaussian sigma to Gaspari-Cohn cutoff distance
    rh = atlas_fieldset()
    afield = geom%afunctionspace%create_field('var',kind=atlas_real(kind_real),levels=0)
    call rh%add(afield)
    call afield%data(real_ptr_1)
    real_ptr_1 = r_base + r_mult*pack(geom%rossby_radius(geom%isc:geom%iec,geom%jsc:geom%jec), .true.)
    ! min based on grid size
    if (r_min_grid .gt. 0.0) then
      real_ptr_1 = max(real_ptr_1,  sqrt(area)*r_min_grid )
    end if
    real_ptr_1 = min(r_max, real_ptr_1)
    real_ptr_1 = max(r_min, real_ptr_1)
    real_ptr_1 = real_ptr_1 * gau2gc ! convert from gaussian sigma to
                                     ! Gaspari-Cohn half width
    call afield%final()

     ! rv
     rv = atlas_fieldset()
     afield = geom%afunctionspace%create_field('var',kind=atlas_real(kind_real),levels=0)
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
