! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Structure holding configuration variables for the 3d error
!! covariance matrices of the SOCA analysis.

module soca_covariance_mod

  use kinds
  implicit none

  !> Fortran derived type to hold configuration data for the SOCA background/model covariance
  type :: soca_3d_covar_config
     real(kind=kind_real) :: Lx=1.0, Ly=.5
     real(kind=kind_real) :: sig_sic
     real(kind=kind_real) :: sig_sit
     real(kind=kind_real) :: sig_ssh
  end type soca_3d_covar_config

#define LISTED_TYPE soca_3d_covar_config

  !> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_3d_cov_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "util/linkedList_c.f"
  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

  !> Setup for the SOCA model's 3d error covariance matrices (B and Q_i)

  !> This routine queries the configuration for the parameters that define the
  !! covariance matrix, and stores the relevant values in the
  !! error covariance structure.

  subroutine soca_3d_covar_setup(c_model, geom, config)

    use soca_constants
    use soca_geom_mod
    use iso_c_binding
    use config_mod
    use fft_mod
    use kinds
    use fckit_log_module, only : fckit_log
    use soca_interph_mod
    
    implicit none
    type(c_ptr), intent(in)   :: c_model  !< The configuration
    type(soca_geom), intent(in) :: geom     !< Geometry
    type(soca_3d_covar_config), intent(inout) :: config !< The covariance structure
    real(kind=kind_real) :: corr_length_scale
    type(soca_hinterp), pointer :: horiz_convol_p
    
    config%sig_sic      = config_get_real(c_model,"sig_sic")
    config%sig_sit      = config_get_real(c_model,"sig_sit")
    config%sig_ssh      = config_get_real(c_model,"sig_ssh")

    call initialize_convolh(geom, horiz_convol_p)
    
  end subroutine soca_3d_covar_setup

  ! ------------------------------------------------------------------------------

  !> Delete for the SOCA model's 3d error covariance matrices
!!$
  subroutine soca_3d_covar_delete(c_key_conf)

    use iso_c_binding

    implicit none
    integer(c_int), intent(inout) :: c_key_conf !< The model covariance structure

    type(soca_3d_covar_config), pointer :: conf !< covar structure

    !call soca_3d_cov_registry%get(c_key_conf, conf)

    !deallocate(conf%sqrt_zonal)
    !deallocate(conf%sqrt_merid)
    !deallocate(conf%sqrt_inv_merid)
    !call soca_3d_cov_registry%remove(c_key_conf)

  end subroutine soca_3d_covar_delete

  !> Multiply by sqrt(C), where C is a 3d covariance matrix

  subroutine soca_3d_covar_sqrt_mult(dx, sqrtCdx, config)
    use iso_c_binding
    use kinds
    use soca_fields
    use soca_interph_mod
    
    implicit none
    type(soca_field), intent(inout)        :: sqrtCdx         !< Full C^1/2 applied to dx
    type(soca_field), intent(in)           :: dx              !< State space increment
    type(soca_3d_covar_config), intent(in) :: config          !< covariance config structure
    type(soca_hinterp), pointer            :: horiz_convol_p  !< pointer to convolution operator
    real(kind=kind_real), allocatable      :: tmp_incr(:)

    call copy(sqrtCdx, dx)
    
    ! sqrtCdx=C.dx
    call initialize_convolh(sqrtCdx%geom, horiz_convol_p)
    allocate(tmp_incr(size(dx%ssh,1)*size(dx%ssh,2)))
    tmp_incr=reshape(dx%ssh,(/size(dx%ssh,1)*size(dx%ssh,2)/))
    call horiz_convol_p%interp_apply(sqrtCdx%ssh, tmp_incr)
    sqrtCdx%ssh=reshape(tmp_incr,(/size(dx%ssh,1),size(dx%ssh,2)/))
    
  end subroutine soca_3d_covar_sqrt_mult

  ! ------------------------------------------------------------------------------

  !> Multiply streamfunction by sqrt(C) - Adjoint

  subroutine soca_3d_covar_sqrt_mult_ad(dx, sqrtCTdx, config)
    use iso_c_binding
    use kinds
    use soca_fields
    use soca_interph_mod
    
    implicit none
    type(soca_field), intent(in)           :: dx
    type(soca_field), intent(inout)        :: sqrtCTdx
    type(soca_3d_covar_config), intent(in) :: config !< covariance config structure
    type(soca_hinterp), pointer            :: horiz_convol_p
    real(kind=kind_real), allocatable      :: tmp_incr(:)
    real(kind=kind_real)      :: crap
    
    call copy(sqrtCTdx, dx)
    
    call initialize_convolh(dx%geom, horiz_convol_p)
    allocate(tmp_incr(size(sqrtCTdx%ssh,1)*size(sqrtCTdx%ssh,2)))
    tmp_incr=reshape(dx%ssh,(/size(dx%ssh,1)*size(dx%ssh,2)/))
    call horiz_convol_p%interpad_apply(sqrtCTdx%ssh,tmp_incr)

  end subroutine soca_3d_covar_sqrt_mult_ad

  ! ------------------------------------------------------------------------------
  
  subroutine soca_3d_covar_D_mult(Ddx, config)
    use iso_c_binding
    use kinds
    use soca_fields
    use soca_interph_mod
    
    implicit none
    type(soca_field), intent(inout)        :: Ddx             !< D applied to dx
    type(soca_3d_covar_config), intent(in) :: config          !< covariance config structure

    Ddx%cicen=config%sig_sic*Ddx%cicen
    Ddx%hicen=config%sig_sit*Ddx%hicen
    Ddx%ssh=config%sig_ssh*Ddx%ssh    
    
  end subroutine soca_3d_covar_D_mult

  ! ------------------------------------------------------------------------------

  subroutine initialize_convolh(geom, horiz_convol_p)
    use ufo_locs_mod  
    use soca_interph_mod
    use soca_geom_mod
    use soca_interph_mod
    
    implicit none

    type(soca_geom), intent(in)              :: geom
    type(soca_hinterp), pointer, intent(out) :: horiz_convol_p

    real(kind=kind_real), allocatable :: lon(:), lat(:)
    logical, save :: convolh_initialized = .false.
    type(soca_hinterp), save, target :: horiz_convol
    integer :: n, nn=4
    character(len=3) :: wgt_type='avg'

    if (.NOT.convolh_initialized) then
       n = size(geom%ocean%lon,1)*size(geom%ocean%lon,2)
       allocate(lon(n),lat(n))
       call horiz_convol%interp_init(n,nn=nn,wgt_type=wgt_type)
       call horiz_convol%interp_compute_weight(geom%ocean%lon,&
            &                                     geom%ocean%lat,&
            &                                     lon,&
            &                                     lat)
       convolh_initialized = .true.
       deallocate(lon,lat)
    end if
    horiz_convol_p => horiz_convol

  end subroutine initialize_convolh

  ! ------------------------------------------------------------------------------

end module soca_covariance_mod
