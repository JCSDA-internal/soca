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
  use hdiag_nicas_mod
implicit none

!> Fortran derived type to hold configuration data for the SOCA background/model covariance
type :: soca_3d_covar_config
  integer :: nx !< Zonal grid dimension
  integer :: ny !< Meridional grid dimension
  real(kind=kind_real)    :: sigma        !< Standard deviation
  real(kind=kind_real)    :: vert_corr    !< Vertical correlation between levels
  !real(kind=kind_real), allocatable :: sqrt_merid(:,:) !< sqrt(meridional correlation matrix)
  !real(kind=kind_real), allocatable :: sqrt_inv_merid(:,:) !< Its inverse
  !real(kind=kind_real), allocatable :: sqrt_zonal(:) !< Spectral weights for sqrt(zonal corr)
  type(hdiag_nicas)    :: nicasB  
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


implicit none
type(c_ptr), intent(in)   :: c_model  !< The configuration
type(soca_geom), intent(in) :: geom     !< Geometry
type(soca_3d_covar_config), intent(inout) :: config !< The covariance structure
!real(kind=kind_real) :: corr_length_scale

!type(hdiag_nicas) :: hdiag



!config%sigma      = config_get_real(c_model,"standard_deviation")
!config%vert_corr  = config_get_real(c_model,"vertical_correlation")
!corr_length_scale = config_get_real(c_model,"horizontal_length_scale")


return
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

! ------------------------------------------------------------------------------

!> Multiply by sqrt(C), where C is a 3d covariance matrix

subroutine soca_3d_covar_sqrt_mult(xincr, xctrl, config)
use iso_c_binding
use kinds
use soca_fields
use type_nam, only: namtype, namcheck
use type_geom, only: geomtype, compute_grid_mesh
use type_bpar, only: bpartype
use type_ndata, only: ndatatype,ndata_dealloc
use type_bdata, only: bdatatype,bdata_alloc
use unstructured_grid_mod
use soca_fields
use nicas_apply_localization, only: apply_localization,apply_localization_sqrt,randomize_localization
use nicas_apply_nicas, only: apply_nicas,apply_nicas_sqrt,apply_nicas_sqrt_ad  

use type_bdata, only: bdatatype
use type_bpar, only: bpartype,bpar_alloc
use model_oops, only: model_oops_coord
use nicas_apply_nicas, only: apply_nicas
use type_randgen,     only: create_randgen
use type_cv, only: cvtype,cv_alloc,cv_random 
use nicas_parameters

implicit none
type(soca_field), intent(inout) :: xincr
type(soca_field), intent(in) :: xctrl
type(soca_3d_covar_config), intent(in) :: config !< covariance config structure

type(namtype)   :: nam                                                            !< Namelist
type(geomtype)  :: geom                                                          !< Geometry
type(bpartype)  :: bpar                                                          !< Block parameters
!type(ndatatype), allocatable :: ndata(:)!(bpar%nb+1)                                             !< NICAS data
!type(ndatatype), allocatable :: ndata(:)
type(ndatatype) :: ndata
type(bdatatype) :: bdata
type(unstructured_grid) :: ug



  

end subroutine soca_3d_covar_sqrt_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C) - Adjoint

subroutine soca_3d_covar_sqrt_mult_ad()

implicit none

! Call to nicas and balance ops ... eventually


end subroutine soca_3d_covar_sqrt_mult_ad

! ------------------------------------------------------------------------------

end module soca_covariance_mod
