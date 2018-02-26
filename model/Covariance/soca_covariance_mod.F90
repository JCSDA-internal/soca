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
  integer :: nx !< Zonal grid dimension
  integer :: ny !< Meridional grid dimension
  real(kind=kind_real)    :: sigma        !< Standard deviation
  real(kind=kind_real)    :: vert_corr    !< Vertical correlation between levels
  real(kind=kind_real), allocatable :: sqrt_merid(:,:) !< sqrt(meridional correlation matrix)
  real(kind=kind_real), allocatable :: sqrt_inv_merid(:,:) !< Its inverse
  real(kind=kind_real), allocatable :: sqrt_zonal(:) !< Spectral weights for sqrt(zonal corr)
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
real(kind=kind_real) :: corr_length_scale

!config%nx         = geom%nx
!config%ny         = geom%ny
config%sigma      = config_get_real(c_model,"standard_deviation")
config%vert_corr  = config_get_real(c_model,"vertical_correlation")
corr_length_scale = config_get_real(c_model,"horizontal_length_scale")

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

!> Multiply by sqrt(C), where C is a 3d covariance matrix

subroutine soca_3d_covar_sqrt_mult(xincr, xctrl, config)
use iso_c_binding
use kinds
use soca_fields
use type_nam, only: namtype, namcheck
use type_geom, only: geomtype, compute_grid_mesh
use type_bpar, only: bpartype
use type_ndata, only: ndatatype,ndata_dealloc
use unstructured_grid_mod
use soca_fields
use nicas_apply_localization, only: apply_localization,apply_localization_from_sqrt,randomize_localization
use tools_const, only: deg2rad
use type_bdata, only: bdatatype
use type_bpar, only: bpartype,bpar_alloc
use model_oops, only: model_oops_coord
use nicas_apply_nicas, only: apply_nicas
use type_randgen,     only: create_randgen

implicit none
type(soca_field), intent(inout) :: xincr
type(soca_field), intent(in) :: xctrl
type(soca_3d_covar_config), intent(in) :: config !< covariance config structure

type(namtype)   :: nam                                                            !< Namelist
type(geomtype)  :: geom                                                          !< Geometry
type(bpartype)  :: bpar                                                          !< Block parameters
!type(ndatatype), allocatable :: ndata(:)!(bpar%nb+1)                                             !< NICAS data
type(ndatatype) :: ndata
type(unstructured_grid) :: ug

!nicas stuff
integer :: nc0a, nl0, nv, nts
real(kind=kind_real), allocatable :: lon(:), lat(:), area(:), vunit(:), rndnum(:)
integer, allocatable :: imask(:,:)

!Grid stuff
integer :: isc, iec, jsc, jec, jjj, jz, il

!Get indices for compute domain (no halo)
isc = xctrl%geom%ocean%G%isc
iec = xctrl%geom%ocean%G%iec    
jsc = xctrl%geom%ocean%G%jsc
jec = xctrl%geom%ocean%G%jec
    
nv = xctrl%geom%ocean%ncat
nl0 = 1
nts = 1
nc0a = (iec - isc + 1) * (jec - jsc + 1 )

allocate( lon(nc0a), lat(nc0a), area(nc0a) )
allocate( vunit(nl0) )
allocate( imask(nc0a, nl0) )    


lon = deg2rad*reshape( xctrl%geom%ocean%lon(isc:iec, jsc:jec), (/nc0a/) )
lat = deg2rad*reshape( xctrl%geom%ocean%lat(isc:iec, jsc:jec), (/nc0a/) ) 

area = reshape( xctrl%geom%ocean%cell_area(isc:iec, jsc:jec), (/nc0a/) )

do jz = 1, nl0       
   vunit(jz) = real(jz)
   imask(1:nc0a,jz) = reshape( xctrl%geom%ocean%mask2d(isc:iec, jsc:jec), (/nc0a/) )
end do

area = 1.0           ! Dummy area
vunit = 1.0          ! Dummy vertical unit
imask = 1            ! Mask

nam%default_seed = .true.
nam%model = 'oops'
nam%mask_type = 'none'
nam%new_hdiag = .false.
nam%displ_diag = .false.
nam%new_param = .false.
nam%new_lct = .false.
nam%mask_check = .false.
nam%new_obsop = .false.
nam%check_dirac = .false.
nam%nc3 = 1
nam%dc = 1.0
nam%nl = nl0
nam%nv = nv
nam%nts = nts
nam%timeslot = 0
nam%datadir = '.'
nam%colorlog = .false.
nam%ens1_ne_offset = 0
do il=1,nam%nl
   nam%levs(il) = il
end do
nam%method = 'cor'
nam%strategy = 'common'
nam%diag_interp = 'natural'
nam%flt_type = 'gc99'

geom%nc0a = nc0a
geom%nl0 = nl0
geom%nlev = nl0

!Initialize random number generator
call create_randgen(nam)
print *,"================="
! Initialize coordinates
call model_oops_coord(geom,lon,lat,area,vunit,imask)
print *,"================="
call bpar_alloc(nam,geom,bpar) 
print *,"================="
call compute_grid_mesh(nam,geom)
print *,"================="
call convert_to_ug(xctrl, ug)
print *,"================="
print *,ndata%nsb
!ndata%nam = nam

call namcheck(nam)

call apply_nicas(geom,ndata,ug%fld)
!call apply_localization_from_sqrt(nam,geom,bpar,ndata,ug%fld)
print *,"================="
!call convert_from_ug(xincr, ug)
print *,"================="

end subroutine soca_3d_covar_sqrt_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C) - Adjoint

subroutine soca_3d_covar_sqrt_mult_ad(kx,ky,xincr,xctl,config)
use iso_c_binding
use kinds
use soca_fields

implicit none
integer(c_int), intent(in)    :: kx            !< Zonal grid spacing
integer(c_int), intent(in)    :: ky            !< Meridional grid spacing
real(c_double), intent(inout) :: xctl(kx,ky,2) !< Result
type(soca_field), intent(in)    :: xincr         !< Streamfunction: psi
type(soca_3d_covar_config), intent(in) :: config !< covar config structure

real(kind=kind_real), allocatable :: xout(:,:,:)
integer :: i, j, k, iri, m
real(kind=kind_real) :: zc, zero, one
real(kind=kind_real) :: zgrid(kx), zfour(kx+2), work(ky)



end subroutine soca_3d_covar_sqrt_mult_ad

! ------------------------------------------------------------------------------

end module soca_covariance_mod
