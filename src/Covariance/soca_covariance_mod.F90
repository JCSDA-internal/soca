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
     real(kind=kind_real) :: ocean_alpha_Lx=1.0      !< Zonal       ] Length scale
     real(kind=kind_real) :: ocean_alpha_Ly=1.0      !< Meridional  ] for
     real(kind=kind_real) :: ocean_alpha_Lz=1.0      !< vertical    ] convolution kernel
     real(kind=kind_real) :: ice_Lx=200.0e3      !< Zonal       ] Length scale
     real(kind=kind_real) :: ice_Lz=1.0      !< vertical    ] convolution kernel     
     character(len=800)   :: D_filename  !< Netcdf file containing
                                         !< the diagonal matrix of standard deviation for
     !< all the fields
  end type soca_3d_covar_config

#define LISTED_TYPE soca_3d_covar_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
 type(registry_t) :: soca_3d_cov_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

  !> Setup for the SOCA model's 3d error covariance matrices (B and Q_i)

  !> This routine queries the configuration for the parameters that define the
  !! covariance matrix, and stores the relevant values in the
  !! error covariance structure.

  subroutine soca_3d_covar_setup(c_conf, geom, config, bkg)

    use soca_geom_mod
    use soca_fields    
    use iso_c_binding
    use config_mod
    use kinds
    use fckit_log_module, only : fckit_log
    use type_bump
    use type_nam
    
    implicit none
    
    type(c_ptr),                   intent(in) :: c_conf   !< The configuration
    type(soca_geom),               intent(in) :: geom     !< Geometry
    type(soca_3d_covar_config), intent(inout) :: config   !< The covariance structure
    type(soca_field),              intent(in) :: bkg      !< Background
    
    call soca_init_lengthscale(c_conf, config, bkg)
    
  end subroutine soca_3d_covar_setup

  ! ------------------------------------------------------------------------------

  !> Delete for the SOCA model's 3d error covariance matrices

  subroutine soca_3d_covar_delete(c_key_conf)

    use iso_c_binding

    implicit none
    integer(c_int), intent(inout) :: c_key_conf !< The model covariance structure

    call soca_3d_cov_registry%remove(c_key_conf)
    !call soca_bump_correlation(destruct=.true.)
    
  end subroutine soca_3d_covar_delete

  ! ------------------------------------------------------------------------------
  
  subroutine soca_init_D(geom, bkg, D_p)
    use soca_geom_mod
    use soca_fields

    implicit none
    
    type(soca_geom),            intent(in) :: geom    !< Geometry
    type(soca_field),           intent(in) :: bkg     !< Background field
    type(soca_field), pointer, intent(out) :: D_p     !< Std of backcround error

    logical,                  save :: D_initialized = .false.
    type(soca_field), save, target :: D  !< Std of backcround error

    if (.not.D_initialized) then
       !call create(D,bkg)
       !call zeros(D)       
    end if
    D_p => D
    
  end subroutine soca_init_D

  subroutine soca_init_lengthscale(c_conf, config, bkg)
    use iso_c_binding
    use config_mod
    use kinds
    use soca_fields
    use fms_io_mod,      only : fms_io_init, fms_io_exit

    implicit none

    type(c_ptr),                 intent(in) :: c_conf   !< The configuration
    type(soca_3d_covar_config), intent(inout) :: config   !< Config parameters for D

    type(soca_field) :: bkg
    type(soca_field) :: lensca    
    integer :: inzo, incat
    character(len=800) :: filename

    ! Get configuration
    config%ocean_alpha_Lx  = config_get_real(c_conf,"ocean_alpha_Lx")
    config%ocean_alpha_Lz  = config_get_real(c_conf,"ocean_alpha_Lz")    
    config%ice_Lx  = config_get_real(c_conf,"ice_Lx")
    config%ice_Lz  = config_get_real(c_conf,"ice_Lz")

    ! Setup a copy of bkg to store rh and rv
    call create_copy(lensca, bkg)

    ! Set rh to a large default value
    call ones(lensca)
    call self_mul(lensca, 9999.9d3)
    
    ! Setup horizontal length scale
    ! Ocean
    call bkg%geom%ocean%get_rossby_radius()    

    do inzo = 1,bkg%geom%ocean%nzo
       lensca%tocn(:,:,inzo) = config%ocean_alpha_Lx*&
            &max(200e3, 5.0*bkg%geom%ocean%rossby_radius)
       lensca%socn(:,:,inzo) = lensca%tocn(:,:,inzo)
       lensca%hocn(:,:,inzo) = lensca%tocn(:,:,inzo)
    end do
    lensca%ssh = lensca%tocn(:,:,1)

    ! Sea-ice
    do incat = 1,bkg%geom%ocean%ncat
       lensca%cicen(:,:,incat+1) = config%ice_Lx
       lensca%hicen(:,:,incat) = config%ice_Lx
    end do

    filename="rh.nc"
    call fld2file(lensca, filename)

    ! Set rv to a large default value
    call ones(lensca)
    call self_mul(lensca, 9999.9d3)
    
    ! Setup vertical length scale
    ! Ocean
    do inzo = 1,bkg%geom%ocean%nzo
       lensca%tocn(:,:,inzo) = config%ocean_alpha_Lz*20
       lensca%socn(:,:,inzo) = config%ocean_alpha_Lz*20
       lensca%hocn(:,:,inzo) = config%ocean_alpha_Lz*20
    end do
    lensca%ssh = config%ocean_alpha_Lz*20

    ! Sea-ice
    do incat = 1,bkg%geom%ocean%ncat
       lensca%cicen(:,:,incat+1) = config%ice_Lx
       lensca%hicen(:,:,incat) = config%ice_Lx
    end do

    filename="rv.nc"
    call fld2file(lensca, filename)   
    
  end subroutine soca_init_lengthscale

  
end module soca_covariance_mod
