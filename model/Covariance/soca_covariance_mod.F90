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
     real(kind=kind_real) :: Lx=1.0      !< Zonal       ] Length scale
     real(kind=kind_real) :: Ly=1.0      !< Meridional  ] for
     real(kind=kind_real) :: Lz=1.0      !< vertical    ] convolution kernel
     real(kind=kind_real) :: sig_sic     !<   
     real(kind=kind_real) :: sig_sit     !< Temporary hack 
     real(kind=kind_real) :: sig_ssh     !<
     real(kind=kind_real) :: sig_tocn     !<
     real(kind=kind_real) :: sig_socn     !<          
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
    config%sig_tocn      = config_get_real(c_model,"sig_tocn")
    config%sig_socn      = config_get_real(c_model,"sig_socn")        

    !call initialize_convolh(geom, horiz_convol_p)
    
  end subroutine soca_3d_covar_setup

  ! ------------------------------------------------------------------------------

  !> Delete for the SOCA model's 3d error covariance matrices

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
    real(kind=kind_real), allocatable      :: tmp_incr(:,:),tmp_incr3d(:,:,:)

    integer :: k,l,m,iter

    call copy(sqrtCdx, dx)
    
    ! Temporary hack while waiting for new interfaces for NICAS
    ! sqrtCdx=C.dx
    allocate(tmp_incr(size(dx%ssh,1),size(dx%ssh,2)))
    !call zeros(sqrtCdx)
!!$!    sqrtCdx%ssh=dx%ssh    
    do iter = 1, 1
       tmp_incr=sqrtCdx%ssh*dx%geom%ocean%mask2d       
       do k = 2, size(dx%ssh,1)-1
          do l = 2, size(dx%ssh,2)-1
             sqrtCdx%ssh(k,l)=(tmp_incr(k+1,l-1)+tmp_incr(k,l-1)+tmp_incr(k-1,l-1)+&
                              &tmp_incr(k+1,l)+tmp_incr(k,l)+tmp_incr(k-1,l)+&
                              &tmp_incr(k+1,l+1)+tmp_incr(k,l+1)+tmp_incr(k-1,l+1))/9.0
          end do
       end do
    end do

    allocate(tmp_incr3d(size(dx%ssh,1),size(dx%ssh,2),dx%geom%ocean%ncat))
    sqrtCdx%hicen=dx%hicen    
    do iter = 1, 1
       tmp_incr3d=sqrtCdx%hicen       
       do k = 2, size(dx%ssh,1)-1
          do l = 2, size(dx%ssh,2)-1
             do m = 1, dx%geom%ocean%ncat
                sqrtCdx%hicen(k,l,m)=(tmp_incr3d(k+1,l-1,m)+tmp_incr3d(k,l-1,m)+tmp_incr3d(k-1,l-1,m)+&
                              &tmp_incr3d(k+1,l,m)+tmp_incr3d(k,l,m)+tmp_incr3d(k-1,l,m)+&
                              &tmp_incr3d(k+1,l+1,m)+tmp_incr3d(k,l+1,m)+tmp_incr3d(k-1,l+1,m))/9.0
             end do
          end do
       end do
    end do

    sqrtCdx%cicen=dx%cicen
    do iter = 1, 2
       tmp_incr3d=sqrtCdx%cicen(:,:,2:dx%geom%ocean%ncat+1)
       do k = 2, size(dx%ssh,1)-1
          do l = 2, size(dx%ssh,2)-1
             do m = 1, dx%geom%ocean%ncat
                sqrtCdx%cicen(k,l,m+1)=(tmp_incr3d(k+1,l-1,m)+tmp_incr3d(k,l-1,m)+tmp_incr3d(k-1,l-1,m)+&
                              &tmp_incr3d(k+1,l,m)+tmp_incr3d(k,l,m)+tmp_incr3d(k-1,l,m)+&
                              &tmp_incr3d(k+1,l+1,m)+tmp_incr3d(k,l+1,m)+tmp_incr3d(k-1,l+1,m))/9.0
             end do
          end do
       end do
    end do

    deallocate(tmp_incr3d)
    allocate(tmp_incr3d(size(dx%ssh,1),size(dx%ssh,2),dx%geom%ocean%nzo))
    sqrtCdx%tocn=dx%tocn    
    do iter = 1, 1
       tmp_incr3d=sqrtCdx%tocn
       do k = 2, size(dx%ssh,1)-1
          do l = 2, size(dx%ssh,2)-1
             do m = 1, dx%geom%ocean%nzo
                sqrtCdx%tocn(k,l,m)=(tmp_incr3d(k+1,l-1,m)+tmp_incr3d(k,l-1,m)+tmp_incr3d(k-1,l-1,m)+&
                              &tmp_incr3d(k+1,l,m)+tmp_incr3d(k,l,m)+tmp_incr3d(k-1,l,m)+&
                              &tmp_incr3d(k+1,l+1,m)+tmp_incr3d(k,l+1,m)+tmp_incr3d(k-1,l+1,m))/9.0
             end do
          end do
       end do
       tmp_incr3d=sqrtCdx%tocn
       do k = 2, size(dx%ssh,1)-1
          do l = 2, size(dx%ssh,2)-1
             do m = 2, dx%geom%ocean%nzo-1
                sqrtCdx%tocn(k,l,m)=(tmp_incr3d(k,l,m-1)+tmp_incr3d(k,l,m)+tmp_incr3d(k,l,m+1))/3.0
             end do
          end do
       end do
    end do

!!$    sqrtCdx%hicen=dx%hicen
!!$    do iter = 1, 1
!!$       tmp_incr3d=sqrtCdx%hicen(:,:,1:dx%geom%ocean%ncat)
!!$       do k = 2, size(dx%ssh,1)-1
!!$          do l = 2, size(dx%ssh,2)-1
!!$             do m = 1, dx%geom%ocean%ncat
!!$                sqrtCdx%hicen(k,l,m)=(tmp_incr3d(k+1,l-1,m)+tmp_incr3d(k,l-1,m)+tmp_incr3d(k-1,l-1,m)+&
!!$                              &tmp_incr3d(k+1,l,m)+tmp_incr3d(k,l,m)+tmp_incr3d(k-1,l,m)+&
!!$                              &tmp_incr3d(k+1,l+1,m)+tmp_incr3d(k,l+1,m)+tmp_incr3d(k-1,l+1,m))/9.0
!!$             end do
!!$          end do
!!$       end do
!!$    end do        
    
    deallocate(tmp_incr, tmp_incr3d)

  end subroutine soca_3d_covar_sqrt_mult

  ! ------------------------------------------------------------------------------

  !> Multiply by sqrt(C) - Adjoint

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
    deallocate(tmp_incr)

  end subroutine soca_3d_covar_sqrt_mult_ad

  ! ------------------------------------------------------------------------------
  
  subroutine soca_3d_covar_D_mult(Ddx, config)
    use iso_c_binding
    use kinds
    use soca_fields
    use soca_interph_mod
    use fms_mod,                 only: read_data, write_data, set_domain
    use fms_io_mod,                only : fms_io_init, fms_io_exit
    
    implicit none
    type(soca_field), intent(inout)        :: Ddx             !< D applied to dx
    type(soca_3d_covar_config), intent(in) :: config          !< covariance config structure
    type(soca_field)    :: sig_cicen             !< D applied to dx

    call create_copy(sig_cicen,Ddx)
    call fms_io_init()
    call read_data('std.nc', 'cicen', sig_cicen%AOGCM%Ice%part_size, domain=Ddx%geom%ocean%G%Domain%mpp_domain)
    call fms_io_exit()

    !Ddx%cicen=1.0Ddx%cicen+0.01

    ! A "bit" of a hack!!!
    !sig_cicen%cicen=exp( -((0.15-sig_cicen%cicen)/0.1)**2 )
    
    !Ddx%cicen=config%sig_sic*sig_cicen%cicen*Ddx%cicen
    Ddx%hicen=config%sig_sit*Ddx%hicen
    Ddx%ssh=config%sig_ssh*Ddx%ssh
    Ddx%tocn=config%sig_tocn*Ddx%tocn
    Ddx%socn=config%sig_socn*Ddx%socn

  end subroutine soca_3d_covar_D_mult

  ! ------------------------------------------------------------------------------

  subroutine initialize_convolh(geom, horiz_convol_p)
    !use ufo_locs_mod  
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
