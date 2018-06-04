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
    use type_bump
    use type_nam
    use tools_const, only: pi,req,deg2rad,rad2deg
    use mpi,             only: mpi_comm_world
    
    implicit none
    type(c_ptr), intent(in)   :: c_model  !< The configuration
    type(soca_geom), intent(in) :: geom     !< Geometry
    type(soca_3d_covar_config), intent(inout) :: config !< The covariance structure
    real(kind=kind_real) :: corr_length_scale
    !type(soca_hinterp), pointer :: horiz_convol_p
    
    config%sig_sic      = config_get_real(c_model,"sig_sic")
    config%sig_sit      = config_get_real(c_model,"sig_sit")
    config%sig_ssh      = config_get_real(c_model,"sig_ssh")
    config%sig_tocn      = config_get_real(c_model,"sig_tocn")
    config%sig_socn      = config_get_real(c_model,"sig_socn")        

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
    tmp_incr=0.0
!!$!    sqrtCdx%ssh=dx%ssh    
    do iter = 1, 3
       tmp_incr=tmp_incr+sqrtCdx%ssh*dx%geom%ocean%mask2d       
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
    do iter = 1, 3
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
    do iter = 1, 3
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
    do iter = 1, 3
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
    !use soca_interph_mod
    use type_bump
    
    implicit none
    type(soca_field), intent(in)           :: dx
    type(soca_field), intent(inout)        :: sqrtCTdx
    type(soca_3d_covar_config), intent(in) :: config !< covariance config structure
    type(bump_type), pointer            :: horiz_convol_p
    real(kind=kind_real), allocatable      :: tmp_incr(:)
    real(kind=kind_real)      :: crap
!!$    
!!$    call copy(sqrtCTdx, dx)
!!$    call initialize_convolh(dx%geom, horiz_convol_p)
!!$    allocate(tmp_incr(size(sqrtCTdx%ssh,1)*size(sqrtCTdx%ssh,2)))
!!$    tmp_incr=reshape(dx%ssh,(/size(dx%ssh,1)*size(dx%ssh,2)/))
!!$    !call horiz_convol_p%interpad_apply(sqrtCTdx%ssh,tmp_incr)
!!$    deallocate(tmp_incr)

  end subroutine soca_3d_covar_sqrt_mult_ad

  ! ------------------------------------------------------------------------------

  !> Multiply by sqrt(C) - Adjoint

  subroutine soca_3d_covar_mult(dx, Cdx, config)
    use iso_c_binding
    use kinds
    use soca_fields
    !use soca_interph_mod
    use type_bump
    
    implicit none
    type(soca_field), intent(in)           :: dx
    type(soca_field), intent(inout)        :: Cdx
    type(soca_3d_covar_config), intent(in) :: config !< covariance config structure
    type(bump_type), pointer            :: horiz_convol_p
    real(kind=kind_real), allocatable      :: tmp_incr(:)
    !Grid stuff
    integer :: isc, iec, jsc, jec, jjj, jz, il, ib

    !--- Initialize geometry to be passed to NICAS
    ! Indices for compute domain (no halo)
    !isc = dx%geom%ocean%G%isc
    !iec = dx%geom%ocean%G%iec
    !jsc = dx%geom%ocean%G%jsc
    !jec = dx%geom%ocean%G%jec
    
    call copy(Cdx, dx)
    print *,'init bump'
    call initialize_convolh(dx%geom, horiz_convol_p)
    print *,'done init bump'    
    allocate(tmp_incr(size(Cdx%ssh,1)*size(Cdx%ssh,2)))
    tmp_incr=reshape(dx%ssh,(/size(dx%ssh,1)*size(dx%ssh,2)/))

    print *,'apply nicas'
    
    call horiz_convol_p%apply_nicas(tmp_incr)
    print *,'apply nicas done'
    Cdx%ssh = reshape(tmp_incr,(/size(dx%ssh,1),size(dx%ssh,2)/))
    
    deallocate(tmp_incr)

  end subroutine soca_3d_covar_mult

  ! ------------------------------------------------------------------------------  
  
  subroutine soca_3d_covar_K_mult_ad(dy, KTdy, traj)
    use iso_c_binding
    use kinds
    use soca_fields
    use soca_balanceop
    
    implicit none
    type(soca_field), intent(inout) :: KTdy     !< K^T dx
    type(soca_field), intent(in)    :: dy       !< Increment        
    type(soca_field), intent(in)    :: traj     !< trajectory    

    !Grid stuff
    integer :: isc, iec, jsc, jec, i, j, k
    real(kind=kind_real) :: t0, s0, dt, ds,deta, z, p, h

    ! Indices for compute domain (no halo)
    isc = traj%geom%ocean%G%isc
    iec = traj%geom%ocean%G%iec
    jsc = traj%geom%ocean%G%jsc
    jec = traj%geom%ocean%G%jec

    call copy(KTdy,dy)

    do i = isc, iec
       do j = jsc, jec
          do k = 1, traj%geom%ocean%nzo
             t0=traj%tocn(i,j,k)
             s0=traj%socn(i,j,k)
             if (k.eq.1) then
                z=traj%hocn(i,j,k)
             else
                z=sum(traj%hocn(i,j,1:k-1))+0.5*traj%hocn(i,j,k)
             end if
             h=traj%hocn(i,j,k)
             p=z
             deta=dy%ssh(i,j)
             call steric_ad(deta, dt, ds, t0, s0, p, h)
             KTdy%tocn(i,j,k)=KTdy%tocn(i,j,k)+dt
             KTdy%socn(i,j,k)=KTdy%socn(i,j,k)+ds
          end do
       end do
    end do

  end subroutine soca_3d_covar_K_mult_ad

  ! ------------------------------------------------------------------------------  
  
  subroutine soca_3d_covar_K_mult(dx, Kdx, traj)
    use iso_c_binding
    use kinds
    use soca_fields
    use soca_balanceop
    
    implicit none
    type(soca_field), intent(inout) :: Kdx      !< K dx
    type(soca_field), intent(in)    :: dx       !< Increment            
    type(soca_field), intent(in)    :: traj     !< trajectory    

    !Grid stuff
    integer :: isc, iec, jsc, jec, i, j, k
    real(kind=kind_real) :: t0, s0, dt, ds,deta, z, p, h

    ! Indices for compute domain (no halo)
    isc = traj%geom%ocean%G%isc
    iec = traj%geom%ocean%G%iec
    jsc = traj%geom%ocean%G%jsc
    jec = traj%geom%ocean%G%jec

    call copy(Kdx,dx)

    do i = isc, iec
       do j = jsc, jec
          do k = 1, traj%geom%ocean%nzo
             t0=traj%tocn(i,j,k)
             s0=traj%socn(i,j,k)
             dt=dx%tocn(i,j,k)
             ds=dx%socn(i,j,k)
             if (k.eq.1) then
                z=traj%hocn(i,j,k)
             else
                z=sum(traj%hocn(i,j,1:k-1))+0.5*traj%hocn(i,j,k)
             end if
             h=traj%hocn(i,j,k)
             p=z
             call steric_tl(deta, dt, ds, t0, s0, p, h)
             Kdx%ssh(i,j)=Kdx%ssh(i,j)+deta    
          end do
       end do
    end do

  end subroutine soca_3d_covar_K_mult

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
    Ddx%tocn=config%sig_tocn*Ddx%tocn
    Ddx%socn=config%sig_socn*Ddx%socn

  end subroutine soca_3d_covar_D_mult

  ! ------------------------------------------------------------------------------

  subroutine initialize_convolh(geom, horiz_convol_p)
    !use ufo_locs_mod  
    use soca_interph_mod
    use soca_geom_mod
    !use soca_interph_mod
    use type_bump
    use type_nam
    use tools_const, only: pi,req,deg2rad,rad2deg
    use mpi,             only: mpi_comm_world
    
    implicit none

    type(soca_geom), intent(in)            :: geom
    type(bump_type), pointer, intent(out)  :: horiz_convol_p

    logical, save                    :: convolh_initialized = .false.
    type(bump_type), save, target    :: horiz_convol

    !Grid stuff
    integer :: isc, iec, jsc, jec, jjj, jz, il, ib
    character(len=1024) :: subr = 'model_write'

    !bump stuff
    integer :: nc0a, nl0, nv, nts
    real(kind=kind_real), allocatable :: lon(:), lat(:), area(:), vunit(:,:)
    logical, allocatable :: lmask(:,:)
    integer, allocatable :: imask(:,:)    
    type(nam_type) :: nam

    real(kind_real), allocatable :: rh(:,:,:,:)     !< Horizontal support radius for covariance (in m)
    real(kind_real), allocatable :: rv(:,:,:,:)     !< Vertical support radius for
    
    if (.NOT.convolh_initialized) then

       !--- Initialize geometry to be passed to NICAS
       ! Indices for compute domain (no halo)
       isc = geom%ocean%G%isc
       iec = geom%ocean%G%iec
       jsc = geom%ocean%G%jsc
       jec = geom%ocean%G%jec

       nv = geom%ocean%ncat + 1                   !< Number of variables
       nl0 = 1                                    !< Number of independent levels
       nts = 1                                    !< Number of time slots
       nc0a = (iec - isc + 1) * (jec - jsc + 1 )  !< Total number of grid cells in the compute domain

       allocate( lon(nc0a), lat(nc0a), area(nc0a) )
       allocate( vunit(nc0a,nl0) )
       allocate( imask(nc0a, nl0), lmask(nc0a, nl0) )
       lon = deg2rad*reshape( geom%ocean%lon(isc:iec, jsc:jec), (/nc0a/) )
       lat = deg2rad*reshape( geom%ocean%lat(isc:iec, jsc:jec), (/nc0a/) ) 
       area = reshape( geom%ocean%cell_area(isc:iec, jsc:jec), (/nc0a/) )
       do jz = 1, nl0       
          vunit(:,jz) = real(jz)
          imask(1:nc0a,jz) = reshape( geom%ocean%mask2d(isc:iec, jsc:jec), (/nc0a/) )
       end do
       vunit = 1.0                      !< Dummy vertical unit

       lmask = .false.
       where (imask.eq.1)
          lmask=.true.
       end where

       print *,'b mat setup ----------------------------'

       horiz_convol%nam%default_seed = .true.
       horiz_convol%nam%new_hdiag = .false.
       horiz_convol%nam%new_param = .true.
       horiz_convol%nam%check_adjoints = .false.
       horiz_convol%nam%check_pos_def = .false.
       horiz_convol%nam%check_sqrt = .false.
       horiz_convol%nam%check_randomization = .false.
       horiz_convol%nam%check_consistency = .false.
       horiz_convol%nam%check_optimality = .false.
       horiz_convol%nam%new_lct = .false.
       horiz_convol%nam%new_obsop = .false.

       horiz_convol%nam%prefix = "soca"
       horiz_convol%nam%method = "cor"
       horiz_convol%nam%strategy = "specific_univariate"
       horiz_convol%nam%sam_write= .true.
       horiz_convol%nam%sam_read= .false.
       horiz_convol%nam%mask_type= "none"
       horiz_convol%nam%mask_check = .false.
       horiz_convol%nam%draw_type = "random_uniform"
       horiz_convol%nam%nc1 = 400
       horiz_convol%nam%ntry = 3
       horiz_convol%nam%nrep =  2
       horiz_convol%nam%nc3 = 10
       horiz_convol%nam%dc = 1000.0e3
       horiz_convol%nam%nl0r = 1
       horiz_convol%nam%ne = 4
       horiz_convol%nam%gau_approx = .false.
       horiz_convol%nam%full_var = .false.
       horiz_convol%nam%local_diag = .false.
       horiz_convol%nam%minim_algo= "hooke"
       horiz_convol%nam%lhomh = .false.
       horiz_convol%nam%lhomv = .false.
       horiz_convol%nam%rvflt =0.0
       horiz_convol%nam%lsqrt = .true.
       horiz_convol%nam%resol =10.0
       horiz_convol%nam%nicas_interp = "bilin"
       horiz_convol%nam%network = .false.
       horiz_convol%nam%mpicom = 2
       horiz_convol%nam%advmode = 0
       horiz_convol%nam%nldwh = 0
       horiz_convol%nam%nldwv = 0
       horiz_convol%nam%diag_rhflt = 0.0
       horiz_convol%nam%diag_interp = "bilin"
       horiz_convol%nam%grid_output = .false.

       allocate(rh(nc0a,nl0,nv,nts))
       allocate(rv(nc0a,nl0,nv,nts))
       rh=500.0e3
       rv=1.0
       call horiz_convol%bump_setup_online(mpi_comm_world,nc0a,nl0,nv,nts,lon,lat,area,vunit,lmask,rh=rh,rv=rv)
       print *,'b mat setup done ----------------------------'
       convolh_initialized = .true.
       deallocate( lon, lat, area, vunit, imask, lmask )       
    end if
    print *,'associate pointer' 
    horiz_convol_p => horiz_convol
    deallocate(rh,rv)
    
  end subroutine initialize_convolh

  ! ------------------------------------------------------------------------------

end module soca_covariance_mod
