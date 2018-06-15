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
    use mpi,              only : mpi_comm_world
    
    implicit none
    
    type(c_ptr), intent(in)   :: c_model  !< The configuration
    type(soca_geom), intent(in) :: geom     !< Geometry
    type(soca_3d_covar_config), intent(inout) :: config !< The covariance structure
    real(kind=kind_real) :: corr_length_scale
    type(bump_type), pointer            :: horiz_convol_p
    
    config%sig_sic      = config_get_real(c_model,"sig_sic")
    config%sig_sit      = config_get_real(c_model,"sig_sit")
    config%sig_ssh      = config_get_real(c_model,"sig_ssh")
    config%sig_tocn      = config_get_real(c_model,"sig_tocn")
    config%sig_socn      = config_get_real(c_model,"sig_socn")        

    call soca_bump_correlation(geom, horiz_convol_p)
    
  end subroutine soca_3d_covar_setup

  ! ------------------------------------------------------------------------------

  !> Delete for the SOCA model's 3d error covariance matrices

  subroutine soca_3d_covar_delete(c_key_conf)

    use iso_c_binding

    implicit none
    integer(c_int), intent(inout) :: c_key_conf !< The model covariance structure

    call soca_3d_cov_registry%remove(c_key_conf)
    call soca_bump_correlation(destruct=.true.)
    
  end subroutine soca_3d_covar_delete

  ! ------------------------------------------------------------------------------
  
  subroutine soca_struct2unstruct(dx_struct, geom, dx_unstruct)
    use soca_geom_mod

    implicit none

    real(kind=kind_real),intent(in)                :: dx_struct(:,:)
    type(soca_geom), intent(in)                    :: geom    
    real(kind=kind_real), allocatable, intent(out) :: dx_unstruct(:)

    integer :: isc, iec, jsc, jec, jjj, jz, il, ib, nc0a

    ! Indices for compute domain (no halo)
    isc = geom%ocean%G%isc
    iec = geom%ocean%G%iec
    jsc = geom%ocean%G%jsc
    jec = geom%ocean%G%jec
    
    nc0a = (iec - isc + 1) * (jec - jsc + 1 )
    allocate(dx_unstruct(nc0a))
    dx_unstruct = reshape( dx_struct(isc:iec, jsc:jec), (/nc0a/) )
    
  end subroutine soca_struct2unstruct

  ! ------------------------------------------------------------------------------
  
  subroutine soca_unstruct2struct(dx_struct, geom, dx_unstruct)
    use soca_geom_mod

    implicit none

    real(kind=kind_real),intent(inout)               :: dx_struct(:,:)
    type(soca_geom), intent(in)                      :: geom    
    real(kind=kind_real), allocatable, intent(inout) :: dx_unstruct(:)

    integer :: isc, iec, jsc, jec, jjj, jz, il, ib, nc0a

    ! Indices for compute domain (no halo)
    isc = geom%ocean%G%isc
    iec = geom%ocean%G%iec
    jsc = geom%ocean%G%jsc
    jec = geom%ocean%G%jec
    
    nc0a = (iec - isc + 1) * (jec - jsc + 1 )

    dx_struct(isc:iec, jsc:jec) = reshape(dx_unstruct,(/size(dx_struct(isc:iec, jsc:jec),1),&
                                                       &size(dx_struct(isc:iec, jsc:jec),2)/))

    deallocate(dx_unstruct)
    
  end subroutine soca_unstruct2struct
  
  ! ------------------------------------------------------------------------------

  subroutine soca_3d_covar_C_mult(dx, Cdx, config)
    use iso_c_binding
    use kinds
    use soca_fields
    use type_bump
    use fms_mod,                 only: read_data, write_data, set_domain
    use fms_io_mod,                only : fms_io_init, fms_io_exit
    
    implicit none
    type(soca_field), intent(in)           :: dx
    type(soca_field), intent(inout)        :: Cdx
    type(soca_3d_covar_config), intent(in) :: config !< covariance config structure
    type(bump_type), pointer            :: horiz_convol_p
    real(kind=kind_real), allocatable      :: tmp_incr(:)
    !Grid stuff
    integer :: isc, iec, jsc, jec, jjj, jz, il, ib, nc0a, icat, izo
    character(len=128) :: filename

    call copy(Cdx, dx)
    print *,'Init bump'
    call soca_bump_correlation(dx%geom, horiz_convol_p)
    print *,'Done init bump'

!!$    call fms_io_init()
!!$    filename='test-cov1.nc'
!!$    call write_data( filename, "ssh", dx%ssh, dx%geom%ocean%G%Domain%mpp_domain)
!!$    call fms_io_exit()
!!$    print *,'wrote cov1 to file'
!!$    read(*,*)

    
    print *,'Apply nicas: ssh'
    call soca_struct2unstruct(dx%ssh, dx%geom, tmp_incr)    
    call horiz_convol_p%apply_nicas(tmp_incr)
    call soca_unstruct2struct(Cdx%ssh, dx%geom, tmp_incr)

!!$    call fms_io_init()
!!$    filename='test-cov2.nc'
!!$    call write_data( filename, "ssh", Cdx%ssh, dx%geom%ocean%G%Domain%mpp_domain)
!!$    call fms_io_exit()
!!$    print *,'wrote cov1 to file'
!!$    read(*,*)

    
    do icat = 1, dx%geom%ocean%ncat
       print *,'Apply nicas: aice, hice, category:',icat       
       call soca_struct2unstruct(dx%cicen(:,:,icat+1), dx%geom, tmp_incr)    
       call horiz_convol_p%apply_nicas(tmp_incr)
       call soca_unstruct2struct(Cdx%cicen(:,:,icat+1), dx%geom, tmp_incr)    

       call soca_struct2unstruct(dx%hicen(:,:,icat), dx%geom, tmp_incr)    
       call horiz_convol_p%apply_nicas(tmp_incr)
       call soca_unstruct2struct(Cdx%hicen(:,:,icat), dx%geom, tmp_incr)    
    end do    

    do izo = 1,1!dx%geom%ocean%nzo
       print *,'Apply nicas: tocn, socn, layer:',izo       
       call soca_struct2unstruct(dx%tocn(:,:,izo), dx%geom, tmp_incr)    
       call horiz_convol_p%apply_nicas(tmp_incr)
       call soca_unstruct2struct(Cdx%tocn(:,:,izo), dx%geom, tmp_incr)    

       call soca_struct2unstruct(dx%socn(:,:,izo), dx%geom, tmp_incr)    
       call horiz_convol_p%apply_nicas(tmp_incr)
       call soca_unstruct2struct(Cdx%socn(:,:,izo), dx%geom, tmp_incr)
    end do    

  end subroutine soca_3d_covar_C_mult

  ! ------------------------------------------------------------------------------  
  
  subroutine soca_3d_covar_K_mult_ad(dy, KTdy, traj)
    use iso_c_binding
    use kinds
    use soca_fields
    use soca_balanceop
    use soca_seaice_balanceop
    
    implicit none
    type(soca_field), intent(inout) :: KTdy     !< K^T dx
    type(soca_field), intent(in)    :: dy       !< Increment        
    type(soca_field), intent(in)    :: traj     !< trajectory    

    !Grid stuff
    integer :: isc, iec, jsc, jec, i, j, k
    real(kind=kind_real) :: tb, sb, dt, ds,deta, z, p, h
    real(kind=kind_real), allocatable :: dcn(:), cnb(:)

    ! Indices for compute domain (no halo)
    isc = traj%geom%ocean%G%isc
    iec = traj%geom%ocean%G%iec
    jsc = traj%geom%ocean%G%jsc
    jec = traj%geom%ocean%G%jec

    call copy(KTdy,dy)

    ! Steric height/density balance
    do i = isc, iec
       do j = jsc, jec
          do k = 1, traj%geom%ocean%nzo
             tb=traj%tocn(i,j,k)
             sb=traj%socn(i,j,k)
             if (k.eq.1) then
                z=traj%hocn(i,j,k)
             else
                z=sum(traj%hocn(i,j,1:k-1))+0.5*traj%hocn(i,j,k)
             end if
             h=traj%hocn(i,j,k)
             p=z
             deta=dy%ssh(i,j)
             call steric_ad(deta, dt, ds, tb, sb, p, h)
             KTdy%tocn(i,j,k)=KTdy%tocn(i,j,k)+dt
             KTdy%socn(i,j,k)=KTdy%socn(i,j,k)+ds
          end do
       end do
    end do

    ! T/C balance
    allocate(dcn(traj%geom%ocean%ncat),cnb(traj%geom%ocean%ncat))
    do i = isc, iec
       do j = jsc, jec
          tb=traj%tocn(i,j,1)
          sb=traj%socn(i,j,1)
          cnb=traj%cicen(i,j,2:)
          call tofc_ad (dt, dcn, tb, sb, cnb)
          do k = 1, traj%geom%ocean%ncat
             KTdy%cicen(i,j,k+1)=KTdy%cicen(i,j,k+1)+dcn(k)
          end do
       end do
    end do

    deallocate(dcn, cnb)
    
  end subroutine soca_3d_covar_K_mult_ad

  ! ------------------------------------------------------------------------------  
  
  subroutine soca_3d_covar_K_mult(dx, Kdx, traj)
    use iso_c_binding
    use kinds
    use soca_fields
    use soca_balanceop
    use soca_seaice_balanceop
    
    implicit none
    type(soca_field), intent(inout) :: Kdx      !< K dx
    type(soca_field), intent(in)    :: dx       !< Increment            
    type(soca_field), intent(in)    :: traj     !< trajectory    

    !Grid stuff
    integer :: isc, iec, jsc, jec, i, j, k
    real(kind=kind_real) :: tb, sb, dt, ds,deta, z, p, h
    real(kind=kind_real), allocatable :: dcn(:), cnb(:)

    ! Indices for compute domain (no halo)
    isc = traj%geom%ocean%G%isc
    iec = traj%geom%ocean%G%iec
    jsc = traj%geom%ocean%G%jsc
    jec = traj%geom%ocean%G%jec

    call copy(Kdx,dx)

    ! Steric height/density balance
    do i = isc, iec
       do j = jsc, jec
          do k = 1, traj%geom%ocean%nzo
             tb=traj%tocn(i,j,k)
             sb=traj%socn(i,j,k)
             dt=dx%tocn(i,j,k)
             ds=dx%socn(i,j,k)
             if (k.eq.1) then
                z=traj%hocn(i,j,k)
             else
                z=sum(traj%hocn(i,j,1:k-1))+0.5*traj%hocn(i,j,k)
             end if
             h=traj%hocn(i,j,k)
             p=z
             call steric_tl(deta, dt, ds, tb, sb, p, h)
             Kdx%ssh(i,j)=Kdx%ssh(i,j)+deta    
          end do
       end do
    end do

    ! T/C balance
    allocate(dcn(traj%geom%ocean%ncat),cnb(traj%geom%ocean%ncat))
    do i = isc, iec
       do j = jsc, jec
          tb=traj%tocn(i,j,1)
          sb=traj%socn(i,j,1)
          cnb=traj%cicen(i,j,2:)
          dcn=dx%cicen(i,j,2:)          
          call tofc_tl (dt, dcn, tb, sb, cnb)
          Kdx%tocn(i,j,1)=Kdx%tocn(i,j,1)+dt
       end do
    end do

    deallocate(dcn)
    
  end subroutine soca_3d_covar_K_mult

  ! ------------------------------------------------------------------------------
  
  subroutine soca_3d_covar_D_mult(Ddx, traj, config)
    use iso_c_binding
    use kinds
    use soca_fields
    use soca_interph_mod
    use fms_mod,                 only: read_data, write_data, set_domain
    use fms_io_mod,                only : fms_io_init, fms_io_exit
    
    implicit none
    type(soca_field), intent(inout)        :: Ddx      !< D applied to dx
    type(soca_field), intent(in)           :: traj     !< trajectory   
    type(soca_3d_covar_config), intent(in) :: config   !< covariance config structure
    type(soca_field)    :: sig_cicen                   !< D applied to dx

    !call create_copy(sig_cicen,Ddx)
    !call fms_io_init()
    !call read_data('std.nc', 'cicen', sig_cicen%AOGCM%Ice%part_size, domain=Ddx%geom%ocean%G%Domain%mpp_domain)
    !call fms_io_exit()

    !Ddx%cicen=1.0Ddx%cicen+0.01

    ! A "bit" of a hack!!!
    !sig_cicen%cicen=exp( -((0.15-traj%cicen)/0.1)**2 )
    !sig_cicen%cicen=exp( -((0.15-traj%cicen)    
    
    Ddx%cicen=config%sig_sic*Ddx%cicen
    Ddx%hicen=config%sig_sit*Ddx%hicen
    Ddx%ssh=config%sig_ssh*Ddx%ssh
    Ddx%tocn=config%sig_tocn*Ddx%tocn
    Ddx%socn=config%sig_socn*Ddx%socn

  end subroutine soca_3d_covar_D_mult

  ! ------------------------------------------------------------------------------

  subroutine soca_bump_correlation(geom, horiz_convol_p, destruct)
    use soca_interph_mod
    use soca_geom_mod
    use type_bump
    use type_nam
    use mpi!,             only: mpi_comm_world
    
    implicit none

    type(soca_geom), optional, intent(in)            :: geom
    type(bump_type), optional, pointer, intent(out)  :: horiz_convol_p
    logical, optional                                :: destruct       ! If true: call bump destructor
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
    integer :: ierr
    real :: start, finish
    
    real(kind_real), allocatable :: rh(:,:,:,:)     !< Horizontal support radius for covariance (in m)
    real(kind_real), allocatable :: rv(:,:,:,:)     !< Vertical support radius for

    ! Desructor
    if (present(destruct)) then
       convolh_initialized = .false.
       call horiz_convol%dealloc()
       return
    end if

    ! Constructor
    if (.NOT.convolh_initialized) then

       !--- Initialize geometry to be passed to NICAS
       ! Indices for compute domain (no halo)
       isc = geom%ocean%G%isc
       iec = geom%ocean%G%iec
       jsc = geom%ocean%G%jsc
       jec = geom%ocean%G%jec

       nv = 1!geom%ocean%ncat + 1                 !< Number of variables
       nl0 = 1                                    !< Number of independent levels
       nts = 1                                    !< Number of time slots
       nc0a = (iec - isc + 1) * (jec - jsc + 1 )  !< Total number of grid cells in the compute domain

       allocate( lon(nc0a), lat(nc0a), area(nc0a) )
       allocate( vunit(nc0a,nl0) )
       allocate( imask(nc0a, nl0), lmask(nc0a, nl0) )

       lon = reshape( geom%ocean%lon(isc:iec, jsc:jec), (/nc0a/) )
       lat = reshape( geom%ocean%lat(isc:iec, jsc:jec), (/nc0a/) )        
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
       !horiz_convol%nam%mask_type= "none"
       horiz_convol%nam%mask_check = .true. !.false.
       horiz_convol%nam%draw_type = "random_uniform"
       horiz_convol%nam%nc1 = 100!0

       horiz_convol%nam%ntry = 3
       horiz_convol%nam%nrep =  2
       horiz_convol%nam%nc3 = 10
       
       horiz_convol%nam%dc = 200.0e3
       horiz_convol%nam%nl0r = 1

       horiz_convol%nam%gau_approx = .false.
       horiz_convol%nam%full_var = .false.
       horiz_convol%nam%local_diag = .false.
       horiz_convol%nam%minim_algo= "hooke"
       horiz_convol%nam%lhomh = .false.
       horiz_convol%nam%lhomv = .false.
       horiz_convol%nam%rvflt =0.0
       horiz_convol%nam%lsqrt = .false.!.true.
       horiz_convol%nam%resol = 10.0 !1000.0!10.0
       horiz_convol%nam%nicas_interp = "bilin"
       horiz_convol%nam%network = .true.!.false.
       horiz_convol%nam%mpicom = 1
       !horiz_convol%nam%advmode = 0
       !horiz_convol%nam%nldwh = 0
       !horiz_convol%nam%nldwv = 0
       !horiz_convol%nam%diag_rhflt = 0.0
       horiz_convol%nam%diag_interp = "bilin"
       horiz_convol%nam%grid_output = .false.

       allocate(rh(nc0a,nl0,nv,nts))
       allocate(rv(nc0a,nl0,nv,nts))
       rh=200.0e3
       rv=1.0

       call cpu_time(start)
       print *,"Time start = ",start," seconds."       
       call horiz_convol%bump_setup_online(mpi_comm_world,nc0a,nl0,nv,nts,lon,lat,area,vunit,lmask,rh=rh,rv=rv)
       call cpu_time(finish)
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       print *,"Time = ",finish-start," seconds."
       convolh_initialized = .true.
       deallocate( lon, lat, area, vunit, imask, lmask )
       deallocate(rh,rv)       
    end if
    horiz_convol_p => horiz_convol

    
  end subroutine soca_bump_correlation

end module soca_covariance_mod
