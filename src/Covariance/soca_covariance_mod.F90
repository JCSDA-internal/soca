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
    
    type(bump_type),  pointer :: horiz_convol_p !< Convolution op from bump
    type(soca_field), pointer :: D_p            !< Std of background error estimated
                                                !< from background
    
    !< Initialize bump
    call soca_bump_correlation(geom, horiz_convol_p, c_conf)
    
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

  subroutine soca_3d_covar_C_mult(dx, config)
    use iso_c_binding
    use kinds
    use soca_fields
    use type_bump
    
    implicit none

    type(soca_field),        intent(inout) :: dx     !< Input: Increment
                                                     !< Output: C dx
    type(soca_3d_covar_config), intent(in) :: config !< covariance config structure

    type(bump_type), pointer          :: horiz_convol_p
    integer :: icat, izo
    
    ! Initialize BUMP and Associate horiz_convol_p 
    call soca_bump_correlation(dx%geom, horiz_convol_p)

    ! Apply convolution to fields
    print *,'Apply nicas: ssh'
    call soca_2d_convol(dx%ssh, horiz_convol_p, dx%geom)
    
    do icat = 1, dx%geom%ocean%ncat
       print *,'Apply nicas: aice, hice, category:',icat
       call soca_2d_convol(dx%cicen(:,:,icat+1), horiz_convol_p, dx%geom)
       call soca_2d_convol(dx%hicen(:,:,icat), horiz_convol_p, dx%geom)       
    end do    

    do izo = 1,dx%geom%ocean%nzo
       print *,'Apply nicas: tocn, socn, layer:',izo
       call soca_2d_convol(dx%tocn(:,:,izo), horiz_convol_p, dx%geom)
       call soca_2d_convol(dx%socn(:,:,izo), horiz_convol_p, dx%geom)       
    end do    

  end subroutine soca_3d_covar_C_mult
  
  ! ------------------------------------------------------------------------------

  subroutine soca_bump_correlation(geom, horiz_convol_p, c_conf, destruct)
    use soca_geom_mod
    use type_bump
    use type_nam
    use mpi!,             only: mpi_comm_world
    use iso_c_binding
    use oobump_mod, only: bump_read_conf
    use fckit_mpi_module, only: fckit_mpi_comm
    
    implicit none

    type(soca_geom),           optional, intent(in) :: geom
    type(bump_type), optional, pointer, intent(out) :: horiz_convol_p
    type(c_ptr),               optional, intent(in) :: c_conf         !< Handle to configuration    
    logical, optional                                :: destruct       ! If true: call bump destructor
    logical, save                    :: convolh_initialized = .false.
    type(bump_type), save, target    :: horiz_convol

    !Grid stuff
    integer :: isc, iec, jsc, jec, jjj, jz, il, ib
    character(len=1024) :: subr = 'model_write'

    !bump stuff
    integer :: nc0a, nl0, nv, nts
    real(kind=kind_real), allocatable :: lon(:), lat(:), area(:), vunit(:,:)
    real(kind=kind_real), allocatable :: rosrad(:)    
    logical, allocatable :: lmask(:,:)
    integer, allocatable :: imask(:,:)    
    type(nam_type) :: nam
    integer :: ierr
    real :: start, finish
    
    real(kind_real), allocatable :: rh(:,:,:,:)     !< Horizontal support radius for covariance (in m)
    real(kind_real), allocatable :: rv(:,:,:,:)     !< Vertical support radius for
    type(fckit_mpi_comm) :: f_comm

    f_comm = fckit_mpi_comm()

    ! Desructor
    if (present(destruct)) then
       convolh_initialized = .false.
       call horiz_convol%dealloc()
       return
    end if

    ! Constructor
    if (.NOT.convolh_initialized) then
       if (.not.(present(c_conf))) then
          call abor1_ftn ("soca_covariance_mod:soca_bump_correlation error, need to specify configuration")
       end if
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

       allocate( lon(nc0a), lat(nc0a), area(nc0a), rosrad(nc0a) )
       allocate( vunit(nc0a,nl0) )
       allocate( imask(nc0a, nl0), lmask(nc0a, nl0) )
       
       lon = reshape( geom%ocean%lon(isc:iec, jsc:jec), (/nc0a/) )
       lat = reshape( geom%ocean%lat(isc:iec, jsc:jec), (/nc0a/) )        
       area = reshape( geom%ocean%cell_area(isc:iec, jsc:jec), (/nc0a/) )
       rosrad = reshape( geom%ocean%rossby_radius(isc:iec, jsc:jec), (/nc0a/) )

       do jz = 1, nl0       
          vunit(:,jz) = real(jz)
          imask(1:nc0a,jz) = reshape( geom%ocean%mask2d(isc:iec, jsc:jec), (/nc0a/) )
       end do
       vunit = 1.0                      !< Dummy vertical unit

       lmask = .false.
       where (imask.eq.1)
          lmask=.true.
       end where

       allocate(rh(nc0a,nl0,nv,nts))
       allocate(rv(nc0a,nl0,nv,nts))

       do jjj=1,nc0a
          rh(jjj,1,1,1)=10.0*rosrad(jjj)
       end do
       where (rh<500e3)
          rh=500e3
       end where
       rv=1.0

       call cpu_time(start)
       print *,"Time start = ",start," seconds."

       ! Compute convolution weight
       call horiz_convol%nam%init() 
       call bump_read_conf(c_conf, horiz_convol)       
       call horiz_convol%setup_online(f_comm%communicator(),nc0a,nl0,nv,nts,lon,lat,area,vunit,lmask)
       call horiz_convol%set_parameter('cor_rh',rh)       
       call horiz_convol%run_drivers()
       call cpu_time(finish)
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       print *,"Time = ",finish-start," seconds."
       convolh_initialized = .true.
       deallocate( lon, lat, area, vunit, imask, lmask )
       deallocate(rh,rv)       
    end if
    horiz_convol_p => horiz_convol
    
  end subroutine soca_bump_correlation

  ! ------------------------------------------------------------------------------

  subroutine soca_2d_convol(dx, horiz_convol_p, geom)

    use soca_geom_mod
    use kinds
    use type_bump
    
    implicit none
    real(kind=kind_real),     intent(inout) :: dx(:,:)
    type(bump_type),             intent(in) :: horiz_convol_p    
    type(soca_geom),             intent(in) :: geom        

    real(kind=kind_real), allocatable :: tmp_incr(:,:,:,:)

    ! Apply 2D convolution
    call soca_struct2unstruct(dx(:,:), geom, tmp_incr)
    call horiz_convol_p%apply_nicas(tmp_incr)
    call soca_unstruct2struct(dx(:,:), geom, tmp_incr)
    
  end subroutine soca_2d_convol

  ! ------------------------------------------------------------------------------
  
  subroutine soca_struct2unstruct(dx_struct, geom, dx_unstruct)
    use soca_geom_mod

    implicit none

    real(kind=kind_real),intent(in)                :: dx_struct(:,:)
    type(soca_geom), intent(in)                    :: geom    
    real(kind=kind_real), allocatable, intent(out) :: dx_unstruct(:,:,:,:)

    integer :: isc, iec, jsc, jec, jjj, jz, il, ib, nc0a

    ! Indices for compute domain (no halo)
    isc = geom%ocean%G%isc
    iec = geom%ocean%G%iec
    jsc = geom%ocean%G%jsc
    jec = geom%ocean%G%jec
    
    nc0a = (iec - isc + 1) * (jec - jsc + 1 )
    allocate(dx_unstruct(nc0a,1,1,1))
    dx_unstruct = reshape( dx_struct(isc:iec, jsc:jec), (/nc0a,1,1,1/) )
    
  end subroutine soca_struct2unstruct

  ! ------------------------------------------------------------------------------
  
  subroutine soca_unstruct2struct(dx_struct, geom, dx_unstruct)
    use soca_geom_mod

    implicit none

    real(kind=kind_real),intent(inout)               :: dx_struct(:,:)
    type(soca_geom), intent(in)                      :: geom    
    real(kind=kind_real), allocatable, intent(inout) :: dx_unstruct(:,:,:,:)

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
  
end module soca_covariance_mod
