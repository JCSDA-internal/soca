!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bumpcorr2d_mod
  use config_mod
  use fckit_mpi_module, only: fckit_mpi_comm  
  use iso_c_binding  
  use kinds
  use oobump_mod, only: bump_read_conf
  use soca_fields
  use soca_geom_mod_c
  use soca_model_geom_type, only : geom_get_domain_indices
  use type_bump
  use type_nam
  
  implicit none
  private

  type, public :: soca_bumpcorr2d_type
     type(bump_type)      :: horizconv     !< bump corr object
     real(kind=kind_real) :: l0     
   contains
     procedure :: initialize => soca_corr_init
     procedure :: apply => soca_corr_apply     
     procedure :: applyad => soca_corrad_apply
     procedure :: finalize => soca_corr_exit
  end type soca_bumpcorr2d_type

contains

  ! ------------------------------------------------------------------------------
  subroutine soca_corr_init(self, geom, c_conf, prefix)
    class(soca_bumpcorr2d_type), intent(inout) :: self   
    type(soca_geom),                intent(in) :: geom   !< Geometry
    type(c_ptr),                    intent(in) :: c_conf    
    character(len=1024),            intent(in) :: prefix !< File prefix for bump
    
    !Grid stuff
    integer :: isc, iec, jsc, jec, jjj, jz, il, ib

    !bump stuff
    integer :: nc0a, nl0, nv, nts
    real(kind=kind_real), allocatable :: lon(:), lat(:), area(:), vunit(:,:)
    real(kind=kind_real), allocatable :: rosrad(:)    
    logical, allocatable :: lmask(:,:)
    integer, allocatable :: imask(:,:)    
    
    real(kind_real), allocatable :: rh(:,:,:,:)     !< Horizontal support radius for covariance (in m)
    real(kind_real), allocatable :: rv(:,:,:,:)     !< Vertical support radius 

    !--- Initialize geometry to be passed to NICAS
    call geom_get_domain_indices(geom%ocean, "compute", isc, iec, jsc, jec)

    nv = 1                                     !< Number of variables
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

    ! Setup land or ice mask
    jz = 1
    imask(1:nc0a,jz) = int(reshape( geom%ocean%mask2d(isc:iec, jsc:jec), (/nc0a/)))

    lmask = .false.
    where (imask.eq.1)
       lmask=.true.
    end where

    ! No vertical convolution, set dummy vertical unit    
    vunit = 1.0d0

    ! Setup horizontal decorrelation length scales
    allocate(rh(nc0a,nl0,nv,nts))
    allocate(rv(nc0a,nl0,nv,nts))
    do jjj=1,nc0a
       rh(jjj,1,1,1)=self%l0 + rosrad(jjj)
    end do
    rv=1.0 ! Vertical scales not used, set to something

    ! Initialize bump namelist/parameters
    call self%horizconv%nam%init()
    call bump_read_conf(c_conf, self%horizconv)
    self%horizconv%nam%prefix = prefix
    print *,'prefix=',self%horizconv%nam%prefix

    ! Compute convolution weight    
    call self%horizconv%setup_online(nc0a,nl0,nv,nts,lon,lat,area,vunit,lmask)
    call self%horizconv%set_parameter('cor_rh',rh)    
    call self%horizconv%run_drivers()

    ! Clean up
    deallocate(lon, lat, area, vunit, imask, lmask, rosrad)
    deallocate(rh,rv)     
    
  end subroutine soca_corr_init

  !--------------------------------------------  
  subroutine soca_corr_apply(self, dx, geom)
    class(soca_bumpcorr2d_type), intent(inout) :: self    
    real(kind=kind_real),        intent(inout) :: dx(:,:)
    type(soca_geom),                intent(in) :: geom        

    real(kind=kind_real), allocatable :: tmp_incr(:,:,:,:)
    real(kind=kind_real), allocatable :: pcv(:)

    integer :: nn

    ! Allocate unstructured tmp_increment and make copy of dx
    call geom%ocean%struct2unstruct(dx(:,:), tmp_incr)

    ! Get control variable size
    call self%horizconv%get_cv_size(nn)
    allocate(pcv(nn))
    pcv = tmp_incr(:,1,1,1)
    print *,'nn=',nn
    
    ! Apply C^1/2
    call self%horizconv%apply_nicas_sqrt(pcv, tmp_incr)

    ! Back to structured grid
    call geom%ocean%unstruct2struct(dx(:,:), tmp_incr)

    ! Clean up
    deallocate(pcv)
    if (allocated(tmp_incr)) deallocate(tmp_incr)

  end subroutine soca_corr_apply

  ! ------------------------------------------------------------------------------

  subroutine soca_corrad_apply(self, dx, geom)
    class(soca_bumpcorr2d_type),    intent(inout) :: self    
    real(kind=kind_real),           intent(inout) :: dx(:,:)
    type(soca_geom),                   intent(in) :: geom        

    real(kind=kind_real), allocatable :: tmp_incr(:,:,:,:)
    real(kind=kind_real), allocatable :: pcv(:)

    integer :: nn

    ! Allocate unstructured tmp_increment and make copy of dx
    call geom%ocean%struct2unstruct(dx(:,:), tmp_incr)

    ! Get control variable size
    call self%horizconv%get_cv_size(nn)
    allocate(pcv(nn))
    pcv = tmp_incr(:,1,1,1)

    ! Apply C^1/2
    call self%horizconv%apply_nicas_sqrt_ad(tmp_incr, pcv)

    ! Back to structured grid
    call geom%ocean%unstruct2struct(dx(:,:), tmp_incr)

    ! Clean up
    deallocate(pcv)
    if (allocated(tmp_incr)) deallocate(tmp_incr)

  end subroutine soca_corrad_apply

 
  !--------------------------------------------  
  subroutine soca_corr_exit(self)
    class(soca_bumpcorr2d_type),    intent(inout) :: self
    
    call self%horizconv%dealloc()

  end subroutine soca_corr_exit

end module soca_bumpcorr2d_mod

