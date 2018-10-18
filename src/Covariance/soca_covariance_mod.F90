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
  use type_bump
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices
  
  implicit none

  !> Fortran derived type to hold configuration data for the SOCA background/model covariance
  type :: soca_cov
     type(bump_type)  :: ocean_conv  !< Ocean convolution op from bump
     type(bump_type)  :: seaice_conv !< Seaice convolution op from bump
     !real(kind=kind_real), allocatable :: seaice_mask(:,:)
     integer, allocatable :: seaice_mask(:,:)          
     type(soca_field), pointer :: bkg         !< Background field (or first guess)
     logical          :: initialized = .false.
  end type soca_cov

#define LISTED_TYPE soca_cov

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
 type(registry_t) :: soca_cov_registry

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

  subroutine soca_cov_setup(self, c_conf, geom, bkg)

    use soca_geom_mod
    use iso_c_binding
    use config_mod
    use type_bump
    use type_nam
    
    implicit none

    type(soca_cov),        intent(inout) :: self   !< The covariance structure    
    type(c_ptr),              intent(in) :: c_conf !< The configuration
    type(soca_geom),          intent(in) :: geom   !< Geometry
    type(soca_field), target, intent(in) :: bkg    !< Background
    
    character(len=3)  :: domain
    integer :: is, ie, js, je, i, j
    
    ! Associate background
    self%bkg => bkg

    ! Define seaice mask from background seaice fraction
    call geom_get_domain_indices(bkg%geom%ocean, "compute", is, ie, js, je)
    allocate(self%seaice_mask(is:ie, js:je))
    self%seaice_mask = 0
    do i = is, ie
       do j = js, je
          if (sum(bkg%cicen(i, j, 2:), 1) * bkg%geom%ocean%mask2d(i, j) .gt. 0.0) then
             self%seaice_mask(i, j) = 1
          else
             self%seaice_mask(i, j) = 0
          end if
       end do
    end do

    ! Initialize ocean bump
    domain = 'ocn'
    call soca_bump_correlation(self, self%ocean_conv, geom, c_conf, domain)

    ! Initialize seaice bump
    domain = 'ice'
    call soca_bump_correlation(self, self%seaice_conv, geom, c_conf, domain)    

    self%initialized = .true.
    
  end subroutine soca_cov_setup

  ! ------------------------------------------------------------------------------

  !> Delete for the SOCA model's 3d error covariance matrices

  subroutine soca_cov_delete(self)

    use iso_c_binding

    implicit none
    type(soca_cov), intent(inout) :: self       !< The covariance structure        

    call self%ocean_conv%dealloc()
    call self%seaice_conv%dealloc()
    nullify(self%bkg)
    deallocate(self%seaice_mask)
    self%initialized = .false.
    
  end subroutine soca_cov_delete

  ! ------------------------------------------------------------------------------

  subroutine soca_cov_C_mult(self, dx)

    use kinds
    use type_bump
    
    implicit none

    type(soca_cov),   intent(inout) :: self !< The covariance structure    
    type(soca_field), intent(inout) :: dx   !< Input: Increment
                                            !< Output: C dx

    integer :: icat, izo
    
    ! Apply convolution to fields
    call soca_2d_convol(dx%ssh, self%ocean_conv, dx%geom)
    
    do icat = 1, dx%geom%ocean%ncat
       print *,'Apply nicas: aice, hice, category:',icat
       call soca_2d_convol(dx%cicen(:,:,icat+1), self%seaice_conv, dx%geom)
       call soca_2d_convol(dx%hicen(:,:,icat), self%seaice_conv, dx%geom)
    end do    

    do izo = 1,dx%geom%ocean%nzo
       print *,'Apply nicas: tocn, socn, layer:',izo
       call soca_2d_convol(dx%tocn(:,:,izo), self%ocean_conv, dx%geom)
       call soca_2d_convol(dx%socn(:,:,izo), self%ocean_conv, dx%geom)       
    end do    

  end subroutine soca_cov_C_mult
  
  ! ------------------------------------------------------------------------------

  subroutine soca_bump_correlation(self, horiz_convol, geom, c_conf, domain)
    use soca_geom_mod
    use type_bump
    use type_nam
    use iso_c_binding
    use oobump_mod, only: bump_read_conf
    use fckit_mpi_module, only: fckit_mpi_comm
    
    implicit none

    type(soca_cov),  intent(inout) :: self   !< The covariance structure
    type(bump_type), intent(inout) :: horiz_convol
    type(soca_geom),    intent(in) :: geom
    type(c_ptr),        intent(in) :: c_conf         !< Handle to configuration    
    character(len=3),   intent(in) :: domain
    
    !Grid stuff
    integer :: isc, iec, jsc, jec, jjj, jz, il, ib
    character(len=1024) :: subr = 'model_write'

    !bump stuff
    integer :: nc0a, nl0, nv, nts
    real(kind=kind_real), allocatable :: lon(:), lat(:), area(:), vunit(:,:)
    real(kind=kind_real), allocatable :: rosrad(:)    
    logical, allocatable :: lmask(:,:)
    integer, allocatable :: imask(:,:)    
    
    real(kind_real), allocatable :: rh(:,:,:,:)     !< Horizontal support radius for covariance (in m)
    real(kind_real), allocatable :: rv(:,:,:,:)     !< Vertical support radius for
    type(fckit_mpi_comm) :: f_comm

    f_comm = fckit_mpi_comm()

    !--- Initialize geometry to be passed to NICAS
    ! Indices for compute domain (no halo)
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

    ! Setup land mask
    jz = 1
    if (domain.eq.'ocn') then
       imask(1:nc0a,jz) = int(reshape( geom%ocean%mask2d(isc:iec, jsc:jec), (/nc0a/)))
    end if
    if (domain.eq.'ice') then
       imask(1:nc0a,jz) = int(reshape( self%seaice_mask(isc:iec, jsc:jec), (/nc0a/)))
    end if

    lmask = .false.
    where (imask.eq.1)
       lmask=.true.
    end where

    ! No vertical convolution, set dummy vertical unit    
    vunit = 1.0d0

    ! Setup horizontal decorrelation length scales
    ! TODO: Length scale should be in configuration
    allocate(rh(nc0a,nl0,nv,nts))
    allocate(rv(nc0a,nl0,nv,nts))
    if (domain.eq.'ocn') then
       do jjj=1,nc0a
          rh(jjj,1,1,1)=20.0d0*rosrad(jjj)
       end do
       where (rh<300d3)
          rh=300d3
       end where
    end if
    if (domain.eq.'ice') then
       rh = 500d3
    end if
    rv=1.0 ! Vertical scales not used, set to something

    ! Initialize bump namelist/parameters
    call horiz_convol%nam%init()
    call bump_read_conf(c_conf, horiz_convol)
    if (domain.eq.'ocn') horiz_convol%nam%prefix = 'ocn'
    if (domain.eq.'ice') horiz_convol%nam%prefix = 'ice'     

    ! Compute convolution weight    
    call horiz_convol%setup_online(f_comm%communicator(),nc0a,nl0,nv,nts,lon,lat,area,vunit,lmask)
    call horiz_convol%set_parameter('cor_rh',rh)       
    call horiz_convol%run_drivers()

    deallocate( lon, lat, area, vunit, imask, lmask )
    deallocate(rh,rv)       

  end subroutine soca_bump_correlation

  ! ------------------------------------------------------------------------------

  subroutine soca_2d_convol(dx, horiz_convol, geom)

    use soca_geom_mod
    use kinds
    use type_bump
    
    implicit none
    real(kind=kind_real), intent(inout) :: dx(:,:)
    type(bump_type),         intent(in) :: horiz_convol    
    type(soca_geom),         intent(in) :: geom        

    real(kind=kind_real), allocatable :: tmp_incr(:,:,:,:)

    ! Apply 2D convolution
    call soca_struct2unstruct(dx(:,:), geom, tmp_incr)
    call horiz_convol%apply_nicas(tmp_incr)
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
