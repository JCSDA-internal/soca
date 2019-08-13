! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_geom_mod

use MOM_domains, only : MOM_domain_type, MOM_infra_init
use MOM_io,      only : io_infra_init
use soca_mom6, only: soca_mom6_config, soca_mom6_init, soca_ice_column, &
                     soca_geomdomain_init
use soca_utils, only: write2pe
use kinds, only: kind_real
use fckit_kdtree_module
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use fms_io_mod, only : fms_io_init, fms_io_exit, &
                       register_restart_field, restart_file_type, &
                       restore_state, query_initialized, &
                       free_restart_type, save_restart
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, &
                            mpp_update_domains
use fms_mod,         only : write_data, read_data
use fms_io_mod,      only : fms_io_init, fms_io_exit

implicit none

private
public :: soca_geom
public :: soca_geom_registry
public :: geom_write, geom_get_domain_indices

!> Geometry data structure
type :: soca_geom
    type(MOM_domain_type), pointer :: Domain !< Ocean model domain
    type(soca_ice_column) :: ice_column !< Sea-ice geometry
    integer :: nx, ny, nzo
    integer :: nzi, nzs, ncat
    integer :: isc, iec, jsc, jec  !< indices of compute domain
    integer :: isd, ied, jsd, jed  !< indices of data domain
    integer :: iscl, iecl, jscl, jecl  !< indices of local compute domain
    integer :: isdl, iedl, jsdl, jedl  !< indices of local data domain
    real(kind=kind_real), allocatable, dimension(:,:) :: lon, lat  !< horizontal grid type
                                                                  !< 2D array of longitude, latitude
    real(kind=kind_real), allocatable, dimension(:,:) :: mask2d    !< 0 = land 1 = ocean
    real(kind=kind_real), allocatable, dimension(:,:) :: shoremask ! Includes shoreline as ocean point (useful for BUMP)
    integer,              allocatable, dimension(:,:) :: ij        ! index of ocean+shore line in compute grid
    real(kind=kind_real), allocatable, dimension(:,:) :: cell_area
    real(kind=kind_real), allocatable, dimension(:,:) :: rossby_radius
    logical :: save_local_domain = .false. ! If true, save the local geometry for each pe.
    contains
    procedure :: init => geom_init
    procedure :: end => geom_end
    procedure :: clone => geom_clone
    procedure :: print => geom_print
    procedure :: get_rossby_radius => geom_rossby_radius
    procedure :: validindex => geom_validindex
    procedure :: gridgen => geom_gridgen
    procedure :: thickness2depth => geom_thickness2depth
    procedure :: struct2unstruct => geom_struct2unstruct
    procedure :: unstruct2struct => geom_unstruct2struct
    procedure :: write => geom_write
end type soca_geom

#define LISTED_TYPE soca_geom

!> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

!> Global registry
type(registry_t) :: soca_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "Utils/linkedList_c.f"

! ------------------------------------------------------------------------------
!> Setup geometry object
subroutine geom_init(self, f_conf)
  class(soca_geom), intent(out) :: self
  type(fckit_configuration), intent(in) :: f_conf

  integer :: isave = 0

  ! Domain decomposition
  call soca_geomdomain_init(self%Domain, self%nzo)

  ! Initialize sea-ice grid
  if ( f_conf%has("num_ice_cat") ) &
      call f_conf%get_or_die("num_ice_cat", self%ice_column%ncat)
  if ( f_conf%has("num_ice_lev") ) &
      call f_conf%get_or_die("num_ice_lev", self%ice_column%nzi)
  if ( f_conf%has("num_sno_lev") ) &
      call f_conf%get_or_die("num_sno_lev", self%ice_column%nzs)

  ! Shortcuts to sea-ice grid size
  self%ncat = self%ice_column%ncat
  self%nzi = self%ice_column%nzi
  self%nzs = self%ice_column%nzs

  ! Allocate geometry arrays
  call geom_allocate(self)

  if ( f_conf%has("read_soca_grid") ) &
      call geom_read(self)

  ! Set output option for local geometry
  if ( f_conf%has("save_local_domain") ) &
      call f_conf%get_or_die("save_local_domain", isave)
  if ( isave == 1 ) self%save_local_domain = .true.

end subroutine geom_init

! ------------------------------------------------------------------------------
!> Geometry destructor
subroutine geom_end(self)
  class(soca_geom), intent(out)  :: self

  if (allocated(self%lon))           deallocate(self%lon)
  if (allocated(self%lat))           deallocate(self%lat)
  if (allocated(self%mask2d))        deallocate(self%mask2d)
  if (allocated(self%shoremask))     deallocate(self%shoremask)
  if (allocated(self%cell_area))     deallocate(self%cell_area)
  if (allocated(self%rossby_radius)) deallocate(self%rossby_radius)
  if (allocated(self%ij))            deallocate(self%ij)
  nullify(self%Domain)

end subroutine geom_end

! ------------------------------------------------------------------------------
!> Clone, self = other
subroutine geom_clone(self, other)
  class(soca_geom), intent( in) :: self
  class(soca_geom), intent(out) :: other

  ! Clone fms domain and vertical levels
  other%Domain => self%Domain
  other%nzo = self%nzo

  ! Clone sea-ice grid
  other%ice_column = self%ice_column
  other%ncat = self%ncat
  other%nzi = self%nzi
  other%nzs = self%nzs

  ! Allocate and clone geometry
  call geom_allocate(other)
  other%lon = self%lon
  other%lat = self%lat
  other%mask2d = self%mask2d
  other%cell_area = self%cell_area
  other%rossby_radius = self%rossby_radius

end subroutine geom_clone

! ------------------------------------------------------------------------------
!>
subroutine geom_gridgen(self)
  class(soca_geom), intent(inout) :: self

  type(soca_mom6_config) :: mom6_config

  ! Generate grid
  call soca_mom6_init(mom6_config, partial_init=.true.)
  self%lon = mom6_config%grid%GeoLonT
  self%lat = mom6_config%grid%GeoLatT
  self%mask2d = mom6_config%grid%mask2dT
  self%cell_area  = mom6_config%grid%areaT

  ! Get Rossby Radius
  call geom_rossby_radius(self)

  ! Output to file
  call geom_write(self)

end subroutine geom_gridgen

! ------------------------------------------------------------------------------
!> Allocate memory and point to mom6 data structure
subroutine geom_allocate(self)
  class(soca_geom), intent(inout) :: self

  integer :: nxny(2), nx, ny
  integer :: nzo, nzi, nzs
  integer :: isd, ied, jsd, jed

  ! Get domain shape (number of levels, indices of data and compute domain)
  call geom_get_domain_indices(self, "compute", self%isc, self%iec, self%jsc, self%jec)
  call geom_get_domain_indices(self, "data", isd, ied, jsd, jed)
  self%isd = isd ;  self%ied = ied ; self%jsd = jsd; self%jed = jed
  call geom_get_domain_indices(self, "compute", self%iscl, self%iecl, self%jscl, self%jecl, local=.true.)
  call geom_get_domain_indices(self, "data", self%isdl, self%iedl, self%jsdl, self%jedl, local=.true.)
  nzo = self%nzo

  ! Allocate arrays on compute domain
  allocate(self%lon(isd:ied,jsd:jed));           self%lon = 0.0_kind_real
  allocate(self%lat(isd:ied,jsd:jed));           self%lat = 0.0_kind_real
  allocate(self%mask2d(isd:ied,jsd:jed));        self%mask2d = 0.0_kind_real
  allocate(self%cell_area(isd:ied,jsd:jed));     self%cell_area = 0.0_kind_real
  allocate(self%rossby_radius(isd:ied,jsd:jed)); self%rossby_radius = 0.0_kind_real
  allocate(self%shoremask(isd:ied,jsd:jed));     self%shoremask = 0.0_kind_real

  ! Extract geometry of interest from model's data structure.
  ! Common to ocean & sea-ice
  ! No halo grid size
  ! TODO: Probably obsolete, remove
  nxny = shape( self%lon )
  nx = nxny(1)
  ny = nxny(2)
  self%nx = nx
  self%ny = ny

  ! Fill halo
  call mpp_update_domains(self%lon, self%Domain%mpp_domain)
  call mpp_update_domains(self%lat, self%Domain%mpp_domain)
  call mpp_update_domains(self%mask2d, self%Domain%mpp_domain)
  call mpp_update_domains(self%cell_area, self%Domain%mpp_domain)

end subroutine geom_allocate

! ------------------------------------------------------------------------------
!> Print geometry info to std output
subroutine geom_print(self)
  class(soca_geom), intent(in) :: self

  print *, 'nx=', self%nx
  print *, 'ny=', self%ny

end subroutine geom_print

! ------------------------------------------------------------------------------
!> Read and store Rossby Radius of deformation
!> TODO: Move out of geometry, use bilinear interp instead of nearest neighbor
subroutine geom_rossby_radius(self)
  class(soca_geom), intent(inout) :: self

  integer :: unit, i, j, n
  real(kind=kind_real) :: dum
  real(kind=kind_real), allocatable :: lon(:),lat(:),rr(:)
  type(kdtree) :: kd
  integer :: isc, iec, jsc, jec
  integer :: index(1), nn, io
  character(len=256) :: geom_output_file = "geom_output.nc"

  unit = 20
  open(unit=unit,file="rossrad.dat",status="old",action="read")
  n = 0
  do
     read(unit,*,iostat=io)
     if (io/=0) exit
     n = n+1
  end do
  rewind(unit)
  allocate(lon(n),lat(n),rr(n))
  do i = 1, n
     read(unit,*) lat(i),lon(i),dum,rr(i)
  end do
  close(unit)

  !--- Initialize kd-tree
  kd = kdtree_create(n, lon, lat)

  isc = self%isc ;  iec = self%iec ; jsc = self%jsc ; jec = self%jec

  !--- Find nearest neighbor
  nn=1 ! Num neighbors
  do i = isc, iec
     do j = jsc, jec
        call kdtree_k_nearest_neighbors(kd,self%lon(i,j),self%lat(i,j),1,index)
        self%rossby_radius(i,j)=rr(index(1))*1e3
     end do
  end do

  ! Release memory
  call kdtree_destroy(kd)

end subroutine geom_rossby_radius

! ------------------------------------------------------------------------------
!> Setup array of "valid index" to inline and pack structured geometry to
!> unstructured geometry
subroutine geom_validindex(self)
  ! Ignores inland mask grid points and
  ! select wet gridpoints and shoreline mask
  class(soca_geom), intent(inout) :: self
  integer :: i, j, ns, cnt
  integer :: isc, iec, jsc, jec
  real(kind=kind_real) :: shoretest

  ! Indices for compute domain (no halo)
  isc = self%isc ; iec = self%iec ; jsc = self%jsc ; jec = self%jec

  ! Extend mask 2 grid point inland TODO:NEED HALO FOR MASK!!!
  self%shoremask = self%mask2d
  do i = isc, iec
     do j = jsc, jec
        self%shoremask(i,j) = self%mask2d(i,j)
     end do
  end do

  ! Get number of valid points
  ns = int(sum(self%shoremask(isc:iec,jsc:jec)))
  allocate(self%ij(2,ns))

!!$    ! Save shoreline + ocean grid point
!!$    cnt = 1
!!$    do i = isc, iec
!!$       do j = jsc, jec
!!$          if (shoretest.gt.0.0d0) then
!!$             self%ij(1, cnt) = i
!!$             self%ij(2, cnt) = j
!!$             cnt = cnt + 1
!!$          end if
!!$       end do
!!$    end do

end subroutine geom_validindex

! ------------------------------------------------------------------------------
!> Write geometry to file
subroutine geom_write(self)
  class(soca_geom), intent(in) :: self

  character(len=256) :: geom_output_file = "soca_gridspec.nc"
  character(len=256) :: geom_output_pe
  integer :: pe
  character(len=8) :: fmt = '(I5.5)'
  character(len=1024) :: strpe
  integer :: ns
  integer :: idr_geom
  type(restart_file_type) :: geom_restart
  type(fckit_mpi_comm) :: f_comm

  ! Setup Communicator
  f_comm = fckit_mpi_comm()

  ! Save global domain
  call fms_io_init()
  idr_geom = register_restart_field(geom_restart, &
                                   &geom_output_file, &
                                   &'lon', &
                                   &self%lon(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &geom_output_file, &
                                   &'lat', &
                                   &self%lat(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &geom_output_file, &
                                   &'area', &
                                   &self%cell_area(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &geom_output_file, &
                                   &'rossby_radius', &
                                   &self%rossby_radius(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &geom_output_file, &
                                   &'mask2d', &
                                   &self%mask2d(:,:), &
                                   domain=self%Domain%mpp_domain)
  call save_restart(geom_restart, directory='')
  call free_restart_type(geom_restart)
  call fms_io_exit()

  if (self%save_local_domain) then
     ! Save local compute grid
     pe = f_comm%rank()

     write (strpe,fmt) pe
     geom_output_pe='geom_output_'//trim(strpe)//'.nc'

     ns = (self%iec - self%isc + 1) * (self%jec - self%jsc + 1 )
     call write2pe(reshape(self%shoremask,(/ns/)),'shoremask',geom_output_pe,.false.)
     call write2pe(reshape(self%mask2d,(/ns/)),'mask',geom_output_pe,.true.)
     call write2pe(reshape(self%lon,(/ns/)),'lon',geom_output_pe,.true.)
     call write2pe(reshape(self%lat,(/ns/)),'lat',geom_output_pe,.true.)
  end if

end subroutine geom_write

! ------------------------------------------------------------------------------
!> Read geometry from file
subroutine geom_read(self)
  class(soca_geom), intent(in) :: self

  character(len=256) :: geom_input_file = "soca_gridspec.nc"
  integer :: idr_geom
  type(restart_file_type) :: geom_restart

  call fms_io_init()
  idr_geom = register_restart_field(geom_restart, &
                                   &geom_input_file, &
                                   &'lon', &
                                   &self%lon(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &geom_input_file, &
                                   &'lat', &
                                   &self%lat(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &geom_input_file, &
                                   &'area', &
                                   &self%cell_area(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &geom_input_file, &
                                   &'rossby_radius', &
                                   &self%rossby_radius(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &geom_input_file, &
                                   &'mask2d', &
                                   &self%mask2d(:,:), &
                                   domain=self%Domain%mpp_domain)
  call restore_state(geom_restart, directory='')
  call free_restart_type(geom_restart)
  call fms_io_exit()

end subroutine geom_read

! ------------------------------------------------------------------------------
!> Get indices for compute or data domain
subroutine geom_get_domain_indices(self, domain_type, is, ie, js, je, local)
  class(soca_geom), intent(in) :: self
  character(len=*),       intent(in) :: domain_type
  integer,               intent(out) :: is, ie, js, je
  logical,      optional, intent(in) :: local

  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed

  call mpp_get_compute_domain(self%Domain%mpp_domain,isc,iec,jsc,jec)
  call mpp_get_data_domain(self%Domain%mpp_domain,isd,ied,jsd,jed)
  if (present(local)) then
     isc = isc - (isd-1) ; iec = iec - (isd-1) ; ied = ied - (isd-1) ; isd = 1
     jsc = jsc - (jsd-1) ; jec = jec - (jsd-1) ; jed = jed - (jsd-1) ; jsd = 1
  end if

  select case (trim(domain_type))
  case ("compute")
     is = isc; ie = iec; js = jsc; je = jec;
  case ("data")
     is = isd; ie = ied; js = jsd; je = jed;

  end select

end subroutine geom_get_domain_indices

! ------------------------------------------------------------------------------
!> Get layer depth from layer thicknesses
subroutine geom_thickness2depth(self, h, z)
  class(soca_geom),     intent(in   ) :: self
  real(kind=kind_real), intent(in   ) :: h(:,:,:) ! Layer thickness
  real(kind=kind_real), intent(inout) :: z(:,:,:) ! Mid-layer depth

  integer :: is, ie, js, je, i, j, k

  ! Should check shape of z
  is = lbound(h,dim=1)
  ie = ubound(h,dim=1)
  js = lbound(h,dim=2)
  je = ubound(h,dim=2)

  !allocate(z(is:ie, js:je, self%nzo))

  do i = is, ie
     do j = js, je
        do k = 1, self%nzo
           if (k.eq.1) then
              z(i,j,k) = 0.5_kind_real*h(i,j,k)
           else
              z(i,j,k) = sum(h(i,j,1:k-1))+0.5_kind_real*h(i,j,k)
           end if
        end do
     end do
  end do
end subroutine geom_thickness2depth

! ------------------------------------------------------------------------------

subroutine geom_struct2unstruct(self, dx_struct, dx_unstruct)
  class(soca_geom),                  intent(in ) :: self
  real(kind=kind_real),              intent(in ) :: dx_struct(:,:)
  real(kind=kind_real), allocatable, intent(out) :: dx_unstruct(:,:,:,:)

  integer :: nc0a

  nc0a = (self%iecl - self%iscl + 1) * (self%jecl - self%jscl + 1 )
  allocate(dx_unstruct(nc0a,1,1,1))
  dx_unstruct = &
  &  reshape(dx_struct(self%iscl:self%iecl, self%jscl:self%jecl), (/nc0a,1,1,1/))

end subroutine geom_struct2unstruct

! ------------------------------------------------------------------------------

subroutine geom_unstruct2struct(self, dx_struct, dx_unstruct)
  class(soca_geom),                  intent(in   ) :: self
  real(kind=kind_real),              intent(inout) :: dx_struct(:,:)
  real(kind=kind_real), allocatable, intent(inout) :: dx_unstruct(:,:,:,:)

  dx_struct(self%iscl:self%iecl, self%jscl:self%jecl) = &
  & reshape(dx_unstruct, (/size(dx_struct(self%iscl:self%iecl, self%jscl:self%jecl),1), &
  &                        size(dx_struct(self%iscl:self%iecl, self%jscl:self%jecl),2)/))

  deallocate(dx_unstruct)

end subroutine geom_unstruct2struct

! ------------------------------------------------------------------------------

end module soca_geom_mod
