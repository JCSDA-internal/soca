! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!> Geometry module
module soca_geom_mod

  ! jedi modules
use atlas_module, only: atlas_functionspace_NodeColumns, atlas_fieldset, &
    atlas_field, atlas_real, atlas_integer, atlas_geometry, atlas_indexkdtree
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use kinds, only: kind_real
use type_fieldset, only: fieldset_type

! mom6 / fms modules
use fms_io_mod, only : fms_io_init, fms_io_exit, &
                       register_restart_field, restart_file_type, &
                       restore_state, free_restart_type, save_restart
use MOM_diag_remap,  only : diag_remap_ctrl, diag_remap_init, diag_remap_configure_axes, &
                            diag_remap_end, diag_remap_update
use MOM_domains, only : MOM_domain_type
use MOM_EOS,         only : EOS_type
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, &
                            mpp_get_global_domain, mpp_update_domains, &
                            CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE                            
! soca modules
use soca_fields_metadata_mod, only : soca_fields_metadata
use soca_mom6, only: soca_mom6_config, soca_mom6_init, soca_geomdomain_init
use soca_utils, only: write2pe, soca_remap_idw


implicit none
private


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

! A default value to indicate that the halo region has not been updated
! during a halo exchange.
! Note, it seems we had been using invalud halo points with lat/lon of 0,0
! oops!
real(kind=kind_real), parameter :: INVALID_HALO = -999_kind_real

! ------------------------------------------------------------------------------
!> Geometry data structure
type, public :: soca_geom
    type(MOM_domain_type), pointer :: Domain !< Ocean model domain
    integer :: nzo, nzo_zstar

    !> \name local domain indices
    !! \{
    integer :: isc, iec, jsc, jec
    !> \}

    !> \name data domain indices
    !! \{
    integer :: isd, ied, jsd, jed
    !> \}

    !> \name global domain indices
    !! \{
    integer :: isg, ieg, jsg, jeg
    !> \}

    !> \name local compute domain indices
    !! \{
    integer :: iscl, iecl, jscl, jecl
    !> \}

    !> \name local data domain indices
    !! \{
    integer :: isdl, iedl, jsdl, jedl
    !> \}

    !> \name iterator dimension
    !! \{
    integer :: iterator_dimension
    !> \}

    !> \name grid latitude/longitude
    !! \{
    real(kind=kind_real), allocatable, dimension(:)   :: lonh !< cell center nominal longitude
    real(kind=kind_real), allocatable, dimension(:)   :: lath !< cell center nominal latitude
    real(kind=kind_real), allocatable, dimension(:)   :: lonq !< cell corner nominal longitude
    real(kind=kind_real), allocatable, dimension(:)   :: latq !< cell corner nominal latitude
    real(kind=kind_real), allocatable, dimension(:,:) :: lon  !< Tracer grid longitude
    real(kind=kind_real), allocatable, dimension(:,:) :: lat  !< Tracer grid latitude
    real(kind=kind_real), allocatable, dimension(:,:) :: lonu !< U grid longitude
    real(kind=kind_real), allocatable, dimension(:,:) :: latu !< U grid latitude
    real(kind=kind_real), allocatable, dimension(:,:) :: lonv !< V grid longitude
    real(kind=kind_real), allocatable, dimension(:,:) :: latv !< V grid latitude
    !> \}

    !> \name ocean/land masks
    !! \{

    !> mask for tracer grid. 0 = land 1 = ocean
    real(kind=kind_real), allocatable, dimension(:,:) :: mask2d
    !> mask for U grid. 0 = land 1 = ocean
    real(kind=kind_real), allocatable, dimension(:,:) :: mask2du
    !> mask for V grid. 0 = land 1 = ocean
    real(kind=kind_real), allocatable, dimension(:,:) :: mask2dv
    !> \}

    !> \name other grid properties
    !! \{
    real(kind=kind_real), allocatable, dimension(:,:) :: sin_rot !< sine of rotation between logical grid north
    real(kind=kind_real), allocatable, dimension(:,:) :: cos_rot !< cosine of rotation between logical grid north
    real(kind=kind_real), allocatable, dimension(:,:) :: dx !< cell x width (m)
    real(kind=kind_real), allocatable, dimension(:,:) :: dy !< cell y width (m)
    real(kind=kind_real), allocatable, dimension(:,:) :: cell_area !< cell area (m^2)
    real(kind=kind_real), allocatable, dimension(:,:) :: rossby_radius !< rossby radius (m) at the gridpoint
    real(kind=kind_real), allocatable, dimension(:,:) :: distance_from_coast !< distance to closest land grid point (m)
    real(kind=kind_real), allocatable, dimension(:,:,:) :: h !< layer thickness (m)
    real(kind=kind_real), allocatable, dimension(:,:,:) :: h_zstar
    !> \}

    !> instance of the metadata that is read in from a config file upon initialization
    type(soca_fields_metadata) :: fields_metadata

    logical, private :: save_local_domain = .false. !< If true, save the local geometry for each pe.
    character(len=:), allocatable :: geom_grid_file !< filename of geometry
    character(len=:), allocatable :: rossby_file !< filename of rossby radius input file (if used)

    type(fckit_mpi_comm) :: f_comm !< MPI communicator

    !> mesh parameters
    type(atlas_functionspace_NodeColumns) :: functionspaceInchalo
    integer, allocatable :: atlas_idx2i(:), atlas_idx2j(:)
    type(atlas_fieldset) :: fieldset !< the geom fields (area, mask, etc)

  contains


    !> \copybrief soca_geom_init \see soca_geom_init
    procedure :: init => soca_geom_init

    !> \copybrief soca_geom_end \see soca_geom_end
    procedure :: end => soca_geom_end

    !> \copybrief soca_geom_fill_atlas_fieldset \see soca_geom_fill_atlas_fieldset
    procedure :: init_fieldset => soca_geom_init_fieldset

    !> \copybrief soca_geom_clone \see soca_geom_clone
    procedure :: clone => soca_geom_clone

    !> \copybrief soca_geom_gridgen \see soca_geom_gridgen
    procedure :: gridgen => soca_geom_gridgen

    !> \copybrief soca_geom_thickness2depth \see soca_geom_thickness2depth
    procedure :: thickness2depth => soca_geom_thickness2depth

    !> \copybrief soca_geom_write \see soca_geom_write
    procedure :: write => soca_geom_write

    ! TODO make private
    procedure :: mesh_valid_nodes_cells => soca_geom_mesh_valid_nodes_cells

    procedure :: atlas_idx2ij => soca_geom_atlas_idx2ij
end type soca_geom

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup geometry object
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_init(self, f_conf, f_comm)
  class(soca_geom),         intent(out) :: self
  type(fckit_configuration), intent(in) :: f_conf
  type(fckit_mpi_comm),   intent(in)    :: f_comm !< MPI communicator for this geometry

  character(len=:), allocatable :: str
  logical :: full_init = .false.

  ! MPI communicator
  self%f_comm = f_comm

  ! Domain decomposition
  call soca_geomdomain_init(self%Domain, self%nzo, f_comm)

  ! User-defined grid filename
  if ( .not. f_conf%get("geom_grid_file", self%geom_grid_file) ) &
     self%geom_grid_file = "soca_gridspec.nc" ! default if not found
  if ( .not. f_conf%get("rossby file", self%rossby_file)) &
    self%rossby_file = "rossrad.dat" ! default if not found

  ! Allocate geometry arrays
  call soca_geom_allocate(self)

  ! Check if a full initialization is required, default to false
  if ( .not. f_conf%get("full_init", full_init) ) full_init = .false.

  ! Read the geometry from file by default,
  ! skip this step if a full init is required
  if ( .not. full_init) call soca_geom_read(self) 

  ! Fill halo
  call mpp_update_domains(self%lon, self%Domain%mpp_domain)
  call mpp_update_domains(self%lat, self%Domain%mpp_domain)
  call mpp_update_domains(self%lonu, self%Domain%mpp_domain)
  call mpp_update_domains(self%latu, self%Domain%mpp_domain)
  call mpp_update_domains(self%lonv, self%Domain%mpp_domain)
  call mpp_update_domains(self%latv, self%Domain%mpp_domain)
  call mpp_update_domains(self%sin_rot, self%Domain%mpp_domain)
  call mpp_update_domains(self%cos_rot, self%Domain%mpp_domain)
  call mpp_update_domains(self%mask2d, self%Domain%mpp_domain)
  call mpp_update_domains(self%mask2du, self%Domain%mpp_domain)
  call mpp_update_domains(self%mask2dv, self%Domain%mpp_domain)
  call mpp_update_domains(self%dx, self%Domain%mpp_domain)
  call mpp_update_domains(self%dy, self%Domain%mpp_domain)
  call mpp_update_domains(self%cell_area, self%Domain%mpp_domain)
  call mpp_update_domains(self%rossby_radius, self%Domain%mpp_domain)
  call mpp_update_domains(self%distance_from_coast, self%Domain%mpp_domain)

  ! Set output option for local geometry
  if ( .not. f_conf%get("save_local_domain", self%save_local_domain) ) &
     self%save_local_domain = .false.

  ! process the fields metadata file
  call f_conf%get_or_die("fields metadata", str)
  call self%fields_metadata%create(str)

  ! retrieve iterator dimension from config
  if ( .not. f_conf%get("iterator dimension", self%iterator_dimension) ) &
      self%iterator_dimension = 2

end subroutine soca_geom_init

! ------------------------------------------------------------------------------
!> Geometry destructor
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_end(self)
  class(soca_geom), intent(out)  :: self

  if (allocated(self%lonh))          deallocate(self%lonh)
  if (allocated(self%lath))          deallocate(self%lath)
  if (allocated(self%lonq))          deallocate(self%lonq)
  if (allocated(self%latq))          deallocate(self%latq)
  if (allocated(self%lon))           deallocate(self%lon)
  if (allocated(self%lat))           deallocate(self%lat)
  if (allocated(self%lonu))          deallocate(self%lonu)
  if (allocated(self%latu))          deallocate(self%latu)
  if (allocated(self%lonv))          deallocate(self%lonv)
  if (allocated(self%latv))          deallocate(self%latv)
  if (allocated(self%sin_rot))       deallocate(self%sin_rot)
  if (allocated(self%cos_rot))       deallocate(self%cos_rot)
  if (allocated(self%mask2d))        deallocate(self%mask2d)
  if (allocated(self%mask2du))       deallocate(self%mask2du)
  if (allocated(self%mask2dv))       deallocate(self%mask2dv)
  if (allocated(self%dx))            deallocate(self%dx)
  if (allocated(self%dy))            deallocate(self%dy)
  if (allocated(self%cell_area))     deallocate(self%cell_area)
  if (allocated(self%rossby_radius)) deallocate(self%rossby_radius)
  if (allocated(self%distance_from_coast)) deallocate(self%distance_from_coast)
  if (allocated(self%h))             deallocate(self%h)
  if (allocated(self%h_zstar))       deallocate(self%h_zstar)
  nullify(self%Domain)
  call self%functionspaceIncHalo%final()

end subroutine soca_geom_end


! --------------------------------------------------------------------------------------------------
!> Fill ATLAS fieldset
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_init_fieldset(self)
  class(soca_geom),  intent(inout) :: self

  integer :: i, j, n, jz
  type(atlas_field) :: fArea, fInterpMask, fVertCoord, fGmask, fOwned, fRossby
  real(kind=kind_real), pointer :: vArea(:,:), vInterpMask(:,:), vVertCoord(:,:), vRossby(:,:)
  integer, pointer :: vGmask(:,:), vOwned(:,:)

  ! create fields, get pointers to their data
  fArea = self%functionspaceInchalo%create_field(name='area', kind=atlas_real(kind_real), levels=1)
  call self%fieldset%add(fArea)
  call fArea%data(vArea)
  
  fInterpMask = self%functionspaceInchalo%create_field(name='interp_mask', kind=atlas_real(kind_real), levels=1)
  call self%fieldset%add(fInterpMask)    
  call fInterpMask%data(vInterpMask)

  fVertCoord = self%functionspaceInchalo%create_field(name='vert_coord', kind=atlas_real(kind_real), levels=self%nzo)
  call self%fieldset%add(fVertCoord)
  call fVertCoord%data(vVertCoord) 

  fGmask = self%functionspaceInchalo%create_field(name='gmask', kind=atlas_integer(kind(0)), levels=self%nzo)
  call self%fieldset%add(fGmask)
  call fGmask%data(vGmask) 

  fOwned = self%functionspaceInchalo%create_field(name='owned', kind=atlas_integer(kind(0)), levels=1)
  call self%fieldset%add(fOwned)
  call fOwned%data(vOwned)

  fRossby = self%functionspaceInchalo%create_field(name='rossby_radius', kind=atlas_real(kind_real), levels=1)
  call self%fieldset%add(fRossby)
  call fRossby%data(vRossby)

  ! set the data
  vOwned = 0 ! need to set to 0, it's the only one that doesn't get halo update
  do n=1,size(self%atlas_idx2i)
    if(.not. self%atlas_idx2ij(n, i, j)) cycle
    vArea(1,n) = self%cell_area(i,j)
    vInterpMask(1,n) = self%mask2d(i,j)
    vGmask(:, n) = int(self%mask2d(i,j))
    vOwned(1, n) = 1
    vRossby(1, n) = self%rossby_radius(i,j)
  end do
  do jz=1,self%nzo
    vVertCoord(jz,:) = real(jz, kind_real)
  end do

  ! halo exchanges for some of the fields 
  call fGmask%halo_exchange()
  call fInterpMask%halo_exchange()
  call fArea%halo_exchange()
  call fVertCoord%halo_exchange()
  call fRossby%halo_exchange()

  ! done, cleanup
  call fArea%final()
  call fInterpMask%final()
  call fVertCoord%final()
  call fGmask%final()
  call fOwned%final()
  call fRossby%final()

end subroutine soca_geom_init_fieldset


! ------------------------------------------------------------------------------
!> Clone, self = other
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_clone(self, other)
  class(soca_geom), intent(inout) :: self
  class(soca_geom), intent(in) :: other

  ! Clone communicator
  self%f_comm = other%f_comm

  ! Clone fms domain and vertical levels
  self%Domain => other%Domain
  self%nzo = other%nzo

  !
  self%geom_grid_file = other%geom_grid_file
  self%rossby_file = other%rossby_file

  self%iterator_dimension = other%iterator_dimension

  ! Allocate and clone geometry
  call soca_geom_allocate(self)
  self%lonh = other%lonh
  self%lath = other%lath
  self%lonq = other%lonq
  self%latq = other%latq
  self%lon = other%lon
  self%lat = other%lat
  self%lonu = other%lonu
  self%latu = other%latu
  self%lonv = other%lonv
  self%latv = other%latv
  self%sin_rot = other%sin_rot
  self%cos_rot = other%cos_rot
  self%mask2d = other%mask2d
  self%mask2du = other%mask2du
  self%mask2dv = other%mask2dv
  self%dx = other%dx
  self%dy = other%dy
  self%cell_area = other%cell_area
  self%rossby_radius = other%rossby_radius
  self%distance_from_coast = other%distance_from_coast
  self%h = other%h
  call self%fields_metadata%clone(other%fields_metadata)
end subroutine soca_geom_clone


! ------------------------------------------------------------------------------
!> Generate the grid with the help of mom6, and save it to a file for use later.
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_gridgen(self)
  class(soca_geom), intent(inout) :: self

  ! allocate variables for regridding to zstar coord
  type(soca_mom6_config) :: mom6_config
  type(diag_remap_ctrl) :: remap_ctrl
  type(EOS_type), pointer :: eqn_of_state
  integer :: k
  real(kind=kind_real), allocatable :: tracer(:,:,:)
  logical :: answers_2018 = .false.

  ! Generate grid
  call soca_mom6_init(mom6_config, partial_init=.true.)
  self%lonh = mom6_config%grid%gridlont
  self%lath = mom6_config%grid%gridlatt
  self%lonq = mom6_config%grid%gridlonb
  self%latq = mom6_config%grid%gridlatb
  self%lon = mom6_config%grid%GeoLonT
  self%lat = mom6_config%grid%GeoLatT
  self%lonu = mom6_config%grid%geoLonCu
  self%latu = mom6_config%grid%geoLatCu
  self%lonv = mom6_config%grid%geoLonCv
  self%latv = mom6_config%grid%geoLatCv

  self%sin_rot = mom6_config%grid%sin_rot
  self%cos_rot = mom6_config%grid%cos_rot

  self%mask2d = mom6_config%grid%mask2dT
  self%mask2du = mom6_config%grid%mask2dCu
  self%mask2dv = mom6_config%grid%mask2dCv
  self%dx = mom6_config%grid%dxT
  self%dy = mom6_config%grid%dyT
  self%cell_area  = mom6_config%grid%areaT
  self%h = mom6_config%MOM_CSp%h

  ! Setup intermediate zstar coordinate
  allocate(tracer(self%isd:self%ied, self%jsd:self%jed, self%nzo))
  tracer = 0.d0 ! dummy tracer
  call diag_remap_init(remap_ctrl, coord_tuple='ZSTAR, ZSTAR, ZSTAR', answers_2018=answers_2018)
  call diag_remap_configure_axes(remap_ctrl, mom6_config%GV, mom6_config%scaling, mom6_config%param_file)
  self%nzo_zstar = remap_ctrl%nz
  if (allocated(self%h_zstar)) deallocate(self%h_zstar)
  allocate(self%h_zstar(self%isd:self%ied, self%jsd:self%jed, 1:remap_ctrl%nz))

  ! Compute intermediate vertical coordinate self%h_zstar
  call diag_remap_update(remap_ctrl, &
                         mom6_config%grid, &
                         mom6_config%GV, &
                         mom6_config%scaling, &
                         self%h, tracer, tracer, eqn_of_state, self%h_zstar)
  call diag_remap_end(remap_ctrl)

  ! Get Rossby Radius
  call soca_geom_rossby_radius(self, self%rossby_file)

  call soca_geom_distance_from_coast(self)

  ! Output to file
  call soca_geom_write(self)

end subroutine soca_geom_gridgen


! ------------------------------------------------------------------------------
!> Allocate memory and point to mom6 data structure
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_allocate(self)
  class(soca_geom), intent(inout) :: self

  integer :: nzo
  integer :: isd, ied, jsd, jed

  nzo = self%nzo

  ! Get domain shape (number of levels, indices of data and compute domain)
  call soca_geom_get_domain_indices(self, "compute", self%isc, self%iec, self%jsc, self%jec)
  call soca_geom_get_domain_indices(self, "data", isd, ied, jsd, jed)
  self%isd = isd ;  self%ied = ied
  self%jsd = jsd; self%jed = jed
  call soca_geom_get_domain_indices(self, "global", self%isg, self%ieg, self%jsg, self%jeg)
  call soca_geom_get_domain_indices(self, "compute", self%iscl, self%iecl, self%jscl, self%jecl, local=.true.)
  call soca_geom_get_domain_indices(self, "data", self%isdl, self%iedl, self%jsdl, self%jedl, local=.true.)

  ! Allocate arrays on compute domain
  ! NOTE, sometimes MOM6 doesn't use all the halo points, we set 
  !  lat/lon to INVALID_HALO so that we can detect this later when deciding
  !  which halo points are invalid/duplicates
  allocate(self%lonh(self%isg:self%ieg));        self%lonh = 0.0_kind_real
  allocate(self%lath(self%jsg:self%jeg));        self%lath = 0.0_kind_real
  allocate(self%lonq(self%isg:self%ieg));        self%lonq = 0.0_kind_real
  allocate(self%latq(self%jsg:self%jeg));        self%latq = 0.0_kind_real
  allocate(self%lon(isd:ied,jsd:jed));           self%lon = INVALID_HALO
  allocate(self%lat(isd:ied,jsd:jed));           self%lat = INVALID_HALO
  allocate(self%lonu(isd:ied,jsd:jed));          self%lonu = 0.0_kind_real
  allocate(self%latu(isd:ied,jsd:jed));          self%latu = 0.0_kind_real
  allocate(self%lonv(isd:ied,jsd:jed));          self%lonv = 0.0_kind_real
  allocate(self%latv(isd:ied,jsd:jed));          self%latv = 0.0_kind_real

  allocate(self%sin_rot(isd:ied,jsd:jed));       self%sin_rot = 0.0_kind_real
  allocate(self%cos_rot(isd:ied,jsd:jed));       self%cos_rot = 0.0_kind_real

  allocate(self%mask2d(isd:ied,jsd:jed));        self%mask2d = 0.0_kind_real
  allocate(self%mask2du(isd:ied,jsd:jed));       self%mask2du = 0.0_kind_real
  allocate(self%mask2dv(isd:ied,jsd:jed));       self%mask2dv = 0.0_kind_real

  allocate(self%dx(isd:ied,jsd:jed));            self%dx = 0.0_kind_real  
  allocate(self%dy(isd:ied,jsd:jed));            self%dy = 0.0_kind_real  
  allocate(self%cell_area(isd:ied,jsd:jed));     self%cell_area = 0.0_kind_real
  allocate(self%rossby_radius(isd:ied,jsd:jed)); self%rossby_radius = 0.0_kind_real
  allocate(self%distance_from_coast(isd:ied,jsd:jed)); self%distance_from_coast = 0.0_kind_real
  allocate(self%h(isd:ied,jsd:jed,1:nzo));       self%h = 0.0_kind_real

end subroutine soca_geom_allocate


! ------------------------------------------------------------------------------
!> Calcuate distance from coast for the ocean points
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_distance_from_coast(self)
  class(soca_geom), intent(inout) :: self

  type(atlas_indexkdtree) :: kd
  type(atlas_geometry) :: ageometry
  integer :: i, j, idx(1)
  integer :: num_land_l, num_land
  integer, allocatable :: rcvcnt(:), displs(:)
  real(kind=kind_real), allocatable :: land_lon(:), land_lat(:)
  real(kind=kind_real), allocatable :: land_lon_l(:), land_lat_l(:)
  real(kind=kind_real) :: closest_lon, closest_lat

  ! collect lat/lon of all land points on all procs
  ! (use the tracer grid and mask for this)
  ! --------------------------------------------------
  allocate(rcvcnt(self%f_comm%size()))
  allocate(displs(self%f_comm%size()))

  num_land_l = count(self%mask2d(self%isc:self%iec, self%jsc:self%jec)==0.0)
  call self%f_comm%allgather(num_land_l, rcvcnt)
  num_land = sum(rcvcnt)

  displs(1) = 0
  do j = 2, self%f_comm%size()
    displs(j) = displs(j-1) + rcvcnt(j-1)
  enddo

  allocate(land_lon_l(num_land_l))
  allocate(land_lat_l(num_land_l))
  land_lon_l = pack(self%lon(self%isc:self%iec, self%jsc:self%jec), &
                  mask=self%mask2d(self%isc:self%iec,self%jsc:self%jec)==0.0)
  land_lat_l = pack(self%lat(self%isc:self%iec, self%jsc:self%jec), &
                  mask=self%mask2d(self%isc:self%iec,self%jsc:self%jec)==0.0)
  allocate(land_lon(num_land))
  allocate(land_lat(num_land))

  call self%f_comm%allgather(land_lon_l, land_lon, num_land_l, rcvcnt, displs)
  call self%f_comm%allgather(land_lat_l, land_lat, num_land_l, rcvcnt, displs)


  ! pass land points to the kd tree
  !---------------------------------------
  ageometry = atlas_geometry("Earth") !< TODO: remove this hardcoded value so
                                      ! we can do DA on Europa at some point.
                                      ! (Next AOP??)
  kd = atlas_indexkdtree(ageometry)
  call kd%reserve(num_land)
  call kd%build(num_land, land_lon, land_lat)


  ! for each point in local domain, lookup distance to nearest land point
  ! ---------------------------------------
  do i = self%isc, self%iec
    do j = self%jsc, self%jec
      call kd%closestPoints( self%lon(i,j), self%lat(i,j), 1, idx )
      self%distance_from_coast(i,j) = ageometry%distance( &
            self%lon(i,j), self%lat(i,j), land_lon(idx(1)), land_lat(idx(1)))
    enddo
  enddo

  ! cleanup
  call kd%final()

end subroutine


! ------------------------------------------------------------------------------
!> Read and store Rossby Radius of deformation
!!
!! Input data is interpolated to the current grid.
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_rossby_radius(self, filename)
  class(soca_geom),           intent(inout) :: self
  character(len=:), allocatable, intent(in) :: filename

  integer :: unit, i, n
  real(kind=kind_real) :: dum
  real(kind=kind_real), allocatable :: lon(:),lat(:),rr(:)
  integer :: isc, iec, jsc, jec
  integer :: io

  ! read in the file
  unit = 20
  open(unit=unit,file=filename,status="old",action="read")
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

  ! convert to meters
  rr = rr * 1e3

  ! remap
  isc = self%isc ;  iec = self%iec ; jsc = self%jsc ; jec = self%jec
  call soca_remap_idw(lon, lat, rr, self%lon(isc:iec,jsc:jec), &
                      self%lat(isc:iec,jsc:jec), self%rossby_radius(isc:iec,jsc:jec) )

end subroutine soca_geom_rossby_radius


! ------------------------------------------------------------------------------
!> Write geometry to file
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_write(self)
  class(soca_geom), intent(in) :: self

  character(len=256) :: geom_output_pe
  integer :: pe
  character(len=8) :: fmt = '(I5.5)'
  character(len=1024) :: strpe
  integer :: ns
  integer :: idr_geom
  type(restart_file_type) :: geom_restart

  ! Save global domain
  call fms_io_init()
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lonh', &
                                   &self%lonh(:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lath', &
                                   &self%lath(:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lonq', &
                                   &self%lonq(:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'latq', &
                                   &self%latq(:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lon', &
                                   &self%lon(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lat', &
                                   &self%lat(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lonu', &
                                   &self%lonu(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'latu', &
                                   &self%latu(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lonv', &
                                   &self%lonv(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'latv', &
                                   &self%latv(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'sin_rot', &
                                   &self%sin_rot(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'cos_rot', &
                                   &self%cos_rot(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'dx', &
                                   &self%dx(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'dy', &
                                   &self%dy(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'area', &
                                   &self%cell_area(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'rossby_radius', &
                                   &self%rossby_radius(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'mask2d', &
                                   &self%mask2d(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'mask2du', &
                                   &self%mask2du(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'mask2dv', &
                                   &self%mask2dv(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'h', &
                                   &self%h(:,:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'nzo_zstar', &
                                   &self%nzo_zstar, &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'h_zstar', &
                                   &self%h_zstar(:,:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'distance_from_coast', &
                                   &self%distance_from_coast(:,:), &
                                   domain=self%Domain%mpp_domain)

  call save_restart(geom_restart, directory='')
  call free_restart_type(geom_restart)
  call fms_io_exit()

  if (self%save_local_domain) then
     ! Save local compute grid
     pe = self%f_comm%rank()

     write (strpe,fmt) pe
     geom_output_pe='geom_output_'//trim(strpe)//'.nc'

     ns = (self%iec - self%isc + 1) * (self%jec - self%jsc + 1 )
     call write2pe(reshape(self%mask2d(self%isc:self%iec,self%jsc:self%jec),(/ns/)),'mask',geom_output_pe,.false.)
     call write2pe(reshape(self%lon(self%isc:self%iec,self%jsc:self%jec),(/ns/)),'lon',geom_output_pe,.true.)
     call write2pe(reshape(self%lat(self%isc:self%iec,self%jsc:self%jec),(/ns/)),'lat',geom_output_pe,.true.)
  end if

end subroutine soca_geom_write


! ------------------------------------------------------------------------------
!> Read geometry from file
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_read(self)
  class(soca_geom), intent(inout) :: self

  integer :: idr_geom
  type(restart_file_type) :: geom_restart

  call fms_io_init()
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lonh', &
                                   &self%lonh(:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lath', &
                                   &self%lath(:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lonq', &
                                   &self%lonq(:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'latq', &
                                   &self%latq(:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lon', &
                                   &self%lon(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lat', &
                                   &self%lat(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lonu', &
                                   &self%lonu(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'latu', &
                                   &self%latu(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'lonv', &
                                   &self%lonv(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'latv', &
                                   &self%latv(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'sin_rot', &
                                   &self%sin_rot(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'cos_rot', &
                                   &self%cos_rot(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'dx', &
                                   &self%dx(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'dy', &
                                   &self%dy(:,:), &
                                   domain=self%Domain%mpp_domain)                                   
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'area', &
                                   &self%cell_area(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'rossby_radius', &
                                   &self%rossby_radius(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'mask2d', &
                                   &self%mask2d(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'mask2du', &
                                   &self%mask2du(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'mask2dv', &
                                   &self%mask2dv(:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'h', &
                                   &self%h(:,:,:), &
                                   domain=self%Domain%mpp_domain)
  idr_geom = register_restart_field(geom_restart, &
                                   &self%geom_grid_file, &
                                   &'distance_from_coast', &
                                   &self%distance_from_coast(:,:), &
                                   domain=self%Domain%mpp_domain)
  call restore_state(geom_restart, directory='')
  call free_restart_type(geom_restart)
  call fms_io_exit()

end subroutine soca_geom_read


! ------------------------------------------------------------------------------
!> Get indices for compute or data domain
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_get_domain_indices(self, domain_type, is, ie, js, je, local)
  class(soca_geom), intent(in) :: self
  character(len=*),       intent(in) :: domain_type !< "compute", "data", or "global"
  integer,               intent(out) :: is, ie, js, je
  logical,      optional, intent(in) :: local

  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: isg, ieg, jsg, jeg

  call mpp_get_compute_domain(self%Domain%mpp_domain,isc,iec,jsc,jec)
  call mpp_get_data_domain(self%Domain%mpp_domain,isd,ied,jsd,jed)
  call mpp_get_global_domain(self%Domain%mpp_domain, isg, ieg, jsg, jeg)
  if (present(local)) then
     isc = isc - (isd-1) ; iec = iec - (isd-1) ; ied = ied - (isd-1) ; isd = 1
     jsc = jsc - (jsd-1) ; jec = jec - (jsd-1) ; jed = jed - (jsd-1) ; jsd = 1
  end if

  select case (trim(domain_type))
  case ("compute")
     is = isc; ie = iec; js = jsc; je = jec;
  case ("data")
     is = isd; ie = ied; js = jsd; je = jed;
  case ("global")
     is = isg; ie = ieg; js = jsg; je = jeg;
  end select

end subroutine soca_geom_get_domain_indices


! ------------------------------------------------------------------------------
!> Get layer depth from layer thicknesses
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_thickness2depth(self, h, z)
  class(soca_geom),     intent(in   ) :: self
  real(kind=kind_real), intent(in   ) :: h(:,:,:) !< Layer thickness
  real(kind=kind_real), intent(inout) :: z(:,:,:) !< Mid-layer depth

  integer :: is, ie, js, je, i, j, k

  ! Should check shape of z
  is = lbound(h,dim=1)
  ie = ubound(h,dim=1)
  js = lbound(h,dim=2)
  je = ubound(h,dim=2)

  ! top layer
  z(:,:,1) = 0.5_kind_real*h(:,:,1)

  ! the rest of the layers
  do i = is, ie
     do j = js, je
        do k = 2, self%nzo
          z(i,j,k) = sum(h(i,j,1:k-1))+0.5_kind_real*h(i,j,k)
        end do
     end do
  end do
end subroutine soca_geom_thickness2depth

! ------------------------------------------------------------------------------
! Get a 2d array of the valid nodes / cells that are to be used on this PE
! for ATLAS mesh generation.
! We assume that owned cells (quads) are generated north and east of the PE's owned nodes (vertices)
subroutine soca_geom_mesh_valid_nodes_cells(self, nodes, cells)
  class(soca_geom),   intent(inout) :: self
  logical, allocatable, intent(out) :: nodes(:,:), cells(:,:)

  integer :: i, j, start_i
  logical :: tripolar, cyclic

  ! TODO, do I need to worry about an extra row of cells on the bottom
  ! when on a bottom PE??
  allocate(nodes(self%isc:self%iec+1, self%jsc:self%jec+1))
  allocate(cells(self%isc:self%iec, self%jsc:self%jec))

  nodes = .true.
  cells = .true.

  tripolar = iand(self%domain%Y_FLAGS, FOLD_NORTH_EDGE) /= 0
  cyclic = iand(self%domain%X_FLAGS, CYCLIC_GLOBAL_DOMAIN) /= 0


  ! -------------------------------------------------------------------------------------
  ! remove halo points that are otherwise goint to be duplicated in the node list for the the PEs.
  
  if(tripolar .and. self%jec == self%domain%NJGLOBAL) then
    ! we are at the tripolar fold, remove some (all?) extra nodes at the northern most row.

    ! (start_i is the x index where we start removing nodes)
    ! Remove all nodes in the eastern half.
    start_i = max(self%domain%NIGLOBAL/2, self%isc)

    ! additionally, this PE crosses the E/W midpoint, remove all extra nodes at the north for this PE
    if (self%isc <= self%domain%NIGLOBAL/2 .and. self%iec > self%domain%NIGLOBAL/2) then      
      start_i = self%isc
    end if
    
    ! remove the nodes
    do i=start_i, self%iec+1
        nodes(i, self%jec+1) = .false.
    end do
  end if

  if (cyclic .and. self%isc == 1 .and. self%iec == self%domain%NIGLOBAL) then  
    ! cyclic boundaries, and this PE spans entire width.
    ! Remove all nodes on the easternmost column.
    do j=self%jsc, self%jec+1
      nodes(self%iec+1, j) = .false.
    end do
  end if

  ! -------------------------------------------------------------------------------------
  ! remove cells on the eastern half, at the top, if this has a tripolar fold
  if (tripolar .and. self%jec == self%domain%NJGLOBAL) then
    ! adjust the cells to be included
    do i=self%isc, self%iec
      if (i >= self%domain%NIGLOBAL/2) then
        cells(i, self%jec) = .false.
      end if
    end do
  end if

end subroutine

! ------------------------------------------------------------------------------

function soca_geom_atlas_idx2ij(self, idx, i, j) result(res)
  class(soca_geom), intent(in) :: self
  integer, intent(in) :: idx
  integer, intent(out) :: i ,j 
  logical :: res

  i = self%atlas_idx2i(idx)
  j = self%atlas_idx2j(idx)
  res = i > 0
end function

! ------------------------------------------------------------------------------

end module soca_geom_mod
