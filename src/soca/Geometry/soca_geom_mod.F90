! (C) Copyright 2017-2024 UCAR
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
use fms_mod, only : fms_init, fms_end
use fms_io_mod, only : fms_io_init, fms_io_exit, &
                       register_restart_field, restart_file_type, &
                       restore_state, free_restart_type, save_restart
use MOM, only : MOM_control_struct, initialize_MOM, MOM_end, get_MOM_state_elements
use MOM_restart, only :MOM_restart_CS ! NOTE remove this when updating MOM6
use MOM_domains, only : MOM_domain_type, MOM_domains_init, MOM_infra_init, MOM_infra_end
use MOM_error_handler, only : MOM_error, MOM_mesg, WARNING, FATAL, is_root_pe
use MOM_file_parser, only : get_param, param_file_type, close_param_file
use MOM_get_input, only : directories, Get_MOM_Input
use MOM_grid, only : ocean_grid_type
use MOM_io, only : io_infra_init, io_infra_end
use MOM_time_manager, only : real_to_time, JULIAN, set_calendar_type
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, &
                            mpp_get_global_domain, mpp_update_domains, &
                            CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
use mpp_mod,only : mpp_init
use time_interp_external_mod, only : time_interp_external_init
use time_manager_mod, only: time_type

! soca modules
use soca_fields_metadata_mod, only : soca_fields_metadata
use soca_utils, only: write2pe

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

    !> \}

    !> instance of the metadata that is read in from a config file upon initialization
    type(soca_fields_metadata) :: fields_metadata

    type(fckit_mpi_comm) :: f_comm !< MPI communicator

    !> mesh parameters
    type(atlas_functionspace_NodeColumns) :: functionspace
    integer, allocatable :: atlas_ij2idx(:,:)
    type(atlas_fieldset) :: fieldset !< the geom fields (area, mask, etc)

  contains


    !> \copybrief soca_geom_init \see soca_geom_init
    procedure :: init => soca_geom_init

    !> \copybrief soca_geom_end \see soca_geom_end
    procedure :: end => soca_geom_end

    !> \copybrief soca_geom_fill_atlas_fieldset \see soca_geom_fill_atlas_fieldset
    procedure :: init_fieldset => soca_geom_init_fieldset

    !> \copybrief soca_geom_write \see soca_geom_write
    procedure :: write => soca_geom_write

    !> \copybrief soca_geom_mesh_valid_nodes_cells \see soca_geom_mesh_valid_nodes_cells
    procedure :: mesh_valid_nodes_cells => soca_geom_mesh_valid_nodes_cells
end type soca_geom

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup geometry object
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_init(self, f_conf, f_comm, gen)
  class(soca_geom),         intent(out) :: self
  type(fckit_configuration), intent(in) :: f_conf
  type(fckit_mpi_comm),   intent(in)    :: f_comm !< MPI communicator for this geometry
  logical,                  intent(in)  :: gen !< if true, we are doing a full init

  ! variables needed by MOM6
  type(param_file_type) :: param_file      ! The structure indicating the file(s)
  type(directories)  :: dirs       !< Relevant dirs/path

  ! variables needed to get grid info from MOM6
  type(time_type) :: Start_time
  type(ocean_grid_type),    pointer :: grid !< Grid metrics
  type(MOM_control_struct)  :: CSp

  ! variables needed for reading gridspec file
  integer :: r
  type(restart_file_type) :: geom_restart
  character(len=:), allocatable :: str

  ! MPI communicator
  self%f_comm = f_comm

  ! use MOM6 to setup domain decomposition
  call mpp_init(localcomm=f_comm%communicator())
  call fms_init()
  call fms_io_init()
  call Get_MOM_Input(param_file, dirs)
  call MOM_domains_init(self%Domain, param_file)
  call get_param(param_file, "soca_mom6", "NK", self%nzo, fail_if_missing=.true.)
  call close_param_file(param_file)

  ! Allocate geometry arrays
  call soca_geom_allocate(self)

  if (gen) then
    ! generate the grid from MOM6
    Start_time = real_to_time(0.0d0)
    call MOM_infra_init(localcomm=self%f_comm%communicator())
    call io_infra_init()
    call set_calendar_type(JULIAN)
    call time_interp_external_init
    call initialize_MOM( Start_time, Start_time, param_file, dirs, CSp )
    call get_MOM_state_elements(CSp, G=grid)

    ! grab the grid properties that we need
    self%lonh = grid%gridlont
    self%lath = grid%gridlatt
    self%lonq = grid%gridlonb
    self%latq = grid%gridlatb
    self%lon =  grid%GeoLonT
    self%lat =  grid%GeoLatT
    self%lonu = grid%geoLonCu
    self%latu = grid%geoLatCu
    self%lonv = grid%geoLonCv
    self%latv = grid%geoLatCv
    self%sin_rot = grid%sin_rot
    self%cos_rot = grid%cos_rot
    self%mask2d = grid%mask2dT
    self%mask2du = grid%mask2dCu
    self%mask2dv = grid%mask2dCv
    self%dx = grid%dxT
    self%dy = grid%dyT
    self%cell_area = grid%areaT

    call MOM_end(CSp)
    call io_infra_end()

  else
    ! Read in the precomputed grid from the soca gridspec file instead
    ! NOTE that we will rerad the gridspec file later for some of the variables
    ! once the altas FunctionSpace has been created.
    if (.not. f_conf%get("geom_grid_file", str)) str = "soca_gridspec.nc"
    r = register_restart_field(geom_restart, str, "lonh",    self%lonh,      self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "lath",    self%lath,      self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "lonq",    self%lonq,      self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "latq",    self%latq,      self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "lon",     self%lon,       self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "lat",     self%lat,       self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "lonu",    self%lonu,      self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "latu",    self%latu,      self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "lonv",    self%lonv,      self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "latv",    self%latv,      self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "sin_rot", self%sin_rot,   self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "cos_rot", self%cos_rot,   self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "dx",      self%dx,        self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "dy",      self%dy,        self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "area",    self%cell_area, self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "mask2d",  self%mask2d,    self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "mask2du", self%mask2du,   self%Domain%mpp_domain)
    r = register_restart_field(geom_restart, str, "mask2dv", self%mask2dv,   self%Domain%mpp_domain)
    call restore_state(geom_restart, directory='')
    call free_restart_type(geom_restart)
  endif

  ! Fill halos
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

  ! process the fields metadata file
  call f_conf%get_or_die("fields metadata", str)
  call self%fields_metadata%create(str)

end subroutine soca_geom_init

! ------------------------------------------------------------------------------
!> Geometry destructor
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_end(self)
  class(soca_geom), intent(out)  :: self

  call fms_io_exit()
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
  nullify(self%Domain)
  call self%functionspace%final()

end subroutine soca_geom_end


! --------------------------------------------------------------------------------------------------
!> Fill ATLAS fieldset
!! NOTE, this is only the set of variables being filled in from the Fortran side,
!! some variables are being set from the C++ side as well.
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_init_fieldset(self, f_conf, gen)
  class(soca_geom),  intent(inout) :: self
  type(fckit_configuration), intent(in) :: f_conf
  logical,                  intent(in)  :: gen !< if true, we are doing a full init

  integer :: i, j, n, jz
  type(atlas_field) :: fArea, fInterpMask, fVertCoord, fGmask, fOwned, fMaskH, fMaskU, fMaskV
  real(kind=kind_real), pointer :: vArea(:,:), vInterpMask(:,:), vVertCoord(:,:), &
                                   vMaskH(:,:), vMaskU(:,:), vMaskV(:,:)
  integer, pointer :: vGmask(:,:), vOwned(:,:)

  ! variables needed for reading in atlas fields from gridspec file
  integer :: v, r
  character(len=:), allocatable :: str
  type(restart_file_type) :: geom_restart
  character(len=20), dimension(2) :: atlasVars
  type(atlas_field) :: aField
  real(kind=kind_real), pointer :: aFieldData(:,:)
  real(kind=kind_real), allocatable :: fieldData(:,:), fieldDataVars(:,:,:)


  ! some fields always have to be generated because they are dependant on the PE layout
  ! -----------------------------------------------------------------------------------------------
  fOwned = self%functionspace%create_field(name='owned', kind=atlas_integer(kind(0)), levels=1)
  call self%fieldset%add(fOwned)
  call fOwned%data(vOwned)
  vOwned = 0 ! need to set to 0, it's the only one that doesn't get halo update

  ! These fields could probably be saved in the gridspec file (I'll clean that up another day!)
  ! For now though they are regenerated
  ! -----------------------------------------------------------------------------------------------
  fArea = self%functionspace%create_field(name='area', kind=atlas_real(kind_real), levels=1)
  call self%fieldset%add(fArea)
  call fArea%data(vArea)

  fInterpMask = self%functionspace%create_field(name='interp_mask', kind=atlas_real(kind_real), levels=1)
  call self%fieldset%add(fInterpMask)
  call fInterpMask%data(vInterpMask)

  fVertCoord = self%functionspace%create_field(name='vert_coord', kind=atlas_real(kind_real), levels=self%nzo)
  call self%fieldset%add(fVertCoord)
  call fVertCoord%data(vVertCoord)

  fGmask = self%functionspace%create_field(name='gmask', kind=atlas_integer(kind(0)), levels=self%nzo)
  call self%fieldset%add(fGmask)
  call fGmask%data(vGmask)

  ! TODO temporary, this needs to be cleaned up
  ! It's here so that answers don't change when moving the gpnorm calculations from fortran
  ! to c++
  fMaskH = self%functionspace%create_field(name='mask_h', kind=atlas_real(kind_real), levels=1)
  call self%fieldset%add(fMaskH)
  call fMaskH%data(vMaskH)

  fMaskU = self%functionspace%create_field(name='mask_u', kind=atlas_real(kind_real), levels=1)
  call self%fieldset%add(fMaskU)
  call fMaskU%data(vMaskU)

  fMaskV = self%functionspace%create_field(name='mask_v', kind=atlas_real(kind_real), levels=1)
  call self%fieldset%add(fMaskV)
  call fMaskV%data(vMaskV)

  ! set the data
  do j=self%jsc,self%jec
    do i=self%isc,self%iec
      n = self%atlas_ij2idx(i,j)
      vArea(1,n) = self%cell_area(i,j)
      vInterpMask(1,n) = self%mask2d(i,j)
      vGmask(:, n) = int(self%mask2d(i,j))
      vOwned(1, n) = 1
      vMaskH(1, n) = self%mask2d(i,j)
      vMaskU(1, n) = self%mask2du(i,j)
      vMaskV(1, n) = self%mask2dv(i,j)
    end do
  end do
  do jz=1,self%nzo
    vVertCoord(jz,:) = real(jz, kind_real)
  end do

  ! halo exchanges for some of the fields
  call fGmask%halo_exchange()
  call fInterpMask%halo_exchange()
  call fArea%halo_exchange()
  call fVertCoord%halo_exchange()

  ! done, cleanup
  call fArea%final()
  call fInterpMask%final()
  call fVertCoord%final()
  call fGmask%final()
  call fOwned%final()
  call fMaskH%final()
  call fMaskU%final()
  call fMaskV%final()

  ! -----------------------------------------------------------------------------------------------
  if (gen) then
    ! some fields are still being generated the fortran side (someday this will be moved to c++)
    call soca_geom_distance_from_coast(self, fieldData)
    aField = self%functionspace%create_field(name='distance_from_coast', kind=atlas_real(kind_real), levels=1)
    call self%fieldset%add(aField)
    call aField%data(aFieldData)
    do j=self%jsc,self%jec
      do i=self%isc,self%iec
        aFieldData(1, self%atlas_ij2idx(i,j)) = fieldData(i, j)
      end do
    end do
    call aField%final()
  else
    ! otherwise, read in files from the gridspec file and place directly into atlas fieldset
    atlasVars = [character(len=20) :: "rossby_radius", "distance_from_coast"]
    if (.not. f_conf%get("geom_grid_file", str)) str = "soca_gridspec.nc"
    allocate(fieldDataVars(self%isd:self%ied, self%jsd:self%jed, size(atlasVars)))

    ! read in from gridspec file
    do v = 1, size(atlasVars)
      r = register_restart_field(geom_restart, str, atlasVars(v), fieldDataVars(:,:,v), self%Domain%mpp_domain)
    end do
    call restore_state(geom_restart, directory='')
    call free_restart_type(geom_restart)

    ! copy from fortran array to atlas field
    do v = 1, size(atlasVars)
      aField = self%functionspace%create_field(atlasVars(v), atlas_real(kind_real), 1)
      call self%fieldset%add(aField)
      call aField%data(aFieldData)
      do j=self%jsc,self%jec
        do i=self%isc,self%iec
          aFieldData(1, self%atlas_ij2idx(i,j)) = fieldDataVars(i,j,v )
        end do
      end do
      call aField%halo_exchange()
    end do
    call aField%final()
  endif
end subroutine soca_geom_init_fieldset


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

  allocate(self%atlas_ij2idx(isd:ied,jsd:jed));  self%atlas_ij2idx = -1

end subroutine soca_geom_allocate


! ------------------------------------------------------------------------------
!> Calcuate distance from coast for the ocean points
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_distance_from_coast(self, dist)
  class(soca_geom), intent(inout) :: self
  real(kind=kind_real), allocatable :: dist(:,:)

  type(atlas_indexkdtree) :: kd
  type(atlas_geometry) :: ageometry
  integer :: i, j, idx(1)
  integer :: num_land_l, num_land
  integer, allocatable :: rcvcnt(:), displs(:)
  real(kind=kind_real), allocatable :: land_lon(:), land_lat(:)
  real(kind=kind_real), allocatable :: land_lon_l(:), land_lat_l(:)
  real(kind=kind_real) :: closest_lon, closest_lat

  allocate(dist(self%isd:self%ied, self%jsd:self%jed))
  dist = 0.0

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
      dist(i,j) = ageometry%distance( &
            self%lon(i,j), self%lat(i,j), land_lon(idx(1)), land_lat(idx(1)))
    enddo
  enddo

  ! cleanup
  call kd%final()

end subroutine

! ------------------------------------------------------------------------------
!> Write geometry to file
!!
!! \related soca_geom_mod::soca_geom
subroutine soca_geom_write(self, f_conf)
  class(soca_geom), intent(in) :: self
  type(fckit_configuration), intent(in) :: f_conf

  ! variables needed for writing restart file
  character(len=:), allocatable :: str
  integer :: r
  type(restart_file_type) :: geom_restart

  ! variables needed for writing atlas fields
  integer :: i, j, v
  character(len=20), dimension(2) :: atlasVars
  type(atlas_field) :: aField
  real(kind=kind_real), pointer :: aFieldData(:,:)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: fieldData

  ! variables for writing local domain
  logical :: save_local
  integer :: ns
  character(len=256) :: geom_output_pe

  ! Save global domain
  if (.not. f_conf%get("geom_grid_file", str)) str = "soca_gridspec.nc"
  r = register_restart_field(geom_restart, str, "lonh",    self%lonh,      self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "lath",    self%lath,      self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "lonq",    self%lonq,      self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "latq",    self%latq,      self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "lon",     self%lon,       self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "lat",     self%lat,       self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "lonu",    self%lonu,      self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "latu",    self%latu,      self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "lonv",    self%lonv,      self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "latv",    self%latv,      self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "sin_rot", self%sin_rot,   self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "cos_rot", self%cos_rot,   self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "dx",      self%dx,        self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "dy",      self%dy,        self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "area",    self%cell_area, self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "mask2d",  self%mask2d,    self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "mask2du", self%mask2du,   self%Domain%mpp_domain)
  r = register_restart_field(geom_restart, str, "mask2dv", self%mask2dv,   self%Domain%mpp_domain)

  ! write out a subset of the fields in the atlas fieldset
  ! (note, someday *most* of the var listed above will be moved to atlas fieldsets)
  atlasVars = [character(len=20) :: "rossby_radius", "distance_from_coast"]
  allocate(fieldData(self%isd:self%ied, self%jsd:self%jed, size(atlasVars)))
  do v = 1, size(atlasVars)
    ! copy from atlas, to our fortran array
    aField = self%fieldset%field(trim(atlasVars(v)))
    call aField%data(aFieldData)
    do j=self%jsc,self%jec
      do i=self%isc,self%iec
        fieldData(i,j,v ) = aFieldData(1, self%atlas_ij2idx(i,j))
      end do
    end do
    call aField%final()

    ! register with FMS
    r = register_restart_field(geom_restart, str, trim(atlasVars(v)), fieldData(:,:, v), self%Domain%mpp_domain)
  end do

  call save_restart(geom_restart, directory='')
  call free_restart_type(geom_restart)

  ! Set output option for local geometry
  if ( .not. f_conf%get("save_local_domain", save_local) ) save_local = .false.
  if (save_local) then
    ! Save local compute grid
    write (geom_output_pe, "(A,I5.5,A)") "geom_output_", self%f_comm%rank(), ".nc"

    ns = (self%iec - self%isc + 1) * (self%jec - self%jsc + 1 )
    call write2pe(reshape(self%mask2d(self%isc:self%iec,self%jsc:self%jec),(/ns/)),'mask',geom_output_pe,.false.)
    call write2pe(reshape(self%lon(self%isc:self%iec,self%jsc:self%jec),(/ns/)),'lon',geom_output_pe,.true.)
    call write2pe(reshape(self%lat(self%isc:self%iec,self%jsc:self%jec),(/ns/)),'lat',geom_output_pe,.true.)
  end if

end subroutine soca_geom_write


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
! Get a 2d array of the valid nodes / cells that are to be used on this PE
! for ATLAS mesh generation.
! We assume that owned cells (quads) are generated north and east of the PE's owned nodes (vertices)
! NOTE: this code assumes that a regional domain does NOT make use of halos on the boundary.
!  We currently have our regional domain surrounded by "virtual land", so this is currently a correct
!  assumption. When we get around to removing that "virtual land", and intend to populate the
!  halo with a boundary conditions, the grid will need to be constructed slightly differently.
subroutine soca_geom_mesh_valid_nodes_cells(self, nodes, cells)
  class(soca_geom),   intent(inout) :: self
  logical, allocatable, intent(out) :: nodes(:,:), cells(:,:)

  integer :: i, j, start_i
  logical :: tripolar, cyclic, regional

  allocate(nodes(self%isc:self%iec+1, self%jsc:self%jec+1))
  allocate(cells(self%isc:self%iec, self%jsc:self%jec))

  nodes = .true.
  cells = .true.

  tripolar = iand(self%domain%Y_FLAGS, FOLD_NORTH_EDGE) /= 0
  cyclic = iand(self%domain%X_FLAGS, CYCLIC_GLOBAL_DOMAIN) /= 0
  regional = iand(self%domain%X_FLAGS, CYCLIC_GLOBAL_DOMAIN) == 0 .and. &
             iand(self%domain%Y_FLAGS, CYCLIC_GLOBAL_DOMAIN) == 0

  ! -------------------------------------------------------------------------------------
  ! we need to mask out cells and/or nodes for 3 special cases.
  ! Remove halo points that are otherwise going to be duplicated in the node list for the the PEs,
  ! or, for regional (for now) remove unowned halo points at the domain boundary.
  ! -------------------------------------------------------------------------------------

  ! The grid is tripolar, and we are a PE at the northern edge
  ! -------------------------------------------------------------------------------------
  if(tripolar .and. self%jec == self%domain%NJGLOBAL) then
      ! (start_i is the x index where we start removing nodes)
    ! Remove all nodes in the eastern half.
    start_i = max(self%domain%NIGLOBAL/2, self%isc)

    ! additionally, this PE crosses the E/W midpoint, remove all extra nodes at the north for this PE
    if (self%isc <= self%domain%NIGLOBAL/2 .and. self%iec > self%domain%NIGLOBAL/2) then
      start_i = self%isc
    end if

    ! remove some nodes
    do i=start_i, self%iec+1
        nodes(i, self%jec+1) = .false.
    end do

    ! remove some cells
    do i=self%isc, self%iec
      if (i >= self%domain%NIGLOBAL/2) then
        cells(i, self%jec) = .false.
      end if
    end do
  end if

  ! the grid is cyclic in the E/W direction and we are a PE at the eastern edge
  ! -------------------------------------------------------------------------------------
  if (cyclic .and. self%isc == 1 .and. self%iec == self%domain%NIGLOBAL) then
    ! cyclic boundaries, and this PE spans entire width.
    ! Remove all nodes on the easternmost column.
    do j=self%jsc, self%jec+1
      nodes(self%iec+1, j) = .false.
    end do
  end if

  ! We are a regional grid
  ! NOTE: this logic should change at some point so that there ARE halos around the
  ! outer boundary
  ! -------------------------------------------------------------------------------------
  if (regional) then
    ! are we are a PE at the top edge?
    if(self%jec == self%domain%NJGLOBAL) then
      ! remove nodes on top
      do i=self%isc,self%iec+1
        nodes(i, self%jec+1) = .false.
      end do
      ! remove cells on top
      do i=self%isc, self%iec
        cells(i, self%jec) = .false.
      end do
    end if

    ! are we a PE at the right edge?
    if(self%iec == self%domain%NIGLOBAL) then
      ! remove nodes on right
      do j=self%jsc,self%jec+1
        nodes(self%iec+1, j) = .false.
      end do
      ! remove cells on right
      do j=self%jsc, self%jec
        cells(self%iec, j) = .false.
      end do
    end if

  end if

end subroutine

! ------------------------------------------------------------------------------

end module soca_geom_mod
