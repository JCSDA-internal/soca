! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_geom_mod

use atlas_module
use MOM_domains, only : MOM_domain_type, MOM_infra_init
use MOM_io,      only : io_infra_init
use soca_mom6, only: soca_mom6_config, soca_mom6_init, soca_ice_column, &
                     soca_geomdomain_init
use soca_utils, only: write2pe, soca_remap_idw
use kinds, only: kind_real
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use fms_io_mod, only : fms_io_init, fms_io_exit, &
                       register_restart_field, restart_file_type, &
                       restore_state, query_initialized, &
                       free_restart_type, save_restart
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, &
                            mpp_get_global_domain, mpp_update_domains
use fms_mod,         only : write_data
use fms_io_mod,      only : fms_io_init, fms_io_exit
use MOM_diag_remap,  only : diag_remap_ctrl, diag_remap_init, diag_remap_configure_axes, &
                            diag_remap_end, diag_remap_update
use MOM_EOS,         only : EOS_type

implicit none

private
public :: soca_geom, &
          geom_write, geom_get_domain_indices

!> Geometry data structure
type :: soca_geom
    type(MOM_domain_type), pointer :: Domain !< Ocean model domain
    type(soca_ice_column) :: ice_column !< Sea-ice geometry
    integer :: nzo, nzo_zstar, nzi, nzs, ncat
    integer :: isc, iec, jsc, jec  !< indices of compute domain
    integer :: isd, ied, jsd, jed  !< indices of data domain
    integer :: isg, ieg, jsg, jeg  !< indices of global domain
    integer :: iscl, iecl, jscl, jecl  !< indices of local compute domain
    integer :: isdl, iedl, jsdl, jedl  !< indices of local data domain
    real(kind=kind_real), allocatable, dimension(:)   :: lonh, lath
    real(kind=kind_real), allocatable, dimension(:)   :: lonq, latq
    real(kind=kind_real), allocatable, dimension(:,:) :: lon, lat !< Tracer point grid
    real(kind=kind_real), allocatable, dimension(:,:) :: lonu, latu !< Zonal velocity grid
    real(kind=kind_real), allocatable, dimension(:,:) :: lonv, latv !< Meridional velocity grid
    real(kind=kind_real), allocatable, dimension(:,:) :: sin_rot, cos_rot !< Rotation between logical grid
                                                                          !< and North
    real(kind=kind_real), allocatable, dimension(:,:) :: mask2d    !< Tracer points. 0 = land 1 = ocean
    real(kind=kind_real), allocatable, dimension(:,:) :: mask2du   !< u        "   . 0 = land 1 = ocean
    real(kind=kind_real), allocatable, dimension(:,:) :: mask2dv   !< v        "   . 0 = land 1 = ocean
    real(kind=kind_real), allocatable, dimension(:,:) :: cell_area
    real(kind=kind_real), allocatable, dimension(:,:) :: rossby_radius
    real(kind=kind_real), allocatable, dimension(:,:,:) :: h
    real(kind=kind_real), allocatable, dimension(:,:,:) :: h_zstar
    logical :: save_local_domain = .false. ! If true, save the local geometry for each pe.
    character(len=:), allocatable :: geom_grid_file
    type(fckit_mpi_comm) :: f_comm
    type(atlas_functionspace) :: afunctionspace
    type(atlas_fieldset) :: afieldset

    contains
    procedure :: init => geom_init
    procedure :: end => geom_end
    procedure :: set_atlas_lonlat => geom_set_atlas_lonlat
    procedure :: fill_atlas_fieldset => geom_fill_atlas_fieldset
    procedure :: clone => geom_clone
    procedure :: get_rossby_radius => geom_rossby_radius
    procedure :: gridgen => geom_gridgen
    procedure :: thickness2depth => geom_thickness2depth
    procedure :: struct2atlas => geom_struct2atlas
    procedure :: atlas2struct => geom_atlas2struct
    procedure :: write => geom_write
end type soca_geom

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup geometry object
subroutine geom_init(self, f_conf, f_comm)
  class(soca_geom),         intent(out) :: self
  type(fckit_configuration), intent(in) :: f_conf
  type(fckit_mpi_comm),   intent(in)    :: f_comm

  logical :: full_init = .false.

  ! MPI communicator
  self%f_comm = f_comm

  ! Domain decomposition
  call soca_geomdomain_init(self%Domain, self%nzo, f_comm)

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

  ! User-defined grid filename
  if ( .not. f_conf%get("geom_grid_file", self%geom_grid_file) ) &
     self%geom_grid_file = "soca_gridspec.nc" ! default if not found

  ! Allocate geometry arrays
  call geom_allocate(self)

  ! Check if a full initialization is required, default to false
  if ( .not. f_conf%get("full_init", full_init) ) full_init = .false.

  ! Read the geometry from file by default,
  ! skip this step if a full init is required
  if ( .not. full_init) call geom_read(self)

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
  call mpp_update_domains(self%cell_area, self%Domain%mpp_domain)

  ! Set output option for local geometry
  if ( .not. f_conf%get("save_local_domain", self%save_local_domain) ) &
     self%save_local_domain = .false.

end subroutine geom_init

! ------------------------------------------------------------------------------
!> Geometry destructor
subroutine geom_end(self)
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
  if (allocated(self%cell_area))     deallocate(self%cell_area)
  if (allocated(self%rossby_radius)) deallocate(self%rossby_radius)
  if (allocated(self%h))             deallocate(self%h)
  if (allocated(self%h_zstar))       deallocate(self%h_zstar)
  nullify(self%Domain)
  call self%afunctionspace%final()
  call self%afieldset%final()

end subroutine geom_end

! --------------------------------------------------------------------------------------------------
!> Set ATLAS lonlat fieldset
subroutine geom_set_atlas_lonlat(self, afieldset)
  class(soca_geom),  intent(inout) :: self
  type(atlas_fieldset), intent(inout) :: afieldset

  real(kind_real), pointer :: real_ptr(:,:)
  type(atlas_field) :: afield

  ! Create lon/lat field
  afield = atlas_field(name="lonlat", kind=atlas_real(kind_real), shape=(/2,(self%iec-self%isc+1)*(self%jec-self%jsc+1)/))
  call afield%data(real_ptr)
  real_ptr(1,:) = pack(self%lon(self%isc:self%iec,self%jsc:self%jec),.true.)
  real_ptr(2,:) = pack(self%lat(self%isc:self%iec,self%jsc:self%jec),.true.)
  call afieldset%add(afield)

end subroutine geom_set_atlas_lonlat

! --------------------------------------------------------------------------------------------------
!> Fill ATLAS fieldset
subroutine geom_fill_atlas_fieldset(self, afieldset)
  class(soca_geom),  intent(inout) :: self
  type(atlas_fieldset), intent(inout) :: afieldset

  integer :: i, jz, n
  integer, pointer :: int_ptr_2(:,:)
  real(kind=kind_real), pointer :: real_ptr_1(:), real_ptr_2(:,:)
  type(atlas_field) :: afield

  ! Add area
  afield = self%afunctionspace%create_field(name='area', kind=atlas_real(kind_real), levels=0)
  call afield%data(real_ptr_1)
  real_ptr_1 = pack(self%cell_area(self%isc:self%iec,self%jsc:self%jec),.true.)
  call afieldset%add(afield)
  call afield%final()

  ! Add vertical unit
  afield = self%afunctionspace%create_field(name='vunit', kind=atlas_real(kind_real), levels=self%nzo)
  call afield%data(real_ptr_2)
  do jz=1,self%nzo
    real_ptr_2(jz,:) = real(jz, kind_real)
  end do
  call afieldset%add(afield)
  call afield%final()

  ! Add geographical mask
  afield = self%afunctionspace%create_field(name='gmask', kind=atlas_integer(kind(0)), levels=self%nzo)
  call afield%data(int_ptr_2)
  do jz=1,self%nzo
    int_ptr_2(jz,:) = int(pack(self%mask2d(self%isc:self%iec,self%jsc:self%jec),.true.))
  end do
  call afieldset%add(afield)
  call afield%final()

end subroutine geom_fill_atlas_fieldset

! ------------------------------------------------------------------------------
!> Clone, self = other
subroutine geom_clone(self, other)
  class(soca_geom), intent( in) :: self
  class(soca_geom), intent(out) :: other

  ! Clone communicator
  other%f_comm = self%f_comm

  ! Clone fms domain and vertical levels
  other%Domain => self%Domain
  other%nzo = self%nzo

  !
  other%geom_grid_file = self%geom_grid_file

  ! Clone sea-ice grid
  other%ice_column = self%ice_column
  other%ncat = self%ncat
  other%nzi = self%nzi
  other%nzs = self%nzs

  ! Allocate and clone geometry
  call geom_allocate(other)
  other%lonh = self%lonh
  other%lath = self%lath
  other%lonq = self%lonq
  other%latq = self%latq
  other%lon = self%lon
  other%lat = self%lat
  other%lonu = self%lonu
  other%latu = self%latu
  other%lonv = self%lonv
  other%latv = self%latv
  other%sin_rot = self%sin_rot
  other%cos_rot = self%cos_rot
  other%mask2d = self%mask2d
  other%mask2du = self%mask2du
  other%mask2dv = self%mask2dv
  other%cell_area = self%cell_area
  other%rossby_radius = self%rossby_radius
  other%h = self%h

end subroutine geom_clone

! ------------------------------------------------------------------------------
!>
subroutine geom_gridgen(self)
  class(soca_geom), intent(inout) :: self

  ! allocate variables for regridding to zstar coord
  type(soca_mom6_config) :: mom6_config
  type(diag_remap_ctrl) :: remap_ctrl
  type(EOS_type), pointer :: eqn_of_state
  integer :: k
  real(kind=kind_real), allocatable :: T(:,:,:), S(:,:,:)

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
  self%cell_area  = mom6_config%grid%areaT
  self%h = mom6_config%MOM_CSp%h

  !Generate new grid based on zstar coordinate
  allocate(T(self%isd:self%ied, self%jsd:self%jed, self%nzo))
  allocate(S(self%isd:self%ied, self%jsd:self%jed, self%nzo))
  T = 0.d0 !fake temp
  S = 0.d0 !fake salinity
  call diag_remap_init(remap_ctrl, coord_tuple='ZSTAR, ZSTAR, ZSTAR')
  call diag_remap_configure_axes(remap_ctrl, mom6_config%GV, mom6_config%scaling, mom6_config%param_file)
  call diag_remap_update(remap_ctrl, mom6_config%grid, mom6_config%GV, mom6_config%scaling, self%h, T, S, eqn_of_state)
  self%nzo_zstar = remap_ctrl%nz
  if (allocated(self%h_zstar)) deallocate(self%h_zstar)
  allocate(self%h_zstar(self%isd:self%ied, self%jsd:self%jed, 1:remap_ctrl%nz))
  self%h_zstar = remap_ctrl%h
  call diag_remap_end(remap_ctrl)

  ! Get Rossby Radius
  call geom_rossby_radius(self)

  ! Output to file
  call geom_write(self)

end subroutine geom_gridgen

! ------------------------------------------------------------------------------
!> Allocate memory and point to mom6 data structure
subroutine geom_allocate(self)
  class(soca_geom), intent(inout) :: self

  integer :: nzo
  integer :: isd, ied, jsd, jed

  ! Get domain shape (number of levels, indices of data and compute domain)
  call geom_get_domain_indices(self, "compute", self%isc, self%iec, self%jsc, self%jec)
  call geom_get_domain_indices(self, "data", isd, ied, jsd, jed)
  self%isd = isd ;  self%ied = ied ; self%jsd = jsd; self%jed = jed
  call geom_get_domain_indices(self, "global", self%isg, self%ieg, self%jsg, self%jeg)
  call geom_get_domain_indices(self, "compute", self%iscl, self%iecl, self%jscl, self%jecl, local=.true.)
  call geom_get_domain_indices(self, "data", self%isdl, self%iedl, self%jsdl, self%jedl, local=.true.)
  nzo = self%nzo

  ! Allocate arrays on compute domain
  allocate(self%lonh(self%isg:self%ieg));        self%lonh = 0.0_kind_real
  allocate(self%lath(self%jsg:self%jeg));        self%lath = 0.0_kind_real
  allocate(self%lonq(self%isg:self%ieg));        self%lonq = 0.0_kind_real
  allocate(self%latq(self%jsg:self%jeg));        self%latq = 0.0_kind_real
  allocate(self%lon(isd:ied,jsd:jed));           self%lon = 0.0_kind_real
  allocate(self%lat(isd:ied,jsd:jed));           self%lat = 0.0_kind_real
  allocate(self%lonu(isd:ied,jsd:jed));          self%lonu = 0.0_kind_real
  allocate(self%latu(isd:ied,jsd:jed));          self%latu = 0.0_kind_real
  allocate(self%lonv(isd:ied,jsd:jed));          self%lonv = 0.0_kind_real
  allocate(self%latv(isd:ied,jsd:jed));          self%latv = 0.0_kind_real

  allocate(self%sin_rot(isd:ied,jsd:jed));       self%sin_rot = 0.0_kind_real
  allocate(self%cos_rot(isd:ied,jsd:jed));       self%cos_rot = 0.0_kind_real

  allocate(self%mask2d(isd:ied,jsd:jed));        self%mask2d = 0.0_kind_real
  allocate(self%mask2du(isd:ied,jsd:jed));       self%mask2du = 0.0_kind_real
  allocate(self%mask2dv(isd:ied,jsd:jed));       self%mask2dv = 0.0_kind_real

  allocate(self%cell_area(isd:ied,jsd:jed));     self%cell_area = 0.0_kind_real
  allocate(self%rossby_radius(isd:ied,jsd:jed)); self%rossby_radius = 0.0_kind_real
  allocate(self%h(isd:ied,jsd:jed,1:nzo));       self%h = 0.0_kind_real

end subroutine geom_allocate

! ------------------------------------------------------------------------------
!> Read and store Rossby Radius of deformation
subroutine geom_rossby_radius(self)
  class(soca_geom), intent(inout) :: self

  integer :: unit, i, n
  real(kind=kind_real) :: dum
  real(kind=kind_real), allocatable :: lon(:),lat(:),rr(:)
  integer :: isc, iec, jsc, jec
  integer :: io

  ! read in the file
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

  ! convert to meters
  rr = rr * 1e3

  ! remap
  isc = self%isc ;  iec = self%iec ; jsc = self%jsc ; jec = self%jec
  call soca_remap_idw(lon, lat, rr, self%lon(isc:iec,jsc:jec), &
                      self%lat(isc:iec,jsc:jec), self%rossby_radius(isc:iec,jsc:jec) )

end subroutine geom_rossby_radius


! ------------------------------------------------------------------------------
!> Write geometry to file
subroutine geom_write(self)
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
  call save_restart(geom_restart, directory='')
  call free_restart_type(geom_restart)
  call fms_io_exit()

  if (self%save_local_domain) then
     ! Save local compute grid
     pe = self%f_comm%rank()

     write (strpe,fmt) pe
     geom_output_pe='geom_output_'//trim(strpe)//'.nc'

     ns = (self%iec - self%isc + 1) * (self%jec - self%jsc + 1 )
     call write2pe(reshape(self%mask2d,(/ns/)),'mask',geom_output_pe,.false.)
     call write2pe(reshape(self%lon,(/ns/)),'lon',geom_output_pe,.true.)
     call write2pe(reshape(self%lat,(/ns/)),'lat',geom_output_pe,.true.)
  end if

end subroutine geom_write

! ------------------------------------------------------------------------------
!> Read geometry from file
subroutine geom_read(self)
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
!> Copy a structured field into an ATLAS fieldset
subroutine geom_struct2atlas(self, dx_struct, dx_atlas)
  class(soca_geom),     intent(in ) :: self
  real(kind=kind_real), intent(in ) :: dx_struct(:,:)
  type(atlas_fieldset), intent(out) :: dx_atlas

  real(kind_real), pointer :: real_ptr(:)
  type(atlas_field) :: afield

  dx_atlas = atlas_fieldset()
  afield = self%afunctionspace%create_field('var_00',kind=atlas_real(kind_real),levels=0)
  call dx_atlas%add(afield)
  call afield%data(real_ptr)
  real_ptr = pack(dx_struct(self%iscl:self%iecl, self%jscl:self%jecl),.true.)
  call afield%final()

end subroutine geom_struct2atlas

! ------------------------------------------------------------------------------
!> Copy a structured field from an ATLAS fieldset
subroutine geom_atlas2struct(self, dx_struct, dx_atlas)
  class(soca_geom),     intent(in   ) :: self
  real(kind=kind_real), intent(inout) :: dx_struct(:,:)
  type(atlas_fieldset), intent(inout) :: dx_atlas

  real(kind_real), pointer :: real_ptr(:)
  logical :: umask(self%iscl:self%iecl,self%jscl:self%jecl)
  type(atlas_field) :: afield

  umask = .true.
  afield = dx_atlas%field('var_00')
  call afield%data(real_ptr)
  dx_struct(self%iscl:self%iecl, self%jscl:self%jecl) = unpack(real_ptr,umask,dx_struct(self%iscl:self%iecl, self%jscl:self%jecl))
  call afield%final()

end subroutine geom_atlas2struct

! ------------------------------------------------------------------------------

end module soca_geom_mod
