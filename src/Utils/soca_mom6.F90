!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! 
! ------------------------------------------------------------------------------
!> Typical geometry/state for ocean/sea-ice models
!>
!>   level 0          ---- Surface                     Tsfc
!>   level 1          ---- Upper snow level            Qsno
!>     .                        .                       .
!>     .                        .                       .
!>     .                        .                       .        
!>   level nzs        ---- Lower snow level            Qsno
!>   level nzs+1      ---- Upper ice level             Qice     
!>     .                        .                       .
!>     .                        .                       .
!>     .                        .                       .
!>   level nzs+nzi    ---- Upper ice level             Qice
!>   level nzs+nzi+1  ---- Ice/Ocean interface level   Tfreeze (from S)
!>   level nzs+nzi+2  ---- Upper level Ocean           SST    
!>
!> 
module soca_mom6

  use ocean_model_mod,         only: ocean_state_type, ocean_public_type
  use time_manager_mod,        only: time_type
  use kinds
  use MOM_grid,                  only : ocean_grid_type
  use MOM_verticalGrid,          only : verticalGrid_type
  use MOM, only : MOM_control_struct
  
  implicit none
  private

  public :: soca_geom_init, soca_geom_end, Coupled, soca_field_init, soca_field_end, soca_ice_column

  type soca_ice_column 
     integer                        :: ncat ! Number of ice categories
     integer                        :: nzi  ! Number of ice levels
     integer                        :: nzs  ! Number of snow levels     
  end type soca_ice_column
  
  type soca_ocn_data_type
     real(kind=kind_real),  pointer :: T(:,:,:)
     real(kind=kind_real),  pointer :: S(:,:,:)
     real(kind=kind_real),  pointer :: U(:,:,:)
     real(kind=kind_real),  pointer :: V(:,:,:)
     real(kind=kind_real),  pointer :: ssh(:,:)
     real(kind=kind_real),  pointer :: H(:,:,:)
  end type soca_ocn_Data_Type

  type soca_ice_data_type 
     type(soca_ice_column) :: ice_column
     real(kind=kind_real),  pointer :: part_size(:,:,:)
     real(kind=kind_real),  pointer :: h_ice(:,:,:)
     real(kind=kind_real),  pointer :: enth_ice(:,:,:,:)
     real(kind=kind_real),  pointer :: sal_ice(:,:,:,:)
     real(kind=kind_real),  pointer :: h_snow(:,:,:)
     real(kind=kind_real),  pointer :: enth_snow(:,:,:,:)
     real(kind=kind_real),  pointer :: T_skin(:,:,:)
  end type soca_ice_data_type
  
  type Coupled
     type (soca_ice_data_type)     :: Ice     
     type (soca_ocn_data_type)     :: Ocn
  end type Coupled

contains

  ! ------------------------------------------------------------------------------
  
  subroutine soca_geom_init(G, GV, ice_column, c_conf)

    use MOM_dyn_horgrid,           only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
    use MOM_transcribe_grid,       only : copy_dyngrid_to_MOM_grid
    use MOM_verticalGrid,          only : verticalGrid_type, verticalGridInit, verticalGridEnd
    use MOM_grid,                  only : MOM_grid_init, ocean_grid_type
    use MOM_file_parser,           only : get_param, param_file_type
    use MOM_get_input,             only : Get_MOM_Input, directories
    use MOM_domains,               only : MOM_domains_init, clone_MOM_domain
    use MOM_hor_index,             only : hor_index_type, hor_index_init
    use fms_io_mod,                only : fms_io_init, fms_io_exit
    use fms_mod,                   only : read_data, write_data
    use MOM_domains,               only : MOM_infra_init, MOM_infra_end
    use MOM_io,                    only : check_nml_error, io_infra_init, io_infra_end
    use MOM_file_parser,           only : open_param_file, close_param_file
    use MOM_fixed_initialization,  only : MOM_initialize_fixed
    use MOM_open_boundary,         only : ocean_OBC_type
    use kinds
    use iso_c_binding
    use config_mod
    
    implicit none

    type(ocean_grid_type),            intent(out) :: G          !< The horizontal grid type (same for ice & ocean)
    type(verticalGrid_type), pointer, intent(out) :: GV         !< Ocean vertical grid
    type(soca_ice_column),            intent(out) :: ice_column !< Ice grid spec
    type(c_ptr),                       intent(in) :: c_conf

    type(dyn_horgrid_type),       pointer  :: dG => NULL()              !< Dynamic grid
    type(hor_index_type)                   :: HI                        !< Horiz index array extents
    type(param_file_type)                  :: param_file                !< Structure to parse for run-time parameters
    type(directories)                      :: dirs                      !< Structure containing several relevant directory paths
    logical                                :: global_indexing = .false. !< If true use global horizontal index DOES NOT WORK
    logical                                :: write_geom_files = .false.!< 
    type(ocean_OBC_type),          pointer :: OBC => NULL()             !< Ocean boundary condition
    integer :: NCat_dflt = 5
    integer :: ncat, nkice, nksnow, km

    ! Initialize fms io
    call fms_io_init()

    ! Parse grid inputs
    call Get_MOM_Input(param_file, dirs)

    ! Domain decomposition/Inintialize mpp domains
    call MOM_domains_init(G%domain, param_file)
    call hor_index_init(G%Domain, HI, param_file, local_indexing=.not.global_indexing)
    call create_dyn_horgrid(dG, HI)
    call clone_MOM_domain(G%Domain, dG%Domain)

    ! Allocate grid arrays
    GV => NULL()
    call verticalGridInit( param_file, GV )
    call MOM_grid_init(G, param_file, HI, global_indexing=global_indexing)

    ! Read/Generate grid
    call MOM_initialize_fixed(dG, OBC, param_file, write_geom_files, dirs%output_directory)
    call copy_dyngrid_to_MOM_grid(dG, G)
    G%ke = GV%ke ; G%g_Earth = GV%g_Earth

    ! Destructor for dynamic grid
    call destroy_dyn_horgrid(dG)
    dG => NULL()
    
    ! Initialize sea-ice grid
    ice_column%ncat = config_get_int(c_conf,"num_ice_cat")
    ice_column%nzi = config_get_int(c_conf,"num_ice_lev")
    ice_column%nzs = config_get_int(c_conf,"num_sno_lev")
    
    call close_param_file(param_file)
    
    call fms_io_exit()

  end subroutine soca_geom_init

  ! ------------------------------------------------------------------------------
  
  subroutine soca_geom_end(G, GV)

    use MOM_grid,                  only : MOM_grid_init, MOM_grid_end, ocean_grid_type
    use MOM_verticalGrid,          only : verticalGrid_type, verticalGridInit, verticalGridEnd    
    
    implicit none

    type(ocean_grid_type),            intent(inout) :: G
    type(verticalGrid_type), pointer, intent(inout) :: GV

    ! Destructors below don't work
    !call MOM_grid_end(G)
    !call verticalGridEnd(GV)

  end subroutine soca_geom_end

  ! ------------------------------------------------------------------------------
  
  subroutine soca_field_init(aogcm, G, GV, ice_column)
    
    use MOM_state_initialization, only : MOM_initialize_state
    use MOM_time_manager,         only : time_type, set_time, time_type_to_real, operator(+)
    use MOM_file_parser,           only : get_param, param_file_type
    use MOM_get_input,             only : Get_MOM_Input, directories
    use MOM, only : MOM_control_struct, initialize_MOM
    use fms_io_mod,                only : fms_io_init, fms_io_exit
    use MOM_file_parser,           only : open_param_file, close_param_file, get_param
    use fms_mod,                   only : read_data
    
    implicit none
    
    type (Coupled),                     intent(out) :: aogcm
    type(ocean_grid_type), intent(inout)           :: G
    type(verticalGrid_type), pointer, intent(inout):: GV
    type(soca_ice_column),              intent(in) :: ice_column
    type(param_file_type) :: param_file  !< structure indicating paramater file to parse
    type(directories)     :: dirs        !< structure with directory paths

    integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nzo, nzi, nzs
    integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
    integer :: ncat, km

    ! Allocate arrays on data domain
    ! Note: Compute domain excludes halo (is, ie, js, je)
    !       Data domain includes halo (isd, ied, jsd, jed)
    is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nzo = G%ke
    Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
    isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
    IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

    ! TODO: Get from config 
    ncat = ice_column%ncat
    nzi = ice_column%nzi
    nzs = ice_column%nzs

    ! Allocate ocean state
    allocate(aogcm%Ocn%T(isd:ied,jsd:jed,nzo))   ; aogcm%Ocn%T(:,:,:) = 0.0_kind_real
    allocate(aogcm%Ocn%S(isd:ied,jsd:jed,nzo))   ; aogcm%Ocn%S(:,:,:) = 0.0_kind_real
    allocate(aogcm%Ocn%ssh(isd:ied,jsd:jed))   ; aogcm%Ocn%ssh(:,:) = 0.0_kind_real
    allocate(aogcm%Ocn%H(isd:ied,jsd:jed,nzo))   ; aogcm%Ocn%H(:,:,:) = 0.0_kind_real

    ! Allocate sea-ice state
    km = ncat + 1
    allocate(aogcm%Ice%part_size(isd:ied, jsd:jed, km)) ; aogcm%Ice%part_size(:,:,:) = 0.0_kind_real
    allocate(aogcm%Ice%h_ice(isd:ied, jsd:jed, ncat)) ; aogcm%Ice%h_ice(:,:,:) = 0.0_kind_real
    allocate(aogcm%Ice%h_snow(isd:ied, jsd:jed, ncat)) ; aogcm%Ice%h_snow(:,:,:) = 0.0_kind_real
    allocate(aogcm%Ice%enth_ice(isd:ied, jsd:jed, ncat, nzi)) ; aogcm%Ice%enth_ice(:,:,:,:) = 0.0_kind_real
    allocate(aogcm%Ice%enth_snow(isd:ied, jsd:jed, ncat, nzs)) ; aogcm%Ice%enth_snow(:,:,:,:) = 0.0_kind_real
    allocate(aogcm%Ice%sal_ice(isd:ied, jsd:jed, ncat, nzi)) ; aogcm%Ice%sal_ice(:,:,:,:) = 0.0_kind_real
    allocate(aogcm%Ice%T_skin(isd:ied, jsd:jed, ncat)) ; aogcm%Ice%T_skin(:,:,:) = 0.0_kind_real

  end subroutine soca_field_init

  ! ------------------------------------------------------------------------------
  
  subroutine soca_field_end(aogcm)
    
    implicit none
    
    type (Coupled),                     intent(inout) :: aogcm

    ! Deallocate ocean state
    deallocate(aogcm%Ocn%T)
    deallocate(aogcm%Ocn%S)
    deallocate(aogcm%Ocn%ssh)
    deallocate(aogcm%Ocn%H)    


    deallocate(aogcm%Ice%part_size)
    deallocate(aogcm%Ice%h_ice)
    deallocate(aogcm%Ice%h_snow)
    deallocate(aogcm%Ice%enth_ice)
    deallocate(aogcm%Ice%enth_snow)
    deallocate(aogcm%Ice%sal_ice)
    deallocate(aogcm%Ice%T_skin)

  end subroutine soca_field_end
  
end module soca_mom6
