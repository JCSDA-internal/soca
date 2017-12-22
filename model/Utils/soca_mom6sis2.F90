module soca_mom6sis2
  use ice_model_mod,           only: ice_data_type
  use ocean_model_mod,         only: ocean_state_type, ocean_public_type
  use time_manager_mod,        only: time_type
  use kinds
  use MOM_grid,                  only : ocean_grid_type
  use MOM_verticalGrid,          only : verticalGrid_type
  
  implicit none
  private

  public :: soca_geom_init, soca_geom_end, Coupled!, soca_model_geom

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to hold geometry data for MOM5 & CICE5 model
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

!!$  type :: soca_model_geom
!!$     type(ocean_grid_type) :: G
!!$     type(VerticalGrid_type), pointer :: GV          
!!$     integer :: nx
!!$     integer :: ny
!!$     integer :: nz
!!$     integer :: ncat
!!$     real(kind=kind_real), pointer :: lon(:,:)       !< The horizontal grid type     !< 2D array of longitude 
!!$     real(kind=kind_real), pointer :: lat(:,:)       !< 2D array of latitude
!!$     real(kind=kind_real), pointer :: z(:)           !<      
!!$     real(kind=kind_real), pointer :: mask2d(:,:)    !< 0 = land 1 = ocean surface mask only
!!$     real(kind=kind_real), pointer :: cell_area(:,:) !<
!!$  end type soca_model_geom

  
  type Coupled
     type   (ice_data_type)             :: Ice
     type (ocean_public_type)           :: Ocean
     type (ocean_state_type),  pointer  :: Ocean_state

     type (time_type)                   :: Time, Time_init, Time_end, Time_step_atmos, Time_step_cpld
     logical                            :: concurrent_ice
     logical                            :: initialized = .false.
     integer                            :: ensemble_id = 1
     integer, allocatable               :: ensemble_pelist(:,:)
  end type Coupled

contains

  subroutine soca_geom_init(G, GV)

    use MOM_dyn_horgrid,           only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
    use MOM_transcribe_grid,       only : copy_dyngrid_to_MOM_grid
    use MOM_verticalGrid,          only : verticalGrid_type, verticalGridInit, verticalGridEnd
    use MOM_grid,                  only : MOM_grid_init, ocean_grid_type
    use MOM_grid,                  only : ocean_grid_type
    use MOM_file_parser,           only : get_param, param_file_type
    use MOM_get_input,             only : Get_MOM_Input, directories
    use MOM_domains,               only : MOM_domains_init, clone_MOM_domain
    use MOM_hor_index,             only : hor_index_type, hor_index_init
    use fms_io_mod,                only : fms_io_init, fms_io_exit
    use MOM_domains,               only : MOM_infra_init, MOM_infra_end
    use MOM_io,                    only : check_nml_error, io_infra_init, io_infra_end, read_data
    use MOM_file_parser,           only : open_param_file, close_param_file
    use MOM_grid_initialize,       only : set_grid_metrics
    use MOM_open_boundary,         only : ocean_OBC_type
    use mpp_mod,                   only : mpp_pe, mpp_npes, mpp_root_pe, mpp_sync
    use MOM_fixed_initialization,  only : MOM_initialize_fixed

    implicit none

    type(ocean_grid_type),            intent(out) :: G          !< The horizontal grid type
    type(verticalGrid_type), pointer, intent(out) :: GV
    
    type(dyn_horgrid_type),       pointer  :: dG => NULL()

    type(hor_index_type)                   :: HI                        !< Horiz index array extents
    type(param_file_type)                  :: param_file                !< Structure to parse for run-time parameters
    type(directories)                      :: dirs                      !< Structure containing several relevant directory paths
    logical                                :: global_indexing = .false. !< If true use global horizontal index DOES NOT WORK
    logical                                :: write_geom_files = .false.
    type(ocean_OBC_type),          pointer :: OBC => NULL()

    ! Initialize fms/mpp io stuff
    call fms_io_init()

    ! Read grid inputs
    call Get_MOM_Input(param_file, dirs)

    ! Inintialize mpp domains
    call MOM_domains_init(G%domain, param_file)
    call hor_index_init(G%Domain, HI, param_file, local_indexing=.not.global_indexing)
    call create_dyn_horgrid(dG, HI)!, bathymetry_at_vel=bathy_at_vel)
    call clone_MOM_domain(G%Domain, dG%Domain)

    ! Allocate grid arrays
    GV => NULL()
    call verticalGridInit( param_file, GV )
    print *,'---------------------------- z=',GV%sLayer        
    call MOM_grid_init(G, param_file, HI, global_indexing=global_indexing)

    ! Read/Generate grid
    call MOM_initialize_fixed(dG, OBC, param_file, write_geom_files, dirs%output_directory)
    call copy_dyngrid_to_MOM_grid(dG, G)
    G%ke = GV%ke ; G%g_Earth = GV%g_Earth
    
    ! Deallocate dynamic grid
    call destroy_dyn_horgrid(dG)

    ! Deassociate pointers to dyn grid and vert grid
    dG => NULL()
    !GV => NULL()    

    call mpp_sync()
    print *,'+++++++++++++++++ out init geom ++++++++++++++++++++++++++'

    call close_param_file(param_file)
    call fms_io_exit()

  end subroutine soca_geom_init

  subroutine soca_geom_end(G, GV)
    use MOM_grid,                  only : MOM_grid_init, MOM_grid_end, ocean_grid_type
    use MOM_verticalGrid,          only : verticalGrid_type, verticalGridInit, verticalGridEnd    
    implicit none
    
    type(ocean_grid_type),            intent(inout) :: G
    type(verticalGrid_type), pointer, intent(inout) :: GV
    
    call MOM_grid_end(G)
    call verticalGridEnd(GV)
    
  end subroutine soca_geom_end
  
end module soca_mom6sis2
