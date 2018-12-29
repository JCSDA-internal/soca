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

  use config_mod
  use fckit_mpi_module, only: fckit_mpi_comm
  
  use fms_io_mod, only : fms_io_init, fms_io_exit
  use fms_mod,    only : read_data, write_data, fms_init, fms_end

  use iso_c_binding  
  use kinds

  use MOM_grid,                  only : ocean_grid_type
  use MOM_verticalGrid,          only : verticalGrid_type
  use MOM, only : MOM_control_struct
  use MOM_forcing_type,    only : forcing, mech_forcing
  use MOM_variables,       only : surface
  use MOM_time_manager,    only : time_type
  use MOM,                 only : MOM_control_struct
  use MOM,                 only : initialize_MOM, step_MOM, MOM_control_struct, MOM_end
  use MOM,                 only : extract_surface_state, finish_MOM_initialization
  use MOM,                 only : get_MOM_state_elements, MOM_state_is_synchronized
  use MOM,                 only : step_offline
  use MOM_domains,         only : MOM_infra_init, MOM_infra_end
  use MOM_domains,         only : MOM_domains_init, clone_MOM_domain
  use MOM_dyn_horgrid,     only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid  
  use MOM_error_handler,   only : MOM_error, MOM_mesg, WARNING, FATAL, is_root_pe
  use MOM_error_handler,   only : callTree_enter, callTree_leave, callTree_waypoint
  use MOM_file_parser,     only : read_param, get_param, log_param, log_version, param_file_type
  use MOM_file_parser,     only : close_param_file
  use MOM_fixed_initialization,  only : MOM_initialize_fixed  
  use MOM_forcing_type,    only : forcing, mech_forcing, forcing_diagnostics
  use MOM_forcing_type,    only : mech_forcing_diags, MOM_forcing_chksum, MOM_mech_forcing_chksum
  use MOM_get_input,       only : directories, Get_MOM_Input, directories  
  use MOM_grid,            only : ocean_grid_type, MOM_grid_init
  use MOM_hor_index,       only : hor_index_type, hor_index_init  
  use MOM_io,              only : open_file, close_file
  use MOM_io,              only : check_nml_error, io_infra_init, io_infra_end
  use MOM_io,              only : APPEND_FILE, ASCII_FILE, READONLY_FILE, SINGLE_FILE
  use MOM_open_boundary,   only : ocean_OBC_type
  use MOM_restart,         only : MOM_restart_CS, save_restart
  use MOM_string_functions,only : uppercase
  use MOM_surface_forcing, only : set_forcing, forcing_save_restart
  use MOM_surface_forcing, only : surface_forcing_init, surface_forcing_CS
  use MOM_time_manager,    only : time_type, set_date, get_date
  use MOM_time_manager,    only : real_to_time, time_type_to_real
  use MOM_time_manager,    only : operator(+), operator(-), operator(*), operator(/)
  use MOM_time_manager,    only : operator(>), operator(<), operator(>=)
  use MOM_time_manager,    only : increment_date, set_calendar_type, month_name
  use MOM_time_manager,    only : JULIAN, GREGORIAN, NOLEAP, THIRTY_DAY_MONTHS
  use MOM_time_manager,    only : NO_CALENDAR
  use MOM_tracer_flow_control, only : tracer_flow_control_CS
  use MOM_transcribe_grid,     only : copy_dyngrid_to_MOM_grid  
  use MOM_variables,       only : surface
  use MOM_verticalGrid,    only : verticalGrid_type
  use MOM_verticalGrid,    only : verticalGrid_type, verticalGridInit, verticalGridEnd  

  use mpp_mod,    only : mpp_init
  
  use time_interp_external_mod, only : time_interp_external_init
  use time_manager_mod,         only: time_type
  
  implicit none
  private

  public :: soca_geom_init, soca_geom_end
  public :: soca_ice_column
  public :: soca_mom6_init, soca_mom6_config, soca_mom6_end

  type soca_ice_column
     integer                        :: ncat ! Number of ice categories
     integer                        :: nzi  ! Number of ice levels
     integer                        :: nzs  ! Number of snow levels
  end type soca_ice_column

  type soca_mom6_config
    type(mech_forcing) :: forces     ! Driving mechanical surface forces
    type(forcing)      :: fluxes     ! Pointers to the thermodynamic forcing fields
                                     ! at the ocean surface.
    type(surface)      :: sfc_state  ! Pointers to the ocean surface state fields.
    type(time_type)    :: Time       ! Model's time before call to step_MOM.
    real               :: dt_forcing ! Coupling time step in seconds.
    type(MOM_control_struct), pointer :: MOM_CSp ! Tracer flow control structure.
 end type soca_mom6_config
  
contains

  ! ------------------------------------------------------------------------------

  subroutine soca_geom_init(G, GV, ice_column, c_conf)

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
    type(fckit_mpi_comm) :: f_comm


    f_comm = fckit_mpi_comm()
    call mpp_init(localcomm=f_comm%communicator())

    ! Initialize fms
    call fms_init()

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
    nullify(GV)

  end subroutine soca_geom_end

  ! ------------------------------------------------------------------------------
  
  subroutine soca_mom6_init(mom6_config)

    use MOM_diag_mediator,   only : diag_ctrl, diag_mediator_close_registration

    implicit none

    type(ocean_grid_type), pointer :: grid
    type(verticalGrid_type), pointer :: GV
    type(directories) :: dirs
    type(time_type), target :: Time       ! A copy of the ocean model's time.
    type(time_type) :: Master_Time        ! The ocean model's master clock. No other
    type(time_type) :: Start_time         ! The start time of the simulation.
    type(time_type) :: Time_step_ocean    ! A time_type version of dt_forcing.
    real :: dt                      ! The baroclinic dynamics time step, in seconds.
    integer :: Restart_control    ! An integer that is bit-tested to determine whether
    ! incremental restart files are saved and whether they
    ! have a time stamped name.  +1 (bit 0) for generic
    ! files and +2 (bit 1) for time-stamped files.  A
    ! restart file is saved at the end of a run segment
    ! unless Restart_control is negative.

    real            :: Time_unit       ! The time unit in seconds for the following input fields.
    type(time_type) :: restint         ! The time between saves of the restart file.

    integer :: date_init(6)=0                ! The start date of the whole simulation.
    integer :: date(6)=-1                    ! Possibly the start date of this run segment.
    integer :: years=0, months=0, days=0     ! These may determine the segment run
    integer :: hours=0, minutes=0, seconds=0 ! length, if read from a namelist.
    integer :: yr, mon, day, hr, mins, sec   ! Temp variables for writing the date.
    type(param_file_type) :: param_file      ! The structure indicating the file(s)
    ! containing all run-time parameters.
    character(len=9)  :: month
    character(len=16) :: calendar = 'julian'
    integer :: calendar_type=-1

    integer :: unit, io_status, ierr
    integer :: ensemble_size, nPEs_per, ensemble_info(6)

    integer, dimension(0) :: atm_PElist, land_PElist, ice_PElist
    integer, dimension(:), allocatable :: ocean_PElist
    logical :: unit_in_use

    logical :: offline_tracer_mode = .false.

    type(tracer_flow_control_CS), pointer :: &
         tracer_flow_CSp => NULL()  !< A pointer to the tracer flow control structure
    type(surface_forcing_CS),  pointer :: surface_forcing_CSp => NULL()
    type(MOM_restart_CS),      pointer :: &
         restart_CSp => NULL()     !< A pointer to the restart control structure
    !! that will be used for MOM restart files.
    type(diag_ctrl), pointer :: diag => NULL() !< Diagnostic structure
    !-----------------------------------------------------------------------

    character(len=4), parameter :: vers_num = 'v2.0'
    ! This include declares and sets the variable "version".
#include "version_variable.h"
    character(len=40)  :: mod_name = "MOM_main (MOM_driver)" ! This module's name.

    integer :: ocean_nthreads = 1
    integer :: ncores_per_node = 36
    logical :: use_hyper_thread = .false.
    integer :: omp_get_num_threads,omp_get_thread_num,get_cpu_affinity,adder,base_cpu
    namelist /ocean_solo_nml/ date_init, calendar, months, days, hours, minutes, seconds,&
         ocean_nthreads, ncores_per_node, use_hyper_thread

    integer :: param_int

    type(soca_mom6_config), intent(out) :: mom6_config
    
    call MOM_infra_init() ; call io_infra_init()

    call MOM_mesg('======== Model being driven by soca ========', 2)
    call callTree_waypoint("soca_mom6::soca_mom6_init")

    ! Provide for namelist specification of the run length and calendar data.
    call open_file(unit, 'input.nml', form=ASCII_FILE, action=READONLY_FILE)
    read(unit, ocean_solo_nml, iostat=io_status)
    call close_file(unit)
    ierr = check_nml_error(io_status,'ocean_solo_nml')
    if (years+months+days+hours+minutes+seconds > 0) then
       if (is_root_pe()) write(*,ocean_solo_nml)
    endif
    calendar = uppercase(calendar)
    if (calendar(1:6) == 'JULIAN') then ;         calendar_type = JULIAN
    elseif (calendar(1:9) == 'GREGORIAN') then ; calendar_type = GREGORIAN
    elseif (calendar(1:6) == 'NOLEAP') then ;    calendar_type = NOLEAP
    elseif (calendar(1:10)=='THIRTY_DAY') then ; calendar_type = THIRTY_DAY_MONTHS
    elseif (calendar(1:11)=='NO_CALENDAR') then; calendar_type = NO_CALENDAR
    elseif (calendar(1:1) /= ' ') then
       call MOM_error(FATAL,'MOM_driver: Invalid namelist value '//trim(calendar)//' for calendar')
    else
       call MOM_error(FATAL,'MOM_driver: No namelist value for calendar')
    endif
    call set_calendar_type(calendar_type)

    Start_time = set_date(date_init(1),date_init(2), date_init(3), &
         date_init(4),date_init(5),date_init(6))

    call time_interp_external_init

    mom6_config%MOM_CSp => NULL()
    Time = Start_time
    call initialize_MOM(Time, Start_time, param_file, dirs, mom6_config%MOM_CSp, restart_CSp, &
         offline_tracer_mode=offline_tracer_mode, diag_ptr=diag, &
         tracer_flow_CSp=tracer_flow_CSp)

    call get_MOM_state_elements(mom6_config%MOM_CSp, G=grid, GV=GV, C_p=mom6_config%fluxes%C_p)
    Master_Time = Time

    call callTree_waypoint("done initialize_MOM")

    call extract_surface_state(mom6_config%MOM_CSp, mom6_config%sfc_state)

    call surface_forcing_init(Time, grid, param_file, diag, &
         surface_forcing_CSp, tracer_flow_CSp)
    call callTree_waypoint("done surface_forcing_init")

    ! Read all relevant parameters and write them to the model log.
    call log_version(param_file, mod_name, version, "")
    call get_param(param_file, mod_name, "DT", param_int, fail_if_missing=.true.)
    dt = real(param_int)
    mom6_config%dt_forcing = dt
    Time_step_ocean = real_to_time(real(mom6_config%dt_forcing, kind=8))

    ! Close the param_file.  No further parsing of input is possible after this.
    call close_param_file(param_file)

    ! Set the forcing for the next steps.
    call set_forcing(mom6_config%sfc_state, mom6_config%forces, mom6_config%fluxes,&
         & Time, Time_step_ocean, grid, &
         surface_forcing_CSp)

    ! Set Time
    mom6_config%Time = Master_Time
    
  end subroutine soca_mom6_init

  ! ------------------------------------------------------------------------------
  
  subroutine soca_mom6_end(mom6_config)
    
    implicit none
    
    type(soca_mom6_config), intent(out) :: mom6_config
    
    call MOM_end(mom6_config%MOM_CSp)

  end subroutine soca_mom6_end
  
end module soca_mom6
