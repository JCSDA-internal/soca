!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
!> Structure holding configuration variables for the  model

module soca_model_mod

  use kinds
  use iso_c_binding
  use soca_geom_mod
  
  implicit none
  
  private
  public :: soca_model
  public :: soca_model_registry
  public :: soca_create
  
  !> Fortran derived type to hold configuration data for the  model
  type :: soca_model
     integer :: nx                !< Zonal grid dimension
     integer :: ny                !< Meridional grid dimension
     real(kind=kind_real) :: dt0  !< dimensional time (seconds)
  end type soca_model

#define LISTED_TYPE soca_model

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_model_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

  subroutine soca_create(self, geom, c_conf)
    
  use MOM_cpu_clock,       only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
  use MOM_cpu_clock,       only : CLOCK_COMPONENT
  use MOM_diag_mediator,   only : enable_averaging, disable_averaging, diag_mediator_end
  use MOM_diag_mediator,   only : diag_ctrl, diag_mediator_close_registration
  use MOM,                 only : initialize_MOM, step_MOM, MOM_control_struct, MOM_end
  use MOM,                 only : extract_surface_state, finish_MOM_initialization
  use MOM,                 only : get_MOM_state_elements, MOM_state_is_synchronized
  use MOM,                 only : step_offline
  use MOM_domains,         only : MOM_infra_init, MOM_infra_end
  use MOM_error_handler,   only : MOM_error, MOM_mesg, WARNING, FATAL, is_root_pe
  use MOM_error_handler,   only : callTree_enter, callTree_leave, callTree_waypoint
  use MOM_file_parser,     only : read_param, get_param, log_param, log_version, param_file_type
  use MOM_file_parser,     only : close_param_file
  use MOM_forcing_type,    only : forcing, mech_forcing, forcing_diagnostics
  use MOM_forcing_type,    only : mech_forcing_diags, MOM_forcing_chksum, MOM_mech_forcing_chksum
  use MOM_get_input,       only : directories
  use MOM_grid,            only : ocean_grid_type
  use MOM_io,              only : file_exists, open_file, close_file
  use MOM_io,              only : check_nml_error, io_infra_init, io_infra_end
  use MOM_io,              only : APPEND_FILE, ASCII_FILE, READONLY_FILE, SINGLE_FILE
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
  use MOM_variables,       only : surface
  use MOM_verticalGrid,    only : verticalGrid_type
  use MOM_write_cputime,   only : write_cputime, MOM_write_cputime_init
  use MOM_write_cputime,   only : write_cputime_start_clock, write_cputime_CS

  use ensemble_manager_mod, only : ensemble_manager_init, get_ensemble_size
  use ensemble_manager_mod, only : ensemble_pelist_setup
  use mpp_mod, only : set_current_pelist => mpp_set_current_pelist
  use time_interp_external_mod, only : time_interp_external_init

  use MOM_ice_shelf, only : initialize_ice_shelf, ice_shelf_end, ice_shelf_CS
  use MOM_ice_shelf, only : shelf_calc_flux, add_shelf_forces, ice_shelf_save_restart
! , add_shelf_flux_forcing, add_shelf_flux_IOB

  use MOM_wave_interface, only: wave_parameters_CS, MOM_wave_interface_init
  use MOM_wave_interface, only: MOM_wave_interface_init_lite, Update_Surface_Waves

  implicit none

  ! A structure with the driving mechanical surface forces
  type(mech_forcing) :: forces
  ! A structure containing pointers to the thermodynamic forcing fields
  ! at the ocean surface.
  type(forcing) :: fluxes

  ! A structure containing pointers to the ocean surface state fields.
  type(surface) :: sfc_state

  ! A pointer to a structure containing metrics and related information.
  type(ocean_grid_type), pointer :: grid
  type(verticalGrid_type), pointer :: GV

  ! If .true., use the ice shelf model for part of the domain.
  logical :: use_ice_shelf

  ! If .true., use surface wave coupling
  logical :: use_waves = .false.

  ! This is .true. if incremental restart files may be saved.
  logical :: permit_incr_restart = .true.

  integer :: ns

  ! nmax is the number of iterations after which to stop so that the
  ! simulation does not exceed its CPU time limit.  nmax is determined by
  ! evaluating the CPU time used between successive calls to write_cputime.
  ! Initially it is set to be very large.
  integer :: nmax=2000000000

  ! A structure containing several relevant directory paths.
  type(directories) :: dirs

  ! A suite of time types for use by MOM
  type(time_type), target :: Time       ! A copy of the ocean model's time.
                                        ! Other modules can set pointers to this and
                                        ! change it to manage diagnostics.
  type(time_type) :: Master_Time        ! The ocean model's master clock. No other
                                        ! modules are ever given access to this.
  type(time_type) :: Time1              ! The value of the ocean model's time at the
                                        ! start of a call to step_MOM.
  type(time_type) :: Start_time         ! The start time of the simulation.
  type(time_type) :: segment_start_time ! The start time of this run segment.
  type(time_type) :: Time_end           ! End time for the segment or experiment.
  type(time_type) :: restart_time       ! The next time to write restart files.
  type(time_type) :: Time_step_ocean    ! A time_type version of dt_forcing.

  real    :: elapsed_time = 0.0   ! Elapsed time in this run in seconds.
  logical :: elapsed_time_master  ! If true, elapsed time is used to set the
                                  ! model's master clock (Time).  This is needed
                                  ! if Time_step_ocean is not an exact
                                  ! representation of dt_forcing.
  real :: dt_forcing              ! The coupling time step in seconds.
  real :: dt                      ! The baroclinic dynamics time step, in seconds.
  real :: dt_off                  ! Offline time step in seconds
  integer :: ntstep               ! The number of baroclinic dynamics time steps
                                  ! within dt_forcing.
  real :: dt_therm
  real :: dt_dyn, dtdia, t_elapsed_seg
  integer :: n, n_max, nts, n_last_thermo
  logical :: diabatic_first, single_step_call
  type(time_type) :: Time2, time_chg

  integer :: Restart_control    ! An integer that is bit-tested to determine whether
                                ! incremental restart files are saved and whether they
                                ! have a time stamped name.  +1 (bit 0) for generic
                                ! files and +2 (bit 1) for time-stamped files.  A
                                ! restart file is saved at the end of a run segment
                                ! unless Restart_control is negative.

  real            :: Time_unit       ! The time unit in seconds for the following input fields.
  type(time_type) :: restint         ! The time between saves of the restart file.
  type(time_type) :: daymax          ! The final day of the simulation.

  integer :: CPU_steps          ! The number of steps between writing CPU time.
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
  integer :: initClock, mainClock, termClock

  logical :: debug               ! If true, write verbose checksums for debugging purposes.
  logical :: offline_tracer_mode ! If false, use the model in prognostic mode where
                                 ! the barotropic and baroclinic dynamics, thermodynamics,
                                 ! etc. are stepped forward integrated in time.
                                 ! If true, then all of the above are bypassed with all
                                 ! fields necessary to integrate only the tracer advection
                                 ! and diffusion equation are read in from files stored from
                                 ! a previous integration of the prognostic model

  type(MOM_control_struct),  pointer :: MOM_CSp => NULL()
  !> A pointer to the tracer flow control structure.
  type(tracer_flow_control_CS), pointer :: &
    tracer_flow_CSp => NULL()  !< A pointer to the tracer flow control structure
  type(surface_forcing_CS),  pointer :: surface_forcing_CSp => NULL()
  type(write_cputime_CS),    pointer :: write_CPU_CSp => NULL()
  type(ice_shelf_CS),        pointer :: ice_shelf_CSp => NULL()
  type(wave_parameters_cs),  pointer :: waves_CSp => NULL()
  type(MOM_restart_CS),      pointer :: &
    restart_CSp => NULL()     !< A pointer to the restart control structure
                              !! that will be used for MOM restart files.
  type(diag_ctrl), pointer :: &
    diag => NULL()            !< A pointer to the diagnostic regulatory structure
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
  
  type(soca_model), intent(inout) :: self
  type(c_ptr),         intent(in) :: c_conf
  type(soca_geom),     intent(in) :: geom

  call write_cputime_start_clock(write_CPU_CSp)

  call MOM_infra_init() ; call io_infra_init()
  ! These clocks are on the global pelist.
  initClock = cpu_clock_id( 'Initialization' )
  mainClock = cpu_clock_id( 'Main loop' )
  termClock = cpu_clock_id( 'Termination' )
  call cpu_clock_begin(initClock)

  call MOM_mesg('======== Model being driven by MOM_driver ========', 2)
  call callTree_waypoint("Program MOM_main, MOM_driver.F90")

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

  Time = Start_time
  call initialize_MOM(Time, Start_time, param_file, dirs, MOM_CSp, restart_CSp, &
                      offline_tracer_mode=offline_tracer_mode, diag_ptr=diag, &
                      tracer_flow_CSp=tracer_flow_CSp)
  
  call get_MOM_state_elements(MOM_CSp, G=grid, GV=GV, C_p=fluxes%C_p)
  Master_Time = Time

  call callTree_waypoint("done initialize_MOM")

  call extract_surface_state(MOM_CSp, sfc_state)

  call surface_forcing_init(Time, grid, param_file, diag, &
                            surface_forcing_CSp, tracer_flow_CSp)
  call callTree_waypoint("done surface_forcing_init")

  segment_start_time = Time
  elapsed_time = 0.0

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod_name, version, "")
  call get_param(param_file, mod_name, "DT", param_int, fail_if_missing=.true.)
  dt = real(param_int)
  dt_forcing = dt
  Time_step_ocean = real_to_time(real(dt_forcing, kind=8))

  ! Close the param_file.  No further parsing of input is possible after this.
  call close_param_file(param_file)
  
  ! Set the forcing for the next steps.
  call set_forcing(sfc_state, forces, fluxes, Time, Time_step_ocean, grid, &
                     surface_forcing_CSp)
  
  ! This call steps the model over a time dt_forcing.
  Time1 = Master_Time ; Time = Master_Time
  call step_MOM(forces, fluxes, sfc_state, Time1, real(dt_forcing, kind=8), MOM_CSp, Waves=Waves_CSP)

end subroutine soca_create

end module soca_model_mod
