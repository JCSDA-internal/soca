module soca_mom6sis2

  !use constants_mod,           only: constants_init

  use time_manager_mod,        only: time_type, set_calendar_type, set_time
  use time_manager_mod,        only: set_date, get_date, days_in_month, month_name
  use time_manager_mod,        only: operator(+), operator(-), operator (<)
  use time_manager_mod,        only: operator (>), operator ( /= ), operator ( / )
  use time_manager_mod,        only: operator (*), THIRTY_DAY_MONTHS, JULIAN
  use time_manager_mod,        only: NOLEAP, NO_CALENDAR, INVALID_CALENDAR
  use time_manager_mod,        only: date_to_string, increment_date
  use time_manager_mod,        only: operator(>=), operator(<=), operator(==)

  use fms_mod,                 only: open_namelist_file, field_exist, file_exist, check_nml_error
  use fms_mod,                 only: uppercase, error_mesg, write_version_number
  use fms_mod,                 only: fms_init, fms_end, stdout
  use fms_mod,                 only: read_data, write_data

  use fms_io_mod,              only: fms_io_init, fms_io_exit
  use fms_io_mod,              only: restart_file_type, register_restart_field, save_restart

  use diag_manager_mod,        only: diag_manager_init, diag_manager_end, diag_grid_end
  use diag_manager_mod,        only: DIAG_OCEAN, DIAG_OTHER, DIAG_ALL, get_base_date
  use diag_manager_mod,        only: diag_manager_set_time_end

  use field_manager_mod,       only: MODEL_ATMOS, MODEL_LAND, MODEL_ICE

  use tracer_manager_mod,      only: tracer_manager_init, get_tracer_index
  use tracer_manager_mod,      only: get_number_tracers, get_tracer_names, NO_TRACER

  use coupler_types_mod,       only: coupler_types_init

  use data_override_mod,       only: data_override_init

  !
  ! model interfaces used to couple the component models:
  !               atmosphere, land, ice, and ocean
  !

  use atmos_model_mod,         only: atmos_model_init, atmos_model_end
  use atmos_model_mod,         only: update_atmos_model_dynamics
  use atmos_model_mod,         only: update_atmos_model_down
  use atmos_model_mod,         only: update_atmos_model_up
  use atmos_model_mod,         only: atmos_data_type
  use atmos_model_mod,         only: land_ice_atmos_boundary_type
  use atmos_model_mod,         only: atmos_data_type_chksum
  use atmos_model_mod,         only: lnd_ice_atm_bnd_type_chksum
  use atmos_model_mod,         only: lnd_atm_bnd_type_chksum
  use atmos_model_mod,         only: ice_atm_bnd_type_chksum
  use atmos_model_mod,         only: atmos_model_restart
  use atmos_model_mod,         only: update_atmos_model_radiation
  use atmos_model_mod,         only: update_atmos_model_state

  use land_model_mod,          only: land_model_init, land_model_end
  use land_model_mod,          only: land_data_type, atmos_land_boundary_type
  use land_model_mod,          only: update_land_model_fast, update_land_model_slow
  use land_model_mod,          only: atm_lnd_bnd_type_chksum
  use land_model_mod,          only: land_data_type_chksum
  use land_model_mod,          only: land_model_restart

  use ice_model_mod,           only: ice_model_init, share_ice_domains, ice_model_end, ice_model_restart
  use ice_model_mod,           only: update_ice_model_fast, set_ice_surface_fields
  use ice_model_mod,           only: ice_data_type, land_ice_boundary_type
  use ice_model_mod,           only: ocean_ice_boundary_type, atmos_ice_boundary_type
  use ice_model_mod,           only: ice_data_type_chksum, ocn_ice_bnd_type_chksum
  use ice_model_mod,           only: atm_ice_bnd_type_chksum, lnd_ice_bnd_type_chksum
  use ice_model_mod,           only: unpack_ocean_ice_boundary, exchange_slow_to_fast_ice
  use ice_model_mod,           only: ice_model_fast_cleanup, unpack_land_ice_boundary
  use ice_model_mod,           only: exchange_fast_to_slow_ice, update_ice_model_slow

  use ocean_model_mod,         only: update_ocean_model, ocean_model_init,  ocean_model_end
  use ocean_model_mod,         only: ocean_public_type, ocean_state_type, ice_ocean_boundary_type
  use ocean_model_mod,         only: ocean_model_restart
  use ocean_model_mod,         only: ocean_public_type_chksum, ice_ocn_bnd_type_chksum
  !
  ! flux_ calls translate information between model grids - see flux_exchange.f90
  !

  use flux_exchange_mod,       only: flux_exchange_init, sfc_boundary_layer
  use flux_exchange_mod,       only: generate_sfc_xgrid, send_ice_mask_sic
  use flux_exchange_mod,       only: flux_down_from_atmos, flux_up_to_atmos
  use flux_exchange_mod,       only: flux_land_to_ice, flux_ice_to_ocean, flux_ocean_to_ice
  use flux_exchange_mod,       only: flux_check_stocks, flux_init_stocks
  use flux_exchange_mod,       only: flux_ocean_from_ice_stocks, flux_ice_to_ocean_stocks
  use flux_exchange_mod,       only: flux_atmos_to_ocean, flux_ex_arrays_dealloc

  use atmos_tracer_driver_mod, only: atmos_tracer_driver_gather_data

  use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_chksum
  use mpp_mod,                 only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_sync
  use mpp_mod,                 only: stderr, stdlog, mpp_error, NOTE, FATAL, WARNING
  use mpp_mod,                 only: mpp_set_current_pelist, mpp_declare_pelist
  use mpp_mod,                 only: input_nml_file

  use mpp_io_mod,              only: mpp_open, mpp_close, mpp_io_clock_on
  use mpp_io_mod,              only: MPP_NATIVE, MPP_RDONLY, MPP_DELETE

  use mpp_domains_mod,         only: mpp_broadcast_domain

  use memutils_mod,            only: print_memuse_stats

  implicit none
  private

  public :: soca_models_init, soca_models_end, Coupled, soca_simple_init

  type Coupled
     type (atmos_data_type)             :: Atm
     type  (land_data_type)             :: Land
     type   (ice_data_type)             :: Ice
     type (ocean_public_type)           :: Ocean
     type (ocean_state_type),  pointer  :: Ocean_state

     type(atmos_land_boundary_type)     :: Atmos_land_boundary
     type(atmos_ice_boundary_type)      :: Atmos_ice_boundary
     type(land_ice_atmos_boundary_type) :: Land_ice_atmos_boundary
     type(land_ice_boundary_type)       :: Land_ice_boundary
     type(ice_ocean_boundary_type)      :: Ice_ocean_boundary
     type(ocean_ice_boundary_type)      :: Ocean_ice_boundary

     type (time_type)                   :: Time, Time_init, Time_end, Time_step_atmos, Time_step_cpld
     logical                            :: concurrent_ice
     logical                            :: initialized = .false.
     integer                            :: ensemble_id = 1
     integer, allocatable               :: ensemble_pelist(:,:)
  end type Coupled

  !-----------------------------------------------------------------------

  character(len=128) :: version = '$Id$'
  character(len=128) :: tag = '$Name$'

  !-----------------------------------------------------------------------
  !---- model defined-types ----

  !  type (atmos_data_type) :: Atm
  !  type  (land_data_type) :: Land
  !  type   (ice_data_type) :: Ice
  !  ! allow members of ocean type to be aliased (ap)
  !  type (ocean_public_type), target :: Ocean
  !  type (ocean_state_type),  pointer :: Ocean_state => NULL()

  !type(atmos_land_boundary_type)     :: Atmos_land_boundary
  !type(atmos_ice_boundary_type)      :: Atmos_ice_boundary
  !type(land_ice_atmos_boundary_type) :: Land_ice_atmos_boundary
  !type(land_ice_boundary_type)       :: Land_ice_boundary
  !type(ice_ocean_boundary_type)      :: Ice_ocean_boundary
  !type(ocean_ice_boundary_type)      :: Ocean_ice_boundary

  !integer :: outunit
  !integer :: ensemble_id = 1
  !integer, allocatable :: ensemble_pelist(:, :)
  !integer :: conc_nthreads = 1
  !integer :: omp_get_thread_num, omp_get_num_threads
  !real :: omp_get_wtime
  !real :: dsec, omp_sec(2)=0.0, imb_sec(2)=0.0


contains

  subroutine soca_simple_init()
    use MOM_domains,        only: MOM_infra_init, num_pes, root_pe, pe_here
    use MOM_file_parser,    only: get_param, log_version, param_file_type
    use MOM_get_input,      only: Get_MOM_Input, directories

    !use ESMF,                only: ESMF_clock, ESMF_time, ESMF_timeInterval, ESMF_TimeInc
    !use ESMF,                only: ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet
    
    implicit none
    integer             :: mpicom_ocn        !< MPI ocn communicator
    integer             :: npes, pe0         !< # of processors and current processor                                  
    type (ocean_public_type), pointer           :: Ocean
    type (ocean_state_type),  pointer  :: Ocean_state => NULL()
    integer :: i

    type(time_type)     :: time_init         !< Start time of coupled model's calendar                                 
    type(time_type)     :: time_in           !< Time at the beginning of the first ocn coupling interval
    
    ! runi-time params                                                                                                 
    type(param_file_type) :: param_file           !< A structure to parse for run-time parameters                      
    type(directories)     :: dirs_tmp             !< A structure containing several relevant directory paths 

    integer :: date_init(6) = (/ 0, 0, 0, 0, 0, 0 /)
    integer :: date(6) = (/ 0, 0, 0, 0, 0, 0 /)    
    integer :: calendar_type = NOLEAP
    integer :: unit, ierr, io
    character(len=256) :: err_msg

    !-----------------------------------------------------------------------
    !------ namelist interface -------

    integer, dimension(6) :: restart_interval = (/ 0, 0, 0, 0, 0, 0/) !< The time interval that write out intermediate restart file.
    !! The format is (yr,mo,day,hr,min,sec).  When restart_interval
    !! is all zero, no intermediate restart file will be written out
    integer, dimension(6) :: current_date     = (/ 0, 0, 0, 0, 0, 0 /) !< The date that the current integration starts with.  (See
    !! force_date_from_namelist.)
    character(len=17) :: calendar = '                 ' !< The calendar type used by the current integration.  Valid values are
    !! consistent with the time_manager module: 'julian', 'noleap', or 'thirty_day'.
    !! The value 'no_calendar' cannot be used because the time_manager's date
    !! functions are used.  All values must be lower case.
    logical :: force_date_from_namelist = .false.  !< Flag that determines whether the namelist variable current_date should override
    !! the date in the restart file `INPUT/coupler.res`.  If the restart file does not
    !! exist then force_date_from_namelist has no effect, the value of current_date
    !! will be used.
    integer :: months=0  !< Number of months the current integration will be run
    integer :: days=0    !< Number of days the current integration will be run
    integer :: hours=0   !< Number of hours the current integration will be run
    integer :: minutes=0 !< Number of minutes the current integration will be run
    integer :: seconds=0 !< Number of seconds the current integration will be run
    integer :: dt_atmos = 0 !< Atmospheric model time step in seconds, including the fat coupling with land and sea ice
    integer :: dt_cpld  = 0 !< Time step in seconds for coupling between ocean and atmospheric models.  This must be an
    !! integral multiple of dt_atmos and dt_ocean.  This is the "slow" timestep.
    integer :: atmos_npes=0 !< The number of MPI tasks to use for the atmosphere
    integer :: ocean_npes=0 !< The number of MPI tasks to use for the ocean
    integer :: ice_npes=0   !< The number of MPI tasks to use for the ice
    integer :: land_npes=0  !< The number of MPI tasks to use for the land
    integer :: atmos_nthreads=1 !< Number of OpenMP threads to use in the atmosphere
    integer :: ocean_nthreads=1 !< Number of OpenMP threads to use in the ocean
    integer :: radiation_nthreads=1 !< Number of threads to use for the radiation.
    logical :: do_atmos =.true. !< Indicates if this component should be executed.  If .FALSE., then execution is skipped.
    !! This is used when ALL the output fields sent by this component to the coupler have been
    !! overridden  using the data_override feature.  This is for advanced users only.
    logical :: do_land =.true. !< See do_atmos
    logical :: do_ice =.true.  !< See do_atmos
    logical :: do_ocean=.true. !< See do_atmos
    logical :: do_flux =.true. !< See do_atmos
    logical :: concurrent=.FALSE. !< If .TRUE., the ocean executes concurrently with the atmosphere-land-ice on a separate
    !! set of PEs.  Concurrent should be .TRUE. if concurrent_ice is .TRUE.
    !! If .FALSE., the execution is serial: call atmos... followed by call ocean...
    logical :: do_concurrent_radiation=.FALSE. !< If .TRUE. then radiation is done concurrently
    logical :: use_lag_fluxes=.TRUE.  !< If .TRUE., the ocean is forced with SBCs from one coupling timestep ago.
    !! If .FALSE., the ocean is forced with most recent SBCs.  For an old leapfrog
    !! MOM4 coupling with dt_cpld=dt_ocean, lag fluxes can be shown to be stable
    !! and current fluxes to be unconditionally unstable.  For dt_cpld>dt_ocean there
    !! is probably sufficient damping for MOM4.  For more modern ocean models (such as
    !! MOM5, GOLD or MOM6) that do not use leapfrog timestepping, use_lag_fluxes=.False.
    !! should be much more stable.
    logical :: concurrent_ice=.FALSE. !< If .TRUE., the slow sea-ice is forced with the fluxes that were used for the
    !! fast ice processes one timestep before.  When used in conjuction with setting
    !! slow_ice_with_ocean=.TRUE., this approach allows the atmosphere and
    !! ocean to run concurrently even if use_lag_fluxes=.FALSE., and it can
    !! be shown to ameliorate or eliminate several ice-ocean coupled instabilities.
    logical :: slow_ice_with_ocean=.FALSE.  !< If true, the slow sea-ice is advanced on the ocean processors.  Otherwise
    !! the slow sea-ice processes are on the same PEs as the fast sea-ice.
    logical :: do_chksum=.FALSE.      !! If .TRUE., do multiple checksums throughout the execution of the model.
    logical :: do_endpoint_chksum=.FALSE.  !< If .TRUE., do checksums of the initial and final states.
    logical :: do_debug=.FALSE.       !< If .TRUE. print additional debugging messages.
    integer :: check_stocks = 0 ! -1: never 0: at end of run only n>0: every n coupled steps
    logical :: use_hyper_thread = .false.
    integer :: ncores_per_node = 0
    logical :: debug_affinity = .false.

    namelist /coupler_nml/ current_date, calendar, force_date_from_namelist,         &
         months, days, hours, minutes, seconds, dt_cpld, dt_atmos, &
         do_atmos, do_land, do_ice, do_ocean, do_flux,             &
         atmos_npes, ocean_npes, ice_npes, land_npes,              &
         atmos_nthreads, ocean_nthreads, radiation_nthreads,       &
         concurrent, do_concurrent_radiation, use_lag_fluxes,      &
         check_stocks, restart_interval, do_debug, do_chksum,      &
         use_hyper_thread, ncores_per_node, debug_affinity,        &
         concurrent_ice, slow_ice_with_ocean, do_endpoint_chksum
    
    call MOM_infra_init(mpicom_ocn)
    npes = num_pes()
    pe0 = root_pe()
    print *,npes, pe0

    allocate(Ocean)
    Ocean_state%is_ocean_PE = .true.
    !allocate(Ocean_state%pelist(npes))
    !Ocean_state%pelist(:) = (/(i,i=pe0,pe0+npes)/)

    unit = open_namelist_file()
    ierr=1;
    do while (ierr /= 0)
       read  (unit, nml=coupler_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'coupler_nml')
    enddo
10  call mpp_close(unit)
    
    call get_MOM_Input(param_file, dirs_tmp)!, check_params=.false.)

    !call ESMF_TimeGet(current_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=seconds, rc=rc)
    !time_init = set_date(year, month, day, hour, minute, seconds, err_msg=err_msg)

    !Balaji: currently written in binary, needs form=MPP_NATIVE
    !call mpp_open( unit, 'INPUT/coupler.res', action=MPP_RDONLY )
    !read( unit,* )calendar_type
    !read( unit,* )date_init
    !read( unit,* )date
    !call mpp_close(unit)

    call set_calendar_type (calendar_type, err_msg)

    date_init = current_date
    date = current_date
    
    time_init = set_date (date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6))
    time_in      = set_date (date(1), date(2), date(3),  &
         date(4), date(5), date(6))


    Ocean_state%Time = Time_in
    !call initialize_MOM(Ocean_state%Time, param_file, Ocean_state%dirs, Ocean_state%MOM_CSp, Time_in)!, &
       !offline_tracer_mode=offline_tracer_mode)

    !call ocean_model_init(Ocean, Ocean_state, time_init, time_in)


    
  end subroutine soca_simple_init
    
  
  !#######################################################################

  !> \brief Initialize all defined exchange grids and all boundary maps
  subroutine soca_models_init(AOGCM) !Atm, Land, Ice, Ocean, Ocean_state)

    use ensemble_manager_mod, only : ensemble_manager_init, get_ensemble_id,ensemble_pelist_setup
    use ensemble_manager_mod, only : get_ensemble_size, get_ensemble_pelist

    implicit none
    !-----------------------------------------------------------------------
    !---- model defined-types ----

    type (Coupled)                    :: AOGCM

    !integer :: ensemble_id = 1
    !integer, allocatable :: ensemble_pelist(:, :)
    integer :: conc_nthreads = 1
    character(len=80) :: text

    integer :: id_atmos_model_init, id_land_model_init, id_ice_model_init
    integer :: id_ocean_model_init, id_flux_exchange_init

    !-----------------------------------------------------------------------
    !------ namelist interface -------

    integer, dimension(6) :: restart_interval = (/ 0, 0, 0, 0, 0, 0/) !< The time interval that write out intermediate restart file.
    !! The format is (yr,mo,day,hr,min,sec).  When restart_interval
    !! is all zero, no intermediate restart file will be written out
    integer, dimension(6) :: current_date     = (/ 0, 0, 0, 0, 0, 0 /) !< The date that the current integration starts with.  (See
    !! force_date_from_namelist.)
    character(len=17) :: calendar = '                 ' !< The calendar type used by the current integration.  Valid values are
    !! consistent with the time_manager module: 'julian', 'noleap', or 'thirty_day'.
    !! The value 'no_calendar' cannot be used because the time_manager's date
    !! functions are used.  All values must be lower case.
    logical :: force_date_from_namelist = .false.  !< Flag that determines whether the namelist variable current_date should override
    !! the date in the restart file `INPUT/coupler.res`.  If the restart file does not
    !! exist then force_date_from_namelist has no effect, the value of current_date
    !! will be used.
    integer :: months=0  !< Number of months the current integration will be run
    integer :: days=0    !< Number of days the current integration will be run
    integer :: hours=0   !< Number of hours the current integration will be run
    integer :: minutes=0 !< Number of minutes the current integration will be run
    integer :: seconds=0 !< Number of seconds the current integration will be run
    integer :: dt_atmos = 0 !< Atmospheric model time step in seconds, including the fat coupling with land and sea ice
    integer :: dt_cpld  = 0 !< Time step in seconds for coupling between ocean and atmospheric models.  This must be an
    !! integral multiple of dt_atmos and dt_ocean.  This is the "slow" timestep.
    integer :: atmos_npes=0 !< The number of MPI tasks to use for the atmosphere
    integer :: ocean_npes=0 !< The number of MPI tasks to use for the ocean
    integer :: ice_npes=0   !< The number of MPI tasks to use for the ice
    integer :: land_npes=0  !< The number of MPI tasks to use for the land
    integer :: atmos_nthreads=1 !< Number of OpenMP threads to use in the atmosphere
    integer :: ocean_nthreads=1 !< Number of OpenMP threads to use in the ocean
    integer :: radiation_nthreads=1 !< Number of threads to use for the radiation.
    logical :: do_atmos =.true. !< Indicates if this component should be executed.  If .FALSE., then execution is skipped.
    !! This is used when ALL the output fields sent by this component to the coupler have been
    !! overridden  using the data_override feature.  This is for advanced users only.
    logical :: do_land =.true. !< See do_atmos
    logical :: do_ice =.true.  !< See do_atmos
    logical :: do_ocean=.true. !< See do_atmos
    logical :: do_flux =.true. !< See do_atmos
    logical :: concurrent=.FALSE. !< If .TRUE., the ocean executes concurrently with the atmosphere-land-ice on a separate
    !! set of PEs.  Concurrent should be .TRUE. if concurrent_ice is .TRUE.
    !! If .FALSE., the execution is serial: call atmos... followed by call ocean...
    logical :: do_concurrent_radiation=.FALSE. !< If .TRUE. then radiation is done concurrently
    logical :: use_lag_fluxes=.TRUE.  !< If .TRUE., the ocean is forced with SBCs from one coupling timestep ago.
    !! If .FALSE., the ocean is forced with most recent SBCs.  For an old leapfrog
    !! MOM4 coupling with dt_cpld=dt_ocean, lag fluxes can be shown to be stable
    !! and current fluxes to be unconditionally unstable.  For dt_cpld>dt_ocean there
    !! is probably sufficient damping for MOM4.  For more modern ocean models (such as
    !! MOM5, GOLD or MOM6) that do not use leapfrog timestepping, use_lag_fluxes=.False.
    !! should be much more stable.
    logical :: concurrent_ice=.FALSE. !< If .TRUE., the slow sea-ice is forced with the fluxes that were used for the
    !! fast ice processes one timestep before.  When used in conjuction with setting
    !! slow_ice_with_ocean=.TRUE., this approach allows the atmosphere and
    !! ocean to run concurrently even if use_lag_fluxes=.FALSE., and it can
    !! be shown to ameliorate or eliminate several ice-ocean coupled instabilities.
    logical :: slow_ice_with_ocean=.FALSE.  !< If true, the slow sea-ice is advanced on the ocean processors.  Otherwise
    !! the slow sea-ice processes are on the same PEs as the fast sea-ice.
    logical :: do_chksum=.FALSE.      !! If .TRUE., do multiple checksums throughout the execution of the model.
    logical :: do_endpoint_chksum=.FALSE.  !< If .TRUE., do checksums of the initial and final states.
    logical :: do_debug=.FALSE.       !< If .TRUE. print additional debugging messages.
    integer :: check_stocks = 0 ! -1: never 0: at end of run only n>0: every n coupled steps
    logical :: use_hyper_thread = .false.
    integer :: ncores_per_node = 0
    logical :: debug_affinity = .false.

    namelist /coupler_nml/ current_date, calendar, force_date_from_namelist,         &
         months, days, hours, minutes, seconds, dt_cpld, dt_atmos, &
         do_atmos, do_land, do_ice, do_ocean, do_flux,             &
         atmos_npes, ocean_npes, ice_npes, land_npes,              &
         atmos_nthreads, ocean_nthreads, radiation_nthreads,       &
         concurrent, do_concurrent_radiation, use_lag_fluxes,      &
         check_stocks, restart_interval, do_debug, do_chksum,      &
         use_hyper_thread, ncores_per_node, debug_affinity,        &
         concurrent_ice, slow_ice_with_ocean, do_endpoint_chksum

    integer :: initClock, mainClock, termClock

    !
    !-----------------------------------------------------------------------
    !     local parameters
    !-----------------------------------------------------------------------
    !
    character(len=48), parameter    :: mod_name = 'coupler_main_mod'
    character(len=64), parameter    :: sub_name = 'soca_models_init'
    character(len=256), parameter   :: error_header =                               &
         '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
    character(len=256), parameter   :: note_header =                                &
         '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    integer :: unit,  ierr, io,    m, i, outunit, logunit, errunit
    integer :: date(6)
    type (time_type) :: Run_length
    character(len=9) :: month
    integer :: pe, npes

    integer :: ens_siz(6), ensemble_size

    integer :: atmos_pe_start=0, atmos_pe_end=0, &
         ocean_pe_start=0, ocean_pe_end=0
    integer :: n
    integer :: diag_model_subset=DIAG_ALL
    logical :: other_fields_exist
    character(len=256) :: err_msg
    integer :: date_restart(6)
    character(len=64)  :: filename, fieldname
    integer :: id_restart, l
    integer :: omp_get_thread_num, omp_get_num_threads
    integer :: get_cpu_affinity, base_cpu, base_cpu_r, adder
    character(len=8)  :: walldate
    character(len=10) :: walltime
    character(len=5)  :: wallzone
    integer           :: wallvalues(8)

    !-----------------------------------------------------------------------
    ! ----- coupled model time -----

    type (time_type) :: Time, Time_init, Time_end, &
         Time_step_atmos, Time_step_cpld
    type(time_type) :: Time_atmos, Time_ocean
    integer :: num_atmos_calls, na
    integer :: num_cpld_calls, nc

    !------ for intermediate restart
    type(restart_file_type), allocatable :: Ice_bc_restart(:), Ocn_bc_restart(:)
    character(len=64),       allocatable :: ice_bc_restart_file(:), ocn_bc_restart_file(:)
    integer                              :: num_ice_bc_restart=0, num_ocn_bc_restart=0
    type(time_type)                      :: Time_restart, Time_restart_current, Time_start
    character(len=32)                    :: timestamp

    ! ----- coupled model initial date -----

    integer :: date_init(6) = (/ 0, 0, 0, 0, 0, 0 /)
    integer :: calendar_type = INVALID_CALENDAR


    

    !-----------------------------------------------------------------------

    print *,'[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[['
    print *,'[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[['
    print *,'[[[[[[[[[[[[[[[[[[[[[[[[[[    START      [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[['
    print *,'[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[['
    
    outunit = stdout()
    errunit = stderr()
    logunit = stdlog()

    if (mpp_pe().EQ.mpp_root_pe()) then
       call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
       write(errunit,*) 'Entering soca_models_init at '&
            //trim(walldate)//' '//trim(walltime)
    endif

    !----- write version to logfile -------
    call write_version_number(version, tag)

    !----- read namelist -------

    unit = open_namelist_file()
    ierr=1;
    do while (ierr /= 0)
       read  (unit, nml=coupler_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'coupler_nml')
    enddo
10  call mpp_close(unit)
    !#endif

    !---- when concurrent is set true and mpp_io_nml io_clock_on is set true, the model
    !---- will crash with error message "MPP_CLOCK_BEGIN: cannot change pelist context of a clock",
    !---- so need to make sure it will not happen
    if (concurrent) then
       if (mpp_io_clock_on()) then
          call error_mesg ('program coupler', 'when coupler_nml variable concurrent is set to true, '// &
               'mpp_io_nml variable io_clock_non can not be set to true.', FATAL )
       endif
    endif

    !---- ncores_per_node must be set when use_hyper_thread = .true.
    if (use_hyper_thread .and. ncores_per_node == 0) then
       call error_mesg ('program copuler', 'coupler_nml ncores_per_node must be set when use_hyper_thread=true', FATAL)
    endif

    !----- read date and calendar type from restart file -----

    force_date_from_namelist = .true.

    !----- use namelist value (either no restart or override flag on) ---

    if ( force_date_from_namelist ) then

       if ( sum(current_date) <= 0 ) then
          call error_mesg ('program coupler',  &
               'no namelist value for base_date or current_date', FATAL)
       else
          date      = current_date
       endif

       !----- override calendar type with namelist value -----

       select case( uppercase(trim(calendar)) )
       case( 'JULIAN' )
          calendar_type = JULIAN
       case( 'NOLEAP' )
          calendar_type = NOLEAP
       case( 'THIRTY_DAY' )
          calendar_type = THIRTY_DAY_MONTHS
       case( 'NO_CALENDAR' )
          calendar_type = NO_CALENDAR
       end select

    endif

    call set_calendar_type (calendar_type, err_msg)
    if (err_msg /= '') then
       call mpp_error(FATAL, 'ERROR in soca_models_init: '//trim(err_msg))
    endif

    if (concurrent .AND. .NOT.(use_lag_fluxes .OR. concurrent_ice) ) &
         call mpp_error( WARNING, 'soca_models_init: you have set concurrent=TRUE, &
         & use_lag_fluxes=FALSE, and concurrent_ice=FALSE &
         & in coupler_nml. When not using lag fluxes, components &
         & will synchronize at two points, and thus run serially.' )
    if (concurrent_ice .AND. .NOT.slow_ice_with_ocean ) call mpp_error(WARNING, &
         'soca_models_init: concurrent_ice is true, but slow ice_with_ocean is &
         & false in coupler_nml.  These two flags should both be true to avoid &
         & effectively serializing the run.' )
    if (use_lag_fluxes .AND. concurrent_ice ) call mpp_error(WARNING, &
         'soca_models_init: use_lag_fluxes and concurrent_ice are both true. &
         & These two coupling options are intended to be exclusive.' )

    !Check with the ensemble_manager module for the size of ensemble
    !and PE counts for each member of the ensemble.
    !
    !NOTE: ensemble_manager_init renames all the output files (restart and diagnostics)
    !      to show which ensemble member they are coming from.
    !      There also need to be restart files for each member of the ensemble in INPUT.
    !
    !NOTE: if the ensemble_size=1 the input/output files will not be renamed.
    !

    !if (mpp_pe().EQ.mpp_root_pe()) then
    !   call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
    !   write(errunit,*) 'Starting initializing ensemble_manager at '&
    !        //trim(walldate)//' '//trim(walltime)
    !endif

    call ensemble_manager_init() ! init pelists for ensembles       

    if (mpp_pe().EQ.mpp_root_pe()) then
       call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
       write(errunit,*) 'Finished initializing ensemble_manager at '&
            //trim(walldate)//' '//trim(walltime)
    endif
    ens_siz = get_ensemble_size()
    ensemble_size = ens_siz(1)
    npes = ens_siz(2)
    !Check for the consistency of PE counts
    if (concurrent) then
       !atmos_npes + ocean_npes must equal npes
       if (atmos_npes.EQ.0 ) atmos_npes = npes - ocean_npes
       if (ocean_npes.EQ.0 ) ocean_npes = npes - atmos_npes
       !both must now be non-zero
       if (atmos_npes.EQ.0 .OR. ocean_npes.EQ.0 ) &
            call mpp_error( FATAL, 'soca_models_init: atmos_npes or ocean_npes must be specified for concurrent coupling.' )
       if (atmos_npes+ocean_npes.NE.npes ) &
            call mpp_error( FATAL, 'soca_models_init: atmos_npes+ocean_npes must equal npes for concurrent coupling.' )
    else                        !serial timestepping
       if ((atmos_npes.EQ.0) .and. (do_atmos .or. do_land .or. do_ice)) atmos_npes = npes
       if ((ocean_npes.EQ.0) .and. (do_ocean)) ocean_npes = npes
       if (max(atmos_npes,ocean_npes).EQ.npes) then !overlapping pelists
          ! do nothing
       else                    !disjoint pelists
          if (atmos_npes+ocean_npes.NE.npes ) call mpp_error( FATAL,  &
               'soca_models_init: atmos_npes+ocean_npes must equal npes for serial coupling on disjoint pelists.' )
       endif
    endif
    if (land_npes == 0 ) land_npes = atmos_npes
    if (land_npes > atmos_npes) call mpp_error(FATAL, 'soca_models_init: land_npes > atmos_npes')

    if (ice_npes  == 0 ) ice_npes  = atmos_npes
    if (ice_npes  > atmos_npes) call mpp_error(FATAL, 'soca_models_init: ice_npes > atmos_npes')
    print *,atmos_npes
    print *,AOGCM%Atm%pelist
    allocate( AOGCM%Atm%pelist  (atmos_npes) )
    allocate( AOGCM%Ocean%pelist(ocean_npes) )
    allocate( AOGCM%Land%pelist (land_npes) )
    allocate( AOGCM%Ice%fast_pelist(ice_npes) )

    !Set up and declare all the needed pelists
    allocate(AOGCM%ensemble_pelist(1:ensemble_size,1:npes))

    call ensemble_pelist_setup(concurrent, atmos_npes, ocean_npes, land_npes, ice_npes, &
         AOGCM%Atm%pelist, AOGCM%Ocean%pelist, AOGCM%Land%pelist, AOGCM%Ice%fast_pelist)
    !set up affinities based on threads

    AOGCM%ensemble_id = get_ensemble_id()
    call get_ensemble_pelist(AOGCM%ensemble_pelist)
    AOGCM%Atm%pe            = ANY(AOGCM%Atm%pelist   .EQ. mpp_pe())
    AOGCM%Ocean%is_ocean_pe = ANY(AOGCM%Ocean%pelist .EQ. mpp_pe())
    AOGCM%Land%pe           = ANY(AOGCM%Land%pelist  .EQ. mpp_pe())

    AOGCM%Ice%shared_slow_fast_PEs = .not.slow_ice_with_ocean
    ! This is where different settings would be applied if the fast and slow
    ! ice occurred on different PEs.
    if (AOGCM%Ice%shared_slow_fast_PEs) then
       ! Fast and slow ice processes occur on the same PEs.
       allocate( AOGCM%Ice%pelist  (ice_npes) )
       AOGCM%Ice%pelist(:) = AOGCM%Ice%fast_pelist(:)
       allocate( AOGCM%Ice%slow_pelist(ice_npes) )
       AOGCM%Ice%slow_pelist(:) = AOGCM%Ice%fast_pelist(:)
    else
       ! Fast ice processes occur a subset of the atmospheric PEs, while
       ! slow ice processes occur on the ocean PEs.
       !### THIS OPTION IS WHOLLY UNTESTED AND I DO NOT KNOW IF IT WILL WORK! - RWH
       allocate( AOGCM%Ice%slow_pelist(ocean_npes) )
       AOGCM%Ice%slow_pelist(:) = AOGCM%Ocean%pelist(:)
       allocate( AOGCM%Ice%pelist  (ice_npes+ocean_npes) )
       ! Set AOGCM%Ice%pelist() to be the union of AOGCM%Ice%fast_pelist and AOGCM%Ice%slow_pelist.
       AOGCM%Ice%pelist(1:ice_npes) = AOGCM%Ice%fast_pelist(:)
       AOGCM%Ice%pelist(ice_npes+1:ice_npes+ocean_npes) = AOGCM%Ocean%pelist(:)
    endif

    AOGCM%Ice%fast_ice_pe = ANY(AOGCM%Ice%fast_pelist(:) .EQ. mpp_pe())
    AOGCM%Ice%slow_ice_pe = ANY(AOGCM%Ice%slow_pelist(:) .EQ. mpp_pe())
    AOGCM%Ice%pe = AOGCM%Ice%fast_ice_pe .OR. AOGCM%Ice%slow_ice_pe

    !Why is the following needed?
    !$  call omp_set_dynamic(.FALSE.)
    !$  call omp_set_nested(.TRUE.)
    if (AOGCM%Atm%pe) then
       call mpp_set_current_pelist( AOGCM%Atm%pelist )
    endif

    !--- initialization clock
    if (AOGCM%Atm%pe) then
       call mpp_set_current_pelist(AOGCM%Atm%pelist)
       id_atmos_model_init = mpp_clock_id( '  Init: atmos_model_init ' )
    endif
    if (AOGCM%Land%pe) then
       call mpp_set_current_pelist(AOGCM%Land%pelist)
       id_land_model_init  = mpp_clock_id( '  Init: land_model_init ' )
    endif
    if (AOGCM%Ice%pe) then
       if (AOGCM%Ice%shared_slow_fast_PEs) then
          call mpp_set_current_pelist(AOGCM%Ice%pelist)
       elseif (AOGCM%Ice%fast_ice_pe) then
          call mpp_set_current_pelist(AOGCM%Ice%fast_pelist)
       elseif (AOGCM%Ice%slow_ice_pe) then
          call mpp_set_current_pelist(AOGCM%Ice%slow_pelist)
       else
          call mpp_error(FATAL, "All Ice%pes must be a part of Ice%fast_ice_pe or Ice%slow_ice_pe")
       endif
       id_ice_model_init   = mpp_clock_id( '  Init: ice_model_init ' )
    endif
    if (AOGCM%Ocean%is_ocean_pe) then
       call mpp_set_current_pelist(AOGCM%Ocean%pelist)
       id_ocean_model_init = mpp_clock_id( '  Init: ocean_model_init ' )
    endif
    call mpp_set_current_pelist(AOGCM%ensemble_pelist(AOGCM%ensemble_id,:))
    id_flux_exchange_init = mpp_clock_id( '  Init: flux_exchange_init' )

    call mpp_set_current_pelist()
    mainClock = mpp_clock_id( 'Main loop' )
    termClock = mpp_clock_id( 'Termination' )

    !Write out messages on root PEs
    if (mpp_pe().EQ.mpp_root_pe()) then
       write( text,'(a,2i6,a,i2.2)' )'Atmos PE range: ', AOGCM%Atm%pelist(1)  , AOGCM%Atm%pelist(atmos_npes)  ,&
            ' ens_', AOGCM%ensemble_id
       call mpp_error( NOTE, 'soca_models_init: '//trim(text) )
       if (ocean_npes .gt. 0) then
          write( text,'(a,2i6,a,i2.2)' )'Ocean PE range: ', AOGCM%Ocean%pelist(1), AOGCM%Ocean%pelist(ocean_npes), &
               ' ens_', AOGCM%ensemble_id
          call mpp_error( NOTE, 'soca_models_init: '//trim(text) )
       else
          write( text,'(a,i2.2)' )'Ocean PE range is not set (do_ocean=.false. and concurrent=.false.) for ens_', &
               AOGCM%ensemble_id
          call mpp_error( NOTE, 'soca_models_init: '//trim(text) )
       endif
       write( text,'(a,2i6,a,i2.2)' )'Land PE range: ', AOGCM%Land%pelist(1)  , AOGCM%Land%pelist(land_npes)  ,&
            ' ens_', AOGCM%ensemble_id
       call mpp_error( NOTE, 'soca_models_init: '//trim(text) )
       write( text,'(a,2i6,a,i2.2)' )'Ice PE range: ', AOGCM%Ice%pelist(1), AOGCM%Ice%pelist(ice_npes), &
            ' ens_', AOGCM%ensemble_id
       call mpp_error( NOTE, 'soca_models_init: '//trim(text) )

       if (concurrent) then
          call mpp_error( NOTE, 'soca_models_init: Running with CONCURRENT coupling.' )

          write( logunit,'(a)' )'Using concurrent coupling...'
          write( logunit,'(a,4i6)' ) &
               'atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end=', &
               AOGCM%Atm%pelist(1)  , AOGCM%Atm%pelist(atmos_npes), AOGCM%Ocean%pelist(1), AOGCM%Ocean%pelist(ocean_npes)
       else
          call mpp_error( NOTE, 'soca_models_init: Running with SERIAL coupling.' )
       endif
       if (use_lag_fluxes) then
          call mpp_error( NOTE, 'soca_models_init: Sending LAG fluxes to ocean.' )
       else
          call mpp_error( NOTE, 'soca_models_init: Sending most recent fluxes to ocean.' )
       endif
       if (concurrent_ice ) call mpp_error( NOTE, &
            'soca_models_init: using lagged slow-ice coupling mode.')
    endif

    !----- write namelist to logfile -----
    if (mpp_pe() == mpp_root_pe() )write( logunit, nml=coupler_nml )

    !----- write current/initial date actually used to logfile file -----

    if ( mpp_pe().EQ.mpp_root_pe() ) &
         write( logunit, 16 )date(1),trim(month_name(date(2))),date(3:6)
16  format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt')

    !-----------------------------------------------------------------------
    !------ initialize diagnostics manager ------

    !jwd Fork here is somewhat dangerous. It relies on "no side effects" from
    !    diag_manager_init. diag_manager_init or this section should be
    !    re-architected to guarantee this or remove this assumption.
    !    For instance, what follows assumes that get_base_date has the same
    !    time for both Atm and Ocean pes. While this should be the case, the
    !    possible error condition needs to be checked

!!$    !diag_model_subset=DIAG_ALL
!!$    if (AOGCM%Atm%pe) then
!!$       call mpp_set_current_pelist(AOGCM%Atm%pelist)
!!$       !if (atmos_npes /= npes) diag_model_subset = DIAG_OTHER  ! change diag_model_subset from DIAG_ALL
!!$    elseif (AOGCM%Ocean%is_ocean_pe) then  ! Error check above for disjoint pelists should catch any problem
!!$       call mpp_set_current_pelist(AOGCM%Ocean%pelist)
!!$       ! The FMS diag manager has a convention that segregates files with "ocean"
!!$       ! in their names from the other files to handle long diag tables.  This
!!$       ! does not work if the ice is on the ocean PEs.
!!$       !if ((ocean_npes /= npes) .and. .not.slow_ice_with_ocean) &
!!$       !     diag_model_subset = DIAG_OCEAN  ! change diag_model_subset from DIAG_ALL
!!$    endif
    !if ( mpp_pe() == mpp_root_pe()) then
    !   call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
    !   write(errunit,*) 'Starting to initialize diag_manager at '&
    !        //trim(walldate)//' '//trim(walltime)
    !endif
    ! initialize diag_manager for processor subset output
    !call diag_manager_init(DIAG_MODEL_SUBSET=diag_model_subset, TIME_INIT=date)
    !call print_memuse_stats( 'diag_manager_init' )
    !if ( mpp_pe() == mpp_root_pe()) then
    !   call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
    !   write(errunit,*) 'Finished initializing diag_manager at '&
    !        //trim(walldate)//' '//trim(walltime)
    !endif
    !-----------------------------------------------------------------------
    !------ reset pelist to "full group" ------

    call mpp_set_current_pelist()
    !----- always override initial/base date with diag_manager value -----

    !call get_base_date ( date_init(1), date_init(2), date_init(3), &
    !     date_init(4), date_init(5), date_init(6)  )

    !----- use current date if no base date ------

    if ( date_init(1) == 0 ) date_init = date

    !----- set initial and current time types ------

    Time_init = set_date (date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6))

    Time      = set_date (date(1), date(2), date(3),  &
         date(4), date(5), date(6))

    Time_start = Time

    !----- compute the ending time -----

    Time_end = Time
    do m=1,months
       Time_end = Time_end + set_time(0,days_in_month(Time_end))
    enddo
    Time_end   = Time_end + set_time(hours*3600+minutes*60+seconds, days)
    !Need to pass Time_end into diag_manager for multiple thread case.
    !call diag_manager_set_time_end(Time_end)

    Run_length = Time_end - Time

    !--- get the time that last intermediate restart file was written out.
    !if (file_exist('INPUT/coupler.intermediate.res')) then
    !   call mpp_open(unit,'INPUT/coupler.intermediate.res',action=MPP_RDONLY)
    !   read(unit,*) date_restart
    !   call mpp_close(unit)
    !else
       date_restart = date
    !endif

    Time_restart_current = Time
    if (ALL(restart_interval ==0)) then
       Time_restart = increment_date(Time_end, 0, 0, 10, 0, 0, 0)   ! no intermediate restart
    else
       Time_restart = set_date(date_restart(1), date_restart(2), date_restart(3),  &
            date_restart(4), date_restart(5), date_restart(6) )
       Time_restart = increment_date(Time_restart, restart_interval(1), restart_interval(2), &
            restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
       if (Time_restart <= Time) call mpp_error(FATAL, &
            '==>Error from program coupler: The first intermediate restart time is no larger than the start time')
    endif

    AOGCM%Time = Time
    AOGCM%Time_init = Time_init
    AOGCM%Time_end = Time_end

    !-----------------------------------------------------------------------
    !----- write time stamps (for start time and end time) ------

    call mpp_open( unit, 'time_stamp.out', nohdrs=.TRUE. )

    month = month_name(date(2))
    if ( mpp_pe().EQ.mpp_root_pe() ) write (unit,20) date, month(1:3)

    call get_date (Time_end, date(1), date(2), date(3),  &
         date(4), date(5), date(6))
    month = month_name(date(2))
    if ( mpp_pe().EQ.mpp_root_pe() ) write (unit,20) date, month(1:3)

    call mpp_close(unit)

20  format (i6,5i4,2x,a3)

    !-----------------------------------------------------------------------
    !----- compute the time steps ------

    Time_step_cpld  = set_time (dt_cpld ,0)
    Time_step_atmos = set_time (dt_atmos,0)

    AOGCM%Time_step_cpld = Time_step_cpld
    AOGCM%Time_step_atmos = Time_step_atmos    
    
    !----- determine maximum number of iterations per loop ------

    num_cpld_calls  = Run_length      / Time_step_cpld
    num_atmos_calls = Time_step_cpld  / Time_step_atmos

    !-----------------------------------------------------------------------
    !------------------- some error checks ---------------------------------

    !----- initial time cannot be greater than current time -------

    if ( Time_init > Time ) call error_mesg ('program coupler',  &
         'initial time is greater than current time', FATAL)

    !----- make sure run length is a multiple of ocean time step ------

    if ( num_cpld_calls * Time_step_cpld  /= Run_length )  &
         call error_mesg ('program coupler',  &
         'run length must be multiple of coupled time step', FATAL)

    ! ---- make sure cpld time step is a multiple of atmos time step ----

    if ( num_atmos_calls * Time_step_atmos /= Time_step_cpld )  &
         call error_mesg ('program coupler',   &
         'cpld time step is not a multiple of the atmos time step', FATAL)

    !
    !       Initialize the tracer manager. This needs to be done on all PEs,
    !       before the individual models are initialized.
    !

    !if (mpp_pe().EQ.mpp_root_pe()) then
    !   call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
    !   write(errunit,*) 'Starting to initialize tracer_manager at '&
    !        //trim(walldate)//' '//trim(walltime)
    !endif
    !call tracer_manager_init
    !if (mpp_pe().EQ.mpp_root_pe()) then
    !   call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
    !   write(errunit,*) 'Finished initializing tracer_manager at '&
    !        //trim(walldate)//' '//trim(walltime)
    !endif
    !
    !       Initialize the coupler types
    !

    if (mpp_pe().EQ.mpp_root_pe()) then
       call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
       write(errunit,*) 'Starting to initialize coupler_types at '&
            //trim(walldate)//' '//trim(walltime)
    endif
    call coupler_types_init
    if (mpp_pe().EQ.mpp_root_pe()) then
       call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
       write(errunit,*) 'Finished initializing coupler_types at '&
            //trim(walldate)//' '//trim(walltime)
    endif

    !-----------------------------------------------------------------------
    !------ initialize component models ------
    !------ grid info now comes from grid_spec file

    if (mpp_pe().EQ.mpp_root_pe()) then
       call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
       write(errunit,*) 'Beginning to initialize component models at '&
            //trim(walldate)//' '//trim(walltime)
    endif
    if (AOGCM%Atm%pe) then
       call mpp_set_current_pelist(AOGCM%Atm%pelist)
       !---- atmosphere ----
       if (mpp_pe().EQ.mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize atmospheric model at '&
               //trim(walldate)//' '//trim(walltime)
       endif

       call mpp_clock_begin(id_atmos_model_init)

       call atmos_model_init( AOGCM%Atm, Time_init, Time, Time_step_atmos, &
            do_concurrent_radiation)

       call mpp_clock_end(id_atmos_model_init)

       if (mpp_pe().EQ.mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing atmospheric model at '&
               //trim(walldate)//' '//trim(walltime)
       endif
       call print_memuse_stats( 'atmos_model_init' )
       call data_override_init(Atm_domain_in = AOGCM%Atm%domain)
    endif
    !---- land ----------
    if (AOGCM%Land%pe) then
       call mpp_set_current_pelist(AOGCM%Land%pelist)
!!$       if (mpp_pe().EQ.mpp_root_pe()) then
!!$          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
!!$          write(errunit,*) 'Starting to initialize land model at '&
!!$               //trim(walldate)//' '//trim(walltime)
!!$       endif
!!$       call mpp_clock_begin(id_land_model_init)
!!$       call land_model_init( AOGCM%Atmos_land_boundary, AOGCM%Land, Time_init, Time, &
!!$            Time_step_atmos, Time_step_cpld )
!!$       call mpp_clock_end(id_land_model_init)
!!$       if (mpp_pe().EQ.mpp_root_pe()) then
!!$          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
!!$          write(errunit,*) 'Finished initializing land model at '&
!!$               //trim(walldate)//' '//trim(walltime)
!!$       endif
!!$       call print_memuse_stats( 'land_model_init' )
       !call data_override_init(Land_domain_in = AOGCM%Land%domain)
       !#ifndef _USE_LEGACY_LAND_
       !        call data_override_init(AOGCM%Land_domainUG_in = AOGCM%Land%ug_domain)
       !#endif
    endif
    !---- ice -----------
    if (AOGCM%Ice%pe) then  ! This occurs for all fast or slow ice PEs.
       if (AOGCM%Ice%fast_ice_pe) then
          call mpp_set_current_pelist(AOGCM%Ice%fast_pelist)
       elseif (AOGCM%Ice%slow_ice_pe) then
          call mpp_set_current_pelist(AOGCM%Ice%slow_pelist)
       else
          call mpp_error(FATAL, "All Ice%pes must be a part of Ice%fast_ice_pe or Ice%slow_ice_pe")
       endif
       if (mpp_pe().EQ.mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize ice model at '&
               //trim(walldate)//' '//trim(walltime)
       endif
       call mpp_clock_begin(id_ice_model_init)
       AOGCM%concurrent_ice=concurrent_ice
       call ice_model_init( AOGCM%Ice, Time_init, Time, Time_step_atmos, &
            Time_step_cpld, Verona_coupler=.false., &
            concurrent_ice=concurrent_ice )
       call mpp_clock_end(id_ice_model_init)

       ! This must be called using the union of the ice PE_lists.
       call mpp_set_current_pelist(AOGCM%Ice%pelist)
       call share_ice_domains(AOGCM%Ice)

       if (mpp_pe().EQ.mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing ice model at '&
               //trim(walldate)//' '//trim(walltime)
       endif
       call print_memuse_stats( 'ice_model_init' )
       if (AOGCM%Ice%fast_ice_pe) then
          call mpp_set_current_pelist(AOGCM%Ice%fast_pelist)
          call data_override_init(Ice_domain_in = AOGCM%Ice%domain)
       endif
    endif
    !---- ocean ---------
    if (AOGCM%Ocean%is_ocean_pe) then
       call mpp_set_current_pelist(AOGCM%Ocean%pelist)
       if (mpp_pe().EQ.mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize ocean model at '&
               //trim(walldate)//' '//trim(walltime)
       endif
       call mpp_clock_begin(id_ocean_model_init)
       AOGCM%Ocean_state => NULL()
       call ocean_model_init( AOGCM%Ocean, AOGCM%Ocean_state, Time_init, Time )
       call mpp_clock_end(id_ocean_model_init)

       if (concurrent) then
          !$           call omp_set_num_threads(ocean_nthreads)
          call mpp_set_current_pelist( AOGCM%Ocean%pelist )
          !$           base_cpu = get_cpu_affinity()
          !$OMP PARALLEL private(adder)
          !$        if (use_hyper_thread) then
          !$          if (mod(omp_get_thread_num(),2) == 0) then
          !$            adder = omp_get_thread_num()/2
          !$          else
          !$            adder = ncores_per_node + omp_get_thread_num()/2
          !$          endif
          !$        else
          !$          adder = omp_get_thread_num()
          !$        endif
          !$       call set_cpu_affinity (base_cpu + adder)
          !$        if (debug_affinity) then
          !$          write(6,*) " ocean  ", get_cpu_affinity(), adder, omp_get_thread_num()
          !$          call flush(6)
          !$        endif
          !$OMP END PARALLEL
       else
          ocean_nthreads = atmos_nthreads
          !$        call omp_set_num_threads(ocean_nthreads)
       endif


       if (mpp_pe().EQ.mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing ocean model at '&
               //trim(walldate)//' '//trim(walltime)
       endif
       call print_memuse_stats( 'ocean_model_init' )
       if (mpp_pe().EQ.mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize data_override at '&
               //trim(walldate)//' '//trim(walltime)
       endif
       call data_override_init(Ocean_domain_in = AOGCM%Ocean%domain )
       if (mpp_pe().EQ.mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing data_override at '&
               //trim(walldate)//' '//trim(walltime)
       endif
    endif

    !---------------------------------------------
    if (mpp_pe().EQ.mpp_root_pe()) then
       call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
       write(errunit,*) 'Finished initializing component models at '&
            //trim(walldate)//' '//trim(walltime)
    endif
    call mpp_set_current_pelist(AOGCM%ensemble_pelist(AOGCM%ensemble_id,:))

    call mpp_broadcast_domain(AOGCM%Ice%domain)
    call mpp_broadcast_domain(AOGCM%Ice%slow_domain_NH)
    call mpp_broadcast_domain(AOGCM%Ocean%domain)
    !-----------------------------------------------------------------------
    !---- initialize flux exchange module ----
    if (mpp_pe().EQ.mpp_root_pe()) then
       call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
       write(errunit,*) 'Starting to initialize flux_exchange at '&
            //trim(walldate)//' '//trim(walltime)
    endif
    call mpp_clock_begin(id_flux_exchange_init)
    !call flux_exchange_init ( Time, AOGCM%Atm, AOGCM%Land, AOGCM%Ice, AOGCM%Ocean, AOGCM%Ocean_state,&
    !     AOGCM%atmos_ice_boundary, AOGCM%land_ice_atmos_boundary, &
    !     AOGCM%land_ice_boundary, AOGCM%ice_ocean_boundary, AOGCM%ocean_ice_boundary, &
    !     do_ocean, dt_atmos=dt_atmos, dt_cpld=dt_cpld)
    call mpp_set_current_pelist(AOGCM%ensemble_pelist(AOGCM%ensemble_id,:))
    call mpp_clock_end(id_flux_exchange_init)
    call mpp_set_current_pelist()
    if (mpp_pe().EQ.mpp_root_pe()) then
       call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
       write(errunit,*) 'Finsihed initializing flux_exchange at '&
            //trim(walldate)//' '//trim(walltime)
    endif

    Time_atmos = Time
    Time_ocean = Time

    !
    !       read in extra fields for the air-sea gas fluxes
    !
    if ( AOGCM%Ice%slow_ice_pe ) then
       call mpp_set_current_pelist(AOGCM%Ice%slow_pelist)
       allocate(Ice_bc_restart(AOGCM%Ice%ocean_fluxes%num_bcs))
       allocate(ice_bc_restart_file(AOGCM%Ice%ocean_fluxes%num_bcs))
       do n = 1, AOGCM%Ice%ocean_fluxes%num_bcs  !{
          if (AOGCM%Ice%ocean_fluxes%bc(n)%num_fields .LE. 0) cycle
          filename = trim(AOGCM%Ice%ocean_fluxes%bc(n)%ice_restart_file)
          do l = 1, num_ice_bc_restart
             if (trim(filename) == ice_bc_restart_file(l)) exit
          enddo
          if (l>num_ice_bc_restart) then
             num_ice_bc_restart = num_ice_bc_restart + 1
             ice_bc_restart_file(l) = trim(filename)
          endif
          filename = 'INPUT/'//trim(filename)
          other_fields_exist = .false.
          do m = 1, AOGCM%Ice%ocean_fluxes%bc(n)%num_fields  !{
             fieldname = trim(AOGCM%Ice%ocean_fluxes%bc(n)%field(m)%name)
             id_restart = register_restart_field(Ice_bc_restart(l), ice_bc_restart_file(l), &
                  fieldname, AOGCM%Ice%ocean_fluxes%bc(n)%field(m)%values, AOGCM%Ice%slow_domain_NH   )
             if ( field_exist(filename, fieldname, AOGCM%Ice%slow_domain_NH) ) then
                other_fields_exist = .true.
                write (outunit,*) trim(note_header), ' Reading restart info for ',         &
                     trim(fieldname), ' from ',  trim(filename)
                call read_data(filename, fieldname, AOGCM%Ice%ocean_fluxes%bc(n)%field(m)%values, AOGCM%Ice%slow_domain_NH)
             elseif (other_fields_exist) then
                call mpp_error(FATAL, trim(error_header) // ' Couldn''t find field ' //     &
                     trim(fieldname) // ' in file ' //trim(filename))
             endif
          enddo  !} m
       enddo  !} n
    endif
    if ( AOGCM%Ocean%is_ocean_pe ) then
       call mpp_set_current_pelist(AOGCM%Ocean%pelist)
       allocate(Ocn_bc_restart(AOGCM%Ocean%fields%num_bcs))
       allocate(ocn_bc_restart_file(AOGCM%Ocean%fields%num_bcs))
       do n = 1, AOGCM%Ocean%fields%num_bcs  !{
          if (AOGCM%Ocean%fields%bc(n)%num_fields .LE. 0) cycle
          filename = trim(AOGCM%Ocean%fields%bc(n)%ocean_restart_file)
          do l = 1, num_ocn_bc_restart
             if (trim(filename) == ocn_bc_restart_file(l)) exit
          enddo
          if (l>num_ocn_bc_restart) then
             num_ocn_bc_restart = num_ocn_bc_restart + 1
             ocn_bc_restart_file(l) = trim(filename)
          endif
          filename = 'INPUT/'//trim(filename)
          other_fields_exist = .false.
          do m = 1, AOGCM%Ocean%fields%bc(n)%num_fields  !{
             fieldname = trim(AOGCM%Ocean%fields%bc(n)%field(m)%name)
             id_restart = register_restart_field(Ocn_bc_restart(l), Ocn_bc_restart_file(l), &
                  fieldname, AOGCM%Ocean%fields%bc(n)%field(m)%values, AOGCM%Ocean%domain    )
             if ( field_exist(filename, fieldname, AOGCM%Ocean%domain) ) then
                other_fields_exist = .true.
                write (outunit,*) trim(note_header), ' Reading restart info for ',         &
                     trim(fieldname), ' from ', trim(filename)
                call read_data(filename, fieldname, AOGCM%Ocean%fields%bc(n)%field(m)%values, AOGCM%Ocean%domain)
             elseif (other_fields_exist) then
                call mpp_error(FATAL, trim(error_header) // ' Couldn''t find field ' //     &
                     trim(fieldname) // ' in file ' //trim(filename))
             endif
          enddo  !} m
       enddo  !} n
    endif

    call mpp_set_current_pelist()

    !-----------------------------------------------------------------------
    !---- open and close dummy file in restart dir to check if dir exists --

    !call mpp_open( unit, 'RESTART/file' )
    !call mpp_close(unit, MPP_DELETE)

    ! Call to daig_grid_end to free up memory used during regional
    ! output setup
    !CALL diag_grid_end()

!!$    !-----------------------------------------------------------------------
!!$    if ( do_endpoint_chksum ) then
!!$       if (AOGCM%Atm%pe) then
!!$          call mpp_set_current_pelist(AOGCM%Atm%pelist)
!!$          call atmos_ice_land_chksum('soca_models_init+', 0, AOGCM%Atm, AOGCM%Land, AOGCM%Ice, &
!!$               AOGCM%Land_ice_atmos_boundary, AOGCM%Atmos_ice_boundary, AOGCM%Atmos_land_boundary)
!!$       endif
!!$       if (AOGCM%Ice%slow_ice_PE) then
!!$          call mpp_set_current_pelist(AOGCM%Ice%slow_pelist)
!!$          call slow_ice_chksum('soca_models_init+', 0, AOGCM%Ice, AOGCM%Ocean_ice_boundary)
!!$       endif
!!$       if (AOGCM%Ocean%is_ocean_pe) then
!!$          call mpp_set_current_pelist(AOGCM%Ocean%pelist)
!!$          call ocean_chksum('soca_models_init+', 0, AOGCM%Ocean, AOGCM%Ice_ocean_boundary)
!!$       endif
!!$    endif

    call mpp_set_current_pelist()
    call print_memuse_stats('soca_models_init')
    
    if (mpp_pe().EQ.mpp_root_pe()) then
       call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
       write(errunit,*) 'Exiting soca_models_init at '&
            //trim(walldate)//' '//trim(walltime)
    endif
    !call fms_io_exit
  end subroutine soca_models_init

  !#######################################################################

  subroutine soca_models_end(AOGCM)!Atm, Land, Ice, Ocean, Ocean_state)
    implicit none
    
    type (Coupled)                    :: AOGCM

    type(atmos_land_boundary_type), pointer     :: Atmos_land_boundary
    type(atmos_ice_boundary_type), pointer      :: Atmos_ice_boundary
    type(land_ice_atmos_boundary_type), pointer :: Land_ice_atmos_boundary
    type(land_ice_boundary_type), pointer       :: Land_ice_boundary
    type(ice_ocean_boundary_type), pointer      :: Ice_ocean_boundary
    type(ocean_ice_boundary_type), pointer      :: Ocean_ice_boundary

    !-----------------------------------------------------------------------

    !if ( do_endpoint_chksum ) then
    !   if (Atm%pe) then
    !      call mpp_set_current_pelist(Atm%pelist)
    !      call atmos_ice_land_chksum('soca_models_end', 0, Atm, AOGCM%Land, Ice, &
    !           AOGCM%Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
    !   endif
    !   if (Ice%slow_ice_PE) then
    !      call mpp_set_current_pelist(Ice%slow_pelist)
    !      call slow_ice_chksum('soca_models_end', 0, Ice, AOGCM%Ocean_ice_boundary)
    !   endif
    !   if (AOGCM%Ocean%is_ocean_pe) then
    !      call mpp_set_current_pelist(AOGCM%Ocean%pelist)
    !      call ocean_chksum('soca_models_end', 0, AOGCM%Ocean, Ice_ocean_boundary)
    !   endif
    !endif

    print *,'[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[['
    print *,'[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[['
    print *,'[[[[[[[[[[[[[[[[[[[[[[[[[[     END       [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[['
    print *,'[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[['
    
    call mpp_set_current_pelist()

    !----- check time versus expected ending time ----

    !if (Time /= Time_end) call error_mesg ('program coupler',  &
    !     'final time does not match expected ending time', WARNING)

    !-----------------------------------------------------------------------
    !the call to fms_io_exit has been moved here
    !this will work for serial code or concurrent (disjoint pelists)
    !but will fail on overlapping but unequal pelists
    if (AOGCM%Ocean%is_ocean_pe) then
       call mpp_set_current_pelist(AOGCM%Ocean%pelist)
       call ocean_model_end (AOGCM%Ocean, AOGCM%Ocean_state, AOGCM%Time)
    endif
    if (AOGCM%Atm%pe) then
       call mpp_set_current_pelist(AOGCM%Atm%pelist)
       call atmos_model_end ( AOGCM%Atm )
    endif
    if (AOGCM%Land%pe) then
       call mpp_set_current_pelist(AOGCM%Land%pelist)
       call land_model_end (Atmos_land_boundary, AOGCM%Land)
    endif
    if (AOGCM%Ice%pe) then  ! This happens on all fast or slow ice PEs.
       if (AOGCM%Ice%slow_ice_PE) then
          call mpp_set_current_pelist(AOGCM%Ice%slow_pelist)
       else ! This must be a fast ice PE.
          call mpp_set_current_pelist(AOGCM%Ice%fast_pelist)
       endif
       call ice_model_end (AOGCM%Ice)
    endif

    !----- write restart file ------
    !call soca_coupler_restart(Time, Time_restart_current, &
    !   Atm, AOGCM%Land, Ice, AOGCM%Ocean, AOGCM%Ocean_state)

    !call diag_manager_end (AOGCM%Time)

    call mpp_set_current_pelist()

    !-----------------------------------------------------------------------

  end subroutine soca_models_end

!!$  !> \brief Writing restart file that contains running time and restart file writing time.
!!$  subroutine soca_coupler_restart(Time_run, Time_res, time_stamp, &
!!$       Atm, Land, Ice, Ocean, Ocean_state)
!!$    type(time_type),   intent(in)           :: Time_run, Time_res
!!$    character(len=*), intent(in),  optional :: time_stamp
!!$    type (atmos_data_type) :: Atm
!!$    type  (land_data_type) :: Land
!!$    type   (ice_data_type) :: Ice
!!$    type (ocean_public_type):: Ocean
!!$    type (ocean_state_type),  pointer :: Ocean_state
!!$
!!$    character(len=128)                      :: file_run, file_res
!!$    integer :: yr, mon, day, hr, min, sec, date(6), unit, n
!!$
!!$    call mpp_set_current_pelist()
!!$
!!$    ! write restart file
!!$    if (present(time_stamp)) then
!!$       file_run = 'RESTART/'//trim(time_stamp)//'.coupler.res'
!!$       file_res = 'RESTART/'//trim(time_stamp)//'.coupler.intermediate.res'
!!$    else
!!$       file_run = 'RESTART/coupler.res'
!!$       file_res = 'RESTART/coupler.intermediate.res'
!!$    endif
!!$
!!$    !----- compute current date ------
!!$    call get_date (Time_run, date(1), date(2), date(3),  &
!!$         date(4), date(5), date(6))
!!$    call mpp_open( unit, file_run, nohdrs=.TRUE. )
!!$    if ( mpp_pe().EQ.mpp_root_pe()) then
!!$       write( unit, '(i6,8x,a)' )calendar_type, &
!!$            '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'
!!$
!!$       write( unit, '(6i6,8x,a)' )date_init, &
!!$            'Model start time:   year, month, day, hour, minute, second'
!!$       write( unit, '(6i6,8x,a)' )date, &
!!$            'Current model time: year, month, day, hour, minute, second'
!!$    endif
!!$    call mpp_close(unit)
!!$
!!$    if (Time_res > Time_start) then
!!$       call mpp_open( unit, file_res, nohdrs=.TRUE. )
!!$       if ( mpp_pe().EQ.mpp_root_pe()) then
!!$          call get_date(Time_res ,yr,mon,day,hr,min,sec)
!!$          write( unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
!!$               'Current intermediate restart time: year, month, day, hour, minute, second'
!!$       endif
!!$       call mpp_close(unit)
!!$    endif
!!$
!!$    if (Ocean%is_ocean_pe) then
!!$       call mpp_set_current_pelist(Ocean%pelist)
!!$       do n = 1, num_ocn_bc_restart
!!$          call save_restart(Ocn_bc_restart(n), time_stamp)
!!$       enddo
!!$    endif
!!$    if (Atm%pe) then
!!$       call mpp_set_current_pelist(Atm%pelist)
!!$       do n = 1, num_ice_bc_restart
!!$          call save_restart(Ice_bc_restart(n), time_stamp)
!!$       enddo
!!$    endif
!!$
!!$
!!$  end subroutine soca_coupler_restart

  !--------------------------------------------------------------------------

  !> \brief Print out checksums for several atm, land and ice variables
  subroutine coupler_chksum(id, timestep, Atm, Land, Ice, Ocean, Ocean_state)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type (atmos_data_type) :: Atm
    type  (land_data_type) :: Land
    type   (ice_data_type) :: Ice
    type (ocean_public_type):: Ocean
    type (ocean_state_type),  pointer :: Ocean_state

    type :: tracer_ind_type
       integer :: atm, ice, lnd ! indices of the tracer in the respective models
    end type tracer_ind_type
    integer                            :: n_atm_tr, n_lnd_tr, n_exch_tr
    integer                            :: n_atm_tr_tot, n_lnd_tr_tot
    integer                            :: i, tr, n, m, outunit
    type(tracer_ind_type), allocatable :: tr_table(:)
    character(32) :: tr_name


    call get_number_tracers (MODEL_ATMOS, num_tracers=n_atm_tr_tot, &
         num_prog=n_atm_tr)
    call get_number_tracers (MODEL_LAND, num_tracers=n_lnd_tr_tot, &
         num_prog=n_lnd_tr)

    ! Assemble the table of tracer number translation by matching names of
    ! prognostic tracers in the atmosphere and surface models; skip all atmos.
    ! tracers that have no corresponding surface tracers.
    allocate(tr_table(n_atm_tr))
    n = 1
    do i = 1,n_atm_tr
       call get_tracer_names( MODEL_ATMOS, i, tr_name )
       tr_table(n)%atm = i
       tr_table(n)%ice = get_tracer_index ( MODEL_ICE,  tr_name )
       tr_table(n)%lnd = get_tracer_index ( MODEL_LAND, tr_name )
       if (tr_table(n)%ice/=NO_TRACER .or. tr_table(n)%lnd/=NO_TRACER) n = n+1
    enddo
    n_exch_tr = n-1

100 FORMAT("CHECKSUM::",A32," = ",Z20)
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

    if (Atm%pe) then
       call mpp_set_current_pelist(Atm%pelist)

       outunit = stdout()
       write(outunit,*) 'BEGIN CHECKSUM(Atm):: ', id, timestep
       write(outunit,100) 'atm%t_bot', mpp_chksum(atm%t_bot)
       write(outunit,100) 'atm%z_bot', mpp_chksum(atm%z_bot)
       write(outunit,100) 'atm%p_bot', mpp_chksum(atm%p_bot)
       write(outunit,100) 'atm%u_bot', mpp_chksum(atm%u_bot)
       write(outunit,100) 'atm%v_bot', mpp_chksum(atm%v_bot)
       write(outunit,100) 'atm%p_surf', mpp_chksum(atm%p_surf)
       write(outunit,100) 'atm%gust', mpp_chksum(atm%gust)
       do tr = 1,n_exch_tr
          n = tr_table(tr)%atm
          if (n /= NO_TRACER) then
             call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
             write(outunit,100) 'atm%'//trim(tr_name), mpp_chksum(Atm%tr_bot(:,:,n))
          endif
       enddo

       write(outunit,100) 'land%t_surf', mpp_chksum(land%t_surf)
       write(outunit,100) 'land%t_ca', mpp_chksum(land%t_ca)
       write(outunit,100) 'land%rough_mom', mpp_chksum(land%rough_mom)
       write(outunit,100) 'land%rough_heat', mpp_chksum(land%rough_heat)
       write(outunit,100) 'land%rough_scale', mpp_chksum(land%rough_scale)
       do tr = 1,n_exch_tr
          n = tr_table(tr)%lnd
          if (n /= NO_TRACER) then
             call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
             !#ifndef _USE_LEGACY_LAND_
             !             write(outunit,100) 'land%'//trim(tr_name), mpp_chksum(Land%tr(:,:,n))
             !#else
             write(outunit,100) 'land%'//trim(tr_name), mpp_chksum(Land%tr(:,:,:,n))
             !#endif
          endif
       enddo

       write(outunit,100) 'ice%t_surf', mpp_chksum(ice%t_surf)
       write(outunit,100) 'ice%rough_mom', mpp_chksum(ice%rough_mom)
       write(outunit,100) 'ice%rough_heat', mpp_chksum(ice%rough_heat)
       write(outunit,100) 'ice%rough_moist', mpp_chksum(ice%rough_moist)
       write(outunit,*) 'STOP CHECKSUM(Atm):: ', id, timestep

       !endif

       !if (Ocean%is_ocean_pe) then
       !call mpp_set_current_pelist(Ocean%pelist)

       write(outunit,*) 'BEGIN CHECKSUM(Ice):: ', id, timestep
       do n = 1, ice%ocean_fields%num_bcs  !{
          do m = 1, ice%ocean_fields%bc(n)%num_fields  !{
             !write(outunit,101) 'ice%', m, n, mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
             write(outunit,101) 'ice%',trim(ice%ocean_fields%bc(n)%name), &
                  trim(ice%ocean_fields%bc(n)%field(m)%name), mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
          enddo  !} m
       enddo  !} n
       write(outunit,*) 'STOP CHECKSUM(Ice):: ', id, timestep

    endif

    deallocate(tr_table)

    call mpp_set_current_pelist()

  end subroutine coupler_chksum

  !#######################################################################

  !> \brief This subroutine calls subroutine that will print out checksums of the elements
  !! of the appropriate type.
  !!
  !! For coupled models typically these types are not defined on all processors.
  !! It is assumed that the appropriate pelist has been set before entering this routine.
  !! This can be achieved in the following way.
  !! ~~~~~~~~~~{.f90}
  !! if (Atm%pe) then
  !!    call mpp_set_current_pelist(Atm%pelist)
  !!    call atmos_ice_land_chksum('MAIN_LOOP-', nc)
  !! endif
  !! ~~~~~~~~~~
  !! If you are on the global pelist before you enter this routine using the above call,
  !! you can return to the global pelist by invoking
  !! ~~~~~~~~~~{.f90}
  !! call mpp_set_current_pelist()
  !! ~~~~~~~~~~
  !! after you exit. This is only necessary if you need to return to the global pelist.
  subroutine atmos_ice_land_chksum(id, timestep, Atm, Land, Ice, &
       Land_ice_atmos_boundary, Atmos_ice_boundary, &
       Atmos_land_boundary)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type (atmos_data_type), intent(in) :: Atm
    type  (land_data_type), intent(in) :: Land
    type   (ice_data_type), intent(in) :: Ice
    type(land_ice_atmos_boundary_type), intent(in) :: Land_ice_atmos_boundary
    type(atmos_ice_boundary_type), intent(in)      :: Atmos_ice_boundary
    type(atmos_land_boundary_type), intent(in)     :: Atmos_land_boundary

    call atmos_data_type_chksum(     id, timestep, Atm)
    call lnd_ice_atm_bnd_type_chksum(id, timestep, Land_ice_atmos_boundary)

    if (Ice%fast_ice_pe) then
       call mpp_set_current_pelist(Ice%fast_pelist)
       call ice_data_type_chksum(   id, timestep, Ice)
       call atm_ice_bnd_type_chksum(id, timestep, Atmos_ice_boundary)
    endif
    if (Land%pe) then
       call mpp_set_current_pelist(Land%pelist)
       call land_data_type_chksum(  id, timestep, Land)
       call atm_lnd_bnd_type_chksum(id, timestep, Atmos_land_boundary)
    endif

    call mpp_set_current_pelist(Atm%pelist)

  end subroutine atmos_ice_land_chksum

  !> \brief This subroutine calls subroutine that will print out checksums of the elements
  !! of the appropriate type.
  !!
  !! For coupled models typically these types are not defined on all processors.
  !! It is assumed that the appropriate pelist has been set before entering this routine.
  !! This can be achieved in the following way.
  !! ~~~~~~~~~~{.f90}
  !! if (Ice%slow_ice_pe) then
  !!    call mpp_set_current_pelist(Ice%slow_pelist)
  !!    call slow_ice_chksum('MAIN_LOOP-', nc)
  !! endif
  !! ~~~~~~~~~~
  !! If you are on the global pelist before you enter this routine using the above call,
  !! you can return to the global pelist by invoking
  !! ~~~~~~~~~~{.f90}
  !! call mpp_set_current_pelist()
  !! ~~~~~~~~~~
  !! after you exit. This is only necessary if you need to return to the global pelist.
  subroutine slow_ice_chksum(id, timestep, Ice, Ocean_ice_boundary)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_data_type), intent(in) :: Ice
    type(ocean_ice_boundary_type), intent(in) :: Ocean_ice_boundary

    call ice_data_type_chksum(    id, timestep, Ice)
    call ocn_ice_bnd_type_chksum( id, timestep, Ocean_ice_boundary)

  end subroutine slow_ice_chksum


  !> \brief This subroutine calls subroutine that will print out checksums of the elements
  !! of the appropriate type.
  !!
  !! For coupled models typically these types are not defined on all processors.
  !! It is assumed that the appropriate pelist has been set before entering this routine.
  !! This can be achieved in the following way.
  !! ~~~~~~~~~~{.f90}
  !! if (Ocean%is_ocean_pe) then
  !!    call mpp_set_current_pelist(Ocean%pelist)
  !!    call ocean_chksum('MAIN_LOOP-', nc)
  !! endif
  !! ~~~~~~~~~~
  !! If you are on the global pelist before you enter this routine using the above call,
  !! you can return to the global pelist by invoking
  !! ~~~~~~~~~~{.f90}
  !! call mpp_set_current_pelist()
  !! ~~~~~~~~~~
  !! after you exit. This is only necessary if you need to return to the global pelist.
  subroutine ocean_chksum(id, timestep, Ocean, Ice_ocean_boundary)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type (ocean_public_type), intent(in) :: Ocean
    type(ice_ocean_boundary_type), intent(in) :: Ice_ocean_boundary

    call ocean_public_type_chksum(id, timestep, Ocean)
    call ice_ocn_bnd_type_chksum( id, timestep, Ice_ocean_boundary)

  end subroutine ocean_chksum


end module soca_mom6sis2
