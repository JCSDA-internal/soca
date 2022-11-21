! (C) Copyright 2020 NOAA
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#define esmf_err_abort(rc) if (esmf_LogFoundError(rc, msg="Aborting UFS", line=__LINE__, file=__FILE__)) call esmf_Finalize(endflag=esmf_END_ABORT)

module soca_ufs_mod

  ! oops
  use datetime_mod
  use duration_mod

  ! fckit
  use fckit_configuration_module, only: fckit_configuration

  ! soca
  use soca_geom_mod,      only: soca_geom
  use soca_state_mod,     only: soca_state
  use soca_fields_mod,     only: soca_field

  ! ufs
  use ESMF
!use ESMF,  only: ESMF_FieldRedistStore, ESMF_FieldRedist, ESMF_RouteHandle
  use NUOPC
  use NUOPC_Driver
  !use module_EARTH_GRID_COMP, only: esmSS => EARTH_REGISTER
  use UFSDriver, only: esmSS => UFSDriver_SS
  use mpp_mod,            only: read_input_nml,mpp_pe


  implicit none
  private

  public :: model_ufs

  !> Fortran derived type to hold model definition
  type :: model_ufs
     type(ESMF_GridComp) :: esmComp
     type(ESMF_State) :: toJedi, fromJedi
     integer :: isc, iec, jsc, jec, npz
     type(esmf_Clock) :: clock
     type(esmf_config) :: cf_main                                         !<-- the configure object
   contains
     procedure :: create
     procedure :: delete
     procedure :: initialize
     procedure :: step
     procedure :: finalize
  end type model_ufs

  character(len=*), parameter :: modname='soca_ufs_mod'

  ! --------------------------------------------------------------------------------------------------

contains

  ! --------------------------------------------------------------------------------------------------

  subroutine create(self, conf, geom)

    implicit none
    class(model_ufs),          intent(inout) :: self
    type(fckit_configuration), intent(in)    :: conf
    type(soca_geom),        intent(in)    :: geom

    integer :: rc, urc, phase, i, cnt
    character(len=20) :: cdate_start, cdate_stop

    type(ESMF_Time)         :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep

    character(len=*),parameter :: subname = modname//' (create)'
    type(ESMF_CplComp),  pointer       :: connectors(:)
    character(len=128) :: name, msg

    ! Initialize ESMF
    call ESMF_Initialize(logkindflag=esmf_LOGKIND_MULTI, &
         defaultCalkind=esmf_CALKIND_GREGORIAN, &
         mpiCommunicator=geom%f_comm%communicator(), rc=rc)
    esmf_err_abort(rc)
    ! Flush log output while debugging
    call ESMF_LogSet(flush=.true., rc=rc)
    esmf_err_abort(rc)
    call ESMF_LogWrite("ESMF Initialized in "//subname, ESMF_LOGMSG_INFO)

    self%isc = geom%isc
    self%iec = geom%iec
    self%jsc = geom%jsc
    self%jec = geom%jec
    self%npz = geom%nzo
    self%cf_main=esmf_configcreate(rc=rc)
    call ESMF_ConfigLoadFile(config=self%cf_main, &
         filename='model_configure', &
         rc=rc)

    ! This call to read_input_nml() seems to be required
    ! for CCPP.  However, it does not belong at this level
    ! but should be handled inside the model itself
    call read_input_nml()
    call ESMF_LogWrite("done reading input nml", ESMF_LOGMSG_INFO)

    ! Create the ESM component
    self%esmComp = ESMF_GridCompCreate(name="esm", rc=rc)
    esmf_err_abort(rc)

    ! SetServices for the ESM component
    call ESMF_GridCompSetServices(self%esmComp, esmSS, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)


    ! Set ESM's Verbosity (High)  - 32513
    call NUOPC_CompAttributeSet(self%esmComp, name="Verbosity", &
         value="32513", rc=rc)
    esmf_err_abort(rc)



    ! Initialize the clock based on contents of model_configure
    ! -------------------------------------------
    call setUFSClock(self,startTime,stopTime)
    call ESMF_GridCompSet(self%esmComp, clock=self%clock, rc=rc)

    ! Create import and export states from
    ! perspective of the exernal system:
    !   toJedi is an IMPORT into Jedi and an EXPORT from ESM
    !   fromJedi is an EXPORT from Jedi and an IMPORT into ESM

    self%toJedi = ESMF_StateCreate(stateintent=ESMF_STATEINTENT_IMPORT, &
         rc=rc)
    esmf_err_abort(rc)


    self%fromJedi = ESMF_StateCreate(stateintent=ESMF_STATEINTENT_EXPORT, &
         rc=rc)
    esmf_err_abort(rc)


    call ESMF_LogWrite("Advertising export from ESM", ESMF_LOGMSG_INFO)
    ! Advertise fields on the exportState, for data coming out of ESM component
    ! Note--only certain fields are available. Check in GFS_surface_generic to see if they are filled
#if 1
    call NUOPC_Advertise(self%toJedi, &
         StandardNames=(/ &
         !              "tocn                                 "/), &   ! Example fields
                        "socn                                 ", &   ! Example fields
                        "tocn                                 ", &   ! Example fields
                        "uocn                                 ", &   ! Example fields
                        "vocn                                 "/), &   ! Example fields
         SharePolicyField="share", &
         TransferOfferGeomObject="cannot provide", rc=rc)
    esmf_err_abort(rc)

#endif
    call ESMF_LogWrite("Advertising imports to ESM", ESMF_LOGMSG_INFO)
    ! Advertise fields on the importState, for data going into ESM componenta

! imports are not yet implemented
#if 0
    call NUOPC_Advertise(self%fromJedi, &
         StandardNames=(/ &
                        "u                                    ", &   ! Example fields
                        "v                                    ", &   ! Example fields
                        "ua                                   ", &   ! Example fields
                        "va                                   ", &   ! Example fields
                        "t                                    ", &   ! Example fields
                        "delp                                 ", &   ! Example fields
                        "sphum                                ", &   ! Example fields
                        "ice_wat                              ", &   ! Example fields
                        "liq_wat                              ", &   ! Example fields
                        "o3mr                                 ", &   ! Example fields
                        "phis                                 ", &   ! Example fields
                        "slmsk                                ", &   ! Example fields
                        "weasd                                ", &   ! Example fields
                        "tsea                                 ", &   ! Example fields
                        "vtype                                ", &   ! Example fields
                        "stype                                ", &   ! Example fields
                        "vfrac                                ", &   ! Example fields
                        "stc                                  ", &   ! Example fields
                        "smc                                  ", &   ! Example fields
                        "snwdph                               ", &   ! Example fields
                        "u_srf                                ", &   ! Example fields
                        "v_srf                                ", &   ! Example fields
                        "f10m                                 "/), &   ! Example fields
         TransferOfferGeomObject="cannot provide", rc=rc)

    esmf_err_abort(rc)
#endif

    call ESMF_StateGet(self%toJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)

    write(msg, "(I2)") cnt
    call ESMF_LogWrite("After filling advertise toJedi state has "//trim(msg)//" items.", &
         ESMF_LOGMSG_INFO)


    ! call ExternalAdvertise phase
    call NUOPC_CompSearchPhaseMap(self%esmComp, &
         methodflag=ESMF_METHOD_INITIALIZE, &
         phaseLabel=label_ExternalAdvertise, phaseIndex=phase, rc=rc)
    esmf_err_abort(rc)

    call ESMF_LogWrite("done with CompSearchPhaseMap ", &
         ESMF_LOGMSG_INFO)

    call ESMF_GridCompInitialize(self%esmComp, phase=phase, &
         importState=self%fromJedi, exportState=self%toJedi, &
         clock=self%clock, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)

    call ESMF_LogWrite("done with GridCompInitialize ", &
         ESMF_LOGMSG_INFO)

    call ESMF_StateGet(self%toJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)

    write(msg, "(I2)") cnt
    call ESMF_LogWrite("After calling advertise toJedi state has "//trim(msg)//" items.", &
         ESMF_LOGMSG_INFO)


    ! Set verbosity flag on connectors
    nullify(connectors);
    call NUOPC_DriverGetComp(self%esmComp, &
         compList=connectors, rc=rc)
    esmf_err_abort(rc)


    call ESMF_LogWrite("About to set connector verbosity", ESMF_LOGMSG_INFO)
    do i=lbound(connectors,1), ubound(connectors,1)
       call ESMF_CplCompGet(connectors(i), name=name, rc=rc)
       esmf_err_abort(rc)

       call NUOPC_CompAttributeSet(connectors(i), name="Verbosity", &
            value="max", rc=rc)
       esmf_err_abort(rc)

       call ESMF_LogWrite(" --> Set verbosity on connector: "//trim(name), &
            ESMF_LOGMSG_INFO)
    enddo

    deallocate(connectors)

    ! call ExternalRealize phase
    call NUOPC_CompSearchPhaseMap(self%esmComp, &
         methodflag=ESMF_METHOD_INITIALIZE, &
         phaseLabel=label_ExternalRealize, phaseIndex=phase, rc=rc)
    esmf_err_abort(rc)

    call ESMF_GridCompInitialize(self%esmComp, phase=phase, &
         importState=self%fromJedi, exportState=self%toJedi, &
         clock=self%clock, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)

    call ESMF_StateGet(self%toJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)

    write(msg, "(I2)") cnt

    call ESMF_LogWrite("Dumping toJedi state with "//trim(msg)//" items", &
         ESMF_LOGMSG_INFO)

    ! call ExternalDataInit phase
    call NUOPC_CompSearchPhaseMap(self%esmComp, &
         methodflag=ESMF_METHOD_INITIALIZE, &
         phaseLabel=label_ExternalDataInit, phaseIndex=phase, rc=rc)
    esmf_err_abort(rc)

    call ESMF_GridCompInitialize(self%esmComp, phase=phase, &
         importState=self%fromJedi, exportState=self%toJedi, &
         clock=self%clock, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)

    call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)

  end subroutine create

! --------------------------------------------------------------------------------------------------

  subroutine initialize(self, state, vdate)

    implicit none

    class(model_ufs),    intent(inout) :: self
    type(soca_state), intent(in)    :: state
    type(datetime),      intent(in)    :: vdate

    character(len=*),parameter :: subname = modname//' (initialize)'

    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)

  end subroutine initialize

! --------------------------------------------------------------------------------------------------

  subroutine step(self, state, vdate_start, vdate_final)

    implicit none

    class(model_ufs),    intent(inout) :: self
    type(soca_state), intent(inout) :: state
    type(datetime),      intent(in)    :: vdate_start
    type(datetime),      intent(in)    :: vdate_final

    ! local variables
    integer :: rc, urc, cnt
    character(len=20) :: strStartTime, strStopTime
    character(len=128) :: name, msg
    type(ESMF_Time) :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep
    integer, save     :: tstep=1
    character(len=80) :: fileName

!-----------------------------------------------------------------------------

    character(len=*),parameter :: subname = modname//' (step)'
    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    call datetime_to_string(vdate_start, strStartTime)
    call datetime_to_string(vdate_final, strStopTime)

    call ESMF_LogWrite(" --> REQUESTED START TIME:"//trim(strStartTime), ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(" --> REQUESTED STOP  TIME:"//trim(strStopTime), ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(self%clock, startTime=startTime, &
         stopTime=stopTime, rc=rc)
    esmf_err_abort(rc)


    call ESMF_TimeSet(startTime, timeString=strStartTime, rc=rc)
    esmf_err_abort(rc)


    call ESMF_TimeSet(stopTime, timeString=strStopTime, rc=rc)
    esmf_err_abort(rc)


    timeStep = stopTime - startTime

    call ESMF_ClockSet(self%clock, startTime=startTime, &
         stopTime=stopTime, currTime=startTime, timeStep=timeStep, rc=rc)
     esmf_err_abort(rc)

    ! step the model forward
    call ESMF_GridCompRun(self%esmComp, &
         importState=self%fromJedi, exportState=self%toJedi, &
         clock=self%clock, userRc=urc, rc=rc)
    esmf_err_abort(rc)
    esmf_err_abort(urc)

    call ESMF_StateGet(self%toJedi, itemCount=cnt, rc=rc)
    esmf_err_abort(rc)
    write(msg, "(I2)") cnt
    call ESMF_LogWrite("after step toJedi state with "//trim(msg)//" items", &
         ESMF_LOGMSG_INFO)
    write(fileName, '("fields_in_esm_import_step",I2.2,".nc")') tstep
    call mom6_to_state(self, state)
    call ESMF_LogWrite("after state write "//trim(msg)//" rc", &
         ESMF_LOGMSG_INFO)

    call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)

  end subroutine step

! --------------------------------------------------------------------------------------------------

  subroutine delete(self)

    implicit none
    class(model_ufs), intent(inout) :: self
    integer :: rc
    character(len=*),parameter :: subname = modname//' (delete)'

    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    call ESMF_GridCompDestroy(self%esmComp, rc=rc)
    esmf_err_abort(rc)


    call ESMF_LogWrite("About to destroy toJedi state "//subname, ESMF_LOGMSG_INFO)

    call ESMF_StateDestroy(self%toJedi, rc=rc)
    esmf_err_abort(rc)


    call ESMF_LogWrite("About to destroy fromJedi state "//subname, ESMF_LOGMSG_INFO)

    call ESMF_StateDestroy(self%fromJedi, rc=rc)
    esmf_err_abort(rc)

    ! Finalize ESMF
    ! -------------
    call ESMF_Finalize(endflag=ESMF_END_KEEPMPI, rc=rc)
    if (rc /= ESMF_SUCCESS) then
       call ESMF_LogWrite("ERROR FINALIZING ESMF "//subname, ESMF_LOGMSG_INFO)
    else
       call ESMF_LogWrite("SUCCESSFULLY FINALIZED ESMF"//subname, ESMF_LOGMSG_INFO)
    endif
  end subroutine delete

  ! --------------------------------------------------------------------------------------------------

  subroutine finalize(self, state, vdate)

    implicit none
    class(model_ufs),    intent(inout) :: self
    type(soca_state), intent(in)    :: state
    type(datetime),intent(in)    :: vdate

    character(len=*),parameter :: subname = modname//' (finalize)'
    ! Clean up is being done in the delete method
    call ESMF_LogWrite("Enter "//subname, ESMF_LOGMSG_INFO)

    call ESMF_LogWrite("Exit "//subname, ESMF_LOGMSG_INFO)

  end subroutine finalize

  subroutine mom6_to_state( self, state )

  implicit none
  type(model_ufs),    intent(in)    :: self
  type(soca_state), intent(inout) :: state

  integer :: num_items, i, rc, rank, lb(3), ub(3), fnpz, ib, nb
  type(ESMF_Field) :: field, destField
  character(len=ESMF_MAXSTR), allocatable :: item_names(:)
  real(kind=ESMF_KIND_R8), pointer :: farrayPtr2(:,:)
  real(kind=ESMF_KIND_R8), pointer :: farrayPtr3(:,:,:)
  character(len=80) :: soca_name
  type(soca_field), pointer :: field_ptr
  type(ESMF_Distgrid)                        :: distgrid
  type(ESMF_Grid)                            :: grid
  type(ESMF_RouteHandle)                      :: routehandle

  real(kind=ESMF_KIND_R8),allocatable,dimension(:,:,:)      :: field_mom6

  distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/360,320/), &
            regDecomp=(/2,4/), &
            rc=rc)
  grid = ESMF_GridCreate(distgrid=distgrid, &
        name="grid", rc=rc)
     if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Array to hold output from UFS in JEDI precision
  ! ------------------------------------------------
  ib = (self%isc-1)*(self%jsc-1)+1
  nb = (self%iec - self%isc + 1) * (self%jec - self%jsc + 1)
  allocate(field_mom6(self%isc:self%iec,self%jsc:self%jec, self%npz))


  ! Get number of items
  ! -------------------
  call ESMF_StateGet(self%toJedi, itemcount = num_items, rc = rc)
  if (rc.ne.0) call abor1_ftn("mom6_to_state: ESMF_StateGet itemcount failed")


  ! Get names of the items
  ! ----------------------
  allocate(item_names(num_items))
  call ESMF_StateGet(self%toJedi, itemnamelist = item_names, rc = rc)
  if (rc.ne.0) call abor1_ftn("mom6_to_state: ESMF_StateGet itemnamelist failed")


  ! Loop over states coming from UFS and convert to JEDI state
  ! -----------------------------------------------------------
  do i = 1, num_items
    ! Create map between UFS name and mom6-jedi name
    ! ----------------------------------------------
    soca_name = trim(item_names(i))
    call ESMF_LogWrite("item name is "//soca_name, ESMF_LOGMSG_INFO)

    ! Only need to extract field from UFS if mom6-jedi needs it
    ! ---------------------------------------------------------
#if 1
!   if (state%has_field(trim(soca_name))) then

      !Get field from the state
      call ESMF_StateGet(self%toJedi, item_names(i), field, rc = rc)
      if (rc.ne.0) call abor1_ftn("mom6_to_state: ESMF_StateGet field failed")

      !Validate the field
      call ESMF_FieldValidate(field, rc = rc)
      if (rc.ne.0) call abor1_ftn("mom6_to_state: ESMF_FieldValidate failed")

      !Get the field rank
      call ESMF_FieldGet(field, rank = rank, rc = rc)
      if (rc.ne.0) call abor1_ftn("mom6_to_state: ESMF_FieldGet rank failed")

      !Convert field to pointer and pointer bounds
      field_mom6 = 0.0_ESMF_KIND_R8
      if (rank == 2) then

        ! this field is going to be on a 1D mesh with an ungridded dimension
        call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr2, totalLBound = lb(1:2), totalUBound = ub(1:2), rc = rc )

        ! create a gridded field to redistribute the data to
        destField = ESMF_FieldCreate(grid=grid, typekind=ESMF_TYPEKIND_R8, &
                 ungriddedLBound=(/1/), ungriddedUBound=(/75/),name='tocn', rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

        ! 1. setup routehandle from source Field to destination Field
        call ESMF_FieldRedistStore(field, destField, routehandle, rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

        ! 2. use precomputed routehandle to redistribute data
        call ESMF_FieldRedist(field, destField, routehandle, rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

        ! verify redist
        call ESMF_FieldGet(destField, localDe=0, farrayPtr=farrayPtr3, totalLBound = lb(1:3), totalUBound = ub(1:3), rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
             
    
        if (rc.ne.0) call abor1_ftn("mom6_to_state: ESMF_FieldGet 2D failed")

        fnpz = ub(3)-lb(3)+1
        field_mom6(self%isc:self%iec,self%jsc:self%jec,1:fnpz) = farrayPtr3(lb(1):ub(1),lb(2):ub(2),:)
        nullify(farrayPtr3)

!     elseif (rank == 3) then
!       call ESMF_FieldGet( field, 0, farrayPtr = farrayPtr3, totalLBound = lb, totalUBound = ub, rc = rc )
!       if (rc.ne.0) call abor1_ftn("mom6_to_state: ESMF_FieldGet 3D failed",rc)

!       fnpz = ub(3)-lb(3)+1
!       field_mom6(self%isc:self%iec,self%jsc:self%jec,1:fnpz) = farrayPtr3(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3))
!       nullify(farrayPtr3)

      else
        deallocate(item_names)
        return 
        call abor1_ftn("mom6_mod: can only handle rank 2 or rank 3 fields from UFS")

      endif
  if(ub(1)<0) then
      ! Check that dimensions match
      if ((ub(1)-lb(1)+1 .ne. self%iec-self%isc+1) .or. (ub(2)-lb(2)+1 .ne. self%jec-self%jsc+1) ) then
        call abor1_ftn("mom6_to_state: dimension mismatch between JEDI and UFS horizontal grid")
      endif

  endif
      ! Get pointer to mom6-jedi side field
      call state%get(trim(soca_name), field_ptr)

      if (field_ptr%nz .ne. fnpz) &
        call abor1_ftn("mom6_to_state: dimension mismatch between JEDI and UFS vertical grid")
      ! Copy from UFS to mom6-jedi
!     field_ptr%val(self%isc:self%iec,self%jsc:self%jec,1:fnpz) = field_mom6(self%isc:self%iec,self%jsc:self%jec,1:fnpz)
!   else
!     call ESMF_LogWrite("Not needed by JEDI is "//soca_name, ESMF_LOGMSG_INFO)
!   endif
#endif
  end do

  deallocate(item_names)

  end subroutine mom6_to_state


  subroutine setUFSClock(self,date_start,date_final)

    class(model_ufs),    intent(inout) :: self
    type(esmf_time),     intent(inout) :: date_start, date_final

    type(ESMF_TimeInterval)            :: runDuration, timestep
    integer :: yy,mm,dd,hh,mns,sec,rc
    !-----------------------------------------------------------------------
    !***  extract the start time from the configuration
    !-----------------------------------------------------------------------
    !
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =YY                            &
                            ,label ='start_year:'                 &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =MM                            &
                            ,label ='start_month:'                &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =DD                            &
                            ,label ='start_day:'                  &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =HH                            &
                            ,label ='start_hour:'                 &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =MNS                           &
                            ,label ='start_minute:'               &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    call ESMF_ConfigGetAttribute(config=self%cf_main                  &
                            ,value =SEC                           &
                            ,label ='start_second:'               &
                            ,rc    =RC)
    esmf_err_abort(rc)

    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    ! Set date_start
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    !
    call ESMF_TimeSet(time=date_start                                 &  !<-- The start time of the forecast (ESMF)
                 ,yy  =YY                                         &  !<-- Year from config file
                 ,mm  =MM                                         &  !<-- Month from config file
                 ,dd  =DD                                         &  !<-- Day from config file
                 ,h   =HH                                         &  !<-- Hour from config file
                 ,m   =MNS                                        &  !<-- Minute from config file
                 ,s   =SEC                                        &  !<-- Second from config file
                ,rc  =RC)
    esmf_err_abort(rc)

    ! Set any time interval here It will be overridden later
    call ESMF_timeintervalset(runduration, d=1, rc=rc)
    call ESMF_timeintervalset(timestep, h=1, rc=rc)
    date_final = date_start + runduration

    self%clock = ESMF_ClockCreate(name="main_clock", &
         timeStep=timestep, startTime=date_start, stopTime=date_final, rc=rc)

    esmf_err_abort(rc)

    call ESMF_ClockPrint(self%clock, options="startTime", &
    preString="Printing startTime to stdout: ", rc=rc)


    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  end subroutine setUFSClock

end module soca_ufs_mod
