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
  use fckit_log_module, only : fckit_log, log
  use soca_geom_mod_c
  use soca_mom6
  use soca_utils
  use soca_fields
  use datetime_mod  
  use MOM,  only : step_MOM
  use MOM_restart, only : save_restart
  use MOM_time_manager,    only : real_to_time, time_type_to_real
  use time_manager_mod, only : time_type, print_time, print_date, set_date
  use MOM_time_manager,    only : operator(+)
  use mpp_domains_mod, only : mpp_update_domains
  
  implicit none

  private
  public :: soca_model
  public :: soca_model_registry
  public :: soca_setup
  public :: soca_initialize_integration
  public :: soca_finalize_integration    
  public :: soca_propagate
  public :: soca_delete
  
  !> Fortran derived type to hold configuration data for the  model
  type :: soca_model
     integer :: nx                !< Zonal grid dimension
     integer :: ny                !< Meridional grid dimension
     real(kind=kind_real) :: dt0  !< dimensional time (seconds)
     integer                :: advance_mom6 !< call mom6 step if true
     type(soca_mom6_config) :: mom6_config  !< MOM6 data structure
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

  ! ------------------------------------------------------------------------------
  !> Initialize model's data structure
  subroutine soca_setup(self)
    type(soca_model), intent(inout) :: self

    call soca_mom6_init(self%mom6_config)

  end subroutine soca_setup

  ! ------------------------------------------------------------------------------
  !> Prepare MOM6 integration
  subroutine soca_initialize_integration(self, flds)
    type(soca_model), intent(inout) :: self
    type(soca_field), intent(inout) :: flds
    
    integer :: isc, iec, jsc, jec
    type(time_type) :: ocean_time   ! The ocean model's clock.
    integer :: year, month, day, hour, minute, second
    character(len=20)  :: strdate
    character(len=1024)  :: buf

    ! Update halo
    call mpp_update_domains(flds%tocn, flds%geom%ocean%G%Domain%mpp_domain)
    call mpp_update_domains(flds%socn, flds%geom%ocean%G%Domain%mpp_domain)

    ! Update MOM's T and S to soca's
    self%mom6_config%MOM_CSp%T = real(flds%tocn, kind=8)
    self%mom6_config%MOM_CSp%S = real(flds%socn, kind=8)

  end subroutine soca_initialize_integration
  
  ! ------------------------------------------------------------------------------
  !> Advance MOM6 one baroclinic time step
  subroutine soca_propagate(self, flds, fldsdate)
    type(soca_model), intent(inout) :: self
    type(soca_field), intent(inout) :: flds
    type(datetime),       intent(in):: fldsdate
    
    integer :: isc, iec, jsc, jec
    type(time_type) :: ocean_time   ! The ocean model's clock.
    integer :: year, month, day, hour, minute, second
    character(len=20)  :: strdate
    character(len=1024)  :: buf

    ! Update halo
    call mpp_update_domains(flds%tocn, flds%geom%ocean%G%Domain%mpp_domain)
    call mpp_update_domains(flds%socn, flds%geom%ocean%G%Domain%mpp_domain)

    ! Update MOM's T and S to soca's
    self%mom6_config%MOM_CSp%T = real(flds%tocn, kind=8)
    self%mom6_config%MOM_CSp%S = real(flds%socn, kind=8)

    ! Set ocean clock
    call datetime_to_string(fldsdate, strdate)
    call soca_str2int(strdate(1:4), year)
    call soca_str2int(strdate(6:7), month)
    call soca_str2int(strdate(9:10), day)
    call soca_str2int(strdate(12:13), hour)
    call soca_str2int(strdate(15:16), minute)    
    call soca_str2int(strdate(18:19), second)    
    self%mom6_config%Time = set_date(year, month, day, hour, minute, second)
    ocean_time = self%mom6_config%Time

    WRITE(buf,*) 'Advancing MOM6 1 time step, starting from: '//&
         &strdate(1:4)//'-'//&
         &strdate(6:7)//'-'//&
         &strdate(9:10)//' '//&
         &strdate(12:13)//':'//&
         &strdate(15:16)
    call log%info(buf,newl=.true.)
    
    if (self%advance_mom6==1) then
       ! Advance MOM in a single step call (advance dyna and thermo)
       call step_MOM(self%mom6_config%forces, &
                     &self%mom6_config%fluxes, &
                     &self%mom6_config%sfc_state, &
                     &self%mom6_config%Time, &
                     &real(self%mom6_config%dt_forcing, kind=8), &
                     &self%mom6_config%MOM_CSp,&
                     &start_cycle=.false.,&
                     &cycle_length=self%mom6_config%MOM_CSp%dt)
    else
       !TODO: Read file
       WRITE(buf,*) 'IO Advance of MOM6: NOT IMPLEMENTED'
       call log%info(buf,newl=.true.)       
    end if
       
    ! Update ocean clock
    ocean_time = ocean_time + real_to_time(self%mom6_config%MOM_CSp%dt)
    self%mom6_config%Time = ocean_time

    ! Update soca fields
    flds%tocn = real(self%mom6_config%MOM_CSp%T, kind=kind_real)
    flds%socn = real(self%mom6_config%MOM_CSp%S, kind=kind_real)
    flds%hocn = real(self%mom6_config%MOM_CSp%h, kind=kind_real)
    flds%ssh = real(self%mom6_config%MOM_CSp%ave_ssh_ibc, kind=kind_real)
    
  end subroutine soca_propagate

  ! ------------------------------------------------------------------------------
  !> Finalize MOM6 integration: Update mom6's state and checkpoint
  subroutine soca_finalize_integration(self, flds)
    type(soca_model), intent(inout) :: self
    type(soca_field), intent(inout) :: flds
    
    integer :: isc, iec, jsc, jec
    type(time_type) :: ocean_time   ! The ocean model's clock.
    integer :: year, month, day, hour, minute, second
    character(len=20)  :: strdate
    character(len=1024)  :: buf

    ! Update halo
    call mpp_update_domains(flds%tocn, flds%geom%ocean%G%Domain%mpp_domain)
    call mpp_update_domains(flds%socn, flds%geom%ocean%G%Domain%mpp_domain)

    ! Update MOM's T and S to soca's
    self%mom6_config%MOM_CSp%T = real(flds%tocn, kind=8)
    self%mom6_config%MOM_CSp%S = real(flds%socn, kind=8)

    ! Checkpoint MOM
    call save_restart(self%mom6_config%dirs%restart_output_dir, &
                     &self%mom6_config%Time,&
                     &self%mom6_config%grid,&
                     &self%mom6_config%restart_CSp,&
                     &GV=self%mom6_config%GV)

  end subroutine soca_finalize_integration
  
  ! ------------------------------------------------------------------------------
  !> Release memory
  subroutine soca_delete(self)
    type(soca_model), intent(inout) :: self

    call soca_mom6_end(self%mom6_config)

  end subroutine soca_delete
  
end module soca_model_mod
