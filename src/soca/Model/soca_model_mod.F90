! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Structure holding configuration variables for the model

module soca_model_mod

use datetime_mod, only: datetime, datetime_to_string
use kinds, only: kind_real

! mom6/fms modules
use fms_io_mod, only : fms_io_init, fms_io_exit
use MOM_restart, only : save_restart
use MOM_surface_forcing, only : set_forcing
use MOM_time_manager, only : operator(+)
use MOM_time_manager, only : real_to_time, time_type_to_real
use MOM, only : step_MOM
use mpp_domains_mod, only : mpp_update_domains
use time_manager_mod, only : time_type, set_date

! soca modules
use soca_fields_mod, only: soca_field
use soca_geom_mod, only: soca_geom
use soca_mom6, only: soca_mom6_config, soca_mom6_init, soca_mom6_end
use soca_state_mod, only: soca_state
use soca_utils, only: soca_str2int

implicit none
private


!> Fortran derived type to hold configuration data for the model
type, public :: soca_model
   integer :: advance_mom6      !< call mom6 step if true
   real(kind=kind_real) :: dt0  !< dimensional time (seconds)
   type(soca_mom6_config) :: mom6_config  !< MOM6 data structure
   real(kind_real), dimension(2) :: tocn_minmax, socn_minmax  !< min, max values

contains

  !> \name constructor/destructor
  !! \{

  !> \copybrief soca_model_setup \see soca_model_setup
  procedure :: setup => soca_model_setup

  !> \copybrief soca_model_delete \see soca_model_delete
  procedure :: delete => soca_model_delete

  !> \}

  !> \name model run steps
  !! \{

  !> \copybrief soca_model_init \see soca_model_init
  procedure :: init => soca_model_init

  !> \copybrief soca_model_propagate \see soca_model_propagate
  procedure :: propagate => soca_model_propagate

  !> \copybrief soca_model_finalize \see soca_model_finalize
  procedure :: finalize => soca_model_finalize

  !> \}

end type soca_model


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Initialize model's data structure
!!
!! \relates soca_model_mod::soca_model
subroutine soca_model_setup(self, geom)
  class(soca_model), intent(inout) :: self
  type(soca_geom),     intent(in) :: geom !< model geometry

  self%mom6_config%f_comm = geom%f_comm
  call soca_mom6_init(self%mom6_config)

end subroutine soca_model_setup


! ------------------------------------------------------------------------------
!> Prepare MOM6 integration
!!
!! \relates soca_model_mod::soca_model
subroutine soca_model_init(self, flds)
  class(soca_model), intent(inout) :: self
  type(soca_state), intent(inout) :: flds !< initial condition

  type(soca_field), pointer :: field
  integer :: i

  ! for each field
  do i=1,size(flds%fields)
    call flds%get(flds%fields(i)%name, field)

    ! Update halos
    call mpp_update_domains(field%val, flds%geom%Domain%mpp_domain)

    ! impose bounds, and set MOM6 state
    select case (field%name)
    case ("tocn")
      if ( self%tocn_minmax(1) /= real(-999., kind=8) ) &
        where( field%val < self%tocn_minmax(1) ) field%val = self%tocn_minmax(1)
      if ( self%tocn_minmax(2) /= real(-999., kind=8) ) &
        where( field%val > self%tocn_minmax(2) ) field%val = self%tocn_minmax(2)
      self%mom6_config%MOM_CSp%T = real(field%val, kind=8)
    case ("socn")
      if ( self%socn_minmax(1) /= real(-999., kind=8) ) &
        where( field%val < self%socn_minmax(1) ) field%val = self%socn_minmax(1)
      if ( self%socn_minmax(2) /= real(-999., kind=8) ) &
        where( field%val > self%socn_minmax(2) ) field%val = self%socn_minmax(2)
      self%mom6_config%MOM_CSp%S = real(field%val, kind=8)
    case ("uocn")
      self%mom6_config%MOM_CSp%u = real(field%val, kind=8)
    case ("vocn")
      self%mom6_config%MOM_CSp%v = real(field%val, kind=8)
    end select

    ! update forcing
    select case(field%name)
    case ("sw")
      field%val(:,:,1) = - real(self%mom6_config%fluxes%sw, kind=kind_real)
    case ("lw")
      field%val(:,:,1) = - real(self%mom6_config%fluxes%lw, kind=kind_real)
    case ("lhf")
      field%val(:,:,1) = - real(self%mom6_config%fluxes%latent, kind=kind_real)
    case ("shf")
      field%val(:,:,1) = - real(self%mom6_config%fluxes%sens, kind=kind_real)
    case ("us")
      field%val(:,:,1) =   real(self%mom6_config%fluxes%ustar, kind=kind_real)
    end select
  end do
end subroutine soca_model_init


! ------------------------------------------------------------------------------
!> Advance MOM6 one baroclinic time step
!!
!! \relates soca_model_mod::soca_model
subroutine soca_model_propagate(self, flds, fldsdate)
  class(soca_model), intent(inout) :: self
  type(soca_state), intent(inout) :: flds
  type(datetime),      intent(in) :: fldsdate

  type(soca_field), pointer :: field
  type(time_type) :: ocean_time  ! The ocean model's clock.
  integer :: year, month, day, hour, minute, second, i
  character(len=20) :: strdate

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

  if (self%advance_mom6==1) then
     ! Set the forcing for the next steps.
     call fms_io_init()
     call set_forcing(self%mom6_config%sfc_state,&
                      self%mom6_config%forces,&
                      self%mom6_config%fluxes,&
                      self%mom6_config%Time,&
                      self%mom6_config%Time_step_ocean,&
                      self%mom6_config%grid, &
                      self%mom6_config%scaling, &
                      self%mom6_config%surface_forcing_CSp)
     call fms_io_exit()

     ! Advance MOM in a single step call (advance dyna and thermo)
     call step_MOM(self%mom6_config%forces, &
                   self%mom6_config%fluxes, &
                   self%mom6_config%sfc_state, &
                   self%mom6_config%Time, &
                   real(self%mom6_config%dt_forcing, kind=8), &
                   self%mom6_config%MOM_CSp,&
                   start_cycle=.false.,&
                   cycle_length=self%mom6_config%MOM_CSp%dt)
  end if

  ! Update ocean clock
  ocean_time = ocean_time + real_to_time(self%mom6_config%MOM_CSp%dt)
  self%mom6_config%Time = ocean_time

  ! Update soca fields
  do i=1,size(flds%fields)
    field => flds%fields(i)
    select case(field%name)
    case ("tocn")
      field%val = real(self%mom6_config%MOM_CSp%T, kind=kind_real)
    case ("socn")
      field%val = real(self%mom6_config%MOM_CSp%S, kind=kind_real)
    case ("hocn")
      field%val = real(self%mom6_config%MOM_CSp%h, kind=kind_real)
    case ("ssh")
      field%val(:,:,1) = real(self%mom6_config%MOM_CSp%ave_ssh_ibc, kind=kind_real)
    case ("uocn")
      field%val = real(self%mom6_config%MOM_CSp%u, kind=kind_real)
    case ("vocn")
      field%val = real(self%mom6_config%MOM_CSp%v, kind=kind_real)
    case ("sw")
      field%val(:,:,1) = - real(self%mom6_config%fluxes%sw, kind=kind_real)
    case ("lw")
      field%val(:,:,1) = - real(self%mom6_config%fluxes%lw, kind=kind_real)
    case ("lhf")
      field%val(:,:,1) = - real(self%mom6_config%fluxes%latent, kind=kind_real)
    case ("shf")
      field%val(:,:,1) = - real(self%mom6_config%fluxes%sens, kind=kind_real)
    case ("us")
      field%val(:,:,1) = real(self%mom6_config%fluxes%ustar, kind=kind_real)
    end select
  end do
end subroutine soca_model_propagate


! ------------------------------------------------------------------------------
!> Finalize MOM6 integration: Update mom6's state and checkpoint
!!
!! \relates soca_model_mod::soca_model
subroutine soca_model_finalize(self, flds)
  class(soca_model), intent(inout) :: self
  type(soca_state), intent(inout) :: flds

  type(soca_field), pointer :: field
  real(kind=8), allocatable :: incr(:,:,:)
  real(kind=8), allocatable :: w(:,:,:), nz(:,:)
  integer :: i, k

  ! for each field
  do i=1,size(flds%fields)
    field => flds%fields(i)

    ! update halos
    call mpp_update_domains(field%val, flds%geom%Domain%mpp_domain)

    ! impose bounds and update MOM6
    select case(field%name)
    case ("tocn")
      if ( self%tocn_minmax(1) /= real(-999., kind=8) ) &
        where( field%val < self%tocn_minmax(1) ) field%val = self%tocn_minmax(1)
      if ( self%tocn_minmax(2) /= real(-999., kind=8) ) &
        where( field%val > self%tocn_minmax(2) ) field%val = self%tocn_minmax(2)
      self%mom6_config%MOM_CSp%T = real(field%val, kind=8)
    case ("socn")
      if ( self%socn_minmax(1) /= real(-999., kind=8) ) &
        where( field%val < self%socn_minmax(1) ) field%val = self%socn_minmax(1)
      if ( self%socn_minmax(2) /= real(-999., kind=8) ) &
        where( field%val > self%socn_minmax(2) ) field%val = self%socn_minmax(2)
      self%mom6_config%MOM_CSp%S = real(field%val, kind=8)
    case ("uocn")
      self%mom6_config%MOM_CSp%u = real(field%val, kind=8)
    case ("vocn")
      self%mom6_config%MOM_CSp%v = real(field%val, kind=8)
    case ("ssh")
      allocate(incr, mold=field%val)
      allocate(w, mold=self%mom6_config%MOM_CSp%h)
      incr = 0.0
      w = 0.0
      where ( self%mom6_config%MOM_CSp%h > 1.0e-3 )
         w = 1.0
      end where
      nz = sum(w, dim=3)
      where ( nz < 1.0 )
         nz = 1e9
      end where
      incr(:,:,1) = real(field%val(:,:,1), kind=8) - self%mom6_config%MOM_CSp%ave_ssh_ibc
      do k = 1, size(w, dim=3)
         self%mom6_config%MOM_CSp%h(:,:,k) = self%mom6_config%MOM_CSp%h(:,:,k) + &
              & field%mask(:,:)*w(:,:,k)*incr(:,:,1)/nz(:,:)
      end do

    end select
  end do

  ! Save MOM restarts with updated SOCA fields
  call save_restart(self%mom6_config%dirs%restart_output_dir, &
                   self%mom6_config%Time, &
                   self%mom6_config%grid, &
                   self%mom6_config%restart_CSp, &
                   GV=self%mom6_config%GV)

end subroutine soca_model_finalize


! ------------------------------------------------------------------------------
!> Release memory
!!
!! \relates soca_model_mod::soca_model
subroutine soca_model_delete(self)
  class(soca_model), intent(inout) :: self

  call soca_mom6_end(self%mom6_config)

end subroutine soca_model_delete

end module soca_model_mod
