!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
!> Setup the model

subroutine c_soca_setup(c_confspec, c_key_geom, c_key_model) bind (c,name='soca_setup_f90')

  use soca_model_mod
  use soca_geom_mod_c
  use iso_c_binding
  use config_mod
  use duration_mod
  use kinds
  use fckit_log_module, only : fckit_log
  use mpi,             only : mpi_comm_world
  use mpp_mod,         only : mpp_init
  use mpp_domains_mod, only : mpp_domains_init
  use fms_io_mod,      only : fms_io_init

  implicit none

  type(c_ptr),       intent(in) :: c_confspec     !< pointer to object of class Config
  integer(c_int),    intent(in) :: c_key_geom     !< Geometry
  integer(c_int), intent(inout) :: c_key_model !< Key to configuration data

  type(soca_model), pointer :: model
  type(soca_geom),  pointer :: geom

  type(duration) :: dtstep
  character(len=20) :: ststep
  character(len=160) :: record
  integer :: int_logical

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_model_registry%init()
  call soca_model_registry%add(c_key_model)
  call soca_model_registry%get(c_key_model, model)

  ! Get local grid size
  model%nx  = geom%ocean%nx
  model%ny  = geom%ocean%ny

  ! Setup time step 
  ststep = config_get_string(c_confspec,len(ststep),"tstep")
  dtstep = trim(ststep)
  model%dt0 = duration_seconds(dtstep)

  ! Setup mom6 advance or identity model
  model%advance_mom6 = config_get_int(c_confspec,"advance_mom6")

  ! Initialize mom6
  if (model%advance_mom6==1) then
     call soca_create(model, geom, c_confspec)
  else
     print *,"Not initializing MOM6"
  end if

  return
end subroutine c_soca_setup

! ------------------------------------------------------------------------------

!> Delete the model

subroutine c_soca_delete(c_key_conf) bind (c,name='soca_delete_f90')

  use soca_model_mod
  use iso_c_binding

  implicit none
  integer(c_int), intent(inout) :: c_key_conf !< Key to configuration structure
  type(soca_model), pointer :: model

  call soca_model_registry%get(c_key_conf, model)
  if (model%advance_mom6==1) then
     call soca_delete(model)
  end if  
  call soca_model_registry%remove(c_key_conf)

  return
end subroutine c_soca_delete

! ------------------------------------------------------------------------------

subroutine c_soca_prepare_integration(c_key_model, c_key_state) &
     & bind(c,name='soca_prepare_integration_f90')

  use iso_c_binding
  use soca_fields
  use soca_model_mod
  use mpi,             only: mpi_comm_world
  use mpp_mod,         only: mpp_init
  use fms_mod,                 only: fms_init
  use fms_io_mod, only : fms_io_init
  use mpp_io_mod,              only: mpp_open, mpp_close
  implicit none
  integer(c_int), intent(in) :: c_key_model  !< Configuration structure
  integer(c_int), intent(in) :: c_key_state !< Model fields

  type(soca_model), pointer :: model
  type(soca_field), pointer  :: flds
  integer :: unit


  call soca_field_registry%get(c_key_state,flds)
  call soca_model_registry%get(c_key_model, model)

end subroutine c_soca_prepare_integration

  ! ------------------------------------------------------------------------------

!> Perform a timestep of the model
subroutine c_soca_propagate(c_key_model, c_key_state) bind(c,name='soca_propagate_f90')

  use iso_c_binding
  use soca_fields
  use soca_model_mod

  implicit none
  integer(c_int), intent(in) :: c_key_model  !< Config structure
  integer(c_int), intent(in) :: c_key_state !< Model fields

  type(soca_model), pointer :: model
  type(soca_field),  pointer :: flds

  call soca_model_registry%get(c_key_model, model)
  call soca_field_registry%get(c_key_state,flds)

  if (model%advance_mom6==1) then
     print *,"Advancing MOM6 1 time step"
     call soca_propagate(model, flds)
  else
     print *,"Not advancing MOM6: Persistence model"
  end if

  return
end subroutine c_soca_propagate
