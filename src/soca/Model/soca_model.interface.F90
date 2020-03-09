! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Setup the model

module soca_model_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod, only: datetime, c_f_datetime
use duration_mod, only: duration, duration_seconds, assignment(=)
use soca_geom_mod, only: soca_geom
use soca_fields_mod, only: soca_fields
use soca_fields_mod_c, only: soca_field_registry
use soca_model_mod, only: soca_model, soca_setup, soca_delete, soca_propagate, &
                          soca_initialize_integration, soca_finalize_integration

implicit none

private
public :: soca_model_registry

#define LISTED_TYPE soca_model

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: soca_model_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

subroutine c_soca_setup(c_conf, c_geom, c_key_model) bind (c,name='soca_setup_f90')

  type(c_ptr),         intent(in) :: c_conf       !< pointer to object of class Config
  type(c_ptr), target, intent(in) :: c_geom   !< Geometry
  integer(c_int),   intent(inout) :: c_key_model  !< Key to configuration data

  type(soca_model), pointer :: model
  type(soca_geom),  pointer :: geom

  type(duration) :: dtstep
  real(c_double), allocatable :: tocn_minmax(:), socn_minmax(:)
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str

  f_conf = fckit_configuration(c_conf)

  call c_f_pointer(c_geom, geom)

  call soca_model_registry%init()
  call soca_model_registry%add(c_key_model)
  call soca_model_registry%get(c_key_model, model)

  ! Get local grid size
  model%nx =  size(geom%lon,1)
  model%ny = size(geom%lon,2)

  ! Setup time step
  call f_conf%get_or_die("tstep", str)
  dtstep = trim(str)
  model%dt0 = duration_seconds(dtstep)

  ! Setup mom6 advance or identity model
  call f_conf%get_or_die("advance_mom6", model%advance_mom6)

  ! Setup defaults for clamping values in the model
  if ( f_conf%has("tocn_minmax") ) then
    call f_conf%get_or_die("tocn_minmax", tocn_minmax)
    model%tocn_minmax = tocn_minmax
  else
    model%tocn_minmax=(/-999., -999./)
  endif
  if ( f_conf%has("socn_minmax") ) then
    call f_conf%get_or_die("socn_minmax", socn_minmax)
    model%socn_minmax = socn_minmax
  else
    model%socn_minmax=(/-999., -999./)
  endif

  ! Initialize mom6
  call soca_setup(model)

  if (allocated(str)) deallocate(str)

  return
end subroutine c_soca_setup

! ------------------------------------------------------------------------------

!> Delete the model
subroutine c_soca_delete(c_key_conf) bind (c,name='soca_delete_f90')

  integer(c_int), intent(inout) :: c_key_conf  !< Key to configuration structure
  type(soca_model), pointer :: model

  call soca_model_registry%get(c_key_conf, model)
  call soca_delete(model)
  call soca_model_registry%remove(c_key_conf)

  return
end subroutine c_soca_delete

! ------------------------------------------------------------------------------
!> Prepare the model or integration
subroutine c_soca_initialize_integration(c_key_model, c_key_state) &
     & bind(c,name='soca_initialize_integration_f90')

  integer(c_int), intent(in) :: c_key_model  !< Configuration structure
  integer(c_int), intent(in) :: c_key_state  !< Model fields

  type(soca_model), pointer :: model
  type(soca_fields),pointer :: flds

  call soca_field_registry%get(c_key_state, flds)
  call soca_model_registry%get(c_key_model, model)

  call soca_initialize_integration(model, flds)

  return
end subroutine c_soca_initialize_integration

! ------------------------------------------------------------------------------

!> Checkpoint model
subroutine c_soca_finalize_integration(c_key_model, c_key_state) &
           bind(c,name='soca_finalize_integration_f90')

  integer(c_int), intent(in) :: c_key_model  !< Configuration structure
  integer(c_int), intent(in) :: c_key_state  !< Model fields

  type(soca_model), pointer :: model
  type(soca_fields),pointer :: flds

  call soca_field_registry%get(c_key_state, flds)
  call soca_model_registry%get(c_key_model, model)

  call soca_finalize_integration(model, flds)

  return
end subroutine c_soca_finalize_integration

! ------------------------------------------------------------------------------

!> Perform a timestep of the model
subroutine c_soca_propagate(c_key_model, c_key_state, c_key_date) bind(c,name='soca_propagate_f90')

  integer(c_int), intent(in) :: c_key_model  !< Config structure
  integer(c_int), intent(in) :: c_key_state  !< Model fields
  type(c_ptr), intent(inout) :: c_key_date   !< DateTime

  type(soca_model), pointer :: model
  type(soca_fields),pointer :: flds
  type(datetime)            :: fldsdate

  call soca_model_registry%get(c_key_model, model)
  call soca_field_registry%get(c_key_state, flds)
  call c_f_datetime(c_key_date, fldsdate)

  call soca_propagate(model, flds, fldsdate)

  return
end subroutine c_soca_propagate

end module soca_model_mod_c
