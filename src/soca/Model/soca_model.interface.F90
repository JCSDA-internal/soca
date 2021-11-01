! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for soca_model_mod::soca_model
module soca_model_mod_c

use datetime_mod, only: datetime, c_f_datetime
use duration_mod, only: duration, duration_seconds, assignment(=)
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding

! soca modules
use soca_geom_mod_c, only: soca_geom_registry
use soca_geom_mod, only: soca_geom
use soca_model_mod, only: soca_model
use soca_state_mod, only: soca_state
use soca_state_reg, only: soca_state_registry

implicit none
private

#define LISTED_TYPE soca_model

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for soca_model instances
type(registry_t), public :: soca_model_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
!> C++ interface for soca_model_mod::soca_model::setup()
subroutine soca_model_setup_c(c_conf, c_key_geom, c_key_model) bind (c,name='soca_model_setup_f90')
  type(c_ptr),       intent(in) :: c_conf       !< pointer to object of class Config
  integer(c_int),    intent(in) :: c_key_geom   !< Geometry
  integer(c_int), intent(inout) :: c_key_model  !< Key to configuration data

  type(soca_model), pointer :: model
  type(soca_geom),  pointer :: geom

  type(duration) :: dtstep
  real(c_double), allocatable :: tocn_minmax(:), socn_minmax(:)
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str

  f_conf = fckit_configuration(c_conf)

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_model_registry%init()
  call soca_model_registry%add(c_key_model)
  call soca_model_registry%get(c_key_model, model)

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
  call model%setup(geom)

  if (allocated(str)) deallocate(str)
end subroutine soca_model_setup_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_model_mod::soca_model::delete()
subroutine soca_model_delete_c(c_key_conf) bind (c,name='soca_model_delete_f90')
  integer(c_int), intent(inout) :: c_key_conf  !< Key to configuration structure

  type(soca_model), pointer :: model

  call soca_model_registry%get(c_key_conf, model)
  call model%delete()
  call soca_model_registry%remove(c_key_conf)
end subroutine soca_model_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_model_mod::soca_model::init()
subroutine soca_model_init_c(c_key_model, c_key_state) &
     & bind(c,name='soca_model_init_f90')
  integer(c_int), intent(in) :: c_key_model  !< Configuration structure
  integer(c_int), intent(in) :: c_key_state  !< Model fields

  type(soca_model), pointer :: model
  type(soca_state),pointer :: flds

  call soca_state_registry%get(c_key_state, flds)
  call soca_model_registry%get(c_key_model, model)

  call model%init(flds)

end subroutine soca_model_init_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_model_mod::soca_model::finalize()
subroutine soca_model_finalize_c(c_key_model, c_key_state) &
           bind(c,name='soca_model_finalize_f90')
  integer(c_int), intent(in) :: c_key_model  !< Configuration structure
  integer(c_int), intent(in) :: c_key_state  !< Model fields

  type(soca_model), pointer :: model
  type(soca_state),pointer :: flds

  call soca_state_registry%get(c_key_state, flds)
  call soca_model_registry%get(c_key_model, model)

  call model%finalize(flds)
end subroutine soca_model_finalize_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_model_mod::soca_model::propagate()
subroutine soca_model_propagate_c(c_key_model, c_key_state, c_key_date) bind(c,name='soca_model_propagate_f90')
  integer(c_int), intent(in) :: c_key_model  !< Config structure
  integer(c_int), intent(in) :: c_key_state  !< Model fields
  type(c_ptr), intent(inout) :: c_key_date   !< DateTime

  type(soca_model), pointer :: model
  type(soca_state),pointer :: flds
  type(datetime)            :: fldsdate

  call soca_model_registry%get(c_key_model, model)
  call soca_state_registry%get(c_key_state, flds)
  call c_f_datetime(c_key_date, fldsdate)

  call model%propagate(flds, fldsdate)
end subroutine soca_model_propagate_c

end module soca_model_mod_c
