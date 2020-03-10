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
use soca_model_mod, only: soca_model, soca_setup, soca_delete, soca_propagate, &
                          soca_initialize_integration, soca_finalize_integration

implicit none
private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine c_soca_setup(c_self, c_conf, c_geom) bind (c,name='soca_setup_f90')
  type(c_ptr),      intent(inout) :: c_self  !< Key to configuration data
  type(c_ptr),         intent(in) :: c_conf  !< pointer to object of class Config
  type(c_ptr), target, intent(in) :: c_geom  !< Geometry

  type(soca_model), pointer :: self
  type(soca_geom),  pointer :: geom

  type(duration) :: dtstep
  real(c_double), allocatable :: tocn_minmax(:), socn_minmax(:)
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str

  f_conf = fckit_configuration(c_conf)

  allocate(self)
  c_self = c_loc(self)
  call c_f_pointer(c_geom, geom)

  ! Get local grid size
  self%nx =  size(geom%lon,1)
  self%ny = size(geom%lon,2)

  ! Setup time step
  call f_conf%get_or_die("tstep", str)
  dtstep = trim(str)
  self%dt0 = duration_seconds(dtstep)

  ! Setup mom6 advance or identity model
  call f_conf%get_or_die("advance_mom6", self%advance_mom6)

  ! Setup defaults for clamping values in the model
  if ( f_conf%has("tocn_minmax") ) then
    call f_conf%get_or_die("tocn_minmax", tocn_minmax)
    self%tocn_minmax = tocn_minmax
  else
    self%tocn_minmax=(/-999., -999./)
  endif
  if ( f_conf%has("socn_minmax") ) then
    call f_conf%get_or_die("socn_minmax", socn_minmax)
    self%socn_minmax = socn_minmax
  else
    self%socn_minmax=(/-999., -999./)
  endif

  ! Initialize mom6
  call soca_setup(self)

  if (allocated(str)) deallocate(str)

  return
end subroutine c_soca_setup

! ------------------------------------------------------------------------------

!> Delete the model
subroutine c_soca_delete(c_self) bind (c,name='soca_delete_f90')
  type(c_ptr), intent(inout) :: c_self

  type(soca_model), pointer :: self

  call c_f_pointer(c_self, self)

  call soca_delete(self)
  deallocate(self)

end subroutine c_soca_delete

! ------------------------------------------------------------------------------
!> Prepare the model or integration
subroutine c_soca_initialize_integration(c_self, c_state) &
  & bind(c,name='soca_initialize_integration_f90')
  type(c_ptr), intent(inout) :: c_self       !< Configuration structure
  type(c_ptr), intent(inout) :: c_state  !< Model fields

  type(soca_model), pointer :: self
  type(soca_fields),pointer :: flds

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_state, flds)

  call soca_initialize_integration(self, flds)

end subroutine c_soca_initialize_integration

! ------------------------------------------------------------------------------

!> Checkpoint model
subroutine c_soca_finalize_integration(c_self, c_state) &
  & bind(c,name='soca_finalize_integration_f90')
  type(c_ptr), intent(inout) :: c_self       !< Configuration structure
  type(c_ptr), intent(inout) :: c_state  !< Model fields

  type(soca_model), pointer :: self
  type(soca_fields),pointer :: flds

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_state, flds)

  call soca_finalize_integration(self, flds)

end subroutine c_soca_finalize_integration

! ------------------------------------------------------------------------------

!> Perform a timestep of the model
subroutine c_soca_propagate(c_self, c_state, c_date) &
  & bind(c,name='soca_propagate_f90')
  type(c_ptr), intent(inout) :: c_self  !< Config structure
  type(c_ptr), intent(inout) :: c_state  !< Model fields
  type(c_ptr), intent(inout) :: c_date   !< DateTime

  type(soca_model), pointer :: self
  type(soca_fields),pointer :: flds
  type(datetime)            :: fldsdate

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_state, flds)
  call c_f_datetime(c_date, fldsdate)

  call soca_propagate(self, flds, fldsdate)

end subroutine c_soca_propagate

end module soca_model_mod_c
