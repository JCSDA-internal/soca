! (C) Copyright 2020-2024 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for soca_state_mod::soca_state
module soca_state_mod_c

use atlas_module, only: atlas_fieldset
use datetime_mod, only: datetime, c_f_datetime
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables

! soca modules
use soca_fields_mod, only: soca_field
use soca_geom_mod_c, only: soca_geom_registry
use soca_geom_mod, only: soca_geom
use soca_increment_mod, only: soca_increment
use soca_increment_reg, only: soca_increment_registry
use soca_state_mod, only: soca_state
use soca_state_reg, only: soca_state_registry
use soca_analytic_mod, only: soca_analytic_state

implicit none
private


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state version of
!! soca_fields_mod::soca_fields::create()
subroutine soca_state_create_c(c_key_self, c_key_geom, c_vars, c_afieldsest) &
  bind(c,name='soca_state_create_f90')
    integer(c_int), intent(inout) :: c_key_self !< Handle to field
    integer(c_int),    intent(in) :: c_key_geom !< Geometry
    type(c_ptr),value, intent(in) :: c_vars     !< List of variables
    type(c_ptr),value, intent(in) :: c_afieldsest

    type(soca_state),pointer :: self
    type(soca_geom),  pointer :: geom
    type(oops_variables)      :: vars
    type(atlas_fieldset)      :: afieldset

    call soca_geom_registry%get(c_key_geom, geom)
    call soca_state_registry%init()
    call soca_state_registry%add(c_key_self)
    call soca_state_registry%get(c_key_self,self)

    vars = oops_variables(c_vars)
    afieldset = atlas_fieldset(c_afieldsest)
    call self%create(geom, vars, afieldset)
    call afieldset%final()

end subroutine soca_state_create_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state version of
!! soca_fields_mod::soca_fields::delete()
subroutine soca_state_delete_c(c_key_self) bind(c,name='soca_state_delete_f90')
    integer(c_int), intent(inout) :: c_key_self

    type(soca_state),    pointer :: self

    call soca_state_registry%get(c_key_self,self)
    call self%delete( )
    call soca_state_registry%remove(c_key_self)

end subroutine soca_state_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state version of
!! soca_fields_mod::soca_fields::read()
subroutine soca_state_read_file_c(c_key_fld, c_conf, c_dt) bind(c,name='soca_state_read_file_f90')
    integer(c_int), intent(in) :: c_key_fld  !< Fields
    type(c_ptr),    intent(in) :: c_conf     !< Configuration
    type(c_ptr), intent(inout) :: c_dt       !< DateTime

    type(soca_state), pointer :: fld
    type(datetime)            :: fdate

    call soca_state_registry%get(c_key_fld,fld)
    call c_f_datetime(c_dt, fdate)
    call fld%read(fckit_configuration(c_conf), fdate)

end subroutine soca_state_read_file_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state version of
!! soca_fields_mod::soca_fields::write_rst()
subroutine soca_state_write_file_c(c_key_fld, c_conf, c_dt) bind(c,name='soca_state_write_file_f90')
    integer(c_int), intent(in) :: c_key_fld  !< Fields
    type(c_ptr),    intent(in) :: c_conf     !< Configuration
    type(c_ptr),    intent(in) :: c_dt       !< DateTime

    type(soca_state), pointer :: fld
    type(datetime)            :: fdate

    call soca_state_registry%get(c_key_fld,fld)
    call c_f_datetime(c_dt, fdate)
    call fld%write_rst(fckit_configuration(c_conf), fdate)

end subroutine soca_state_write_file_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::tohgrid()
subroutine soca_state_tohgrid_c(c_key_self) bind(c,name='soca_state_tohgrid_f90')
  integer(c_int),     intent(in) :: c_key_self

  type(soca_state), pointer :: self

  call soca_state_registry%get(c_key_self,self)
  call self%tohpoints()

end subroutine soca_state_tohgrid_c


! ------------------------------------------------------------------------------
subroutine scoa_state_analytic_c(c_key_self, c_conf, c_dt) &
    bind(c,name='soca_state_analytic_f90')
  integer (c_int),     intent(in   ) :: c_key_self
  TYPE (c_ptr), value, intent(in   ) :: c_conf
  TYPE (c_ptr),        intent(inout) :: c_dt

  type(soca_state), pointer :: self
  type(datetime) :: fdate

  call soca_state_registry%get(c_key_self,self)
  call c_f_datetime (c_dt, fdate)
  call soca_analytic_state(self)

end subroutine scoa_state_analytic_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state version of
!! soca_fields_mod::soca_fields::update_fields()
subroutine soca_state_update_fields_c(c_key_self, c_vars) &
           bind (c,name='soca_state_update_fields_f90')

integer(c_int),     intent(in) :: c_key_self
type(c_ptr), value, intent(in) :: c_vars

type(soca_state), pointer :: f_self
type(oops_variables)      :: f_vars

! LinkedList
! ----------
call soca_state_registry%get(c_key_self, f_self)

! Fortrain APIs
! -------------
f_vars = oops_variables(c_vars)

! Call implementation
! -------------------
call f_self%update_fields(f_vars)

end subroutine soca_state_update_fields_c

end module soca_state_mod_c
