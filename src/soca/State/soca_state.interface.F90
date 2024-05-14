! (C) Copyright 2020-2024 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for soca_state_mod::soca_state
module soca_state_mod_c

use atlas_module, only: atlas_fieldset, atlas_field, atlas_real, atlas_metadata
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
    call self%sync_to_atlas()

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
!! soca_fields_mod::soca_fields::copy()
subroutine soca_state_copy_c(c_key_self,c_key_rhs) bind(c,name='soca_state_copy_f90')
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_rhs

    type(soca_state), pointer :: self
    type(soca_state), pointer :: rhs

    call soca_state_registry%get(c_key_self,self)
    call soca_state_registry%get(c_key_rhs,rhs)

    call rhs%sync_from_atlas()
    call self%copy(rhs)
    call self%sync_to_atlas()

end subroutine soca_state_copy_c


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
    call fld%sync_to_atlas()

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
    call fld%sync_from_atlas()
    call c_f_datetime(c_dt, fdate)
    call fld%write_rst(fckit_configuration(c_conf), fdate)

end subroutine soca_state_write_file_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::rotate()
subroutine soca_state_rotate2grid_c(c_key_self, c_uvars, c_vvars) bind(c,name='soca_state_rotate2grid_f90')
  integer(c_int), intent(in)     :: c_key_self
  type(c_ptr), value, intent(in) :: c_uvars
  type(c_ptr), value, intent(in) :: c_vvars

  type(soca_state), pointer :: self
  type(oops_variables) :: uvars, vvars

  uvars = oops_variables(c_uvars)
  vvars = oops_variables(c_vvars)

  call soca_state_registry%get(c_key_self,self)
  call self%rotate(coordinate="grid", uvars=uvars, vvars=vvars)
  call self%sync_to_atlas()

end subroutine soca_state_rotate2grid_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::rotate()
subroutine soca_state_rotate2north_c(c_key_self, c_uvars, c_vvars) bind(c,name='soca_state_rotate2north_f90')
  integer(c_int),     intent(in) :: c_key_self
  type(c_ptr), value, intent(in) :: c_uvars
  type(c_ptr), value, intent(in) :: c_vvars

  type(soca_state), pointer :: self
  type(oops_variables) :: uvars, vvars

  uvars = oops_variables(c_uvars)
  vvars = oops_variables(c_vvars)

  call soca_state_registry%get(c_key_self,self)
  call self%rotate(coordinate="north", uvars=uvars, vvars=vvars)
  call self%sync_to_atlas()

end subroutine soca_state_rotate2north_c

! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::tohgrid()
subroutine soca_state_tohgrid_c(c_key_self) bind(c,name='soca_state_tohgrid_f90')
  integer(c_int),     intent(in) :: c_key_self

  type(soca_state), pointer :: self

  call soca_state_registry%get(c_key_self,self)
  call self%tohpoints()
  call self%sync_to_atlas()

end subroutine soca_state_tohgrid_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::convert()
subroutine soca_state_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='soca_state_change_resol_f90')
    integer(c_int), intent(in) :: c_key_fld
    integer(c_int), intent(in) :: c_key_rhs

    type(soca_state), pointer :: fld, rhs

    call soca_state_registry%get(c_key_fld,fld)
    call soca_state_registry%get(c_key_rhs,rhs)

    call rhs%sync_from_atlas()

    ! TODO (Guillaume or Travis) implement == in geometry or something to that effect.
    if (size(fld%geom%lon,1)==size(rhs%geom%lon,1) .and. size(fld%geom%lat,2)==size(rhs%geom%lat,2) .and. &
      fld%geom%nzo==rhs%geom%nzo ) then
      call fld%copy(rhs)
    else
      call fld%convert(rhs)
    endif
    call fld%sync_to_atlas()

end subroutine soca_state_change_resol_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::logexpon()
subroutine soca_state_logtrans_c(c_key_self, c_trvars) bind(c,name='soca_state_logtrans_f90')
  integer(c_int), intent(in)     :: c_key_self
  type(c_ptr), value, intent(in) :: c_trvars

  type(soca_state), pointer :: self
  type(oops_variables) :: trvars

  trvars = oops_variables(c_trvars)

  call soca_state_registry%get(c_key_self,self)
  call self%sync_from_atlas()
  call self%logexpon(transfunc="log", trvars=trvars)
  call self%sync_to_atlas()

end subroutine soca_state_logtrans_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::logexpon()
subroutine soca_state_expontrans_c(c_key_self, c_trvars) bind(c,name='soca_state_expontrans_f90')
  integer(c_int), intent(in)     :: c_key_self
  type(c_ptr), value, intent(in) :: c_trvars

  type(soca_state), pointer :: self
  type(oops_variables) :: trvars

  trvars = oops_variables(c_trvars)

  call soca_state_registry%get(c_key_self,self)
  call self%logexpon(transfunc="expon", trvars=trvars)
  call self%sync_to_atlas()

end subroutine soca_state_expontrans_c


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
  call self%sync_to_atlas()

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
call f_self%sync_from_atlas()

! Fortrain APIs
! -------------
f_vars = oops_variables(c_vars)

! Call implementation
! -------------------
call f_self%update_fields(f_vars)
call f_self%sync_to_atlas()

end subroutine soca_state_update_fields_c

end module soca_state_mod_c
