! (C) Copyright 2020-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for soca_increment_mod::soca_increment
module soca_increment_mod_c

use atlas_module, only: atlas_fieldset
use datetime_mod, only: datetime, c_f_datetime
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds, only: kind_real
use oops_variables_mod, only : oops_variables

! soca modules
use soca_geom_iter_mod_c, only: soca_geom_iter_registry
use soca_geom_iter_mod, only : soca_geom_iter
use soca_geom_mod_c, only: soca_geom_registry
use soca_geom_mod, only: soca_geom
use soca_increment_mod, only : soca_increment
use soca_increment_reg, only : soca_increment_registry
use soca_state_mod, only : soca_state
use soca_state_reg, only : soca_state_registry

implicit none
private

contains


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment version of
!! soca_fields_mod::soca_fields::create()
subroutine soca_increment_create_c(c_key_self, c_key_geom, c_vars) bind(c,name='soca_increment_create_f90')
  integer(c_int), intent(inout) :: c_key_self !< Handle to field
  integer(c_int),    intent(in) :: c_key_geom !< Geometry
  type(c_ptr),value, intent(in) :: c_vars     !< List of variables

  type(soca_increment),pointer :: self
  type(soca_geom),  pointer :: geom
  type(oops_variables)      :: vars

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_increment_registry%init()
  call soca_increment_registry%add(c_key_self)
  call soca_increment_registry%get(c_key_self,self)

  vars = oops_variables(c_vars)
  call self%create(geom, vars)

end subroutine soca_increment_create_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment version of
!! soca_fields_mod::soca_fields::delete()
subroutine soca_increment_delete_c(c_key_self) bind(c,name='soca_increment_delete_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(soca_increment),    pointer :: self

  call soca_increment_registry%get(c_key_self,self)
  call self%delete( )
  call soca_increment_registry%remove(c_key_self)

end subroutine soca_increment_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment::dirac()
subroutine soca_increment_dirac_c(c_key_self,c_conf) bind(c,name='soca_increment_dirac_f90')
  integer(c_int), intent(in) :: c_key_self
  type(c_ptr),    intent(in) :: c_conf !< Configuration

  type(soca_increment), pointer :: self

  call soca_increment_registry%get(c_key_self,self)
  call self%dirac(fckit_configuration(c_conf))

end subroutine soca_increment_dirac_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment::random()
subroutine soca_increment_random_c(c_key_self) bind(c,name='soca_increment_random_f90')
  integer(c_int), intent(in) :: c_key_self

  type(soca_increment), pointer :: self

  call soca_increment_registry%get(c_key_self,self)
  call self%random()

end subroutine soca_increment_random_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment version of
!! soca_fields_mod::soca_fields::copy()
subroutine soca_increment_copy_c(c_key_self,c_key_rhs) bind(c,name='soca_increment_copy_f90')
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(soca_increment), pointer :: self
  type(soca_increment), pointer :: rhs

  call soca_increment_registry%get(c_key_self,self)
  call soca_increment_registry%get(c_key_rhs,rhs)

  call self%copy(rhs)

end subroutine soca_increment_copy_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment::convert()
subroutine soca_increment_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='soca_increment_change_resol_f90')
  integer(c_int), intent(in) :: c_key_fld
  integer(c_int), intent(in) :: c_key_rhs

  type(soca_increment), pointer :: fld, rhs

  call soca_increment_registry%get(c_key_fld,fld)
  call soca_increment_registry%get(c_key_rhs,rhs)

  ! TODO (Guillaume or Travis) implement == in geometry or something to that effect.
  if ( size(fld%geom%lon,1)==size(rhs%geom%lon,1) .and. &
        size(fld%geom%lat,2)==size(rhs%geom%lat,2) .and. &
        fld%geom%nzo==rhs%geom%nzo ) then
      call fld%copy(rhs)
  else
      call fld%convert(rhs)
  endif

end subroutine soca_increment_change_resol_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment::to_atlas()
subroutine soca_increment_to_fieldset_c(c_key_self, c_vars, c_afieldset) &
     bind (c,name='soca_increment_to_fieldset_f90')
  integer(c_int),         intent(in) :: c_key_self
  type(c_ptr),     value, intent(in) :: c_vars
  type(c_ptr),     value, intent(in) :: c_afieldset

  type(soca_increment), pointer :: self
  type(oops_variables) :: vars
  type(atlas_fieldset) :: afieldset

  call soca_increment_registry%get(c_key_self,self)
  vars = oops_variables(c_vars)
  afieldset = atlas_fieldset(c_afieldset)

  call self%to_fieldset(vars, afieldset)

end subroutine


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment::from_fieldset()
subroutine soca_increment_from_fieldset_c(c_key_self, c_vars, c_afieldset) &
    bind (c,name='soca_increment_from_fieldset_f90')
  integer(c_int),         intent(in) :: c_key_self
  type(c_ptr),     value, intent(in) :: c_vars
  type(c_ptr),     value, intent(in) :: c_afieldset

  type(soca_increment), pointer :: self
  type(oops_variables)          :: vars
  type(atlas_fieldset)          :: afieldset

  call soca_increment_registry%get(c_key_self, self)
  vars = oops_variables(c_vars)
  afieldset = atlas_fieldset(c_afieldset)

  call self%from_fieldset(vars, afieldset)

end subroutine


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment version of
!! soca_fields_mod::soca_fields::read()
subroutine soca_increment_read_file_c(c_key_fld, c_conf, c_dt) bind(c,name='soca_increment_read_file_f90')
  integer(c_int), intent(in) :: c_key_fld  !< Fields
  type(c_ptr),    intent(in) :: c_conf     !< Configuration
  type(c_ptr), intent(inout) :: c_dt       !< DateTime

  type(soca_increment), pointer :: fld
  type(datetime)            :: fdate

  call soca_increment_registry%get(c_key_fld,fld)
  call c_f_datetime(c_dt, fdate)
  call fld%read(fckit_configuration(c_conf), fdate)

end subroutine soca_increment_read_file_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment version of
!! soca_fields_mod::soca_fields::write()
subroutine soca_increment_write_file_c(c_key_fld, c_conf, c_dt) bind(c,name='soca_increment_write_file_f90')
  integer(c_int), intent(in) :: c_key_fld  !< Fields
  type(c_ptr),    intent(in) :: c_conf     !< Configuration
  type(c_ptr),    intent(in) :: c_dt       !< DateTime

  type(soca_increment), pointer :: fld
  type(datetime)            :: fdate

  call soca_increment_registry%get(c_key_fld,fld)
  call c_f_datetime(c_dt, fdate)
  call fld%write_rst(fckit_configuration(c_conf), fdate)

end subroutine soca_increment_write_file_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment RMS
subroutine soca_increment_rms_c(c_key_fld, prms) bind(c,name='soca_increment_rms_f90')
  integer(c_int),    intent(in) :: c_key_fld
  real(c_double), intent(inout) :: prms

  type(soca_increment), pointer :: fld
  real(kind=kind_real)      :: zz

  call soca_increment_registry%get(c_key_fld,fld)

  call fld%dot_prod(fld, zz)
  prms = sqrt(zz)

end subroutine soca_increment_rms_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment::getpoint()
subroutine soca_increment_getpoint_c(c_key_fld,c_key_iter,values, values_len) bind(c,name='soca_increment_getpoint_f90')
  integer(c_int), intent(in) :: c_key_fld
  integer(c_int), intent(in) :: c_key_iter
  integer(c_int), intent(in) :: values_len
  real(c_double), intent(inout) :: values(values_len)

  type(soca_increment),      pointer :: fld
  type(soca_geom_iter), pointer :: iter

  call soca_increment_registry%get(c_key_fld,fld)
  call soca_geom_iter_registry%get(c_key_iter,iter)

  call fld%getpoint(iter, values)

end subroutine soca_increment_getpoint_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment::setpoint()
subroutine soca_increment_setpoint_c(c_key_fld,c_key_iter,values, values_len) bind(c,name='soca_increment_setpoint_f90')
  integer(c_int), intent(inout) :: c_key_fld
  integer(c_int), intent(in) :: c_key_iter
  integer(c_int), intent(in) :: values_len
  real(c_double), intent(in) :: values(values_len)

  type(soca_increment),      pointer :: fld
  type(soca_geom_iter), pointer :: iter

  call soca_increment_registry%get(c_key_fld,fld)
  call soca_geom_iter_registry%get(c_key_iter,iter)

  call fld%setpoint(iter, values)

end subroutine soca_increment_setpoint_c


! ------------------------------------------------------------------------------
!> C++ interface to get soca_increment_mod::soca_increment dimension sizes
subroutine soca_incrementnum_c(c_key_fld, nx, ny, nzo, nf) bind(c,name='soca_increment_sizes_f90')
  integer(c_int),         intent(in) :: c_key_fld
  integer(kind=c_int), intent(inout) :: nx, ny, nzo, nf

  type(soca_increment), pointer :: fld

  call soca_increment_registry%get(c_key_fld,fld)

  nx = size(fld%geom%lon,1)
  ny = size(fld%geom%lon,2)
  nzo = fld%geom%nzo
  nf = size(fld%fields)

end subroutine soca_incrementnum_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment version of
!! soca_fields_mod::soca_fields::update_fields()
subroutine soca_increment_update_fields_c(c_key_self, c_vars) &
           bind (c,name='soca_increment_update_fields_f90')

integer(c_int),     intent(in) :: c_key_self
type(c_ptr), value, intent(in) :: c_vars

type(soca_increment), pointer :: f_self
type(oops_variables)          :: f_vars

! LinkedList
! ----------
call soca_increment_registry%get(c_key_self, f_self)

! Fortrain APIs
! -------------
f_vars = oops_variables(c_vars)

! Call implementation
! -------------------
call f_self%update_fields(f_vars)

end subroutine soca_increment_update_fields_c



! ------------------------------------------------------------------------------
!> C++ interface for computing horizontal decorrelation length scales
!> soca_increment_mod::soca_increment
!! using  soca_increment_mod::soca_increment::horiz_scales()
subroutine soca_increment_horiz_scales_c(c_key_self, c_conf) &
     bind(c,name='soca_increment_horiz_scales_f90')
  integer(c_int), intent(in) :: c_key_self
  type(c_ptr),    intent(in) :: c_conf !< Configuration

  type(soca_increment), pointer :: f_self


  call soca_increment_registry%get(c_key_self, f_self)
  call f_self%horiz_scales(fckit_configuration(c_conf))

end subroutine soca_increment_horiz_scales_c



! ------------------------------------------------------------------------------
!> C++ interface for computing vertical decorrelation length scales
!> soca_increment_mod::soca_increment
!! using  soca_increment_mod::soca_increment::vert_scales()
subroutine soca_increment_vert_scales_c(c_key_self, c_vert) bind(c,name='soca_increment_vert_scales_f90')
  integer(c_int), intent(in) :: c_key_self
  real(c_double), intent(in) :: c_vert

  type(soca_increment), pointer :: f_self
  real(kind=kind_real) :: vert

  call soca_increment_registry%get(c_key_self, f_self)
  vert = c_vert
  call f_self%vert_scales(c_vert)

end subroutine soca_increment_vert_scales_c

! ------------------------------------------------------------------------------
end module
