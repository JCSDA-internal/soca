! (C) Copyright 2020-2021 UCAR
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
subroutine soca_state_create_c(c_key_self, c_key_geom, c_vars) bind(c,name='soca_state_create_f90')
    integer(c_int), intent(inout) :: c_key_self !< Handle to field
    integer(c_int),    intent(in) :: c_key_geom !< Geometry
    type(c_ptr),value, intent(in) :: c_vars     !< List of variables

    type(soca_state),pointer :: self
    type(soca_geom),  pointer :: geom
    type(oops_variables)      :: vars

    call soca_geom_registry%get(c_key_geom, geom)
    call soca_state_registry%init()
    call soca_state_registry%add(c_key_self)
    call soca_state_registry%get(c_key_self,self)

    vars = oops_variables(c_vars)
    call self%create(geom, vars)

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
!! soca_fields_mod::soca_fields::zeros()
subroutine soca_state_zero_c(c_key_self) bind(c,name='soca_state_zero_f90')
    integer(c_int), intent(in) :: c_key_self

    type(soca_state), pointer :: self

    call soca_state_registry%get(c_key_self,self)
    call self%zeros()

end subroutine soca_state_zero_c


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

    call self%copy(rhs)

end subroutine soca_state_copy_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state version of
!! soca_fields_mod::soca_fields::axpy()
subroutine soca_state_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='soca_state_axpy_f90')
    integer(c_int), intent(in) :: c_key_self
    real(c_double), intent(in) :: c_zz
    integer(c_int), intent(in) :: c_key_rhs

    type(soca_state), pointer :: self
    real(kind=kind_real)      :: zz
    type(soca_state), pointer :: rhs

    call soca_state_registry%get(c_key_self,self)
    call soca_state_registry%get(c_key_rhs,rhs)
    zz = c_zz

    call self%axpy(zz,rhs)

end subroutine soca_state_axpy_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::add_incr()
subroutine soca_state_add_incr_c(c_key_self,c_key_rhs) bind(c,name='soca_state_add_incr_f90')
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_rhs

    type(soca_state),     pointer :: self
    type(soca_increment), pointer :: rhs

    call soca_state_registry%get(c_key_self,self)
    call soca_increment_registry%get(c_key_rhs,rhs)

    call self%add_incr(rhs)

end subroutine soca_state_add_incr_c


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
!> C++ interface for soca_state_mod::soca_state version of
!! soca_fields_mod::soca_fields::gpnorm()
subroutine soca_state_gpnorm_c(c_key_fld, kf, pstat) bind(c,name='soca_state_gpnorm_f90')
    integer(c_int),    intent(in) :: c_key_fld
    integer(c_int),    intent(in) :: kf
    real(c_double), intent(inout) :: pstat(3*kf)

    type(soca_state), pointer :: fld
    real(kind=kind_real)      :: zstat(3, kf)
    integer :: jj, js, jf

    call soca_state_registry%get(c_key_fld,fld)

    call fld%gpnorm(kf, zstat)
    jj=0
    do jf = 1, kf
        do js = 1, 3
        jj=jj+1
        pstat(jj) = zstat(js,jf)
        enddo
    enddo

end subroutine soca_state_gpnorm_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state RMS
subroutine soca_state_rms_c(c_key_fld, prms) bind(c,name='soca_state_rms_f90')
    integer(c_int),    intent(in) :: c_key_fld
    real(c_double), intent(inout) :: prms

    type(soca_state), pointer :: fld
    real(kind=kind_real)      :: zz

    call soca_state_registry%get(c_key_fld,fld)

    call fld%dot_prod(fld, zz)
    prms = sqrt(zz)

end subroutine soca_state_rms_c


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

end subroutine soca_state_rotate2north_c

! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::tohgrid()
subroutine soca_state_tohgrid_c(c_key_self) bind(c,name='soca_state_tohgrid_f90')
  integer(c_int),     intent(in) :: c_key_self

  type(soca_state), pointer :: self

  call soca_state_registry%get(c_key_self,self)
  call self%tohpoints()

end subroutine soca_state_tohgrid_c

! ------------------------------------------------------------------------------
!> C++ interface to get soca_state_mod::soca_state dimensions sizes
subroutine soca_state_sizes_c(c_key_fld, nx, ny, nzo, nf) bind(c,name='soca_state_sizes_f90')
    integer(c_int),         intent(in) :: c_key_fld
    integer(kind=c_int), intent(inout) :: nx, ny, nzo, nf

    type(soca_state), pointer :: fld

    call soca_state_registry%get(c_key_fld,fld)

    nx = size(fld%geom%lon,1)
    ny = size(fld%geom%lon,2)
    nzo = fld%geom%nzo
    nf = size(fld%fields)

end subroutine soca_state_sizes_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::convert()
subroutine soca_state_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='soca_state_change_resol_f90')
    integer(c_int), intent(in) :: c_key_fld
    integer(c_int), intent(in) :: c_key_rhs

    type(soca_state), pointer :: fld, rhs

    call soca_state_registry%get(c_key_fld,fld)
    call soca_state_registry%get(c_key_rhs,rhs)

    ! TODO (Guillaume or Travis) implement == in geometry or something to that effect.
    if (size(fld%geom%lon,1)==size(rhs%geom%lon,1) .and. size(fld%geom%lat,2)==size(rhs%geom%lat,2) .and. &
      fld%geom%nzo==rhs%geom%nzo ) then
      call fld%copy(rhs)
    else
      call fld%convert(rhs)
    endif

end subroutine soca_state_change_resol_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state version of
!! soca_fields_mod::soca_fields::serial_size()
subroutine soca_state_serial_size_c(c_key_self,c_key_geom,c_vec_size) bind (c,name='soca_state_serial_size_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  integer(c_size_t), intent(out) :: c_vec_size

  type(soca_state), pointer :: self
  type(soca_geom),  pointer :: geom
  integer :: vec_size

  call soca_state_registry%get(c_key_self,self)
  call soca_geom_registry%get(c_key_geom,geom)

  call self%serial_size(geom, vec_size)
  c_vec_size = vec_size

end subroutine soca_state_serial_size_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state version of
!! soca_fields_mod::soca_fields::serialize()
subroutine soca_state_serialize_c(c_key_self,c_key_geom,c_vec_size,c_vec) bind (c,name='soca_state_serialize_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  integer(c_size_t), intent(in) :: c_vec_size
  real(c_double), intent(out) :: c_vec(c_vec_size)

  type(soca_state), pointer :: self
  type(soca_geom),  pointer :: geom

  integer :: vec_size

  vec_size = c_vec_size
  call soca_state_registry%get(c_key_self,self)
  call soca_geom_registry%get(c_key_geom,geom)

  call self%serialize(geom, vec_size, c_vec)

end subroutine soca_state_serialize_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state version of
!! soca_fields_mod::soca_fields::deserialize()
subroutine soca_state_deserialize_c(c_key_self,c_key_geom,c_vec_size,c_vec,c_index) bind (c,name='soca_state_deserialize_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  integer(c_size_t), intent(in) :: c_vec_size
  real(c_double), intent(in) :: c_vec(c_vec_size)
  integer(c_size_t), intent(inout) :: c_index

  type(soca_state), pointer :: self
  type(soca_geom),  pointer :: geom
  integer :: vec_size, idx

  vec_size = c_vec_size
  idx = c_index
  call soca_state_registry%get(c_key_self,self)
  call soca_geom_registry%get(c_key_geom,geom)

  call self%deserialize(geom,vec_size,c_vec, idx)
  c_index=idx

end subroutine soca_state_deserialize_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_state_mod::soca_state::logexpon()
subroutine soca_state_logtrans_c(c_key_self, c_trvars) bind(c,name='soca_state_logtrans_f90')
  integer(c_int), intent(in)     :: c_key_self
  type(c_ptr), value, intent(in) :: c_trvars

  type(soca_state), pointer :: self
  type(oops_variables) :: trvars

  trvars = oops_variables(c_trvars)

  call soca_state_registry%get(c_key_self,self)
  call self%logexpon(transfunc="log", trvars=trvars)

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


! ------------------------------------------------------------------------------
!> C++ interface for State version of soca_field_mod::soca_field::to_fieldset()
subroutine soca_state_to_fieldset_c(c_key_self, c_vars, c_fieldset) &
    bind (c, name='soca_state_to_fieldset_f90')
  integer(c_int),       intent(in) :: c_key_self
  type(c_ptr), value,   intent(in) :: c_vars
  type(c_ptr), value,   intent(in) :: c_fieldset

  type(soca_state), pointer :: self
  type(oops_variables)      :: vars
  type(atlas_fieldset)      :: afieldset

  call soca_state_registry%get(c_key_self, self)
  vars = oops_variables(c_vars)
  afieldset = atlas_fieldset(c_fieldset)

  call self%to_fieldset(vars, afieldset)
end subroutine

! ------------------------------------------------------------------------------
!> C++ interface for soca_increment_mod::soca_increment::from_fieldset()
subroutine soca_state_from_fieldset_c(c_key_self, c_vars, c_afieldset) &
  bind (c,name='soca_state_from_fieldset_f90')
integer(c_int),         intent(in) :: c_key_self
type(c_ptr),     value, intent(in) :: c_vars
type(c_ptr),     value, intent(in) :: c_afieldset

type(soca_state), pointer :: self
type(oops_variables)          :: vars
type(atlas_fieldset)          :: afieldset

call soca_state_registry%get(c_key_self, self)
vars = oops_variables(c_vars)
afieldset = atlas_fieldset(c_afieldset)

call self%from_fieldset(vars, afieldset)

end subroutine

end module soca_state_mod_c
