! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Interfaces to be called from C++ for Fortran handling of model fields

! ------------------------------------------------------------------------------
module soca_state_mod_c

use iso_c_binding

use datetime_mod, only: datetime, c_f_datetime
use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use oops_variables_mod
use soca_geom_mod_c, only: soca_geom_registry
use soca_geom_mod, only: soca_geom
use soca_increment_mod
use soca_increment_reg
use soca_state_mod
use soca_state_reg
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod, only: ufo_geovals

implicit none
private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


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

subroutine soca_state_delete_c(c_key_self) bind(c,name='soca_state_delete_f90')
    integer(c_int), intent(inout) :: c_key_self

    type(soca_state),    pointer :: self

    call soca_state_registry%get(c_key_self,self)
    call self%delete( )
    call soca_state_registry%remove(c_key_self)

end subroutine soca_state_delete_c

! ------------------------------------------------------------------------------

subroutine soca_state_zero_c(c_key_self) bind(c,name='soca_state_zero_f90')
    integer(c_int), intent(in) :: c_key_self

    type(soca_state), pointer :: self

    call soca_state_registry%get(c_key_self,self)
    call self%zeros()

end subroutine soca_state_zero_c


! ------------------------------------------------------------------------------

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
subroutine soca_state_expontrans_c(c_key_self, c_trvars) bind(c,name='soca_state_expontrans_f90')
  integer(c_int), intent(in)     :: c_key_self
  type(c_ptr), value, intent(in) :: c_trvars

  type(soca_state), pointer :: self
  type(oops_variables) :: trvars

  trvars = oops_variables(c_trvars)

  call soca_state_registry%get(c_key_self,self)
  call self%logexpon(transfunc="expon", trvars=trvars)

end subroutine soca_state_expontrans_c

end module soca_state_mod_c
