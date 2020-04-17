! (C) Copyright 2020-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_increment_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use unstructured_grid_mod, only: unstructured_grid, unstructured_grid_registry
use datetime_mod, only: datetime, c_f_datetime
use oops_variables_mod
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_locs_mod, only: ufo_locs
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod, only: ufo_geovals
use soca_geom_mod, only: soca_geom
use soca_geom_mod_c, only: soca_geom_registry
use soca_geom_iter_mod, only: soca_geom_iter, soca_geom_iter_registry
use soca_increment_mod
use soca_increment_reg
use soca_state_mod
use soca_state_reg

implicit none
private

contains

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

  subroutine soca_increment_delete_c(c_key_self) bind(c,name='soca_increment_delete_f90')
    integer(c_int), intent(inout) :: c_key_self

    type(soca_increment),    pointer :: self

    call soca_increment_registry%get(c_key_self,self)
    call self%delete( )
    call soca_increment_registry%remove(c_key_self)

  end subroutine soca_increment_delete_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_zero_c(c_key_self) bind(c,name='soca_increment_zero_f90')
    integer(c_int), intent(in) :: c_key_self

    type(soca_increment), pointer :: self

    call soca_increment_registry%get(c_key_self,self)
    call self%zeros()

  end subroutine soca_increment_zero_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_dirac_c(c_key_self,c_conf) bind(c,name='soca_increment_dirac_f90')
    integer(c_int), intent(in) :: c_key_self
    type(c_ptr),    intent(in) :: c_conf !< Configuration

    type(soca_increment), pointer :: self

    call soca_increment_registry%get(c_key_self,self)
    call self%dirac(fckit_configuration(c_conf))

  end subroutine soca_increment_dirac_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_random_c(c_key_self) bind(c,name='soca_increment_random_f90')
    integer(c_int), intent(in) :: c_key_self

    type(soca_increment), pointer :: self

    call soca_increment_registry%get(c_key_self,self)
    call self%random()

  end subroutine soca_increment_random_c

  ! ------------------------------------------------------------------------------

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

  subroutine soca_increment_self_add_c(c_key_self,c_key_rhs) bind(c,name='soca_increment_self_add_f90')
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_rhs

    type(soca_increment), pointer :: self
    type(soca_increment), pointer :: rhs

    call soca_increment_registry%get(c_key_self,self)
    call soca_increment_registry%get(c_key_rhs,rhs)

    call self%add(rhs)

  end subroutine soca_increment_self_add_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_self_schur_c(c_key_self,c_key_rhs) bind(c,name='soca_increment_self_schur_f90')
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_rhs

    type(soca_increment), pointer :: self
    type(soca_increment), pointer :: rhs

    call soca_increment_registry%get(c_key_self,self)
    call soca_increment_registry%get(c_key_rhs,rhs)

    call self%schur(rhs)

  end subroutine soca_increment_self_schur_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_self_sub_c(c_key_self,c_key_rhs) bind(c,name='soca_increment_self_sub_f90')
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_rhs

    type(soca_increment), pointer :: self
    type(soca_increment), pointer :: rhs

    call soca_increment_registry%get(c_key_self,self)
    call soca_increment_registry%get(c_key_rhs,rhs)

    call self%sub(rhs)

  end subroutine soca_increment_self_sub_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_self_mul_c(c_key_self,c_zz) bind(c,name='soca_increment_self_mul_f90')
    integer(c_int), intent(in) :: c_key_self
    real(c_double), intent(in) :: c_zz

    type(soca_increment), pointer :: self
    real(kind=kind_real) :: zz

    call soca_increment_registry%get(c_key_self,self)
    zz = c_zz

    call self%mul(zz)

  end subroutine soca_increment_self_mul_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='soca_increment_axpy_f90')
    integer(c_int), intent(in) :: c_key_self
    real(c_double), intent(in) :: c_zz
    integer(c_int), intent(in) :: c_key_rhs

    type(soca_increment), pointer :: self
    real(kind=kind_real)      :: zz
    type(soca_increment), pointer :: rhs

    call soca_increment_registry%get(c_key_self,self)
    call soca_increment_registry%get(c_key_rhs,rhs)
    zz = c_zz

    call self%axpy(zz,rhs)

  end subroutine soca_increment_axpy_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_dot_prod_c(c_key_fld1,c_key_fld2,c_prod) bind(c,name='soca_increment_dot_prod_f90')
    integer(c_int),    intent(in) :: c_key_fld1, c_key_fld2
    real(c_double), intent(inout) :: c_prod

    type(soca_increment), pointer :: fld1, fld2
    real(kind=kind_real) :: zz

    call soca_increment_registry%get(c_key_fld1,fld1)
    call soca_increment_registry%get(c_key_fld2,fld2)

    call fld1%dot_prod(fld2,zz)

    c_prod = zz

  end subroutine soca_increment_dot_prod_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_diff_incr_c(c_key_lhs,c_key_x1,c_key_x2) bind(c,name='soca_increment_diff_incr_f90')
    integer(c_int), intent(in) :: c_key_lhs
    integer(c_int), intent(in) :: c_key_x1
    integer(c_int), intent(in) :: c_key_x2

    type(soca_increment), pointer :: lhs
    type(soca_state),     pointer :: x1
    type(soca_state),     pointer :: x2

    call soca_increment_registry%get(c_key_lhs,lhs)
    call soca_state_registry%get(c_key_x1,x1)
    call soca_state_registry%get(c_key_x2,x2)
    call x1%diff_incr(x2, lhs)

  end subroutine soca_increment_diff_incr_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='soca_increment_change_resol_f90')
    integer(c_int), intent(in) :: c_key_fld
    integer(c_int), intent(in) :: c_key_rhs

    type(soca_increment), pointer :: fld, rhs

    call soca_increment_registry%get(c_key_fld,fld)
    call soca_increment_registry%get(c_key_rhs,rhs)

    ! TODO implement a proper change of resolution, just copying for now
    call fld%copy(rhs)

  end subroutine soca_increment_change_resol_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_ug_coord_c(c_key_fld, c_key_ug) bind (c,name='soca_increment_ug_coord_f90')
    integer(c_int), intent(in) :: c_key_fld
    integer(c_int), intent(in) :: c_key_ug

    type(soca_increment), pointer :: fld
    type(unstructured_grid), pointer :: ug

    call soca_increment_registry%get(c_key_fld,fld)
    call unstructured_grid_registry%get(c_key_ug,ug)

    call fld%ug_coord(ug)

  end subroutine soca_increment_ug_coord_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_field_to_ug_c(c_key_fld, c_key_ug, c_its) bind (c,name='soca_increment_field_to_ug_f90')
    integer(c_int), intent(in) :: c_key_fld
    integer(c_int), intent(in) :: c_key_ug
    integer(c_int), intent(in) :: c_its

    type(soca_increment),       pointer :: fld
    type(unstructured_grid), pointer :: ug
    integer                          :: its

    call soca_increment_registry%get(c_key_fld,fld)
    call unstructured_grid_registry%get(c_key_ug,ug)
    its = c_its+1

    call fld%to_ug(ug, its)

  end subroutine soca_increment_field_to_ug_c

  ! ------------------------------------------------------------------------------

  subroutine soca_increment_field_from_ug_c(c_key_fld, c_key_ug, c_its) bind (c,name='soca_increment_field_from_ug_f90')
    integer(c_int), intent(in) :: c_key_fld
    integer(c_int), intent(in) :: c_key_ug
    integer(c_int), intent(in) :: c_its

    type(soca_increment),       pointer :: fld
    type(unstructured_grid), pointer :: ug
    integer                          :: its

    call soca_increment_registry%get(c_key_fld,fld)
    call unstructured_grid_registry%get(c_key_ug,ug)
    its = c_its+1

    call fld%from_ug(ug, its)

  end subroutine soca_increment_field_from_ug_c

  ! ------------------------------------------------------------------------------

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

  subroutine soca_increment_gpnorm_c(c_key_fld, kf, pstat) bind(c,name='soca_increment_gpnorm_f90')
    integer(c_int),    intent(in) :: c_key_fld
    integer(c_int),    intent(in) :: kf
    real(c_double), intent(inout) :: pstat(3*kf)

    type(soca_increment), pointer :: fld
    real(kind=kind_real)      :: zstat(3, kf)
    integer :: jj, js, jf

    call soca_increment_registry%get(c_key_fld,fld)

    call fld%gpnorm(kf, zstat)
    jj=0
    do jf = 1, kf
       do js = 1, 3
          jj=jj+1
          pstat(jj) = zstat(js,jf)
       enddo
    enddo

  end subroutine soca_increment_gpnorm_c

  ! ------------------------------------------------------------------------------

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

  subroutine soca_incrementnum_c(c_key_fld, nx, ny, nzo, nzi, ncat, nf) bind(c,name='soca_increment_sizes_f90')
    integer(c_int),         intent(in) :: c_key_fld
    integer(kind=c_int), intent(inout) :: nx, ny, nzo, nzi, ncat, nf

    type(soca_increment), pointer :: fld

    call soca_increment_registry%get(c_key_fld,fld)

    nx = size(fld%geom%lon,1)
    ny = size(fld%geom%lon,2)
    nzo = fld%geom%nzo
    nzi = fld%geom%nzi
    ncat = fld%geom%ncat
    nf = size(fld%fields)

  end subroutine soca_incrementnum_c

end module
