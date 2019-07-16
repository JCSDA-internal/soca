! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Interfaces to be called from C++ for Fortran handling of model fields

! ------------------------------------------------------------------------------
module soca_fields_mod_c
  use iso_c_binding
  use soca_fields
  use kinds
  use unstructured_grid_mod
  use datetime_mod
  use soca_interpfields_mod
  use soca_geom_mod_c
  use soca_geom_mod, only: soca_geom
  use ufo_locs_mod_c
  use ufo_locs_mod
  use ufo_geovals_mod_c
  use ufo_geovals_mod
  use variables_mod
  use soca_getvaltraj_mod, only: soca_getvaltraj, soca_getvaltraj_registry

  implicit none

  public :: soca_field_registry
  
#define LISTED_TYPE soca_field

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_field_registry

contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

subroutine soca_field_create_c(c_key_self, c_key_geom, c_vars) bind(c,name='soca_field_create_f90')
  integer(c_int), intent(inout) :: c_key_self !< Handle to field
  integer(c_int),    intent(in) :: c_key_geom !< Geometry
  type(c_ptr),       intent(in) :: c_vars     !< List of variables

  type(soca_field), pointer :: self
  type(soca_geom),  pointer :: geom
  type(oops_vars)           :: vars

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_field_registry%init()
  call soca_field_registry%add(c_key_self)
  call soca_field_registry%get(c_key_self,self)

  call oops_vars_create(c_vars, vars)
  call create(self, geom, vars)

end subroutine soca_field_create_c

! ------------------------------------------------------------------------------

subroutine soca_field_delete_c(c_key_self) bind(c,name='soca_field_delete_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(soca_field),     pointer :: self

  call soca_field_registry%get(c_key_self,self)
  call delete(self)
  call soca_field_registry%remove(c_key_self)

end subroutine soca_field_delete_c

! ------------------------------------------------------------------------------

subroutine soca_field_zero_c(c_key_self) bind(c,name='soca_field_zero_f90')
  integer(c_int), intent(in) :: c_key_self

  type(soca_field),  pointer :: self

  call soca_field_registry%get(c_key_self,self)
  call zeros(self)

end subroutine soca_field_zero_c

! ------------------------------------------------------------------------------

subroutine soca_field_dirac_c(c_key_self,c_conf) bind(c,name='soca_field_dirac_f90')
  integer(c_int), intent(in) :: c_key_self
  type(c_ptr),    intent(in) :: c_conf !< Configuration

  type(soca_field),  pointer :: self

  call soca_field_registry%get(c_key_self,self)
  call dirac(self,c_conf)

end subroutine soca_field_dirac_c

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

subroutine soca_field_random_c(c_key_self) bind(c,name='soca_field_random_f90')
  integer(c_int), intent(in) :: c_key_self

  type(soca_field),  pointer :: self

  call soca_field_registry%get(c_key_self,self)
  call random(self)

end subroutine soca_field_random_c

! ------------------------------------------------------------------------------

subroutine soca_field_copy_c(c_key_self,c_key_rhs) bind(c,name='soca_field_copy_f90')
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(soca_field), pointer :: self
  type(soca_field), pointer :: rhs

  call soca_field_registry%get(c_key_self,self)
  call soca_field_registry%get(c_key_rhs,rhs)

  call copy(self, rhs)

end subroutine soca_field_copy_c

! ------------------------------------------------------------------------------

subroutine soca_field_self_add_c(c_key_self,c_key_rhs) bind(c,name='soca_field_self_add_f90')
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(soca_field), pointer :: self
  type(soca_field), pointer :: rhs

  call soca_field_registry%get(c_key_self,self)
  call soca_field_registry%get(c_key_rhs,rhs)

  call self_add(self,rhs)

end subroutine soca_field_self_add_c

! ------------------------------------------------------------------------------

subroutine soca_field_self_schur_c(c_key_self,c_key_rhs) bind(c,name='soca_field_self_schur_f90')
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(soca_field), pointer :: self
  type(soca_field), pointer :: rhs

  call soca_field_registry%get(c_key_self,self)
  call soca_field_registry%get(c_key_rhs,rhs)

  call self_schur(self,rhs)

end subroutine soca_field_self_schur_c

! ------------------------------------------------------------------------------

subroutine soca_field_self_sub_c(c_key_self,c_key_rhs) bind(c,name='soca_field_self_sub_f90')
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(soca_field), pointer :: self
  type(soca_field), pointer :: rhs

  call soca_field_registry%get(c_key_self,self)
  call soca_field_registry%get(c_key_rhs,rhs)

  call self_sub(self,rhs)

end subroutine soca_field_self_sub_c

! ------------------------------------------------------------------------------

subroutine soca_field_self_mul_c(c_key_self,c_zz) bind(c,name='soca_field_self_mul_f90')
  integer(c_int), intent(in) :: c_key_self
  real(c_double), intent(in) :: c_zz

  type(soca_field), pointer :: self
  real(kind=kind_real) :: zz

  call soca_field_registry%get(c_key_self,self)
  zz = c_zz

  call self_mul(self,zz)

end subroutine soca_field_self_mul_c

! ------------------------------------------------------------------------------

subroutine soca_field_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='soca_field_axpy_f90')
  integer(c_int), intent(in) :: c_key_self
  real(c_double), intent(in) :: c_zz
  integer(c_int), intent(in) :: c_key_rhs

  type(soca_field), pointer :: self
  real(kind=kind_real)      :: zz
  type(soca_field), pointer :: rhs

  call soca_field_registry%get(c_key_self,self)
  call soca_field_registry%get(c_key_rhs,rhs)
  zz = c_zz

  call axpy(self,zz,rhs)

end subroutine soca_field_axpy_c

! ------------------------------------------------------------------------------

subroutine soca_field_dot_prod_c(c_key_fld1,c_key_fld2,c_prod) bind(c,name='soca_field_dot_prod_f90')
  integer(c_int),    intent(in) :: c_key_fld1, c_key_fld2
  real(c_double), intent(inout) :: c_prod

  type(soca_field), pointer :: fld1, fld2
  real(kind=kind_real) :: zz

  call soca_field_registry%get(c_key_fld1,fld1)
  call soca_field_registry%get(c_key_fld2,fld2)

  call dot_prod(fld1,fld2,zz)

  c_prod = zz

end subroutine soca_field_dot_prod_c

! ------------------------------------------------------------------------------

subroutine soca_field_add_incr_c(c_key_self,c_key_rhs) bind(c,name='soca_field_add_incr_f90')
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(soca_field), pointer :: self
  type(soca_field), pointer :: rhs

  call soca_field_registry%get(c_key_self,self)
  call soca_field_registry%get(c_key_rhs,rhs)

  call add_incr(self,rhs)

end subroutine soca_field_add_incr_c

! ------------------------------------------------------------------------------

subroutine soca_field_diff_incr_c(c_key_lhs,c_key_x1,c_key_x2) bind(c,name='soca_field_diff_incr_f90')
  integer(c_int), intent(in) :: c_key_lhs
  integer(c_int), intent(in) :: c_key_x1
  integer(c_int), intent(in) :: c_key_x2

  type(soca_field), pointer :: lhs
  type(soca_field), pointer :: x1
  type(soca_field), pointer :: x2

  call soca_field_registry%get(c_key_lhs,lhs)
  call soca_field_registry%get(c_key_x1,x1)
  call soca_field_registry%get(c_key_x2,x2)

  call diff_incr(lhs,x1,x2)

end subroutine soca_field_diff_incr_c

! ------------------------------------------------------------------------------

subroutine soca_field_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='soca_field_change_resol_f90')
  integer(c_int), intent(in) :: c_key_fld
  integer(c_int), intent(in) :: c_key_rhs

  type(soca_field), pointer :: fld, rhs

  call soca_field_registry%get(c_key_fld,fld)
  call soca_field_registry%get(c_key_rhs,rhs)

  call change_resol(fld,rhs)

end subroutine soca_field_change_resol_c

! ------------------------------------------------------------------------------

subroutine soca_field_ug_coord_c(c_key_fld, c_key_ug) bind (c,name='soca_field_ug_coord_f90')
  integer(c_int), intent(in) :: c_key_fld
  integer(c_int), intent(in) :: c_key_ug

  type(soca_field), pointer :: fld
  type(unstructured_grid), pointer :: ug

  call soca_field_registry%get(c_key_fld,fld)
  call unstructured_grid_registry%get(c_key_ug,ug)

  call ug_coord(fld, ug)

end subroutine soca_field_ug_coord_c

! ------------------------------------------------------------------------------

subroutine soca_field_field_to_ug_c(c_key_fld, c_key_ug, c_its) bind (c,name='soca_field_field_to_ug_f90')
  integer(c_int), intent(in) :: c_key_fld
  integer(c_int), intent(in) :: c_key_ug
  integer(c_int), intent(in) :: c_its

  type(soca_field),        pointer :: fld
  type(unstructured_grid), pointer :: ug
  integer                          :: its

  call soca_field_registry%get(c_key_fld,fld)
  call unstructured_grid_registry%get(c_key_ug,ug)
  its = c_its+1

  call field_to_ug(fld, ug, its)

end subroutine soca_field_field_to_ug_c

! ------------------------------------------------------------------------------

subroutine soca_field_field_from_ug_c(c_key_fld, c_key_ug, c_its) bind (c,name='soca_field_field_from_ug_f90')
  integer(c_int), intent(in) :: c_key_fld
  integer(c_int), intent(in) :: c_key_ug
  integer(c_int), intent(in) :: c_its

  type(soca_field),        pointer :: fld
  type(unstructured_grid), pointer :: ug
  integer                          :: its

  call soca_field_registry%get(c_key_fld,fld)
  call unstructured_grid_registry%get(c_key_ug,ug)
  its = c_its+1

  call field_from_ug(fld, ug, its)

end subroutine soca_field_field_from_ug_c

! ------------------------------------------------------------------------------

subroutine soca_field_read_file_c(c_key_fld, c_conf, c_dt) bind(c,name='soca_field_read_file_f90')
  integer(c_int), intent(in) :: c_key_fld  !< Fields
  type(c_ptr),    intent(in) :: c_conf     !< Configuration
  type(c_ptr), intent(inout) :: c_dt       !< DateTime

  type(soca_field), pointer :: fld
  type(datetime)            :: fdate

  call soca_field_registry%get(c_key_fld,fld)
  call c_f_datetime(c_dt, fdate)
  call read_file(fld, c_conf, fdate)

end subroutine soca_field_read_file_c

! ------------------------------------------------------------------------------

subroutine soca_field_write_file_c(c_key_fld, c_conf, c_dt) bind(c,name='soca_field_write_file_f90')
  integer(c_int), intent(in) :: c_key_fld  !< Fields
  type(c_ptr),    intent(in) :: c_conf     !< Configuration
  type(c_ptr),    intent(in) :: c_dt       !< DateTime

  type(soca_field), pointer :: fld
  type(datetime)            :: fdate

  call soca_field_registry%get(c_key_fld,fld)
  call c_f_datetime(c_dt, fdate)
  call write_file(fld, c_conf, fdate)

end subroutine soca_field_write_file_c

! ------------------------------------------------------------------------------

subroutine soca_field_gpnorm_c(c_key_fld, kf, pstat) bind(c,name='soca_field_gpnorm_f90')
  integer(c_int),    intent(in) :: c_key_fld
  integer(c_int),    intent(in) :: kf
  real(c_double), intent(inout) :: pstat(3*kf)

  type(soca_field), pointer :: fld
  real(kind=kind_real)      :: zstat(3, kf)
  integer :: jj, js, jf

  call soca_field_registry%get(c_key_fld,fld)

  call gpnorm(fld, kf, zstat)
  jj=0
  do jf = 1, kf
     do js = 1, 3
        jj=jj+1
        pstat(jj) = zstat(js,jf)
     enddo
  enddo

end subroutine soca_field_gpnorm_c

! ------------------------------------------------------------------------------

subroutine soca_field_rms_c(c_key_fld, prms) bind(c,name='soca_field_rms_f90')
  integer(c_int),    intent(in) :: c_key_fld
  real(c_double), intent(inout) :: prms

  type(soca_field), pointer :: fld
  real(kind=kind_real)      :: zz

  call soca_field_registry%get(c_key_fld,fld)

  call fldrms(fld, zz)

  prms = zz

end subroutine soca_field_rms_c

! ------------------------------------------------------------------------------

subroutine soca_field_interp_nl_c(c_key_fld,c_key_loc,c_vars,c_key_gom) bind(c,name='soca_field_interp_nl_f90')
  integer(c_int), intent(in) :: c_key_fld
  integer(c_int), intent(in) :: c_key_loc
  type(c_ptr),    intent(in) :: c_vars     !< List of requested variables
  integer(c_int), intent(in) :: c_key_gom

  type(soca_field),  pointer :: fld
  type(ufo_locs),    pointer :: locs
  type(ufo_geovals), pointer :: gom
  type(oops_vars)            :: vars

  call oops_vars_create(c_vars, vars)

  call soca_field_registry%get(c_key_fld,fld)
  call ufo_locs_registry%get(c_key_loc,locs)
  call ufo_geovals_registry%get(c_key_gom,gom)

  call getvalues(fld, locs, vars, gom)

end subroutine soca_field_interp_nl_c

! ------------------------------------------------------------------------------

subroutine soca_field_interp_nl_traj_c(c_key_fld,c_key_loc,c_vars,c_key_gom,c_key_traj) bind(c,name='soca_field_interp_nl_traj_f90')
  integer(c_int),           intent(in) :: c_key_fld
  integer(c_int),           intent(in) :: c_key_loc
  type(c_ptr),              intent(in) :: c_vars     !< List of requested variables
  integer(c_int),           intent(in) :: c_key_gom
  integer(c_int), optional, intent(in) :: c_key_traj !< Trajectory for interpolation/transforms

  type(soca_field),      pointer :: fld
  type(ufo_locs),        pointer :: locs
  type(oops_vars)                :: vars
  type(ufo_geovals),     pointer :: gom
  type(soca_getvaltraj), pointer :: traj

  call oops_vars_create(c_vars, vars)

  call soca_field_registry%get(c_key_fld,fld)
  call ufo_locs_registry%get(c_key_loc,locs)
  call ufo_geovals_registry%get(c_key_gom,gom)
  call soca_getvaltraj_registry%get(c_key_traj,traj)

  call getvalues(fld, locs, vars, gom, traj, interp_type='nl')

end subroutine soca_field_interp_nl_traj_c

! ------------------------------------------------------------------------------

subroutine soca_field_interp_tl_c(c_key_fld,c_key_loc,c_vars,c_key_gom,c_key_traj) bind(c,name='soca_field_interp_tl_f90')
  integer(c_int),           intent(in) :: c_key_fld
  integer(c_int),           intent(in) :: c_key_loc
  type(c_ptr),              intent(in) :: c_vars     !< List of requested variables
  integer(c_int),           intent(in) :: c_key_gom
  integer(c_int), optional, intent(in) :: c_key_traj !< Trajectory for interpolation/transforms

  type(soca_field),      pointer :: fld
  type(ufo_locs),        pointer :: locs
  type(oops_vars)                :: vars
  type(ufo_geovals),     pointer :: gom
  type(soca_getvaltraj), pointer :: traj

  call oops_vars_create(c_vars, vars)

  call soca_field_registry%get(c_key_fld,fld)
  call ufo_locs_registry%get(c_key_loc,locs)
  call ufo_geovals_registry%get(c_key_gom,gom)
  call soca_getvaltraj_registry%get(c_key_traj,traj)

  call getvalues(fld, locs, vars, gom, traj, interp_type='tl')

end subroutine soca_field_interp_tl_c

! ------------------------------------------------------------------------------

subroutine soca_field_interp_ad_c(c_key_fld,c_key_loc,c_vars,c_key_gom,c_key_traj) bind(c,name='soca_field_interp_ad_f90')
  integer(c_int), intent(in) :: c_key_fld
  integer(c_int), intent(in) :: c_key_loc
  type(c_ptr),    intent(in) :: c_vars     !< List of requested variables
  integer(c_int), intent(in) :: c_key_gom
  integer(c_int), intent(in) :: c_key_traj !< Trajectory for interpolation/transforms

  type(soca_field),      pointer :: fld
  type(ufo_locs),        pointer :: locs
  type(oops_vars)                :: vars
  type(ufo_geovals),     pointer :: gom
  type(soca_getvaltraj), pointer :: traj

  call oops_vars_create(c_vars, vars)

  call soca_field_registry%get(c_key_fld,fld)
  call ufo_locs_registry%get(c_key_loc,locs)
  call ufo_geovals_registry%get(c_key_gom,gom)
  call soca_getvaltraj_registry%get(c_key_traj,traj)

  call getvalues_ad(fld, locs, vars, gom, traj)

end subroutine soca_field_interp_ad_c

! ------------------------------------------------------------------------------

subroutine soca_fieldnum_c(c_key_fld, nx, ny, nzo, nzi, ncat, nf) bind(c,name='soca_field_sizes_f90')
  integer(c_int),         intent(in) :: c_key_fld
  integer(kind=c_int), intent(inout) :: nx, ny, nzo, nzi, ncat, nf

  type(soca_field), pointer :: fld

  call soca_field_registry%get(c_key_fld,fld)

  nx = fld%geom%nx
  ny = fld%geom%ny
  nzo = fld%geom%nzo
  nzi = fld%geom%nzi
  ncat = fld%geom%ncat
  nf = fld%nf

end subroutine soca_fieldnum_c

end module soca_fields_mod_c
