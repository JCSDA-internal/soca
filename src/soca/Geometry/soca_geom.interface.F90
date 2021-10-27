! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_geom_mod_c

use atlas_module, only: atlas_fieldset, atlas_functionspace_pointcloud
use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module,           only: fckit_mpi_comm
use soca_geom_mod, only: soca_geom
use oops_variables_mod
use soca_fields_metadata_mod

implicit none

private
public :: soca_geom_registry

#define LISTED_TYPE soca_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: soca_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
!> Setup geometry object
subroutine c_soca_geo_setup(c_key_self, c_conf, c_comm) bind(c,name='soca_geo_setup_f90')

  integer(c_int),  intent(inout) :: c_key_self
  type(c_ptr),        intent(in) :: c_conf
  type(c_ptr), value, intent(in) :: c_comm

  type(soca_geom), pointer :: self

  call soca_geom_registry%init()
  call soca_geom_registry%add(c_key_self)
  call soca_geom_registry%get(c_key_self,self)

  call self%init(fckit_configuration(c_conf), fckit_mpi_comm(c_comm) )

end subroutine c_soca_geo_setup

! --------------------------------------------------------------------------------------------------
!> Set ATLAS lonlat fieldset
subroutine c_soca_geo_set_atlas_lonlat(c_key_self, c_afieldset)  bind(c,name='soca_geo_set_atlas_lonlat_f90')

  integer(c_int), intent(in) :: c_key_self
  type(c_ptr), intent(in), value :: c_afieldset

  type(soca_geom), pointer :: self
  type(atlas_fieldset) :: afieldset

  call soca_geom_registry%get(c_key_self,self)
  afieldset = atlas_fieldset(c_afieldset)

  call self%set_atlas_lonlat(afieldset)

end subroutine c_soca_geo_set_atlas_lonlat

! --------------------------------------------------------------------------------------------------
!> Set ATLAS functionspace pointer
subroutine c_soca_geo_set_atlas_functionspace_pointer(c_key_self,c_afunctionspace) &
 & bind(c,name='soca_geo_set_atlas_functionspace_pointer_f90')

  integer(c_int), intent(in)     :: c_key_self
  type(c_ptr), intent(in), value :: c_afunctionspace

  type(soca_geom),pointer :: self

  call soca_geom_registry%get(c_key_self,self)

  self%afunctionspace = atlas_functionspace_pointcloud(c_afunctionspace)

end subroutine c_soca_geo_set_atlas_functionspace_pointer

! --------------------------------------------------------------------------------------------------
!> Fill ATLAS fieldset
subroutine c_soca_geo_fill_atlas_fieldset(c_key_self, c_afieldset) &
 & bind(c,name='soca_geo_fill_atlas_fieldset_f90')

  integer(c_int),     intent(in) :: c_key_self
  type(c_ptr), value, intent(in) :: c_afieldset

  type(soca_geom), pointer :: self
  type(atlas_fieldset) :: afieldset

  call soca_geom_registry%get(c_key_self,self)
  afieldset = atlas_fieldset(c_afieldset)

  call self%fill_atlas_fieldset(afieldset)

end subroutine c_soca_geo_fill_atlas_fieldset

! ------------------------------------------------------------------------------
!> Clone geometry object
subroutine c_soca_geo_clone(c_key_self, c_key_other) bind(c,name='soca_geo_clone_f90')

  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in)    :: c_key_other

  type(soca_geom), pointer :: self, other

  call soca_geom_registry%add(c_key_self)
  call soca_geom_registry%get(c_key_self, self)
  call soca_geom_registry%get(c_key_other, other )

  call self%clone(other)

end subroutine c_soca_geo_clone

! ------------------------------------------------------------------------------
!> Generate grid
subroutine c_soca_geo_gridgen(c_key_self) bind(c,name='soca_geo_gridgen_f90')

  integer(c_int), intent(inout) :: c_key_self

  type(soca_geom), pointer :: self

  call soca_geom_registry%get(c_key_self, self)

  call self%gridgen()

end subroutine c_soca_geo_gridgen

! ------------------------------------------------------------------------------
!> Geometry destructor
subroutine c_soca_geo_delete(c_key_self) bind(c,name='soca_geo_delete_f90')

  integer(c_int), intent(inout) :: c_key_self

  type(soca_geom), pointer :: self

  call soca_geom_registry%get(c_key_self, self)
  call self%end()
  call soca_geom_registry%remove(c_key_self)

end subroutine c_soca_geo_delete

! ------------------------------------------------------------------------------
!> return begin and end of local geometry
subroutine c_soca_geo_start_end(c_key_self, ist, iend, jst, jend, kst, kend) bind(c, name='soca_geo_start_end_f90')

  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: ist, iend, jst, jend, kst, kend

  type(soca_geom), pointer :: self
  call soca_geom_registry%get(c_key_self, self)

  ist  = self%isc
  iend = self%iec
  jst  = self%jsc
  jend = self%jec
  kst  = self%ksc
  kend = self%kec

end subroutine c_soca_geo_start_end

! ------------------------------------------------------------------------------
!> return dimension of the GeometryIterator
subroutine c_soca_geo_iterator_dimension(c_key_self, itd) bind(c, name='soca_geo_iterator_dimension_f90')
  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: itd ! iterator dimension

  type(soca_geom), pointer :: self
  call soca_geom_registry%get(c_key_self, self)

  itd = self%iterator_dimension
end subroutine

! ------------------------------------------------------------------------------
subroutine c_soca_geo_get_num_levels(c_key_self, c_vars, c_levels_size, c_levels) &
           bind(c, name='soca_geo_get_num_levels_f90')
  integer(c_int),     intent(in)  :: c_key_self
  type(c_ptr), value, intent(in)  :: c_vars
  integer(c_size_t),  intent(in)  :: c_levels_size
  integer(c_size_t),  intent(out) :: c_levels(c_levels_size)

  type(soca_geom), pointer :: self
  type(oops_variables)     :: vars
  integer :: i
  character(len=:), allocatable :: field_name
  type(soca_field_metadata) :: field

  call soca_geom_registry%get(c_key_self, self)
  vars = oops_variables(c_vars)

  do i = 1,vars%nvars()
    field_name = vars%variable(i)
    field = self%fields_metadata%get(field_name)
    select case(field%levels)
    case ("1")
      c_levels(i) = 1
    case ("full_ocn")
      if (field_name == field%getval_name_surface) then
        c_levels(i) = 1
      else
        c_levels(i) = self%nzo
      end if
    case default
      call abor1_ftn('ERROR in c_soca_geo_get_num_levels, unknown "levels" '//field%levels)
    end select
  end do

end subroutine c_soca_geo_get_num_levels

! ------------------------------------------------------------------------------

end module soca_geom_mod_c
