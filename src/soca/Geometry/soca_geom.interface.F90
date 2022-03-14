! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


!> C++ interfaces for soca_geom_mod::soca_geom
module soca_geom_mod_c

use atlas_module, only: atlas_fieldset, atlas_functionspace_pointcloud
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module,           only: fckit_mpi_comm
use iso_c_binding
use oops_variables_mod, only: oops_variables
use kinds, only: kind_real

! soca modules
use soca_geom_mod, only: soca_geom
use soca_fields_metadata_mod, only : soca_field_metadata


implicit none
private


#define LISTED_TYPE soca_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for soca_geom instances
type(registry_t), public :: soca_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"


! ------------------------------------------------------------------------------
!> C++ interface for soca_geom_mod::soca_geom::init()
subroutine soca_geo_setup_c(c_key_self, c_conf, c_comm) bind(c,name='soca_geo_setup_f90')
  integer(c_int),  intent(inout) :: c_key_self
  type(c_ptr),        intent(in) :: c_conf
  type(c_ptr), value, intent(in) :: c_comm

  type(soca_geom), pointer :: self

  call soca_geom_registry%init()
  call soca_geom_registry%add(c_key_self)
  call soca_geom_registry%get(c_key_self,self)

  call self%init(fckit_configuration(c_conf), fckit_mpi_comm(c_comm) )
end subroutine soca_geo_setup_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for soca_geom_mod::soca_geom::set_atlas_lonlat()
subroutine soca_geo_set_atlas_lonlat_c(c_key_self, c_afieldset)  bind(c,name='soca_geo_set_atlas_lonlat_f90')
  integer(c_int), intent(in) :: c_key_self
  type(c_ptr), intent(in), value :: c_afieldset

  type(soca_geom), pointer :: self
  type(atlas_fieldset) :: afieldset

  call soca_geom_registry%get(c_key_self,self)
  afieldset = atlas_fieldset(c_afieldset)

  call self%set_atlas_lonlat(afieldset)
end subroutine soca_geo_set_atlas_lonlat_c


! --------------------------------------------------------------------------------------------------
!> C++ interface to get atlas functionspace pointr from  soca_geom_mod::soca_geom
subroutine soca_geo_set_atlas_functionspace_pointer_c(c_key_self,c_afunctionspace) &
  bind(c,name='soca_geo_set_atlas_functionspace_pointer_f90')
  integer(c_int), intent(in)     :: c_key_self
  type(c_ptr), intent(in), value :: c_afunctionspace

  type(soca_geom),pointer :: self

  call soca_geom_registry%get(c_key_self,self)

  self%afunctionspace = atlas_functionspace_pointcloud(c_afunctionspace)
end subroutine soca_geo_set_atlas_functionspace_pointer_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for soca_geom_mod::soca_geom::fill_atlas_fieldset()
subroutine soca_geo_fill_atlas_fieldset_c(c_key_self, c_afieldset) &
 & bind(c,name='soca_geo_fill_atlas_fieldset_f90')

  integer(c_int),     intent(in) :: c_key_self
  type(c_ptr), value, intent(in) :: c_afieldset

  type(soca_geom), pointer :: self
  type(atlas_fieldset) :: afieldset

  call soca_geom_registry%get(c_key_self,self)
  afieldset = atlas_fieldset(c_afieldset)

  call self%fill_atlas_fieldset(afieldset)
end subroutine soca_geo_fill_atlas_fieldset_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_geom_mod::soca_geom::clone()
subroutine soca_geo_clone_c(c_key_self, c_key_other) bind(c,name='soca_geo_clone_f90')
  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in)    :: c_key_other

  type(soca_geom), pointer :: self, other

  call soca_geom_registry%add(c_key_self)
  call soca_geom_registry%get(c_key_self, self)
  call soca_geom_registry%get(c_key_other, other )

  call self%clone(other)
end subroutine soca_geo_clone_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_geom_mod::soca_geom::gridgen()
subroutine soca_geo_gridgen_c(c_key_self) bind(c,name='soca_geo_gridgen_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(soca_geom), pointer :: self

  call soca_geom_registry%get(c_key_self, self)

  call self%gridgen()
end subroutine soca_geo_gridgen_c


! ------------------------------------------------------------------------------
!> C++ interface for soca_geom_mod::soca_geom::end()
subroutine soca_geo_delete_c(c_key_self) bind(c,name='soca_geo_delete_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(soca_geom), pointer :: self

  call soca_geom_registry%get(c_key_self, self)
  call self%end()
  call soca_geom_registry%remove(c_key_self)
end subroutine soca_geo_delete_c


! ------------------------------------------------------------------------------
!> C++ interface to return begin and end of local geometry in soca_geom
subroutine soca_geo_start_end(c_key_self, ist, iend, jst, jend, kst, kend) bind(c, name='soca_geo_start_end_f90')
  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: ist, iend, jst, jend, kst, kend

  type(soca_geom), pointer :: self
  call soca_geom_registry%get(c_key_self, self)

  ist  = self%isc
  iend = self%iec
  jst  = self%jsc
  jend = self%jec
  kst  = 1
  kend = self%nzo
end subroutine soca_geo_start_end


! ------------------------------------------------------------------------------
!> C++ interface to get dimension of the GeometryIterator
subroutine soca_geo_iterator_dimension(c_key_self, itd) bind(c, name='soca_geo_iterator_dimension_f90')
  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: itd ! iterator dimension

  type(soca_geom), pointer :: self
  call soca_geom_registry%get(c_key_self, self)

  itd = self%iterator_dimension
end subroutine soca_geo_iterator_dimension

! ------------------------------------------------------------------------------
!> C++ interface to get number of levels for soca_geom
subroutine soca_geo_get_num_levels_c(c_key_self, c_vars, c_levels_size, c_levels) &
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
end subroutine soca_geo_get_num_levels_c


! ------------------------------------------------------------------------------
!> Get the number of points valid in the local grid.
subroutine soca_geo_gridsize_c(c_key_self, c_grid, c_masked, c_halo, c_size) &
           bind(c, name='soca_geo_gridsize_f90')
  integer(c_int),    intent(in) :: c_key_self
  character(c_char), intent(in) :: c_grid
  logical(c_bool),   intent(in) :: c_masked, c_halo
  integer(c_int),    intent(out) :: c_size

  type(soca_geom), pointer :: self
  real(kind=kind_real),     pointer :: mask(:,:) => null() !< field mask
  real(kind=kind_real),     pointer :: lon(:,:) => null()  !< field lon
  real(kind=kind_real),     pointer :: lat(:,:) => null()  !< field lat
  integer :: is, ie, js, je

  call soca_geom_registry%get(c_key_self, self)

  ! get the correct grid and mask
  select case(c_grid)
  case ('h')
    lon => self%lon
    lat => self%lat
    mask => self%mask2d
  case ('u')
    lon => self%lonu
    lat => self%latu
    mask => self%mask2du
  case ('v')
    lon => self%lonv
    lat => self%latv
    mask => self%mask2dv
  case default
    call abor1_ftn('error in soca_geo_gridsize_c. grid: '//c_grid)
  end select

  ! get the starting/ending index based on whether we need the halo
  if (c_halo) then
    is = self%isd; ie = self%ied
    js = self%jsd; je = self%jed
  else
    is = self%isc; ie = self%iec
    js = self%jsc; je = self%jec
  end if

  ! count number of point depending on whether masked
  if (c_masked) then
    c_size = count(mask(is:ie, js:je) /= 0)
  else
    c_size = (ie - is + 1) * (je - js + 1)
  end if
end subroutine soca_geo_gridsize_c


! ------------------------------------------------------------------------------
subroutine soca_geo_gridlatlon_c(c_key_self, c_grid, c_masked, c_halo, c_size, &
    c_lat, c_lon) bind(c, name='soca_geo_gridlatlon_f90')

  integer(c_int),    intent(in) :: c_key_self
  character(c_char), intent(in) :: c_grid
  logical(c_bool),   intent(in) :: c_masked, c_halo
  integer(c_int),    intent(in) :: c_size
  real(c_double),    intent(inout)  :: c_lat(c_size), c_lon(c_size)

  type(soca_geom), pointer :: self
  real(kind=kind_real),     pointer :: mask(:,:) => null() !< field mask
  real(kind=kind_real),     pointer :: lon(:,:) => null()  !< field lon
  real(kind=kind_real),     pointer :: lat(:,:) => null()  !< field lat
  integer :: is, ie, js, je, idx, i, j

  call soca_geom_registry%get(c_key_self, self)

  ! get the correct grid and mask
  select case(c_grid)
  case ('h')
    lon => self%lon
    lat => self%lat
    mask => self%mask2d
  case ('u')
    lon => self%lonu
    lat => self%latu
    mask => self%mask2du
  case ('v')
    lon => self%lonv
    lat => self%latv
    mask => self%mask2dv
  case default
    call abor1_ftn('error in soca_geo_gridlatlon_c. grid: '//c_grid)
  end select

  ! get the starting/ending index based on whether we need the halo
  if (c_halo) then
    is = self%isd; ie = self%ied
    js = self%jsd; je = self%jed
  else
    is = self%isc; ie = self%iec
    js = self%jsc; je = self%jec
  end if

  if (c_masked)   then
    c_lat(:) = pack(lat(is:ie,js:je), mask=mask(is:ie,js:je)/=0)
    c_lon(:) = pack(lon(is:ie,js:je), mask=mask(is:ie,js:je)/=0)
  else
    c_lat(:) = pack(lat(is:ie,js:je), mask=.true.)
    c_lon(:) = pack(lon(is:ie,js:je), mask=.true.)
  end if
end subroutine

! ------------------------------------------------------------------------------
subroutine soca_geo_getvargrid_c(c_key_self, c_vars, c_grid, c_masked) &
  bind(c, name="soca_geo_getvargrid_f90")
integer(c_int),     intent(in)  :: c_key_self
type(c_ptr), value, intent(in)  :: c_vars
character(c_char),  intent(inout) :: c_grid
logical(c_bool),    intent(inout) :: c_masked

type(soca_geom), pointer :: self
type(oops_variables)     :: vars
type(soca_field_metadata) :: metadata
character(len=:), allocatable :: var_name

call soca_geom_registry%get(c_key_self, self)
vars = oops_variables(c_vars)

! TODO cleanup
if (vars%nvars() /= 1) then
  call abor1_ftn('error in soca_geo_getvargrid_c. Wrong number of vars')
end if

var_name = trim(vars%variable(1))
metadata = self%fields_metadata%get(var_name)
c_grid = metadata%grid
c_masked = metadata%masked

end subroutine soca_geo_getvargrid_c

! ------------------------------------------------------------------------------
end module soca_geom_mod_c
