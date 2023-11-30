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

use mpp_domains_mod, only : mpp_update_domains

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
subroutine soca_geo_lonlat_c(c_key_self, c_afieldset)  bind(c,name='soca_geo_lonlat_f90')
  integer(c_int), intent(in) :: c_key_self
  type(c_ptr), intent(in), value :: c_afieldset

  type(soca_geom), pointer :: self
  type(atlas_fieldset) :: afieldset

  call soca_geom_registry%get(c_key_self,self)
  afieldset = atlas_fieldset(c_afieldset)

  call self%lonlat(afieldset)
end subroutine soca_geo_lonlat_c


! --------------------------------------------------------------------------------------------------
!> C++ interface to get atlas functionspace pointr from  soca_geom_mod::soca_geom
subroutine soca_geo_set_atlas_functionspace_pointer_c(c_key_self, c_functionspaceIncHalo) &
  bind(c,name='soca_geo_set_atlas_functionspace_pointer_f90')
  integer(c_int), intent(in)     :: c_key_self
  type(c_ptr), intent(in), value :: c_functionspaceIncHalo

  type(soca_geom),pointer :: self

  call soca_geom_registry%get(c_key_self,self)

  self%functionspaceIncHalo = atlas_functionspace_pointcloud(c_functionspaceIncHalo)
end subroutine soca_geo_set_atlas_functionspace_pointer_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for soca_geom_mod::soca_geom::fill_atlas_fieldset()
subroutine soca_geo_to_fieldset_c(c_key_self, c_afieldset) &
 & bind(c,name='soca_geo_to_fieldset_f90')

  integer(c_int),     intent(in) :: c_key_self
  type(c_ptr), value, intent(in) :: c_afieldset

  type(soca_geom), pointer :: self
  type(atlas_fieldset) :: afieldset

  call soca_geom_registry%get(c_key_self,self)
  afieldset = atlas_fieldset(c_afieldset)

  call self%to_fieldset(afieldset)
end subroutine soca_geo_to_fieldset_c


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
!> Get the number of points valid in the local grid (skipping masked points)
subroutine soca_geo_gridsize_c(c_key_self, c_halo, c_size) &
           bind(c, name='soca_geo_gridsize_f90')
  integer(c_int),    intent(in) :: c_key_self
  logical(c_bool),   intent(in) :: c_halo    !< true if halo should be included in number of points
  integer(c_int),    intent(out):: c_size    !< the resulting number of gridpoints

  type(soca_geom), pointer :: self

  call soca_geom_registry%get(c_key_self, self)
  if (c_halo) then
    c_size = self%ngrid_halo_valid
  else
    c_size = self%ngrid
  end if
end subroutine soca_geo_gridsize_c

! ------------------------------------------------------------------------------
!> Get the lat/lons for the h grid
!> If c_halo is true, then only the "valid" halo points are used
!> (i.e. duplicate halo points and invalid halo points are skipped)
subroutine soca_geo_gridlatlon_c(c_key_self, c_halo, c_size, &
    c_lat, c_lon) bind(c, name='soca_geo_gridlatlon_f90')

  integer(c_int),    intent(in) :: c_key_self
  logical(c_bool),   intent(in) :: c_halo
  integer(c_int),    intent(in) :: c_size
  real(c_double),    intent(inout)  :: c_lat(c_size), c_lon(c_size)

  type(soca_geom), pointer :: self
  real(kind=kind_real),     pointer :: lon(:,:) => null()  !< field lon
  real(kind=kind_real),     pointer :: lat(:,:) => null()  !< field lat

  call soca_geom_registry%get(c_key_self, self)
  
  ! TODO: re-enable the ability to get lat/lon for different grids?
  ! For now just use h grid
  lon => self%lon
  lat => self%lat
  
  if (c_halo) then
    ! get all lat/lon points, except for duplicate or invalid halo points
    c_lat(:) = pack(self%lat, mask=self%valid_halo_mask)
    c_lon(:) = pack(self%lon, mask=self%valid_halo_mask)
  else
    ! get all non-halo lat/lon points
    c_lat(:) = pack(self%lat(self%isc:self%iec, self%jsc:self%jec), .true.)
    c_lon(:) = pack(self%lon(self%isc:self%iec, self%jsc:self%jec), .true.)
  end if
end subroutine

! ------------------------------------------------------------------------------
subroutine soca_geo_get_mesh_size_c(c_key_self, c_nodes, c_quads) bind(c, name='soca_geo_get_mesh_size_f90')
  integer(c_int), intent(in)  :: c_key_self
  integer(c_int), intent(out) :: c_nodes, c_quads

  type(soca_geom), pointer :: self
  logical, allocatable :: valid_nodes(:,:), valid_cells(:,:)

  call soca_geom_registry%get(c_key_self, self)
  
  call self%mesh_valid_nodes_cells(valid_nodes, valid_cells)
  c_nodes = count(valid_nodes)
  c_quads = count(valid_cells)

end subroutine

! ------------------------------------------------------------------------------
! TODO move this into the soca_geom_mod file?
subroutine soca_geo_get_mesh_c(c_key_self, c_nodes, c_lon, c_lat, c_ghosts, c_global_idx, c_remote_idx, c_partition, &
    c_quads) bind(c, name='soca_geo_get_mesh_f90')
  integer(c_int), intent(in)    :: c_key_self
  integer(c_int), intent(in)    :: c_nodes, c_quads
  integer(c_int), intent(inout) :: c_ghosts(c_nodes), c_global_idx(c_nodes), c_remote_idx(c_nodes), c_partition(c_nodes)
  real(c_double), intent(inout) :: c_lon(c_nodes), c_lat(c_nodes)

  !logical :: N_tripolar, EW_cyclic
  integer :: idx, i, j
  integer :: nx, ny
  integer, allocatable :: global_idx(:,:), local_idx(:,:), partition(:,:)
  logical, allocatable :: valid_nodes(:,:), valid_cells(:,:)
  
  type(soca_geom), pointer :: self

  call soca_geom_registry%get(c_key_self, self)

  ! parameters pulled from grid
  nx=self%domain%NIGLOBAL
  ny=self%domain%NJGLOBAL

  ! find global/local indexes, and the PE owner of each grid point
  ! (requiring a halo communications)
  allocate(global_idx(self%isd:self%ied, self%jsd:self%jed))
  allocate(local_idx(self%isd:self%ied, self%jsd:self%jed))
  allocate(partition(self%isd:self%ied, self%jsd:self%jed))
  idx=0 ! local index
  global_idx = -999
  local_idx = -999
  partition = -999
  do j=self%jsc,self%jec
    do i=self%isc, self%iec
      idx = idx + 1
      local_idx(i,j) = idx
      global_idx(i,j) = (j-1)*nx + i
      partition(i,j) = self%f_comm%rank()
    end do
  end do
  call mpp_update_domains(local_idx, self%Domain%mpp_domain)
  call mpp_update_domains(global_idx, self%Domain%mpp_domain)
  call mpp_update_domains(partition, self%Domain%mpp_domain)

  ! find which quads / vertices are we going to skip (in case of non cyclic or tripolar_fold special cases)
  call self%mesh_valid_nodes_cells(valid_nodes, valid_cells)

  ! fill in the arrays
  c_ghosts = 1
  idx=1
  do j=self%jsc,self%jec+1
    do i=self%isc, self%iec+1
      if (.not. valid_nodes(i,j)) cycle

      c_lon(idx) = self%lon(i,j)
      c_lat(idx) = self%lat(i,j)
      if(j .le. self%jec .and. i .le. self%iec) then
        c_ghosts(idx) = 0
      end if

      c_global_idx(idx) = global_idx(i,j)
      c_remote_idx(idx) = local_idx(i,j)
      c_partition(idx) = partition(i,j)

      idx = idx + 1
    end do
  end do

  if (c_nodes /= idx-1) then
    ! TODO do a proper assert / error
    stop 42
  end if
  
end subroutine
! ------------------------------------------------------------------------------
end module soca_geom_mod_c
