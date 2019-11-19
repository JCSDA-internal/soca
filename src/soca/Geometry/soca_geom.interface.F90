! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_geom_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use soca_geom_mod, only: soca_geom

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
subroutine c_soca_geo_setup(c_key_self, c_conf) bind(c,name='soca_geo_setup_f90')

  integer(c_int), intent(inout) :: c_key_self
  type(c_ptr),       intent(in) :: c_conf

  type(soca_geom), pointer :: self

  call soca_geom_registry%init()
  call soca_geom_registry%add(c_key_self)
  call soca_geom_registry%get(c_key_self,self)

  call self%init(fckit_configuration(c_conf))

end subroutine c_soca_geo_setup

! ------------------------------------------------------------------------------
!> Clone geometry object
subroutine c_soca_geo_clone(c_key_self, c_key_other) bind(c,name='soca_geo_clone_f90')

  integer(c_int), intent(in   ) :: c_key_self
  integer(c_int), intent(inout) :: c_key_other

  type(soca_geom), pointer :: self, other

  call soca_geom_registry%add(c_key_other)
  call soca_geom_registry%get(c_key_other, other)
  call soca_geom_registry%get(c_key_self , self )

  call self%clone(other)

end subroutine c_soca_geo_clone

! ------------------------------------------------------------------------------
!> Generate grid
subroutine c_soca_geo_gridgen(c_key_self, c_conf) bind(c,name='soca_geo_gridgen_f90')

  integer(c_int), intent(inout) :: c_key_self
  type(c_ptr),       intent(in) :: c_conf

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
subroutine c_soca_geo_start_end(c_key_self, ist, iend, jst, jend) bind(c, name='soca_geo_start_end_f90')

  implicit none

  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: ist, iend, jst, jend 

  type(soca_geom), pointer :: self
  call soca_geom_registry%get(c_key_self, self)

  ist  = self%isc
  iend = self%iec
  jst  = self%jsc
  jend = self%jec

end subroutine c_soca_geo_start_end

! ------------------------------------------------------------------------------
!> return nx, ny, nzo for global grid
subroutine c_soca_geo_global_grid_size(c_key_self, nx, ny, nzo) bind(c, name='soca_geo_global_grid_size_f90')

  implicit none

  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: nx, ny, nzo

  type(soca_geom), pointer :: self
  call soca_geom_registry%get(c_key_self, self)

  nx   = self%nx
  ny   = self%ny
  nzo  = self%nzo

end subroutine c_soca_geo_global_grid_size 

! ------------------------------------------------------------------------------
!> Dump basic geometry info in file and std output
subroutine c_soca_geo_info(c_key_self) bind(c,name='soca_geo_info_f90')

  integer(c_int), intent(in   ) :: c_key_self

  type(soca_geom), pointer :: self

  call soca_geom_registry%get(c_key_self , self)
  call self%print()

end subroutine c_soca_geo_info

! ------------------------------------------------------------------------------

end module soca_geom_mod_c
