! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_geom_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use soca_geom_mod, only: soca_geom

implicit none
private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup geometry object
subroutine c_soca_geo_setup(c_self, c_conf) bind(c,name='soca_geo_setup_f90')
  type(c_ptr), intent(inout) :: c_self
  type(c_ptr),    intent(in) :: c_conf

  type(soca_geom), pointer :: self

  allocate(self)
  c_self = c_loc(self)

  call self%init(fckit_configuration(c_conf))

end subroutine c_soca_geo_setup

! ------------------------------------------------------------------------------
!> Clone geometry object
subroutine c_soca_geo_clone(c_self, c_other) bind(c,name='soca_geo_clone_f90')
  type(c_ptr), intent(inout) :: c_self
  type(c_ptr),    intent(in) :: c_other

  type(soca_geom), pointer :: self, other

  allocate(self)
  c_self = c_loc(self)
  call c_f_pointer(c_other, other)

  call self%clone(other)

end subroutine c_soca_geo_clone

! ------------------------------------------------------------------------------
!> Generate grid
subroutine c_soca_geo_gridgen(c_self) bind(c,name='soca_geo_gridgen_f90')
  type(c_ptr), intent(in) :: c_self

  type(soca_geom), pointer :: self
  call c_f_pointer(c_self, self)

  call self%gridgen()

end subroutine c_soca_geo_gridgen

! ------------------------------------------------------------------------------
!> Geometry destructor
subroutine c_soca_geo_delete(c_self) bind(c,name='soca_geo_delete_f90')
  type(c_ptr), intent(inout) :: c_self

  type(soca_geom), pointer :: self

  call c_f_pointer(c_self, self)

  call self%end()

  deallocate(self)

end subroutine c_soca_geo_delete

! ------------------------------------------------------------------------------
!> return begin and end of local geometry
subroutine c_soca_geo_start_end(c_self, ist, iend, jst, jend) bind(c, name='soca_geo_start_end_f90')
  type(c_ptr),      intent(in):: c_self
  integer(c_int), intent(out) :: ist, iend, jst, jend

  type(soca_geom), pointer :: self

  call c_f_pointer(c_self, self)

  ist  = self%isc
  iend = self%iec
  jst  = self%jsc
  jend = self%jec

end subroutine c_soca_geo_start_end


! ------------------------------------------------------------------------------

end module soca_geom_mod_c
