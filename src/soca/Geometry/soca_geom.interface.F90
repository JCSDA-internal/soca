! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

! ------------------------------------------------------------------------------
!> Setup geometry object
subroutine c_soca_geo_setup(c_key_self, c_conf) bind(c,name='soca_geo_setup_f90')
  use iso_c_binding
  use soca_geom_mod, only: soca_geom, soca_geom_registry
  implicit none
  integer(c_int), intent(inout) :: c_key_self
  type(c_ptr),       intent(in) :: c_conf

  type(soca_geom), pointer :: self

  call soca_geom_registry%init()
  call soca_geom_registry%add(c_key_self)
  call soca_geom_registry%get(c_key_self,self)

  call self%init(c_conf)

end subroutine c_soca_geo_setup

! ------------------------------------------------------------------------------
!> Clone geometry object
subroutine c_soca_geo_clone(c_key_self, c_key_other) bind(c,name='soca_geo_clone_f90')
  use iso_c_binding
  use soca_geom_mod, only: soca_geom, soca_geom_registry
  implicit none
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
  use iso_c_binding
  use soca_geom_mod, only: soca_geom, soca_geom_registry
  implicit none
  integer(c_int), intent(inout) :: c_key_self
  type(c_ptr),       intent(in) :: c_conf

  type(soca_geom), pointer :: self

  call soca_geom_registry%get(c_key_self,self)

  print *,'geom key = ',c_key_self
  call self%gridgen()

end subroutine c_soca_geo_gridgen

! ------------------------------------------------------------------------------
!> Geometry destructor
subroutine c_soca_geo_delete(c_key_self) bind(c,name='soca_geo_delete_f90')
  use iso_c_binding
  use soca_geom_mod, only: soca_geom, soca_geom_registry
  implicit none
  integer(c_int), intent(inout) :: c_key_self

  type(soca_geom), pointer :: self

  call soca_geom_registry%get(c_key_self, self)
  call self%end()
  call soca_geom_registry%remove(c_key_self)

end subroutine c_soca_geo_delete

! ------------------------------------------------------------------------------
!> Dump basic geometry info in file and std output
subroutine c_soca_geo_info(c_key_self) bind(c,name='soca_geo_info_f90')
  use iso_c_binding
  use soca_geom_mod, only: soca_geom, soca_geom_registry
  implicit none
  integer(c_int), intent(in   ) :: c_key_self

  type(soca_geom), pointer :: self

  call soca_geom_registry%get(c_key_self , self)
  call self%print()

end subroutine c_soca_geo_info

! ------------------------------------------------------------------------------
