!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!> Fortran module handling geometry for MOM6

module soca_geom_mod

  use iso_c_binding
  use config_mod
  use kinds
  use soca_model_geom_type
  
  implicit none
  private
  public :: soca_geom_registry, soca_geom

  type :: soca_geom
     type( soca_model_geom ) :: ocean ! Includes sea-ice
  end type soca_geom
  
#define LISTED_TYPE soca_geom

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_geom_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

  ! ------------------------------------------------------------------------------
  subroutine c_soca_geo_setup(c_key_self, c_conf) bind(c,name='soca_geo_setup_f90')
    use netcdf
    use soca_mom6

    implicit none

    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr),       intent(in) :: c_conf
    type(soca_geom),      pointer :: self
    
    call soca_geom_registry%init()
    call soca_geom_registry%add(c_key_self)
    call soca_geom_registry%get(c_key_self,self)

    call self%ocean%init(c_conf)
    call self%ocean%get_rossby_radius()
    call self%ocean%validindex() !BUG: Needs a halo of 2 to work
    call self%ocean%infotofile()
    
  end subroutine c_soca_geo_setup

  ! ------------------------------------------------------------------------------

  subroutine c_soca_geo_clone(c_key_self, c_key_other) bind(c,name='soca_geo_clone_f90')
    use soca_mom6
    implicit none
    integer(c_int), intent(in   ) :: c_key_self
    integer(c_int), intent(inout) :: c_key_other

    type(soca_geom), pointer :: self, other

    call soca_geom_registry%add(c_key_other)
    call soca_geom_registry%get(c_key_other, other)
    call soca_geom_registry%get(c_key_self , self )

    call self%ocean%clone(other%ocean)
    
  end subroutine c_soca_geo_clone

  ! ------------------------------------------------------------------------------

  subroutine c_soca_geo_delete(c_key_self) bind(c,name='soca_geo_delete_f90')
    use soca_mom6

    use fms_io_mod,      only: fms_io_init, fms_io_exit    
    implicit none
    integer(c_int), intent(inout) :: c_key_self     
    type(soca_geom), pointer :: self

    call soca_geom_registry%get(c_key_self, self)
    call soca_geom_end(self%ocean%G, self%ocean%GV)
    call self%ocean%end()
    call soca_geom_registry%remove(c_key_self)

  end subroutine c_soca_geo_delete

  ! ------------------------------------------------------------------------------

  subroutine c_soca_geo_info(c_key_self) bind(c,name='soca_geo_info_f90')

    implicit none
    integer(c_int), intent(in   ) :: c_key_self
    type(soca_geom), pointer :: self

    call soca_geom_registry%get(c_key_self , self)
    call self%ocean%print()
    call self%ocean%infotofile()    
    
  end subroutine c_soca_geo_info

  ! ------------------------------------------------------------------------------

end module soca_geom_mod
