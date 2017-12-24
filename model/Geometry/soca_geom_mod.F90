!> Fortran module handling geometry for MOM5 & CICE5 model.
!! @ Todo: Remove sea-ice mask (maybe?), add nx0, ny0 to .json
module soca_geom_mod

  use iso_c_binding
  use config_mod
  use kinds
  !use soca_mom6sis2,             only : soca_model_geom
  use soca_model_geom_type
  
  implicit none
  private
  public :: soca_geom_registry, soca_geom

  type :: soca_geom
     type( soca_model_geom ) :: ocean
     type( soca_model_geom ) :: ice
     type( soca_model_geom ) :: snow
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
    use interface_ncread_fld, only: ncread_fld
    use soca_mom6sis2

    implicit none

    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr),       intent(in) :: c_conf
    type(soca_geom),      pointer :: self
    integer                       :: nxny(2), nx, ny
    
    call soca_geom_registry%init()
    call soca_geom_registry%add(c_key_self)
    call soca_geom_registry%get(c_key_self,self)

    call soca_geom_init(self%ocean%G, self%ocean%GV)
    
    nxny = shape( self%ocean%G%GeoLonT )
    nx = nxny(1)
    ny = nxny(2)    

    ! Ocean
    self%ocean%lon => self%ocean%G%GeoLonT
    self%ocean%lat => self%ocean%G%GeoLatT    
    self%ocean%mask2d => self%ocean%G%mask2dT
    self%ocean%cell_area => self%ocean%G%areaT
    self%ocean%z => self%ocean%GV%sLayer   ! pb here: Not initialized ...    
    self%ocean%nx = nx
    self%ocean%ny = ny    
    self%ocean%nz = self%ocean%G%ke
    self%ocean%ncat = 0    

    ! Assume Same horizontal grid for Ocean and Sea-ice 
    ! Sea-ice
    self%ice%lon => self%ocean%G%GeoLonT
    self%ice%lat => self%ocean%G%GeoLatT    
    self%ice%mask2d => self%ocean%G%mask2dT
    self%ice%cell_area => self%ocean%G%areaT
    self%ice%nx = nx
    self%ice%ny = ny
    self%ice%nz = 7      ! HARDCODED ... NEED TO FIX
    self%ice%ncat = 5    ! HARDCODED ... NEED TO FIX    
    
    ! Snow over Sea-ice
    self%snow%lon => self%ocean%G%GeoLonT
    self%snow%lat => self%ocean%G%GeoLatT    
    self%snow%mask2d => self%ocean%G%mask2dT
    self%snow%cell_area => self%ocean%G%areaT    
    self%snow%nx = nx
    self%snow%ny = ny
    self%snow%nz = 1      ! HARDCODED ... NEED TO FIX
    self%snow%ncat = self%ice%ncat    ! HARDCODED ... NEED TO FIX

    call self%ocean%infotofile()
    
  end subroutine c_soca_geo_setup

  ! ------------------------------------------------------------------------------

  subroutine c_soca_geo_clone(c_key_self, c_key_other) bind(c,name='soca_geo_clone_f90')
    use soca_mom6sis2
    implicit none
    integer(c_int), intent(in   ) :: c_key_self
    integer(c_int), intent(inout) :: c_key_other

    type(soca_geom), pointer :: self, other

    call soca_geom_registry%add(c_key_other)
    call soca_geom_registry%get(c_key_other, other)
    call soca_geom_registry%get(c_key_self , self )

    call self%ocean%clone(other%ocean)
    call self%ice%clone(other%ice)
    call self%snow%clone(other%snow)        

  end subroutine c_soca_geo_clone

  ! ------------------------------------------------------------------------------

  subroutine c_soca_geo_delete(c_key_self) bind(c,name='soca_geo_delete_f90')
    use soca_mom6sis2

    use fms_io_mod,      only: fms_io_init, fms_io_exit    
    implicit none
    integer(c_int), intent(inout) :: c_key_self     
    type(soca_geom), pointer :: self

    call soca_geom_registry%get(c_key_self, self)
    call soca_geom_end(self%ocean%G, self%ocean%GV)
    call soca_geom_registry%remove(c_key_self)

    !self => NULL()
    
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
