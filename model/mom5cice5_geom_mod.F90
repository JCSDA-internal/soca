
!> Fortran module handling geometry for MOM5 & CICE5 model

module mom5cice5_geom_mod

  use iso_c_binding
  use config_mod
  use kinds
  !use netcdf
  !use ncutils

  implicit none
  private
  public :: mom5cice5_geom
  public :: mom5cice5_geom_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to hold geometry data for MOM5 & CICE5 model
  type :: mom5cice5_geom
     integer :: nx
     integer :: ny
     integer :: nzo
     integer :: nzi
     integer :: ncat
     character(len=128) :: filename
     real(kind=kind_real), allocatable :: lon(:,:)
     real(kind=kind_real), allocatable :: lat(:,:)     
  end type mom5cice5_geom

#define LISTED_TYPE mom5cice5_geom

  !> Linked list interface - defines registry_t type
#include "linkedList_i.f"

  !> Global registry
  type(registry_t) :: mom5cice5_geom_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_geo_setup(c_key_self, c_conf) bind(c,name='mom5cice5_geo_setup_f90')
    !use netcdf
    !use ncutils
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr), intent(in)    :: c_conf
    type(mom5cice5_geom), pointer :: self
    integer :: varid, ncid, nxdimid, nydimid

    call mom5cice5_geom_registry%init()
    call mom5cice5_geom_registry%add(c_key_self)
    call mom5cice5_geom_registry%get(c_key_self,self)

    self%nx = config_get_int(c_conf, "nx")
    self%ny = config_get_int(c_conf, "ny")
    self%nzo = config_get_int(c_conf, "nzo")
    self%nzi = config_get_int(c_conf, "nzi")
    self%ncat = config_get_int(c_conf, "ncat")
    self%filename = config_get_string(c_conf, len(self%filename), "filename")
    allocate(self%lon(self%nx,self%ny))
    allocate(self%lat(self%nx,self%ny))

    !call nccheck(nf90_open(self%filename, nf90_nowrite,ncid))
    !Get the size of the state
    !call nccheck(nf90_inq_dimid(ncid, 'grid_x_T', nxdimid))
    !call nccheck(nf90_inquire_dimension(ncid, nxdimid, len = self%nx))
    !call nccheck(nf90_inq_dimid(ncid, 'grid_y_T', nydimid))
    !call nccheck(nf90_inquire_dimension(ncid, nydimid, len = self%ny))

    !call nccheck(nf90_inq_varid(ncid, 'x_T', varid))
    !call nccheck(nf90_get_var(ncid, varid, self%lon))
    !call nccheck(nf90_inq_varid(ncid, 'y_T', varid))
    !call nccheck(nf90_get_var(ncid, varid, self%lat))

    !call check(nf90_close(ncid))

  end subroutine c_mom5cice5_geo_setup

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_geo_clone(c_key_self, c_key_other) bind(c,name='mom5cice5_geo_clone_f90')
    implicit none
    integer(c_int), intent(in   ) :: c_key_self
    integer(c_int), intent(inout) :: c_key_other

    type(mom5cice5_geom), pointer :: self, other

    call mom5cice5_geom_registry%add(c_key_other)
    call mom5cice5_geom_registry%get(c_key_other, other)
    call mom5cice5_geom_registry%get(c_key_self , self )
    other%nx = self%nx
    other%ny = self%ny
    other%nzo = self%nzo
    other%nzi = self%nzi
    other%ncat = self%ncat
    other%filename = self%filename
    other%lon = self%lon
    other%lat = self%lat

  end subroutine c_mom5cice5_geo_clone

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_geo_delete(c_key_self) bind(c,name='mom5cice5_geo_delete_f90')

    implicit none
    integer(c_int), intent(inout) :: c_key_self     
    type(mom5cice5_geom), pointer :: self

    call mom5cice5_geom_registry%get(c_key_self , self )
    if (allocated(self%lon)) deallocate(self%lon)
    if (allocated(self%lat)) deallocate(self%lat)
    call mom5cice5_geom_registry%remove(c_key_self)

  end subroutine c_mom5cice5_geo_delete

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_geo_info(c_key_self, c_nx, c_ny, c_nzo, c_nzi, c_ncat) bind(c,name='mom5cice5_geo_info_f90')
    implicit none
    integer(c_int), intent(in   ) :: c_key_self
    integer(c_int), intent(inout) :: c_nx
    integer(c_int), intent(inout) :: c_ny
    integer(c_int), intent(inout) :: c_nzo
    integer(c_int), intent(inout) :: c_nzi
    integer(c_int), intent(inout) :: c_ncat
    type(mom5cice5_geom), pointer :: self

    call mom5cice5_geom_registry%get(c_key_self , self )
    c_nx = self%nx
    c_ny = self%ny
    c_nzo = self%nzo
    c_nzi = self%nzi
    c_ncat = self%ncat

  end subroutine c_mom5cice5_geo_info

  ! ------------------------------------------------------------------------------

end module mom5cice5_geom_mod
