
!> Fortran module handling geometry for MOM5 & CICE5 model

module mom5cice5_geom_mod

  use iso_c_binding
  use config_mod
  use kinds
  use netcdf
  use ncutils

  implicit none
  private
  public :: mom5cice5_geom
  public :: mom5cice5_geom_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to hold geometry data for MOM5 & CICE5 model
  type :: mom5cice5_geom
     integer :: nx
     integer :: ny
     integer :: nzo                                 !< Number of ocean levels
     integer :: nzi                                 !< Number of ice levels
     integer :: ncat                                !< Number of ice thickness categories
     character(len=128) :: gridfname                !< Name of file containing the grid specs
     real(kind=kind_real), allocatable :: lon(:,:)  !< 2D array of longitude 
     real(kind=kind_real), allocatable :: lat(:,:)  !< 2D array of latitude
     !real(kind=kind_real), allocatable :: zi(:,:,:) !< 3D array of depth bellow ice/snow surface    
     real(kind=kind_real), allocatable :: mask(:,:) !< 0 = land 1 = ocean NEED TO MERGE WITH ICE MASK
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
    use netcdf
    use interface_ncread_fld, only: ncread_fld

    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr), intent(in)    :: c_conf
    type(mom5cice5_geom), pointer :: self
    integer :: varid, ncid, nxdimid, nydimid
    character(len=128)  :: varname

    call mom5cice5_geom_registry%init()
    call mom5cice5_geom_registry%add(c_key_self)
    call mom5cice5_geom_registry%get(c_key_self,self)

    self%nx = config_get_int(c_conf, "nx")
    self%ny = config_get_int(c_conf, "ny")
    self%nzo = config_get_int(c_conf, "nzo")
    self%nzi = config_get_int(c_conf, "nzi")
    self%ncat = config_get_int(c_conf, "ncat")
    self%gridfname = config_get_string(c_conf, len(self%gridfname), "gridfname")

    allocate(self%lon(self%nx,self%ny))
    allocate(self%lat(self%nx,self%ny))
    allocate(self%mask(self%nx,self%ny))

    print *,self%gridfname
    varname='x_T'; call ncread_fld(self%gridfname, varname, self%lon, self%nx, self%ny)
    varname='y_T'; call ncread_fld(self%gridfname, varname, self%lat, self%nx, self%ny)
    varname='wet'; call ncread_fld(self%gridfname, varname, self%mask, self%nx, self%ny)

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
    other%gridfname = self%gridfname
    other%lon = self%lon
    other%lat = self%lat
    other%mask = self%mask

  end subroutine c_mom5cice5_geo_clone

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_geo_delete(c_key_self) bind(c,name='mom5cice5_geo_delete_f90')

    implicit none
    integer(c_int), intent(inout) :: c_key_self     
    type(mom5cice5_geom), pointer :: self

    call mom5cice5_geom_registry%get(c_key_self , self )
    if (allocated(self%lon)) deallocate(self%lon)
    if (allocated(self%lat)) deallocate(self%lat)
    if (allocated(self%mask)) deallocate(self%mask)
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

end module mom5cice5_geom_mod
