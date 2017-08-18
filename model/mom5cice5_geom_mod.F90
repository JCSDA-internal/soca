
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
  !>
  !>   level 0          ---- Surface                     Tsfc
  !>   level 1          ---- Upper snow level            Qsno
  !>     .                        .                       .
  !>     .                        .                       .
  !>     .                        .                       .        
  !>   level nzs        ---- Lower snow level            Qsno
  !>   level nzs+1      ---- Upper ice level             Qice     
  !>     .                        .                       .
  !>     .                        .                       .
  !>     .                        .                       .
  !>   level nzs+nzi    ---- Upper ice level             Qice
  !>   level nzs+nzi+1  ---- Ice/Ocean interface level   Tfreeze (from S)
  !>   level nzs+nzi+2  ---- Upper level Ocean           SST    
  !>
  
  type :: mom5cice5_geom
     integer :: nx
     integer :: ny
     integer :: nzs                                      !< Number of snow levels (!!!! Fields hard coded for 1 level !!!!)
     integer :: nzi                                      !< Number of ice levels     
     integer :: nzo                                      !< Number of ocean levels
     integer :: ncat                                     !< Number of ice thickness categories
     character(len=128) :: gridfname                     !< Name of file containing the grid specs
     character(len=128) :: icemaskfname                  !< Name of file containing the grid specs
     real(kind=kind_real), allocatable :: lon(:,:)       !< 2D array of longitude 
     real(kind=kind_real), allocatable :: lat(:,:)       !< 2D array of latitude
     integer,              allocatable :: level(:)       !< 1D array of levels. See beginning of script for def.
     real(kind=kind_real), allocatable :: mask(:,:)      !< 0 = land 1 = ocean surface mask only
     real(kind=kind_real), allocatable :: icemask(:,:)   !< 0 = land/liquid ocean 1 = some ice
     real(kind=kind_real), allocatable :: cell_area(:,:) !< 
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
    integer :: jj
    integer :: start2(2), count2(2)
    
    call mom5cice5_geom_registry%init()
    call mom5cice5_geom_registry%add(c_key_self)
    call mom5cice5_geom_registry%get(c_key_self,self)

    self%nx = config_get_int(c_conf, "nx")
    self%ny = config_get_int(c_conf, "ny")
    self%nzo = config_get_int(c_conf, "nzo")      ! Number of ocean levels
    self%nzi = config_get_int(c_conf, "nzi")      ! Number of sea-ice levels
    self%nzs = config_get_int(c_conf, "nzs")      ! Number of snow levels
    self%ncat = config_get_int(c_conf, "ncat")
    self%gridfname = config_get_string(c_conf, len(self%gridfname), "gridfname")
    self%icemaskfname = config_get_string(c_conf, len(self%icemaskfname), "icemaskfname")

    allocate(self%lon(self%nx,self%ny))
    allocate(self%lat(self%nx,self%ny))
    allocate(self%level(self%nzs+self%nzi+self%nzo+2))    
    allocate(self%mask(self%nx,self%ny))
    allocate(self%icemask(self%nx,self%ny))
    allocate(self%cell_area(self%nx,self%ny))

    !do jj = 0,size(self%level, 1)
    do jj = 1,size(self%level, 1)       
       self%level(jj) = jj
    end do
    start2 = (/1,1/)
    count2 = (/self%nx,self%ny/)    
    print *,'ocean grid file name:',self%gridfname
    print *,'ice mask file name:',self%icemaskfname
    varname='x_T'; call ncread_fld(self%gridfname, varname, self%lon, start2, count2)
    varname='y_T'; call ncread_fld(self%gridfname, varname, self%lat, start2, count2)
    varname='wet'; call ncread_fld(self%gridfname, varname, self%mask, start2, count2)
    varname='area_T'; call ncread_fld(self%gridfname, varname, self%cell_area, start2, count2)
    varname='iceumask'; call ncread_fld(self%icemaskfname, varname, self%icemask, start2, count2)

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
    other%nzs = self%nzs
    other%ncat = self%ncat
    other%gridfname = self%gridfname
    other%lon = self%lon
    other%lat = self%lat
    other%level = self%level
    other%mask = self%mask
    other%icemask = self%icemask
    other%cell_area = self%cell_area

  end subroutine c_mom5cice5_geo_clone

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_geo_delete(c_key_self) bind(c,name='mom5cice5_geo_delete_f90')

    implicit none
    integer(c_int), intent(inout) :: c_key_self     
    type(mom5cice5_geom), pointer :: self

    call mom5cice5_geom_registry%get(c_key_self , self )
    if (allocated(self%lon)) deallocate(self%lon)
    if (allocated(self%lat)) deallocate(self%lat)
    if (allocated(self%level)) deallocate(self%level)    
    if (allocated(self%mask)) deallocate(self%mask)
    if (allocated(self%icemask)) deallocate(self%icemask)
    if (allocated(self%cell_area)) deallocate(self%cell_area)
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

  subroutine c_mom5cice5_geo_getgeofld(c_key_self, geoloc, geofld, level) bind(c,name='mom5cice5_geo_getgeofld_f90')

    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: geofld     !< 0: lon, 1: lat, 2: mask, 3: cell_area
    integer(c_int), intent(in) :: level    
    type(mom5cice5_geom), pointer :: self
    real(kind=kind_real), allocatable :: geoloc(:) 
    integer :: jx, jy, jj

    call mom5cice5_geom_registry%get(c_key_self , self )

    ! Extract geometry elements from mom5cice5 geom object
    allocate(geoloc(self%nx*self%ny))
    if (geofld==0) geoloc=reshape(self%lon, (/self%nx*self%ny/))       
    if (geofld==1) geoloc=reshape(self%lat, (/self%nx*self%ny/))
    if (geofld==2) then
       if (level==self%nzs+self%nzi+2) then ! Last level is the ocean mask
          geoloc=reshape(self%mask, (/self%nx*self%ny/))
       else
          geoloc=reshape(self%icemask, (/self%nx*self%ny/))
       end if
    end if
    if (geofld==3) geoloc=reshape(self%cell_area, (/self%nx*self%ny/))
  end subroutine c_mom5cice5_geo_getgeofld

end module mom5cice5_geom_mod
