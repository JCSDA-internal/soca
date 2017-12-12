!> Fortran module handling geometry for MOM5 & CICE5 model.
!! @ Todo: Remove sea-ice mask (maybe?), add nx0, ny0 to .json
module soca_geom_mod

  use iso_c_binding
  use config_mod
  use kinds
  !use netcdf
  use ncutils

  !use MOM_grid
  
  implicit none
  private
  public :: soca_geom
  public :: soca_geom_registry

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
  
  type :: soca_geom
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
     real(kind=kind_real), allocatable :: icemask(:,:)   !< 0 = land/liquid ocean; 1 = some ice
     real(kind=kind_real), allocatable :: cell_area(:,:) !<
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
    use constants_mod,           only: constants_init    
    use ocean_model_mod,         only: ocean_public_type, ocean_state_type, ice_ocean_boundary_type
    use ice_model_mod,           only: ice_data_type
    
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr), intent(in)    :: c_conf
    type(soca_geom), pointer :: self
    integer :: varid, ncid, nxdimid, nydimid
    character(len=128)  :: varname
    integer :: jj, nx0, ny0
    integer :: start2(2), count2(2), pe
    type (ocean_public_type) :: Ocean_
    type (ocean_state_type), pointer :: Ocean_state_ => NULL()
    type   (ice_data_type) :: Ice_
    
    call soca_geom_registry%init()
    call soca_geom_registry%add(c_key_self)
    call soca_geom_registry%get(c_key_self,self)

    call constants_init
    call soca_models_init(Ocean_, Ocean_state_, Ice_) ! switch ocean_state_type to public in
                                                  ! ocean_model_MOM.F90

    print *,'Ocean_state_%grid%geoLonT=',Ocean_state_%grid%geoLonT(1,1)
    print *,'Ocean_state_%grid%geoLonT=',Ice_%part_size(1,1,1)
    
    
    print *,'============================================================='
    print *,'============================================================='
    print *,'============================================================='
    print *,'===================out of coupler init======================='
    print *,'============================================================='
    print *,'============================================================='
    print *,'============================================================='
    
    !call ocean_model_restart(Ocean_state_)!, timestamp)
    !call ice_model_restart(Ice_)

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
    nx0=1 !20
    ny0=1 !60
    start2 = (/nx0,ny0/)
    count2 = (/self%nx,self%ny/)

    varname='x_T'; call ncread_fld(self%gridfname, varname, self%lon, start2, count2)
    varname='y_T'; call ncread_fld(self%gridfname, varname, self%lat, start2, count2)
    varname='wet'; call ncread_fld(self%gridfname, varname, self%mask, start2, count2)
    varname='area_T'; call ncread_fld(self%gridfname, varname, self%cell_area, start2, count2)
    varname='iceumask'; call ncread_fld(self%icemaskfname, varname, self%icemask, start2, count2)

  end subroutine c_soca_geo_setup

  ! ------------------------------------------------------------------------------

  subroutine c_soca_geo_clone(c_key_self, c_key_other) bind(c,name='soca_geo_clone_f90')

    implicit none
    integer(c_int), intent(in   ) :: c_key_self
    integer(c_int), intent(inout) :: c_key_other

    type(soca_geom), pointer :: self, other

    call soca_geom_registry%add(c_key_other)
    call soca_geom_registry%get(c_key_other, other)
    call soca_geom_registry%get(c_key_self , self )
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

  end subroutine c_soca_geo_clone

  ! ------------------------------------------------------------------------------

  subroutine c_soca_geo_delete(c_key_self) bind(c,name='soca_geo_delete_f90')
    use soca_mom6sis2
    use mpp_mod,         only: mpp_exit    
    use fms_io_mod,      only: fms_io_init, fms_io_exit    
    implicit none
    integer(c_int), intent(inout) :: c_key_self     
    type(soca_geom), pointer :: self

    !call fms_io_exit
    print *,'================ THE END ====================='


    call soca_models_end()

    call soca_geom_registry%get(c_key_self, self)
    if (allocated(self%lon)) deallocate(self%lon)
    if (allocated(self%lat)) deallocate(self%lat)
    if (allocated(self%level)) deallocate(self%level)    
    if (allocated(self%mask)) deallocate(self%mask)
    if (allocated(self%icemask)) deallocate(self%icemask)
    if (allocated(self%cell_area)) deallocate(self%cell_area)
    call soca_geom_registry%remove(c_key_self)


!call mpp_domains_exit
!call mpp_exit

    
  end subroutine c_soca_geo_delete

  ! ------------------------------------------------------------------------------

  subroutine c_soca_geo_info(c_key_self, c_nx, c_ny, c_nzo, c_nzi, c_ncat) bind(c,name='soca_geo_info_f90')

    implicit none
    integer(c_int), intent(in   ) :: c_key_self
    integer(c_int), intent(inout) :: c_nx
    integer(c_int), intent(inout) :: c_ny
    integer(c_int), intent(inout) :: c_nzo
    integer(c_int), intent(inout) :: c_nzi
    integer(c_int), intent(inout) :: c_ncat
    type(soca_geom), pointer :: self

    call soca_geom_registry%get(c_key_self , self )
    c_nx = self%nx
    c_ny = self%ny
    c_nzo = self%nzo
    c_nzi = self%nzi
    c_ncat = self%ncat

  end subroutine c_soca_geo_info

  ! ------------------------------------------------------------------------------

end module soca_geom_mod
