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

  type :: soca_model_geom
     integer :: nx
     integer :: ny
     integer :: nz
     integer :: ncat
     real(kind=kind_real), pointer :: lon(:,:)       !< 2D array of longitude 
     real(kind=kind_real), pointer :: lat(:,:)       !< 2D array of latitude
     real(kind=kind_real), pointer :: z(:)         !<      
     real(kind=kind_real), pointer :: mask2d(:,:)      !< 0 = land 1 = ocean surface mask only
     real(kind=kind_real), pointer :: cell_area(:,:) !<
  end type soca_model_geom

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
    use MOM_grid,                  only : MOM_grid_init, MOM_grid_end, ocean_grid_type
    use MOM_file_parser,    only: get_param, log_version, param_file_type
    use MOM_get_input,      only: Get_MOM_Input, directories
    use MOM_domains,              only : MOM_domains_init, clone_MOM_domain
    use MOM_io,                   only : MOM_io_init
    use MOM_hor_index,             only : hor_index_type, hor_index_init    
    !use constants_mod,           only: constants_init    
    !use ocean_model_mod,         only: ocean_public_type, ocean_state_type, ice_ocean_boundary_type
    !use ice_model_mod,           only: ice_data_type
    
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr), intent(in)    :: c_conf
    type(soca_geom), pointer :: self
    integer :: varid, ncid, nxdimid, nydimid
    character(len=128)  :: varname
    integer :: jj, nx0, ny0
    integer :: start2(2), count2(2), pe

    type(ocean_grid_type) :: G          !< The horizontal grid type
    type(param_file_type) :: param_file           !< A structure to parse for run-time parameters
    type(directories)     :: dirs_tmp             !< A structure containing several relevant directory paths
    type(hor_index_type)            :: HI  !  A hor_index_type for array extents
    logical :: global_indexing   ! If true use global horizontal index values instead
    
    call soca_geom_registry%init()
    call soca_geom_registry%add(c_key_self)
    call soca_geom_registry%get(c_key_self,self)

    !call get_MOM_Input(param_file, dirs_tmp)!, check_params=.false.)
    !call MOM_domains_init(G%domain, param_file)!, symmetric=symmetric)
    !!call MOM_debugging_init(param_file)
    !!call diag_mediator_infrastructure_init()
    !call MOM_io_init(param_file)

    !call hor_index_init(G%Domain, HI, param_file, &
     !                 local_indexing=.not.global_indexing)
    
    !print *,'calling MOM_grid_init'
    !call MOM_grid_init(G, param_file)
    !print *,'=============isc=',G%isc,G%isd_global,maxval(G%mask2dT),G%geoLonT
    !call soca_simple_init()

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

!!$    other%nx = self%nx
!!$    other%ny = self%ny
!!$    other%nzo = self%nzo
!!$    other%nzi = self%nzi
!!$    other%nzs = self%nzs
!!$    other%ncat = self%ncat
!!$    other%gridfname = self%gridfname
!!$    other%lon = self%lon
!!$    other%lat = self%lat
!!$    other%level = self%level
!!$    other%mask = self%mask
!!$    other%icemask = self%icemask
!!$    other%cell_area = self%cell_area

  end subroutine c_soca_geo_clone

  ! ------------------------------------------------------------------------------

  subroutine c_soca_geo_delete(c_key_self) bind(c,name='soca_geo_delete_f90')
    !use soca_mom6sis2
    !use mpp_mod,         only: mpp_exit    
    use fms_io_mod,      only: fms_io_init, fms_io_exit    
    implicit none
    integer(c_int), intent(inout) :: c_key_self     
    type(soca_geom), pointer :: self

    call soca_geom_registry%get(c_key_self, self)
    !if (allocated(self%lon)) deallocate(self%lon)
    !if (allocated(self%lat)) deallocate(self%lat)
    !if (allocated(self%level)) deallocate(self%level)    
    !if (allocated(self%mask)) deallocate(self%mask)
    !if (allocated(self%icemask)) deallocate(self%icemask)
    !if (allocated(self%cell_area)) deallocate(self%cell_area)
    call soca_geom_registry%remove(c_key_self)


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
!!$    c_nx = self%nx
!!$    c_ny = self%ny
!!$    c_nzo = self%nzo
!!$    c_nzi = self%nzi
!!$    c_ncat = self%ncat

  end subroutine c_soca_geo_info

  ! ------------------------------------------------------------------------------

end module soca_geom_mod
