module soca_model_geom_type

  use MOM_grid,                  only : ocean_grid_type
  use MOM_verticalGrid,          only : verticalGrid_type
  use ice_grid,                  only : ice_grid_type  
  use kinds

  implicit none
  private

  type, public :: soca_model_geom
     type(ocean_grid_type)            :: G     !< Ocean/sea-ice horizontal grid
     type(VerticalGrid_type), pointer :: GV    !< Ocean vertical grid
     type(ocean_grid_type)            :: seaice_G     !< Ocean/sea-ice horizontal grid     
     type(ice_grid_type)              :: IG    !< Ice grid
     ! Short-cut variables and convenience pointers
     integer :: nx
     integer :: ny
     integer :: nzo
     integer :: nzi
     integer :: nzs     
     integer :: ncat
     real(kind=kind_real), allocatable, dimension(:,:) :: lon      !< The horizontal grid type     !< 2D array of longitude 
     real(kind=kind_real), allocatable, dimension(:,:) :: lat       !< 2D array of latitude
     real(kind=kind_real), allocatable, dimension(:)   :: z           !<      
     real(kind=kind_real), allocatable, dimension(:,:) :: mask2d    !< 0 = land 1 = ocean surface mask only
     real(kind=kind_real), allocatable, dimension(:,:) :: cell_area !<
   contains
     procedure :: init => geom_init
     procedure :: end => geom_end     
     procedure :: shortcuts => geom_associate     
     procedure :: clone => geom_clone
     procedure :: print => geom_print
     procedure :: infotofile => geom_infotofile
  end type soca_model_geom

contains

  subroutine geom_init(self)
    
    use kinds
    use soca_mom6sis2, only : soca_geom_init !, soca_ice_geom_init
    
    implicit none

    class(soca_model_geom),   intent(out)  :: self    

    call soca_geom_init(self%G, self%GV, self%IG)
    !call soca_ice_geom_init(self%seaice_G, self%IG)    
    call geom_associate(self)
    
  end subroutine geom_init

  subroutine geom_end(self)
    
    implicit none

    class(soca_model_geom),   intent(out)  :: self    

    if (allocated(self%lon)) deallocate(self%lon)
    if (allocated(self%lat)) deallocate(self%lat)
    if (allocated(self%mask2d)) deallocate(self%mask2d)
    if (allocated(self%cell_area)) deallocate(self%cell_area)    
    
  end subroutine geom_end
    
  subroutine geom_clone(self, other)

    implicit none

    class(soca_model_geom), intent(in)  :: self
    class(soca_model_geom), intent(out) :: other

    other%G = self%G
    other%IG = self%IG
    call geom_associate(other)    
    
  end subroutine geom_clone

  subroutine geom_associate(self)

    implicit none

    class(soca_model_geom), intent(inout)  :: self
    integer                       :: nxny(2), nx, ny
    
    nxny = shape( self%G%GeoLonT )
    nx = nxny(1)
    ny = nxny(2)

    ! Extract geometry of interest from model's data structure.
    ! Common to ocean & sea-ice
    self%nx = nx
    self%ny = ny

    ! Can't point to data structure, so allocating ...
    allocate(self%lon(nx, ny))
    allocate(self%lat(nx, ny))
    allocate(self%mask2d(nx, ny))
    allocate(self%cell_area(nx, ny))
    
    self%lon = self%G%GeoLonT
    self%lat = self%G%GeoLatT
    self%mask2d = self%G%mask2dT
    self%cell_area = self%G%areaT
    
    ! Ocean
    self%nzo = self%G%ke
    !self%z => self%GV%sLayer

    ! Sea-ice
    self%ncat = self%IG%CatIce
    self%nzi = self%IG%NkIce
    self%nzs = self%IG%NkSnow    
    
  end subroutine geom_associate

  subroutine geom_print(self)

    implicit none

    class(soca_model_geom), intent(in) :: self

    print *, 'nx=', self%nx
    print *, 'ny=', self%ny    

  end subroutine geom_print
  
  subroutine geom_infotofile(self)

    use mpp_mod,          only : mpp_pe, mpp_npes, mpp_root_pe, mpp_sync
    use fms_mod,         only : get_mosaic_tile_grid, write_data, set_domain
    use fms_io_mod,      only : fms_io_init, fms_io_exit
    
    implicit none
    class(soca_model_geom), intent(in) :: self
    integer :: unit
    character(len=256) :: geom_output_file = "geom_output.nc"
    character(len=256) :: geom_field_name  = "none"

    call fms_io_init()
    !call set_domain( self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "lon", self%lon, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "lat", self%lat, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "z", self%z, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "area", self%cell_area, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "mask2d", self%mask2d, self%G%Domain%mpp_domain)        
    call fms_io_exit()    
    
  end subroutine geom_infotofile

end module soca_model_geom_type
