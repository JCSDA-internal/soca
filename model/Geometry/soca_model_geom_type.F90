module soca_model_geom_type

  use MOM_grid,                  only : ocean_grid_type
  use MOM_verticalGrid,          only : verticalGrid_type
  use kinds

  implicit none
  private

  type, public :: soca_model_geom
     type(ocean_grid_type) :: G
     type(VerticalGrid_type), pointer :: GV          
     integer :: nx
     integer :: ny
     integer :: nz
     integer :: ncat
     real(kind=kind_real), pointer :: lon(:,:)       !< The horizontal grid type     !< 2D array of longitude 
     real(kind=kind_real), pointer :: lat(:,:)       !< 2D array of latitude
     real(kind=kind_real), pointer :: z(:)           !<      
     real(kind=kind_real), pointer :: mask2d(:,:)    !< 0 = land 1 = ocean surface mask only
     real(kind=kind_real), pointer :: cell_area(:,:) !<
   contains
     procedure :: clone => geom_clone
     procedure :: print => geom_print
     procedure :: infotofile => geom_infotofile     
  end type soca_model_geom

contains

    subroutine geom_clone(self, other)
    implicit none
    class(soca_model_geom), intent(in)  :: self
    class(soca_model_geom), intent(out) :: other

    other%nx = self%nx
    other%ny = self%ny    
    other%nz = self%nz
    other%ncat = self%ncat
    other%lon => self%lon
    other%lat => self%lat
    other%z => self%z
    other%mask2d => self%mask2d
    other%cell_area => self%cell_area

  end subroutine geom_clone

  subroutine geom_print(self)

    implicit none

    class(soca_model_geom), intent(in) :: self

    print *, 'nx=', self%nx
    print *, 'ny=', self%ny    

  end subroutine geom_print
  
  subroutine geom_infotofile(self)

    use mpp_mod,                   only : mpp_pe, mpp_npes, mpp_root_pe, mpp_sync
    use fms_mod,         only : get_mosaic_tile_grid, write_data, set_domain
    use fms_io_mod,                only : fms_io_init, fms_io_exit
    
    implicit none
    class(soca_model_geom), intent(in) :: self
    integer :: unit
    character(len=256) :: geom_output_file = "geom_output.nc"
    character(len=256) :: geom_field_name  = "none"

    !unit = 20 + mpp_pe()
    !write(unit,*)'pe=', mpp_pe(), mpp_npes()
    !write(unit,*)'nx, ny = ', self%nx, self%ny
    !write(unit,*)'lon(:,1) = ', self%lon(:,1)

    call fms_io_init()
    call write_data( geom_output_file, "lon", self%lon, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "lat", self%lat, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "z", self%z, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "area", self%cell_area, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "mask2d", self%mask2d, self%G%Domain%mpp_domain)        
    call fms_io_exit()    
    
  end subroutine geom_infotofile

end module soca_model_geom_type
