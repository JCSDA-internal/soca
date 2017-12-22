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
  end type soca_model_geom

contains

  subroutine geom_clone(self, other)
    class(soca_model_geom), intent(in)  :: self
    class(soca_model_geom), intent(out) :: other

    other%nx = self%nx
    other%ny = self%ny    
    other%nz = self%nz
    other%ncat = self%ncat
    other%lon = self%lon
    other%lat = self%lat
    other%z = self%z
    other%mask2d = self%mask2d
    other%cell_area = self%cell_area

  end subroutine geom_clone

  subroutine geom_print(this)
    class(soca_model_geom), intent(in) :: this
    print *,'[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]'
    print *, 'nx=', this%nx
    print *, 'ny=', this%ny    
    print *,'[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]'    
  end subroutine geom_print

end module soca_model_geom_type
