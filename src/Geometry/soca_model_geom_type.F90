!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_model_geom_type

  use MOM_grid,                  only : ocean_grid_type
  use MOM_verticalGrid,          only : verticalGrid_type
  use ice_grid,                  only : ice_grid_type  
  use kinds

  implicit none
  private
  public :: geom_infotofile, geom_get_domain_indices
  
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
     real(kind=kind_real), allocatable, dimension(:,:) :: mask2d    !< 0 = land 1 = ocean
     real(kind=kind_real), allocatable, dimension(:,:) :: obsmask   !< 0 = land and halo 1 = ocean     
     real(kind=kind_real), allocatable, dimension(:,:) :: cell_area !<
     real(kind=kind_real), allocatable, dimension(:,:) :: rossby_radius !<     
   contains
     procedure :: init => geom_init
     procedure :: end => geom_end     
     procedure :: shortcuts => geom_associate     
     procedure :: clone => geom_clone
     procedure :: print => geom_print
     procedure :: get_rossby_radius => geom_rossby_radius
     procedure :: infotofile => geom_infotofile
  end type soca_model_geom

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

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
    if (allocated(self%obsmask)) deallocate(self%obsmask)    
    if (allocated(self%cell_area)) deallocate(self%cell_area)
    if (allocated(self%rossby_radius)) deallocate(self%rossby_radius)    
    
  end subroutine geom_end

  ! ------------------------------------------------------------------------------
  
  subroutine geom_clone(self, other)

    implicit none

    class(soca_model_geom), intent(in)  :: self
    class(soca_model_geom), intent(out) :: other

    other%G = self%G
    other%IG = self%IG
    call geom_associate(other)    
    
  end subroutine geom_clone

  ! ------------------------------------------------------------------------------
  
  subroutine geom_associate(self)

    implicit none

    class(soca_model_geom), intent(inout)  :: self
    integer                       :: nxny(2), nx, ny

    integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nzo, nzi, nzs
    integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

    ! Allocate arrays on data domain
    ! Note: Compute domain excludes halo (is, ie, js, je)
    !       Data domain includes halo (isd, ied, jsd, jed)
    ! Compute
    is   = self%G%isc  ; ie   = self%G%iec  ; js   = self%G%jsc  ; je   = self%G%jec ; nzo = self%G%ke

    ! Data
    isd  = self%G%isd  ; ied  = self%G%ied  ; jsd  = self%G%jsd  ; jed  = self%G%jed

    !Isq  = self%G%IscB ; Ieq  = self%G%IecB ; Jsq  = self%G%JscB ; Jeq  = self%G%JecB    
    !IsdB = self%G%IsdB ; IedB = self%G%IedB ; JsdB = self%G%JsdB ; JedB = self%G%JedB
   
    nxny = shape( self%G%GeoLonT )
    nx = nxny(1)
    ny = nxny(2)

    ! Extract geometry of interest from model's data structure.
    ! Common to ocean & sea-ice
    self%nx = nx
    self%ny = ny

    ! Can't point to data structure, so allocating ...
    allocate(self%lon(isd:ied,jsd:jed))
    allocate(self%lat(isd:ied,jsd:jed))    
    allocate(self%mask2d(isd:ied,jsd:jed))
    allocate(self%obsmask(isd:ied,jsd:jed))    
    allocate(self%cell_area(isd:ied,jsd:jed))
    allocate(self%rossby_radius(isd:ied,jsd:jed))
    
    self%lon = self%G%GeoLonT
    self%lat = self%G%GeoLatT
    self%mask2d = self%G%mask2dT
    self%cell_area = self%G%areaT

    ! Setting up mask used to qc out observation that are on land or out
    ! of the compute domain
    self%obsmask = self%G%mask2dT        

!!$    self%obsmask(isd:is-1,:)=0.0
!!$    self%obsmask(ie+1:,:)=0.0
!!$    self%obsmask(:,jsd:js-1)=0.0
!!$    self%obsmask(:,je+1:)=0.0    
    
    self%obsmask(isd:is,:)=0.0
    self%obsmask(ie:,:)=0.0
    self%obsmask(:,jsd:js)=0.0
    self%obsmask(:,je:)=0.0    

    ! Ocean
    self%nzo = self%G%ke
    !self%z => self%GV%sLayer

    ! Sea-ice
    self%ncat = self%IG%CatIce
    self%nzi = self%IG%NkIce
    self%nzs = self%IG%NkSnow    
    
  end subroutine geom_associate

  ! ------------------------------------------------------------------------------

  subroutine geom_print(self)

    implicit none

    class(soca_model_geom), intent(in) :: self

    print *, 'nx=', self%nx
    print *, 'ny=', self%ny    

  end subroutine geom_print

  ! ------------------------------------------------------------------------------

  subroutine geom_rossby_radius(self)
    use kinds
    use type_kdtree, only: kdtree_type
    use type_mpl    
    use tools_const, only: pi,req,deg2rad,rad2deg
    use fms_mod,         only : get_mosaic_tile_grid, write_data, set_domain, read_data
    use fms_io_mod,      only : fms_io_init, fms_io_exit
    use mpi
    use fckit_mpi_module, only: fckit_mpi_comm
    
    implicit none

    class(soca_model_geom), intent(inout) :: self

    integer :: unit, i, j, n
    real(kind=kind_real), allocatable :: lon(:),lat(:),rr(:)
    logical, allocatable :: mask(:)    
    type(kdtree_type) :: kdtree
    type(mpl_type) :: mpl    
    real(kind=kind_real) :: dum, dist(1),lonm(1),latm(1)
    integer :: isc, iec, jsc, jec
    integer :: index(1), nn, io
    character(len=256) :: geom_output_file = "geom_output.nc"
    type(fckit_mpi_comm) :: f_comm

    f_comm = fckit_mpi_comm()
        
    unit = 20
    open(unit=unit,file="rossrad.dat",status="old",action="read")
    n = 0
    do
       read(unit,*,iostat=io)
       if (io/=0) exit
       n = n+1
    end do
    rewind(unit)
    allocate(lon(n),lat(n),rr(n),mask(n))
    do i = 1, n
       read(unit,*) lat(i),lon(i),dum,rr(i)
    end do
    close(unit)

    !--- Initialize kd-tree
    mask=.true.
    where (lon>180.0)
       lon=lon-360.0
    end where
    lon=deg2rad*reshape(lon,(/n/))
    lat=deg2rad*reshape(lat,(/n/))

    call mpl%init(f_comm%communicator())
    call kdtree%create(mpl, n,lon,lat,mask)

    !--- Find nearest neighbor    
    isc   = self%G%isc
    iec   = self%G%iec
    jsc   = self%G%jsc
    jec   = self%G%jec

    nn=1 ! Num neighbors
    do i = isc, iec
       do j = jsc, jec
          lonm=self%lon(i,j)
          if (lonm(1)>180.0) lonm=lonm-360.0
          lonm=deg2rad*lonm
          latm(1)=deg2rad*self%lat(i,j)          
          call kdtree%find_nearest_neighbors(lonm(1),&
                                            &latm(1),&
                                            &nn,index,dist)
          self%rossby_radius(i,j)=rr(index(1))*1e3
       end do
    end do

    call fms_io_init()
    call write_data( geom_output_file, "rossby_radius", self%rossby_radius*self%mask2d, self%G%Domain%mpp_domain)    
    !call write_data( geom_output_file, "mask2d", self%mask2d, self%G%Domain%mpp_domain)
!!$    print *,'rosbby:',self%rossby_radius
!!$    print *,'============================================'
!!$    read(*,*)    
!!$    call read_data(geom_output_file,"rossby_radius",self%rossby_radius,domain=self%G%Domain%mpp_domain)
!!$    print *,'rosbby:',self%rossby_radius
!!$    read(*,*)
    call fms_io_exit()
    
  end subroutine geom_rossby_radius

  ! ------------------------------------------------------------------------------
  
  subroutine geom_infotofile(self)

    use mpp_mod,          only : mpp_pe, mpp_npes, mpp_root_pe, mpp_sync
    use fms_mod,         only : get_mosaic_tile_grid, write_data, set_domain
    use fms_io_mod,      only : fms_io_init, fms_io_exit
    use soca_utils
    
    implicit none
    class(soca_model_geom), intent(in) :: self
    character(len=256) :: geom_output_file = "geom_output.nc"
    character(len=256) :: geom_output_pe,varname
    integer :: pe
    character(len=8) :: fmt = '(I5.5)'
    character(len=1024) :: strpe
    integer :: isc,iec,jsc,jec 
    
    call fms_io_init()
    ! Save full domain
    call write_data( geom_output_file, "lon", self%lon, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "lat", self%lat, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "area", self%cell_area, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "rossby_radius", self%rossby_radius, self%G%Domain%mpp_domain)
    call write_data( geom_output_file, "mask2d", self%mask2d, self%G%Domain%mpp_domain)
    call fms_io_exit()    

    ! Save local compute grid
    ! Get compute domain
    isc = self%G%isc    
    iec = self%G%iec
    jsc = self%G%jsc    
    jec = self%G%jec    

    pe = mpp_pe()
    write (strpe,fmt) pe
    geom_output_pe='geom_output_'//trim(strpe)//'.nc'
    varname='obsmask'
    call write2pe(reshape(self%obsmask,(/self%nx*self%ny/)),varname,geom_output_pe,.false.)
    varname='lon'
    call write2pe(reshape(self%lon,(/self%nx*self%ny/)),varname,geom_output_pe,.true.)
    varname='lat'
    call write2pe(reshape(self%lat,(/self%nx*self%ny/)),varname,geom_output_pe,.true.)        

    print *,'--------------------- ',pe,mpp_npes(),strpe,geom_output_pe

  end subroutine geom_infotofile

  ! ------------------------------------------------------------------------------
  
  subroutine geom_get_domain_indices(self, domain_type, is, ie, js, je)

    implicit none

    class(soca_model_geom), intent(in)  :: self
    character(7),            intent(in) :: domain_type
    integer,                intent(out) :: is, ie, js, je

    select case (trim(domain_type))
       case ("compute")
          is = self%G%isc
          ie = self%G%iec
          js = self%G%jsc
          je = self%G%jec

       case ("data")
          is = self%G%isd
          ie = self%G%ied
          js = self%G%jsd
          je = self%G%jed

       end select

  end subroutine geom_get_domain_indices
  
end module soca_model_geom_type
