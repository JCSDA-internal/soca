!
! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_omb_stats_mod

  use kinds
  use netcdf
  use soca_utils
  use type_mpl
  use tools_const, only: deg2rad,rad2deg
  use type_kdtree, only: kdtree_type

  implicit none
  private

  type, public :: soca_omb_stats
     integer                           :: nlocs     
     real(kind=kind_real), allocatable :: lon(:)       
     real(kind=kind_real), allocatable :: lat(:)       
     real(kind=kind_real), allocatable :: bgerr(:)
     real(kind=kind_real), allocatable :: bgerr_model(:,:)     
   contains
     procedure :: init => soca_omb_stats_init
     procedure :: bin => soca_omb_stats_bin
     procedure :: exit => soca_omb_stats_exit
  end type soca_omb_stats

contains

  ! ------------------------------------------------------------------------------  
  subroutine soca_omb_stats_init(self)
    class(soca_omb_stats), intent(inout) :: self

    integer(kind=4) :: ncid
    integer(kind=4) :: dimid
    integer(kind=4) :: varid

    ! TODO: reading from master only, then broadcast
    call nc_check(nf90_open('godas_sst_bgerr.nc', nf90_nowrite,ncid))

    ! Get the size of the horizontal grid
    call nc_check(nf90_inq_dimid(ncid, 'nlocs', dimid))
    call nc_check(nf90_inquire_dimension(ncid, dimid, len = self%nlocs))

    allocate(self%lon(self%nlocs), self%lat(self%nlocs), self%bgerr(self%nlocs))

    ! Get longitude    
    call nc_check(nf90_inq_varid(ncid,'longitude',varid))
    call nc_check(nf90_get_var(ncid,varid,self%lon))

    ! Get latitude
    call nc_check(nf90_inq_varid(ncid,'latitude',varid))
    call nc_check(nf90_get_var(ncid,varid,self%lat))

    ! Get omb stats
    call nc_check(nf90_inq_varid(ncid,'sst_bgerr',varid))
    call nc_check(nf90_get_var(ncid,varid,self%bgerr))        

    ! Close netcdf file
    call nc_check(nf90_close(ncid))

    ! Rotate longitude
    where (self%lon>180.0_kind_real)
       self%lon=self%lon-360.0_kind_real
    end where

    
  end subroutine soca_omb_stats_init

  ! ------------------------------------------------------------------------------
  subroutine soca_omb_stats_bin(self, lon, lat)
    class(soca_omb_stats), intent(inout) :: self
    real(kind=kind_real),     intent(in) :: lon(:,:)       
    real(kind=kind_real),     intent(in) :: lat(:,:) 

    integer :: is, ie, js, je, i, j
    type(mpl_type) :: mpl
    type(kdtree_type) :: kdtree
    integer :: index(1), nn=1
    real(kind=kind_real) :: lonm(1), latm(1), dist(1)
    
    ! Get local array bounds
    is = lbound(lon,dim=1)
    ie = ubound(lon,dim=1)
    js = lbound(lon,dim=2)
    je = ubound(lon,dim=2)        

    allocate(self%bgerr_model(is:ie,js:je))
    self%bgerr_model = 0.5_kind_real ! Set to dummy value for now
                                     ! TODO: call to kdtree

    !--- Initialize kd-tree
    call mpl%init()
    call kdtree%alloc(mpl, self%nlocs)
    call kdtree%init(mpl, deg2rad*self%lon, deg2rad*self%lat)

    do i = is, ie
       do j = js, je
          lonm=lon(i,j)
          if (lonm(1)>180.0_kind_real) lonm=lonm-360.0_kind_real
          lonm=deg2rad*lonm
          latm(1)=deg2rad*lat(i,j)
          call kdtree%find_nearest_neighbors(mpl,&
                                            &lonm(1),&
                                            &latm(1),&
                                            &nn,index,dist)
          self%bgerr_model(i,j)=self%bgerr(index(1))
       end do
    end do
    
    ! Release memory
    call kdtree%dealloc()

    
  end subroutine soca_omb_stats_bin

  ! ------------------------------------------------------------------------------
  subroutine soca_omb_stats_exit(self)
    class(soca_omb_stats), intent(inout) :: self
    
    deallocate(self%lon, self%lat, self%bgerr, self%bgerr_model)
    
  end subroutine soca_omb_stats_exit

end module soca_omb_stats_mod
