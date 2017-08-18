
!> Handle fields for the QG model

module mom5cice5_fields

  use config_mod
  use mom5cice5_geom_mod
  use mom5cice5_vars_mod
  use kinds

  implicit none
  private

  public :: mom5cice5_field, &
       & create, delete, zeros, dirac, random, copy, &
       & self_add, self_schur, self_sub, self_mul, axpy, &
       & dot_prod, add_incr, diff_incr, &
       & read_file, write_file, gpnorm, fldrms, &
       & change_resol, convert_to_ug, convert_from_ug
  public :: mom5cice5_field_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to hold QG fields
  type :: mom5cice5_field
     integer :: nx                     !< Zonal grid dimension
     integer :: ny                     !< Meridional grid dimension
     integer :: nzo                    !< Number of z levels in the ocean
     integer :: nzi                    !< Number of levels in sea-ice
     integer :: nzs                    !< Number of levels in snow       !!!!!!!!!!! HARD CODED TO 1 !!!!!!!!!!!!!!!
     integer :: ncat                   !< Number of sea-ice thickness categories    
     integer :: nf                     !< Number of fields
     character(len=128) :: gridfname   !< Grid file name
     character(len=128) :: cicefname   !< Fields file name for cice
     character(len=128) :: momfname    !< Fields file name for mom
     real(kind=kind_real), allocatable :: lat(:,:)            !< ... from geom
     real(kind=kind_real), allocatable :: lon(:,:)            !< ... from geom  !!!! NEED TO FIND BETTER SOLUTION !!!!!!
     integer,              allocatable :: level(:)            !< ... from geom
     real(kind=kind_real), allocatable :: icemask(:,:)        !< ... from geom Sea-ice mask
     real(kind=kind_real), allocatable :: mask(:,:)           !< ... from geom ocean mask !!!!! ONLY SURFACE !!!!!!        
     real(kind=kind_real), allocatable :: cicen(:,:,:)        !< Sea-ice fraction
     real(kind=kind_real), allocatable :: hicen(:,:,:)        !< Sea-ice thickness
     real(kind=kind_real), allocatable :: vicen(:,:,:)        !< Sea-ice volume
     real(kind=kind_real), allocatable :: hsnon(:,:,:)        !< Snow depth over sea-ice
     real(kind=kind_real), allocatable :: vsnon(:,:,:)        !< Snow volume over sea-ice  
     real(kind=kind_real), allocatable :: tsfcn(:,:,:)        !< Temperature over sea-ice or snow
     real(kind=kind_real), allocatable :: qsnon(:,:,:)        !< Enthalpy of snow
     real(kind=kind_real), allocatable :: sicnk(:,:,:,:)      !< Salin_wity of sea-ice
     real(kind=kind_real), allocatable :: sssoc(:,:)          !< Ocean (surface) Salinity
     real(kind=kind_real), allocatable :: qicnk(:,:,:,:)      !< Enthalpy of sea-ice
     real(kind=kind_real), allocatable :: tlioc(:,:)          !< Liquid ocean temperature 
     real(kind=kind_real), allocatable :: sstoc(:,:)          !< Average temperature of grid cell 
     character(len=5), allocatable :: fldnames(:)             !< Variable identifiers
  end type mom5cice5_field

#define LISTED_TYPE mom5cice5_field

  !> Linked list interface - defines registry_t type
#include "linkedList_i.f"

  !> Global registry
  type(registry_t) :: mom5cice5_field_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine create(self, geom, vars)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    type(mom5cice5_geom),  intent(in)    :: geom
    type(mom5cice5_vars),  intent(in)    :: vars

    self%nx   = geom%nx
    self%ny   = geom%ny
    self%nzo  = geom%nzo
    self%nzi  = geom%nzi
    self%nzs  = geom%nzs
    self%ncat = geom%ncat
    self%gridfname = geom%gridfname
    self%nf   = vars%nv

    allocate(self%lon(self%nx,self%ny))
    allocate(self%lat(self%nx,self%ny))
    allocate(self%mask(self%nx,self%ny))
    allocate(self%icemask(self%nx,self%ny))
    allocate(self%level(self%nzs+self%nzi+self%nzo+2))    

    self%lat = geom%lat
    self%lon = geom%lon
    self%mask = geom%mask
    self%icemask = geom%icemask
    self%level = geom%level

    allocate(self%cicen(self%nx,self%ny,self%ncat))
    allocate(self%hicen(self%nx,self%ny,self%ncat))
    allocate(self%vicen(self%nx,self%ny,self%ncat))
    allocate(self%hsnon(self%nx,self%ny,self%ncat))
    allocate(self%vsnon(self%nx,self%ny,self%ncat))
    allocate(self%tsfcn(self%nx,self%ny,self%ncat))
    allocate(self%qsnon(self%nx,self%ny,self%ncat))
    allocate(self%sicnk(self%nx,self%ny,self%ncat,self%nzi))
    allocate(self%sssoc(self%nx,self%ny))
    allocate(self%qicnk(self%nx,self%ny,self%ncat,self%nzi))
    allocate(self%tlioc(self%nx,self%ny))
    allocate(self%sstoc(self%nx,self%ny))    

    self%cicen=0.0_kind_real
    self%hicen=0.0_kind_real
    self%vicen=0.0_kind_real        
    self%hsnon=0.0_kind_real
    self%vsnon=0.0_kind_real
    self%tsfcn=0.0_kind_real
    self%qsnon=0.0_kind_real
    self%sicnk=0.0_kind_real
    self%sssoc=0.0_kind_real    
    self%qicnk=0.0_kind_real
    self%tlioc=0.0_kind_real    
    self%sstoc=0.0_kind_real

    if (self%nf>12) then
       call abor1_ftn ("mom5cice5_fields:create error number of fields")       
    endif
    allocate(self%fldnames(self%nf))
    self%fldnames(:)=vars%fldnames(:)

    call check(self)

  end subroutine create

  ! ------------------------------------------------------------------------------

  subroutine delete(self)
    implicit none
    type(mom5cice5_field), intent(inout) :: self

    call check(self)

    if (allocated(self%cicen)) deallocate(self%cicen)
    if (allocated(self%hicen)) deallocate(self%hicen)
    if (allocated(self%vicen)) deallocate(self%vicen)        
    if (allocated(self%hsnon)) deallocate(self%hsnon)
    if (allocated(self%vsnon)) deallocate(self%vsnon)
    if (allocated(self%tsfcn)) deallocate(self%tsfcn)    
    if (allocated(self%qsnon)) deallocate(self%qsnon)
    if (allocated(self%sicnk)) deallocate(self%sicnk)
    if (allocated(self%sssoc)) deallocate(self%sssoc)
    if (allocated(self%qicnk)) deallocate(self%qicnk)
    if (allocated(self%tlioc)) deallocate(self%tlioc)
    if (allocated(self%sstoc)) deallocate(self%sstoc)    
    if (allocated(self%fldnames)) deallocate(self%fldnames)

  end subroutine delete

  ! ------------------------------------------------------------------------------

  subroutine zeros(self)
    implicit none
    type(mom5cice5_field), intent(inout) :: self

    call check(self)

    self%cicen=0.0_kind_real
    self%hicen=0.0_kind_real
    self%vicen=0.0_kind_real        
    self%hsnon=0.0_kind_real
    self%vsnon=0.0_kind_real
    self%tsfcn=0.0_kind_real
    self%qsnon=0.0_kind_real
    self%sicnk=0.0_kind_real
    self%sssoc=0.0_kind_real    
    self%qicnk=0.0_kind_real
    self%tlioc=0.0_kind_real    
    self%sstoc=0.0_kind_real

  end subroutine zeros

  ! ------------------------------------------------------------------------------

  subroutine dirac(self, c_conf)
    use iso_c_binding
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    type(c_ptr), intent(in)       :: c_conf   !< Configuration
    integer :: ndir,idir,ildir,ifdir,ioff
    integer,allocatable :: ixdir(:),iydir(:)
    character(len=3) :: idirchar

    call check(self)

    ! Get Diracs positions
    ndir = config_get_int(c_conf,"ndir")
    allocate(ixdir(ndir))
    allocate(iydir(ndir))
    do idir=1,ndir
       write(idirchar,'(i3)') idir
       ixdir(idir) = config_get_int(c_conf,"ixdir("//trim(adjustl(idirchar))//")")
       iydir(idir) = config_get_int(c_conf,"iydir("//trim(adjustl(idirchar))//")")
    end do
    !ildir = config_get_int(c_conf,"ildir")
    !ifdir = config_get_int(c_conf,"ifdir")

    ! Check 
    !if (ndir<1) call abor1_ftn("qg_fields:dirac non-positive ndir")
    !if (any(ixdir<1).or.any(ixdir>self%nx)) call abor1_ftn("qg_fields:dirac invalid ixdir")
    !if (any(iydir<1).or.any(iydir>self%ny)) call abor1_ftn("qg_fields:dirac invalid iydir")
    !if ((ildir<1).or.(ildir>self%nl)) call abor1_ftn("qg_fields:dirac invalid ildir")
    !if ((ifdir<1).or.(ifdir>self%nf)) call abor1_ftn("qg_fields:dirac invalid ifdir")

    ! Setup Diracs
    call zeros(self)
    !ioff = (ifdir-1)*self%nl
    do idir=1,ndir
       self%qicnk(ixdir(idir),iydir(idir),1,4) = 1.0 ! Surface temp incr for cat 1
    end do

  end subroutine dirac

  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

  subroutine random(self)

    implicit none
    type(mom5cice5_field), intent(inout) :: self

    call check(self)
    call random_number(self%cicen); self%cicen=self%cicen-sum(self%cicen)

  end subroutine random

  ! ------------------------------------------------------------------------------

  subroutine copy(self,rhs)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    type(mom5cice5_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%lat = rhs%lat
    self%lon = rhs%lon
    self%level = rhs%level
    self%cicen = rhs%cicen
    self%hicen = rhs%hicen
    self%vicen = rhs%vicen
    self%hsnon = rhs%hsnon
    self%vsnon = rhs%vsnon
    self%tsfcn = rhs%tsfcn
    self%qsnon = rhs%qsnon
    self%sicnk = rhs%sicnk
    self%sssoc = rhs%sssoc
    self%qicnk = rhs%qicnk
    self%tlioc = rhs%tlioc
    self%sstoc = rhs%sstoc

    return
  end subroutine copy

  ! ------------------------------------------------------------------------------

  subroutine self_add(self,rhs)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    type(mom5cice5_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen=self%cicen+rhs%cicen
    self%hicen=self%hicen+rhs%hicen
    self%vicen=self%vicen+rhs%vicen
    self%hsnon=self%hsnon+rhs%hsnon
    self%vsnon=self%vsnon+rhs%vsnon
    self%tsfcn=self%tsfcn+rhs%tsfcn
    self%qsnon=self%qsnon+rhs%qsnon
    self%sicnk=self%sicnk+rhs%sicnk
    self%sssoc=self%sssoc+rhs%sssoc
    self%qicnk=self%qicnk+rhs%qicnk
    self%tlioc=self%tlioc+rhs%tlioc
    self%sstoc=self%sstoc+rhs%sstoc

    return
  end subroutine self_add

  ! ------------------------------------------------------------------------------

  subroutine self_schur(self,rhs)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    type(mom5cice5_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen=self%cicen*rhs%cicen
    self%hicen=self%hicen*rhs%hicen
    self%vicen=self%vicen*rhs%vicen
    self%hsnon=self%hsnon*rhs%hsnon
    self%vsnon=self%vsnon*rhs%vsnon
    self%tsfcn=self%tsfcn*rhs%tsfcn
    self%qsnon=self%qsnon*rhs%qsnon
    self%sicnk=self%sicnk*rhs%sicnk
    self%sssoc=self%sssoc*rhs%sssoc
    self%qicnk=self%qicnk*rhs%qicnk
    self%tlioc=self%tlioc*rhs%tlioc
    self%sstoc=self%sstoc*rhs%sstoc

    return
  end subroutine self_schur

  ! ------------------------------------------------------------------------------

  subroutine self_sub(self,rhs)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    type(mom5cice5_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen=self%cicen-rhs%cicen
    self%hicen=self%hicen-rhs%hicen
    self%vicen=self%vicen-rhs%vicen
    self%hsnon=self%hsnon-rhs%hsnon
    self%vsnon=self%vsnon-rhs%vsnon
    self%tsfcn=self%tsfcn-rhs%tsfcn
    self%qsnon=self%qsnon-rhs%qsnon
    self%sicnk=self%sicnk-rhs%sicnk
    self%sssoc=self%sssoc-rhs%sssoc
    self%qicnk=self%qicnk-rhs%qicnk
    self%tlioc=self%tlioc-rhs%tlioc
    self%sstoc=self%sstoc-rhs%sstoc

    return
  end subroutine self_sub

  ! ------------------------------------------------------------------------------

  subroutine self_mul(self,zz)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    real(kind=kind_real), intent(in) :: zz

    call check(self)

    self%cicen = zz * self%cicen
    self%hicen = zz * self%hicen
    self%vicen = zz * self%vicen
    self%hsnon = zz * self%hsnon
    self%vsnon = zz * self%vsnon
    self%tsfcn = zz * self%tsfcn
    self%qsnon = zz * self%qsnon
    self%sicnk = zz * self%sicnk
    self%sssoc = zz * self%sssoc
    self%qicnk = zz * self%qicnk
    self%tlioc = zz * self%tlioc
    self%sstoc = zz * self%sstoc

    return
  end subroutine self_mul

  ! ------------------------------------------------------------------------------

  subroutine axpy(self,zz,rhs)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    real(kind=kind_real), intent(in) :: zz
    type(mom5cice5_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen=self%cicen + zz * rhs%cicen
    self%hicen=self%hicen + zz * rhs%hicen
    self%vicen=self%vicen + zz * rhs%vicen
    self%hsnon=self%hsnon + zz * rhs%hsnon
    self%vsnon=self%vsnon + zz * rhs%vsnon
    self%tsfcn=self%tsfcn + zz * rhs%tsfcn
    self%qsnon=self%qsnon + zz * rhs%qsnon
    self%sicnk=self%sicnk + zz * rhs%sicnk
    self%sssoc=self%sssoc + zz * rhs%sssoc
    self%qicnk=self%qicnk + zz * rhs%qicnk
    self%tlioc=self%tlioc + zz * rhs%tlioc
    self%sstoc=self%sstoc + zz * rhs%sstoc        

    return
  end subroutine axpy

  ! ------------------------------------------------------------------------------

  subroutine dot_prod(fld1,fld2,zprod)
    implicit none
    type(mom5cice5_field), intent(in) :: fld1, fld2
    real(kind=kind_real), intent(out) :: zprod
    integer :: jj, kk
    call check_resolution(fld1, fld2)
    if (fld1%nf /= fld2%nf .or. fld1%nzi /= fld2%nzi) then
       call abor1_ftn("mom5cice5_fields:field_prod error number of fields")
    endif

    zprod = 0.0_kind_real
    do jj = 1, fld1%ncat
       zprod=sum(fld1%cicen(:,:,jj)*fld2%cicen(:,:,jj)*fld1%icemask) + &
         sum(fld1%hicen(:,:,jj)*fld2%hicen(:,:,jj)*fld1%icemask) + &
         sum(fld1%vicen(:,:,jj)*fld2%vicen(:,:,jj)*fld1%icemask) + &
         sum(fld1%hsnon(:,:,jj)*fld2%hsnon(:,:,jj)*fld1%icemask) + &
         sum(fld1%vsnon(:,:,jj)*fld2%vsnon(:,:,jj)*fld1%icemask) + &
         sum(fld1%tsfcn(:,:,jj)*fld2%tsfcn(:,:,jj)*fld1%icemask) + &
         sum(fld1%qsnon(:,:,jj)*fld2%qsnon(:,:,jj)*fld1%icemask)
    end do

    do jj = 1, fld1%ncat
       do kk = 1,fld1%nzi
          zprod = zprod + &
         sum(fld1%sicnk(:,:,jj,kk)*fld2%sicnk(:,:,jj,kk)*fld1%icemask) + &
         sum(fld1%qicnk(:,:,jj,kk)*fld2%qicnk(:,:,jj,kk)*fld1%icemask)
       end do
    end do
    zprod = zprod + sum(fld1%sssoc*fld2%sssoc*fld1%mask) + &
                    sum(fld1%tlioc*fld2%tlioc*fld1%mask) + &
                    sum(fld1%sstoc*fld2%sstoc*fld1%mask)

    return
  end subroutine dot_prod

  ! ------------------------------------------------------------------------------

  subroutine add_incr(self,rhs)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    type(mom5cice5_field), intent(in)    :: rhs

    call check(self)
    call check(rhs)

    call self_add(self,rhs)

    return
  end subroutine add_incr

  ! ------------------------------------------------------------------------------

  subroutine diff_incr(lhs,x1,x2)
    implicit none
    type(mom5cice5_field), intent(inout) :: lhs
    type(mom5cice5_field), intent(in)    :: x1
    type(mom5cice5_field), intent(in)    :: x2

    call check(lhs)
    call check(x1)
    call check(x2)

    call zeros(lhs)


    if (x1%nx==x2%nx .and. x1%ny==x2%ny) then
       if (lhs%nx==x1%nx .and. lhs%ny==x1%ny) then
          lhs%cicen = x1%cicen - x2%cicen
          lhs%hicen = x1%hicen - x2%hicen
          lhs%vicen = x1%vicen - x2%vicen
          lhs%hsnon = x1%hsnon - x2%hsnon
          lhs%vsnon = x1%vsnon - x2%vsnon
          lhs%tsfcn = x1%tsfcn - x2%tsfcn
          lhs%qsnon = x1%qsnon - x2%qsnon
          lhs%sicnk = x1%sicnk - x2%sicnk
          lhs%sssoc = x1%sssoc - x2%sssoc
          lhs%qicnk = x1%qicnk - x2%qicnk
          lhs%tlioc = x1%tlioc - x2%tlioc
          lhs%sstoc = x1%sstoc - x2%sstoc         
       else
          call abor1_ftn("mom5cice5_fields:diff_incr: not coded for low res increment yet")
       endif
    else
       call abor1_ftn("mom5cice5_fields:diff_incr: states not at same resolution")
    endif

    return
  end subroutine diff_incr

  ! ------------------------------------------------------------------------------

  subroutine change_resol(fld,rhs)
    implicit none
    type(mom5cice5_field), intent(inout) :: fld
    type(mom5cice5_field), intent(in)    :: rhs
    real(kind=kind_real), allocatable :: ztmp(:,:)
    real(kind=kind_real) :: dy1, dy2, ya, yb, dx1, dx2, xa, xb
    integer :: jx, jy, jf, iy, ia, ib

    call check(fld)
    call check(rhs)
    call copy(fld,rhs) !!!!!!!!! NOT IMPLEMENTED YET !!!!!!!!!!!!!
    !call abor1_ftn("mom5cice5_fields:field_resol: untested code")

    return
  end subroutine change_resol

  ! ------------------------------------------------------------------------------

  subroutine read_file(fld, c_conf, vdate)
    ! Needs more interface clean-up here...
    use iso_c_binding
    use datetime_mod
    use fckit_log_module, only : log
    use netcdf
    use ncutils
    use interface_ncread_fld, only: ncread_fld
    implicit none
    type(mom5cice5_field), intent(inout) :: fld      !< Fields
    type(c_ptr), intent(in)       :: c_conf   !< Configuration
    type(datetime), intent(inout) :: vdate    !< DateTime

    integer, parameter :: iunit=10
    integer, parameter :: max_string_length=800 ! Yuk!
    character(len=max_string_length+50) :: record
    character(len=max_string_length) :: filename
    character(len=20) :: sdate, fmtn
    character(len=4)  :: cnx
    character(len=11) :: fmt1='(X,ES24.16)'
    character(len=1024)  :: buf
    character(len=128)  :: varname, basename
    integer :: ic, iy, il, ix, is, jx, jy, jf, iread, nf, level
    real(kind=kind_real), allocatable :: zz(:), var3d(:,:,:)

    integer :: nx0, ny0
    integer :: nx, ny, varid, ncid
    integer :: start2(2), count2(2)
    integer :: start3(3), count3(3)
    integer :: start4(4), count4(4)    

    nx0=1 !20
    ny0=350 !60

    print *,fld%nx, fld%ny
    print *,shape(fld%cicen)
    print *,fld%lon(10,10)
    iread = 0
    if (config_element_exists(c_conf,"read_from_file")) then
       iread = config_get_int(c_conf,"read_from_file")
    endif
    print *,'IREAD=',iread
    if (iread==0) then
       call log%warning("mom5cice5_fields:read_file: Inventing State")
       call invent_state(fld,c_conf)
       sdate = config_get_string(c_conf,len(sdate),"date")
       WRITE(buf,*) 'validity date is: '//sdate
       call log%info(buf)
       call datetime_set(sdate, vdate)
    else
       !iread = 0
       sdate = config_get_string(c_conf,len(sdate),"date")
       WRITE(buf,*) 'validity date is: '//sdate
       call log%info(buf)
       call datetime_set(sdate, vdate)       

       ! Read Sea-Ice
       fld%cicefname = config_get_string(c_conf, len(fld%cicefname), "cicefname")
       !WRITE(buf,*) 'cice fname:',fld%cicefname
       print *, 'cice fname:',fld%cicefname
       start3 = (/nx0,ny0,1/)
       count3 = (/fld%nx,fld%ny,fld%ncat/)
       varname='aicen'; call ncread_fld(fld%cicefname, varname, fld%cicen, fld%nx, fld%ny, fld%ncat, start3, count3)
       varname='vicen'; call ncread_fld(fld%cicefname, varname, fld%vicen, fld%nx, fld%ny, fld%ncat, start3, count3)
       varname='vsnon'; call ncread_fld(fld%cicefname, varname, fld%vicen, fld%nx, fld%ny, fld%ncat, start3, count3)
       varname='Tsfcn'; call ncread_fld(fld%cicefname, varname, fld%tsfcn, fld%nx, fld%ny, fld%ncat, start3, count3)
       allocate(var3d(fld%nx,fld%ny,fld%ncat))
       do level=1,fld%nzi
          basename='qice'; call fld_name_int2str(basename, level, varname)
          call ncread_fld(fld%cicefname, varname, var3d, fld%nx, fld%ny, fld%ncat, start3, count3)
          fld%qicnk(:,:,:,level)=var3d
          basename='sice'; call fld_name_int2str(basename, level, varname)
          call ncread_fld(fld%cicefname, varname, var3d, fld%nx, fld%ny, fld%ncat, start3, count3)
          fld%sicnk(:,:,:,level)=var3d
       end do
       deallocate(var3d)
       !do level=1,fld%nzs !!!!!!! CURRENTLY HARD CODED FOR NZS=1 !!!!!!!!!!!!!!!!
       basename='qsno'; call fld_name_int2str(basename, 1, varname)
       print *, varname
       call ncread_fld(fld%cicefname, varname, fld%qsnon, fld%nx, fld%ny, fld%ncat, start3, count3)

       ! Read Ocean
       fld%momfname = config_get_string(c_conf, len(fld%momfname), "momfname")
       print *,'mom fname:',fld%momfname
       start4 = (/nx0,ny0,1,1/)
       count4 = (/fld%nx,fld%ny,fld%nzo,1/)
       varname='temp'
       call ncread_fld(fld%momfname, varname, fld%sstoc, fld%nx, fld%ny, start4, count4)
       varname='salt'
       call ncread_fld(fld%momfname, varname, fld%sssoc, fld%nx, fld%ny, start4, count4)

    endif

    call check(fld)

    return
  end subroutine read_file

  ! ------------------------------------------------------------------------------

  subroutine write_file(fld, c_conf, vdate)
    use iso_c_binding
    use datetime_mod
    use fckit_log_module, only : log
    use netcdf
    use ncutils

    implicit none
    type(mom5cice5_field), intent(inout) :: fld    !< Fields
    type(c_ptr), intent(in)    :: c_conf           !< Configuration
    type(datetime), intent(inout) :: vdate         !< DateTime
    integer, parameter :: max_string_length=800    ! Yuk!
    character(len=max_string_length) :: filename
    character(len=128)  :: varname
    character(len=20) :: sdate
    real(kind=8) :: missing=-999d0

    integer :: ncid, varid, dimids2d(2), dimids4d(3)
    integer :: jx, jy, varid_lon, varid_lat, varid_sst
    integer :: x_dimid, y_dimid, z_dimid, cat_dimid, status
    integer :: catnum

    catnum=1 !!!!!!!!!!!! HARD CODED CATEGORY !!!!!!!!!!!!!!!!!!!!

    call check(fld)

    print *,'writing stuf ................'
    
    filename = config_get_string(c_conf, len(filename), "filename")
    varname = config_get_string(c_conf, len(varname), "varname")

    call nccheck( nf90_create(filename, nf90_clobber, ncid) )
    call nccheck( nf90_def_dim(ncid, "xaxis_1", fld%nx, x_dimid) )
    call nccheck( nf90_def_dim(ncid, "yaxis_1", fld%ny, y_dimid) )
    call nccheck( nf90_def_dim(ncid, "zaxis_1", fld%nzi, z_dimid) )
    !call nccheck( nf90_def_dim(ncid, "cataxis_1", fld%ncat, cat_dimid) )
    !dimids4d =  (/ x_dimid, y_dimid, cat_dimid, z_dimid /)
    dimids4d =  (/ x_dimid, y_dimid, z_dimid /)
    dimids2d =  (/ x_dimid, y_dimid /)

    print *,'writing ....', fld%nx, fld%ny
    
    do jx=1,fld%nx
       do jy=1,fld%ny
          if (nint(fld%icemask(jx,jy))==0) fld%qicnk(jx,jy,:,:) = missing
          if (nint(fld%mask(jx,jy))==0) fld%sstoc(jx,jy) = missing
       end do
    end do

    call nccheck( nf90_def_var(ncid, 'qicnk', nf90_double, dimids4d, varid) )
    call nccheck( nf90_put_att(ncid, varid, '_FillValue', missing) )
    call nccheck( nf90_def_var(ncid, 'sstoc', nf90_double, dimids2d, varid_sst) )
    call nccheck( nf90_put_att(ncid, varid_sst, '_FillValue', missing) )
    call nccheck( nf90_def_var(ncid, 'lat', nf90_double, dimids2d, varid_lat) )
    call nccheck( nf90_def_var(ncid, 'lon', nf90_double, dimids2d, varid_lon) )        
    call nccheck( nf90_enddef(ncid) )
    !call nccheck( nf90_put_var(ncid, varid, fld%qicnk(:,:,catnum,:)))
    call nccheck( nf90_put_var(ncid, varid, fld%qicnk(:,:,catnum,:)))    
    call nccheck( nf90_put_var(ncid, varid_sst, fld%sstoc))
    call nccheck( nf90_put_var(ncid, varid_lat, fld%lat))
    call nccheck( nf90_put_var(ncid, varid_lon, fld%lon))        
    !print *,'doneriting ....'
    !call nccheck( nf90_def_var(ncid, varname, nf90_double, dimids2d, varid) )
    !call nccheck( nf90_enddef(ncid) )
    !call nccheck( nf90_put_var(ncid, varid, fld%sstoc ) )
    call nccheck( nf90_close(ncid) )

    print *,'in write field ...............'
    
    call datetime_to_string(vdate, sdate)

    return
  end subroutine write_file

  ! ------------------------------------------------------------------------------

  subroutine gpnorm(fld, nf, pstat) 
    implicit none
    type(mom5cice5_field), intent(in) :: fld
    integer, intent(in) :: nf
    real(kind=kind_real), intent(inout) :: pstat(3, nf) !> [average, min, max]
    real(kind=kind_real) :: zz
    integer :: jj,joff

    call check(fld)
    
    !pstat(1,:)=minval(fld%cicen)
    !pstat(2,:)=maxval(fld%cicen)
    !pstat(3,:)=abs(maxval(fld%cicen)-minval(fld%cicen))
    
    !call abor1_ftn("mom5cice5_fields_gpnorm: error not implemented")
    !print *,'pstat=',pstat
    call dot_prod(fld,fld,zz)    

    pstat = sqrt(zz)
    !call random_number(pstat)
    
    return
  end subroutine gpnorm

  ! ------------------------------------------------------------------------------

  subroutine fldrms(fld, prms)
    implicit none
    type(mom5cice5_field), intent(in) :: fld
    real(kind=kind_real), intent(out) :: prms
    integer :: jf,jy,jx,ii
    real(kind=kind_real) :: zz

    call check(fld)

    !zz = 0.0_kind_real

    !do jf=1,fld%nl*fld%nf
    !   do jy=1,fld%ny
    !      do jx=1,fld%nx
    !         zz = zz + fld%gfld3d(jx,jy,jf)*fld%gfld3d(jx,jy,jf)
    !      enddo
    !   enddo
    !enddo

    !ii = fld%nl*fld%nf*fld%ny*fld%nx
    call dot_prod(fld,fld,prms)
    
    !prms = sqrt(zz/real(ii,kind_real))
    !prms = 1.0_kind_real

  end subroutine fldrms

  ! ------------------------------------------------------------------------------

  subroutine lin_weights(kk,delta1,delta2,k1,k2,w1,w2)
    implicit none
    integer, intent(in)  :: kk
    real(kind=kind_real), intent(in)     :: delta1,delta2
    integer, intent(out) :: k1,k2
    real(kind=kind_real), intent(out)    :: w1,w2

    integer :: ii
    real(kind=kind_real) :: zz

    zz=real(kk-1,kind_real)*delta1
    zz=zz/delta2
    ii=int(zz)
    w1=zz-real(ii,kind_real)
    w2=1.0_kind_real-w1
    k1=ii+1
    k2=ii+2

    return
  end subroutine lin_weights

  ! ------------------------------------------------------------------------------

  subroutine convert_to_ug(self, ug)
    use unstructured_grid_mod
    use mom5cice5_thermo
    
    implicit none
    type(mom5cice5_field), intent(in) :: self
    type(unstructured_grid), intent(inout) :: ug
    real(kind=kind_real), allocatable :: zz(:)
    real(kind=kind_real), allocatable :: vv(:)    
    integer, allocatable :: cmask(:)
    integer :: jx,jy,jz,jk
    integer :: nz_total     ! Total number of levels in the 3D fields
    integer :: n_vars       ! Number of 3D variables 
    integer :: n_surf_vars  ! Number of surf vars (sould be 0 for ocean/ice)
    integer :: cat_num      ! !!!!!!!!! only doing 1 category for now !!!!!!!!!!!

    cat_num = 1
    nz_total = size(self%level)
    allocate(zz(nz_total))
    allocate(vv(nz_total))
    allocate(cmask(nz_total))
    do jz = 1,nz_total
       zz(jz) = real(self%level(jz))
    end do
    call create_unstructured_grid(ug, nz_total, zz)

    n_vars = 1      !!!!! START WITH ONLY ONE VAR !!!!!!!!! 
    n_surf_vars = 0 !!!!! NO SURFACE VAR !!!!!!!!! 

    do jy=1,self%ny
       do jx=1,self%nx
          jk = 1
          cmask(jk) = int(self%icemask(jx,jy))       ! Surface T
          vv(jk) = self%tsfcn(jx,jy,cat_num)
          jk = jk + 1
          do jz = 1,self%nzs                         ! Snow T
             cmask(jk) = int(self%icemask(jx,jy))    !
             !vv(jk) = Ts_nl(self%qsnon(jx,jy,cat_num))
             vv(jk) = self%qsnon(jx,jy,cat_num)             
             jk = jk + 1
          end do
          do jz = 1,self%nzi                         ! Ice T
             cmask(jk) = int(self%icemask(jx,jy))    !
             !vv(jk) = Ti_nl(self%qicnk(jx,jy,cat_num,jz),self%sicnk(jx,jy,cat_num,jz))
             vv(jk) = self%qicnk(jx,jy,cat_num,jz)
             jk = jk + 1
          end do
          cmask(jk) = int(self%icemask(jx,jy))       ! Ice/Ocean interface
          !vv(jk) = Tm(self%sssoc(jx,jy))             ! Tf = -mu * S
          vv(jk) = self%sssoc(jx,jy)             ! Tf = -mu * S          
          jk = jk + 1          
          do jz = 1,self%nzo                         ! Ocean
             !cmask(jk) = int(self%mask(jx,jy))       !
             cmask(jk) = int(self%icemask(jx,jy))   !             
             vv(jk) = self%sstoc(jx,jy)              !
             jk = jk + 1
          end do

          cmask(:) = int(self%mask(jx,jy))           ! A corriger !!!
          !print *,'cmask=',cmask
          !if (self%icemask(jx,jy)>0.0) then
             !print *,vv(:)
             !read(*,*)
          !end if
          !if (cmask(1)==1) read(*,*)
          
          call add_column(ug, self%lat(jx,jy), self%lon(jx,jy), &
               nz_total, &
               n_vars, &
               n_surf_vars, &
               cmask, &
               0)
          ug%last%column%cols(:) = vv(:)
       enddo
    enddo
  end subroutine convert_to_ug

  ! ------------------------------------------------------------------------------

  subroutine convert_from_ug(self, ug)
    use unstructured_grid_mod
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    type(unstructured_grid), intent(in) :: ug
    type(column_element), pointer :: current
    real(kind=kind_real) :: dx, dy
    integer :: jx,jy,jz,jk
    integer :: n_vars       ! Number of 3D variables 
    integer :: n_surf_vars  ! Number of surf vars (sould be 0 for ocean/ice)
    integer :: cat_num      ! !!!!!!!!! only doing 1 category for now !!!!!!!!!!!

!!!!!!!!! code inverse of convert_to_ug !!!!!!!!!!!!!!!

    print *,'in convert from ug ............'

    
    current => ug%head
    cat_num = 1 
    n_vars = 1      !!!!! START WITH ONLY ONE VAR !!!!!!!!! 
    n_surf_vars = 0 !!!!! NO SURFACE VAR !!!!!!!!! 

    do jy=1,self%ny
       do jx=1,self%nx
          jk = 1
          self%tsfcn(jx,jy,cat_num) = current%column%cols(jk)         ! Tsfcs
          jk = jk + 1
          do jz = 1,self%nzs                                          ! Q Snow
             self%qsnon(jx,jy,cat_num) = current%column%cols(jk)
             jk = jk + 1
          end do
          do jz = 1,self%nzi                                          ! Q Ice
             self%qicnk(jx,jy,cat_num,jz) = current%column%cols(jk)
             jk = jk + 1
          end do
          self%sssoc(jx,jy) = current%column%cols(jk)                 ! Ice/Ocean interface,
                                                                      ! Tf = -mu * S          
          jk = jk + 1          
          do jz = 1,self%nzo                                          ! Ocean SST
             self%sstoc(jx,jy) = current%column%cols(jk)              !
             jk = jk + 1
          end do          
          current => current%next
       enddo
       print*, 'from',jy
    enddo

  end subroutine convert_from_ug

  ! ------------------------------------------------------------------------------

  function common_vars(x1, x2)

    implicit none
    type(mom5cice5_field), intent(in) :: x1, x2
    integer :: common_vars
    integer :: jf

    ! We assume here that one set of fields is a subset of the other,
    ! that fields are always in the same order starting with x,
    ! and that the common fields are the first ones.

    common_vars = min(x1%nf, x2%nf)
    do jf = 1, common_vars
       if (x1%fldnames(jf)/=x2%fldnames(jf)) &
            & call abor1_ftn("common_vars: fields do not match")
    enddo
    if (x1%nzi /= x2%nzi) call abor1_ftn("common_vars: error number of levels")
    common_vars = x1%nzi * common_vars

  end function common_vars

  ! ------------------------------------------------------------------------------

  subroutine check_resolution(x1, x2)

    implicit none
    type(mom5cice5_field), intent(in) :: x1, x2

    ! NEEDS WORK !!!
    if (x1%nx /= x2%nx .or.  x1%ny /= x2%ny ) then
       call abor1_ftn ("mom5cice5_fields: resolution error")
    endif
    call check(x1)
    call check(x2)

  end subroutine check_resolution

  ! ------------------------------------------------------------------------------

  subroutine check(self)
    implicit none
    type(mom5cice5_field), intent(in) :: self
    logical :: bad

    bad = .false.
    bad = bad .or. (size(self%cicen, 1) /= self%nx)

    ! add more test here ...

    if (bad) then
       write(0,*)'nx, ny, nf, nzi, nzo = ',self%nx,self%ny,self%nf,self%nzi,self%nzo
       call abor1_ftn ("mom5cice5_fields: field not consistent")
    endif

  end subroutine check

  ! ------------------------------------------------------------------------------

end module mom5cice5_fields
