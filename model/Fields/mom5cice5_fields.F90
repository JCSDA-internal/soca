
!> Handle fields for the  model

module mom5cice5_fields

  use config_mod
  use mom5cice5_geom_mod
  use mom5cice5_goms_mod
  use mom5cice5_locs_mod  
  use mom5cice5_vars_mod
  use kinds

  implicit none
  private

  public :: mom5cice5_field, &
       & create, delete, zeros, dirac, random, copy, &
       & self_add, self_schur, self_sub, self_mul, axpy, &
       & dot_prod, add_incr, diff_incr, &
       & read_file, write_file, gpnorm, fldrms, &
       & change_resol, interp_tl, interp_ad, convert_to_ug, convert_from_ug
  public :: mom5cice5_field_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to hold fields
  type :: mom5cice5_field
     type(mom5cice5_geom), pointer     :: geom                  !< MOM5 & CICE5 Geometry
     integer                           :: nf                    !< Number of fields
     character(len=128)                :: gridfname             !< Grid file name
     character(len=128)                :: cicefname             !< Fields file name for cice
     character(len=128)                :: momfname              !< Fields file name for mom
     real(kind=kind_real), allocatable :: cicen(:,:,:)          !< Sea-ice fraction                 (nx,ny,ncat)
     real(kind=kind_real), allocatable :: hicen(:,:,:)          !< Sea-ice thickness                (nx,ny,ncat)
     real(kind=kind_real), allocatable :: vicen(:,:,:)          !< Sea-ice volume                   (nx,ny,ncat)
     real(kind=kind_real), allocatable :: hsnon(:,:,:)          !< Snow depth over sea-ice          (nx,ny,ncat)
     real(kind=kind_real), allocatable :: vsnon(:,:,:)          !< Snow volume over sea-ice         (nx,ny,ncat) 
     real(kind=kind_real), allocatable :: tsfcn(:,:,:)          !< Temperature over sea-ice or snow (nx,ny,ncat)
     real(kind=kind_real), allocatable :: qsnon(:,:,:)          !< Enthalpy of snow                 (nx,ny,ncat)
     real(kind=kind_real), allocatable :: sicnk(:,:,:,:)        !< Salin_wity of sea-ice            (nx,ny,ncat,nzi)
     real(kind=kind_real), allocatable :: sssoc(:,:)            !< Ocean (surface) Salinity         (nx,ny,nzo)
     real(kind=kind_real), allocatable :: qicnk(:,:,:,:)        !< Enthalpy of sea-ice              (nx,ny,ncat,nzi)
     real(kind=kind_real), allocatable :: sstoc(:,:)            !< Average temperature of grid cell (nx,ny,nzo)
     character(len=5), allocatable     :: fldnames(:)           !< Variable identifiers             (nf)
     integer, allocatable              :: numfld_per_fldname(:) !< Number of 2d fields for each     (nf) 
                                                                !< element of fldnames 
  end type mom5cice5_field

#define LISTED_TYPE mom5cice5_field

  !> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: mom5cice5_field_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "util/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine create(self, geom, vars)
    implicit none
    type(mom5cice5_field), intent(inout)          :: self
    type(mom5cice5_geom),  pointer, intent(in)    :: geom
    type(mom5cice5_vars),  intent(in)          :: vars        

    integer :: ivar

    self%geom => geom
    self%gridfname = geom%gridfname
    self%nf   = vars%nv

    allocate(self%numfld_per_fldname(vars%nv))
    
    do ivar=1,vars%nv
       select case(vars%fldnames(ivar))
       case ('cicen','hicen','vicen','hsnon','vsnon','tsfcn')
          self%numfld_per_fldname(ivar)=geom%ncat
       case ('sicnk','qicnk')
          self%numfld_per_fldname(ivar)=geom%ncat*geom%nzi
       case ('qsnon')
          self%numfld_per_fldname(ivar)=geom%ncat*geom%nzs
       case ('sssoc','sstoc')
         self%numfld_per_fldname(ivar)=geom%nzo
       case default
          call abor1_ftn("c_mom5cice5_fields: undefined variables")
       end select
    end do

    
    allocate(self%cicen(self%geom%nx,self%geom%ny,self%geom%ncat))
    allocate(self%hicen(self%geom%nx,self%geom%ny,self%geom%ncat))
    allocate(self%vicen(self%geom%nx,self%geom%ny,self%geom%ncat))
    allocate(self%hsnon(self%geom%nx,self%geom%ny,self%geom%ncat))
    allocate(self%vsnon(self%geom%nx,self%geom%ny,self%geom%ncat))
    allocate(self%tsfcn(self%geom%nx,self%geom%ny,self%geom%ncat))
    allocate(self%qsnon(self%geom%nx,self%geom%ny,self%geom%ncat))
    allocate(self%sicnk(self%geom%nx,self%geom%ny,self%geom%ncat,self%geom%nzi))
    allocate(self%sssoc(self%geom%nx,self%geom%ny))
    allocate(self%qicnk(self%geom%nx,self%geom%ny,self%geom%ncat,self%geom%nzi))
    allocate(self%sstoc(self%geom%nx,self%geom%ny))    

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
    self%sstoc=0.0_kind_real

    if (self%nf>11) then
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
    !if (any(ixdir<1).or.any(ixdir>self%geom%nx)) call abor1_ftn("qg_fields:dirac invalid ixdir")
    !if (any(iydir<1).or.any(iydir>self%geom%ny)) call abor1_ftn("qg_fields:dirac invalid iydir")
    !if ((ildir<1).or.(ildir>self%nl)) call abor1_ftn("qg_fields:dirac invalid ildir")
    !if ((ifdir<1).or.(ifdir>self%nf)) call abor1_ftn("qg_fields:dirac invalid ifdir")

    ! Setup Diracs
    call zeros(self)
    !ioff = (ifdir-1)*self%nl
    do idir=1,ndir
       !self%qicnk(ixdir(idir),iydir(idir),1,4) = 1.0 ! Surface temp incr for cat 1
       !self%tsfcn(ixdir(idir),iydir(idir),1) = 1.0 ! Surface temp incr for cat 1
       self%sstoc(ixdir(idir),iydir(idir)) = 1.0 ! Surface temp incr for cat 1
       !self%cicen(ixdir(idir),iydir(idir),3) = 1.0 ! Surface temp incr for cat 1
    end do

  end subroutine dirac

  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

  subroutine random(self)

    implicit none
    type(mom5cice5_field), intent(inout) :: self

    call check(self)
    call random_number(self%cicen); self%cicen=self%cicen-sum(self%cicen) !<--- NO GOOD !!!!

  end subroutine random

  ! ------------------------------------------------------------------------------

  subroutine copy(self,rhs)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    type(mom5cice5_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

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
    if (fld1%nf /= fld2%nf .or. fld1%geom%nzi /= fld2%geom%nzi) then
       call abor1_ftn("mom5cice5_fields:field_prod error number of fields")
    endif

    zprod = 0.0_kind_real
    do jj = 1, fld1%geom%ncat
       zprod=sum(fld1%cicen(:,:,jj)*fld2%cicen(:,:,jj)*fld1%geom%icemask) + &
            sum(fld1%hicen(:,:,jj)*fld2%hicen(:,:,jj)*fld1%geom%icemask) + &
            sum(fld1%vicen(:,:,jj)*fld2%vicen(:,:,jj)*fld1%geom%icemask) + &
            sum(fld1%hsnon(:,:,jj)*fld2%hsnon(:,:,jj)*fld1%geom%icemask) + &
            sum(fld1%vsnon(:,:,jj)*fld2%vsnon(:,:,jj)*fld1%geom%icemask) + &
            sum(fld1%tsfcn(:,:,jj)*fld2%tsfcn(:,:,jj)*fld1%geom%icemask) + &
            sum(fld1%qsnon(:,:,jj)*fld2%qsnon(:,:,jj)*fld1%geom%icemask)
    end do

    do jj = 1, fld1%geom%ncat
       do kk = 1,fld1%geom%nzi
          zprod = zprod + &
               sum(fld1%sicnk(:,:,jj,kk)*fld2%sicnk(:,:,jj,kk)*fld1%geom%icemask) + &
               sum(fld1%qicnk(:,:,jj,kk)*fld2%qicnk(:,:,jj,kk)*fld1%geom%icemask)
       end do
    end do
    zprod = zprod + sum(fld1%sssoc*fld2%sssoc*fld1%geom%mask) + &
         sum(fld1%sstoc*fld2%sstoc*fld1%geom%mask)
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


    if (x1%geom%nx==x2%geom%nx .and. x1%geom%ny==x2%geom%ny) then
       if (lhs%geom%nx==x1%geom%nx .and. lhs%geom%ny==x1%geom%ny) then
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
    call copy(fld,rhs)
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
    use mom5cice5_thermo

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
    integer :: ic, iy, il, ix, is, jx, jy, jk, jcat, jf, iread, nf, level
    real(kind=kind_real), allocatable :: zz(:), var3d(:,:,:)

    integer :: nx0, ny0
    integer :: nx, ny, varid, ncid
    integer :: start2(2), count2(2)
    integer :: start3(3), count3(3)
    integer :: start4(4), count4(4)    

    nx0=1 !20
    ny0=1 !60

    iread = 0
    if (config_element_exists(c_conf,"read_from_file")) then
       iread = config_get_int(c_conf,"read_from_file")
    endif
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
       count3 = (/fld%geom%nx,fld%geom%ny,fld%geom%ncat/)
       varname='aicen'; call ncread_fld(fld%cicefname, varname, fld%cicen, fld%geom%nx, fld%geom%ny, fld%geom%ncat, start3, count3)
       varname='vicen'; call ncread_fld(fld%cicefname, varname, fld%vicen, fld%geom%nx, fld%geom%ny, fld%geom%ncat, start3, count3)
       varname='vsnon'; call ncread_fld(fld%cicefname, varname, fld%vicen, fld%geom%nx, fld%geom%ny, fld%geom%ncat, start3, count3)
       varname='Tsfcn'; call ncread_fld(fld%cicefname, varname, fld%tsfcn, fld%geom%nx, fld%geom%ny, fld%geom%ncat, start3, count3)
       allocate(var3d(fld%geom%nx,fld%geom%ny,fld%geom%ncat))
       do level=1,fld%geom%nzi
          basename='qice'; call fld_name_int2str(basename, level, varname)
          call ncread_fld(fld%cicefname, varname, var3d, fld%geom%nx, fld%geom%ny, fld%geom%ncat, start3, count3)
          fld%qicnk(:,:,:,level)=var3d
          basename='sice'; call fld_name_int2str(basename, level, varname)
          call ncread_fld(fld%cicefname, varname, var3d, fld%geom%nx, fld%geom%ny, fld%geom%ncat, start3, count3)
          fld%sicnk(:,:,:,level)=var3d
       end do

       deallocate(var3d)
       !do level=1,fld%geom%nzs !!!!!!! CURRENTLY HARD CODED FOR NZS=1 !!!!!!!!!!!!!!!!
       basename='qsno'; call fld_name_int2str(basename, 1, varname)
       print *, varname
       call ncread_fld(fld%cicefname, varname, fld%qsnon, fld%geom%nx, fld%geom%ny, fld%geom%ncat, start3, count3)

       ! Read Ocean
       fld%momfname = config_get_string(c_conf, len(fld%momfname), "momfname")
       print *,'mom fname:',fld%momfname
       start4 = (/nx0,ny0,1,1/)
       count4 = (/fld%geom%nx,fld%geom%ny,fld%geom%nzo,1/)
       varname='temp'
       call ncread_fld(fld%momfname, varname, fld%sstoc, fld%geom%nx, fld%geom%ny, start4, count4)
       varname='salt'
       call ncread_fld(fld%momfname, varname, fld%sssoc, fld%geom%nx, fld%geom%ny, start4, count4)

       print *, 'OUT OF READ_FILE'
    endif

    ! Enthalpy to temperature
    do jx = 1,fld%geom%nx
       do jy = 1,fld%geom%ny
          do jk = 1,fld%geom%nzi
             do jcat = 1,fld%geom%ncat             
                fld%qicnk(jx,jy,jcat,jk) = Ti_nl(fld%qicnk(jx,jy,jcat,jk),fld%sicnk(jx,jy,jcat,jk))
             end do
          end do
       end do
    end do

    call check(fld)

    return
  end subroutine read_file

  ! ------------------------------------------------------------------------------

  subroutine write_file(fld, c_conf, vdate)
    use iso_c_binding
    use datetime_mod
    use fckit_log_module, only : fckit_log
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
    character(len=1024):: buf

    integer :: ncid, varid, dimids2d(2), dimids3d(3), dimids4d(3)
    integer :: jx, jy, varid_lon, varid_lat, varid_sst, varid_cicen, varid_tsfcn
    integer :: x_dimid, y_dimid, z_dimid, cat_dimid, status
    integer :: catnum

    catnum=1 !!!!!!!!!!!! HARD CODED CATEGORY !!!!!!!!!!!!!!!!!!!!

    print *,'============== IN WRITE FILE =================='

    call check(fld)

    !filename = config_get_string(c_conf, len(filename), "filename")
    !varname = config_get_string(c_conf, len(varname), "varname")

    filename = genfilename(c_conf,max_string_length,vdate)
    WRITE(buf,*) 'field:write_file: writing '//filename
    call fckit_log%info(buf)

    call nccheck( nf90_create(filename, nf90_clobber, ncid) )
    call nccheck( nf90_def_dim(ncid, "xaxis_1", fld%geom%nx, x_dimid) )
    call nccheck( nf90_def_dim(ncid, "yaxis_1", fld%geom%ny, y_dimid) )
    call nccheck( nf90_def_dim(ncid, "zaxis_1", fld%geom%nzi, z_dimid) )
    call nccheck( nf90_def_dim(ncid, "cataxis_1", fld%geom%ncat, cat_dimid) )

    !call nccheck( nf90_def_dim(ncid, "cataxis_1", fld%geom%ncat, cat_dimid) )
    !dimids4d =  (/ x_dimid, y_dimid, cat_dimid, z_dimid /)
    dimids4d =  (/ x_dimid, y_dimid, z_dimid /)
    dimids3d =  (/ x_dimid, y_dimid, cat_dimid /)
    dimids2d =  (/ x_dimid, y_dimid /)

    do jx=1,fld%geom%nx
       do jy=1,fld%geom%ny
          if (nint(fld%geom%mask(jx,jy))==0) fld%qicnk(jx,jy,:,:) = missing          
          if (nint(fld%geom%mask(jx,jy))==0) fld%sstoc(jx,jy) = missing
       end do
    end do

    call nccheck( nf90_def_var(ncid, 'qicnk', nf90_double, dimids4d, varid) )
    call nccheck( nf90_put_att(ncid, varid, '_FillValue', missing) )
    call nccheck( nf90_def_var(ncid, 'cicen', nf90_double, dimids3d, varid_cicen) )
    call nccheck( nf90_put_att(ncid, varid_cicen, '_FillValue', missing) )
    call nccheck( nf90_def_var(ncid, 'tsfcn', nf90_double, dimids3d, varid_tsfcn) )
    call nccheck( nf90_put_att(ncid, varid_tsfcn, '_FillValue', missing) )        
    call nccheck( nf90_def_var(ncid, 'sstoc', nf90_double, dimids2d, varid_sst) )
    call nccheck( nf90_put_att(ncid, varid_sst, '_FillValue', missing) )
    call nccheck( nf90_def_var(ncid, 'lat', nf90_double, dimids2d, varid_lat) )
    call nccheck( nf90_def_var(ncid, 'lon', nf90_double, dimids2d, varid_lon) )        
    call nccheck( nf90_enddef(ncid) )
    call nccheck( nf90_put_var(ncid, varid_tsfcn, fld%tsfcn))
    call nccheck( nf90_put_var(ncid, varid, fld%qicnk(:,:,catnum,:)))    
    call nccheck( nf90_put_var(ncid, varid_sst, fld%sstoc))
    call nccheck( nf90_put_var(ncid, varid_cicen, fld%cicen))
    call nccheck( nf90_put_var(ncid, varid_lat, fld%geom%lat))
    call nccheck( nf90_put_var(ncid, varid_lon, fld%geom%lon))        

    !call nccheck( nf90_def_var(ncid, varname, nf90_double, dimids2d, varid) )
    !call nccheck( nf90_enddef(ncid) )
    !call nccheck( nf90_put_var(ncid, varid, fld%sstoc ) )
    call nccheck( nf90_close(ncid) )

    call datetime_to_string(vdate, sdate)

    return
  end subroutine write_file

  ! ------------------------------------------------------------------------------

  function genfilename (c_conf,length,vdate)
    use iso_c_binding
    use datetime_mod
    use duration_mod
    type(c_ptr), intent(in)    :: c_conf  !< Configuration
    integer, intent(in) :: length
    character(len=length) :: genfilename
    type(datetime), intent(in) :: vdate

    character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
         & prefix, mmb
    type(datetime) :: rdate
    type(duration) :: step
    integer lenfn

    ! here we should query the length and then allocate "string".
    ! But Fortran 90 does not allow variable-length allocatable strings.
    ! config_get_string checks the string length and aborts if too short.
    fdbdir = config_get_string(c_conf,len(fdbdir),"datadir")
    expver = config_get_string(c_conf,len(expver),"exp")
    typ    = config_get_string(c_conf,len(typ)   ,"type")

    if (typ=="ens") then
       mmb = config_get_string(c_conf, len(mmb), "member")
       lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
       prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ) // "." // TRIM(mmb)
    else
       lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
       prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ)
    endif

    if (typ=="fc" .or. typ=="ens") then
       referencedate = config_get_string(c_conf,len(referencedate),"date")
       call datetime_to_string(vdate, validitydate)
       call datetime_create(TRIM(referencedate), rdate)
       call datetime_diff(vdate, rdate, step)
       call duration_to_string(step, sstep)
       lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
       genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
    endif

    if (typ=="an") then
       call datetime_to_string(vdate, validitydate)
       lenfn = lenfn + 1 + LEN_TRIM(validitydate)
       genfilename = TRIM(prefix) // "." // TRIM(validitydate)
    endif

    if (lenfn>length) &
         & call abor1_ftn("fields:genfilename: filename too long")

  end function genfilename

  ! ------------------------------------------------------------------------------
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

    !print *,'pstat=',pstat

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

    call dot_prod(fld,fld,prms)
    prms = sqrt(prms)

  end subroutine fldrms

  ! ------------------------------------------------------------------------------

  subroutine interp_tl(fld, locs, gom)
    implicit none
    type(mom5cice5_field), intent(in)   :: fld
    type(mom5cice5_locs), intent(in)    :: locs
    type(mom5cice5_goms), intent(inout) :: gom
    character(2)                        :: op_type='TL'

    call check(fld)

    call nicas_interph(fld, locs, gom, op_type)

  end subroutine interp_tl

  ! ------------------------------------------------------------------------------

  subroutine interp_ad(fld, locs, gom)
    implicit none
    type(mom5cice5_field), intent(inout) :: fld
    type(mom5cice5_locs), intent(in)     :: locs
    type(mom5cice5_goms), intent(inout)  :: gom
    character(2)                        :: op_type='AD'

    call check(fld)
    !call nicas_interph(fld, locs, gom, op_type)

  end subroutine interp_ad

  ! ------------------------------------------------------------------------------

  subroutine nicas_interph(fld, locs, gom, op_type)

    use type_linop
    use tools_interp, only: interp_horiz
    use type_randgen, only: rng,initialize_sampling,create_randgen !randgentype
    use module_namelist, only: namtype    
    use tools_const, only : deg2rad

    type(mom5cice5_field), intent(in)    :: fld
    type(mom5cice5_locs), intent(in)     :: locs
    type(mom5cice5_goms), intent(inout)  :: gom
    character(2), intent(in)             :: op_type !('TL' or 'AD')

    integer :: Nc, No, var_type_index, Ncat
    integer :: ivar, gom_dim1, cnt_fld
    character(len=1024)  :: buf
    logical,allocatable :: mask(:), masko(:)               ! < mask (ncells, nlevels)
    real(kind=kind_real), allocatable :: lon(:), lat(:), lono(:), lato(:), fld_src(:), fld_dst(:)
    type(namtype) :: nam !< Namelist variables

    Nc = fld%geom%nx*fld%geom%ny
    No = locs%nloc
    Ncat = fld%geom%ncat
    if (No>0) then
       allocate(lon(Nc), lat(Nc), mask(Nc), fld_src(Nc))    ! <--- Not memory efficient ...
       allocate(masko(No), fld_dst(No), lono(No), lato(No)) ! <--- use pointers?

       masko = .true. ! Figured out what's the use for masko????
       !Some issues with the mask, FIX IT!!!!
       !where(reshape(fld%geom%mask,(/Nc/)).eq.0.0_kind_real)
       mask = .true.
       !end where
       if (.not.(gom%hinterp_initialized)) then
          print *,'INITIALIZE INTERP'
          rng = create_randgen(nam)
          lono = deg2rad*locs%xyz(1,:)
          lato = deg2rad*locs%xyz(2,:)
          lon = deg2rad*reshape(fld%geom%lon, (/Nc/))     ! Inline grid, structured to un-structured
          lat = deg2rad*reshape(fld%geom%lat, (/Nc/))     ! and change to SI Units
          call interp_horiz(rng, Nc, lon,  lat,  mask, &
            No, lono, lato, masko, &
            gom%hinterp_op)
          gom%hinterp_initialized = .true.
       end if

       !Finish Initializing gom
       if (.not.allocated(gom%values)) then
          gom_dim1=sum(fld%numfld_per_fldname(1:gom%nvar)) ! WILL CREATE ISSUES:
                                                           ! Assume the order of var type is preserved
                                                           ! [cicen, hicen, ...]
          allocate(gom%values(gom_dim1,gom%nobs))
       end if
       if (.not.allocated(gom%numfld_per_fldname)) then
          allocate(gom%numfld_per_fldname(gom%nvar))
          gom%numfld_per_fldname=fld%numfld_per_fldname ! Will be used in obs oper          
       end if !probably need to assert shape of gom%values==(gom_dim1,gom%nobs)
          
       select case (op_type)
       case ('TL')
          if (.not.allocated(fld_src)) allocate(fld_src(Nc))
          cnt_fld=0
          !Loop through variable types
          do var_type_index=1,gom%nvar
             !Loop through variable fields
             do ivar=1,fld%numfld_per_fldname(var_type_index)
                cnt_fld=cnt_fld+1
                print *,'Apply interp op to variable:',gom%variables(var_type_index),' field num:',cnt_fld
                fld_src = reshape(fld%cicen(:,:,ivar), (/Nc/))
                call apply_linop(gom%hinterp_op, fld_src, fld_dst)
                gom%values(cnt_fld,gom%used:gom%used+No-1)=fld_dst(1:No)
             end do
          end do
          deallocate(fld_src)
       case ('AD')
          call abor1_ftn("nicas_interph: Wrapper for adjoint not implemented yet")
          !call apply_linop_ad(hinterp_op,fld_dst,fld_src)
          !CHECK: WHERE DO WE TRIGGER THE ADJOINT TESTS? 
          !put fld_src
       end select
       gom%used=gom%used+locs%nloc
    end if
  end subroutine nicas_interph

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
    nz_total = size(self%geom%level)
    allocate(zz(nz_total))
    allocate(vv(nz_total))
    allocate(cmask(nz_total))
    do jz = 1,nz_total
       zz(jz) = real(self%geom%level(jz))
    end do
    call create_unstructured_grid(ug, nz_total, zz)

    n_vars = 1      !!!!! START WITH ONLY ONE VAR !!!!!!!!! 
    n_surf_vars = 0 !!!!! NO SURFACE VAR !!!!!!!!! 

    do jy=1,self%geom%ny
       do jx=1,self%geom%nx
          jk = 1
          cmask(jk) = int(self%geom%mask(jx,jy))       ! Surface T
          vv(jk) = self%tsfcn(jx,jy,cat_num)
          jk = jk + 1
          do jz = 1,self%geom%nzs                              ! Snow T
             cmask(jk) = int(self%geom%mask(jx,jy))    !
             !vv(jk) = Ts_nl(self%qsnon(jx,jy,cat_num))
             vv(jk) = self%qsnon(jx,jy,cat_num)             
             jk = jk + 1
          end do
          do jz = 1,self%geom%nzi                              ! Ice T
             cmask(jk) = int(self%geom%mask(jx,jy))    !
             !vv(jk) = Ti_nl(self%qicnk(jx,jy,cat_num,jz),self%sicnk(jx,jy,cat_num,jz))
             vv(jk) = self%qicnk(jx,jy,cat_num,jz)
             jk = jk + 1
          end do
          cmask(jk) = int(self%geom%mask(jx,jy))       ! Ice/Ocean interface
          !vv(jk) = Tm(self%sssoc(jx,jy))                 ! Tf = -mu * S
          vv(jk) = self%sssoc(jx,jy)                      ! Tf = -mu * S          
          jk = jk + 1          
          do jz = 1,self%geom%nzo                              ! Ocean
             cmask(jk) = int(self%geom%mask(jx,jy))    !             
             vv(jk) = self%sstoc(jx,jy)                   !
             jk = jk + 1
          end do

          !cmask(:) = int(self%geom%mask(jx,jy))           ! Some issues with the mask
          !print *,'cmask=',cmask
          !if (self%icemask(jx,jy)>0.0) then
          !print *,vv(:)
          !read(*,*)
          !end if
          !if (cmask(1)==1) read(*,*)

          call add_column(ug, self%geom%lat(jx,jy), self%geom%lon(jx,jy), self%geom%cell_area(jx, jy), &
               nz_total, &
               n_vars, &
               n_surf_vars, &
               cmask, &
               0)
          ug%last%column%fld3d(:) = vv(:)
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

    current => ug%head
    cat_num = 1 
    n_vars = 1      !!!!! START WITH ONLY ONE VAR !!!!!!!!! 
    n_surf_vars = 0 !!!!! NO SURFACE VAR !!!!!!!!! 

    do jy=1,self%geom%ny
       do jx=1,self%geom%nx
          jk = 1
          self%tsfcn(jx,jy,cat_num) = current%column%fld3d(jk)         ! Tsfcs
          jk = jk + 1
          do jz = 1,self%geom%nzs                                          ! Q Snow
             self%qsnon(jx,jy,cat_num) = current%column%fld3d(jk)
             jk = jk + 1
          end do
          do jz = 1,self%geom%nzi                                          ! Q Ice
             self%qicnk(jx,jy,cat_num,jz) = current%column%fld3d(jk)
             jk = jk + 1
          end do
          self%sssoc(jx,jy) = current%column%fld3d(jk)                 ! Ice/Ocean interface,
          ! Tf = -mu * S          
          jk = jk + 1          
          do jz = 1,self%geom%nzo                                          ! Ocean SST
             self%sstoc(jx,jy) = current%column%fld3d(jk)              !
             jk = jk + 1
          end do
          current => current%next
       enddo
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
    if (x1%geom%nzi /= x2%geom%nzi) call abor1_ftn("common_vars: error number of levels")
    common_vars = x1%geom%nzi * common_vars

  end function common_vars

  ! ------------------------------------------------------------------------------

  subroutine check_resolution(x1, x2)

    implicit none
    type(mom5cice5_field), intent(in) :: x1, x2

    ! NEEDS WORK !!!
    if (x1%geom%nx /= x2%geom%nx .or.  x1%geom%ny /= x2%geom%ny ) then
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
    bad = bad .or. (size(self%cicen, 1) /= self%geom%nx)

    ! add more test here ...

    if (bad) then
       write(0,*)'nx, ny, nf, nzi, nzo = ',self%geom%nx,self%geom%ny,self%nf,self%geom%nzi,self%geom%nzo
       call abor1_ftn ("mom5cice5_fields: field not consistent")
    endif

  end subroutine check

  ! ------------------------------------------------------------------------------

end module mom5cice5_fields
