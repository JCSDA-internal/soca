
!> Handle fields for the QG model

module mom5cice5_fields

  use config_mod
  use mom5cice5_geom_mod
  use mom5cice5_vars_mod
  use kinds

  implicit none
  private

  public :: mom5cice5_field, &
       & create, delete, zeros, random, copy, &
       & self_add, self_schur, self_sub, self_mul, axpy, &
       & dot_prod, add_incr, diff_incr, &
       & read_file, write_file, gpnorm, fldrms, &
       & change_resol
  public :: mom5cice5_field_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to hold QG fields
  type :: mom5cice5_field
     integer :: nx                     !< Zonal grid dimension
     integer :: ny                     !< Meridional grid dimension
     integer :: nzo                    !< Number of z levels in the ocean
     integer :: nzi                    !< Number of levels in sea-ice  
     integer :: ncat                   !< Number of sea-ice thickness categories    
     integer :: nf                     !< Number of fields

     real(kind=kind_real), allocatable :: cn(:,:,:)        !< Sea-ice fraction
     real(kind=kind_real), allocatable :: hicen(:,:,:)     !< Sea-ice thickness
     real(kind=kind_real), allocatable :: vicen(:,:,:)     !< Sea-ice volume
     real(kind=kind_real), allocatable :: hsnon(:,:,:)     !< Snow depth over sea-ice
     real(kind=kind_real), allocatable :: vsnon(:,:,:)     !< Snow volume over sea-ice  
     real(kind=kind_real), allocatable :: tsfcn(:,:,:)     !< Temperature over sea-ice or snow
     real(kind=kind_real), allocatable :: qsnon(:,:,:)     !< Enthalpy of snow
     real(kind=kind_real), allocatable :: sicenk(:,:,:,:)  !< Salinity of sea-ice
     real(kind=kind_real), allocatable :: so(:,:)          !< Ocean (surface) Salinity
     real(kind=kind_real), allocatable :: qicenk(:,:,:,:)  !< Enthalpy of sea-ice
     real(kind=kind_real), allocatable :: to(:,:)          !< Liquid ocean temperature
     real(kind=kind_real), allocatable :: tsst(:,:)        !< Average temperature 
     character(len=1), allocatable :: fldnames(:)      !< Variable identifiers
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
    integer :: ioff

    self%nx   = geom%nx
    self%ny   = geom%ny
    self%nzo  = geom%nzo
    self%nzi  = geom%nzi
    self%ncat = geom%ncat
    self%nf   = vars%nv

    allocate(self%cn(self%nx,self%ny,self%ncat))
    allocate(self%hicen(self%nx,self%ny,self%ncat))
    allocate(self%vicen(self%nx,self%ny,self%ncat))
    allocate(self%hsnon(self%nx,self%ny,self%ncat))
    allocate(self%vsnon(self%nx,self%ny,self%ncat))
    allocate(self%tsfcn(self%nx,self%ny,self%ncat))
    allocate(self%qsnon(self%nx,self%ny,self%ncat))
    allocate(self%sicenk(self%nx,self%ny,self%ncat,self%nzi))
    allocate(self%so(self%nx,self%ny))
    allocate(self%qicenk(self%nx,self%ny,self%ncat,self%nzi))
    allocate(self%to(self%nx,self%ny))
    allocate(self%tsst(self%nx,self%ny))    

    self%cn(:,:,:)=0.0_kind_real
    self%hicen(:,:,:)=0.0_kind_real
    self%vicen(:,:,:)=0.0_kind_real        
    self%hsnon(:,:,:)=0.0_kind_real
    self%vsnon(:,:,:)=0.0_kind_real
    self%tsfcn(:,:,:)=0.0_kind_real
    self%qsnon(:,:,:)=0.0_kind_real
    self%sicenk(:,:,:,:)=0.0_kind_real
    self%so(:,:)=0.0_kind_real    
    self%qicenk(:,:,:,:)=0.0_kind_real
    self%to(:,:)=0.0_kind_real    
    self%tsst(:,:)=0.0_kind_real

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

    if (allocated(self%cn)) deallocate(self%cn)
    if (allocated(self%hicen)) deallocate(self%hicen)
    if (allocated(self%vicen)) deallocate(self%vicen)        
    if (allocated(self%hsnon)) deallocate(self%hsnon)
    if (allocated(self%vsnon)) deallocate(self%vsnon)
    if (allocated(self%tsfcn)) deallocate(self%tsfcn)    
    if (allocated(self%qsnon)) deallocate(self%qsnon)
    if (allocated(self%sicenk)) deallocate(self%sicenk)
    if (allocated(self%so)) deallocate(self%so)
    if (allocated(self%qicenk)) deallocate(self%qicenk)
    if (allocated(self%to)) deallocate(self%to)
    if (allocated(self%tsst)) deallocate(self%tsst)    
    if (allocated(self%fldnames)) deallocate(self%fldnames)

  end subroutine delete

  ! ------------------------------------------------------------------------------

  subroutine zeros(self)
    implicit none
    type(mom5cice5_field), intent(inout) :: self

    call check(self)

    self%cn=0.0_kind_real
    self%hicen=0.0_kind_real
    self%vicen=0.0_kind_real        
    self%hsnon=0.0_kind_real
    self%vsnon=0.0_kind_real
    self%tsfcn=0.0_kind_real
    self%qsnon=0.0_kind_real
    self%sicenk=0.0_kind_real
    self%so=0.0_kind_real    
    self%qicenk=0.0_kind_real
    self%to=0.0_kind_real    
    self%tsst=0.0_kind_real

  end subroutine zeros

  ! ------------------------------------------------------------------------------

  subroutine random(self)
    use random_vectors_gauss_mod
    implicit none
    type(mom5cice5_field), intent(inout) :: self

    call check(self)

    call random_vector_gauss(self%gfld3d(:,:,:))

  end subroutine random

  ! ------------------------------------------------------------------------------

  subroutine copy(self,rhs)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    type(mom5cice5_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cn=rhs%cn
    self%hicen=rhs%hicen
    self%vicen=rhs%vicen
    self%hsnon=rhs%hsnon
    self%vsnon=rhs%vsnon
    self%tsfcn=rhs%tsfcn
    self%qsnon=rhs%qsnon
    self%sicenk=rhs%sicenk
    self%so=rhs%so
    self%qicenk=rhs%qicenk
    self%to=rhs%to
    self%tsst=rhs%tsst

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

    self%cn=self%cn+rhs%cn
    self%hicen=self%hicen+rhs%hicen
    self%vicen=self%vicen+rhs%vicen
    self%hsnon=self%hsnon+rhs%hsnon
    self%vsnon=self%vsnon+rhs%vsnon
    self%tsfcn=self%tsfcn+rhs%tsfcn
    self%qsnon=self%qsnon+rhs%qsnon
    self%sicenk=self%sicenk+rhs%sicenk
    self%so=self%so+rhs%so
    self%qicenk=self%qicenk+rhs%qicenk
    self%to=self%to+rhs%to
    self%tsst=self%tsst+rhs%tsst

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

    self%cn=self%cn*rhs%cn
    self%hicen=self%hicen*rhs%hicen
    self%vicen=self%vicen*rhs%vicen
    self%hsnon=self%hsnon*rhs%hsnon
    self%vsnon=self%vsnon*rhs%vsnon
    self%tsfcn=self%tsfcn*rhs%tsfcn
    self%qsnon=self%qsnon*rhs%qsnon
    self%sicenk=self%sicenk*rhs%sicenk
    self%so=self%so*rhs%so
    self%qicenk=self%qicenk*rhs%qicenk
    self%to=self%to*rhs%to
    self%tsst=self%tsst*rhs%tsst

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

    self%cn=self%cn-rhs%cn
    self%hicen=self%hicen-rhs%hicen
    self%vicen=self%vicen-rhs%vicen
    self%hsnon=self%hsnon-rhs%hsnon
    self%vsnon=self%vsnon-rhs%vsnon
    self%tsfcn=self%tsfcn-rhs%tsfcn
    self%qsnon=self%qsnon-rhs%qsnon
    self%sicenk=self%sicenk-rhs%sicenk
    self%so=self%so-rhs%so
    self%qicenk=self%qicenk-rhs%qicenk
    self%to=self%to-rhs%to
    self%tsst=self%tsst-rhs%tsst

    return
  end subroutine self_sub

  ! ------------------------------------------------------------------------------

  subroutine self_mul(self,zz)
    implicit none
    type(mom5cice5_field), intent(inout) :: self
    real(kind=kind_real), intent(in) :: zz

    call check(self)

    self%cn = zz * self%cn
    self%hicen = zz * self%hicen
    self%vicen = zz * self%vicen
    self%hsnon = zz * self%hsnon
    self%vsnon = zz * self%vsnon
    self%tsfcn = zz * self%tsfcn
    self%qsnon = zz * self%qsnon
    self%sicenk = zz * self%sicenk
    self%so = zz * self%so
    self%qicenk = zz * self%qicenk
    self%to = zz * self%to
    self%tsst = zz * self%tsst

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

    self%cn=self%cn + zz * rhs%cn
    self%hicen=self%hicen + zz * rhs%hicen
    self%vicen=self%vicen + zz * rhs%vicen
    self%hsnon=self%hsnon + zz * rhs%hsnon
    self%vsnon=self%vsnon + zz * rhs%vsnon
    self%tsfcn=self%tsfcn + zz * rhs%tsfcn
    self%qsnon=self%qsnon + zz * rhs%qsnon
    self%sicenk=self%sicenk + zz * rhs%sicenk
    self%so=self%so + zz * rhs%so
    self%qicenk=self%qicenk + zz * rhs%qicenk
    self%to=self%to + zz * rhs%to
    self%tsst=self%tsst + zz * rhs%tsst        

    return
  end subroutine axpy

  ! ------------------------------------------------------------------------------

  subroutine dot_prod(fld1,fld2,zprod)
    implicit none
    type(mom5cice5_field), intent(in) :: fld1, fld2

    call check_resolution(fld1, fld2)
    if (fld1%nf /= fld2%nf .or. fld1%nl /= fld2%nl) then
       call abor1_ftn("mom5cice5_fields:field_prod error number of fields")
    endif
    call abor1_ftn("mom5cice5_fields:field_prod should never dot_product full state")    

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
          lhs%cn = x1%cn - x2%cn
          lhs%hicen = x1%hicen - x2%hicen
          lhs%vicen = x1%vicen - x2%vicen
          lhs%hsnon = x1%hsnon - x2%hsnon
          lhs%vsnon = x1%vsnon - x2%vsnon
          lhs%tsfcn = x1%tsfcn - x2%tsfcn
          lhs%qsnon = x1%qsnon - x2%qsnon
          lhs%sicenk = x1%sicenk -x2%sicenk
          lhs%so = x1%so - x2%so
          lhs%qicenk = x1%qicenk - x2%qicenk
          lhs%to = x1%to - x2%to
          lhs%tsst = x1%tsst - x2%tsst         
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

    call abor1_ftn("mom5cice5_fields:field_resol: untested code")

    return
  end subroutine change_resol

  ! ------------------------------------------------------------------------------

  subroutine read_file(fld, c_conf, vdate)
    ! Needs more interface clean-up here...
    use iso_c_binding
    use datetime_mod
    use fckit_log_module, only : log

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
    integer :: ic, iy, il, ix, is, jx, jy, jf, iread, nf
    real(kind=kind_real), allocatable :: zz(:)

    iread = 1
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
       call zeros(fld)
       filename = config_get_string(c_conf,len(filename),"filename")
       WRITE(buf,*) 'mom5cice5_field:read_file: opening '//filename
       call log%info(buf)
       open(unit=iunit, file=trim(filename), form='formatted', action='read')

       read(iunit,*) ix, iy, il, ic, is
       if (ix /= fld%nx .or. iy /= fld%ny .or. il /= fld%nl) then
          write (record,*) "mom5cice5_fields:read_file: ", &
               & "input fields have wrong dimensions: ",ix,iy,il
          call log%error(record)
          write (record,*) "mom5cice5_fields:read_file: expected: ",fld%nx,fld%ny,fld%nl
          call log%error(record)
          call abor1_ftn("mom5cice5_fields:read_file: input fields have wrong dimensions")
       endif

       read(iunit,*) sdate
       WRITE(buf,*) 'validity date is: '//sdate
       call log%info(buf)
       call datetime_set(sdate, vdate)

       if (fld%nx>9999)  call abor1_ftn("Format too small")
       write(cnx,'(I4)')fld%nx
       fmtn='('//trim(cnx)//fmt1//')'

       nf = min(fld%nf, ic)
       do jf=1,il*nf
          do jy=1,fld%ny
             read(iunit,fmtn) (fld%gfld3d(jx,jy,jf), jx=1,fld%nx)
          enddo
       enddo
       ! Skip un-necessary data from file if any
       allocate(zz(fld%nx))
       do jf=nf*il+1, ic*il
          do jy=1,fld%ny
             read(iunit,fmtn) (zz(jx), jx=1,fld%nx)
          enddo
       enddo
       deallocate(zz)

       if (fld%lbc) then
          do jf=1,4
             read(iunit,fmt1) fld%xbound(jf)
          enddo
          do jf=1,4
             read(iunit,fmtn) (fld%qbound(jx,jf), jx=1,fld%nx)
          enddo
       endif

       close(iunit)
    endif

    call check(fld)

    return
  end subroutine read_file

  ! ------------------------------------------------------------------------------

  subroutine write_file(fld, c_conf, vdate)
    use iso_c_binding
    use datetime_mod
    use fckit_log_module, only : log

    implicit none
    type(mom5cice5_field), intent(in) :: fld    !< Fields
    type(c_ptr), intent(in)    :: c_conf !< Configuration
    type(datetime), intent(in) :: vdate  !< DateTime

    integer, parameter :: iunit=11
    integer, parameter :: max_string_length=800 ! Yuk!
    character(len=max_string_length+50) :: record
    character(len=max_string_length) :: filename
    character(len=20) :: sdate, fmtn
    character(len=4)  :: cnx
    character(len=11) :: fmt1='(X,ES24.16)'
    character(len=1024):: buf
    integer :: jf, jy, jx, is

    call check(fld)

    filename = genfilename(c_conf,max_string_length,vdate)
    WRITE(buf,*) 'mom5cice5_field:write_file: writing '//filename
    call log%info(buf)
    open(unit=iunit, file=trim(filename), form='formatted', action='write')

    is=0
    if (fld%lbc) is=1

    write(iunit,*) fld%nx, fld%ny, fld%nl, fld%nf, is

    call datetime_to_string(vdate, sdate)
    write(iunit,*) sdate

    if (fld%nx>9999)  call abor1_ftn("Format too small")
    write(cnx,'(I4)')fld%nx
    fmtn='('//trim(cnx)//fmt1//')'

    do jf=1,fld%nl*fld%nf
       do jy=1,fld%ny
          write(iunit,fmtn) (fld%gfld3d(jx,jy,jf), jx=1,fld%nx)
       enddo
    enddo

    if (fld%lbc) then
       do jf=1,4
          write(iunit,fmt1) fld%xbound(jf)
       enddo
       do jf=1,4
          write(iunit,fmtn) (fld%qbound(jx,jf), jx=1,fld%nx)
       enddo
    endif

    close(iunit)

    return
  end subroutine write_file

  ! ------------------------------------------------------------------------------

  subroutine gpnorm(fld, nf, pstat)
    implicit none
    type(mom5cice5_field), intent(in) :: fld
    integer, intent(in) :: nf
    real(kind=kind_real), intent(inout) :: pstat(3, nf)
    integer :: jj,joff

    call check(fld)

    do jj=1,fld%nf
       joff=(jj-1)*fld%nl
       pstat(1,jj)=minval(fld%gfld3d(:,:,joff+1:joff+fld%nl))
       pstat(2,jj)=maxval(fld%gfld3d(:,:,joff+1:joff+fld%nl))
       pstat(3,jj)=sqrt(sum(fld%gfld3d(:,:,joff+1:joff+fld%nl)**2) &
            & /real(fld%nl*fld%nx*fld%ny,kind_real))
    enddo
    jj=jj-1

    if (fld%lbc) then
       jj=jj+1
       pstat(1,jj)=minval(fld%xbound(:))
       pstat(2,jj)=maxval(fld%xbound(:))
       pstat(3,jj)=sqrt(sum(fld%xbound(:)**2)/real(4,kind_real))

       jj=jj+1
       pstat(1,jj)=minval(fld%qbound(:,:))
       pstat(2,jj)=maxval(fld%qbound(:,:))
       pstat(3,jj)=sqrt(sum(fld%qbound(:,:)**2)/real(4*fld%nx,kind_real))
    endif

    if (jj /= nf) call abor1_ftn("mom5cice5_fields_gpnorm: error number of fields")

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

    zz = 0.0

    do jf=1,fld%nl*fld%nf
       do jy=1,fld%ny
          do jx=1,fld%nx
             zz = zz + fld%gfld3d(jx,jy,jf)*fld%gfld3d(jx,jy,jf)
          enddo
       enddo
    enddo

    ii = fld%nl*fld%nf*fld%ny*fld%nx

    prms = sqrt(zz/real(ii,kind_real))

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
         & call abor1_ftn("mom5cice5_fields:genfilename: filename too long")

  end function genfilename

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
    if (x1%nl /= x2%nl) call abor1_ftn("common_vars: error number of levels")
    common_vars = x1%nl * common_vars

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
    bad = bad .or. (size(self%cn, 1) /= self%nx)

    ! add more test here ...

    if (bad) then
       write(0,*)'nx, ny, nf, nl, lbc = ',self%nx,self%ny,self%nf,self%nl,self%lbc
       call abor1_ftn ("mom5cice5_fields: field not consistent")
    endif

  end subroutine check

  ! ------------------------------------------------------------------------------

end module mom5cice5_fields
