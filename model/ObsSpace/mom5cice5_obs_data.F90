
!> Handle observations for the MOM5CICE5 model

module mom5cice5_obs_data

  use iso_c_binding
  use string_f_c_mod
  use config_mod
  use datetime_mod
  use duration_mod
  use mom5cice5_goms_mod
  use mom5cice5_locs_mod
  use mom5cice5_obs_vectors
  use mom5cice5_obsoper_mod
  use mom5cice5_vars_mod
  use fckit_log_module, only : fckit_log
  use kinds

  implicit none
  private

  public :: obs_data, obs_setup, obs_delete, obs_get, obs_put, obs_count, max_string
  public :: obs_data_registry

  ! ------------------------------------------------------------------------------
  integer, parameter :: max_string=800
  ! ------------------------------------------------------------------------------

  !> A type to represent observation data
  type obs_data
     integer                   :: ngrp
     character(len=max_string) :: filein
     character(len=max_string) :: fileout     
     type(group_data), pointer :: grphead => null()
  end type obs_data

#define LISTED_TYPE obs_data

  !> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: obs_data_registry

  ! ------------------------------------------------------------------------------

  !> A type to represent a linked list of observation group data
  type group_data
     character(len=50)           :: grpname
     type(group_data), pointer   :: next => null()
     integer                     :: nobs
     integer, allocatable        :: seqnos(:)
     type(datetime), allocatable :: times(:)
     type(column_data), pointer  :: colhead => null()
  end type group_data

  ! ------------------------------------------------------------------------------

  !> A type to represent a linked list of observation columns
  type column_data
     character(len=50)                 :: colname
     type(column_data), pointer        :: next => null()
     integer                           :: ncol
     real(kind=kind_real), allocatable :: values(:,:)
  end type column_data

  ! ------------------------------------------------------------------------------

  !> Fortran generic
  interface obs_count
     module procedure obs_count_time, obs_count_all, obs_count_indx
  end interface obs_count

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "util/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine obs_setup(fin, fout, self)
    implicit none
    type(obs_data), intent(inout) :: self
    character(len=*), intent(in) :: fin, fout

    self%ngrp=0
    self%filein =fin
    self%fileout=fout

    if (self%filein/="") call obs_read(self)
    call fckit_log%debug("TRACE: mom5cice5_obs_data:obs_setup: done")

  end subroutine obs_setup

  ! ------------------------------------------------------------------------------

  subroutine obs_delete(self)
    implicit none
    type(obs_data), intent(inout) :: self
    type(group_data), pointer :: jgrp
    type(column_data), pointer :: jcol
    integer :: jo

    if (self%fileout/="") call obs_write(self)

    do while (associated(self%grphead))
       jgrp=>self%grphead
       self%grphead=>jgrp%next
       do jo=1,jgrp%nobs
          call datetime_delete(jgrp%times(jo))
       enddo
       deallocate(jgrp%times)
       do while (associated(jgrp%colhead))
          jcol=>jgrp%colhead
          jgrp%colhead=>jcol%next
          deallocate(jcol%values)
          deallocate(jcol)
       enddo
       deallocate(jgrp)
    enddo

  end subroutine obs_delete

  ! ------------------------------------------------------------------------------

  subroutine obs_get(self, req, col, ovec)
    implicit none
    type(obs_data), intent(in) :: self
    character(len=*), intent(in) :: req, col
    type(obs_vect), intent(inout) :: ovec

    type(group_data),  pointer :: jgrp
    type(column_data), pointer :: jcol
    integer :: jo, jc
    character(len=250) :: record

    write(record,*)"obs_get req=",req
    call fckit_log%info(record)
    write(record,*)"obs_get col=",col
    call fckit_log%info(record)

    ! Find obs group !!!! DEFINE OBS GROUP !!!!!!!!!
    call findgroup(self,req,jgrp)
    if (.not.associated(jgrp)) then
       jgrp=>self%grphead
       do while (associated(jgrp))
          write(record,*)"Group ",jgrp%grpname," exists."
          call fckit_log%info(record)
          jgrp=>jgrp%next
       enddo
       write(record,*)"Cannot find ",req," ."
       call fckit_log%error(record)
       call abor1_ftn("mom5cice5_obs_get: obs group not found")
    endif

    ! Find obs column
    call findcolumn(jgrp,col,jcol)
    if (.not.associated(jcol)) call abor1_ftn("mom5cice5_obs_get: obs column not found")

    ! Get data
    if (allocated(ovec%values)) deallocate(ovec%values)
    ovec%nobs=jgrp%nobs
    ovec%ncol=jcol%ncol
    allocate(ovec%values(ovec%ncol,ovec%nobs))

    do jo=1,jgrp%nobs
       do jc=1,jcol%ncol
          ovec%values(jc,jo)=jcol%values(jc,jo)
       enddo
    enddo

    write(record,*)"obs_get nobs, ncol=",jgrp%nobs,jcol%ncol
    call fckit_log%debug("TRACE: " // record)

  end subroutine obs_get

  ! ------------------------------------------------------------------------------

  subroutine obs_put(self, req, col, ovec)
    implicit none
    type(obs_data), intent(inout) :: self
    character(len=*), intent(in) :: req, col
    type(obs_vect), intent(in) :: ovec

    type(group_data),  pointer :: jgrp
    type(column_data), pointer :: jcol
    integer :: jo, jc
    character(len=250) :: record

    ! Find obs group
    call findgroup(self,req,jgrp)
    if (.not.associated(jgrp)) then
       jgrp=>self%grphead
       do while (associated(jgrp))
          write(record,*)"Group ",jgrp%grpname," exists."
          call fckit_log%info(record)
          jgrp=>jgrp%next
       enddo
       write(record,*)"Cannot find ",req," ."
       call fckit_log%error(record)
       call abor1_ftn("mom5cice5_obs_put: obs group not found")
    endif

    ! Find obs column (and add it if not there)
    call findcolumn(jgrp,col,jcol)
    if (.not.associated(jcol)) then
       if (.not.associated(jgrp%colhead)) call abor1_ftn("mom5cice5_obs_put: no locations")
       jcol=>jgrp%colhead
       do while (associated(jcol%next))
          jcol=>jcol%next
       enddo
       allocate(jcol%next)
       jcol=>jcol%next

       jcol%colname=col
       jcol%ncol=ovec%ncol
       allocate(jcol%values(jcol%ncol,jgrp%nobs))
    endif

    ! Put data
    if (ovec%nobs/=jgrp%nobs) call abor1_ftn("mom5cice5_obs_put: error obs number")
    if (ovec%ncol/=jcol%ncol) call abor1_ftn("mom5cice5_obs_put: error col number")
    do jo=1,jgrp%nobs
       do jc=1,jcol%ncol
          jcol%values(jc,jo)=ovec%values(jc,jo)
       enddo
    enddo

  end subroutine obs_put

  ! ------------------------------------------------------------------------------

  subroutine obs_locations(c_key_self, lreq, c_req, c_t1, c_t2, c_key_locs) bind(c,name='mom5cice5_obsdb_locations_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: lreq
    character(kind=c_char,len=1), intent(in) :: c_req(lreq+1)
    type(c_ptr), intent(in) :: c_t1, c_t2
    integer(c_int), intent(inout) :: c_key_locs

    type(obs_data), pointer :: self
    character(len=lreq) :: req
    type(datetime) :: t1, t2
    type(mom5cice5_locs), pointer :: locs
    type(obs_vect) :: ovec
    character(len=8) :: col="Location"

    call obs_data_registry%get(c_key_self, self)
    call c_f_string(c_req, req)
    call c_f_datetime(c_t1, t1)
    call c_f_datetime(c_t2, t2)

    call obs_time_get(self, req, col, t1, t2, ovec)

    call mom5cice5_locs_registry%init()
    call mom5cice5_locs_registry%add(c_key_locs)
    call mom5cice5_locs_registry%get(c_key_locs,locs)

    call mom5cice5_loc_setup(locs, ovec)

    deallocate(ovec%values)

  end subroutine obs_locations

  ! ------------------------------------------------------------------------------

  subroutine obs_getgom(c_key_self, lreq, c_req, c_key_vars, c_t1, c_t2, c_key_gom)&
       &bind(c,name='mom5cice5_obsdb_getgom_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: lreq
    character(kind=c_char,len=1), intent(in) :: c_req(lreq+1)
    integer(c_int), intent(in) :: c_key_vars
    type(c_ptr), intent(in) :: c_t1, c_t2
    integer(c_int), intent(inout) :: c_key_gom

    type(obs_data), pointer :: self
    character(len=lreq) :: req
    type(mom5cice5_vars), pointer :: vars
    type(datetime) :: t1, t2
    type(mom5cice5_goms), pointer :: gom

    integer :: nobs
    integer, allocatable :: mobs(:)

    character(len=21) :: t1str, t2str, tstr
    
    call obs_data_registry%get(c_key_self, self)
    call c_f_string(c_req, req)
    call mom5cice5_vars_registry%get(c_key_vars, vars)
    call c_f_datetime(c_t1, t1)
    call c_f_datetime(c_t2, t2)

    call obs_count(self, req, t1, t2, nobs)
    allocate(mobs(nobs))
    
    call obs_count(self, req, t1, t2, mobs)

    allocate(gom)
    call mom5cice5_goms_registry%init()
    call mom5cice5_goms_registry%add(c_key_gom)
    call mom5cice5_goms_registry%get(c_key_gom,gom)

    call gom_setup(gom, vars, mobs)
    deallocate(mobs)

  end subroutine obs_getgom

  ! ------------------------------------------------------------------------------

  subroutine obs_err_generate(c_key_self, c_key_type, perr) bind(c,name='mom5cice5_obsdb_seterr_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_type
    real(c_double), intent(in) :: perr

    type(obs_data), pointer :: self
    type(mom5cice5_obsoper), pointer :: otyp
    type(obs_vect) :: obserr
    integer :: nobs

    call obs_data_registry%get(c_key_self, self)
    call mom5cice5_obsoper_registry%get(c_key_type, otyp)

    call obs_count(self, otyp%request, nobs)
    call obsvec_setup(obserr,otyp%ncol,nobs)
    obserr%values(:,:)=perr
    call obs_put(self, trim(otyp%request), "ObsErr", obserr)
    deallocate(obserr%values)

  end subroutine obs_err_generate

  ! ------------------------------------------------------------------------------

  subroutine obs_generate(c_key_self, lreq, c_req, c_conf, c_bgn, c_step, ktimes, kobs) bind(c,name='mom5cice5_obsdb_generate_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: lreq
    character(kind=c_char,len=1), intent(in) :: c_req(lreq+1)
    type(c_ptr), intent(in)    :: c_conf
    type(c_ptr), intent(in)    :: c_bgn
    type(c_ptr), intent(in)    :: c_step
    integer(c_int), intent(in)  :: ktimes
    integer(c_int), intent(inout) :: kobs
    type(obs_data), pointer :: self
    character(len=lreq) :: req
    type(datetime) :: bgn
    type(duration) :: step

    integer :: nlocs
    type(datetime), allocatable :: times(:)
    type(obs_vect) :: obsloc

    character(len=21) :: tstr
    
    call obs_data_registry%get(c_key_self, self)
    call c_f_string(c_req, req)
    call c_f_datetime(c_bgn, bgn)
    call c_f_duration(c_step, step)

    nlocs  = config_get_int(c_conf, "obs_density");
    kobs=nlocs*ktimes;
    !call datetime_to_string(bgn, tstr)
    !print *,'=============== obs bgn time=',tstr    
    !print *,'=============== obs_density=',nlocs
    
    allocate(times(kobs))

    call generate_locations(c_conf, nlocs, ktimes, bgn, step, times, obsloc)
    call obs_create(self, trim(req), times, obsloc)

    deallocate(times)
    deallocate(obsloc%values)

  end subroutine obs_generate

  ! ------------------------------------------------------------------------------

  subroutine obs_nobs(c_key_self, lreq, c_req, kobs) bind(c,name='mom5cice5_obsdb_nobs_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: lreq
    character(kind=c_char,len=1), intent(in) :: c_req(lreq+1)
    integer(c_int), intent(inout) :: kobs

    type(obs_data), pointer :: self
    character(len=lreq) :: req
    integer :: iobs

    call obs_data_registry%get(c_key_self, self)
    call c_f_string(c_req, req)

    call obs_count(self, req, iobs)
    kobs = iobs

  end subroutine obs_nobs

  ! ------------------------------------------------------------------------------

  subroutine obs_time_get(self, req, col, t1, t2, ovec)
    implicit none
    type(obs_data), intent(in)    :: self
    character(len=*), intent(in)  :: req, col
    type(datetime), intent(in)    :: t1, t2
    type(obs_vect), intent(inout) :: ovec

    type(group_data),  pointer :: jgrp
    type(column_data), pointer :: jcol
    integer :: jo, jc, iobs

    character(len=21) :: t1str, t2str, tstr
    !character(kind=c_char,len=1) :: cstring(21)

    ! Find obs group
    call findgroup(self,req,jgrp)
    if (.not.associated(jgrp)) call abor1_ftn("obs_time_get: obs group not found")

    ! Find obs column
    call findcolumn(jgrp,col,jcol)
    if (.not.associated(jcol)) call abor1_ftn("obs_time_get: obs column not found")

    ! Time selection
    iobs=0
    do jo=1,jgrp%nobs
       if (t1<jgrp%times(jo) .and. jgrp%times(jo)<=t2) then
          iobs=iobs+1
          !print *,'iobs++!! = ',iobs
          !call datetime_to_string(t1, t1str)
          !call datetime_to_string(t2, t2str)
          !call datetime_to_string(jgrp%times(jo), tstr)
          !print *,'t1=',t1str,' t=',tstr,' t2=',t2str
          !read(*,*)
       end if
    enddo

    ! Get data
    if (ovec%nobs/=iobs .or. ovec%ncol/=jcol%ncol) then
       if (allocated(ovec%values)) deallocate(ovec%values)
       ovec%nobs=iobs
       ovec%ncol=jcol%ncol
       allocate(ovec%values(ovec%ncol,ovec%nobs))
    endif

    iobs=0
    do jo=1,jgrp%nobs
       if (t1<jgrp%times(jo) .and. jgrp%times(jo)<=t2) then
          iobs=iobs+1
          do jc=1,jcol%ncol
             ovec%values(jc,iobs)=jcol%values(jc,jo)
          enddo
       endif
    enddo
    
  end subroutine obs_time_get

  ! ------------------------------------------------------------------------------

  subroutine obs_count_time(self, req, t1, t2, kobs)
    implicit none
    type(obs_data), intent(in)   :: self
    character(len=*), intent(in) :: req
    type(datetime), intent(in)   :: t1, t2
    integer, intent(inout)       :: kobs

    type(group_data),  pointer :: jgrp
    integer :: jo

    ! Find obs group
    call findgroup(self,req,jgrp)
    if (.not.associated(jgrp)) call abor1_ftn("obs_count: obs group not found")

    ! Time selection
    kobs=0
    do jo=1,jgrp%nobs
       if (t1<jgrp%times(jo) .and. jgrp%times(jo)<=t2) kobs=kobs+1
    enddo

  end subroutine obs_count_time

  ! ------------------------------------------------------------------------------

  subroutine obs_count_indx(self, req, t1, t2, kobs)
    implicit none
    type(obs_data), intent(in)   :: self
    character(len=*), intent(in) :: req
    type(datetime), intent(in)   :: t1, t2
    integer, intent(inout)       :: kobs(:)

    type(group_data),  pointer :: jgrp
    integer :: jo, io

    ! Find obs group
    call findgroup(self,req,jgrp)
    if (.not.associated(jgrp)) call abor1_ftn("obs_count: obs group not found")

    ! Time selection
    io=0
    do jo=1,jgrp%nobs
       if (t1<jgrp%times(jo) .and. jgrp%times(jo)<=t2) then
          io=io+1
          kobs(io)=jo
       endif
    enddo

  end subroutine obs_count_indx

  ! ------------------------------------------------------------------------------

  subroutine obs_count_all(self, req, kobs)
    implicit none
    type(obs_data), intent(in) :: self
    character(len=*), intent(in) :: req
    integer, intent(inout) :: kobs
    type(group_data),  pointer :: jgrp

    call findgroup(self,req,jgrp)

    if (associated(jgrp)) then
       kobs=jgrp%nobs
    else
       kobs=0
    endif

  end subroutine obs_count_all

  ! ------------------------------------------------------------------------------

  subroutine obs_create(self, req, times, locs)
    implicit none
    type(obs_data), intent(inout) :: self
    character(len=*), intent(in) :: req
    type(datetime), intent(in) :: times(:)
    type(obs_vect), intent(in) :: locs
    type(group_data), pointer :: igrp
    integer :: jo, jc

    call findgroup(self,req,igrp)
    if (associated(igrp)) call abor1_ftn("obs_create: obs group already exists")

    if (associated(self%grphead)) then
       igrp=>self%grphead
       do while (associated(igrp%next))
          igrp=>igrp%next
       enddo
       allocate(igrp%next)
       igrp=>igrp%next
    else
       allocate(self%grphead)
       igrp=>self%grphead
    endif

    igrp%grpname=req
    igrp%nobs=size(times)
    allocate(igrp%times(igrp%nobs))
    igrp%times(:)=times(:)

    allocate(igrp%colhead)
    igrp%colhead%colname="Location"
    igrp%colhead%ncol=3
    allocate(igrp%colhead%values(3,igrp%nobs))
    if (locs%ncol/=3) call abor1_ftn("obs_create: error locations not 3D")
    if (locs%nobs/=igrp%nobs) call abor1_ftn("obs_create: error locations number")
    do jo=1,igrp%nobs
       do jc=1,3
          igrp%colhead%values(jc,jo)=locs%values(jc,jo)
       enddo
    enddo

    self%ngrp=self%ngrp+1

  end subroutine obs_create

  ! ------------------------------------------------------------------------------
  !  Private
  ! ------------------------------------------------------------------------------

  subroutine obs_read(self)
    implicit none
    type(obs_data), intent(inout) :: self
    integer :: iin, icol, jo, jc, jg, ncol
    type(group_data), pointer :: jgrp
    type(column_data), pointer :: jcol
    real(kind=kind_real), allocatable :: ztmp(:)
    character(len=20) :: stime
    character(len=max_string+50) :: record

    character(len=21) :: tstr

    iin=90
    write(record,*)'obs_read: opening ',trim(self%filein)
    call fckit_log%info(record)
    open(unit=iin, file=trim(self%filein), form='formatted', action='read')

    read(iin,*)self%ngrp
    do jg=1,self%ngrp
       if (jg==1) then
          allocate(self%grphead)
          jgrp=>self%grphead
       else
          allocate(jgrp%next)
          jgrp=>jgrp%next
       endif
       read(iin,*)jgrp%grpname
       read(iin,*)jgrp%nobs
       write(record,*)'obs_read: reading ',jgrp%nobs,' ',jgrp%grpname,' observations.'
       call fckit_log%info(record)
       allocate(jgrp%times(jgrp%nobs))

       read(iin,*)ncol
       icol=0
       do jc=1,ncol
          if (jc==1) then
             allocate(jgrp%colhead)
             jcol=>jgrp%colhead
          else
             allocate(jcol%next)
             jcol=>jcol%next
          endif
          read(iin,*)jcol%colname, jcol%ncol
          icol=icol+jcol%ncol
          allocate(jcol%values(jcol%ncol,jgrp%nobs))
       enddo

       allocate(ztmp(icol))
       do jo=1,jgrp%nobs
          read(iin,*)stime,ztmp(:)
          call datetime_create(stime,jgrp%times(jo))
          call datetime_to_string(jgrp%times(jo), tstr)
          icol=0
          jcol=>jgrp%colhead
          do while (associated(jcol))
             do jc=1,jcol%ncol
                icol=icol+1
                jcol%values(jc,jo)=ztmp(icol)
             enddo
             jcol=>jcol%next
          enddo
       enddo
       deallocate(ztmp)
    enddo

    close(iin)

  end subroutine obs_read

  ! ------------------------------------------------------------------------------

  subroutine obs_write(self)
    implicit none
    type(obs_data), intent(in) :: self
    integer :: iout, icol, jc, jo
    type(group_data), pointer :: jgrp
    type(column_data), pointer :: jcol
    real(kind=kind_real), allocatable :: ztmp(:)
    character(len=20) :: stime

    iout=91
    open(unit=iout, file=trim(self%fileout), form='formatted', action='write')    
    write(iout,*)self%ngrp
    jgrp=>self%grphead
    do while (associated(jgrp))
       write(iout,*)jgrp%grpname
       write(iout,*)jgrp%nobs

       icol=0
       jcol=>jgrp%colhead
       do while (associated(jcol))
          icol=icol+1
          jcol=>jcol%next
       enddo
       write(iout,*)icol

       icol=0
       jcol=>jgrp%colhead
       do while (associated(jcol))
          write(iout,*)jcol%colname, jcol%ncol
          icol=icol+jcol%ncol
          jcol=>jcol%next
       enddo

       allocate(ztmp(icol))
       do jo=1,jgrp%nobs
          icol=0
          jcol=>jgrp%colhead
          do while (associated(jcol))
             do jc=1,jcol%ncol
                icol=icol+1
                ztmp(icol)=jcol%values(jc,jo)
             enddo
             jcol=>jcol%next
          enddo
          call datetime_to_string(jgrp%times(jo),stime)
          write(iout,*)stime,ztmp(:)
       enddo
       deallocate(ztmp)

       jgrp=>jgrp%next
    enddo

    close(iout)

  end subroutine obs_write

  ! ------------------------------------------------------------------------------

  subroutine findgroup(self,req,find)
    type(obs_data), intent(in) :: self
    character(len=*), intent(in) :: req
    type(group_data), pointer, intent(inout) :: find

    find=>self%grphead
    do while (associated(find))
       if (find%grpname==req) exit
       find=>find%next
    enddo

  end subroutine findgroup

  ! ------------------------------------------------------------------------------

  subroutine findcolumn(grp,col,find)
    type(group_data), intent(in) :: grp
    character(len=*), intent(in) :: col
    type(column_data), pointer, intent(inout) :: find

    find=>grp%colhead
    do while (associated(find))
       if (find%colname==col) exit
       find=>find%next
    enddo

  end subroutine findcolumn

  ! ------------------------------------------------------------------------------

end module mom5cice5_obs_data
