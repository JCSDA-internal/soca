
!> Fortran module handling interpolated (to obs locations) model variables

module mom5cice5_goms_mod

  use iso_c_binding
  use mom5cice5_geom_mod
  use mom5cice5_vars_mod
  use kinds

  implicit none
  private
  public :: mom5cice5_goms, gom_setup
  public :: mom5cice5_goms_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to hold interpolated fields required by the obs operators
  type :: mom5cice5_goms
     !type(mom5cice5_geom), pointer :: geom !< MOM5 & CICE5 Geometry     
     integer :: nobs                                   ! Number of obs (time and loc) in DA window (?)
     integer :: nvar                                   ! Number of variables in gom
     integer :: used
     integer, allocatable :: indx(:)
     integer, allocatable :: numfld_per_fldname(:)     ! Cloned from field in interp
     real(kind=kind_real), allocatable :: values(:,:)  ! nvar x nobs
     ! values(:,i)=[cicen(1:ncat),
     !              hicen(1:ncat),
     !              vicen(1:ncat),
     !              hsnon(1:ncat),     
     !              ...]
     character(len=5), allocatable :: variables(:)
     logical :: lalloc
  end type mom5cice5_goms

#define LISTED_TYPE mom5cice5_goms

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: mom5cice5_goms_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_gom_create(c_key_self) bind(c,name='mom5cice5_gom_create_f90')

    implicit none
    integer(c_int), intent(inout) :: c_key_self

    type(mom5cice5_goms), pointer :: self

    print *,'&&&&&&&&&&&&&&&&&&&&&&& IN GOM CREATE &&&&&&&&&&&&&&&&&&&'
    
    call mom5cice5_goms_registry%init()
    call mom5cice5_goms_registry%add(c_key_self)
    call mom5cice5_goms_registry%get(c_key_self, self)

    self%lalloc = .false.

  end subroutine c_mom5cice5_gom_create

  ! ------------------------------------------------------------------------------

  subroutine gom_setup(self, vars, kobs)
    implicit none
    type(mom5cice5_goms), intent(inout) :: self
    type(mom5cice5_vars), intent(in) :: vars
    integer, intent(in) :: kobs(:)     ! Obs indices

    print *,'@@@@@@@@@@@@@@@@@@@@@ IN GOM SETUP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    
    self%nobs=size(kobs)
    self%nvar=vars%nv
    self%used=1
    allocate(self%indx(self%nobs))
    self%indx(:)=kobs(:)
    allocate(self%variables(self%nvar))
    self%variables(1:self%nvar)=vars%fldnames(1:vars%nv)
    
    if (.not.allocated(self%numfld_per_fldname)) then
       print *,'in gom setup'
       print *,vars%fldnames
       read(*,*)
       allocate(self%numfld_per_fldname(vars%nv))  
    end if
    self%lalloc = .true.

  end subroutine gom_setup

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_gom_delete(c_key_self) bind(c,name='mom5cice5_gom_delete_f90')

    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(mom5cice5_goms), pointer :: self

    call mom5cice5_goms_registry%get(c_key_self, self)
    if (self%lalloc) then
       if (allocated(self%values)) deallocate(self%values) !Not allocated in constructor
       deallocate(self%indx)
       deallocate(self%variables)
    endif
    call mom5cice5_goms_registry%remove(c_key_self)

  end subroutine c_mom5cice5_gom_delete

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_gom_zero(c_key_self) bind(c,name='mom5cice5_gom_zero_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    type(mom5cice5_goms), pointer :: self
    call mom5cice5_goms_registry%get(c_key_self, self)
    self%values(:,:)=0.0_kind_real
  end subroutine c_mom5cice5_gom_zero

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_gom_random(c_key_self) bind(c,name='mom5cice5_gom_random_f90')
    use random_vectors_mod
    implicit none
    integer(c_int), intent(in) :: c_key_self
    type(mom5cice5_goms), pointer :: self
    call mom5cice5_goms_registry%get(c_key_self, self)
    call random_vector(self%values(:,:))
  end subroutine c_mom5cice5_gom_random

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_gom_mult(c_key_self, zz) bind(c,name='mom5cice5_gom_mult_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    real(c_double), intent(in) :: zz
    type(mom5cice5_goms), pointer :: self
    integer :: jo, jv

    call mom5cice5_goms_registry%get(c_key_self, self)
    do jo=1,self%nobs
       do jv=1,self%nvar
          self%values(jv,jo) = zz * self%values(jv,jo)
       enddo
    enddo

  end subroutine c_mom5cice5_gom_mult

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_gom_dotprod(c_key_self, c_key_other, prod) bind(c,name='mom5cice5_gom_dotprod_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self, c_key_other
    real(c_double), intent(inout) :: prod
    type(mom5cice5_goms), pointer :: self, other
    integer :: jo, jv

    call mom5cice5_goms_registry%get(c_key_self, self)
    call mom5cice5_goms_registry%get(c_key_other, other)
    prod=0.0_kind_real
    do jo=1,self%nobs
       do jv=1,self%nvar
          prod=prod+self%values(jv,jo)*other%values(jv,jo)
       enddo
    enddo

  end subroutine c_mom5cice5_gom_dotprod

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_gom_minmaxavg(c_key_self, kobs, pmin, pmax, prms) bind(c,name='mom5cice5_gom_minmaxavg_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(inout) :: kobs
    real(c_double), intent(inout) :: pmin, pmax, prms
    type(mom5cice5_goms), pointer :: self

    call mom5cice5_goms_registry%get(c_key_self, self)

    kobs = self%nobs
    pmin=minval(self%values(:,:))
    pmax=maxval(self%values(:,:))
    prms=sqrt(sum(self%values(:,:)**2)/real(self%nobs*self%nvar,kind_real))

  end subroutine c_mom5cice5_gom_minmaxavg

  ! ------------------------------------------------------------------------------
  subroutine mom5cice5_gom_read_file_c(c_key_self, c_conf) bind(c,name='mom5cice5_gom_read_file_f90')
    !use nc_diag_write_mod!, only: nc_diag_init
    use nc_diag_write_mod, only: nc_diag_init, nc_diag_metadata, nc_diag_write    
    use config_mod
    use fckit_log_module, only : fckit_log
    implicit none
    integer(c_int), intent(in) :: c_key_self
    type(c_ptr), intent(in)    :: c_conf
    type(mom5cice5_goms), pointer :: self

    integer, parameter :: iunit=10
    integer, parameter :: max_string_length=250 ! Yuk!
    character(len=max_string_length) :: filename, record
    character(len=4)  :: cnx
    character(len=17) :: fmtn
    character(len=11) :: fmt1='(X,ES24.16)'
    integer :: jj, jo, jv

    print *,'@@@@@@@@@@@@@@@@@@@ IN GOM READ @@@@@@@@@@@@@@@@@@@@@@@@@@' 

    call mom5cice5_goms_registry%get(c_key_self, self)
    if (self%lalloc) call abor1_ftn("mom5cice5_gom_read_file gom alredy allocated")

    filename = config_get_string(c_conf,len(filename),"filename")
    write(record,*)'mom5cice5_gom_read_file: opening '//trim(filename)
    call fckit_log%info(record)
    open(unit=iunit, file=trim(filename), form='formatted', action='read')

    ! open netcdf file and read dimensions
    !call nc_diag_init(filename)!, iunit)
    !call write_split_conv_diag_nc(fn, hdr, mass, wind, append_suffix)
    !call nc_diag_metadata("Station_ID",              conv_wind(i)%Station_ID)x    

    read(iunit,*) self%nobs, self%nvar, self%used
    allocate(self%indx(self%nobs))
    allocate(self%variables(self%nvar))
    allocate(self%values(self%nvar,self%nobs))

    read(iunit,*) self%indx(:)
    do jv=1,self%nvar
       read(iunit,*) self%variables(jv)
    enddo

    if (self%nvar>9999)  call abor1_ftn("Format too small")
    write(cnx,'(I4)')self%nvar
    fmtn='('//trim(cnx)//fmt1//')'

    do jo=1,self%nobs
       read(iunit,fmtn) (self%values(jj,jo), jj=1,self%nvar)
    enddo

    close(iunit)
!!$    !call nc_diag_write_close(filename)
    self%lalloc = .true.

  end subroutine mom5cice5_gom_read_file_c

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_gom_write_file_c(c_key_self, c_conf) bind(c,name='mom5cice5_gom_write_file_f90')
    use nc_diag_write_mod!, only: nc_diag_init 
    use config_mod
    use fckit_log_module, only : fckit_log
    implicit none
    integer(c_int), intent(in) :: c_key_self
    type(c_ptr), intent(in) :: c_conf
    type(mom5cice5_goms), pointer :: self

    integer, parameter :: iunit=10
    integer, parameter :: max_string_length=250 ! Yuk!
    character(len=max_string_length) :: filename, record
    character(len=4)  :: cnx
    character(len=17) :: fmtn
    character(len=11) :: fmt1='(X,ES24.16)'
    integer :: jj, jo, jv

!!!!!!!!     FOR EXAMLE: LOOK AT gsidiag_conv_bin2nc4.f90 IN IODA REPO  !!!!!!!!!!!!!

    call mom5cice5_goms_registry%get(c_key_self, self)
    if (.not.self%lalloc) call abor1_ftn("mom5cice5_gom_write_file gom not allocated")

    filename = config_get_string(c_conf,len(filename),"filename")
    write(record,*)'mom5cice5_gom_write_file: opening '//trim(filename)
    call fckit_log%info(record)
    !call nc_diag_init(filename)!, iunit)

    open(unit=iunit, file=trim(filename), form='formatted', action='write')

    write(iunit,*) self%nobs, self%nvar, self%used
    write(iunit,*) self%indx(:)
    do jv=1,self%nvar
       write(iunit,*) self%variables(jv)
    enddo

    if (self%nvar>9999) call abor1_ftn("Format too small")
    write(cnx,'(I4)')self%nvar
    fmtn='('//trim(cnx)//fmt1//')'

    do jo=1,self%nobs
       write(iunit,fmtn) (self%values(jj,jo), jj=1,self%nvar)
    enddo

    close(iunit)

  end subroutine mom5cice5_gom_write_file_c

  ! ------------------------------------------------------------------------------  

end module mom5cice5_goms_mod
