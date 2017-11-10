
!> Fortran module handling interpolated (to obs locations) model variables

module mom5cice5_goms_mod

  use iso_c_binding
  use mom5cice5_geom_mod
  use mom5cice5_vars_mod
  use kinds
  !Interpolation related modules
  use type_linop
  use tools_interp, only: interp_horiz
  !use type_randgen, only: rng,initialize_sampling,create_randgen
  use module_namelist, only: namtype
  
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
     integer, allocatable :: numfld_per_fldname(:)        ! Cloned from field in interp
     real(kind=kind_real), allocatable :: values(:,:)  ! nvar x nobs
                                                       ! values(:,i)=[cicen(1:ncat),
                                                       !              hicen(1:ncat),
                                                       !              vicen(1:ncat),
                                                       !              hsnon(1:ncat),     
                                                       !              ...]
     character(len=5), allocatable :: variables(:)
     logical :: lalloc
     type(linoptype) :: hinterp_op
     logical :: hinterp_initialized !True:  hinterp_op has been initialized
                                    !False: hinterp_op not initialized
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

  subroutine gom_setup(self, vars, kobs)
    implicit none
    type(mom5cice5_goms), intent(inout) :: self
    type(mom5cice5_vars), intent(in) :: vars
    integer, intent(in) :: kobs(:)     ! Obs indices

    self%nobs=size(kobs)
    self%nvar=vars%nv
    self%used=1
    allocate(self%indx(self%nobs))
    self%indx(:)=kobs(:)
    allocate(self%variables(self%nvar))
    self%variables(1:self%nvar)=vars%fldnames(1:vars%nv)
    ! self%values is allocated after the initialization
    ! of the interpolation (in fields)
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

end module mom5cice5_goms_mod
