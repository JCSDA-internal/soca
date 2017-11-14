
!> Fortran module handling observation vectors

module mom5cice5_obs_vectors

  use iso_c_binding
  use random_vectors_mod
  use kinds

  implicit none
  private
  public :: obs_vect, obsvec_setup
  public :: mom5cice5_obs_vect_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to represent an observation vector
  type obs_vect
     integer :: nobs=0
     integer :: ncol=0
     real(kind=kind_real), allocatable :: values(:,:)  ! ncol x nobs
  end type obs_vect

#define LISTED_TYPE obs_vect

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: mom5cice5_obs_vect_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_obsvec_setup(c_key_self, ncol, nobs) bind(c,name='mom5cice5_obsvec_setup_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    integer(c_int), intent(in) :: ncol, nobs
    type(obs_vect), pointer :: self

    call mom5cice5_obs_vect_registry%init()
    call mom5cice5_obs_vect_registry%add(c_key_self)
    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    call obsvec_setup(self, ncol, nobs)
    
  end subroutine c_mom5cice5_obsvec_setup

  ! ------------------------------------------------------------------------------

  subroutine obsvec_setup(self, nc, no)
    implicit none
    type(obs_vect), intent(inout) :: self
    integer(c_int), intent(in)    :: nc, no

    self%ncol=nc
    self%nobs=no
    if (allocated(self%values)) deallocate(self%values)
    allocate(self%values(self%ncol,self%nobs))
    
  end subroutine obsvec_setup

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_obsvec_clone(c_key_self, c_key_other) bind(c,name='mom5cice5_obsvec_clone_f90')
    implicit none
    integer(c_int), intent(in)    :: c_key_self
    integer(c_int), intent(inout) :: c_key_other
    type(obs_vect), pointer :: self, other

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    call mom5cice5_obs_vect_registry%init()
    call mom5cice5_obs_vect_registry%add(c_key_other)
    call mom5cice5_obs_vect_registry%get(c_key_other,other)
    other%ncol=self%ncol
    other%nobs=self%nobs
    allocate(other%values(other%ncol,other%nobs))

  end subroutine c_mom5cice5_obsvec_clone

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_obsvec_delete(c_key_self) bind(c,name='mom5cice5_obsvec_delete_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(obs_vect), pointer :: self

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    deallocate(self%values)
    call mom5cice5_obs_vect_registry%remove(c_key_self)

  end subroutine c_mom5cice5_obsvec_delete

  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_assign(c_key_self, c_key_rhs) bind(c,name='mom5cice5_obsvec_assign_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_rhs
    type(obs_vect), pointer :: self, rhs

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    call mom5cice5_obs_vect_registry%get(c_key_rhs,rhs)
    if (rhs%ncol/=self%ncol .or. rhs%nobs/=self%nobs) then
       deallocate(self%values)
       self%ncol=rhs%ncol
       self%nobs=rhs%nobs
       allocate(self%values(self%ncol,self%nobs))
    endif
    self%values(:,:)=rhs%values(:,:)

  end subroutine c_mom5cice5_obsvec_assign
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_zero(c_key_self) bind(c,name='mom5cice5_obsvec_zero_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    type(obs_vect), pointer :: self

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    self%values(:,:)=0.0_kind_real

  end subroutine c_mom5cice5_obsvec_zero
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_mul_scal(c_key_self, zz) bind(c,name='mom5cice5_obsvec_mul_scal_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    real(c_double), intent(in) :: zz
    type(obs_vect), pointer :: self

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    self%values(:,:)=zz*self%values(:,:)

  end subroutine c_mom5cice5_obsvec_mul_scal
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_add(c_key_self, c_key_other) bind(c,name='mom5cice5_obsvec_add_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_other
    type(obs_vect), pointer :: self, other

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    call mom5cice5_obs_vect_registry%get(c_key_other,other)
    self%values(:,:)=self%values(:,:)+other%values(:,:)

  end subroutine c_mom5cice5_obsvec_add
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_sub(c_key_self, c_key_other) bind(c,name='mom5cice5_obsvec_sub_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_other
    type(obs_vect), pointer :: self, other

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    call mom5cice5_obs_vect_registry%get(c_key_other,other)
    self%values(:,:)=self%values(:,:)-other%values(:,:)

  end subroutine c_mom5cice5_obsvec_sub
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_mul(c_key_self, c_key_other) bind(c,name='mom5cice5_obsvec_mul_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_other
    type(obs_vect), pointer :: self, other

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    call mom5cice5_obs_vect_registry%get(c_key_other,other)
    self%values(:,:)=self%values(:,:)*other%values(:,:)

  end subroutine c_mom5cice5_obsvec_mul
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_div(c_key_self, c_key_other) bind(c,name='mom5cice5_obsvec_div_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_other
    type(obs_vect), pointer :: self, other

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    call mom5cice5_obs_vect_registry%get(c_key_other,other)
    self%values(:,:)=self%values(:,:)/other%values(:,:)

  end subroutine c_mom5cice5_obsvec_div
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_axpy(c_key_self, zz, c_key_other) bind(c,name='mom5cice5_obsvec_axpy_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    real(c_double), intent(in) :: zz
    integer(c_int), intent(in) :: c_key_other
    type(obs_vect), pointer :: self, other

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    call mom5cice5_obs_vect_registry%get(c_key_other,other)
    self%values(:,:)=self%values(:,:)+zz*other%values(:,:)

  end subroutine c_mom5cice5_obsvec_axpy
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_invert(c_key_self) bind(c,name='mom5cice5_obsvec_invert_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    type(obs_vect), pointer :: self

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    self%values(:,:)=1.0_kind_real/self%values(:,:)

  end subroutine c_mom5cice5_obsvec_invert
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_random(c_key_self) bind(c,name='mom5cice5_obsvec_random_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    type(obs_vect), pointer :: self

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    call random_vector(self%values)

  end subroutine c_mom5cice5_obsvec_random
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_dotprod(c_key_self, c_key_other, zz) bind(c,name='mom5cice5_obsvec_dotprod_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_other
    real(c_double), intent(inout) :: zz
    type(obs_vect), pointer :: self, other
    integer :: jc, jo

    zz=0.0_kind_real
    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    call mom5cice5_obs_vect_registry%get(c_key_other,other)
    do jo=1,self%nobs
       do jc=1,self%ncol
          zz = zz + self%values(jc,jo)*other%values(jc,jo)
       enddo
    enddo

  end subroutine c_mom5cice5_obsvec_dotprod
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_minmaxavg(c_key_self, zmin, zmax, zavg) bind(c,name='mom5cice5_obsvec_minmaxavg_f90')
    implicit none
    integer(c_int), intent(in)    :: c_key_self
    real(c_double), intent(inout) :: zmin, zmax, zavg
    type(obs_vect), pointer :: self
    integer :: jc, jo

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    if (self%nobs>0.and.self%ncol>0) then
       if (.not.allocated(self%values)) call abor1_ftn("obsvec_minmax: obs vector not allocated")
       zmin=self%values(1,1)
       zmax=self%values(1,1)
       zavg=0.0_kind_real
       do jo=1,self%nobs
          do jc=1,self%ncol
             zmin=MIN(self%values(jc,jo),zmin)
             zmax=MAX(self%values(jc,jo),zmax)
             zavg=zavg+self%values(jc,jo)
          enddo
       enddo
       zavg=zavg/real(self%ncol*self%nobs,kind_real)
    else
       zmin=0.0_kind_real
       zmax=0.0_kind_real
       zavg=0.0_kind_real
    endif

  end subroutine c_mom5cice5_obsvec_minmaxavg
  ! ------------------------------------------------------------------------------
  subroutine c_mom5cice5_obsvec_nobs(c_key_self, kobs) bind(c,name='mom5cice5_obsvec_nobs_f90')
    implicit none
    integer(c_int), intent(in)    :: c_key_self
    integer(c_int), intent(inout) :: kobs
    type(obs_vect), pointer :: self

    call mom5cice5_obs_vect_registry%get(c_key_self,self)
    kobs=self%nobs*self%ncol

  end subroutine c_mom5cice5_obsvec_nobs
  ! ------------------------------------------------------------------------------

end module mom5cice5_obs_vectors
