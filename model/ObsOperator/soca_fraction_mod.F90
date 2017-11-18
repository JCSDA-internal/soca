!> Fortran module for the simulation of sea-ice fraction observations
!! @todo Write detailed documentation.
module soca_fraction_mod

  use iso_c_binding
  use config_mod
  use duration_mod
  use soca_obs_data
  use soca_obs_vectors
  use soca_obsoper_mod
  use soca_vars_mod
  use soca_locs_mod
  use soca_goms_mod
  use kinds

  implicit none
  private

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------

  subroutine c_soca_fraction_setup(c_key_self, c_conf) bind(c,name='soca_fraction_setup_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr), intent(in)       :: c_conf

    type(soca_obsoper), pointer :: self
    character(len=5) :: svars(5) = (/"cicen","hicen","vicen","hsnon","vsnon"/)
    integer :: ncol
    
    print *,'============ In fraction setup ===== ',svars
    call soca_obsoper_registry%init()
    call soca_obsoper_registry%add(c_key_self)
    call soca_obsoper_registry%get(c_key_self, self)

    ncol = 1
    call soca_oper_setup(self, c_conf, svars(:), ncol)

  end subroutine c_soca_fraction_setup

  ! ------------------------------------------------------------------------------

  subroutine c_soca_fraction_delete(c_key_self) bind(c,name='soca_fraction_delete_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self

    type(soca_obsoper), pointer :: self

    call soca_obsoper_registry%get(c_key_self, self)
    call soca_oper_delete(self)
    call soca_obsoper_registry%remove(c_key_self)

  end subroutine c_soca_fraction_delete

  ! ------------------------------------------------------------------------------

  subroutine soca_fraction_equiv(c_key_gom, c_key_hofx, c_bias) &
       & bind(c,name='soca_fraction_equiv_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_gom
    integer(c_int), intent(in) :: c_key_hofx
    real(c_double), intent(in) :: c_bias
    type(soca_goms), pointer  :: gom
    type(obs_vect), pointer :: hofx
    integer :: io, jo, ncat

    print *,'In fraction equiv fortran .... registering keys ...'
    !read(*,*)
    call soca_goms_registry%get(c_key_gom, gom) 
    call soca_obs_vect_registry%get(c_key_hofx,hofx)

    print *,'In fraction equiv fortran .... DONE registering keys ...'
    !read(*,*)
    print *,'allocated numfld_per...:',allocated(gom%numfld_per_fldname),gom%numfld_per_fldname(1)
    ncat=gom%numfld_per_fldname(1)
    print *,'ncat=',ncat,' nobs=',gom%nobs
    !read(*,*)        
    do jo=1,gom%nobs
       io=gom%indx(jo)
       hofx%values(1,io)=sum(gom%values(1:ncat,jo))! + c_bias
    enddo
    print *,'end nl'
    
  end subroutine soca_fraction_equiv

  ! ------------------------------------------------------------------------------

  subroutine soca_fraction_equiv_tl(c_key_gom, c_key_hofx, c_bias) &
       & bind(c,name='soca_fraction_equiv_tl_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_gom
    integer(c_int), intent(in) :: c_key_hofx
    real(c_double), intent(in) :: c_bias
    type(soca_goms), pointer  :: gom
    type(obs_vect), pointer :: hofx
    integer :: io, jo, ncat

    call soca_goms_registry%get(c_key_gom, gom)
    call soca_obs_vect_registry%get(c_key_hofx,hofx)

    ncat=gom%numfld_per_fldname(1)
    do jo=1,gom%nobs
       io=gom%indx(jo)
       hofx%values(1,io)=sum(gom%values(1:ncat,jo))! + c_bias
    enddo

  end subroutine soca_fraction_equiv_tl

  ! ------------------------------------------------------------------------------

  subroutine soca_fraction_equiv_ad(c_key_gom, c_key_hofx, c_bias) &
       & bind(c,name='soca_fraction_equiv_ad_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_gom
    integer(c_int), intent(in) :: c_key_hofx
    real(c_double), intent(inout) :: c_bias
    type(soca_goms), pointer  :: gom
    type(obs_vect), pointer :: hofx
    integer :: io, jo, nco, ncat

    call soca_goms_registry%get(c_key_gom, gom)
    call soca_obs_vect_registry%get(c_key_hofx,hofx)

    print *,' in adjoint .....shape of gom%values:',shape(gom%values)
    print *,'shape of hofx%values:',shape(hofx%values)

    ncat=gom%numfld_per_fldname(1)
    do jo=1,gom%nobs
       io=gom%indx(jo)
       do nco=1,ncat
          gom%values(nco,jo)=hofx%values(1,io)
       end do
    enddo
    print *,'end ad'
    !read(*,*)       

  end subroutine soca_fraction_equiv_ad

  ! ------------------------------------------------------------------------------

end module soca_fraction_mod
