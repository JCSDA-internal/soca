!> Fortran module for the simulation of sea-ice fraction observations
!! @todo Write detailed documentation.
module mom5cice5_fraction_mod

  use iso_c_binding
  use config_mod
  use duration_mod
  use mom5cice5_obs_data
  use mom5cice5_obs_vectors
  use mom5cice5_obsoper_mod
  use mom5cice5_vars_mod
  use mom5cice5_locs_mod
  use mom5cice5_goms_mod
  use kinds

  implicit none
  private

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_fraction_setup(c_key_self, c_conf) bind(c,name='mom5cice5_fraction_setup_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr), intent(in)       :: c_conf

    type(mom5cice5_obsoper), pointer :: self
    character(len=5) :: svars(5) = (/"cicen","hicen","vicen","hsnon","vsnon"/)
    integer :: ncol
    
    print *,'============ In fraction setup ===== ',svars
    call mom5cice5_obsoper_registry%init()
    call mom5cice5_obsoper_registry%add(c_key_self)
    call mom5cice5_obsoper_registry%get(c_key_self, self)

    ncol = 1
    call mom5cice5_oper_setup(self, c_conf, svars(:), ncol)

  end subroutine c_mom5cice5_fraction_setup

  ! ------------------------------------------------------------------------------

  subroutine c_mom5cice5_fraction_delete(c_key_self) bind(c,name='mom5cice5_fraction_delete_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self

    type(mom5cice5_obsoper), pointer :: self

    call mom5cice5_obsoper_registry%get(c_key_self, self)
    call mom5cice5_oper_delete(self)
    call mom5cice5_obsoper_registry%remove(c_key_self)

  end subroutine c_mom5cice5_fraction_delete

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_fraction_equiv(c_key_gom, c_key_hofx, c_bias) &
       & bind(c,name='mom5cice5_fraction_equiv_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_gom
    integer(c_int), intent(in) :: c_key_hofx
    real(c_double), intent(in) :: c_bias
    type(mom5cice5_goms), pointer  :: gom
    type(obs_vect), pointer :: hofx
    integer :: io, jo, ncat

    print *,'In fraction equiv fortran .... registering keys ...'
    call mom5cice5_goms_registry%get(c_key_gom, gom) 
    call mom5cice5_obs_vect_registry%get(c_key_hofx,hofx)

    print *,'In fraction equiv fortran .... DONE registering keys ...'
    !read(*,*)
    print *,'allocated numfld_per...:',allocated(gom%numfld_per_fldname)
    !read(*,*)    
    ncat=gom%numfld_per_fldname(1)
    print *,'ncat=',ncat
    do jo=1,gom%nobs
       io=gom%indx(jo)
       hofx%values(1,io)=sum(gom%values(1:ncat,jo)) + c_bias
    enddo

  end subroutine mom5cice5_fraction_equiv

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_fraction_equiv_tl(c_key_gom, c_key_hofx, c_bias) &
       & bind(c,name='mom5cice5_fraction_equiv_tl_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_gom
    integer(c_int), intent(in) :: c_key_hofx
    real(c_double), intent(in) :: c_bias
    type(mom5cice5_goms), pointer  :: gom
    type(obs_vect), pointer :: hofx
    integer :: io, jo, ncat

    call mom5cice5_goms_registry%get(c_key_gom, gom)
    call mom5cice5_obs_vect_registry%get(c_key_hofx,hofx)

    print *,' in tl .....'
    read(*,*)
    
    ncat=gom%numfld_per_fldname(1)
    do jo=1,gom%nobs
       io=gom%indx(jo)
       hofx%values(1,io)=sum(gom%values(1:ncat,jo)) + c_bias
    enddo

  end subroutine mom5cice5_fraction_equiv_tl

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_fraction_equiv_ad(c_key_gom, c_key_hofx, c_bias) &
       & bind(c,name='mom5cice5_fraction_equiv_ad_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_gom
    integer(c_int), intent(in) :: c_key_hofx
    real(c_double), intent(inout) :: c_bias
    type(mom5cice5_goms), pointer  :: gom
    type(obs_vect), pointer :: hofx
    integer :: io, jo, nco, ncat

    call mom5cice5_goms_registry%get(c_key_gom, gom)
    call mom5cice5_obs_vect_registry%get(c_key_hofx,hofx)

    print *,' in adjoint .....'
    read(*,*)
    ncat=gom%numfld_per_fldname(1)
    do jo=1,gom%nobs
       io=gom%indx(jo)
       do nco=1,ncat          
          gom%values(nco,jo)=hofx%values(1,io)
          !gom%values(nco,jo)+
       end do
       !c_bias = c_bias + hofx%values(1,io)

       !gom%values(1,jo)=hofx%values(1,io)
       !c_bias = c_bias + hofx%values(1,io)
    enddo

  end subroutine mom5cice5_fraction_equiv_ad

  ! ------------------------------------------------------------------------------

end module mom5cice5_fraction_mod
