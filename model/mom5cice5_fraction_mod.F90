
!> Fortran module for fraction observations for the QG model
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
    character(len=1) :: svars(1) = (/"x"/)

    call mom5cice5_obsoper_registry%init()
    call mom5cice5_obsoper_registry%add(c_key_self)
    call mom5cice5_obsoper_registry%get(c_key_self, self)

    call mom5cice5_oper_setup(self, c_conf, svars, 1)

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
    integer :: io, jo

    call mom5cice5_goms_registry%get(c_key_gom, gom) 
    call mom5cice5_obs_vect_registry%get(c_key_hofx,hofx)

    do jo=1,gom%nobs
       io=gom%indx(jo)
       hofx%values(1,io)=sum(gom%values(1:gom%geom%ncat,jo)) + c_bias
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
    integer :: io, jo

    call mom5cice5_goms_registry%get(c_key_gom, gom)
    call mom5cice5_obs_vect_registry%get(c_key_hofx,hofx)

    do jo=1,gom%nobs
       io=gom%indx(jo)
       hofx%values(1,io)=sum(gom%values(1:gom%geom%ncat,jo)) + c_bias       
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
    integer :: io, jo, nco

    call mom5cice5_goms_registry%get(c_key_gom, gom)
    call mom5cice5_obs_vect_registry%get(c_key_hofx,hofx)

    do jo=1,gom%nobs
       io=gom%indx(jo)
       do nco=1,gom%geom%ncat
          gom%values(nco,jo)=gom%values(nco,jo)+hofx%values(1,io)
       end do
       c_bias = c_bias + hofx%values(1,io)
       ! CODE ADJOINT BELOW
       !gom%values(1,jo)=hofx%values(1,io)
       !c_bias = c_bias + hofx%values(1,io)
    enddo

  end subroutine mom5cice5_fraction_equiv_ad

  ! ------------------------------------------------------------------------------

end module mom5cice5_fraction_mod
