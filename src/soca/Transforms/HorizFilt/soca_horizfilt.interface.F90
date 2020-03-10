! (C) Copyright 2017-2020 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module c_soca_horizfilt_mod
  use iso_c_binding
  use fckit_configuration_module, only: fckit_configuration
  use soca_horizfilt_mod
  use soca_geom_mod, only : soca_geom
  use soca_fields_mod
  use oops_variables_mod

  implicit none
  private

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------
  !> Setup for the filtering operator

  subroutine c_soca_horizfilt_setup(c_self, c_conf, c_geom, c_traj, c_vars) &
    & bind (c,name='soca_horizfilt_setup_f90')
    type(c_ptr),      intent(inout) :: c_self   !< The filtering structure
    type(c_ptr),         intent(in) :: c_conf       !< The configuration
    type(c_ptr), target, intent(in) :: c_geom   !< Geometry
    type(c_ptr),         intent(in) :: c_traj   !< Trajectory
    type(c_ptr),  value, intent(in) :: c_vars       !< List of variables

    type(soca_horizfilt_type), pointer :: self
    type(soca_geom),           pointer :: geom
    type(soca_fields),         pointer :: traj
    type(oops_variables)               :: vars

    allocate(self)
    c_self = c_loc(self)
    call c_f_pointer(c_geom, geom)
    call c_f_pointer(c_traj, traj)
    vars = oops_variables(c_vars)

    call soca_horizfilt_setup(self, fckit_configuration(c_conf), geom, traj, vars)

  end subroutine c_soca_horizfilt_setup

  ! ------------------------------------------------------------------------------
  !> Delete filtering operator

  subroutine c_soca_horizfilt_delete(c_self) &
    & bind (c,name='soca_horizfilt_delete_f90')
    type(c_ptr), intent(inout) :: c_self  !< The filtering structure

    type(soca_horizfilt_type), pointer :: self

    call c_f_pointer(c_self, self)

    call soca_horizfilt_delete(self)
    deallocate(self)

  end subroutine c_soca_horizfilt_delete

  ! ------------------------------------------------------------------------------
  !> Multiply

  subroutine c_soca_horizfilt_mult(c_self, c_in, c_out) &
    & bind(c,name='soca_horizfilt_mult_f90')
    type(c_ptr),    intent(in) :: c_self  !< The filtering structure
    type(c_ptr),    intent(in) :: c_in    !<    "   to Increment in
    type(c_ptr), intent(inout) :: c_out   !<    "   to Increment out

    type(soca_horizfilt_type), pointer :: self
    type(soca_fields),         pointer :: xin, xout

    call c_f_pointer(c_self, self)
    call c_f_pointer(c_in, xin)
    call c_f_pointer(c_out, xout)

    call soca_horizfilt_mult(self, xin, xout) !< xout = C.xout

  end subroutine c_soca_horizfilt_mult

  ! ------------------------------------------------------------------------------
  !> Multiply adjoint

  subroutine c_soca_horizfilt_mult_ad(c_self, c_in, c_out) &
    & bind(c,name='soca_horizfilt_multad_f90')
    type(c_ptr),    intent(in) :: c_self  !< The filtering structure
    type(c_ptr),    intent(in) :: c_in    !<    "   to Increment in
    type(c_ptr), intent(inout) :: c_out   !<    "   to Increment out

    type(soca_horizfilt_type),  pointer :: self
    type(soca_fields),          pointer :: xin, xout

    call c_f_pointer(c_self, self)
    call c_f_pointer(c_in, xin)
    call c_f_pointer(c_out, xout)

    call soca_horizfilt_multad(self, xin, xout) !< xout = C^T.xout

  end subroutine c_soca_horizfilt_mult_ad

end module c_soca_horizfilt_mod
