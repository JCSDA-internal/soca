! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_vertconv_mod_c

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use kinds, only: kind_real
use soca_fields_mod, only: soca_fields
use soca_vertconv_mod, only: soca_vertconv, soca_conv_setup, &
                             soca_conv, soca_conv_ad

implicit none
private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Constructor for Vertconv
subroutine c_soca_vertconv_setup(c_self, c_conf, c_traj, c_bkg) &
    bind(c,name='soca_vertconv_setup_f90')
  type(c_ptr), intent(inout) :: c_self   !< The Vertconv structure
  type(c_ptr),    intent(in) :: c_conf       !< The configuration
  type(c_ptr),    intent(in) :: c_traj   !< trajectory
  type(c_ptr),    intent(in) :: c_bkg    !< background

  type(soca_vertconv), pointer :: self
  type(soca_fields),   pointer :: traj
  type(soca_fields),   pointer :: bkg

  allocate(self)
  c_self = c_loc(self)
  call c_f_pointer(c_traj, traj)
  call c_f_pointer(c_bkg, bkg)

  call soca_conv_setup (self, bkg, traj, fckit_configuration(c_conf))

end subroutine c_soca_vertconv_setup

! ------------------------------------------------------------------------------
!> Destructor for Vertconv
subroutine c_soca_vertconv_delete(c_self) bind(c,name='soca_vertconv_delete_f90')
  type(c_ptr), intent(inout) :: c_self  !< The background covariance structure

  type(soca_vertconv), pointer :: self

  ! Deallocate trajectory and backgroun
  ! TODO
  ! Deallocate ocean depth array
  ! TODO

  call c_f_pointer(c_self, self)

  if (associated(self%traj)) nullify(self%traj)
  if (associated(self%bkg)) nullify(self%bkg)
  deallocate(self)

end subroutine c_soca_vertconv_delete

! ------------------------------------------------------------------------------
!> Multiplication
subroutine c_soca_vertconv_mult_f90(c_self, c_a, c_m)&
    bind(c,name='soca_vertconv_mult_f90')
  type(c_ptr), intent(in) :: c_self  !< config
  type(c_ptr), intent(in) :: c_a     !< Increment in
  type(c_ptr), intent(in) :: c_m     !< Increment out

  type(soca_vertconv), pointer :: self
  type(soca_fields),   pointer :: dxa  ! in
  type(soca_fields),   pointer :: dxm  ! out
  type(soca_fields),   pointer :: traj


  call c_f_pointer(c_self, self)
  call c_f_pointer(c_a, dxa)
  call c_f_pointer(c_m, dxm)

  !< Computes dxm = Vertconv dxa
  ! dxm = dxa
  call dxm%copy( dxa)

  ! Apply forward convolution operator to T & S
  call soca_conv(self, dxm, dxa)

end subroutine c_soca_vertconv_mult_f90

! ------------------------------------------------------------------------------
!> Multiplication adjoint
subroutine c_soca_vertconv_multad_f90(c_self, c_m, c_a)&
    bind(c,name='soca_vertconv_multad_f90')
  type(c_ptr), intent(in) :: c_self  !< config
  type(c_ptr), intent(in) :: c_a     !< Increment out
  type(c_ptr), intent(in) :: c_m     !< Increment in

  type(soca_vertconv),    pointer :: self
  type(soca_fields),      pointer :: dxa
  type(soca_fields),      pointer :: dxm

  call c_f_pointer(c_self, self)
  call c_f_pointer(c_a, dxa)
  call c_f_pointer(c_m, dxm)

  ! dxa = dxm
  call dxa%copy(dxm)

  ! Apply adjoint of convolution operator
  call soca_conv_ad(self, dxm, dxa)

end subroutine c_soca_vertconv_multad_f90

end module soca_vertconv_mod_c
