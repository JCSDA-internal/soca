
!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

! ------------------------------------------------------------------------------
!> Constructor for D (standard deviation of background error)
subroutine c_soca_lengthscale_setup(c_key_self, c_conf) bind(c,name='soca_lengthscale_setup_f90')
  use iso_c_binding
  use soca_lengthscale_mod

  integer(c_int), intent(inout) :: c_key_self   !< The D structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration

  type(soca_lengthscale_config), pointer :: self

  call soca_lengthscale_registry%init()
  call soca_lengthscale_registry%add(c_key_self)
  call soca_lengthscale_registry%get(c_key_self, self)

  call soca_lengthscale_setup(c_conf, self)

end subroutine c_soca_lengthscale_setup

! ------------------------------------------------------------------------------
!> Destructor for D
subroutine c_soca_lengthscale_delete(c_key_self) bind(c,name='soca_lengthscale_delete_f90')
  use iso_c_binding
  use soca_lengthscale_mod

  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure

  call soca_lengthscale_registry%remove(c_key_self)

end subroutine c_soca_lengthscale_delete


