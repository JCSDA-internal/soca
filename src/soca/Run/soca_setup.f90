!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
!#include <fms_platform.h>

! ------------------------------------------------------------------------------

subroutine soca_setup(c_conf) bind(c,name='soca_setup_f')
  use iso_c_binding
  use config_mod

  implicit none

  type(c_ptr), intent(in) :: c_conf

end subroutine soca_setup

! ------------------------------------------------------------------------------

subroutine soca_finalize() bind(c,name='soca_finalize_f')

  implicit none

  return

end subroutine soca_finalize

! ------------------------------------------------------------------------------
