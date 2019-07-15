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
  use mpp_mod,         only: mpp_init
  use fms_mod,                 only: fms_init
  use fckit_mpi_module, only: fckit_mpi_comm

  implicit none

  type(c_ptr), intent(in) :: c_conf
  type(fckit_mpi_comm) :: f_comm

  f_comm = fckit_mpi_comm()
  
  !call mpp_init(localcomm=f_comm%communicator())
  !call mpp_domains_init()  
  !call fms_init()

end subroutine soca_setup

! ------------------------------------------------------------------------------

subroutine soca_finalize() bind(c,name='soca_finalize_f')

  implicit none

  return

end subroutine soca_finalize

! ------------------------------------------------------------------------------
