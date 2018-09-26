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
  use mpi,             only: mpi_comm_world
  use mpp_mod,         only: mpp_init
  use fms_mod,                 only: fms_init
  use fms_io_mod,              only: fms_io_init
  use mpp_io_mod,              only: mpp_open, mpp_close, MPP_DELETE
  use mpp_domains_mod, only: mpp_domains_init
  use fckit_mpi_module, only: fckit_mpi_comm
  
  implicit none

  type(c_ptr), intent(in) :: c_conf
  type(fckit_mpi_comm) :: f_comm

  f_comm = fckit_mpi_comm()
  
  call mpp_init(localcomm=f_comm%communicator())
  call fms_init()

end subroutine soca_setup

! ------------------------------------------------------------------------------

subroutine soca_finalize() bind(c,name='soca_finalize_f')

  use fms_io_mod,      only: fms_io_exit
  use mpp_mod,         only: mpp_exit, mpp_sync
  use mpp_io_mod,              only: mpp_open, mpp_close, MPP_DELETE
  use fms_mod,                 only: fms_end
  use fms_io_mod,              only: fms_io_exit
  use mpi
  implicit none

  !call mpi_finalize(ierr)
  !call mpp_exit()
  call fms_io_exit()
  !call mpp_sync()

end subroutine soca_finalize

! ------------------------------------------------------------------------------
