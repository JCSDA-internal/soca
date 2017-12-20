
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
implicit none

type(c_ptr), intent(in) :: c_conf
integer :: stackmax = 4000000, unit

call mpp_init(localcomm=mpi_comm_world)
!call mpp_domains_init
call fms_io_init

!call mpp_open( unit, 'EGRESS' )
!call mpp_close(unit, MPP_DELETE)
!if (config_element_exists(c_conf,"stackmax")) stackmax = config_get_int(c_conf,"stackmax")
!call mpp_domains_set_stack_size(stackmax)

end subroutine soca_setup

! ------------------------------------------------------------------------------

subroutine soca_finalize() bind(c,name='soca_finalize_f')

use fms_io_mod,      only: fms_io_exit
  !use soca_mom6sis2
use mpp_mod,         only: mpp_exit
    use mpp_io_mod,              only: mpp_open, mpp_close, MPP_DELETE
    use fms_mod,                 only: fms_end
    use fms_io_mod,              only: fms_io_exit
implicit none
integer :: unit
    !call mpp_open( unit, 'EGRESS' )
    !call mpp_close(unit, MPP_DELETE)
!call fms_io_exit()
    print *,'============================================================='
    print *,'============================================================='
    print *,'============================================================='
    print *,'===================calling mpp exit======================='
!call mpp_exit()
!call fms_io_exit()

end subroutine soca_finalize

! ------------------------------------------------------------------------------
