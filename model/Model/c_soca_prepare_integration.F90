
subroutine c_soca_prepare_integration(c_key_conf, c_key_state) &
     & bind(c,name='soca_prepare_integration_f90')

  use iso_c_binding
  use soca_fields
  use soca_configs
  use soca_constants
  use soca_mom6sis2, only : soca_models_init, soca_models_end, Coupled
use mpi,             only: mpi_comm_world
use mpp_mod,         only: mpp_init
use fms_mod,                 only: fms_init
use fms_io_mod, only : fms_io_init
  use mpp_io_mod,              only: mpp_open, mpp_close
  implicit none
  integer(c_int), intent(in) :: c_key_conf  !< Configuration structure
  integer(c_int), intent(in) :: c_key_state !< Model fields

  type(soca_config), pointer :: conf
  type(soca_field), pointer  :: flds
integer :: unit
  ! ------------------------------------------------------------------------------

  call soca_field_registry%get(c_key_state,flds)
  call soca_config_registry%get(c_key_conf, conf)

!  call mpp_init(localcomm=mpi_comm_world)
!  call fms_init
!  call fms_io_init
  
!!$  inquire(file="EGRESS", exist=flds%AOGCM%initialized)
!!$  if ( .not. flds%AOGCM%initialized ) then
!!$     print *,'---------------in prepare_integration--------------------'     
!!$     call soca_models_init(flds%AOGCM)
!!$     call mpp_open( unit, 'EGRESS' )
!!$     call mpp_close(unit)
!!$  end if
!!$
!!$  print *,'ice stuff:',flds%AOGCM%Ice%part_size

  !call soca_models_init(flds%AOGCM)
  ! -- Do stuff below

end subroutine c_soca_prepare_integration
