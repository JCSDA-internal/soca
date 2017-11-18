!> Setup the model

subroutine c_soca_setup(c_confspec, c_key_geom, c_key_confdata) bind (c,name='soca_setup_f90')

  use soca_constants
  use soca_configs
  use soca_geom_mod
  use iso_c_binding
  use config_mod
  use duration_mod
  use kinds
  use fckit_log_module, only : fckit_log

  implicit none
  type(c_ptr), intent(in)    :: c_confspec         !< pointer to object of class Config
  integer(c_int), intent(in) :: c_key_geom         !< Geometry
  integer(c_int), intent(inout) :: c_key_confdata  !< Key to configuration data

  type(soca_config), pointer :: config
  type(soca_geom), pointer :: geom

  integer :: icentre, jcentre, ii, jj
  real(kind=kind_real) :: distx, disty
  type(duration) :: dtstep
  character(len=20) :: ststep
  character(len=160) :: record

  ! ------------------------------------------------------------------------------

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_config_registry%init()
  call soca_config_registry%add(c_key_confdata)
  call soca_config_registry%get(c_key_confdata, config)

  config%nx  = geom%nx
  config%ny  = geom%ny
  write(record,*)'c_soca_setup: nx, ny=',config%nx,config%ny
  call fckit_log%info(record)

  ststep = config_get_string(c_confspec,len(ststep),"tstep")
  dtstep = trim(ststep)
  config%dt0 = duration_seconds(dtstep)
  write(record,*)'c_soca_setup: dt0=',config%dt0
  call fckit_log%info(record)

  ! ------------------------------------------------------------------------------
  return
end subroutine c_soca_setup
