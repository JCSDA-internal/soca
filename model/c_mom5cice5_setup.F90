!> Setup the model

subroutine c_mom5cice5_setup(c_confspec, c_key_geom, c_key_confdata) bind (c,name='mom5cice5_setup_f90')

  use mom5cice5_constants
  use mom5cice5_configs
  use mom5cice5_geom_mod
  use iso_c_binding
  use config_mod
  use duration_mod
  use kinds
  use fckit_log_module, only : fckit_log

  implicit none
  type(c_ptr), intent(in)    :: c_confspec         !< pointer to object of class Config
  integer(c_int), intent(in) :: c_key_geom         !< Geometry
  integer(c_int), intent(inout) :: c_key_confdata  !< Key to configuration data

  type(mom5cice5_config), pointer :: config
  type(mom5cice5_geom), pointer :: geom

  integer :: icentre, jcentre, ii, jj
  real(kind=kind_real) :: distx, disty
  type(duration) :: dtstep
  character(len=20) :: ststep
  character(len=160) :: record

  ! ------------------------------------------------------------------------------

  call mom5cice5_geom_registry%get(c_key_geom, geom)
  call mom5cice5_config_registry%init()
  call mom5cice5_config_registry%add(c_key_confdata)
  call mom5cice5_config_registry%get(c_key_confdata, config)

  config%nx  = geom%nx
  config%ny  = geom%ny
  write(record,*)'c_mom5cice5_setup: nx, ny=',config%nx,config%ny
  call fckit_log%info(record)

  ststep = config_get_string(c_confspec,len(ststep),"tstep")
  dtstep = trim(ststep)
  config%dt0 = duration_seconds(dtstep)
  write(record,*)'c_mom5cice5_setup: dt0=',config%dt0
  call fckit_log%info(record)

  ! ------------------------------------------------------------------------------
  return
end subroutine c_mom5cice5_setup
