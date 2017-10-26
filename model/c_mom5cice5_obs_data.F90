
!> Fortran module handling observations for the QG model

module c_mom5cice5_obs_data

  use mom5cice5_obs_data
  use mom5cice5_obs_vectors
  use mom5cice5_obsoper_mod
  use iso_c_binding
  use config_mod
  use fckit_log_module, only : fckit_log
  use string_f_c_mod

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_obs_setup(c_key_self, c_conf) bind(c,name='mom5cice5_obsdb_setup_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(c_ptr), intent(in)    :: c_conf !< configuration

    type(obs_data), pointer :: self
    character(len=max_string) :: fin, fout
    character(len=max_string+30) :: record

    print *,'=================IN OBS SETUP ==================='
    
    if (config_element_exists(c_conf,"ObsData.ObsDataIn")) then
       fin  = config_get_string(c_conf,max_string,"ObsData.ObsDataIn.obsfile")       
    else
       fin  = ""
    endif
    write(record,*)'mom5cice5_obs_setup: file in =',trim(fin)
    call fckit_log%info(record)

    fout = config_get_string(c_conf,max_string, "ObsData.ObsDataOut.obsfile")
    write(record,*)'mom5cice5_obs_setup: file out=',trim(fout)
    call fckit_log%info(record)

    call obs_data_registry%init()
    call obs_data_registry%add(c_key_self)
    call obs_data_registry%get(c_key_self, self)
    call obs_setup(trim(fin), trim(fout), self)

  end subroutine mom5cice5_obs_setup

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_obs_delete(c_key_self) bind(c,name='mom5cice5_obsdb_delete_f90')
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    type(obs_data), pointer :: self

    call obs_data_registry%get(c_key_self, self)
    call obs_delete(self)
    call obs_data_registry%remove(c_key_self)

  end subroutine mom5cice5_obs_delete

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_obs_get(c_key_self, lreq, c_req, lcol, c_col, c_key_ovec) bind(c,name='mom5cice5_obsdb_get_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: lreq, lcol
    character(kind=c_char,len=1), intent(in) :: c_req(lreq+1), c_col(lcol+1)
    integer(c_int), intent(in) :: c_key_ovec

    type(obs_data), pointer :: self
    type(obs_vect), pointer :: ovec
    character(len=lreq) :: req
    character(len=lcol) :: col

    call obs_data_registry%get(c_key_self, self)
    call mom5cice5_obs_vect_registry%get(c_key_ovec,ovec)
    call c_f_string(c_req, req)
    call c_f_string(c_col, col)

    call obs_get(self, trim(req), trim(col), ovec)

  end subroutine mom5cice5_obs_get

  ! ------------------------------------------------------------------------------

  subroutine mom5cice5_obs_put(c_key_self, lreq, c_req, lcol, c_col, c_key_ovec) bind(c,name='mom5cice5_obsdb_put_f90')
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: lreq, lcol
    character(kind=c_char,len=1), intent(in) :: c_req(lreq+1), c_col(lcol+1)
    integer(c_int), intent(in) :: c_key_ovec

    type(obs_data), pointer :: self
    type(obs_vect), pointer :: ovec
    character(len=lreq) :: req
    character(len=lcol) :: col

    call obs_data_registry%get(c_key_self, self)
    call mom5cice5_obs_vect_registry%get(c_key_ovec,ovec)
    call c_f_string(c_req, req)
    call c_f_string(c_col, col)

    call obs_put(self, trim(req), trim(col), ovec)

  end subroutine mom5cice5_obs_put

  ! ------------------------------------------------------------------------------

end module c_mom5cice5_obs_data
