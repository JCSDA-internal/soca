!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!> Perform a timestep of the model

subroutine c_soca_propagate(c_key_conf, c_key_state) bind(c,name='soca_propagate_f90')

use iso_c_binding
use soca_fields
use soca_configs

implicit none
integer(c_int), intent(in) :: c_key_conf  !< Config structure
integer(c_int), intent(in) :: c_key_state !< Model fields

type(soca_config), pointer :: conf
type(soca_field),  pointer :: flds

! ------------------------------------------------------------------------------

call soca_config_registry%get(c_key_conf, conf)
call soca_field_registry%get(c_key_state,flds)

call propagate(flds, conf)

! ------------------------------------------------------------------------------
return
end subroutine c_soca_propagate
