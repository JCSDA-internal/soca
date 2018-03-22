!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

!> Perform a timestep of the model

subroutine propagate(flds,config)

  use soca_fields
  use soca_configs
  use soca_constants
  use kinds

  implicit none
  type(soca_field),  intent(inout) :: flds
  type(soca_config), intent(in)    :: config

  print *,'Advance model (not!).'

  ! ------------------------------------------------------------------------------
  return
end subroutine propagate
