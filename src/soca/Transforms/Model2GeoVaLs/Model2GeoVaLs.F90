! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interface for converting model variables to geovals (mostly identity function)
module soca_model2geovals_mod_c

use iso_c_binding
use kinds, only: kind_real

! soca modules
use soca_fields_mod, only: soca_field
use soca_geom_mod_c, only: soca_geom_registry
use soca_geom_mod, only: soca_geom
use soca_increment_mod, only: soca_increment
use soca_increment_reg, only: soca_increment_registry
use soca_state_mod, only: soca_state
use soca_state_reg, only: soca_state_registry

implicit none
private


!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> C++ interface for linear change of variables from model to geovals
!!
!! Only the identity operator is needed for the linear variables
!! \throws abor1_ftn aborts if the field name cannot be in the "getval_name*"
!! section of the variable metadata
subroutine soca_model2geovals_linear_changevar_f90(c_key_geom, c_key_dxin, c_key_dxout) &
  bind(c,name='soca_model2geovals_linear_changevar_f90')
  integer(c_int), intent(in) :: c_key_geom, c_key_dxin, c_key_dxout

  type(soca_geom),  pointer :: geom
  type(soca_increment), pointer :: dxin, dxout
  type(soca_field), pointer :: field
  integer :: i

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_increment_registry%get(c_key_dxin, dxin)
  call soca_increment_registry%get(c_key_dxout, dxout)

  ! identity operators
  do i=1, size(dxout%fields)
    call dxin%get(dxout%fields(i)%metadata%name, field)

    if (field%metadata%getval_name == dxout%fields(i)%name) then
      dxout%fields(i)%val(:,:,:) =  field%val(:,:,:)  !< full field
    elseif (field%metadata%getval_name_surface == dxout%fields(i)%name) then
      dxout%fields(i)%val(:,:,1) = field%val(:,:,1) !< surface only of a 3D field
    else
      call abor1_ftn( 'error in soca_model2geovals_linear_changevar_f90 processing ' &
                       // dxout%fields(i)%name )
    endif

  end do
end subroutine


!-------------------------------------------------------------------------------
!> C++ interface for linear change of variables from geovals to model
!!
!! Only the identity operator is need for the linear variables.
!! \throws abor1_ftn aborts if the field name cannot be in the "getval_name*"
!! section of the variable metadata
subroutine soca_model2geovals_linear_changevarAD_f90(c_key_geom, c_key_dxin, c_key_dxout) &
  bind(c,name='soca_model2geovals_linear_changevarAD_f90')
  integer(c_int), intent(in) :: c_key_geom, c_key_dxin, c_key_dxout

  type(soca_geom),  pointer :: geom
  type(soca_increment), pointer :: dxin, dxout
  type(soca_field), pointer :: field
  integer :: i

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_increment_registry%get(c_key_dxin, dxin)
  call soca_increment_registry%get(c_key_dxout, dxout)

  ! identity operators
  do i=1, size(dxin%fields)
    call dxout%get(dxin%fields(i)%metadata%name, field)

    if(field%metadata%getval_name == dxin%fields(i)%name) then
      field%val = field%val + dxin%fields(i)%val !< full field
    elseif(field%metadata%getval_name_surface == dxin%fields(i)%name) then
      field%val(:,:,1) = field%val(:,:,1) + dxin%fields(i)%val(:,:,1) !< surface only
    else
      call abor1_ftn( 'error in soca_model2geovals_linear_changevarAD_f90 processing ' &
                       // dxin%fields(i)%name )
    end if

  end do
end subroutine


!-------------------------------------------------------------------------------
!> C++ interface for the non-linear change of variables from model to geovals
!!
!! This is *mostly* an identity operator, except for a small number of derived variables
!! that are to be calculated here ("distance_from_coast", "sea_area_fraction", etc.)
!! \throws abor1_ftn aborts if field name is not handled.
subroutine soca_model2geovals_changevar_f90(c_key_geom, c_key_xin, c_key_xout) &
  bind(c,name='soca_model2geovals_changevar_f90')
  integer(c_int), intent(in) :: c_key_geom, c_key_xin, c_key_xout

  type(soca_geom),  pointer :: geom
  type(soca_state), pointer :: xin, xout
  type(soca_field), pointer :: field
  integer :: i

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_state_registry%get(c_key_xin, xin)
  call soca_state_registry%get(c_key_xout, xout)
!
  do i=1, size(xout%fields)

    ! special cases
    select case (xout%fields(i)%name)

    ! fields that are obtained from geometry
    case ('distance_from_coast')
      xout%fields(i)%val(:,:,1) = real(geom%distance_from_coast, kind=kind_real)

    case ('sea_area_fraction')
      xout%fields(i)%val(:,:,1) = real(geom%mask2d, kind=kind_real)

    case ('mesoscale_representation_error')
      ! Representation errors: dx/R
      ! TODO, why is the halo left to 0 for RR ??
      xout%fields(i)%val(geom%isc:geom%iec, geom%jsc:geom%jec, 1) = &
          geom%mask2d(geom%isc:geom%iec, geom%jsc:geom%jec) * &
          sqrt(geom%cell_area(geom%isc:geom%iec, geom%jsc:geom%jec) / &
               geom%rossby_radius(geom%isc:geom%iec, geom%jsc:geom%jec))

    ! special derived state variables
    case ('surface_temperature_where_sea')
      call xin%get('tocn', field)
      xout%fields(i)%val(:,:,1) = field%val(:,:,1) + 273.15_kind_real

    case ('sea_floor_depth_below_sea_surface')
      call xin%get('hocn', field)
      xout%fields(i)%val(:,:,1) = sum(field%val, dim=3)

    ! identity operators
    case default
      call xin%get(xout%fields(i)%metadata%name, field)
      if (field%metadata%getval_name == xout%fields(i)%name) then
        xout%fields(i)%val(:,:,:) =  field%val(:,:,:) !< full field
      elseif (field%metadata%getval_name_surface == xout%fields(i)%name) then
        xout%fields(i)%val(:,:,1) = field%val(:,:,1) !< surface only of a 3D field
      else
        call abor1_ftn( 'error in soca_model2geovals_changevar_f90 processing ' &
                        // xout%fields(i)%name )
      endif

    end select

  end do
end subroutine

!-------------------------------------------------------------------------------

end module
