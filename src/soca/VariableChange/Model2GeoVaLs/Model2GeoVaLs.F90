module soca_model2geovals_mod

use iso_c_binding

use kinds, only: kind_real

use soca_fields_mod
use soca_fieldspec_mod
use soca_geom_mod
use soca_geom_mod_c, only: soca_geom_registry
use soca_state_mod
use soca_state_reg
use soca_increment_mod
use soca_increment_reg

use fckit_configuration_module, only: fckit_configuration

implicit none
private

contains

!-------------------------------------------------------------------------------
subroutine soca_model2geovals_linear_changevar_f90(c_key_geom, c_key_dxin, c_key_dxout) &
  bind(c,name='soca_model2geovals_linear_changevar_f90')
  integer(c_int), intent(in) :: c_key_geom, c_key_dxin, c_key_dxout

  type(soca_geom),  pointer :: geom
  type(soca_increment), pointer :: dxin, dxout
  type(soca_field), pointer :: field

  type(soca_fieldspec) :: fieldspec
  integer :: i

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_increment_registry%get(c_key_dxin, dxin)
  call soca_increment_registry%get(c_key_dxout, dxout)
!
  do i=1, size(dxout%fields)
    ! identity operators
    fieldspec =  geom%fieldspecs%get_cf(dxout%fields(i)%name)
    if (fieldspec%getval_name == dxout%fields(i)%name) then
      ! full field
      call dxin%get(fieldspec%name, field)
      dxout%fields(i)%val(:,:,:) =  field%val(:,:,:)

    elseif (fieldspec%getval_name_surface == dxout%fields(i)%name) then
      ! surface only of a 3D field
      call dxin%get(fieldspec%name,field)
      dxout%fields(i)%val(:,:,1) = field%val(:,:,1)
    endif
  end do
    ! end select
end subroutine


!-------------------------------------------------------------------------------
subroutine soca_model2geovals_linear_changevarAD_f90(c_key_geom, c_key_dxin, c_key_dxout) &
  bind(c,name='soca_model2geovals_linear_changevarAD_f90')
  integer(c_int), intent(in) :: c_key_geom, c_key_dxin, c_key_dxout

  type(soca_geom),  pointer :: geom
  type(soca_increment), pointer :: dxin, dxout
  type(soca_field), pointer :: field

  type(soca_fieldspec) :: fieldspec
  integer :: i
  logical :: b

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_increment_registry%get(c_key_dxin, dxin)
  call soca_increment_registry%get(c_key_dxout, dxout)

  do i=1, size(dxin%fields)
    fieldspec = geom%fieldspecs%get_cf(dxin%fields(i)%name)
    call dxout%get(fieldspec%name, field)

    if(fieldspec%getval_name == dxin%fields(i)%name) then
      field%val = field%val + dxin%fields(i)%val
    elseif(fieldspec%getval_name_surface == dxin%fields(i)%name) then
      field%val(:,:,1) = field%val(:,:,1) + dxin%fields(i)%val(:,:,1)
    else
      stop 51
    end if
  end do
end subroutine
!-------------------------------------------------------------------------------


subroutine soca_model2geovals_changevar_f90(c_key_geom, c_key_xin, c_key_xout) &
  bind(c,name='soca_model2geovals_changevar_f90')
  integer(c_int), intent(in) :: c_key_geom, c_key_xin, c_key_xout

  type(soca_geom),  pointer :: geom
  type(soca_state), pointer :: xin, xout
  type(soca_field), pointer :: field

  type(soca_fieldspec) :: fieldspec
  integer :: i


  call soca_geom_registry%get(c_key_geom, geom)
  call soca_state_registry%get(c_key_xin, xin)
  call soca_state_registry%get(c_key_xout, xout)
!
  do i=1, size(xout%fields)
    select case (xout%fields(i)%name)

    case ('distance_from_coast')
      xout%fields(i)%val(:,:,1) = real(geom%distance_from_coast, kind=kind_real)

    ! case ('surface_temperature_where_sea')
    !   call xin%get('tocn', field)
    !   xout%fields(i)%val = field%val

    case ('sea_floor_depth_below_sea_surface')
      call xin%get('hocn', field)
      xout%fields(i)%val(:,:,1) = sum(field%val, dim=3)

    case ('sea_area_fraction')
      xout%fields(i)%val(:,:,1) = real(geom%mask2d, kind=kind_real)

    case ('mesoscale_representation_error')
      ! Representation errors: dx/R
      ! TODO, why is the halo left to 0 for RR ??
      print *, "DBG RR ", minval(geom%rossby_radius), maxval(geom%rossby_radius)
      xout%fields(i)%val(geom%isc:geom%iec, geom%jsc:geom%jec, 1) = &
          geom%mask2d(geom%isc:geom%iec, geom%jsc:geom%jec) * &
          sqrt(geom%cell_area(geom%isc:geom%iec, geom%jsc:geom%jec) / &
               geom%rossby_radius(geom%isc:geom%iec, geom%jsc:geom%jec))


    case default
      ! identity operators
      fieldspec =  geom%fieldspecs%get_cf(xout%fields(i)%name)
      if (fieldspec%getval_name == xout%fields(i)%name) then
        ! full field
        call xin%get(fieldspec%name, field)
        xout%fields(i)%val(:,:,:) =  field%val(:,:,:)

      elseif (fieldspec%getval_name_surface == xout%fields(i)%name) then
        ! surface only of a 3D field
        call xin%get(fieldspec%name,field)
        xout%fields(i)%val(:,:,1) = field%val(:,:,1)
      endif
    end select

  end do
end subroutine

end module