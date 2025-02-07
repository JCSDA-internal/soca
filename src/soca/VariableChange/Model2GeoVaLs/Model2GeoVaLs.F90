! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interface for converting model variables to geovals (mostly identity function)
module soca_model2geovals_mod_c

use iso_c_binding
use kinds, only: kind_real
use atlas_module, only: atlas_field

! soca modules
use soca_fields_mod, only: soca_field
use soca_geom_mod_c, only: soca_geom_registry
use soca_geom_mod, only: soca_geom
use soca_state_mod, only: soca_state
use soca_state_reg, only: soca_state_registry

implicit none
private


!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> C++ interface for the non-linear change of variables from model to geovals
!!
!! This is *mostly* an identity operator, except for a small number of derived variables
!! that are to be calculated here ("distance_from_coast", "sea_area_fraction", etc.)
!! \throws abor1_ftn aborts if field name is not handled.
!!
!! TODO: can this be moved to the pure C++ side?? Probably yes
subroutine soca_model2geovals_changevar_f90(c_key_geom, c_key_xin, c_key_xout) &
  bind(c,name='soca_model2geovals_changevar_f90')
  integer(c_int), intent(in) :: c_key_geom, c_key_xin, c_key_xout

  type(soca_geom),  pointer :: geom
  type(soca_state), pointer :: xin, xout
  type(soca_field), pointer :: soca_fieldin, soca_fieldout
  integer :: i, ii, jj, kk, idx

  type(atlas_field) :: field_out, field2
  real(kind=kind_real), pointer :: data_out(:,:), data2(:,:)

  call soca_geom_registry%get(c_key_geom, geom)
  call soca_state_registry%get(c_key_xin, xin)
  call soca_state_registry%get(c_key_xout, xout)

  do i=1, xout%afieldset%size()
    field_out = xout%afieldset%field(i)
    call field_out%data(data_out)

    ! special cases
    select case (field_out%name())

    ! fields that are obtained from geometry
    case ('latitude')
      do jj=geom%jsc,geom%jec;
        do ii=geom%isc,geom%iec
          data_out(1, geom%atlas_ij2idx(ii,jj)) = real(geom%lat(ii,jj), kind=kind_real)
        end do
      end do
      call field_out%set_dirty()

    case ('longitude')
      do jj=geom%jsc,geom%jec
        do ii=geom%isc,geom%iec
          data_out(1, geom%atlas_ij2idx(ii,jj)) = real(geom%lon(ii,jj), kind=kind_real)
        end do
      end do
      call field_out%set_dirty()

    case ('sea_water_depth')
      field2 = xin%afieldset%field("sea_water_cell_thickness");
      call field2%data(data2)
      do ii = 1, field2%shape(2)
        data_out(1, ii) = 0.5 * data2(1, ii)
        do kk =2, field2%shape(1)
          data_out(kk, ii) = data_out(kk-1,ii) +  0.5 * (data2(kk, ii) + data2(kk-1, ii))
        end do
      end do
      call field2%final()

    case ('distance_from_coast')
      field2 = geom%fieldset%field("distance_from_coast")
      call field2%data(data2)
      do ii = 1, field2%shape(2)
        data_out(1, ii) = data2(1, ii)
      end do
      call field2%final()

    case ('sea_area_fraction')
      do jj=geom%jsc,geom%jec
        do ii=geom%isc,geom%iec
          data_out(1, geom%atlas_ij2idx(ii,jj)) = real(geom%mask2d(ii,jj), kind=kind_real)
        end do
      end do
      call field_out%set_dirty()

    case ('mesoscale_representation_error')
      ! Representation errors: dx/R
      field2 = geom%fieldset%field("rossby_radius")
      call field2%data(data2)
      do jj=geom%jsc,geom%jec
        do ii=geom%isc,geom%iec
          idx = geom%atlas_ij2idx(ii,jj)
          data_out(1, idx) = geom%mask2d(ii,jj) * &
              sqrt(geom%cell_area(ii, jj)) / &
              data2(1, idx)
        end do
      end do
      call field2%final()
      call field_out%set_dirty()

    ! special derived state variables
    case ('skin_temperature_at_surface_where_sea')
      field2 = xin%afieldset%field("sea_water_potential_temperature")
      call field2%data(data2)
      do ii = 1, field2%shape(2)
        do kk = 1, field2%shape(1)
          ! TODO: should we be skipping land??
          data_out(kk, ii) = data2(kk, ii) + 273.15_kind_real
        end do
      end do
      call field2%final()

    case ('sea_floor_depth_below_sea_surface')
      field2 = xin%afieldset%field("sea_water_cell_thickness")
      call field2%data(data2)
      do ii = 1, field2%shape(2)
        data_out(1, ii) = data2(1, ii)
        do kk = 2, field2%shape(1)
          data_out(1, ii) = data_out(1, ii) + data2(kk, ii)
        end do
      end do
      call field2%final()

    ! identity operators
    case default
      ! TODO remove the dependency on the fields structure (requires
      ! chaning how metadata is stored)
      call xout%get(field_out%name(), soca_fieldout)
      call xin%get(soca_fieldout%metadata%name, soca_fieldin)
      field2 = xin%afieldset%field(soca_fieldin%name)
      call field2%data(data2)
      if (field_out%name() == soca_fieldin%metadata%name ) then
        ! full 3D field
        do ii = 1, field2%shape(2)
          do kk = 1, field2%shape(1)
            data_out(kk, ii) = data2(kk, ii)
          end do
        end do
      elseif (soca_fieldin%metadata%name_surface == field_out%name()) then
        ! surface only of a 3D field
        do ii = 1, field2%shape(2)
          data_out(1, ii) = data2(1, ii)
        end do
      else
        call abor1_ftn( 'error in soca_model2geovals_changevar_f90 processing ' &
                        // field_out%name() )
      end if
      call field2%final()
    end select

    call field_out%final()
  end do
end subroutine

!-------------------------------------------------------------------------------

end module
