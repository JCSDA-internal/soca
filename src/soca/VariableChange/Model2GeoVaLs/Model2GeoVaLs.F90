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
subroutine soca_model2geovals_changevar_f90(c_key_geom, c_key_xin, c_key_xout) &
  bind(c,name='soca_model2geovals_changevar_f90')
  integer(c_int), intent(in) :: c_key_geom, c_key_xin, c_key_xout

  type(soca_geom),  pointer :: geom
  type(soca_state), pointer :: xin, xout
  type(soca_field), pointer :: field
  integer :: i, ii, jj, kk

  type(atlas_field) :: aField
  real(kind=kind_real), pointer :: aFieldPtr(:,:)


  call soca_geom_registry%get(c_key_geom, geom)
  call soca_state_registry%get(c_key_xin, xin)
  call soca_state_registry%get(c_key_xout, xout)

  call xin%sync_from_atlas()
!
  do i=1, size(xout%fields)

    ! special cases
    select case (xout%fields(i)%name)

    ! fields that are obtained from geometry
    case ('latitude')
      xout%fields(i)%val(:,:,1) = real(geom%lat, kind=kind_real)

    case ('longitude')
      xout%fields(i)%val(:,:,1) = real(geom%lon, kind=kind_real)

    case ('sea_water_depth')
      call xin%get('hocn', field)
        xout%fields(i)%val = 0.5 * field%val
        do kk = 2, field%nz
          xout%fields(i)%val(:,:,kk) = xout%fields(i)%val(:,:,kk) + &
                                      sum(field%val(:,:,1:kk-1), dim=3)
        end do

    case ('distance_from_coast')
      aField = geom%fieldset%field("distance_from_coast")
      call aField%data(aFieldPtr)
      do jj=geom%jsc,geom%jec
        do ii=geom%isc,geom%iec
          xout%fields(i)%val(ii,jj,1) = aFieldPtr(1, geom%atlas_ij2idx(ii,jj))
        end do
      end do
      call aField%final()

    case ('sea_area_fraction')
      xout%fields(i)%val(:,:,1) = real(geom%mask2d, kind=kind_real)

    case ('mesoscale_representation_error')
      ! Representation errors: dx/R
      ! TODO, why is the halo left to 0 for RR ??
      aField = geom%fieldset%field("rossby_radius")
      call aField%data(aFieldPtr)
      do jj=geom%jsc,geom%jec
        do ii=geom%isc,geom%iec
          xout%fields(i)%val(ii,jj,1) = &
            geom%mask2d(ii, jj) * &
            sqrt(geom%cell_area(ii, jj)) / &
            aFieldPtr(1, geom%atlas_ij2idx(ii,jj))
        end do
      end do
      call aField%final()

    ! special derived state variables
    case ('skin_temperature_at_surface_where_sea')
      call xin%get('tocn', field)
      xout%fields(i)%val(:,:,1) = field%val(:,:,1) + 273.15_kind_real

    case ('sea_floor_depth_below_sea_surface')
      call xin%get('hocn', field)
      xout%fields(i)%val(:,:,1) = sum(field%val, dim=3)

    ! identity operators
    case default
      call xin%get(xout%fields(i)%metadata%name, field)
      if (xout%fields(i)%name == field%metadata%name .or. &
          xout%fields(i)%name == field%metadata%getval_name ) then
        xout%fields(i)%val(:,:,:) =  field%val(:,:,:) !< full field
      elseif (field%metadata%getval_name_surface == xout%fields(i)%name) then
        xout%fields(i)%val(:,:,1) = field%val(:,:,1) !< surface only of a 3D field
      else
        call abor1_ftn( 'error in soca_model2geovals_changevar_f90 processing ' &
                        // xout%fields(i)%name )
      endif

    end select

  end do

  call xout%sync_to_atlas()
end subroutine

!-------------------------------------------------------------------------------

end module
