!
! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bkgerrfilt_mod
  use config_mod
  use datetime_mod  
  use iso_c_binding
  use kinds
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices
  use soca_utils
  use soca_omb_stats_mod
  use fckit_mpi_module
  
  implicit none

  !> Fortran derived type to hold configuration D
  type :: soca_bkgerrfilt_config
     type(soca_field),    pointer :: bkg
     type(soca_field)             :: filt
     real(kind=kind_real)         :: efold_z           ! E-folding scale
     real(kind=kind_real)         :: scale             ! Rescaling factor
     real(kind=kind_real)         :: ocn_depth_min     ! Minimum depth
     integer                      :: isc, iec, jsc, jec
  end type soca_bkgerrfilt_config

#define LISTED_TYPE soca_bkgerrfilt_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_bkgerrfilt_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------
  !> Setup the static background error
  subroutine soca_bkgerrfilt_setup(c_conf, self, bkg)
    type(soca_bkgerrfilt_config), intent(inout) :: self
    type(soca_field),    target, intent(in) :: bkg
    type(c_ptr),                 intent(in) :: c_conf

    integer :: isc, iec, jsc, jec, i, j, k, nl
    real(kind=kind_real), allocatable :: dvdz(:), v(:), h(:)
    real(kind=kind_real) :: dt, ds, t0, s0, p, lon, lat
    real(kind=kind_real) :: detas, efold
    type(datetime) :: vdate
    character(len=800) :: fname = 'soca_bkgerrfilt.nc'
    logical :: read_from_file = .false.

    ! Get number of ocean levels
    nl = size(bkg%hocn,3)

    ! Allocate memory for bkgerrfiltor and set to zero
    call create_copy(self%filt, bkg)
    call zeros(self%filt)
    
    ! Read parameters from config
    self%ocn_depth_min = config_get_real(c_conf,"ocean_depth_min")
    self%scale         = config_get_real(c_conf,"rescale_bkgerr")
    self%efold_z       = config_get_real(c_conf,"efold_z")
    
    ! Associate background
    self%bkg => bkg

    ! Setup rescaling and masks
    call geom_get_domain_indices(bkg%geom%ocean, "compute", isc, iec, jsc, jec)
    self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec
    do i = isc, iec
       do j = jsc, jec
          if (sum(bkg%hocn(i,j,:)).gt.self%ocn_depth_min) then             
             self%filt%ssh(i,j) = self%scale
             do k = 1, nl
                if (bkg%hocn(i,j,k).gt.1e-3_kind_real) then
                   ! Only apply if layer is thick enough
                   efold = self%scale*exp(-self%bkg%layer_depth(i,j,k)/self%efold_z)
                else
                   ! Set to zero if layer is too thin
                   efold = 0.0_kind_real
                end if
                self%filt%tocn(i,j,k) = efold
                self%filt%socn(i,j,k) = efold
             end do
          else
             ! Set to zero if ocean is too shallow
             self%filt%ssh(i,j)    = 0.0_kind_real
             self%filt%tocn(i,j,:) = 0.0_kind_real
             self%filt%socn(i,j,:) = 0.0_kind_real
          end if
          
          ! Do nothing for sea-ice
          self%filt%seaice%cicen(i,j,:) =  1.0_kind_real
          self%filt%seaice%hicen(i,j,:) =  1.0_kind_real
       end do
    end do
    
    ! Save filtered background error
    call soca_fld2file(self%filt, fname)
    
  end subroutine soca_bkgerrfilt_setup

  ! ------------------------------------------------------------------------------
  !> Apply background error: dxm = D dxa
  subroutine soca_bkgerrfilt_mult(self, dxa, dxm)
    type(soca_bkgerrfilt_config),    intent(in) :: self    
    type(soca_field),            intent(in) :: dxa
    type(soca_field),         intent(inout) :: dxm

    integer :: i, j

    do i = self%isc, self%iec
       do j = self%jsc, self%jec
          if (self%bkg%geom%ocean%mask2d(i,j).eq.1) then
             dxm%ssh(i,j) = self%filt%ssh(i,j) * dxa%ssh(i,j)
             dxm%tocn(i,j,:) = self%filt%tocn(i,j,:) * dxa%tocn(i,j,:)
             dxm%socn(i,j,:) = self%filt%socn(i,j,:)  * dxa%socn(i,j,:)

             dxm%seaice%cicen(i,j,:) =  self%filt%seaice%cicen(i,j,:) * dxa%seaice%cicen(i,j,:)
             dxm%seaice%hicen(i,j,:) =  self%filt%seaice%hicen(i,j,:) * dxa%seaice%hicen(i,j,:)
          else
             dxm%ssh(i,j) = 0.0_kind_real
             dxm%tocn(i,j,:) = 0.0_kind_real
             dxm%socn(i,j,:) = 0.0_kind_real

             dxm%seaice%cicen(i,j,:) = 0.0_kind_real
             dxm%seaice%hicen(i,j,:) = 0.0_kind_real 
          end if
       end do
    end do
    ! Surface fields
    call dxm%ocnsfc%copy(dxa%ocnsfc)

  end subroutine soca_bkgerrfilt_mult

  ! ------------------------------------------------------------------------------
  
end module soca_bkgerrfilt_mod


