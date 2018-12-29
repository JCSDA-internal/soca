!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
!> Structure holding configuration variables for the  model

module soca_model_mod

  use kinds
  use iso_c_binding
  use soca_geom_mod
  use soca_mom6
  use soca_fields  
  use MOM,  only : step_MOM
  use mpp_domains_mod, only : mpp_update_domains
  
  implicit none

  private
  public :: soca_model
  public :: soca_model_registry
  public :: soca_create
  public :: soca_propagate
  public :: soca_delete
  
  !> Fortran derived type to hold configuration data for the  model
  type :: soca_model
     integer :: nx                !< Zonal grid dimension
     integer :: ny                !< Meridional grid dimension
     real(kind=kind_real) :: dt0  !< dimensional time (seconds)
     integer                :: advance_mom6 ! call mom6 step if true
     type(soca_mom6_config) :: mom6_config
  end type soca_model

#define LISTED_TYPE soca_model

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_model_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

  subroutine soca_create(self, geom, c_conf)

    implicit none
    
    type(soca_model), intent(inout) :: self
    type(c_ptr),         intent(in) :: c_conf
    type(soca_geom),     intent(in) :: geom

    call soca_mom6_init(self%mom6_config)

  end subroutine soca_create

  ! ------------------------------------------------------------------------------
  
  subroutine soca_propagate(self, flds)

    implicit none
    
    type(soca_model), intent(inout) :: self
    type(soca_field), intent(inout) :: flds

    integer :: isc, iec, jsc, jec
    
    ! Update halo
    call mpp_update_domains(flds%tocn, flds%geom%ocean%G%Domain%mpp_domain)
    call mpp_update_domains(flds%socn, flds%geom%ocean%G%Domain%mpp_domain)    

    ! Update MOM's T and S to soca's
    self%mom6_config%MOM_CSp%T = real(flds%tocn, kind=8)
    self%mom6_config%MOM_CSp%S = real(flds%socn, kind=8)    

    ! Advance MOM 1 baroclinic time step    
    call step_MOM(self%mom6_config%forces, &
                 &self%mom6_config%fluxes, &
                 &self%mom6_config%sfc_state, &
                 &self%mom6_config%Time, &
                 &real(self%mom6_config%dt_forcing, kind=8), &
                 &self%mom6_config%MOM_CSp)

    ! Update soca fields
    flds%tocn = real(self%mom6_config%MOM_CSp%T, kind=kind_real)
    flds%socn = real(self%mom6_config%MOM_CSp%S, kind=kind_real)
    flds%hocn = real(self%mom6_config%MOM_CSp%h, kind=kind_real)
    flds%ssh = real(self%mom6_config%MOM_CSp%ave_ssh_ibc, kind=kind_real)

  end subroutine soca_propagate

  subroutine soca_delete(self)

    implicit none
    
    type(soca_model), intent(inout) :: self

    call soca_mom6_end(self%mom6_config)

  end subroutine soca_delete
  
  
end module soca_model_mod
