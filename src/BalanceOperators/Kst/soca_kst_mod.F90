!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_kst_mod

  use kinds
  implicit none

  !> Fortran derived type to hold configuration Kst
  type :: soca_kst_config
     real(kind=kind_real) :: dsdtmax !> 1.0 [psu/K]
     real(kind=kind_real) :: dsdzmin !> 3.0e-3 [psu/m] 
     real(kind=kind_real) :: dtdzmin !> 1.0e-3 [K/m]
  end type soca_kst_config

#define LISTED_TYPE soca_kst_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_kst_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

  subroutine soca_kst_setup(c_conf, config)
    use iso_c_binding
    use config_mod
    use kinds

    implicit none

    type(c_ptr),              intent(in) :: c_conf   !< The configuration
    type(soca_kst_config), intent(inout) :: config   !< Config parameters for Kst

    config%dsdtmax      = config_get_real(c_conf,"dsdtmax")
    config%dsdzmin      = config_get_real(c_conf,"dsdzmin")
    config%dtdzmin      = config_get_real(c_conf,"dtdzmin")

  end subroutine soca_kst_setup

  !==========================================================================
  subroutine soca_soft_jacobian (jac, t, s, h)
    !==========================================================================
    !
    ! Jacobian of Sb=S(T) relative to the reference state t, s. jac=dS/dT at (t,s)
    !
    ! Input:
    ! ------
    ! s      : Background practical salinity                   [g/kg]
    ! t      : Background potential Temperature                [deg C]
    ! h      : Layer thickness                                 [m]
    !
    ! Output:
    ! -------
    ! jac    : Jacobian [ds1/dt1, ...,dsN/dtN];                [psu/deg C]
    !    
    !--------------------------------------------------------------------------

    use kinds

    implicit none

    real(kind=kind_real),               intent(in) :: t(:), s(:), h(:)
    real(kind=kind_real), allocatable, intent(out) :: jac(:)

    real(kind=kind_real), allocatable :: dtdz(:), dsdz(:)
    real(kind=kind_real) :: dsdtmax = 1.0    !> Need to go in conf
    real(kind=kind_real) :: dsdzmin = 3.0e-3 !> Need to go in conf
    real(kind=kind_real) :: dtdzmin = 1.0e-3 !> Need to go in conf

    integer :: nl

    ! Allocate
    nl = size(t,1)
    allocate(dtdz(nl),dsdz(nl),jac(nl))

    call soca_diff(dtdz,t,h)
    call soca_diff(dsdz,s,h)
    jac=0.0
    where (abs(dtdz)>dtdzmin)
       jac=dsdz/dtdz
    end where
    where (abs(jac)>dsdtmax)
       jac=0.0
    end where

  end subroutine soca_soft_jacobian

  !==========================================================================
  subroutine soca_soft_tl (ds, dt, t0, s0, h)
    !==========================================================================
    !
    ! Tangent of soft relative to the reference state t0, s0
    !
    ! Input:
    ! ------
    ! dt     : Potential temperature                           [deg C]
    ! s0     : Salinity trajectory                             [psu]
    ! t0     : Potential Temperature trajectory                [deg C]
    ! h      : Layer thickness                                 [m]
    !
    ! Output:
    ! -------
    ! ds  : salinity                                   [psu]
    !    
    !--------------------------------------------------------------------------

    use kinds

    implicit none

    real(kind=kind_real), intent(in)  :: dt(:), t0(:), s0(:), h(:)
    real(kind=kind_real), intent(out) :: ds(:)

    real(kind=kind_real), allocatable :: jac(:)    !< Mid-layer depth

    integer :: nl !< Number of layers

    ! Allocate
    nl = size(dt,1)
    allocate(jac(nl))

    ! Compute Jacobian
    call soca_soft_jacobian (jac, t0, s0, h)

    ! TLM
    ds = jac * dt

  end subroutine soca_soft_tl

  !==========================================================================
  subroutine soca_soft_ad (ds, dt, t0, s0, h)
    !==========================================================================
    !
    ! Adjoint of tangent of soft relative to the reference state t0, s0
    !
    ! Input:
    ! ------
    ! ds  : salinity                                           [psu]
    ! s0     : Salinity trajectory                             [psu]
    ! t0     : Potential Temperature trajectory                [deg C]
    ! h      : Layer thickness                                 [m]
    !
    ! Output:
    ! -------
    ! dt     : Potential temperature                           [deg C]    
    !    
    !--------------------------------------------------------------------------

    use kinds

    implicit none

    real(kind=kind_real), intent(in)  :: ds(:), t0(:), s0(:), h(:)
    real(kind=kind_real), intent(out) :: dt(:)

    real(kind=kind_real), allocatable :: jac(:)    !< Mid-layer depth

    integer :: nl !< Number of layers

    ! Allocate
    nl = size(dt,1)
    allocate(jac(nl))

    ! Compute Jacobian
    call soca_soft_jacobian (jac, t0, s0, h)

    ! TLM
    dt = jac * ds

  end subroutine soca_soft_ad

  subroutine soca_diff(dvdz,v,h)
    use kinds

    implicit none

    real(kind=kind_real), intent(in)  :: v(:), h(:)
    real(kind=kind_real), intent(out) :: dvdz(:)

    integer :: k, ik

    k = size(v,1)
    do ik = 2, k-1
       dvdz(ik) = (v(ik+1)-v(ik-1))/(h(ik)+0.5*h(ik+1)+h(ik-1))
    end do
    dvdz(1) = dvdz(2)
    dvdz(k) = dvdz(k-1)    

  end subroutine soca_diff

end module soca_kst_mod
