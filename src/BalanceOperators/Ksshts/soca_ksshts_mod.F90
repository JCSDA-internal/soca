!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_ksshts_mod

  use kinds
  implicit none

  !> Fortran derived type to hold configuration Ksshts
  type :: soca_ksshts_config
     real(kind=kind_real) :: dsdtmax !> 1.0 [psu/K]
     real(kind=kind_real) :: dsdzmin !> 3.0e-3 [psu/m] 
     real(kind=kind_real) :: dtdzmin !> 1.0e-3 [K/m]
  end type soca_ksshts_config

#define LISTED_TYPE soca_ksshts_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_ksshts_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

  subroutine soca_ksshts_setup(c_conf, config)
    use iso_c_binding
    use config_mod
    use kinds

    implicit none

    type(c_ptr),              intent(in) :: c_conf   !< The configuration
    type(soca_ksshts_config), intent(inout) :: config   !< Config parameters for Ksshts

  end subroutine soca_ksshts_setup

  !==========================================================================
  subroutine soca_steric_jacobian (jac, t, s, p, h, lon, lat)
    !==========================================================================
    !
    ! Jacobian of stericnl relative to the reference state t0, s0
    !
    ! Input:
    ! ------
    ! s      : Ref. Absolute Salinity                          [g/kg]
    ! t      : Ref. Conservative Temperature                   [deg C]
    ! p      : Pressure                                        [dbar]
    ! h      : Layer thickness                                 [m]
    ! lon    : Longitude                                       [DEG E]
    ! lat    : Latitude                                        [DEG N]      
    !
    ! Output:
    ! -------
    ! jac    : Jacobian [detas/dt1, ...,detas/dtN;             [m/deg C]  
    !                    detas/ds1, ...,detas/dsN]             [m/(g/kg)]
    !    
    !--------------------------------------------------------------------------

    use gsw_mod_toolbox, only : gsw_rho, gsw_rho_first_derivatives,gsw_sa_from_sp, gsw_pt_from_ct
    use gsw_mod_kinds
    use kinds

    implicit none

    real(kind=kind_real), intent(in)  :: t, s, p, h, lon, lat
    real(kind=kind_real), intent(out) :: jac(2)
    real(kind=kind_real) :: rho0, sa, ct, lon_rot
    real(kind=kind_real) :: drhods, drhodt, drhodp, eps=1.0e-8

    ! Insitu density
    rho0 = soca_rho(s, t, p, lon, lat)
    drhodt = (soca_rho(s, t+eps, p, lon, lat)-rho0)/eps
    drhods = (soca_rho(s+eps, t, p, lon, lat)-rho0)/eps    

    jac(1)=-h*drhodt/rho0
    jac(2)=-h*drhods/rho0
    
  end subroutine soca_steric_jacobian

  !==========================================================================
  subroutine soca_steric_tl (detas, dt, ds, t0, s0, p, h, lon, lat)
    !==========================================================================
    !
    ! Tangent of stericnl relative to the reference state t0, s0
    !
    ! Input:
    ! ------
    ! ds     : Absolute Salinity                               [g/kg]
    ! dt     : Conservative Temperature                        [deg C]
    ! s0     : Ref. Absolute Salinity                          [g/kg]
    ! t0     : Ref. Conservative Temperature                   [deg C]
    ! p      : Pressure                                        [dbar]
    ! h      : Layer thickness                                 [m]
    ! lon    : Longitude                                       [DEG E]
    ! lat    : Latitude                                        [DEG N]    
    !
    ! Output:
    ! -------
    ! detas  : steric height                                   [m]
    !    
    !--------------------------------------------------------------------------

    use kinds

    implicit none

    real(kind=kind_real), intent(in) :: dt, ds, t0, s0, p, h, lon, lat
    real(kind=kind_real), intent(out) :: detas
    real(kind=kind_real) :: jac(2) 

    call soca_steric_jacobian (jac, t0, s0, p, h, lon, lat)
    detas = jac(1)*dt + jac(2)*ds

  end subroutine soca_steric_tl

  !==========================================================================
  subroutine soca_steric_ad (detas, dt_ad, ds_ad, t0, s0, p, h, lon, lat)
    !==========================================================================
    !
    ! Adjoint of sterictl relative to the trajectory t0, s0
    !
    ! Input:
    ! ------
    ! detas  : Steric Height                                   [m]    
    ! s0     : Traj. for Absolute Salinity                     [g/kg]
    ! t0     : Traj. for Conservative Temperature              [deg C]
    ! p      : Pressure                                        [dbar]
    ! h      : Layer thickness                                 [m]
    ! lon    : Longitude                                       [DEG E]
    ! lat    : Latitude                                        [DEG N]    
    !
    ! Output:
    ! -------
    !
    ! ds_ad  : Adjoint var for Absolute Salinity               [g/kg]
    ! dt_ad  : Adjoint var for Conservative Temperature        [deg C]
    !
    !--------------------------------------------------------------------------

    use kinds

    implicit none

    real(kind=kind_real),  intent(in) :: t0, s0, p, h, detas, lon, lat
    real(kind=kind_real), intent(out) :: dt_ad, ds_ad
    real(kind=kind_real) :: jac(2) 

    call soca_steric_jacobian (jac, t0, s0, p, h, lon, lat)    
    dt_ad = jac(1)*detas
    ds_ad = jac(2)*detas

  end subroutine soca_steric_ad

  elemental function soca_rho(sp, pt, p, lon, lat)
    use kinds
    use gsw_mod_toolbox, only : gsw_rho, gsw_rho_first_derivatives,gsw_sa_from_sp, gsw_ct_from_pt    
    real(kind=kind_real), intent(in)  :: pt, sp, p, lon, lat
    real(kind=kind_real) :: sa, ct, lon_rot, soca_rho

    !Rotate longitude if necessary
    lon_rot = lon
    if (lon<-180.0) lon_rot=lon+360.0
    if (lon>180.0) lon_rot=lon-360.0
    
    ! Convert practical salinity to absolute salinity    
    sa = gsw_sa_from_sp (sp, p, lon_rot, lat)

    ! Convert potential temperature to concervative temperature
    ct = gsw_ct_from_pt (sa, pt)

    ! Insitu density
    soca_rho = gsw_rho(sa,ct,p)

    return
  end function soca_rho
  
end module soca_ksshts_mod
