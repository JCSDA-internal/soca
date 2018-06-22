!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_balanceop

  implicit none
  
contains
  !==========================================================================
  subroutine soca_steric_jacobian (jac, t, s, p, h)
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
    !
    ! Output:
    ! -------
    ! jac    : Jacobian [detas/dt1, ...,detas/dtN;             [m/deg C]  
    !                    detas/ds1, ...,detas/dsN]             [m/(g/kg)]
    !    
    !--------------------------------------------------------------------------

    use gsw_mod_toolbox, only : gsw_rho, gsw_rho_first_derivatives
    use gsw_mod_kinds
    use kinds

    implicit none

    real(kind=kind_real), intent(in)  :: t, s, p, h
    real(kind=kind_real), intent(out) :: jac(2)
    real(kind=kind_real) :: rho0
    real(kind=kind_real) :: drhods, drhodt, drhodp

    rho0 = gsw_rho(s,t,p)
    call gsw_rho_first_derivatives(s,t,p,drhods, drhodt, drhodp)
    jac(1)=-h*drhodt/rho0
    jac(2)=-h*drhods/rho0
     !jac(3)=(rho-rho0)/rho0 !=detas/dh 
    
  end subroutine soca_steric_jacobian

  !==========================================================================
  subroutine soca_steric_tl (detas, dt, ds, t0, s0, p, h)
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
    !
    ! Output:
    ! -------
    ! detas  : steric height                                   [m]
    !    
    !--------------------------------------------------------------------------

    use kinds

    implicit none

    real(kind=kind_real), intent(in) :: dt, ds, t0, s0, p, h
    real(kind=kind_real), intent(out) :: detas
    real(kind=kind_real) :: jac(2) 

    call soca_steric_jacobian (jac, t0, s0, p, h)
    detas = jac(1)*dt + jac(2)*ds

  end subroutine soca_steric_tl

  !==========================================================================
  subroutine soca_steric_ad (detas, dt_ad, ds_ad, t0, s0, p, h)
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

    real(kind=kind_real), intent(in) :: t0, s0, p, h, detas
    real(kind=kind_real),            intent(out) :: dt_ad, ds_ad
    real(kind=kind_real) :: jac(2) 

    call soca_steric_jacobian (jac, t0, s0, p, h)    
    dt_ad = jac(1)*detas
    ds_ad = jac(2)*detas

  end subroutine soca_steric_ad

end module soca_balanceop
