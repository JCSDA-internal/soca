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
end module soca_balanceop
