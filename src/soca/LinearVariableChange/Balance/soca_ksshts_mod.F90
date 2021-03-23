! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_ksshts_mod

use kinds, only: kind_real
use soca_utils, only: soca_rho
use gsw_mod_toolbox, only : gsw_rho, gsw_rho_first_derivatives,gsw_sa_from_sp, gsw_pt_from_ct

implicit none

private
public :: soca_ksshts, soca_steric_jacobian

!> Fortran derived type to hold configuration Ksshts
type :: soca_ksshts
   real(kind=kind_real), allocatable :: kssht(:,:,:) !> deta(i,j)/dT(i,j,k)
   real(kind=kind_real), allocatable :: ksshs(:,:,:) !> deta(i,j)/dS(i,j,k)
end type soca_ksshts

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

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

  real(kind=kind_real), intent(in)  :: t, s, p, h, lon, lat
  real(kind=kind_real), intent(out) :: jac(2)
  real(kind=kind_real) :: rho0
  real(kind=kind_real) :: drhods, drhodt, eps=1.0e-8

  ! Insitu density
  rho0 = soca_rho(s, t, p, lon, lat)
  drhodt = (soca_rho(s, t+eps, p, lon, lat)-rho0)/eps
  drhods = (soca_rho(s+eps, t, p, lon, lat)-rho0)/eps

  jac(1)=-h*drhodt/rho0
  jac(2)=-h*drhods/rho0

end subroutine soca_steric_jacobian

end module soca_ksshts_mod
