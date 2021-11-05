! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> variable transform: SSH balance
module soca_ksshts_mod

use kinds, only: kind_real
use soca_utils, only: soca_rho
use gsw_mod_toolbox, only : gsw_rho, gsw_rho_first_derivatives,gsw_sa_from_sp, gsw_pt_from_ct

implicit none
private
public :: soca_steric_jacobian

!> Hold the jacobians for the ssh to s/t balance transform
!!
!! should be populated by calls to soca_ksshts_mod::soca_ksshts::soca_steric_jacobian
!! \see soca_balance_mod::soca_balance_setup
type, public :: soca_ksshts
   real(kind=kind_real), allocatable :: kssht(:,:,:) !< deta(i,j)/dT(i,j,k)
   real(kind=kind_real), allocatable :: ksshs(:,:,:) !< deta(i,j)/dS(i,j,k)
end type soca_ksshts


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Jacobian of stericnl relative to the reference state t0, s0
!!
!! \relates soca_ksshts_mod::soca_ksshts
!! \param s: ref. Absolute Salinity [g/kg]
!! \param t: Ref. Conservative Temperature [deg C]
!! \param p: Pressure [dbar]
!! \param h: Layer thickness [m]
!! \param lon: Longitude [DEG E]
!! \param lat: Latitude[DEG N]
!! \param jac: Jacobian
!!   - [0] = detas/dt1, ...,detas/dtN; [m/deg C]
!!   - [1] = detas/ds1, ...,detas/dsN; [m/(g/kg)]
subroutine soca_steric_jacobian (jac, t, s, p, h, lon, lat)
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
