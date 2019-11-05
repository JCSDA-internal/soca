! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_kst_mod

use kinds, only: kind_real
use soca_fields_mod, only: soca_field
use soca_utils, only: soca_diff

implicit none

private
public :: soca_kst, soca_soft_jacobian

!> Fortran derived type to hold the setup for Kst
type :: soca_kst
   real(kind=kind_real) :: dsdtmax !> 1.0    [psu/K]
   real(kind=kind_real) :: dsdzmin !> 3.0e-3 [psu/m]
   real(kind=kind_real) :: dtdzmin !> 1.0e-3 [K/m]
   integer              :: nlayers
   real(kind=kind_real), allocatable :: jacobian(:,:,:) !> dS/dT(i,j,k)
end type soca_kst

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!==========================================================================
subroutine soca_soft_jacobian (jac, t, s, h, dsdtmax, dsdzmin, dtdzmin)
  !==========================================================================
  !
  ! Jacobian of Sb=S(T) relative to the reference state t, s. jac=dS/dT at (t,s)
  !
  ! Input:
  ! ------
  ! s      : Background practical salinity                   [g/kg]
  ! t      : Background potential Temperature                [deg C]
  ! h      : Layer thickness                                 [m]
  ! config : Configuration for soft
  !
  ! Output:
  ! -------
  ! jac    : Jacobian [ds1/dt1, ...,dsN/dtN];                [psu/deg C]
  !
  !--------------------------------------------------------------------------
  real(kind=kind_real),                 intent(in) :: t(:), s(:), h(:)
  real(kind=kind_real),                 intent(in) :: dsdtmax, dsdzmin, dtdzmin
  real(kind=kind_real), allocatable, intent(inout) :: jac(:) ! jac=ds/dt

  real(kind=kind_real), allocatable :: dtdz(:), dsdz(:)
  integer :: nl, z
  real(kind=kind_real) :: j

  ! Allocate
  nl = size(t,1)
  allocate(dtdz(nl),dsdz(nl))
  if (.not.allocated(jac)) allocate(jac(nl))

  call soca_diff(dtdz,t,h)
  call soca_diff(dsdz,s,h)

  jac = 0.0
  do z=1,nl
    jac(z) = 0.0

    ! Limit application of soft according to configuration
    if ( abs(dtdz(z)) < dtdzmin ) cycle
    if ( abs(dsdz(z)) < dsdzmin ) cycle

    ! Jacobian of soft
    j=dsdz(z)/dtdz(z)

    ! Limit application of soft according to configuration
    if ( abs(j) > dsdtmax ) cycle

    ! if we reach this point in the code, the jacobian is usable
    jac(z) = j;
  end do

end subroutine soca_soft_jacobian

end module soca_kst_mod
