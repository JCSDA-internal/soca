!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_kst_mod

  use kinds
  use soca_fields
  
  implicit none

  !> Fortran derived type to hold the setup for Kst
  type :: soca_kst
     real(kind=kind_real) :: dsdtmax !> 1.0    [psu/K]
     real(kind=kind_real) :: dsdzmin !> 3.0e-3 [psu/m] 
     real(kind=kind_real) :: dtdzmin !> 1.0e-3 [K/m]
     real(kind=kind_real), allocatable :: jacobian(:,:,:) !> dS/dT(i,j,k)
  end type soca_kst

#define LISTED_TYPE soca_kst

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

    use kinds

    implicit none

    real(kind=kind_real),                 intent(in) :: t(:), s(:), h(:)
    real(kind=kind_real),                 intent(in) :: dsdtmax, dsdzmin, dtdzmin
    real(kind=kind_real), allocatable, intent(inout) :: jac(:) ! jac=ds/dt

    real(kind=kind_real), allocatable :: dtdz(:), dsdz(:)
    !real(kind=kind_real) :: dsdtmax = 1.0    !> Need to go in conf
    !real(kind=kind_real) :: dsdzmin = 3.0e-3 !> Need to go in conf
    !real(kind=kind_real) :: dtdzmin = 1.0e-3 !> Need to go in conf

    integer :: nl

    ! Allocate
    nl = size(t,1)
    allocate(dtdz(nl),dsdz(nl))
    if (.not.allocated(jac)) allocate(jac(nl))

    call soca_diff(dtdz,t,h)
    call soca_diff(dsdz,s,h)

    ! Jacobian of soft
    jac=dsdz/dtdz

    ! Limit application of soft according to configuration
    where (abs(dtdz)<dtdzmin)
       jac=0.0
    end where
    where (abs(dsdz)<dsdzmin)
       jac=0.0
    end where    
    where (abs(jac)>dsdtmax)
       jac=0.0
    end where

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !jac(1:20) = 0.0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  end subroutine soca_soft_jacobian

  !==========================================================================
  subroutine soca_soft_tl (ds, dt, jac)
    !==========================================================================
    !
    ! Tangent of soft relative to the reference state t0, s0
    !
    ! Input:
    ! ------
    ! dt     : Potential temperature                           [deg C]
    ! jac    : ds/dt                                           [psu/K]
    !
    ! Output:
    ! -------
    ! ds  : salinity                                   [psu]
    !    
    !--------------------------------------------------------------------------

    use kinds

    implicit none

    real(kind=kind_real), intent(in)  :: dt(:), jac(:)
    real(kind=kind_real), intent(out) :: ds(:)

    integer :: nl !< Number of layers

    ! TLM
    ds = jac * dt

  end subroutine soca_soft_tl

  !==========================================================================
  subroutine soca_soft_ad (ds, dt, jac)
    !==========================================================================
    !
    ! Adjoint of tangent of soft relative to the reference state t0, s0
    !
    ! Input:
    ! ------
    ! ds  : salinity                                           [psu]
    ! jac    : ds/dt                                           [psu/K]
    !
    ! Output:
    ! -------
    ! dt     : Potential temperature                           [deg C]    
    !    
    !--------------------------------------------------------------------------

    use kinds

    implicit none

    real(kind=kind_real), intent(in)  :: ds(:), jac(:)
    real(kind=kind_real), intent(out) :: dt(:)

    ! AD of TLM
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
