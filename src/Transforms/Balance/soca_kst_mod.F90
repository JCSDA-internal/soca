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

  ! ------------------------------------------------------------------------------
contains

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

  end subroutine soca_soft_jacobian

  ! ------------------------------------------------------------------------------
  
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
