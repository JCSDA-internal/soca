!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_ktc_mod

  use kinds
  implicit none

  !> Fortran derived type to hold configuration Ktc
  type :: soca_ktc_config

  end type soca_ktc_config

#define LISTED_TYPE soca_ktc_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_ktc_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

  subroutine soca_ktc_setup(c_conf, config)
    use iso_c_binding
    use config_mod
    use kinds

    implicit none

    type(c_ptr),              intent(in) :: c_conf   !< The configuration
    type(soca_ktc_config), intent(inout) :: config   !< Config parameters for Ktc


  end subroutine soca_ktc_setup

  !==========================================================================
  subroutine soca_tofc_tl (dt, dcn, tb, sb, cnb)
    !==========================================================================
    !
    !
    ! Input:
    ! ------
    ! dcn    : Seaice fraction for each category               [1]        
    ! ds     : Practical Salinity                              [psu]    
    ! tb     : Ref. Potential Temperature                      [deg C]
    ! sb     : Ref. Practical Salinity                         [psu]
    ! cnb    : Ref. sea-ice concentration                      [psu]        
    !
    ! Output:
    ! -------
    ! dt     : Potential Temperature                           [deg C]    
    !    
    !--------------------------------------------------------------------------

    use kinds

    implicit none

    real(kind=kind_real), intent(in) :: dcn(:), tb, sb, cnb(:)
    real(kind=kind_real), intent(out) :: dt

    real(kind=kind_real) :: dc, a, b, cb, alpha=2.0, mu=0.054, tm
    integer :: ncat

    ncat=size(dcn)
    dc=sum(dcn)
    cb=sum(cnb)
    tm=-mu*sb

    if (cb.lt.1.0) then    
       a=-alpha!(tm-tb)/(1.0-cb)
       !a=(tm-tb)/(1.0-cb)
    else
       a=0
    end if

    dt = a*dc
    
  end subroutine soca_tofc_tl

  !==========================================================================
  subroutine soca_tofc_ad (dt, dcn, tb, sb, cnb)
    !==========================================================================
    !
    ! Adjoint
    !
    ! Input:
    ! ------
    !
    ! Output:
    ! -------
    !
    !
    !--------------------------------------------------------------------------
    use kinds
    implicit none    

    real(kind=kind_real), intent(in) :: dt, tb, sb, cnb(:)
    real(kind=kind_real), allocatable, intent(out) :: dcn(:)

    real(kind=kind_real) :: dc, a, b, cb, alpha=2.0, mu=0.054, tm
    integer :: ncat

    ncat=size(dcn)
    allocate(dcn(ncat))
    dc=sum(dcn)
    cb=sum(cnb)
    tm=-mu*sb

    if (cb.lt.1.0) then    
       a=-alpha!(tm-tb)/(1.0-cb)
       !a=(tm-tb)/(1.0-cb)       
    else
       a=0.0
    end if
    
    dcn = a*dt
    
  end subroutine soca_tofc_ad

  
end module soca_ktc_mod
