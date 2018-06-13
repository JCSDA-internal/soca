!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_seaice_balanceop

  implicit none
  
contains
  
  !==========================================================================
  subroutine tofc_tl (dt, dcn, tb, sb, cnb)
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
    
  end subroutine tofc_tl

  !==========================================================================
  subroutine tofc_ad (dt, dcn, tb, sb, cnb)
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
    
  end subroutine tofc_ad
  
end module soca_seaice_balanceop
