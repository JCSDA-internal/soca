!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_vertconv_mod

  use kinds
  use soca_fields

  implicit none

  !> Fortran derived type to hold the setup for Vertconv
  type :: soca_vertconv
     real(kind=kind_real) :: lz                    !> Vertical decorrelation
     !> scale [m]
     type(soca_field)     :: traj                  !> Trajectory
     type(soca_field)     :: bkg                   !> Background
     real(kind=kind_real), allocatable :: z(:,:,:) !> Ocean Depth [m]
     integer              :: isc, iec, jsc, jec    !> Compute domain 
  end type soca_vertconv

#define LISTED_TYPE soca_vertconv

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_vertconv_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------

  subroutine soca_conv (self, convdx, dx)

    use kinds

    implicit none

    type(soca_vertconv), intent(in) :: self
    type(soca_field),    intent(in) :: dx
    type(soca_field),   intent(inout) :: convdx

    real(kind=kind_real), allocatable :: z(:), zp(:)
    real(kind=kind_real) :: lz2, dist2, coef, lz
    integer :: nl, j, k, id, jd

    lz = self%lz
    nl = size(self%z,3)
    lz2 = lz**2

    allocate(z(nl), zp(nl))
    do id = self%isc, self%iec
       do jd = self%jsc, self%jec
          if (self%bkg%geom%ocean%mask2d(id,jd).eq.1) then
             z(:) = self%z(id,jd,:)
             zp = z
             do j = 1, nl
                convdx%tocn(id,jd,j) = 0.0d0
                convdx%socn(id,jd,j) = 0.0d0             
                do k = 1,nl
                   dist2 = (z(j)-zp(k))**2
                   coef = exp(-dist2/lz2)
                   convdx%tocn(id,jd,j) = convdx%tocn(id,jd,j) &
                        &+ dx%tocn(id,jd,k)*coef
                   convdx%socn(id,jd,j) = convdx%socn(id,jd,j) &
                        &+ dx%socn(id,jd,k)*coef
                end do
             end do
          end if
       end do
    end do
    deallocate(z, zp)

  end subroutine soca_conv

  ! ------------------------------------------------------------------------------  
  subroutine soca_conv_ad (self, convdx, dx)

    use kinds

    implicit none

    type(soca_vertconv), intent(in) :: self
    type(soca_field), intent(inout) :: dx     ! OUT
    type(soca_field),    intent(in) :: convdx ! IN

    real(kind=kind_real), allocatable :: z(:), zp(:)
    real(kind=kind_real) :: lz2, dist2, coef, lz
    integer :: nl, j, k, id, jd

    lz = self%lz
    nl = size(self%z,3)
    lz2 = lz**2

    allocate(z(nl), zp(nl))
    do id = self%isc, self%iec
       do jd = self%jsc, self%jec
          z(:) = self%z(id,jd,:)
          zp = z
          if (self%bkg%geom%ocean%mask2d(id,jd).eq.1) then
             dx%tocn(id,jd,:) = 0.0d0
             dx%socn(id,jd,:) = 0.0d0                

             do j = nl, 1, -1
                do k = nl, 1, -1
                   dist2 = (z(j)-zp(k))**2
                   coef = exp(-dist2/lz2)
                   dx%tocn(id,jd,k) = dx%tocn(id,jd,k) + coef*convdx%tocn(id,jd,j)
                   dx%socn(id,jd,k) = dx%socn(id,jd,k) + coef*convdx%socn(id,jd,j)
                end do
             end do
          end if
       end do
    end do
    deallocate(z, zp)

  end subroutine soca_conv_ad

end module soca_vertconv_mod
