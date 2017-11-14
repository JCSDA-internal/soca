! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module for generating random vectors
module random_vectors_mod

use, intrinsic :: iso_c_binding
use kinds

implicit none
private
public :: random_vector

! ------------------------------------------------------------------------------
!> Fortran generic for generating random 1d, 2d and 3d arrays
interface random_vector
module procedure random_vector_1, random_vector_2, random_vector_3
end interface
! ------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
interface
subroutine random_c(kk, pp) bind(C,name='random_f')
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), intent(in) :: kk
  real(kind=c_double), intent(out) :: pp
! The line below:
! real(kind=c_double), intent(out) :: pp(:)
! would not work because then fortran passes an array descriptor
! (somebody could pass a non-contiguous array using array syntax)
! instead of the address of the first element of the array. YT
end subroutine random_c
end interface
!-------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Generate a random 1d array of reals
subroutine random_vector_1(xx)
implicit none
real(kind=kind_real), intent(inout) :: xx(:)
real(kind=c_double) :: zz(size(xx))
integer(c_int) :: nn

nn = size(xx)
call random_c(nn, zz(1))
xx(:)=zz(:)

end subroutine random_vector_1

! ------------------------------------------------------------------------------

!> Generate a random 2d array of reals
subroutine random_vector_2(xx)
implicit none
real(kind=kind_real), intent(inout) :: xx(:,:)
real(kind=c_double) :: zz(size(xx))
integer(c_int) :: nn

nn = size(xx)
call random_c(nn, zz(1))
xx = reshape(zz, shape(xx))

end subroutine random_vector_2

! ------------------------------------------------------------------------------

!> Generate a random 3d array of reals
subroutine random_vector_3(xx)
implicit none
real(kind=kind_real), intent(inout) :: xx(:,:,:)
real(kind=c_double) :: zz(size(xx))
integer(c_int) :: nn

nn = size(xx)
call random_c(nn, zz(1))
xx = reshape(zz, shape(xx))

end subroutine random_vector_3

! ------------------------------------------------------------------------------

end module random_vectors_mod
