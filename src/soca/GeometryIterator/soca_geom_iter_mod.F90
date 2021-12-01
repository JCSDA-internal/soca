! (C) Copyright 2019-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Geometry iterator
module soca_geom_iter_mod

use kinds, only : kind_real
use soca_geom_mod, only: soca_geom

implicit none
private


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Geometry iterator
!!
!! When initialized, the iterator points to the first valid local grid cell.
!! Calls to soca_geom_iter::next() moves the iterator forward, and calls to
!! soca_geom_iter::current() retrieves the lat/lon of the current grid cell.
!! The iterator is mainly used by soca_increment_mod::soca_increment::getpoint()
!! and soca_increment_mod::soca_increment::setpoint()
type, public :: soca_geom_iter
  type(soca_geom), pointer :: geom => null() !< Geometry

  integer :: iindex = 1  !< i index of current grid point
  integer :: jindex = 1  !< j index of current grid point
  integer :: kindex = 1  !< k index of current grid point 

contains

  !> \copybrief soca_geom_iter_setup \see soca_geom_iter_setup
  procedure :: setup => soca_geom_iter_setup

  !> \copybrief soca_geom_iter_clone \see soca_geom_iter_clone
  procedure :: clone => soca_geom_iter_clone

  !> \copybrief soca_geom_iter_equals \see soca_geom_iter_equals
  procedure :: equals => soca_geom_iter_equals

  !> \copybrief soca_geom_iter_current \see soca_geom_iter_current
  procedure :: current => soca_geom_iter_current

  !> \copybrief soca_geom_iter_next \see soca_geom_iter_next
  procedure :: next => soca_geom_iter_next

end type soca_geom_iter


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Setup for the geometry iterator
!!
!! \relates soca_geom_iter_mod::soca_geom_iter
subroutine soca_geom_iter_setup(self, geom, iindex, jindex, kindex)
  class(soca_geom_iter),    intent(inout) :: self
  type(soca_geom), pointer, intent(   in) :: geom !< Pointer to geometry
  integer,         intent(   in) :: iindex, jindex, kindex  !< starting index

  ! Associate geometry
  self%geom => geom

  ! Define iindex/jindex/kindex for local tile
  self%iindex = iindex
  self%jindex = jindex
  self%kindex = kindex

end subroutine soca_geom_iter_setup


! ------------------------------------------------------------------------------
!> Clone for the geometry iterator from \p other to \p self
!!
!! \relates soca_geom_iter_mod::soca_geom_iter
subroutine soca_geom_iter_clone(self, other)
  class(soca_geom_iter), intent(inout) :: self
  type(soca_geom_iter),  intent(   in) :: other !< Other geometry iterator to clone from

  ! Associate geometry
  self%geom => other%geom

  ! Copy iindex/jindex/kindex
  self%iindex = other%iindex
  self%jindex = other%jindex
  self%kindex = other%kindex

end subroutine soca_geom_iter_clone


! ------------------------------------------------------------------------------
!> Check for the geometry iterator equality (pointing to same i/j location)
!!
!! \relates soca_geom_iter_mod::soca_geom_iter
subroutine soca_geom_iter_equals(self, other, equals)
  class(soca_geom_iter), intent( in) :: self
  type(soca_geom_iter),  intent( in) :: other  !< Other geometry iterator
  integer,               intent(out) :: equals !< Equality flag

  ! Initialization
  equals = 0

  ! Check equality
  if (associated(self%geom, other%geom)) then
    if (self%geom%iterator_dimension .eq. 2) then
      if ((self%iindex==other%iindex) .and. (self%jindex==other%jindex)) equals = 1
    elseif (self%geom%iterator_dimension .eq. 3) then
      if ((self%iindex==other%iindex) .and. (self%jindex==other%jindex) .and. &
          (self%kindex==other%kindex) ) equals = 1
    else
      call abor1_ftn('soca_geom_iter_equals: unknown geom%iterator_dimension')
    endif
  endif

end subroutine soca_geom_iter_equals


! ------------------------------------------------------------------------------
!> Get geometry iterator current lat/lon
!!
!! \throws abor1_ftn aborts if iterator is out of bounds
!! \relates soca_geom_iter_mod::soca_geom_iter
subroutine soca_geom_iter_current(self, lon, lat, depth)

  ! Passed variables
  class(soca_geom_iter), intent( in) :: self !< Geometry iterator
  real(kind_real),    intent(out) :: lat  !< Latitude
  real(kind_real),    intent(out) :: lon  !< Longitude
  real(kind_real),    intent(out) :: depth!< Depth

  real(kind_real) :: depth1d(1,1,self%geom%nzo), h1d(1,1,self%geom%nzo)

  ! Check iindex/jindex
  if (self%iindex == -1 .AND. self%jindex == -1) then
    ! special case of {-1,-1} means end of the grid
    lat = self%geom%lat(self%geom%iec,self%geom%jec)
    lon = self%geom%lon(self%geom%iec,self%geom%jec)
  elseif (self%iindex < self%geom%isc .OR. self%iindex > self%geom%iec .OR. &
          self%jindex < self%geom%jsc .OR. self%jindex > self%geom%jec) then
    ! outside of the grid
    call abor1_ftn('soca_geom_iter_current: lat/lon iterator out of bounds')
  else
    ! inside of the grid
    lat = self%geom%lat(self%iindex,self%jindex)
    lon = self%geom%lon(self%iindex,self%jindex)
  endif

  ! check kindex
  if (self%geom%iterator_dimension .eq. 2) then
    depth = -99999
  elseif (self%geom%iterator_dimension .eq. 3) then
    h1d(1,1,:) = self%geom%h(self%iindex,self%jindex,:)
    call self%geom%thickness2depth(h1d, depth1d)
    if (self%kindex == -1) then
      ! special case of {-1} means end of the grid
      depth = depth1d(1,1,self%geom%nzo)
    elseif (self%kindex == 0) then
      ! special case of the surface fields
      depth = 0;
    elseif (self%kindex < 0 .OR. self%kindex > self%geom%nzo) then
      ! out of range
      call abor1_ftn('soca_geom_iter_current: depth iterator out of bounds')
    else
      ! inside of the 3D grid
      depth = depth1d(1,1,self%kindex)
    endif
  else
    call abor1_ftn('soca_geom_iter_current: unknown geom%iterator_dimension')
  endif

end subroutine soca_geom_iter_current


! ------------------------------------------------------------------------------
!> Update geometry iterator to next point
!!
!! \todo skip over masked points
!! \relates soca_geom_iter_mod::soca_geom_iter
subroutine soca_geom_iter_next(self)
  class(soca_geom_iter), intent(inout) :: self
  integer :: iindex, jindex, kindex

  iindex = self%iindex
  jindex = self%jindex
  kindex = self%kindex

  ! increment by 1
  if (self%geom%iterator_dimension .eq. 2) then
    if (iindex.lt.self%geom%iec) then
      iindex = iindex + 1
    elseif (iindex.eq.self%geom%iec) then
      iindex = self%geom%isc
      jindex = jindex + 1
    end if

    if (jindex > self%geom%jec) then
      iindex=-1
      jindex=-1
    end if

  elseif (self%geom%iterator_dimension .eq. 3) then
    if (iindex.lt.self%geom%iec) then
      iindex = iindex + 1
    elseif (iindex.eq.self%geom%iec) then
      iindex = self%geom%isc
      if (jindex.lt.self%geom%jec) then
        jindex = jindex + 1
      elseif (jindex.eq.self%geom%jec) then
        jindex = self%geom%jsc
        kindex = kindex + 1
      end if !j loop
    end if !iloop

    if (kindex > self%geom%nzo) then
      iindex=-1
      jindex=-1
      kindex=-1
    end if !kloop
  else
    call abor1_ftn('soca_geom_iter_next: unknown geom%iterator_dimension')
  endif

  self%iindex = iindex
  self%jindex = jindex
  self%kindex = kindex

end subroutine soca_geom_iter_next
! ------------------------------------------------------------------------------

end module soca_geom_iter_mod
