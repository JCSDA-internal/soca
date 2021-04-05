! (C) Copyright 2021-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_fieldspec_mod

use fckit_configuration_module, only: fckit_configuration, fckit_yamlconfiguration
use fckit_pathname_module, only : fckit_pathname

implicit none
private

public :: soca_fieldspecs, soca_fieldspec


type :: soca_fieldspec
  character(len=:),  allocatable :: name !< internal name used only by soca code
  character(len=:),  allocatable :: getval_name !< variable name used by UFO
  character(len=:),  allocatable :: getval_name_surface  ! internal name used by UFO for the suface (if this is a 3D field)
  character(len=1)               :: grid     !< "h", "u" or "v"
  logical                        :: masked   !< should use land mask when interpolating
  logical                        :: dummy_atm !< a meaningless dummy field, for the CRTM hacks
  character(len=:),  allocatable :: io_file  !< the restart file domain (ocn, sfc, ice)
  character(len=:),  allocatable :: io_name  !< the name use in the restart IO
  character(len=:),  allocatable :: levels   !< "1", or "full_ocn"

end type

! ------------------------------------------------------------------------------

type :: soca_fieldspecs

private
  type(soca_fieldspec), allocatable :: fieldspecs(:)

contains
  procedure :: create => soca_fieldspecs_create
  procedure :: clone  => soca_fieldspecs_clone
  procedure :: get    => soca_fieldspecs_get
end type

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

contains

subroutine soca_fieldspecs_create(self, filename)
  class(soca_fieldspecs), intent(inout) :: self
  character(len=:), allocatable :: filename

  type(fckit_configuration)  :: conf
  type(fckit_Configuration), allocatable :: conf_list(:)

  integer :: i, j
  logical :: bool
  character(len=:), allocatable :: str

  conf = fckit_yamlconfiguration( fckit_pathname(filename))
  call conf%get_or_die("", conf_list)
  allocate(self%fieldspecs(size(conf_list)))
  do i=1,size(self%fieldspecs)

    call conf_list(i)%get_or_die("name", self%fieldspecs(i)%name)

    if(.not. conf_list(i)%get("grid", str)) str = 'h'
    self%fieldspecs(i)%grid = str

    if(.not. conf_list(i)%get("levels", str)) str = "1"
    self%fieldspecs(i)%levels = str

    if(.not. conf_list(i)%get("masked", bool)) bool = .true.
    self%fieldspecs(i)%masked = bool

    if(.not. conf_list(i)%get("getval name", str)) str=self%fieldspecs(i)%name
    self%fieldspecs(i)%getval_name = str

    if(.not. conf_list(i)%get("getval name surface", str)) str=""
    self%fieldspecs(i)%getval_name_surface = str

    if(.not. conf_list(i)%get("io name", str)) str = ""
    self%fieldspecs(i)%io_name = str

    if(.not. conf_list(i)%get("io file", str)) str = ""
    self%fieldspecs(i)%io_file = str

    if(.not. conf_list(i)%get("dummy_atm", bool)) bool = .false.
    self%fieldspecs(i)%dummy_atm = bool
  end do

  ! check for duplicates
  do i=1,size(self%fieldspecs)
    do j=i+1,size(self%fieldspecs)
      if ( self%fieldspecs(i)%name == self%fieldspecs(j)%name .or. &
           self%fieldspecs(i)%name == self%fieldspecs(j)%getval_name .or. &
           self%fieldspecs(i)%name == self%fieldspecs(j)%getval_name_surface .or. &
           self%fieldspecs(i)%getval_name == self%fieldspecs(j)%name .or. &
           self%fieldspecs(i)%getval_name == self%fieldspecs(j)%getval_name .or. &
           self%fieldspecs(i)%getval_name == self%fieldspecs(j)%getval_name_surface .or. &
           ( self%fieldspecs(i)%getval_name_surface /=  "" .and. ( &
             self%fieldspecs(i)%getval_name_surface == self%fieldspecs(j)%name .or. &
             self%fieldspecs(i)%getval_name_surface == self%fieldspecs(j)%getval_name ))) then
        print *, i, self%fieldspecs(i)%name, j, self%fieldspecs(j)%name
        stop 100
      end if
    end do
  end do



end subroutine

! ------------------------------------------------------------------------------

subroutine soca_fieldspecs_clone(self, other)
  class(soca_fieldspecs), intent(in) :: self
  class(soca_fieldspecs), intent(out) ::  other

  integer :: i
  allocate(other%fieldspecs(size(self%fieldspecs)))
  do i = 1, size(self%fieldspecs)
    other%fieldspecs(i)%name = self%fieldspecs(i)%name
    other%fieldspecs(i)%getval_name = self%fieldspecs(i)%getval_name
    other%fieldspecs(i)%getval_name_surface = self%fieldspecs(i)%getval_name_surface
    other%fieldspecs(i)%io_file = self%fieldspecs(i)%io_file
    other%fieldspecs(i)%io_name = self%fieldspecs(i)%io_name
    other%fieldspecs(i)%grid = self%fieldspecs(i)%grid
    other%fieldspecs(i)%levels = self%fieldspecs(i)%levels
    other%fieldspecs(i)%masked = self%fieldspecs(i)%masked
    other%fieldspecs(i)%dummy_atm = self%fieldspecs(i)%dummy_atm
  end do

end subroutine


! ------------------------------------------------------------------------------
function soca_fieldspecs_get(self, name) result(field)
  class(soca_fieldspecs), intent(in) :: self
  character(len=:), allocatable :: name
  type(soca_fieldspec) :: field

  integer :: i

  do i=1,size(self%fieldspecs)
    if( trim(self%fieldspecs(i)%name) == trim(name) .or. &
        trim(self%fieldspecs(i)%getval_name) == trim(name) .or. &
        trim(self%fieldspecs(i)%getval_name_surface) == trim(name) ) then
      field = self%fieldspecs(i)
      return
    endif
  enddo

  print *, "DBG nope, couldnt find ", name
  stop 1

end function


end module