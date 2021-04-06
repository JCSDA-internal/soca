! (C) Copyright 2021-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_fields_metadata_mod

use fckit_configuration_module, only: fckit_configuration, fckit_yamlconfiguration
use fckit_pathname_module, only : fckit_pathname

implicit none
private

public :: soca_field_metadata, soca_fields_metadata


type :: soca_field_metadata
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

type :: soca_fields_metadata

private
  type(soca_field_metadata), allocatable :: metadata(:)

contains
  procedure :: create => soca_fields_metadata_create
  procedure :: clone  => soca_fields_metadata_clone
  procedure :: get    => soca_fields_metadata_get
end type

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

contains

subroutine soca_fields_metadata_create(self, filename)
  class(soca_fields_metadata), intent(inout) :: self
  character(len=:), allocatable :: filename

  type(fckit_configuration)  :: conf
  type(fckit_Configuration), allocatable :: conf_list(:)

  integer :: i, j
  logical :: bool
  character(len=:), allocatable :: str

  conf = fckit_yamlconfiguration( fckit_pathname(filename))
  call conf%get_or_die("", conf_list)
  allocate(self%metadata(size(conf_list)))
  do i=1,size(self%metadata)

    call conf_list(i)%get_or_die("name", self%metadata(i)%name)

    if(.not. conf_list(i)%get("grid", str)) str = 'h'
    self%metadata(i)%grid = str

    if(.not. conf_list(i)%get("levels", str)) str = "1"
    self%metadata(i)%levels = str

    if(.not. conf_list(i)%get("masked", bool)) bool = .true.
    self%metadata(i)%masked = bool

    if(.not. conf_list(i)%get("getval name", str)) str=self%metadata(i)%name
    self%metadata(i)%getval_name = str

    if(.not. conf_list(i)%get("getval name surface", str)) str=""
    self%metadata(i)%getval_name_surface = str

    if(.not. conf_list(i)%get("io name", str)) str = ""
    self%metadata(i)%io_name = str

    if(.not. conf_list(i)%get("io file", str)) str = ""
    self%metadata(i)%io_file = str

    if(.not. conf_list(i)%get("dummy_atm", bool)) bool = .false.
    self%metadata(i)%dummy_atm = bool
  end do

  ! check for duplicates
  do i=1,size(self%metadata)
    do j=i+1,size(self%metadata)
      if ( self%metadata(i)%name == self%metadata(j)%name .or. &
           self%metadata(i)%name == self%metadata(j)%getval_name .or. &
           self%metadata(i)%name == self%metadata(j)%getval_name_surface .or. &
           self%metadata(i)%getval_name == self%metadata(j)%name .or. &
           self%metadata(i)%getval_name == self%metadata(j)%getval_name .or. &
           self%metadata(i)%getval_name == self%metadata(j)%getval_name_surface .or. &
           ( self%metadata(i)%getval_name_surface /=  "" .and. ( &
             self%metadata(i)%getval_name_surface == self%metadata(j)%name .or. &
             self%metadata(i)%getval_name_surface == self%metadata(j)%getval_name ))) then
        print *, i, self%metadata(i)%name, j, self%metadata(j)%name
        stop 100
      end if
    end do
  end do



end subroutine

! ------------------------------------------------------------------------------

subroutine soca_fields_metadata_clone(self, other)
  class(soca_fields_metadata), intent(in) :: self
  class(soca_fields_metadata), intent(out) ::  other

  integer :: i
  allocate(other%metadata(size(self%metadata)))
  do i = 1, size(self%metadata)
    other%metadata(i)%name = self%metadata(i)%name
    other%metadata(i)%getval_name = self%metadata(i)%getval_name
    other%metadata(i)%getval_name_surface = self%metadata(i)%getval_name_surface
    other%metadata(i)%io_file = self%metadata(i)%io_file
    other%metadata(i)%io_name = self%metadata(i)%io_name
    other%metadata(i)%grid = self%metadata(i)%grid
    other%metadata(i)%levels = self%metadata(i)%levels
    other%metadata(i)%masked = self%metadata(i)%masked
    other%metadata(i)%dummy_atm = self%metadata(i)%dummy_atm
  end do

end subroutine


! ------------------------------------------------------------------------------
function soca_fields_metadata_get(self, name) result(field)
  class(soca_fields_metadata), intent(in) :: self
  character(len=:), allocatable :: name
  type(soca_field_metadata) :: field

  integer :: i

  do i=1,size(self%metadata)
    if( trim(self%metadata(i)%name) == trim(name) .or. &
        trim(self%metadata(i)%getval_name) == trim(name) .or. &
        trim(self%metadata(i)%getval_name_surface) == trim(name) ) then
      field = self%metadata(i)
      return
    endif
  enddo

  print *, "DBG nope, couldnt find ", name
  stop 1

end function


end module