! (C) Copyright 2021-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Metadata for soca_fields
module soca_fields_metadata_mod

use fckit_configuration_module, only: fckit_configuration, fckit_yamlconfiguration
use fckit_pathname_module, only : fckit_pathname

implicit none
private


! ------------------------------------------------------------------------------
!> Holds all of the user configurable meta data associated with a single field
!!
!! Instances of these types are to be held by soca_fields_metadata
type, public :: soca_field_metadata
  character(len=:),  allocatable :: name     !< internal name used only by soca code
  character(len=1)               :: grid     !< "h", "u" or "v"
  logical                        :: masked   !< should use land mask when interpolating
  character(len=:),  allocatable :: levels   !< "1", or "full_ocn"
  character(len=:),  allocatable :: getval_name !< variable name used by UFO
  character(len=:),  allocatable :: getval_name_surface  ! name used by UFO for the surface (if this is a 3D field)
  character(len=:),  allocatable :: io_file  !< the restart file domain (ocn, sfc, ice)
  character(len=:),  allocatable :: io_name  !< the name use in the restart IO
  character(len=:),  allocatable :: property  !< physical property of the field, "none" or "positive_definite"
  logical                        :: vert_interp   !< true if the field can be vertically interpolated
end type


! ------------------------------------------------------------------------------
!> A collection of soca_field_metadata types representing ALL possible fields
!! (state, increment, other derived) in soca. These are read in from a configuration file.
type, public :: soca_fields_metadata

  type(soca_field_metadata), private, allocatable :: metadata(:)

contains

  !> \copybrief soca_fields_metadata_create \see soca_fields_metadata_create
  procedure :: create => soca_fields_metadata_create

  !> \copybrief soca_fields_metadata_clone \see soca_fields_metadata_clone
  procedure :: clone  => soca_fields_metadata_clone

  !> \copybrief soca_fields_metadata_get \see soca_fields_metadata_get
  procedure :: get    => soca_fields_metadata_get
end type


! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Create the main soca_fields_metadata instance by reading in parameters from a
!! yaml file.
!!
!! See the members of soca_field_metadata for a list of valid options
!!
!! \throws abor1_ftn aborts if there are duplicate fields
!! \relates soca_fields_metadata_mod::soca_fields_metadata
subroutine soca_fields_metadata_create(self, filename)
  class(soca_fields_metadata), intent(inout) :: self
  character(len=:), allocatable, intent(in) :: filename !< filename of the yaml configuration

  type(fckit_configuration)  :: conf
  type(fckit_Configuration), allocatable :: conf_list(:)

  integer :: i, j
  logical :: bool
  character(len=:), allocatable :: str

  ! parse all the metadata from a yaml configuration file
  conf = fckit_yamlconfiguration( fckit_pathname(filename))
  call conf%get_or_die("", conf_list)
  allocate(self%metadata(size(conf_list)))
  do i=1,size(self%metadata)

    call conf_list(i)%get_or_die("name", self%metadata(i)%name)

    if(.not. conf_list(i)%get("grid", str)) str = 'h'
    self%metadata(i)%grid = str

    if(.not. conf_list(i)%get("masked", bool)) bool = .true.
    self%metadata(i)%masked = bool

    if(.not. conf_list(i)%get("levels", str)) str = "1"
    self%metadata(i)%levels = str

    if(.not. conf_list(i)%get("getval name", str)) str=self%metadata(i)%name
    self%metadata(i)%getval_name = str

    if(.not. conf_list(i)%get("getval name surface", str)) str=""
    self%metadata(i)%getval_name_surface = str

    if(.not. conf_list(i)%get("io name", str)) str = ""
    self%metadata(i)%io_name = str

    if(.not. conf_list(i)%get("io file", str)) str = ""
    self%metadata(i)%io_file = str

    if(.not. conf_list(i)%get("property", str)) str = "none"
    self%metadata(i)%property = str

    if(.not. conf_list(i)%get("vert interp", bool)) then
       if (self%metadata(i)%levels == "1" ) then
          bool = .false.
       else
          bool = .true.
       end if
    end if
    self%metadata(i)%vert_interp = bool

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
        str=repeat(" ",1024)
        write(str, *) "Duplicate field metadata: ", i, self%metadata(i)%name, &
                                                    j, self%metadata(j)%name
        call abor1_ftn(str)
      end if
    end do
  end do

end subroutine


! ------------------------------------------------------------------------------
!> Make a copy from \rhs to \p self
!!
!! \relates soca_fields_metadata_mod::soca_fields_metadata
subroutine soca_fields_metadata_clone(self, rhs)
  class(soca_fields_metadata), intent(inout) :: self
  class(soca_fields_metadata), intent(in) :: rhs !< metadata to clone \b from

  self%metadata = rhs%metadata

end subroutine


! ------------------------------------------------------------------------------
!> Get the metadata for the field with the given name
!!
!! The \p name can match any of \c name, \c getval_name, or \c getval_name_surface
function soca_fields_metadata_get(self, name) result(field)
  class(soca_fields_metadata), intent(in) :: self
  character(len=:), allocatable, intent(in) :: name !< the name to search for
  type(soca_field_metadata) :: field

  integer :: i

  ! find the field by any of its internal or getval names
  do i=1,size(self%metadata)
    if( trim(self%metadata(i)%name) == trim(name) .or. &
        trim(self%metadata(i)%getval_name) == trim(name) .or. &
        trim(self%metadata(i)%getval_name_surface) == trim(name) ) then
      field = self%metadata(i)
      return
    endif
  enddo

  call abor1_ftn("Unable to find field metadata for: " // name)

end function

end module
