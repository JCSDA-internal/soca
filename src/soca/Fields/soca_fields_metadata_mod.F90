! (C) Copyright 2021-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Metadata for soca_fields
module soca_fields_metadata_mod

use fckit_configuration_module, only: fckit_configuration, fckit_yamlconfiguration
use fckit_pathname_module, only : fckit_pathname
use kinds, only: kind_real

implicit none
private


! ------------------------------------------------------------------------------
!> Holds all of the user configurable meta data associated with a single field
!!
!! Instances of these types are to be held by soca_fields_metadata
type, public :: soca_field_metadata
  character(len=:),  allocatable :: name     !< name used by soca and JEDI
  character(len=:),  allocatable :: name_surface  !< name used by UFO for the surface (if this is a 3D field)
  character(len=1)               :: grid     !< "h", "u" or "v"
  logical                        :: masked   !< should use land mask when interpolating
  character(len=:),  allocatable :: levels   !< "1", or "full_ocn"
  character(len=:),  allocatable :: io_file  !< the restart file domain (ocn, sfc, ice). Or if "CONSTANT" use the value in "constant_value"
  character(len=:),  allocatable :: io_name  !< the name use in the restart IO
  character(len=:),  allocatable :: io_sup_name  !< the IO name of the super set variable
  character(len=:),  allocatable :: property  !< physical property of the field, "none" or "positive_definite"
  integer                        :: categories  !< number of seaice categories
  integer                        :: category    !< category index of the seaice field
  real(kind=kind_real)           :: fillvalue
  logical                        :: vert_interp   !< true if the field can be vertically interpolated
  real(kind=kind_real)           :: constant_value !< An optional value to use globally for the field
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
!> helper function to replace a substring in a string
function templateStr(str_in, find, replace) result(str_out)
  character(len=*), intent(in) :: str_in
  character(len=*), intent(in) :: find
  character(len=*), intent(in) :: replace
  character(len=:), allocatable :: str_out
  integer :: pos

  pos = index(str_in, find)
  if (pos > 0) then
    str_out = str_in(:pos-1) // trim(adjustl(replace)) // str_in(pos+len(find):)
  end if
end function

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

  integer :: i, j, k, l
  logical :: bool
  real(kind=kind_real) :: r
  character(len=:), allocatable :: str
  character(len=10) :: str10
  real(kind=kind_real) :: val
  type(soca_field_metadata), allocatable :: metadata_tmp(:)

  ! parse all the metadata from a yaml configuration file,
  ! and temporarily store it in metadata_tmp
  conf = fckit_yamlconfiguration( fckit_pathname(filename))
  call conf%get_or_die("", conf_list)
  allocate(metadata_tmp(size(conf_list)))
  do i=1,size(metadata_tmp)

    call conf_list(i)%get_or_die("name", metadata_tmp(i)%name)

    if(.not. conf_list(i)%get("name surface", str)) str=""
    metadata_tmp(i)%name_surface = str

    if(.not. conf_list(i)%get("grid", str)) str = 'h'
    metadata_tmp(i)%grid = str

    if(.not. conf_list(i)%get("masked", bool)) bool = .true.
    metadata_tmp(i)%masked = bool

    if(.not. conf_list(i)%get("levels", str)) str = "1"
    metadata_tmp(i)%levels = str

    if(.not. conf_list(i)%get("io name", str)) str = ""
    metadata_tmp(i)%io_name = str

    if(.not. conf_list(i)%get("io sup name", str)) str = ""
    metadata_tmp(i)%io_sup_name = str

    if(.not. conf_list(i)%get("io file", str)) str = ""
    metadata_tmp(i)%io_file = str

    if(.not. conf_list(i)%get("property", str)) str = "none"
    metadata_tmp(i)%property = str

    ! if(.not. conf_list(i)%get("category", val)) val = -1
    ! metadata_tmp(i)%category = val

    if(.not. conf_list(i)%get("categories", val)) val = -1
    metadata_tmp(i)%categories = val

    if(.not. conf_list(i)%get("fill value", val)) val = 0.0
    metadata_tmp(i)%fillvalue = val

    if(.not. conf_list(i)%get("vert interp", bool)) then
       if (metadata_tmp(i)%levels == "1" ) then
          bool = .false.
       else
          bool = .true.
       end if
    end if
    metadata_tmp(i)%vert_interp = bool

    if(conf_list(i)%get("constant value", r)) then
      if (.not. metadata_tmp(i)%io_file == "") then
        str=repeat(" ", 1024)
        write(str, *) "error in field metadata file for '", metadata_tmp(i)%name, &
          "' :  'io file' cannot be set if 'constant value' is given"
        call abor1_ftn(str)
      end if
      metadata_tmp(i)%constant_value = r
      metadata_tmp(i)%io_file = "CONSTANT"
    end if
  end do

  ! break out any templated category fields
  ! 1. determine final array size
  j = 0
  do i=1,size(metadata_tmp)
    j = j + merge(metadata_tmp(i)%categories, 1, metadata_tmp(i)%categories > 0)
  end do
  allocate(self%metadata(j))
  ! 2. create final metadata with categories expanded
  j = 1
  do i=1,size(metadata_tmp)
    if (metadata_tmp(i)%categories <= 0) then
      self%metadata(j) = metadata_tmp(i)
      j = j + 1
    else
      do k=1,metadata_tmp(i)%categories
        write(str10, '(I10)') k
        self%metadata(j) = metadata_tmp(i)
        self%metadata(j)%category = k
        self%metadata(j)%name = templateStr(metadata_tmp(i)%name, "<CATEGORY>", str10)
        self%metadata(j)%io_name = templateStr(metadata_tmp(i)%io_name, "<CATEGORY>", str10)
        self%metadata(j)%name_surface = templateStr(metadata_tmp(i)%name_surface, "<CATEGORY>", str10)
        j = j + 1
      end do
    end if
  end do

  ! check for duplicates
  do i=1,size(self%metadata)
    do j=i+1,size(self%metadata)
      if ( self%metadata(i)%name == self%metadata(j)%name .or. &
           self%metadata(i)%name == self%metadata(j)%name_surface .or. &
           ( self%metadata(i)%name_surface /=  "" .and. &
             self%metadata(i)%name_surface == self%metadata(j)%name)) then
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
  character(len=*),            intent(in) :: name !< name of field to find
  type(soca_field_metadata) :: field

  integer :: i

  ! find the field by any of its internal or getval names
  do i=1,size(self%metadata)
    if( trim(self%metadata(i)%name) == trim(name) .or. &
        trim(self%metadata(i)%name_surface) == trim(name) ) then
      field = self%metadata(i)
      return
    endif
  enddo

  call abor1_ftn("Unable to find field metadata for: " // name)

end function

end module
