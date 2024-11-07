! (C) Copyright 2017-2024 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


!> Handle fields for the model.
!!
!! soca_fields represents a state or increment, and contains one or more
!! soca_field instances for each of the fields. The metadata associated
!! with a given field is stored in soca_fields_metadata_mod::soca_fields_metadata
module soca_fields_mod

use atlas_module, only: atlas_fieldset, atlas_field, atlas_real, atlas_metadata

! JEDI modules
use datetime_mod, only: datetime, datetime_set, datetime_to_string, datetime_to_string_io, &
                        datetime_create, datetime_diff
use duration_mod, only: duration, duration_to_string
use fckit_configuration_module, only: fckit_configuration
use logger_mod
use fckit_mpi_module, only: fckit_mpi_min, fckit_mpi_max, fckit_mpi_sum
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables
use tools_const, only: deg2rad

! MOM6 / FMS modules
use fms_io_mod, only: fms_io_init, fms_io_exit, register_restart_field, &
                      restart_file_type, restore_state, free_restart_type, save_restart, &
                      file_exist, field_exist
use fms_mod,    only: write_data, set_domain
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h, &
                          end_remapping
use mpp_domains_mod, only : mpp_update_domains

! SOCA modules
use soca_fields_metadata_mod, only : soca_field_metadata
use soca_geom_mod, only : soca_geom
use soca_utils, only: soca_mld
use soca_utils, only: soca_stencil_interp, soca_stencil_neighbors

implicit none
private


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Holds all data and metadata related to a single field variable.
!!
!! Instances of these types are to be held by soca_fields.
!! The members soca_field::mask can remain \c null, in which it is assumed that
!! no mask is used.
type, public :: soca_field

  !> The internally used name of the field.
  character(len=:),     allocatable :: name

  !> The number of vertical levels.
  integer                           :: nz

  !> The actual field data.
  real(kind=kind_real), allocatable :: val(:,:,:)

  !> Pointer to the relevant mask in soca_geom_mod::soca_geom
  !!
  !! If \c null, it is assumed that no mask is present
  real(kind=kind_real),     pointer :: mask(:,:) => null()!!

  !> Pointer to the relevant longitudes in soca_geom_mod::soca_geom
  !!
  !! \note This should never remain \c null() after initialization of the class.
  real(kind=kind_real),     pointer :: lon(:,:) => null()

  !> Pointer to the relevant latitudes in soca_geom_mod::soca_geom
  !!
  !! \note This should never remain \c null() after initialization of the class.
  real(kind=kind_real),     pointer :: lat(:,:) => null()

  !> Parameters for the field as determined by the configuration yaml.
  !!
  !! see soca_fields_metadata_mod::soca_field_metadata
  type(soca_field_metadata)         :: metadata

contains

  !>\copybrief soca_field_copy \see soca_field_copy
  procedure :: copy            => soca_field_copy

  !>\copybrief soca_field_delete \see soca_field_delete
  procedure :: delete          => soca_field_delete

  !>\copybrief soca_field_check_congruent \see soca_field_check_congruent
  procedure :: check_congruent => soca_field_check_congruent

  !>\copybrief soca_field_update_halo \see soca_field_update_halo
  procedure :: update_halo     => soca_field_update_halo

  !>\copybrief soca_field_stencil_interp \see soca_field_stencil_interp
  procedure :: stencil_interp  => soca_field_stencil_interp

  !>\copybrief soca_field_fill_masked \see soca_field_fill_masked
  procedure :: fill_masked     => soca_field_fill_masked

end type soca_field


! ------------------------------------------------------------------------------
!> A collection of soca_field types representing a collective state or increment.
!!
!! The base class for soca_increment_mod::soca_increment and soca_state_mod::soca_state
type, public :: soca_fields

  !> Pointer to the relevant soca_geom_mod::soca_geom
  !!
  !! \note This should never remain \c null() after initialization of the class.
  type(soca_geom),  pointer :: geom => null()

  !> The soca_field instances that make up the fields
  type(soca_field), allocatable :: fields(:)

  type(atlas_fieldset) :: aFieldset

contains
  !> \name constructors / destructors
  !! \{

  !> \copybrief soca_fields_create \see soca_fields_create
  procedure :: create => soca_fields_create

  !> \copybrief soca_fields_copy \see soca_fields_copy
  procedure :: copy   => soca_fields_copy

  !> \copybrief soca_fields_delete \see soca_fields_delete
  procedure :: delete => soca_fields_delete

  !> \}

  !> \name field getters/checkers
  !! \{

  !> \copybrief soca_fields_get \see soca_fields_get
  procedure :: get    => soca_fields_get

  !> \copybrief soca_fields_has \see soca_fields_has
  procedure :: has    => soca_fields_has

  !> \copybrief soca_fields_check_congruent \see soca_fields_check_congruent
  procedure :: check_congruent => soca_fields_check_congruent

  !> \copybrief soca_fields_check_subset \see soca_fields_check_subset
  procedure :: check_subset    => soca_fields_check_subset

  !> \}

  !> \name math operators
  !! \{

  !> \copybrief soca_fields_ones \see soca_fields_ones
  procedure :: ones     => soca_fields_ones

  !> \copybrief soca_fields_zeros \see soca_fields_zeros
  procedure :: zeros    => soca_fields_zeros

  !> \}

  !> \name I/O
  !! \{

  !> \copybrief soca_fields_read \see soca_fields_read
  procedure :: read      => soca_fields_read
  procedure, private :: read_seaice => soca_fields_read_seaice

  !> \copybrief soca_fields_write_file \see soca_fields_write_file
  procedure :: write_file=> soca_fields_write_file

  !> \copybrief soca_fields_write_rst \see soca_fields_write_rst
  procedure :: write_rst => soca_fields_write_rst

  !> \}

  !> \name misc
  !! \{

  !> \copybrief soca_fields_update_halos \see soca_fields_update_halos
  procedure :: update_halos => soca_fields_update_halos

  !> \copybrief soca_fields_tohpoints \see soca_fields_tohpoints
  procedure :: tohpoints  => soca_fields_tohpoints
  !> \}


  !> \copybrief soca_fields_update_fields \see soca_fields_update_fields
  procedure :: update_fields => soca_fields_update_fields
  procedure :: update_metadata => soca_fields_update_metadata

  !> \name Temporary sync between C++ atlas and our fortran array.
  !> This will go away once the transition to atlas is complete
  !! \{
  procedure :: sync_to_atlas => soca_fields_sync_to_atlas
  procedure :: sync_from_atlas => soca_fields_sync_from_atlas
  !> \}

end type soca_fields


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
contains


! ------------------------------------------------------------------------------
! soca_field subroutines
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Copy a field from \p rhs to \p self.
!!
!! If the fields are not congruent, this subroutine will throw an error.
!! \p self must be allocated first.
!! \relates soca_fields_mod::soca_field
subroutine soca_field_copy(self, rhs)
  class(soca_field), intent(inout) :: self !< The field to copy \b to
  type(soca_field),  intent(in)    :: rhs !< The field to copy \b from

  call self%check_congruent(rhs)

  ! the only variable that should be different is %val
  self%val = rhs%val

  ! NOTE: the pointers (mask, lat, lon) will be different, but should NOT
  ! be changed to point to rhs pointers. Bad things happen
end subroutine soca_field_copy


! ------------------------------------------------------------------------------
!> Update the data in the halo region of the field.
!!
!! \relates soca_fields_mod::soca_field
!! \todo have field keep a pointer to its relevant sections of soca_geom?
subroutine soca_field_update_halo(self, geom)
  class(soca_field),     intent(inout) :: self
  type(soca_geom), pointer, intent(in) :: geom !< soca_geom from soca_fields

  call mpp_update_domains(self%val, geom%Domain%mpp_domain)
end subroutine soca_field_update_halo

! ------------------------------------------------------------------------------
!> Make sure the two fields are the same in terms of name, size, shape.
!!
!! \throws abor1_ftn Halts program if fields are not congruent
!! \relates soca_fields_mod::soca_field
subroutine soca_field_check_congruent(self, rhs)
  class(soca_field), intent(in) :: self
  type(soca_field),  intent(in) :: rhs !< other field to check for congruency
  integer :: i

  if ( self%nz /= rhs%nz ) call abor1_ftn("soca_field:  self%nz /= rhs%nz")
  if ( self%name /= rhs%name ) call abor1_ftn("soca_field:  self%name /= rhs%name")
  if ( size(shape(self%val)) /= size(shape(rhs%val)) ) &
    call abor1_ftn("soca_field: shape of self%val /= rhs%val")
  do i =1, size(shape(self%val))
    if (size(self%val, dim=i) /= size(rhs%val, dim=i)) &
      call abor1_ftn("soca_field: shape of self%val /= rhs%val")
  end do
end subroutine soca_field_check_congruent

! ------------------------------------------------------------------------------
!> Perform spatial interpolation between adjacent grid point in the same stencil
!!
!! Interpolation used is inverse distance weidghted, taking into
!! consideration the mask and using at most 6 neighbors.
subroutine soca_field_stencil_interp(self, geom, fromto)
  class(soca_field), intent(inout) :: self
  class(soca_geom),    intent(in) :: geom   !< geometry
  character(len=4),     intent(in) :: fromto !< "u2h", "v2h"

  integer :: i, j
  real(kind=kind_real), allocatable :: val_tmp(:,:,:)
  real(kind=kind_real) :: val_max = 9e8_kind_real
  integer :: ij(2,6), sti, nn
  real(kind_real) :: lon_src(6), lat_src(6)
  real(kind=kind_real), allocatable :: val(:,:)
  real(kind=kind_real), allocatable :: lonsrc_local(:,:), latsrc_local(:,:)
  real(kind=kind_real), allocatable :: londst_local(:,:), latdst_local(:,:)
  real(kind=kind_real), allocatable :: masksrc_local(:,:), maskdst_local(:,:)

  ! Initialize temporary arrays
  allocate(val_tmp, mold=self%val)
  val_tmp = 0_kind_real

  ! Identify source and destination grids
  select case(fromto)
  case("vtoh")
     ! Horizontal interpolation: v-points to h-points
     allocate(lonsrc_local, mold=geom%lonv); lonsrc_local = geom%lonv
     allocate(latsrc_local, mold=geom%latv); latsrc_local = geom%latv
     allocate(masksrc_local, mold=geom%mask2dv);  masksrc_local = geom%mask2dv
     allocate(londst_local, mold=geom%lon);  londst_local = geom%lon
     allocate(latdst_local, mold=geom%lat);  latdst_local = geom%lat
     allocate(maskdst_local, mold=geom%mask2d);  maskdst_local = geom%mask2d

  case("utoh")
     ! Horizontal interpolation: u-points to h-points
     allocate(lonsrc_local, mold=geom%lonu); lonsrc_local = geom%lonu
     allocate(latsrc_local, mold=geom%latu); latsrc_local = geom%latu
     allocate(masksrc_local, mold=geom%mask2du);  masksrc_local = geom%mask2du
     allocate(londst_local, mold=geom%lon);  londst_local = geom%lon
     allocate(latdst_local, mold=geom%lat);  latdst_local = geom%lat
     allocate(maskdst_local, mold=geom%mask2d);  maskdst_local = geom%mask2d

  case default
     call abor1_ftn('soca_field::stencil_interp, option '//fromto//&
                    ' not implemented yet')

  end select

  ! Interpolate
  allocate(val(6,self%nz))
  do j = geom%jsc, geom%jec
     do i = geom%isc, geom%iec
        ! destination on land, skip
        if (maskdst_local(i,j) == 0_kind_real) cycle

        ! get the 6 or less src-point neighbors surrounding the (i,j) dst-point
        call soca_stencil_neighbors(fromto, i, j, ij)
        nn = 1
        val = 0_kind_real
        do sti = 1, 6
           ! source point on land, skip
           if (masksrc_local(ij(1,sti), ij(2,sti)) == 0_kind_real) cycle

           ! outcroping of layers, skip
           if (abs(self%val(ij(1,sti), ij(2,sti),1)) > val_max) cycle

           ! store the valid neighbors
           lon_src(nn) = lonsrc_local(ij(1,sti), ij(2,sti))
           lat_src(nn) = latsrc_local(ij(1,sti), ij(2,sti))
           val(nn,:) = self%val(ij(1,sti), ij(2,sti),:)
           nn = nn + 1
        end do
        nn = nn - 1

        ! val_tmp: interpolated val at (i,j) dst-point along layers
        if ( nn >=1 ) then
           call soca_stencil_interp(lon_src, lat_src, &
                                    londst_local(i,j), latdst_local(i,j), &
                                    val, val_tmp(i,j,:), nn)
        end if
     end do
  end do
  self%val = val_tmp

end subroutine soca_field_stencil_interp

! ------------------------------------------------------------------------------
!> Fill masked values
!!
!! Needed when reading fms history which can contain NaN's over land
subroutine soca_field_fill_masked(self, geom)
  class(soca_field), intent(inout) :: self
  type(soca_geom),      intent(in) :: geom

  integer :: i, j

  if (.not. associated(self%mask)) return
  do j = geom%jsc, geom%jec
    do i = geom%isc, geom%iec
      if (self%mask(i,j)==0) self%val(i,j,:) = self%metadata%fillvalue
    end do
  end do

end subroutine soca_field_fill_masked

! ------------------------------------------------------------------------------
!> Delete the soca_field object.
!!
!! \relates soca_fields_mod::soca_field
subroutine soca_field_delete(self)
  class(soca_field), intent(inout) :: self

  if (allocated(self%val)) deallocate(self%val)
end subroutine


! ------------------------------------------------------------------------------
! soca_fields subroutines
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> For a given list of field names, initialize the properties of those fields
!!
!! \param[in] vars: List of variables to initialize. They must be present in the
!!   configuration file used to create soca_fields_metadata_mod::soca_fields_metadata
!!
!! \throws abor1_ftn aborts if illegal grid or levels specified
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_init_vars(self, vars)
  class(soca_fields),         intent(inout) :: self
  character(len=:), allocatable, intent(in) :: vars(:)

  integer :: i, nz

  allocate(self%fields(size(vars)))
  do i=1,size(vars)
    self%fields(i)%name = trim(vars(i))

    ! get the field metadata parameters that are read in from a config file
    self%fields(i)%metadata = self%geom%fields_metadata%get(self%fields(i)%name)
    ! Set grid location and masks
    select case(self%fields(i)%metadata%grid)
    case ('h')
      self%fields(i)%lon => self%geom%lon
      self%fields(i)%lat => self%geom%lat
      if (self%fields(i)%metadata%masked) &
        self%fields(i)%mask => self%geom%mask2d
    case ('u')
      self%fields(i)%lon => self%geom%lonu
      self%fields(i)%lat => self%geom%latu
      if (self%fields(i)%metadata%masked) &
        self%fields(i)%mask => self%geom%mask2du
    case ('v')
        self%fields(i)%lon => self%geom%lonv
        self%fields(i)%lat => self%geom%latv
        if (self%fields(i)%metadata%masked) &
          self%fields(i)%mask => self%geom%mask2dv
    case default
      call abor1_ftn('soca_fields::create(): Illegal grid '// &
                     self%fields(i)%metadata%grid // &
                     ' given for ' // self%fields(i)%name)
    end select

    ! determine number of levels
    if (self%fields(i)%name == self%fields(i)%metadata%name_surface) then
      ! if this field is a surface getval, override the number of levels with 1
      nz = 1
    else
      select case(self%fields(i)%metadata%levels)
      case ('full_ocn')
        nz = self%geom%nzo
      case default
        read(self%fields(i)%metadata%levels, *) nz
      end select
    endif

    ! allocate space
    self%fields(i)%nz = nz
    allocate(self%fields(i)%val(&
      self%geom%isd:self%geom%ied, &
      self%geom%jsd:self%geom%jed, &
      nz ))

  end do
end subroutine


! ------------------------------------------------------------------------------
!> Create a new set of fields, allocate space for them, and initialize to zero
!!
!! \see soca_fields_init_vars
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_create(self, geom, vars, aFieldset)
  class(soca_fields),        intent(inout) :: self
  type(soca_geom),  pointer, intent(inout) :: geom !< geometry to associate with the fields
  type(oops_variables),      intent(in) :: vars !< list of field names to create
  type(atlas_fieldset),      intent(in) :: aFieldset

  character(len=:), allocatable :: vars_str(:)
  integer :: i

  self%afieldset = aFieldset

  ! make sure current object has not already been allocated
  if (allocated(self%fields)) &
    call abor1_ftn("soca_fields::create(): object already allocated")

  ! associate geometry
  self%geom => geom

  ! initialize the variable parameters
  allocate(character(len=1024) :: vars_str(vars%nvars()))
  do i=1,vars%nvars()
    vars_str(i) = trim(vars%variable(i))
  end do
  call soca_fields_init_vars(self, vars_str)

  ! set everything to zero
  call self%zeros()

end subroutine soca_fields_create


! ------------------------------------------------------------------------------
!> delete all the fields
!!
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_delete(self)
  class(soca_fields), intent(inout) :: self
  integer :: i

  ! clear the fields and nullify pointers
  nullify(self%geom)
  do i = 1, size(self%fields)
    call self%fields(i)%delete()
  end do
  deallocate(self%fields)

end subroutine


! ------------------------------------------------------------------------------
!> Copy the contents of \p rhs to \p self.
!!
!! \p self will be initialized with the variable names in \p rhs if
!! not already initialized.
!!
!! \see soca_fields_init_vars
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_copy(self, rhs)
  class(soca_fields), intent(inout) :: self
  class(soca_fields),  intent(in)    :: rhs !< fields to copy from

  character(len=:), allocatable :: vars_str(:)
  integer :: i
  type(soca_field), pointer :: rhs_fld

  ! initialize the variables based on the names in rhs
  if (.not. allocated(self%fields)) then
    self%geom => rhs%geom
    allocate(character(len=1024) :: vars_str(size(rhs%fields)))
    do i=1, size(vars_str)
      vars_str(i) = rhs%fields(i)%name
    end do
    call soca_fields_init_vars(self, vars_str)
  end if

  ! copy values from rhs to self, only if the variable exists
  !  in self
  do i=1,size(self%fields)
    call rhs%get(self%fields(i)%name, rhs_fld)
    call self%fields(i)%copy(rhs_fld)
  end do

end subroutine


! ------------------------------------------------------------------------------
!> Get a pointer to the soca_field with the given name.
!!
!! \note use soca_fields::has() if you need to check for optional fields
!! \throws abor1_ftn If no field exists with that name, the prorgam aborts
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_get(self, name, field)
  class(soca_fields), target, intent(in)  :: self
  character(len=*),           intent(in)  :: name !< name of field to find
  type(soca_field), pointer,  intent(out) :: field  !< a pointer to the resulting field

  integer :: i

  ! find the field with the given name
  do i=1,size(self%fields)
    if (trim(name) == self%fields(i)%name) then
      field => self%fields(i)
      return
    end if
  end do

  ! oops, the field was not found
  call abor1_ftn("soca_fields::get():  cannot find field "//trim(name))
end subroutine


! ------------------------------------------------------------------------------
!> Returns whether a field with the given name exists
!!
!! \relates soca_fields_mod::soca_fields
function soca_fields_has(self, name) result(res)
  class(soca_fields), intent(in) :: self
  character(len=*),   intent(in) :: name !< name of field to find

  logical :: res
  integer :: i

  res = .false.
  do i=1,size(self%fields)
    if (trim(name) == self%fields(i)%name) then
      res = .true.
      return
    end if
  end do
end function


! ------------------------------------------------------------------------------
!> Update the halo region of all fields.
!!
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_update_halos(self)
  class(soca_fields), intent(inout) :: self
  integer :: i

  do i=1,size(self%fields)
    call self%fields(i)%update_halo(self%geom)
  end do
end subroutine soca_fields_update_halos


! ------------------------------------------------------------------------------
!> Set the value of all fields to one.
!!
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_ones(self)
  class(soca_fields), intent(inout) :: self
  integer :: i

  do i = 1, size(self%fields)
    self%fields(i)%val = 1.0_kind_real
  end do

end subroutine soca_fields_ones


! ------------------------------------------------------------------------------
!> Reset the value of all fields to zero.
!!
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_zeros(self)
  class(soca_fields), intent(inout) :: self
  integer :: i

  do i = 1, size(self%fields)
    self%fields(i)%val = 0.0_kind_real
  end do

end subroutine soca_fields_zeros


! ------------------------------------------------------------------------------
!> read a set of fields from a file
!!
!! \param[in] f_conf : Configuration with the following parameters
!!    - "read_from_file" :
!!      - 0 = Invent the state
!!      - 1 = read state
!!      - 2 = (nothing??)
!!      - 3 = read increment
!!    - "remap_filename" : (optional) the filename containing "h" to perform the
!!      vertical remapping of these fields after they are loaded.
!!    - "date" : (required if read_from_file == 0)
!!    - "basename" : The common part of the path prepended to the following
!!       \c *_filename parameters
!!    - "ocn_filename" : ocean filename
!!    - "sfc_filename" : (optional) surface field filename
!!    - "ice_filename" : (optional) ice field filename
!!    - "wav_filename" : (optoinal) wave field filename
!!    - "bio_filename" : (optoinal) biochemistry field filename
!! \param[inout] vdate : If fields are being invented (read_from_file == 0),
!!    the \p vdate is used as the valid date of the fields. If the fields are
!!    being read in as a state (read_from_file == 1), \p vdate is set the the
!!    date from the files
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_read(self, f_conf, vdate)
  class(soca_fields), target, intent(inout) :: self
  type(fckit_configuration),  intent(in)    :: f_conf
  type(datetime),             intent(inout) :: vdate

  integer, parameter :: max_string_length=800
  character(len=max_string_length) :: ocn_filename, sfc_filename, ice_filename, wav_filename, bio_filename, filename
  character(len=:), allocatable :: basename, incr_filename
  integer :: iread = 0, id
  integer :: ii
  logical :: vert_remap=.false.
  character(len=max_string_length) :: remap_filename
  real(kind=kind_real), allocatable :: h_common(:,:,:)    !< layer thickness to remap to
  type(restart_file_type), target :: ocean_restart, sfc_restart, ice_restart, wav_restart, bio_restart
  type(restart_file_type) :: ocean_remap_restart
  type(restart_file_type), pointer :: restart
  integer :: idr
  integer :: isd, ied, jsd, jed
  integer :: isc, iec, jsc, jec
  integer :: i, j, k, nz, n
  type(remapping_CS)  :: remapCS
  character(len=:), allocatable :: str
  real(kind=kind_real), allocatable :: h_common_ij(:), hocn_ij(:), varocn_ij(:), varocn2_ij(:)
  logical :: read_sfc, read_ice, read_wav, read_bio, soca_io
  type(soca_field), pointer :: field, field2, hocn, mld, layer_depth
  type(oops_variables) :: seaice_categories_vars

  if ( f_conf%has("read_from_file") ) &
      call f_conf%get_or_die("read_from_file", iread)

  ! Get Indices for data domain and allocate common layer depth array
  isd = self%geom%isd ; ied = self%geom%ied
  jsd = self%geom%jsd ; jed = self%geom%jed

  ! Check if vertical remapping needs to be applied
  nz = self%geom%nzo
  if ( f_conf%has("remap_filename") ) then
     vert_remap = .true.
     call f_conf%get_or_die("remap_filename", str)
     remap_filename = str
     allocate(h_common(isd:ied,jsd:jed,nz))
     h_common = 0.0_kind_real

     ! Read common vertical coordinate from file
     call fms_io_init()
     idr = register_restart_field(ocean_remap_restart, remap_filename, 'h', h_common, &
          domain=self%geom%Domain%mpp_domain)
     call restore_state(ocean_remap_restart, directory='')
     call free_restart_type(ocean_remap_restart)
     call fms_io_exit()
  end if

  ! Create unit increment
  if ( f_conf%has("Identity") ) then
     call f_conf%get_or_die("Identity", id)
     if ( id==1 ) call self%ones()
     call f_conf%get_or_die("date", str)
     call datetime_set(str, vdate)
  end if

  ! TODO redo this to be generic

  ! iread = 1 (state) or 3 (increment): Read restart file
  if ((iread==1).or.(iread==3)) then
    if (self%has("sea_water_cell_thickness")) call self%get("sea_water_cell_thickness", hocn)
    ! filename for ocean
    call f_conf%get_or_die("basename", str)
    basename = str
    call f_conf%get_or_die("ocn_filename", str)
    ocn_filename = trim(basename) // trim(str)

    ! filename for ocn sfc
    read_sfc = .false.
    sfc_filename=""
    if ( f_conf%has("sfc_filename") ) then
      call f_conf%get_or_die("basename", str)
      basename = str
      call f_conf%get_or_die("sfc_filename", str)
      sfc_filename = trim(basename)//trim(str)
    end if

    ! filename for ice
    read_ice = .false.
    ice_filename=""
    if ( f_conf%has("ice_filename") ) then
      call f_conf%get_or_die("basename", str)
      basename = str
      call f_conf%get_or_die("ice_filename", str)
      ice_filename = trim(basename)//trim(str)
    end if

    ! filename for wav
    read_wav = .false.
    wav_filename=""
    if ( f_conf%has("wav_filename") ) then
      call f_conf%get_or_die("basename", str)
      basename = str
      call f_conf%get_or_die("wav_filename", str)
      wav_filename = trim(basename)//trim(str)
    end if

    ! filename for bio
    read_bio = .false.
    bio_filename=""
    if ( f_conf%has("bio_filename") ) then
      call f_conf%get_or_die("basename", str)
      basename = str
      call f_conf%get_or_die("bio_filename", str)
      bio_filename = trim(basename)//trim(str)
    end if

    call fms_io_init()

    seaice_categories_vars = oops_variables()
    ! built-in variables
    do i=1,size(self%fields)

      if(self%fields(i)%metadata%io_file == "CONSTANT") then
        self%fields(i)%val(:,:,:) = self%fields(i)%metadata%constant_value

      else if(self%fields(i)%metadata%io_file /= "") then
        ! which file are we reading from?
        select case(self%fields(i)%metadata%io_file)
        case ('ocn')
          filename = ocn_filename
          restart => ocean_restart
        case ('sfc')
          if (sfc_filename == "") cycle ! we have sfc fields, but no file to read from
          filename = sfc_filename
          restart => sfc_restart
          read_sfc = .true.
        case ('ice')
          filename = ice_filename
          restart => ice_restart
          ! note, read_ice is set below
        case ('wav')
          filename = wav_filename
          restart => wav_restart
          read_wav = .true.
        case ('bio')
          filename = bio_filename
          restart => bio_restart
          read_bio = .true.
        case default
          call abor1_ftn('read_file(): illegal io_file: '//self%fields(i)%metadata%io_file)
        end select

        ! check if the field has a category dimension and skip if it does
        if (self%fields(i)%metadata%io_file == "ice") then
          if (self%fields(i)%metadata%categories > 0 ) then
            ! check if the file was constructed by soca or comes from the CICE history
            ! The CICE history aggregates the category and level in 1 array
            ! The SOCA io considers categories to be separate variables and will index the naming
            if(file_exist(filename) .and. field_exist(filename, self%fields(i)%metadata%io_name)) then
              read_ice = .true.
            else
              call seaice_categories_vars%push_back(self%fields(i)%name)
              cycle
            end if
          else
            ! not a category field, treat normally
            read_ice = .true.
          end if
        end if

        ! setup to read
        if (self%fields(i)%nz == 1) then
          idr = register_restart_field(restart, filename, self%fields(i)%metadata%io_name, &
              self%fields(i)%val(:,:,1), domain=self%geom%Domain%mpp_domain)
        else
          idr = register_restart_field(restart, filename, self%fields(i)%metadata%io_name, &
              self%fields(i)%val(:,:,:), domain=self%geom%Domain%mpp_domain)
        end if
      end if
    end do

    call restore_state(ocean_restart, directory='')
    call free_restart_type(ocean_restart)
    if (read_sfc) then
      call restore_state(sfc_restart, directory='')
      call free_restart_type(sfc_restart)
    end if
    if (read_ice) then
      call restore_state(ice_restart, directory='')
      call free_restart_type(ice_restart)
    end if
    if (read_wav) then
      call restore_state(wav_restart, directory='')
      call free_restart_type(wav_restart)
    end if
    if (read_bio) then
      call restore_state(bio_restart, directory='')
      call free_restart_type(bio_restart)
    end if

    call fms_io_exit()

    ! read sea ice variables with categoriy and/or levels dimensions
    if (seaice_categories_vars%nvars() > 0) then
      call self%read_seaice(ice_filename, seaice_categories_vars)
    end if

    ! Change masked values
    do n=1,size(self%fields)
       field => self%fields(n)
       call field%fill_masked(self%geom)
    end do

    ! Update halo and return if reading increment
    if (iread==3) then !
       do n=1,size(self%fields)
         field => self%fields(n)
         call mpp_update_domains(field%val, self%geom%Domain%mpp_domain)
      end do
      return
   end if

    ! Indices for compute domain
    isc = self%geom%isc ; iec = self%geom%iec
    jsc = self%geom%jsc ; jec = self%geom%jec

    ! Remap layers if needed
    if (vert_remap) then

      ! output log of  what fields are going to be interpolated vertically
      if ( self%geom%f_comm%rank() == 0 ) then
        do n=1,size(self%fields)
          if (.not. self%fields(n)%metadata%vert_interp) cycle
          call oops_log%info("vertically remapping "//trim(self%fields(n)%name))
        end do
      end if

      allocate(h_common_ij(nz), hocn_ij(nz), varocn_ij(nz), varocn2_ij(nz))
      call initialize_remapping(remapCS,'PCM')
      do i = isc, iec
        do j = jsc, jec
          h_common_ij = h_common(i,j,:)
          hocn_ij = hocn%val(i,j,:)

          do n=1,size(self%fields)
            field => self%fields(n)
            ! TODO Vertical remapping is only valid if the field is on the tracer grid point.
            if (.not. field%metadata%vert_interp) cycle
            if (associated(field%mask) .and. field%mask(i,j).eq.1) then
               varocn_ij = field%val(i,j,:)
               call remapping_core_h(remapCS, nz, h_common_ij, varocn_ij,&
                      &nz, hocn_ij, varocn2_ij)
               field%val(i,j,:) = varocn2_ij
            else
               field%val(i,j,:) = 0.0_kind_real
            end if
          end do
        end do
      end do
      hocn%val = h_common
      deallocate(h_common_ij, hocn_ij, varocn_ij, varocn2_ij)
      call end_remapping(remapCS)
    end if

    ! Initialize mid-layer depth from layer thickness
    ! TODO, this shouldn't live here, it should be part of the variable change class only
    if (self%has("sea_water_depth")) then
      call self%get("sea_water_depth", layer_depth)
        layer_depth%val = 0.5 * hocn%val
        do k = 2, hocn%nz
          layer_depth%val(:,:,k) = layer_depth%val(:,:,k) + sum(hocn%val(:,:,1:k-1), dim=3)
        end do
    end if

    ! Compute mixed layer depth TODO: Move somewhere else ...
    if (self%has("ocean_mixed_layer_thickness") .and. self%has("sea_water_depth")) then
      call self%get("sea_water_potential_temperature", field)
      call self%get("sea_water_salinity", field2)
      call self%get("ocean_mixed_layer_thickness", mld)
      mld%val = 0.0
      do i = isc, iec
        do j = jsc, jec
            if (self%geom%mask2d(i,j)==0) cycle

            mld%val(i,j,1) = soca_mld(&
                &field2%val(i,j,:),&
                &field%val(i,j,:),&
                &layer_depth%val(i,j,:),&
                &self%geom%lon(i,j),&
                &self%geom%lat(i,j))
        end do
      end do
    end if

    ! Update halo
    do n=1,size(self%fields)
      field => self%fields(n)
      call mpp_update_domains(field%val, self%geom%Domain%mpp_domain)
    end do

    ! Set vdate if reading state
    if (iread==1) then
      call f_conf%get_or_die("date", str)
      call datetime_set(str, vdate)
    end if

    return
  end if

end subroutine soca_fields_read

! Populate an empty oop_variable instance with the unique CICE variables
subroutine get_cice_vars(self, cice_vars, ncat, nlev, cice_vars_type)
  type(soca_fields), intent(inout) :: self
  type(oops_variables), intent(in) :: cice_vars
  integer, intent(out) :: ncat, nlev
  character(len=5), intent(in) :: cice_vars_type

  integer :: i, levels

  select case (trim(cice_vars_type))
    case ("dynam")
      ! get the variables with a category dimension only (dynamic variables)
      nlev = 1
      do i=1,size(self%fields)
        if (self%fields(i)%metadata%io_file == "ice" .and.&
            & .not. cice_vars%has(self%fields(i)%metadata%io_sup_name)) then
          if (self%fields(i)%metadata%levels == '1' .and. self%fields(i)%metadata%categories > 0) then
            call cice_vars%push_back(self%fields(i)%metadata%io_sup_name)
          end if
          ncat = self%fields(i)%metadata%categories
        end if
      end do
    case ("therm")
      ! get the variables with category and level dimensions (thermodynamic variables)
      nlev = -1
      do i=1,size(self%fields)
        if (self%fields(i)%metadata%io_file == "ice" .and.&
            & .not. cice_vars%has(self%fields(i)%metadata%io_sup_name)) then
          read(self%fields(i)%metadata%levels, *) levels
          if (levels > 1 .and. self%fields(i)%metadata%categories > 0) then
            call cice_vars%push_back(self%fields(i)%metadata%io_sup_name)
            ncat = self%fields(i)%metadata%categories
            nlev = levels
          end if
        end if
      end do

    case default
      ! abort here
  end select

end subroutine get_cice_vars


subroutine soca_fields_read_seaice(self, filename, seaice_categories_vars)
  class(soca_fields), intent(inout) :: self
  character(800), intent(in) :: filename  !TODO: there's probably a better way to do this
  type(oops_variables), intent(in) :: seaice_categories_vars

  type(oops_variables) :: cice_vars_cats, cice_vars_cats_levs
  type(restart_file_type) :: restart

  integer :: i, ncat, icelevs, snowlevs, idr, cnt, io_index
  real(kind=kind_real), allocatable :: tmp3d(:,:,:,:), tmp4d(:,:,:,:,:)

  ! check what cice variables with category dimension need to be read
  cice_vars_cats = oops_variables()  ! used to store the unique cice io variables with a category dimension
  call get_cice_vars(self, cice_vars_cats, ncat, icelevs, "dynam")

  ! read the cice variables with category dimension only
  if (cice_vars_cats%nvars() > 0) then
    allocate(tmp3d(self%geom%isd:self%geom%ied,self%geom%jsd:self%geom%jed,ncat,cice_vars_cats%nvars()))
    tmp3d = 0.0_kind_real
    call fms_io_init()
    do i=1,cice_vars_cats%nvars()
      idr = register_restart_field(restart, filename, cice_vars_cats%variable(i), &
                         tmp3d(:,:,:,i), domain=self%geom%Domain%mpp_domain)
    end do
    call restore_state(restart, directory='')
    call free_restart_type(restart)
    call fms_io_exit()

    ! copy the variable into the corresponding field
    cnt = 1
    do i = 1, size(self%fields)
      if (self%fields(i)%metadata%io_file == "ice" .and.&
         &self%fields(i)%metadata%levels == '1' .and.&
         &self%fields(i)%metadata%categories > 0) then

        ! get the index of cice_vars that correspond to the io_sup_name
        io_index = cice_vars_cats%find(self%fields(i)%metadata%io_sup_name)
        self%fields(i)%val(:,:,1) = tmp3d(:,:,self%fields(i)%metadata%category,io_index)
      end if
    end do
  end if

  ! check what cice variables with category and level dimension need to be read
  cice_vars_cats_levs = oops_variables()  ! used to store the unique cice io variables with category and level dimensions
  call get_cice_vars(self, cice_vars_cats_levs, ncat, icelevs, "therm")

  ! read the seaice (not snow) cice variables with category and level dimensions
  if (cice_vars_cats_levs%nvars() > 0) then
    allocate(tmp4d(self%geom%isd:self%geom%ied,self%geom%jsd:self%geom%jed,icelevs,&
    &ncat,cice_vars_cats_levs%nvars()))
    tmp4d = 0.0_kind_real
    call fms_io_init()
    do i=1,cice_vars_cats_levs%nvars()
      idr = register_restart_field(restart, filename, cice_vars_cats_levs%variable(i), &
                         tmp4d(:,:,:,:,i), domain=self%geom%Domain%mpp_domain)
    end do
    call restore_state(restart, directory='')
    call free_restart_type(restart)
    call fms_io_exit()

    ! copy the variable into the corresponding field
    cnt = 1
    do i = 1, size(self%fields)
      if (self%fields(i)%metadata%io_file == "ice" .and.&
         &self%fields(i)%nz > 1 .and.&
         &self%fields(i)%metadata%categories > 0) then

        ! get the index of cice_vars that correspond to the io_sup_name
        io_index = cice_vars_cats_levs%find(self%fields(i)%metadata%io_sup_name)
        self%fields(i)%val(:,:,:) = tmp4d(:,:,:,self%fields(i)%metadata%category,io_index)
      end if
    end do
  end if
end subroutine soca_fields_read_seaice

! ------------------------------------------------------------------------------
!> Make sure two sets of fields are the same shape (same variables, same resolution)
!!
!! \throws abor1_ftn aborts if two fields are not congruent.
!!
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_check_congruent(self, rhs)
  class(soca_fields), intent(in) :: self
  class(soca_fields), intent(in) :: rhs !< other fields to check for congruency

  integer :: i, j

  ! number of fields should be the same
  if (size(self%fields) /= size(rhs%fields)) &
    call abor1_ftn("soca_fields: contains different number of fields")

  ! each field should match (name, size, shape)
  do i=1,size(self%fields)
    if (self%fields(i)%name /= rhs%fields(i)%name) &
      call abor1_ftn("soca_fields: field have different names")
    do j = 1, size(shape(self%fields(i)%val))
      if (size(self%fields(i)%val, dim=j) /= size(rhs%fields(i)%val, dim=j) ) then
        call abor1_ftn("soca_fields: field '"//self%fields(i)%name//"' has different dimensions")
      end if
    end do
  end do
end subroutine soca_fields_check_congruent


! ------------------------------------------------------------------------------
!> make sure two sets of fields are the same shape for fields they have in common
!!
!! \throws abor1_ftn aborts if \p self is not a subset of \p rhs
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_check_subset(self, rhs)
  class(soca_fields), intent(in) :: self
  class(soca_fields), intent(in) :: rhs !< other field that \p self should be subset of

  type(soca_field), pointer :: fld
  integer :: i, j

  ! each field should match (name, size, shape)
  do i=1,size(self%fields)
    if (.not. rhs%has(self%fields(i)%name)) &
      call abor1_ftn("soca_fields: self is not a subset of rhs")
    call rhs%get(self%fields(i)%name, fld)
    do j = 1, size(shape(fld%val))
      if (size(self%fields(i)%val, dim=j) /= size(fld%val, dim=j) ) then
        call abor1_ftn("soca_fields: field '"//self%fields(i)%name//"' has different dimensions")
      end if
    end do
  end do
end subroutine soca_fields_check_subset


! ------------------------------------------------------------------------------
!> Save soca fields to file using fms write_data
!!
!! \param[in] filename : The name of the file to save to
!!
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_write_file(self, filename)
  class(soca_fields),  intent(in) :: self
  character(len=*),   intent(in) :: filename

  integer :: ii

  call fms_io_init()
  call set_domain( self%geom%Domain%mpp_domain )

  ! write out all fields
  do ii = 1, size(self%fields)
    call write_data( filename, self%fields(ii)%name, self%fields(ii)%val(:,:,:), self%geom%Domain%mpp_domain)
  end do


  call fms_io_exit()
end subroutine soca_fields_write_file


! ------------------------------------------------------------------------------
!> Save soca fields in a restart format
!!
!! TODO this can be generalized even more
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_write_rst(self, f_conf, vdate)
  class(soca_fields), target, intent(inout) :: self      !< Fields
  type(fckit_configuration),  intent(in)    :: f_conf   !< Configuration
  type(datetime),             intent(inout) :: vdate    !< DateTime

  integer, parameter :: max_string_length=800
  character(len=max_string_length) :: ocn_filename, sfc_filename, ice_filename, wav_filename, bio_filename, filename
  type(restart_file_type), target :: ocean_restart, sfc_restart, ice_restart, wav_restart, bio_restart
  type(restart_file_type), pointer :: restart
  integer :: idr, i
  type(soca_field), pointer :: field
  logical :: write_sfc, write_ice, write_wav, write_bio, date_cols

  write_ice = .false.
  write_sfc = .false.
  write_wav = .false.
  write_bio = .false.
  call fms_io_init()

  ! Get date IO format (colons or not?)
  date_cols = .true.
  if (f_conf%has("date colons")) then
    call f_conf%get_or_die("date colons", date_cols)
  end if

  ! filenames
  ocn_filename = soca_genfilename(f_conf,max_string_length,vdate,date_cols,"ocn")
  sfc_filename = soca_genfilename(f_conf,max_string_length,vdate,date_cols,"sfc")
  ice_filename = soca_genfilename(f_conf, max_string_length,vdate,date_cols,"ice")
  wav_filename = soca_genfilename(f_conf, max_string_length,vdate,date_cols,"wav")
  bio_filename = soca_genfilename(f_conf, max_string_length,vdate,date_cols,"bio")

  ! built in variables
  do i=1,size(self%fields)
    field => self%fields(i)
    call field%fill_masked(self%geom)
    if (len_trim(field%metadata%io_file) /= 0) then
      ! which file are we writing to
      select case(field%metadata%io_file)
      case ('ocn')
        filename = ocn_filename
        restart => ocean_restart
      case ('sfc')
        filename = sfc_filename
        restart => sfc_restart
        write_sfc = .true.
      case ('ice')
        filename = ice_filename
        restart => ice_restart
        write_ice = .true.
      case ('wav')
        filename = wav_filename
        restart => wav_restart
        write_wav = .true.
      case ('bio')
        filename = bio_filename
        restart => bio_restart
        write_bio = .true.
      case default
        call abor1_ftn('soca_write_restart(): illegal io_file: '//field%metadata%io_file)
      end select

      ! write
      if (field%nz == 1) then
        idr = register_restart_field( restart, filename, field%metadata%io_name, &
          field%val(:,:,1), domain=self%geom%Domain%mpp_domain)
      else
        idr = register_restart_field( restart, filename, field%metadata%io_name, &
        field%val(:,:,:), domain=self%geom%Domain%mpp_domain)
      end if
    end if
  end do

  ! write out and cleanup
  call save_restart(ocean_restart, directory='')
  call free_restart_type(ocean_restart)
  if (write_sfc) then
    call save_restart(sfc_restart, directory='')
    call free_restart_type(sfc_restart)
  end if
  if (write_ice) then
    call save_restart(ice_restart, directory='')
    call free_restart_type(ice_restart)
  end if
  if (write_wav) then
    call save_restart(wav_restart, directory='')
    call free_restart_type(wav_restart)
  end if
  if (write_bio) then
    call save_restart(bio_restart, directory='')
    call free_restart_type(bio_restart)
 endif
  call fms_io_exit()

end subroutine soca_fields_write_rst

! ------------------------------------------------------------------------------
!> Interpolates from uv-points location to h-points.
!!
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_tohpoints(self)
  class(soca_fields), intent(inout) :: self !< self

  integer :: i
  real(kind=kind_real), allocatable :: val(:,:,:)
  real(kind=kind_real), pointer :: lon_out(:,:) => null()
  real(kind=kind_real), pointer :: lat_out(:,:) => null()
  character(len=4) :: fromto

  ! Associate lon_out and lat_out with the h-grid
  lon_out => self%geom%lon
  lat_out => self%geom%lat

  ! Apply interpolation to all fields, when necessary
  do i=1,size(self%fields)
    ! Check if already on h-points
    if (self%fields(i)%metadata%grid == 'h') cycle

    ! Interpolate to different location of the stencil
    fromto = self%fields(i)%metadata%grid//'toh'
    call self%fields(i)%stencil_interp(self%geom, fromto)
    call self%fields(i)%update_halo(self%geom)

    ! Update grid location to h-points
    self%fields(i)%metadata%grid = 'h'
    self%fields(i)%lon => self%geom%lon
    self%fields(i)%lat => self%geom%lat
 end do

end subroutine soca_fields_tohpoints

! ------------------------------------------------------------------------------
!> update fields, using list of variables the method removes fields not in the
!! list and allocates fields in the list but not allocated
!!
!! \see soca_fields_serialize
!! \relates soca_fields_mod::soca_fields

subroutine soca_fields_update_fields(self, vars)

  class(soca_fields),   intent(inout) :: self
  type(oops_variables), intent(in)    :: vars  ! New variable the field should have

  type(soca_fields) :: tmp_fields
  type(soca_field), pointer :: field
  integer :: f

  ! create new fields
  call tmp_fields%create(self%geom, vars, self%aFieldset)

  ! copy over where already existing
  do f = 1, size(tmp_fields%fields)
    if (self%has(tmp_fields%fields(f)%name)) then
      call self%get(tmp_fields%fields(f)%name, field)
      call tmp_fields%fields(f)%copy(field)
    end if
  end do

  ! move ownership of fields from tmp to self
  call move_alloc(tmp_fields%fields, self%fields)

  call self%update_metadata()

end subroutine soca_fields_update_fields

! ------------------------------------------------------------------------------
! Internal module functions/subroutines
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Generate filename (based on oops/qg)
!!
!! The configuration \p f_conf is expected to provide the following
!! - "datadir" : the directory the filenames should be prefixed with
!! - "exp" : experiment name
!! - "type" : one of "fc", "an", "incr", "ens"
!! - "member" : required only if "type == ens"
function soca_genfilename(f_conf,length,vdate,date_cols,domain_type)
  type(fckit_configuration),  intent(in) :: f_conf
  integer,                    intent(in) :: length
  type(datetime),             intent(in) :: vdate
  logical,                    intent(in) :: date_cols  !< Date written with colons or not
  character(len=3), optional, intent(in) :: domain_type

  character(len=length)                  :: soca_genfilename
  character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
       & prefix, mmb
  type(datetime) :: rdate
  type(duration) :: step
  integer lenfn
  character(len=:), allocatable :: str

  call f_conf%get_or_die("datadir", str)
  fdbdir = str
  call f_conf%get_or_die("exp", str)
  expver = str
  call f_conf%get_or_die("type", str)
  typ = str

  if (present(domain_type)) then
     expver = trim(domain_type)//"."//expver
  else
     expver = "ocn.ice."//expver
  end if
  if (typ=="ens") then
     call f_conf%get_or_die("member", str)
     mmb = str
     lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
     prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ) // "." // TRIM(mmb)
  else
     lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
     prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ)
  endif

  if (typ=="fc" .or. typ=="ens") then
     call f_conf%get_or_die("date", str)
     referencedate = str
     if (date_cols) then
       call datetime_create(trim(referencedate),rdate)
       call datetime_diff(vdate,rdate,step)
       call duration_to_string(step,sstep)
     else
       call datetime_create(trim(referencedate),rdate)
       call datetime_to_string_io(rdate,referencedate)
       call datetime_diff(vdate,rdate,step)
       call duration_to_string(step,sstep)
     endif
     lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
     soca_genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
  endif

  if (typ=="an" .or. typ=="incr") then
     if (date_cols) then
       call datetime_to_string(vdate,validitydate)
     else
       call datetime_to_string_io(vdate,validitydate)
     endif
     lenfn = lenfn + 1 + LEN_TRIM(validitydate)
     soca_genfilename = TRIM(prefix) // "." // TRIM(validitydate)
  endif

  if (lenfn>length) &
       & call abor1_ftn("fields:genfilename: filename too long")

   if ( allocated(str) ) deallocate(str)

end function soca_genfilename

! ------------------------------------------------------------------------------
! Copy the data from the internal "fields" to the atlas fieldset.
! We also need to make sure the metadata on the fieldset is set correctly
subroutine soca_fields_sync_to_atlas(self)
  class(soca_fields),   intent(inout)    :: self

  type(atlas_field) :: afield
  integer :: v, n, i, j
  type(atlas_metadata) :: meta
  real(kind=kind_real), pointer :: real_ptr(:,:)

  ! copy from internal fortran fields to atlas fieldset
  do v=1,size(self%fields)
    ! get/create field
    if (self%afieldset%has_field( self%fields(v)%name)) then
      afield = self%afieldset%field( self%fields(v)%name)
    else
      afield = self%geom%functionspace%create_field( &
        name= self%fields(v)%name, kind=atlas_real(kind_real), levels= self%fields(v)%nz)
      call self%afieldset%add(afield)
    end if

    ! create and fill field
    call afield%data(real_ptr)
    real_ptr = 0.0_kind_real  ! set all points to zero, overwrite owned values below
    do j=self%geom%jsc,self%geom%jec
      do i=self%geom%isc,self%geom%iec
        real_ptr(:, self%geom%atlas_ij2idx(i,j)) =  self%fields(v)%val(i,j,:)
      end do
    end do
    call afield%set_dirty(.true.)  ! indicate halo values are out-of-date

    call afield%final()
  end do

  ! update that atlas fieldset metadata, if needed
  call self%update_metadata()
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_fields_sync_from_atlas(self)
  class(soca_fields),   intent(inout)    :: self

  integer :: v, i, j
  type(atlas_field) :: afield
  real(kind=kind_real), pointer :: real_ptr(:,:)
  type(soca_field), pointer :: field

  ! TODO, remove fields that no longer exist?

  do v = 1, self%afieldset%size()
    afield = self%afieldset%field(v)

    if (.not. self%has(afield%name())) then
      cycle
      ! IS THIS A BUG, why are there fields in the internal field that are not part of the atlas fieldset????
      ! call abor1_ftn('fields_sync_from_atlas: variable '//trim(afield%name())//' not found in fields')
    endif

    call self%get(afield%name(), field)
    call afield%data(real_ptr)
    do j = self%geom%jsc, self%geom%jec
      do i = self%geom%isc, self%geom%iec
        field%val(i,j,:) = real_ptr(:, self%geom%atlas_ij2idx(i,j))
      end do
    end do
    call field%update_halo(self%geom)
  call afield%final()
  end do


  call self%update_metadata()
end subroutine


! ------------------------------------------------------------------------------
! update the metadata in the atlas fieldset based on fields metadata
! TODO this should probably just be combined with the update fields method
subroutine soca_fields_update_metadata(self)
  class(soca_fields), intent(inout) :: self
  integer :: n
  type(atlas_field) :: afield
  type(atlas_metadata) :: ameta
  type(soca_field_metadata) :: metadata

  do n =1, self%afieldset%size()
    afield = self%afieldset%field(n)
    ameta = afield%metadata()
    metadata = self%geom%fields_metadata%get(afield%name())

    call ameta%set('interp_type', 'default')
    if (metadata%masked) then
      call ameta%set('mask', 'interp_mask')
    end if
  end do
end subroutine

end module soca_fields_mod
