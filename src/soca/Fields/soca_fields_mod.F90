! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


!> Handle fields for the model.
!!
!! soca_fields represents a state or increment, and contains one or more
!! soca_field instances for each of the fields. The metadata associated
!! with a given field is stored in soca_fields_metadata_mod::soca_fields_metadata
module soca_fields_mod

! JEDI modules
use datetime_mod, only: datetime, datetime_set, datetime_to_string, &
                        datetime_create, datetime_diff
use duration_mod, only: duration, duration_to_string
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_min, fckit_mpi_max, fckit_mpi_sum
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables
use tools_const, only: deg2rad

! MOM6 / FMS modules
use fms_io_mod, only: fms_io_init, fms_io_exit, register_restart_field, &
                      restart_file_type, restore_state, free_restart_type, save_restart
use fms_mod,    only: write_data, set_domain
use horiz_interp_mod, only : horiz_interp_type
use horiz_interp_spherical_mod, only : horiz_interp_spherical, horiz_interp_spherical_del, &
                                       horiz_interp_spherical_new
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h, &
                          end_remapping
use mpp_domains_mod, only : mpp_update_domains

! SOCA modules
use soca_fields_metadata_mod, only : soca_field_metadata
use soca_geom_mod, only : soca_geom
use soca_utils, only: soca_mld

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

  !> \copybrief soca_fields_add \see soca_fields_add
  procedure :: add      => soca_fields_add

  !> \copybrief soca_fields_axpy \see soca_fields_axpy
  procedure :: axpy     => soca_fields_axpy

  !> \copybrief soca_fields_dotprod \see soca_fields_dotprod
  procedure :: dot_prod => soca_fields_dotprod

  !> \copybrief soca_fields_gpnorm \see soca_fields_gpnorm
  procedure :: gpnorm   => soca_fields_gpnorm

  !> \copybrief soca_fields_mul \see soca_fields_mul
  procedure :: mul      => soca_fields_mul

  !> \copybrief soca_fields_sub \see soca_fields_sub
  procedure :: sub      => soca_fields_sub

  !> \copybrief soca_fields_ones \see soca_fields_ones
  procedure :: ones     => soca_fields_ones

  !> \copybrief soca_fields_zeros \see soca_fields_zeros
  procedure :: zeros    => soca_fields_zeros

  !> \}

  !> \name I/O
  !! \{

  !> \copybrief soca_fields_read \see soca_fields_read
  procedure :: read      => soca_fields_read

  !> \copybrief soca_fields_write_file \see soca_fields_write_file
  procedure :: write_file=> soca_fields_write_file

  !> \copybrief soca_fields_write_rst \see soca_fields_write_rst
  procedure :: write_rst => soca_fields_write_rst

  !> \}

  !> \name misc
  !! \{

  !> \copybrief soca_fields_update_halos \see soca_fields_update_halos
  procedure :: update_halos => soca_fields_update_halos

  !> \copybrief soca_fields_colocate \see soca_fields_colocate
  procedure :: colocate  => soca_fields_colocate
  !> \}

  !> \name serialization
  !! \{

  !> \copybrief soca_fields_serial_size \see soca_fields_serial_size
  procedure :: serial_size => soca_fields_serial_size

  !> \copybrief soca_fields_serialize \see soca_fields_serialize
  procedure :: serialize   => soca_fields_serialize

  !> \copybrief soca_fields_deserialize \see soca_fields_deserialize
  procedure :: deserialize => soca_fields_deserialize

  !> \}

  !> \copybrief soca_fields_update_fields \see soca_fields_update_fields
  procedure :: update_fields => soca_fields_update_fields

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
!> Perform spatial interpolation between two grids.
!!
!! Interpolation used is inverse distance weidghted, taking into
!! consideration the mask.
!! \param[in] geom: The geometry to interpolate to
!! \param[in] interp2d: interpolation object created by calling
!!     \c horiz_interp_spherical_new() in FMS
!!
!! \relates soca_fields_mod::soca_field
subroutine soca_field_stencil_interp(self, geom, interp2d)
  class(soca_field),     intent(inout) :: self
  type(soca_geom), pointer, intent(in) :: geom
  type(horiz_interp_type),  intent(in) :: interp2d

  integer :: k
  real(kind=kind_real), allocatable :: val(:,:,:)

  allocate(val, mold=self%val)
  val = self%val
  do k = 1, self%nz
     call horiz_interp_spherical(interp2d, &
          & val(geom%isd:geom%ied, geom%jsd:geom%jed,k), &
          & self%val(geom%isc:geom%iec, geom%jsc:geom%jec,k))
  end do
  call self%update_halo(geom)
end subroutine soca_field_stencil_interp


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
    if (self%fields(i)%name == self%fields(i)%metadata%getval_name_surface) then
      ! if this field is a surface getval, override the number of levels with 1
      nz = 1
    else
      select case(self%fields(i)%metadata%levels)
      case ('full_ocn')
        nz = self%geom%nzo
      case ('1') ! TODO, generalize to work with any number?
        nz = 1
      case default
        call abor1_ftn('soca_fields::create(): Illegal levels '//self%fields(i)%metadata%levels// &
                       ' given for ' // self%fields(i)%name)
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
subroutine soca_fields_create(self, geom, vars)
  class(soca_fields),        intent(inout) :: self
  type(soca_geom),  pointer, intent(inout) :: geom !< geometry to associate with the fields
  type(oops_variables),      intent(in) :: vars !< list of field names to create

  character(len=:), allocatable :: vars_str(:)
  integer :: i

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
  write(6,*) 'size of fields is ',size(self%fields)
  do i=1,size(self%fields)
    write(6,*) 'name of field ',i,' is ',trim(self%fields(i)%name)
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
!> Add two sets of fields together
!!
!! \f$ self = self + rhs \f$
!!
!! \throws abor1_ftn aborts if two fields are not congruent
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_add(self, rhs)
  class(soca_fields), intent(inout) :: self
  class(soca_fields),     intent(in) :: rhs !< other field to add
  integer :: i

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! add
  do i=1,size(self%fields)
    self%fields(i)%val = self%fields(i)%val + rhs%fields(i)%val
  end do
end subroutine soca_fields_add


! ------------------------------------------------------------------------------
!> subtract two sets of fields
!!
!! \f$ self = self - rhs \f$
!!
!! \throws abor1_ftn aborts if two fields are not congruent
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_sub(self, rhs)
  class(soca_fields), intent(inout) :: self
  class(soca_fields),     intent(in) :: rhs !< other field to subtract
  integer :: i

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! subtract
  do i=1,size(self%fields)
    self%fields(i)%val = self%fields(i)%val - rhs%fields(i)%val
  end do
end subroutine soca_fields_sub


! ------------------------------------------------------------------------------
!> Multiply a set of fields by a constant.
!!
!! \f$ self = zz * self \f$
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_mul(self, zz)
  class(soca_fields), intent(inout) :: self
  real(kind=kind_real),  intent(in) :: zz !< the constant by which to multipy the field
  integer :: i

  do i=1,size(self%fields)
    self%fields(i)%val = zz * self%fields(i)%val
  end do
end subroutine soca_fields_mul


! ------------------------------------------------------------------------------
!> Add two fields (multiplying the rhs first)
!!
!! \f$self = self + zz * rhs\f$
!!
!! \throws abor1_ftn aborts if \p is not a subset of \rhs
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_axpy(self, zz, rhs)
  class(soca_fields), target, intent(inout) :: self
  real(kind=kind_real),       intent(in)    :: zz !< constant by which to multiply other rhs
  class(soca_fields),         intent(in)    :: rhs !< other field to add

  type(soca_field), pointer :: f_rhs, f_lhs
  integer :: i

  ! make sure fields are correct shape
  call self%check_subset(rhs)

  do i=1,size(self%fields)
    f_lhs => self%fields(i)
    if (.not. rhs%has(f_lhs%name)) cycle
    call rhs%get(f_lhs%name, f_rhs)
    f_lhs%val = f_lhs%val + zz *f_rhs%val
  end do
end subroutine soca_fields_axpy


! ------------------------------------------------------------------------------
!> Calculate the global dot product of two sets of fields.
!!
!! \throws abor1_ftn aborts if two fields are not congruent
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_dotprod(self, rhs, zprod)
  class(soca_fields), target, intent(in)  :: self
  class(soca_fields), target, intent(in)  :: rhs !< field 2 of dot product
  real(kind=kind_real),       intent(out) :: zprod !< The resulting dot product

  real(kind=kind_real) :: local_zprod
  integer :: ii, jj, kk, n
  type(soca_field), pointer :: field1, field2

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! loop over (almost) all fields
  local_zprod = 0.0_kind_real
  do n=1,size(self%fields)
    field1 => self%fields(n)
    field2 => rhs%fields(n)

    ! add the given field to the dot product (only using the compute domain)
    do ii = self%geom%isc, self%geom%iec
      do jj = self%geom%jsc, self%geom%jec
        ! masking
        if (associated(field1%mask)) then
          if (field1%mask(ii,jj) < 1) cycle
        endif

        ! add to dot product
        do kk=1,field1%nz
          local_zprod = local_zprod + field1%val(ii,jj,kk) * field2%val(ii,jj,kk)
        end do
      end do
    end do
  end do

  ! Get global dot product
  call self%geom%f_comm%allreduce(local_zprod, zprod, fckit_mpi_sum())
end subroutine soca_fields_dotprod


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
  character(len=max_string_length) :: ocn_filename, sfc_filename, ice_filename, wav_filename, filename
  character(len=:), allocatable :: basename, incr_filename
  integer :: iread = 0
  integer :: ii
  logical :: vert_remap=.false.
  character(len=max_string_length) :: remap_filename
  real(kind=kind_real), allocatable :: h_common(:,:,:)    !< layer thickness to remap to
  type(restart_file_type), target :: ocean_restart, sfc_restart, ice_restart, wav_restart
  type(restart_file_type) :: ocean_remap_restart
  type(restart_file_type), pointer :: restart
  integer :: idr
  integer :: isd, ied, jsd, jed
  integer :: isc, iec, jsc, jec
  integer :: i, j, nz, n
  type(remapping_CS)  :: remapCS
  character(len=:), allocatable :: str
  real(kind=kind_real), allocatable :: h_common_ij(:), hocn_ij(:), varocn_ij(:), varocn2_ij(:)
  logical :: read_sfc, read_ice, read_wav
  type(soca_field), pointer :: field, field2, hocn, mld, layer_depth

  if ( f_conf%has("read_from_file") ) &
      call f_conf%get_or_die("read_from_file", iread)

  call self%get("hocn", hocn)

  ! Get Indices for data domain and allocate common layer depth array
  isd = self%geom%isd ; ied = self%geom%ied
  jsd = self%geom%jsd ; jed = self%geom%jed

  ! Check if vertical remapping needs to be applied
  nz = hocn%nz
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

  ! iread = 0: Invent state
  if (iread==0) then
     call self%zeros()
     call f_conf%get_or_die("date", str)
     call datetime_set(str, vdate)
  end if

  ! TODO redo this to be generic

  ! iread = 1 (state) or 3 (increment): Read restart file
  if ((iread==1).or.(iread==3)) then

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

    call fms_io_init()

    ! built-in variables
    do i=1,size(self%fields)
      if(self%fields(i)%metadata%io_name /= "") then
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
          read_ice = .true.
        case ('wav')
          filename = wav_filename
          restart => wav_restart
          read_wav = .true.
        case default
          call abor1_ftn('read_file(): illegal io_file: '//self%fields(i)%metadata%io_file)
        end select

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

    call fms_io_exit()

    ! Indices for compute domain
    isc = self%geom%isc ; iec = self%geom%iec
    jsc = self%geom%jsc ; jec = self%geom%jec

    ! Initialize mid-layer depth from layer thickness
    if (self%has("layer_depth")) then
      call self%get("layer_depth", layer_depth)
      call self%geom%thickness2depth(hocn%val, layer_depth%val)
    end if

    ! Compute mixed layer depth TODO: Move somewhere else ...
    if (self%has("mld") .and. self%has("layer_depth")) then
      call self%get("tocn", field)
      call self%get("socn", field2)
      call self%get("mld", mld)
      do i = isc, iec
        do j = jsc, jec
            mld%val(i,j,1) = soca_mld(&
                &field2%val(i,j,:),&
                &field%val(i,j,:),&
                &layer_depth%val(i,j,:),&
                &self%geom%lon(i,j),&
                &self%geom%lat(i,j))
        end do
      end do
    end if

    ! Remap layers if needed
    if (vert_remap) then
      allocate(h_common_ij(nz), hocn_ij(nz), varocn_ij(nz), varocn2_ij(nz))
      call initialize_remapping(remapCS,'PCM')
      do i = isc, iec
        do j = jsc, jec
          h_common_ij = h_common(i,j,:)
          hocn_ij = hocn%val(i,j,:)

          do n=1,size(self%fields)
            field => self%fields(n)
            select case(field%name)
            ! TODO remove hardcoded variable names here
            ! TODO Add u and v. Remapping u and v will require interpolating h
            case ('tocn','socn')
              if (associated(field%mask) .and. field%mask(i,j).eq.1) then
                varocn_ij = field%val(i,j,:)
                call remapping_core_h(remapCS, nz, h_common_ij, varocn_ij,&
                      &nz, hocn_ij, varocn2_ij)
                field%val(i,j,:) = varocn2_ij
              else
                field%val(i,j,:) = 0.0_kind_real
              end if
            end select
          end do
        end do
      end do
      hocn%val = h_common
      deallocate(h_common_ij, hocn_ij, varocn_ij, varocn2_ij)
      call end_remapping(remapCS)
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


! ------------------------------------------------------------------------------
!> calculate global statistics for each field (min, max, average)
!!
!! \param[in] nf: The number of fields, should be equal to the size of
!!     soca_fields::fields
!! \param[out] pstat: a 2D array with shape (i,j). For each field index
!!     i is set as 0 = min, 1 = max, 2 = average, for j number of fields.
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_gpnorm(self, nf, pstat)
  class(soca_fields),      intent(in) :: self
  integer,                 intent(in) :: nf
  real(kind=kind_real),   intent(out) :: pstat(3, nf)

  logical :: mask(self%geom%isc:self%geom%iec, self%geom%jsc:self%geom%jec)
  real(kind=kind_real) :: ocn_count, local_ocn_count, tmp(3)
  integer :: jj, isc, iec, jsc, jec
  type(soca_field), pointer :: field

  ! Indices for compute domain
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  ! calculate global min, max, mean for each field
  do jj=1, size(self%fields)
    call self%get(self%fields(jj)%name, field)

    ! get the mask and the total number of grid cells
    if (.not. associated(field%mask)) then
       mask = .true.
     else
       mask = field%mask(isc:iec, jsc:jec) > 0.0
     end if
    local_ocn_count = count(mask)
    call self%geom%f_comm%allreduce(local_ocn_count, ocn_count, fckit_mpi_sum())

    ! calculate global min/max/mean
    call fldinfo(field%val(isc:iec,jsc:jec,:), mask, tmp)
    call self%geom%f_comm%allreduce(tmp(1), pstat(1,jj), fckit_mpi_min())
    call self%geom%f_comm%allreduce(tmp(2), pstat(2,jj), fckit_mpi_max())
    call self%geom%f_comm%allreduce(tmp(3), pstat(3,jj), fckit_mpi_sum())
    pstat(3,jj) = pstat(3,jj)/ocn_count
  end do
end subroutine soca_fields_gpnorm


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

  ! some other derived fields that should be written out
  call write_data( filename, "rossby_radius", self%geom%rossby_radius, self%geom%Domain%mpp_domain)

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
  character(len=max_string_length) :: ocn_filename, sfc_filename, ice_filename, wav_filename, filename
  type(restart_file_type), target :: ocean_restart, sfc_restart, ice_restart, wav_restart
  type(restart_file_type), pointer :: restart
  integer :: idr, i
  type(soca_field), pointer :: field
  logical :: write_sfc, write_ice, write_wav

  write_ice = .false.
  write_sfc = .false.
  write_wav = .false.
  call fms_io_init()

  ! filenames
  ocn_filename = soca_genfilename(f_conf,max_string_length,vdate,"ocn")
  sfc_filename = soca_genfilename(f_conf,max_string_length,vdate,"sfc")
  ice_filename = soca_genfilename(f_conf, max_string_length,vdate,"ice")
  wav_filename = soca_genfilename(f_conf, max_string_length,vdate,"wav")

  ! built in variables
  do i=1,size(self%fields)
    field => self%fields(i)
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
  call fms_io_exit()

end subroutine soca_fields_write_rst

! ------------------------------------------------------------------------------
!> Colocate by interpolating from one c-grid location to another.
!!
!! \warning only works on the "h" grid currently (not the "u" or "v" grid)
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_colocate(self, cgridlocout)
  class(soca_fields),    intent(inout) :: self !< self
  character(len=1),         intent(in) :: cgridlocout !< colocate to cgridloc (u, v or h)

  integer :: i, k
  real(kind=kind_real), allocatable :: val(:,:,:)
  real(kind=kind_real), pointer :: lon_out(:,:) => null()
  real(kind=kind_real), pointer :: lat_out(:,:) => null()
  type(soca_geom),  pointer :: g => null()
  type(horiz_interp_type) :: interp2d

  ! Associate lon_out and lat_out according to cgridlocout
  select case(cgridlocout)
  ! TODO: Test colocation to u and v grid
  !case ('u')
  !  lon_out => self%geom%lonu
  !  lat_out => self%geom%latu
  !case ('v')
  !  lon_out => self%geom%lonv
  !  lat_out => self%geom%latv
  case ('h')
    lon_out => self%geom%lon
    lat_out => self%geom%lat
  case default
    call abor1_ftn('soca_fields::colocate(): unknown c-grid location '// cgridlocout)
  end select

  ! Apply interpolation to all fields, when necessary
  do i=1,size(self%fields)

    ! Check if already colocated
    if (self%fields(i)%metadata%grid == cgridlocout) cycle

    ! Initialize fms spherical idw interpolation
     g => self%geom
     call horiz_interp_spherical_new(interp2d, &
       & real(deg2rad*self%fields(i)%lon(g%isd:g%ied,g%jsd:g%jed), 8), &
       & real(deg2rad*self%fields(i)%lat(g%isd:g%ied,g%jsd:g%jed), 8), &
       & real(deg2rad*lon_out(g%isc:g%iec,g%jsc:g%jec), 8), &
       & real(deg2rad*lat_out(g%isc:g%iec,g%jsc:g%jec), 8))

    ! Make a temporary copy of field
    if (allocated(val)) deallocate(val)
    allocate(val, mold=self%fields(i)%val)
    val = self%fields(i)%val

    ! Interpolate all levels
    do k = 1, self%fields(i)%nz
      call self%fields(i)%stencil_interp(self%geom, interp2d)
    end do

    ! Update c-grid location
    self%fields(i)%metadata%grid = cgridlocout
    select case(cgridlocout)
    ! TODO: Test colocation to u and v grid
    !case ('u')
    !  self%fields(i)%lon => self%geom%lonu
    !  self%fields(i)%lat => self%geom%latu
    !case ('v')
    !  self%fields(i)%lon => self%geom%lonv
    !  self%fields(i)%lat => self%geom%latv
    case ('h')
      self%fields(i)%lon => self%geom%lon
      self%fields(i)%lat => self%geom%lat
    end select

 end do
 call horiz_interp_spherical_del(interp2d)

end subroutine soca_fields_colocate


! ------------------------------------------------------------------------------
!> Number of elements to return in the serialized array
!!
!! \see soca_fields_serialize
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_serial_size(self, geom, vec_size)
  class(soca_fields),    intent(in)  :: self
  type(soca_geom),       intent(in)  :: geom !< todo remove, not needed?
  integer,               intent(out) :: vec_size !< resulting size of vector

  integer :: i

  ! Loop over fields
  vec_size = 0
  do i=1,size(self%fields)
    vec_size = vec_size + size(self%fields(i)%val)
  end do

end subroutine soca_fields_serial_size


! ------------------------------------------------------------------------------
!> Return the fields as a serialized array
!!
!! \see soca_fields_serial_size
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_serialize(self, geom, vec_size, vec)
  class(soca_fields),    intent(in)  :: self
  type(soca_geom),       intent(in)  :: geom  !< todo remove this, not needed?
  integer,               intent(in)  :: vec_size !< size of vector to return
  real(kind=kind_real),  intent(out) :: vec(vec_size) !< fields as a serialized vector

  integer :: index, i, nn

  ! Loop over fields, levels and horizontal points
  index = 1
  do i=1,size(self%fields)
    nn = size(self%fields(i)%val)
    vec(index:index+nn-1) = reshape(self%fields(i)%val, (/ nn /) )
    index = index + nn
  end do

end subroutine soca_fields_serialize

! ------------------------------------------------------------------------------
!> Deserialize, creating fields from a single serialized array
!!
!! \see soca_fields_serialize
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_deserialize(self, geom, vec_size, vec, index)
  class(soca_fields), intent(inout) :: self
  type(soca_geom),       intent(in)    :: geom !< todo remove this, not needed?
  integer,               intent(in)    :: vec_size !< size of \p vec
  real(kind=kind_real),  intent(in)    :: vec(vec_size) !< vector to deserialize
  integer,               intent(inout) :: index !< index in \p vec at which to start deserializing

  integer :: i, nn

  ! Loop over fields, levels and horizontal points
  do i=1,size(self%fields)
    nn = size(self%fields(i)%val)
    self%fields(i)%val = reshape(vec(index+1:index+1+nn), shape(self%fields(i)%val))
    index = index + nn
  end do

end subroutine soca_fields_deserialize

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
  call tmp_fields%create(self%geom, vars)

  ! copy over where already existing
  do f = 1, size(tmp_fields%fields)
    if (self%has(tmp_fields%fields(f)%name)) then
      call self%get(tmp_fields%fields(f)%name, field)
      call tmp_fields%fields(f)%copy(field)
    end if
  end do

  ! move ownership of fields from tmp to self
  call move_alloc(tmp_fields%fields, self%fields)

end subroutine soca_fields_update_fields

! ------------------------------------------------------------------------------
! Internal module functions/subroutines
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Calculate min/max/mean statistics for a given field, using a mask.
!!
!! \param[in] fld : the field to calculate the statistics on
!! \param[in] mask : statistics are only calculated where \p mask is \c .true.
!! \param[out] info : [0] = min, [1] = max, [2] = average
subroutine fldinfo(fld, mask, info)
  real(kind=kind_real),  intent(in) :: fld(:,:,:)
  logical,               intent(in) :: mask(:,:)
  real(kind=kind_real), intent(out) :: info(3)

  integer :: z
  real(kind=kind_real) :: tmp(3,size(fld, dim=3))

  ! calculate the min/max/sum separately for each masked level
  do z = 1, size(tmp, dim=2)
     tmp(1,z) = minval(fld(:,:,z), mask=mask)
     tmp(2,z) = maxval(fld(:,:,z), mask=mask)
     tmp(3,z) = sum(   fld(:,:,z), mask=mask) / size(fld, dim=3)
  end do

  ! then combine the min/max/sum over all levels
  info(1) = minval(tmp(1,:))
  info(2) = maxval(tmp(2,:))
  info(3) = sum(   tmp(3,:))
end subroutine fldinfo


! ------------------------------------------------------------------------------
!> Generate filename (based on oops/qg)
!!
!! The configuration \p f_conf is expected to provide the following
!! - "datadir" : the directory the filenames should be prefixed with
!! - "exp" : experiment name
!! - "type" : one of "fc", "an", "incr", "ens"
!! - "member" : required only if "type == ens"
function soca_genfilename (f_conf,length,vdate,domain_type)
  type(fckit_configuration),  intent(in) :: f_conf
  integer,                    intent(in) :: length
  type(datetime),             intent(in) :: vdate
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
     call datetime_to_string(vdate, validitydate)
     call datetime_create(TRIM(referencedate), rdate)
     call datetime_diff(vdate, rdate, step)
     call duration_to_string(step, sstep)
     lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
     soca_genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
  endif

  if (typ=="an" .or. typ=="incr") then
     call datetime_to_string(vdate, validitydate)
     lenfn = lenfn + 1 + LEN_TRIM(validitydate)
     soca_genfilename = TRIM(prefix) // "." // TRIM(validitydate)
  endif

  if (lenfn>length) &
       & call abor1_ftn("fields:genfilename: filename too long")

   if ( allocated(str) ) deallocate(str)

end function soca_genfilename


end module soca_fields_mod
