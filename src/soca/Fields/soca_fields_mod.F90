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
use iso_c_binding, only: c_ptr
use logger_mod
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables

! MOM6 / FMS modules
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h, &
                          end_remapping
use mpp_domains_mod, only : mpp_update_domains

! SOCA I/O
use soca_io_mod, only : soca_io_reader, soca_io_writer, &
                        soca_io_file_exists, soca_io_var_exists, &
                        soca_io_writers_commit_ensemble, &
                        soca_io_readers_commit_ensemble, &
                        soca_io_ensemble_root_pe, &
                        soca_io_ensemble_write_parallel

! SOCA modules
use soca_fields_metadata_mod, only : soca_field_metadata
use soca_geom_mod, only : soca_geom
use soca_utils, only: soca_mld
use soca_utils, only: soca_stencil_interp, soca_stencil_neighbors

implicit none
private

public :: soca_fields_write_ensemble
public :: soca_fields_read_ensemble


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

  !> \copybrief soca_fields_delete \see soca_fields_delete
  procedure :: delete => soca_fields_delete

  !> \}

  !> \name field getters/checkers
  !! \{

  !> \copybrief soca_fields_get \see soca_fields_get
  procedure :: get    => soca_fields_get

  !> \copybrief soca_fields_has \see soca_fields_has
  procedure :: has    => soca_fields_has

  !> \}

  !> \name math operators
  !! \{

  !> \copybrief soca_fields_ones \see soca_fields_ones
  procedure :: ones     => soca_fields_ones

  !> \}

  !> \name I/O
  !! \{

  !> \copybrief soca_fields_read \see soca_fields_read
  procedure :: read      => soca_fields_read
  procedure, private :: read_seaice => soca_fields_read_seaice

  !> \copybrief soca_fields_write_rst \see soca_fields_write_rst
  procedure :: write_rst => soca_fields_write_rst

  !> \}

  !> \name misc
  !! \{

  !> \copybrief soca_fields_tohpoints \see soca_fields_tohpoints
  procedure :: tohpoints  => soca_fields_tohpoints
  !> \}


  !> \copybrief soca_fields_update_fields \see soca_fields_update_fields
  procedure :: update_fields => soca_fields_update_fields

end type soca_fields

! ------------------------------------------------------------------------------
! Used to hold info when processing an atlas field
type varwrapper
  type(atlas_field) :: afield
  real(kind=kind_real), pointer :: adata(:,:)
  type(soca_field), pointer :: field
  real(kind=kind_real), allocatable :: data(:,:,:)
end type varwrapper

! ------------------------------------------------------------------------------
! Polymorphic pointer wrapper so soca_fields_*_ensemble can take heterogeneous
! arrays of soca_state and soca_increment (both extend soca_fields). Used only
! by the bulk I/O entrypoints; built on the Fortran side from arrays of C++
! registry keys.
type, public :: soca_fields_ptr_t
  class(soca_fields), pointer :: p => null()
end type soca_fields_ptr_t

! ------------------------------------------------------------------------------
! One writer's worth of varwrapper state in the ensemble write path. Per-writer
! arrays are nested under this type so a single 1-D allocatable indexed by
! "writer index" carries all the per-(member, domain) vars buffers.
type :: ensemble_vars_t
  type(varwrapper), allocatable :: vars(:)
end type ensemble_vars_t

! ------------------------------------------------------------------------------
! One member's worth of read state, carried across the soca_fields_read split
! into prepare (config + queue) and finalize (copy-to-atlas + post-processing).
! For the ensemble path, an array of plans is built up across members so that
! all per-domain readers can be flushed in one bulk soca_io_readers_commit_ensemble
! call, then each member's plan is finalized to drive its post-processing.
type :: read_plan_t
  integer :: iread = 0
  ! per-active-domain reader + scratch buffers. readers(ri) reads into
  ! domain_vars(ri)%vars(n)%data; same ri indexes both. Size 0 when
  ! iread is not 1 or 3 (the no-read paths).
  type(soca_io_reader),   allocatable :: readers(:)
  type(ensemble_vars_t),  allocatable :: domain_vars(:)
  ! vertical remap target (allocatable so its absence flags "no remap").
  ! target-ness comes from declaring plan instances with the target attribute.
  real(kind=kind_real),  allocatable :: h_common(:,:,:)
  ! CICE category-vars deferred-read state (handled in finalize via read_seaice).
  ! seaice_categories_vars is only initialized (oops_variables() empty_ctor)
  ! when iread==1 or 3; finalize matches that guard before destructing.
  type(oops_variables) :: seaice_categories_vars
  character(len=:),       allocatable :: ice_filename
  ! Probe results that drive post-processing
  logical :: compute_icethickness = .false.
  logical :: compute_snowthickness = .false.
end type read_plan_t


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
contains


! ------------------------------------------------------------------------------
! soca_field subroutines
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Perform spatial interpolation between adjacent grid point in the same stencil
!!
!! Interpolation used is inverse distance weidghted, taking into
!! consideration the mask and using at most 6 neighbors.
subroutine soca_field_stencil_interp(field, geom, fromto)
  real(kind=kind_real), allocatable, intent(inout) :: field(:,:,:)
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
  allocate(val_tmp, mold=field)
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
  allocate(val(6,size(field, 3)))
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
           if (abs(field(ij(1,sti), ij(2,sti),1)) > val_max) cycle

           ! store the valid neighbors
           lon_src(nn) = lonsrc_local(ij(1,sti), ij(2,sti))
           lat_src(nn) = latsrc_local(ij(1,sti), ij(2,sti))
           val(nn,:) = field(ij(1,sti), ij(2,sti),:)
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
  field = val_tmp

end subroutine soca_field_stencil_interp


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

  self%afieldset = aFieldset
  self%geom => geom

  ! make sure current object has not already been allocated
  if (allocated(self%fields)) &
    call abor1_ftn("soca_fields::create(): object already allocated")

  ! initialize the variables
  call self%update_fields(vars)

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
  deallocate(self%fields)
  call self%afieldset%final()

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
!> Set the value of all fields to one.
!!
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_ones(self)
  class(soca_fields), intent(inout) :: self
  type(atlas_field) :: field
  real(kind=kind_real), pointer :: fdata(:,:)
  integer :: i

  do i = 1, self%afieldset%size()
    field = self%afieldset%field(i)
    call field%data(fdata)
    fdata(:,:) = 1.0_kind_real
    call field%set_dirty(.false.)
  end do
  call field%final()

end subroutine soca_fields_ones


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

  type(read_plan_t), target :: plan
  integer :: r

  call soca_fields_read_prepare(self, f_conf, vdate, plan)

  ! Drain the per-domain readers serially in the single-state path.
  ! soca_fields_read_ensemble flushes all members' readers together via
  ! soca_io_readers_commit_ensemble instead.
  if (allocated(plan%readers)) then
    do r = 1, size(plan%readers)
      call plan%readers(r)%commit()
    end do
  end if

  call soca_fields_read_finalize(self, plan)
end subroutine soca_fields_read


! ------------------------------------------------------------------------------
!> Build a read_plan_t for one soca_fields member: handle pre-read config
!! (Identity, vdate from "date", h_common for vertical remap), open per-domain
!! readers, and enqueue reads onto them. The readers are NOT committed here --
!! callers do that, either per-reader (single-state path) or in bulk via
!! soca_io_readers_commit_ensemble (ensemble path).
!!
!! \param[inout] self  fields to populate
!! \param[in]    f_conf config block (basename, *_filename, read_from_file, ...)
!! \param[inout] vdate  set from f_conf%"date" if Identity or iread==1
!! \param[out]   plan   per-domain readers + scratch + post-processing flags
!! \param[in]    reader_pe optional reader_pe for all readers built here
!!               (h_common + per-domain). Ensemble path passes
!!               soca_io_ensemble_root_pe(m, nmembers) so different members'
!!               readers run on different PEs concurrently.
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_read_prepare(self, f_conf, vdate, plan, reader_pe)
  class(soca_fields), target,    intent(inout) :: self
  type(fckit_configuration),     intent(in)    :: f_conf
  type(datetime),                intent(inout) :: vdate
  type(read_plan_t), target,     intent(inout) :: plan
  integer, optional,             intent(in)    :: reader_pe

  character(len=:), allocatable :: str, basename, filename
  type(soca_io_reader) :: hreader
  integer :: d, f, i, n, idx, nreaders, ri
  type(soca_field_metadata) :: field_meta
  character(len=3), dimension(5) :: domains

  domains = [character(len=3) :: "ocn", "sfc", "ice", "wav", "bio"]
  ! Default iread to 1 when 'basename' is set: the original soca_fields_read had
  ! `integer :: iread = 0` which (Fortran gotcha) gave it implicit SAVE, so
  ! callers like saber's "model file" path (basename + ocn_filename but no
  ! read_from_file) silently relied on iread being inherited from a prior call.
  ! Inferring iread=1 from the presence of basename keeps that behavior
  ! deterministic without forcing a YAML schema change.
  plan%iread = 0
  if (f_conf%has("basename")) plan%iread = 1
  if (f_conf%has("read_from_file")) call f_conf%get_or_die("read_from_file", plan%iread)

  ! Vertical remap target: read h_common immediately on its own reader.
  ! Kept serial-per-member in v1; lives outside plan%readers so the bulk
  ! ensemble commit isn't blocked on this small single-field read.
  if (f_conf%has("remap_filename")) then
    call f_conf%get_or_die("remap_filename", str)
    allocate(plan%h_common(self%geom%isd:self%geom%ied, &
                           self%geom%jsd:self%geom%jed, self%geom%nzo))
    plan%h_common = 0.0_kind_real
    if (present(reader_pe)) then
      call hreader%init(self%geom%Domain%mpp_domain, str, reader_pe=reader_pe)
    else
      call hreader%init(self%geom%Domain%mpp_domain, str)
    end if
    call hreader%enqueue('h', plan%h_common)
    call hreader%commit()
  end if

  ! Identity (unit increment); independent of iread. Sets vdate from "date".
  if (f_conf%has("Identity")) then
    call f_conf%get_or_die("Identity", i)
    if (i == 1) call self%ones()
    call f_conf%get_or_die("date", str)
    call datetime_set(str, vdate)
  end if

  if (plan%iread /= 1 .and. plan%iread /= 3) then
    allocate(plan%readers(0), plan%domain_vars(0))
    return
  end if

  plan%seaice_categories_vars = oops_variables()

  if (plan%iread == 1) then
    call f_conf%get_or_die("date", str)
    call datetime_set(str, vdate)
  end if
  call f_conf%get_or_die("basename", basename)

  ! Probe for compute_ice/snow thickness from the ice file (if any).
  if (f_conf%get("ice_filename", str)) then
    filename = trim(basename) // trim(str)
    plan%ice_filename = filename
    field_meta = self%geom%fields_metadata%get("sea_ice_thickness")
    if (.not. soca_io_var_exists(filename, field_meta%io_name)) then
      plan%compute_icethickness = .true.
    end if
    field_meta = self%geom%fields_metadata%get("sea_ice_snow_thickness")
    if (.not. soca_io_var_exists(filename, field_meta%io_name)) then
      plan%compute_snowthickness = .true.
    end if
  end if

  ! Pass 1: count domains that have both a filename in config and matching fields.
  nreaders = 0
  do d = 1, size(domains)
    if (.not. f_conf%get(domains(d)//"_filename", str)) cycle
    n = 0
    do i = 1, size(self%fields)
      if (self%fields(i)%metadata%io_file == domains(d)) n = n + 1
    end do
    if (n > 0) nreaders = nreaders + 1
  end do

  allocate(plan%readers(nreaders))
  allocate(plan%domain_vars(nreaders))

  ! Pass 2: init each reader, enqueue all matching vars. Vars whose data are
  ! not allocated (the deferred-CICE-category case) keep allocated(data)=.false.
  ! and finalize's copy-to-atlas loop skips them.
  ri = 0
  do d = 1, size(domains)
    if (.not. f_conf%get(domains(d)//"_filename", str)) cycle
    filename = trim(basename) // trim(str)

    n = 0
    do i = 1, size(self%fields)
      if (self%fields(i)%metadata%io_file == domains(d)) n = n + 1
    end do
    if (n == 0) cycle

    ri = ri + 1
    allocate(plan%domain_vars(ri)%vars(n))
    if (present(reader_pe)) then
      call plan%readers(ri)%init(self%geom%Domain%mpp_domain, filename, reader_pe=reader_pe)
    else
      call plan%readers(ri)%init(self%geom%Domain%mpp_domain, filename)
    end if

    n = 0
    do f = 1, size(self%fields)
      if (self%fields(f)%metadata%io_file /= domains(d)) cycle
      if (domains(d) == "ice" .and. self%fields(f)%metadata%categories > 0) then
        ! CICE-history-style aggregated cat+lev files take a separate read path
        ! in finalize via self%read_seaice. SOCA-written files have each
        ! category as a separate var and read through here normally.
        if (soca_io_file_exists(filename) .and. &
            soca_io_var_exists(filename, self%fields(f)%metadata%io_name)) then
        else
          call plan%seaice_categories_vars%push_back(self%fields(f)%name)
          cycle
        end if
      end if
      n = n + 1

      associate (v => plan%domain_vars(ri)%vars(n))
        v%field => self%fields(f)
        v%afield = self%afieldset%field(v%field%name)
        call v%afield%data(v%adata)
        allocate(v%data(self%geom%isd:self%geom%ied, &
                        self%geom%jsd:self%geom%jed, v%field%nz))
        if (v%field%nz == 1) then
          if ((self%fields(f)%metadata%name == "sea_ice_thickness") &
              .and. plan%compute_icethickness) then
            field_meta = self%geom%fields_metadata%get("sea_ice_volume")
            call plan%readers(ri)%enqueue(field_meta%io_name, v%data(:,:,1))
          else if ((self%fields(f)%metadata%name == "sea_ice_snow_thickness") &
                   .and. plan%compute_snowthickness) then
            field_meta = self%geom%fields_metadata%get("sea_ice_snow_volume")
            call plan%readers(ri)%enqueue(field_meta%io_name, v%data(:,:,1))
          else
            call plan%readers(ri)%enqueue(v%field%metadata%io_name, v%data(:,:,1))
          end if
        else
          call plan%readers(ri)%enqueue(v%field%metadata%io_name, v%data(:,:,:))
        end if
      end associate
    end do
    ! unused slots at tail of plan%domain_vars(ri)%vars (skipped CICE cases)
    ! keep allocated(v%data)=.false.; finalize handles them.
  end do
end subroutine soca_fields_read_prepare


! ------------------------------------------------------------------------------
!> Drain a populated read_plan_t into self: copy per-domain scratch into atlas
!! fields with land-mask handling, run the existing post-processing chain
!! (CICE-category deferred read, ice/snow thickness compute, mid-layer depth,
!! mixed-layer depth, optional vertical remap), and release scratch.
!!
!! Assumes readers in plan%readers have already been committed.
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_read_finalize(self, plan)
  class(soca_fields), target,    intent(inout) :: self
  type(read_plan_t), target,     intent(inout) :: plan

  integer :: d, f, i, j, k, n, idx, ri
  type(remapping_CS) :: remapCS
  type(atlas_field) :: afield1, afield2, afield3, afield4
  real(kind=kind_real), pointer :: adata1(:,:), adata2(:,:), adata3(:,:), adata4(:,:)
  real(kind=kind_real), allocatable :: h_common_ij(:), hocn_ij(:), varocn_ij(:), varocn2_ij(:)
  logical :: remap_this_point

  ! Identity-only / no-read path still needs h_common dealloc; nothing else to do.
  if (plan%iread /= 1 .and. plan%iread /= 3) then
    if (allocated(plan%h_common)) deallocate(plan%h_common)
    return
  end if

  ! Set constant-valued atlas fields (independent of any read).
  do f = 1, size(self%fields)
    if (self%fields(f)%metadata%io_file == "CONSTANT") then
      afield1 = self%afieldset%field(self%fields(f)%name)
      call afield1%data(adata1)
      adata1(:,:) = self%fields(f)%metadata%constant_value
    end if
  end do

  ! Copy per-domain scratch buffers into atlas fields with land-mask handling.
  do ri = 1, size(plan%readers)
    associate (vars => plan%domain_vars(ri)%vars)
      do n = 1, size(vars)
        if (.not. allocated(vars(n)%data)) cycle  ! deferred CICE-category slot
        vars(n)%adata(:,:) = 0.0
        do j = self%geom%jsc, self%geom%jec
          do i = self%geom%isc, self%geom%iec
            idx = self%geom%atlas_ij2idx(i,j)
            if (associated(vars(n)%field%mask)) then
              if (vars(n)%field%mask(i,j) == 0) then
                vars(n)%adata(:, idx) = vars(n)%field%metadata%fillvalue
              else
                vars(n)%adata(:, idx) = vars(n)%data(i,j,:)
              end if
            else
              vars(n)%adata(:, idx) = vars(n)%data(i,j,:)
            end if
          end do
        end do
        call vars(n)%afield%set_dirty()
      end do
      do n = 1, size(vars)
        if (.not. allocated(vars(n)%data)) cycle
        deallocate(vars(n)%data)
        call vars(n)%afield%final()
      end do
    end associate
    deallocate(plan%domain_vars(ri)%vars)
  end do

  ! Deferred-CICE-category read (its own internal I/O cycle; kept serial in v1).
  if (plan%seaice_categories_vars%nvars() > 0) then
    call self%read_seaice(plan%ice_filename, plan%seaice_categories_vars)
  end if

  ! Compute ice thickness from ice_volume / ice_area_fraction.
  if (plan%compute_icethickness .and. self%afieldset%has("sea_ice_thickness")) then
    afield1 = self%afieldset%field("sea_ice_thickness")
    afield2 = self%afieldset%field("sea_ice_area_fraction")
    call afield1%data(adata1)
    call afield2%data(adata2)
    do j = self%geom%jsc, self%geom%jec
      do i = self%geom%isc, self%geom%iec
        idx = self%geom%atlas_ij2idx(i,j)
        if (adata2(1,idx) > 0.0) then
          adata1(1,idx) = adata1(1,idx) / adata2(1,idx)
        else
          adata1(1,idx) = 0.0_kind_real
        end if
      end do
    end do
    call afield1%set_dirty()
  end if

  ! Compute snow thickness from snow_volume / ice_area_fraction.
  if (plan%compute_snowthickness .and. self%afieldset%has("sea_ice_snow_thickness")) then
    afield1 = self%afieldset%field("sea_ice_snow_thickness")
    afield2 = self%afieldset%field("sea_ice_area_fraction")
    call afield1%data(adata1)
    call afield2%data(adata2)
    do j = self%geom%jsc, self%geom%jec
      do i = self%geom%isc, self%geom%iec
        idx = self%geom%atlas_ij2idx(i,j)
        if (adata2(1,idx) > 0.0) then
          adata1(1,idx) = adata1(1,idx) / adata2(1,idx)
        else
          adata1(1,idx) = 0.0_kind_real
        end if
      end do
    end do
    call afield1%set_dirty()
  end if

  ! TODO: this should live in a variable-change class, not here.
  ! Mid-layer depth from layer thickness (running half-thickness sum).
  if (self%afieldset%has("sea_water_depth")) then
    afield1 = self%afieldset%field("sea_water_depth")
    afield2 = self%afieldset%field("sea_water_cell_thickness")
    call afield1%data(adata1)
    call afield2%data(adata2)
    do j = self%geom%jsc, self%geom%jec
      do i = self%geom%isc, self%geom%iec
        idx = self%geom%atlas_ij2idx(i,j)
        adata1(:,idx) = 0.5 * adata2(:,idx)
        do k = 2, afield1%shape(1)
          adata1(k,idx) = adata1(k,idx) + sum(adata2(1:k-1,idx))
        end do
      end do
    end do
    call afield1%set_dirty()
  end if

  ! TODO: move out. Mixed-layer depth from T/S/depth.
  if (self%afieldset%has("ocean_mixed_layer_thickness")) then
    afield1 = self%afieldset%field("ocean_mixed_layer_thickness")
    afield2 = self%afieldset%field("sea_water_salinity")
    afield3 = self%afieldset%field("sea_water_potential_temperature")
    afield4 = self%afieldset%field("sea_water_depth")
    call afield1%data(adata1)
    call afield2%data(adata2)
    call afield3%data(adata3)
    call afield4%data(adata4)
    adata1(:,:) = 0.0_kind_real
    do j = self%geom%jsc, self%geom%jec
      do i = self%geom%isc, self%geom%iec
        if (self%geom%mask2d(i,j) == 0) cycle
        idx = self%geom%atlas_ij2idx(i,j)
        adata1(1,idx) = soca_mld(adata2(:,idx), adata3(:,idx), adata4(:,idx), &
          self%geom%lon(i,j), self%geom%lat(i,j))
      end do
    end do
    call afield1%set_dirty()
  end if

  ! Optional vertical remap onto h_common.
  if (allocated(plan%h_common)) then
    call initialize_remapping(remapCS, 'PCM')
    allocate(h_common_ij(self%geom%nzo), hocn_ij(self%geom%nzo), &
             varocn_ij(self%geom%nzo), varocn2_ij(self%geom%nzo))
    afield1 = self%afieldset%field("sea_water_cell_thickness")
    call afield1%data(adata1)
    do n = 1, size(self%fields)
      if (.not. self%fields(n)%metadata%vert_interp) cycle
      if (self%geom%f_comm%rank() == 0) then
        call oops_log%info("vertically remapping "//trim(self%fields(n)%name))
      end if
      afield2 = self%afieldset%field(self%fields(n)%name)
      call afield2%data(adata2)
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          remap_this_point = .not. associated(self%fields(n)%mask)
          if (associated(self%fields(n)%mask)) then
            remap_this_point = self%fields(n)%mask(i,j) .gt. 0.0
          end if
          if (remap_this_point) then
            h_common_ij(:) = plan%h_common(i,j,:)
            hocn_ij(:) = adata1(:, idx)
            varocn_ij(:) = adata2(:, idx)
            call remapping_core_h(remapCS, self%geom%nzo, h_common_ij, varocn_ij, &
                                  self%geom%nzo, hocn_ij, varocn2_ij)
            adata2(:, idx) = varocn2_ij
          else
            adata2(:, idx) = 0.0_kind_real
          end if
        end do
      end do
      call afield2%set_dirty()
    end do
    call end_remapping(remapCS)
    deallocate(h_common_ij, hocn_ij, varocn_ij, varocn2_ij)
  end if

  ! Release the empty-ctor'd seaice_categories_vars handle from prepare so
  ! the C++-side variables_ allocation doesn't leak per-member.
  call plan%seaice_categories_vars%destruct()

  if (allocated(plan%h_common)) deallocate(plan%h_common)
  call afield1%final()
  call afield2%final()
  call afield3%final()
  call afield4%final()
end subroutine soca_fields_read_finalize


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
  character(len=*), intent(in) :: filename
  type(oops_variables), intent(in) :: seaice_categories_vars

  type(oops_variables) :: cice_vars_cats, cice_vars_cats_levs
  type(soca_io_reader) :: reader
  type(atlas_field) :: afield
  real(kind=kind_real), pointer :: adata(:,:)

  integer :: i, j, k, f, ncat, icelevs, snowlevs, cnt, io_index, idx
  ! `target` so soca_io_reader can pointer-associate to slices of these
  real(kind=kind_real), allocatable, target :: tmp3d(:,:,:,:), tmp4d(:,:,:,:,:)

  ! check what cice variables with category dimension need to be read
  cice_vars_cats = oops_variables()  ! used to store the unique cice io variables with a category dimension
  call get_cice_vars(self, cice_vars_cats, ncat, icelevs, "dynam")

  ! read the cice variables with category dimension only
  if (cice_vars_cats%nvars() > 0) then
    allocate(tmp3d(self%geom%isd:self%geom%ied,self%geom%jsd:self%geom%jed,ncat,cice_vars_cats%nvars()))
    tmp3d = 0.0_kind_real
    call reader%init(self%geom%Domain%mpp_domain, filename)
    do i=1,cice_vars_cats%nvars()
      call reader%enqueue(cice_vars_cats%variable(i), tmp3d(:,:,:,i))
    end do
    call reader%commit()

    ! copy the variable into the corresponding field
    cnt = 1
    do f = 1, size(self%fields)
      if (self%fields(f)%metadata%io_file == "ice" .and.&
         &self%fields(f)%metadata%levels == '1' .and.&
         &self%fields(f)%metadata%categories > 0) then

        ! get the index of cice_vars that correspond to the io_sup_name
        io_index = cice_vars_cats%find(self%fields(f)%metadata%io_sup_name)

        afield = self%afieldset%field(self%fields(f)%name)
        call afield%data(adata)
        do j=self%geom%jsc, self%geom%jec
          do i=self%geom%isc, self%geom%iec
            idx = self%geom%atlas_ij2idx(i,j)
            adata(1,idx) = tmp3d(i,j,self%fields(f)%metadata%category,io_index)
          end do
        end do
        call afield%set_dirty()
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
    call reader%init(self%geom%Domain%mpp_domain, filename)
    do i=1,cice_vars_cats_levs%nvars()
      call reader%enqueue(cice_vars_cats_levs%variable(i), tmp4d(:,:,:,:,i))
    end do
    call reader%commit()

    ! copy the variable into the corresponding field
    cnt = 1
    do f = 1, size(self%fields)
      if (self%fields(f)%metadata%io_file == "ice" .and.&
         &self%fields(f)%nz > 1 .and.&
         &self%fields(f)%metadata%categories > 0) then

        ! get the index of cice_vars that correspond to the io_sup_name
        io_index = cice_vars_cats_levs%find(self%fields(f)%metadata%io_sup_name)

        afield = self%afieldset%field(self%fields(f)%name)
        call afield%data(adata)
        do j=self%geom%jsc, self%geom%jec
          do i=self%geom%isc, self%geom%iec
            idx = self%geom%atlas_ij2idx(i,j)
            adata(:,idx) = tmp4d(i,j,:,self%fields(f)%metadata%category,io_index)
          end do
        end do
        call afield%set_dirty()
      end if
    end do
  end if
  call afield%final()
end subroutine soca_fields_read_seaice


! ------------------------------------------------------------------------------
!> Build one per-domain writer for `fld`: enqueue every field whose io_file
!! matches `dom`, copying each atlas field's compute-domain slice into a fresh
!! varwrapper buffer (masked cells set to the field's fillvalue). `vars` carries
!! the buffers and has the TARGET attribute because the writer holds pointers
!! into them; the caller must keep `vars` alive (and unmutated) until the
!! writer is committed. Leaves `vars` unallocated and `writer` untouched when no
!! field maps to `dom`. Shared by soca_fields_write_rst (single state) and
!! soca_fields_write_ensemble (bulk), so the atlas->compute-domain copy and the
!! enqueue logic live in exactly one place.
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_enqueue_domain(fld, dom, filename, writer, vars, root_pe)
  class(soca_fields),            intent(inout) :: fld
  character(len=*),              intent(in)    :: dom
  character(len=*),              intent(in)    :: filename
  type(soca_io_writer),          intent(inout) :: writer
  type(varwrapper), allocatable, target, intent(out) :: vars(:)
  integer, optional,             intent(in)    :: root_pe

  integer :: i, j, f, n, idx

  ! How many fields land in this domain's file?
  n = 0
  do f = 1, size(fld%fields)
    if (fld%fields(f)%metadata%io_file == dom) n = n + 1
  end do
  if (n == 0) return

  if (present(root_pe)) then
    call writer%init(fld%geom%Domain%mpp_domain, filename, root_pe=root_pe)
  else
    call writer%init(fld%geom%Domain%mpp_domain, filename)
  end if

  allocate(vars(n))
  n = 0
  do f = 1, size(fld%fields)
    if (fld%fields(f)%metadata%io_file /= dom) cycle
    n = n + 1

    ! copy atlas field into a per-PE compute-domain temporary
    vars(n)%afield = fld%aFieldset%field(fld%fields(f)%name)
    call vars(n)%afield%data(vars(n)%adata)
    allocate(vars(n)%data(fld%geom%isd:fld%geom%ied, &
                          fld%geom%jsd:fld%geom%jed, vars(n)%afield%shape(1)))

    ! copy, setting masked values to fillvalue
    if (associated(fld%fields(f)%mask)) vars(n)%data = fld%fields(f)%metadata%fillvalue
    do j = fld%geom%jsc, fld%geom%jec
      do i = fld%geom%isc, fld%geom%iec
        idx = fld%geom%atlas_ij2idx(i, j)
        if (associated(fld%fields(f)%mask)) then
          if (fld%fields(f)%mask(i, j) == 0) cycle
        end if
        vars(n)%data(i, j, :) = vars(n)%adata(:, idx)
      end do
    end do

    ! enqueue with the writer
    if (vars(n)%afield%shape(1) == 1) then
      call writer%enqueue(fld%fields(f)%metadata%io_name, vars(n)%data(:,:,1))
    else
      call writer%enqueue(fld%fields(f)%metadata%io_name, vars(n)%data(:,:,:))
    end if
  end do
end subroutine soca_fields_enqueue_domain


! ------------------------------------------------------------------------------
!> Save soca fields in a restart format
!!
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_write_rst(self, f_conf, vdate)
  class(soca_fields), target, intent(inout) :: self      !< Fields
  type(fckit_configuration),  intent(in)    :: f_conf   !< Configuration
  type(datetime),             intent(inout) :: vdate    !< DateTime

  integer, parameter :: max_string_length=800
  type(soca_io_writer) :: writer
  integer :: d, n
  logical :: date_cols

  character(len=3), allocatable :: domains(:)
  character(len=:), allocatable :: domain_filename

  ! target so vars(n)%data sections are valid pointer targets for soca_io enqueue
  type(varwrapper), allocatable, target :: vars(:)

  ! Get date IO format (colons or not?)
  date_cols = .true.
  if (f_conf%has("date colons")) then
    call f_conf%get_or_die("date colons", date_cols)
  end if

  ! Set up domain info
  domains = [character(len=3) :: "ocn", "sfc", "ice", "wav", "bio"]

  ! for each domain, get the fields to be written out and write them
  do d=1,size(domains)
    domain_filename = soca_genfilename(f_conf,max_string_length,vdate,date_cols,domains(d))

    call soca_fields_enqueue_domain(self, domains(d), domain_filename, writer, vars)
    if (.not. allocated(vars)) cycle  ! no fields in this domain

    ! write the file
    call writer%commit()

    ! cleanup
    do n=1,size(vars)
      deallocate(vars(n)%data)
      call vars(n)%afield%final()
    end do
    deallocate(vars)
  end do
end subroutine soca_fields_write_rst

! ------------------------------------------------------------------------------
!> Bulk write multiple states/increments in one ensemble pass. Each member's
!! per-domain files are written from a rotated writer PE (assigned by
!! soca_io_ensemble_root_pe) so disk writes happen concurrently across the
!! world communicator. Falls back to sequential per-member commit when
!! geometry.io.'ensemble write' is 'sequential'.
!!
!! Behavior per-member is identical to soca_fields_write_rst: both build each
!! per-domain writer via soca_fields_enqueue_domain (same domain set, filename,
!! mask/fillvalue handling and byte-level netCDF output). Only the dispatch
!! (rotated root + ensemble orchestrator) changes.
!!
!! \param fields_ptrs polymorphic pointers to soca_state or soca_increment
!! \param f_confs    per-member configurations (must align 1:1 with fields_ptrs)
!! \param vdates     per-member valid dates
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_write_ensemble(fields_ptrs, c_confs, vdates)
  type(soca_fields_ptr_t),    intent(inout) :: fields_ptrs(:)
  type(c_ptr),                intent(in)    :: c_confs(:)
  type(datetime),             intent(inout) :: vdates(:)

  integer, parameter :: max_string_length = 800
  integer :: nmembers, m, d, f, n, wi, nwriters, root_pe
  character(len=3), allocatable :: domains(:)
  character(len=:), allocatable :: domain_filename
  logical :: date_cols
  class(soca_fields), pointer :: fld
  type(soca_io_writer), allocatable, target :: writers(:)
  type(ensemble_vars_t), allocatable, target :: vars_arr(:)
  type(fckit_configuration) :: f_conf

  nmembers = size(fields_ptrs)
  if (nmembers == 0) return
  if (size(c_confs) /= nmembers .or. size(vdates) /= nmembers) then
    call abor1_ftn("soca_fields_write_ensemble: input array size mismatch")
  end if

  domains = [character(len=3) :: "ocn", "sfc", "ice", "wav", "bio"]

  ! Pass 1: count writers (one per (member, domain) with vars).
  nwriters = 0
  do m = 1, nmembers
    fld => fields_ptrs(m)%p
    do d = 1, size(domains)
      n = 0
      do f = 1, size(fld%fields)
        if (fld%fields(f)%metadata%io_file == domains(d)) n = n + 1
      end do
      if (n > 0) nwriters = nwriters + 1
    end do
  end do
  if (nwriters == 0) return

  allocate(writers(nwriters))
  allocate(vars_arr(nwriters))

  ! Pass 2: build writers and enqueue. Construct fckit_configuration on the fly
  ! per member (array assignment of fckit_configuration drops the underlying
  ! c_ptr in some compiler/version combos, so we materialize a fresh handle
  ! each iteration -- same pattern as the single-state interface).
  wi = 0
  do m = 1, nmembers
    fld => fields_ptrs(m)%p
    f_conf = fckit_configuration(c_confs(m))
    date_cols = .true.
    if (f_conf%has("date colons")) then
      call f_conf%get_or_die("date colons", date_cols)
    end if

    root_pe = soca_io_ensemble_root_pe(m, nmembers)

    do d = 1, size(domains)
      ! count vars for this (member, domain); must stay in lockstep with the
      ! pass-1 nwriters count so wi indexes the preallocated writers array.
      n = 0
      do f = 1, size(fld%fields)
        if (fld%fields(f)%metadata%io_file == domains(d)) n = n + 1
      end do
      if (n == 0) cycle

      wi = wi + 1
      domain_filename = soca_genfilename(f_conf, max_string_length, vdates(m), date_cols, domains(d))
      call soca_fields_enqueue_domain(fld, domains(d), domain_filename, &
                                      writers(wi), vars_arr(wi)%vars, root_pe=root_pe)
    end do
  end do

  ! Commit. The 'sequential' fallback is kept for A/B; it does not benefit
  ! from rotated root_pes (each writer still gathers to its assigned PE), so
  ! reproducibility holds.
  if (soca_io_ensemble_write_parallel()) then
    call soca_io_writers_commit_ensemble(writers)
  else
    do wi = 1, nwriters
      call writers(wi)%commit()
    end do
  end if

  ! Cleanup. Atlas fields use refcounted handles -> need explicit final().
  do wi = 1, nwriters
    if (.not. allocated(vars_arr(wi)%vars)) cycle
    do n = 1, size(vars_arr(wi)%vars)
      if (allocated(vars_arr(wi)%vars(n)%data)) deallocate(vars_arr(wi)%vars(n)%data)
      call vars_arr(wi)%vars(n)%afield%final()
    end do
    deallocate(vars_arr(wi)%vars)
  end do
  deallocate(writers, vars_arr)
end subroutine soca_fields_write_ensemble


! ------------------------------------------------------------------------------
!> Bulk read multiple states/increments. Mirrors soca_fields_write_ensemble:
!! prepare per-member read_plan_t with a rotated reader_pe, flatten all
!! per-(member, domain) readers into a single array, drive the bulk I/O via
!! soca_io_readers_commit_ensemble (strided or scatter, sync or async per the
!! geometry.io knobs), then finalize each member to copy scratch into atlas
!! fields and run post-read transforms.
!!
!! \param fields_ptrs polymorphic pointers to soca_state or soca_increment
!! \param c_confs     per-member configurations (c_ptrs; materialized inline)
!! \param vdates      per-member valid dates
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_read_ensemble(fields_ptrs, c_confs, vdates)
  type(soca_fields_ptr_t),    intent(inout) :: fields_ptrs(:)
  type(c_ptr),                intent(in)    :: c_confs(:)
  type(datetime),             intent(inout) :: vdates(:)

  integer :: nmembers, m, ri, nreaders_total, root_pe, k
  type(read_plan_t), allocatable, target :: plans(:)
  type(soca_io_reader), allocatable :: flat_readers(:)
  type(fckit_configuration) :: f_conf

  nmembers = size(fields_ptrs)
  if (nmembers == 0) return
  if (size(c_confs) /= nmembers .or. size(vdates) /= nmembers) then
    call abor1_ftn("soca_fields_read_ensemble: input array size mismatch")
  end if

  allocate(plans(nmembers))

  ! Pass 1: prepare each member. Rotated reader_pe spreads the bulk reads
  ! across PEs so they hit different files concurrently in the scatter path.
  ! Construct fckit_configuration on the fly per member -- array-of-handles
  ! through the c-binding has dropped the c_ptr in some compiler/version
  ! combos; same pattern as soca_fields_write_ensemble.
  do m = 1, nmembers
    f_conf = fckit_configuration(c_confs(m))
    root_pe = soca_io_ensemble_root_pe(m, nmembers)
    call soca_fields_read_prepare(fields_ptrs(m)%p, f_conf, vdates(m), plans(m), &
                                  reader_pe=root_pe)
  end do

  ! Pass 2: flatten all per-domain readers into one array for the bulk commit.
  ! Intrinsic assignment deep-copies the read_entry vars(:) but shallow-copies
  ! their pointer components (data1d/2d/3d/4d -> plans(m)%domain_vars(...)%data
  ! scratch buffers) and the domain pointer. The read_entry gbuf_* allocatables
  ! are unallocated at copy time (only populated during reader_stage_read on
  ! flat_readers), so the deep-copy doesn't duplicate any payload. The bulk
  ! scatter delivers data straight into the buffers finalize reads from.
  nreaders_total = 0
  do m = 1, nmembers
    nreaders_total = nreaders_total + size(plans(m)%readers)
  end do
  if (nreaders_total > 0) then
    allocate(flat_readers(nreaders_total))
    k = 0
    do m = 1, nmembers
      do ri = 1, size(plans(m)%readers)
        k = k + 1
        flat_readers(k) = plans(m)%readers(ri)
      end do
    end do
    call soca_io_readers_commit_ensemble(flat_readers)
    deallocate(flat_readers)
  end if

  ! Pass 3: per-member finalize -- copy scratch -> atlas, run post-processing
  ! (CICE category deferred read, ice/snow thickness, mid-layer depth, MLD,
  ! vertical remap), release scratch.
  do m = 1, nmembers
    call soca_fields_read_finalize(fields_ptrs(m)%p, plans(m))
  end do

  deallocate(plans)
end subroutine soca_fields_read_ensemble


! ------------------------------------------------------------------------------
!> Interpolates from uv-points location to h-points.
!!
!! \relates soca_fields_mod::soca_fields
subroutine soca_fields_tohpoints(self)
  class(soca_fields), intent(inout) :: self !< self

  integer :: i,j,k,n,idx
  character(len=4) :: fromto

  type(atlas_field) :: afield
  real(kind=kind_real), pointer :: adata(:,:)
  real(kind=kind_real), allocatable :: fdata(:,:,:)

  ! Apply interpolation to all fields, when necessary
  do n=1,size(self%fields)
    ! Check if already on h-points
    if (self%fields(n)%metadata%grid == 'h') cycle

    ! Interpolate to different location of the stencil
    fromto = self%fields(n)%metadata%grid//'toh'

    ! convert from atlas to 3d fortran field...
    ! because I don't want to fully refactor stencil interpolation
    allocate(fdata(self%geom%isd:self%geom%ied, self%geom%jsd:self%geom%jed, self%fields(n)%nz))
    afield = self%aFieldset%field(self%fields(n)%name)
    call afield%data(adata)
    do j=self%geom%jsc, self%geom%jec
      do i=self%geom%isc, self%geom%iec
        idx = self%geom%atlas_ij2idx(i,j)
        do k=1,afield%shape(1)
          fdata(i,j,k) = adata(k, idx)
        end do
      end do
    end do
    call mpp_update_domains(fdata, self%geom%Domain%mpp_domain)

    ! interp
    call soca_field_stencil_interp(fdata, self%geom, fromto)

    !copy back to atlas
    do j=self%geom%jsc, self%geom%jec
      do i=self%geom%isc, self%geom%iec
        idx = self%geom%atlas_ij2idx(i,j)
        do k=1,afield%shape(1)
          adata(k, idx) = fdata(i,j,k)
        end do
      end do
    end do
    deallocate(fdata)

    call afield%set_dirty()

    ! Update grid location to h-points
    self%fields(n)%metadata%grid = 'h'
    self%fields(n)%lon => self%geom%lon
    self%fields(n)%lat => self%geom%lat
 end do
 call afield%final()

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
  type(atlas_field) :: afield
  real(kind=kind_real), pointer :: adata(:,:)
  integer :: f, i, j, idx
  character(len=:), allocatable :: vars_str(:)

  type(atlas_metadata) :: ameta
  type(soca_field_metadata) :: metadata

  ! reinitialize variable parameters
  if (allocated(self%fields)) deallocate(self%fields)
  allocate(character(len=100) :: vars_str(vars%nvars()))
  do i=1,vars%nvars()
    vars_str(i) = vars%variable(i)
  end do
  call soca_fields_init_vars(self, vars_str)

  ! create new atlas fields
  do f=1,size(self%fields)
    if (.not. self%aFieldset%has(self%fields(f)%name)) then
      afield = self%geom%functionspace%create_field( &
        name=self%fields(f)%name, kind=atlas_real(kind_real), &
        levels=self%fields(f)%nz)
      call afield%data(adata)
      call self%afieldset%add(afield)
      adata(:,:) = 0.0_kind_real

      ! set metadata
      ameta = afield%metadata()
      metadata = self%geom%fields_metadata%get(afield%name())
      call ameta%set('interp_type', 'default')
      if (metadata%masked) then
        call ameta%set('mask', 'interp_mask')
      end if
      call ameta%set('nearest 3d level', 'top')
    end if
  end do
  call ameta%final()
  call afield%final()
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
     lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep) + 3
     soca_genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep) // ".nc"
  endif

  if (typ=="an" .or. typ=="incr") then
     if (date_cols) then
       call datetime_to_string(vdate,validitydate)
     else
       call datetime_to_string_io(vdate,validitydate)
     endif
     lenfn = lenfn + 1 + LEN_TRIM(validitydate) + 3
     soca_genfilename = TRIM(prefix) // "." // TRIM(validitydate) // ".nc"
  endif

  if (lenfn>length) &
       & call abor1_ftn("fields:genfilename: filename too long")

   if ( allocated(str) ) deallocate(str)

end function soca_genfilename

end module soca_fields_mod
