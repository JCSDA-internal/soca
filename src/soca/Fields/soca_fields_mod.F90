! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Handle fields for the model

module soca_fields_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: log, fckit_log
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_min, fckit_mpi_max, &
                            fckit_mpi_sum
use unstructured_grid_mod, only: unstructured_grid, &
                                 allocate_unstructured_grid_coord, &
                                 allocate_unstructured_grid_field
use datetime_mod, only: datetime, datetime_set
use duration_mod, only: duration
use kinds, only: kind_real
use oops_variables_mod
use fms_mod,    only: read_data, write_data, set_domain
use fms_io_mod, only: fms_io_init, fms_io_exit, &
                      register_restart_field, restart_file_type, &
                      restore_state, query_initialized, &
                      free_restart_type, save_restart
use mpp_domains_mod, only : mpp_update_domains
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h, end_remapping
use random_mod, only: normal_distribution
use soca_geom_mod, only : soca_geom
use soca_geom_iter_mod, only : soca_geom_iter
use soca_fieldsutils_mod, only: soca_genfilename, fldinfo
use soca_utils, only: soca_mld
use soca_geostrophy_mod

implicit none

private
public :: soca_fields, soca_field

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

!> Holds all data and metadata related to a single field variable
type :: soca_field
  character(len=:),     allocatable :: name       !< the internally used name of the field
  integer                           :: nz         !< the number of levels
  real(kind=kind_real), allocatable :: val(:,:,:) !< the actual data
  real(kind=kind_real),     pointer :: mask(:,:) => null() !< field mask
  character(len=:),     allocatable :: cf_name    !< the (optional) name needed by UFO
  character(len=:),     allocatable :: io_name    !< the (optional) name use in the restart IO
  character(len=:),     allocatable :: io_file    !< the (optional) restart file domain
                                                  ! (ocn, sfc, ice)
contains
  procedure :: delete  => soca_field_delete
  procedure :: copy    => soca_field_copy
  procedure :: check_congruent => soca_field_check_congruent
end type soca_field

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

!> Holds a collection of soca_field types, and the public suroutines
!> to manipulate them. Represents all the fields of a given state
!> or increment
type :: soca_fields
   type(soca_geom),  pointer :: geom           !< MOM6 Geometry
   type(soca_field), pointer :: fields(:) => null()

   ! Ocean diagnostics (TODO, generalize)
   real(kind=kind_real), allocatable :: mld(:,:)           !< Mixed layer depth (nx,ny)
   real(kind=kind_real), allocatable :: layer_depth(:,:,:) !< Mid-layer depth (nx,ny,nz0)

contains
  ! constructors / destructors
  procedure :: copy   => soca_fields_copy
  procedure :: create => soca_fields_create
  procedure :: delete => soca_fields_delete

  ! field getters
  procedure :: get    => soca_fields_get
  procedure :: has    => soca_fields_has
  procedure :: check_congruent => soca_fields_check_congruent

  ! math
  procedure :: add      => soca_fields_add
  procedure :: add_incr => soca_fields_add_incr
  procedure :: axpy     => soca_fields_axpy
  procedure :: diff_incr=> soca_fields_diff_incr
  procedure :: dirac    => soca_fields_dirac
  procedure :: dot_prod => soca_fields_dotprod
  procedure :: gpnorm   => soca_fields_gpnorm
  procedure :: mul      => soca_fields_mul
  procedure :: random   => soca_fields_random
  procedure :: schur    => soca_fields_schur
  procedure :: sub      => soca_fields_sub
  procedure :: zeros    => soca_fields_zeros

  ! IO
  procedure :: from_ug   => soca_fields_from_ug
  procedure :: to_ug     => soca_fields_to_ug
  procedure :: ug_coord  => soca_fields_ug_coord

  procedure :: getpoint  => soca_fields_getpoint
  procedure :: setpoint  => soca_fields_setpoint

  procedure :: read      => soca_fields_read
  procedure :: write_file=> soca_fields_write_file
  procedure :: write_rst => soca_fields_write_rst
end type soca_fields


contains


! ------------------------------------------------------------------------------
! soca_field subroutines
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> copy a field from rhs to self. Self must be allocated first
subroutine soca_field_copy(self, rhs)
  class(soca_field), intent(inout) :: self
  type(soca_field),  intent(in)    :: rhs

  call self%check_congruent(rhs)

  ! copy fields
  self%name = rhs%name
  self%nz = rhs%nz
  self%val = rhs%val
  self%cf_name = rhs%cf_name
  self%io_name = rhs%io_name
  self%io_file = rhs%io_file
  self%mask => rhs%mask
end subroutine soca_field_copy


! ------------------------------------------------------------------------------
! make sure the two fields are the same in terms of name, size, shape
subroutine soca_field_check_congruent(self, rhs)
  class(soca_field), intent(in) :: self
  type(soca_field),  intent(in) :: rhs
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
!> delete the soca_field object
subroutine soca_field_delete(self)
  class(soca_field), intent(inout) :: self

  deallocate(self%val)
end subroutine


! ------------------------------------------------------------------------------
! soca_fields subroutines
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> for a given list of field names, initialize the properties of those fields
! NOTE: this information should be moved into a yaml file
! TODO, allocate space for derived variables
subroutine soca_fields_init_vars(self, vars)
  type(soca_fields),          intent(inout) :: self
  character(len=:), allocatable, intent(in) :: vars(:)

  integer :: i, nz

  allocate(self%fields(size(vars)))
  do i=1,size(vars)
    self%fields(i)%name = trim(vars(i))

    ! determine number of levels, and if masked
    select case(self%fields(i)%name)
    case ('tocn','socn','hocn')
      nz = self%geom%nzo
      self%fields(i)%mask => self%geom%mask2d
    case ('uocn')
      nz = self%geom%nzo
      self%fields(i)%mask => self%geom%mask2du
    case ('vocn')
      nz = self%geom%nzo
      self%fields(i)%mask => self%geom%mask2dv
    case ('hicen','hsnon', 'cicen')
      nz = self%geom%ncat
      self%fields(i)%mask => self%geom%mask2d
    case ('ssh')
      nz = 1
      self%fields(i)%mask => self%geom%mask2d
    case ('sw', 'lhf', 'shf', 'lw', 'us')
      nz = 1
    case default
      call abor1_ftn('soca_fields::create(): unknown field '// self%fields(i)%name)
    end select

    ! allocate space
    self%fields(i)%nz = nz
    allocate(self%fields(i)%val(&
      self%geom%isd:self%geom%ied, &
      self%geom%jsd:self%geom%jed, &
      nz ))

    ! set other variables associated with each field
    self%fields(i)%cf_name = ""
    self%fields(i)%io_name = ""
    select case(self%fields(i)%name)
    case ('tocn')
      self%fields(i)%cf_name = "sea_water_potential_temperature"
      self%fields(i)%io_file = "ocn"
      self%fields(i)%io_name = "Temp"
    case ('socn')
      self%fields(i)%cf_name = "sea_water_practical_salinity"
      self%fields(i)%io_file = "ocn"
      self%fields(i)%io_name = "Salt"
    case ('ssh')
      self%fields(i)%cf_name = "sea_surface_height_above_geoid"
      self%fields(i)%io_file = "ocn"
      self%fields(i)%io_name = "ave_ssh"
    case ('hocn')
      self%fields(i)%cf_name = "sea_water_cell_thickness"
      self%fields(i)%io_file = "ocn"
      self%fields(i)%io_name = "h"
   case ('uocn')
      self%fields(i)%cf_name = "sea_water_zonal_current"
      self%fields(i)%io_file = "ocn"
      self%fields(i)%io_name = "u"
   case ('vocn')
      self%fields(i)%cf_name = "sea_water_meridional_current"
      self%fields(i)%io_file = "ocn"
      self%fields(i)%io_name = "v"
    case ('hicen')
      self%fields(i)%cf_name = "sea_ice_category_thickness"
      self%fields(i)%io_file = "ice"
    case ('cicen')
      self%fields(i)%cf_name = "sea_ice_category_area_fraction"
      self%fields(i)%io_file = "ice"
    case ('hsnon')
      self%fields(i)%io_file = "ice"
    case ('sw')
      self%fields(i)%cf_name = "net_downwelling_shortwave_radiation"
      self%fields(i)%io_file = "sfc"
      self%fields(i)%io_name = "sw_rad"
    case ('lw')
      self%fields(i)%cf_name = "net_downwelling_longwave_radiation"
      self%fields(i)%io_file = "sfc"
      self%fields(i)%io_name = "lw_rad"
    case ('lhf')
      self%fields(i)%cf_name = "upward_latent_heat_flux_in_air"
      self%fields(i)%io_file = "sfc"
      self%fields(i)%io_name = "latent_heat"
    case ('shf')
      self%fields(i)%cf_name = "upward_sensible_heat_flux_in_air"
      self%fields(i)%io_file = "sfc"
      self%fields(i)%io_name = "sens_heat"
    case ('us')
      self%fields(i)%cf_name = "friction_velocity_over_water"
      self%fields(i)%io_file = "sfc"
      self%fields(i)%io_name = "fric_vel"
    end select

  end do
end subroutine


! ------------------------------------------------------------------------------
!> Create a new set of fields, allocate space for them, and initialize to zero
subroutine soca_fields_create(self, geom, vars)
  class(soca_fields),        intent(inout) :: self
  type(soca_geom),  pointer, intent(inout) :: geom
  type(oops_variables),         intent(in) :: vars  !< list of field names to create

  character(len=:), allocatable :: vars_str(:)
  integer :: i

  ! make sure current object has not already been allocated
  if (associated(self%fields)) &
    call abor1_ftn("soca_fields::create(): object already allocated")

  ! associate geometry
  self%geom => geom

  ! initialize the variable parameters
  allocate(character(len=1024) :: vars_str(vars%nvars()))
  do i=1,vars%nvars()
    vars_str(i) = trim(vars%variable(i))
  end do
  call soca_fields_init_vars(self, vars_str)

  ! TODO: generalize derived type fields as well (mld, layer_depth, etc) ?
  allocate(self%mld(geom%isd:geom%ied,geom%jsd:geom%jed))
  allocate(self%layer_depth(geom%isd:geom%ied,geom%jsd:geom%jed,geom%nzo))

  ! set everything to zero
  call self%zeros()
end subroutine soca_fields_create


! ------------------------------------------------------------------------------
!> delete all the fields
subroutine soca_fields_delete(self)
  class(soca_fields), intent(inout) :: self
  integer :: i

  ! clear the fields and nullify pointers
  nullify(self%geom)
  do i = 1, size(self%fields)
    call self%fields(i)%delete()
  end do
  deallocate(self%fields)
  nullify(self%fields)

  ! delete derived variables
  ! TODO generalize
  deallocate(self%mld)
  deallocate(self%layer_depth)
end subroutine


! ------------------------------------------------------------------------------
!> Copy the contents of rhs to self. Self will be initialized with the variable
!> names in rhs if not already initialized
subroutine soca_fields_copy(self, rhs)
  class(soca_fields), intent(inout) :: self
  type(soca_fields),  intent(in)    :: rhs

  character(len=:), allocatable :: vars_str(:)
  integer :: i

  ! initialize the variables based on the names in rhs
  if (.not. associated(self%fields)) then
    self%geom => rhs%geom
    allocate(character(len=1024) :: vars_str(size(rhs%fields)))
    do i=1, size(vars_str)
      vars_str(i) = rhs%fields(i)%name
    end do
    call soca_fields_init_vars(self, vars_str)
  end if

  ! copy values from rhs to self
  do i=1,size(self%fields)
    call self%fields(i)%copy(rhs%fields(i))
  end do

  ! copy (and allocate) the derived variables
  ! TODO generalize
  if (.not. allocated(self%mld)) then
    allocate(self%mld(self%geom%isd:self%geom%ied,self%geom%jsd:self%geom%jed))
    allocate(self%layer_depth(self%geom%isd:self%geom%ied,self%geom%jsd:self%geom%jed,self%geom%nzo))
  end if
  self%mld = rhs%mld
  self%layer_depth = rhs%layer_depth
end subroutine


! ------------------------------------------------------------------------------
!> get a pointer to the soca_field with the given name.
!> If no field exists with that name, the prorgam aborts
!> (use soca_fields%has() if you need to check for optional fields)
subroutine soca_fields_get(self, name, field)
  class(soca_fields),         intent(in) :: self
  character(len=*),           intent(in) :: name   !< name of field to find
  type(soca_field), pointer, intent(out) :: field  !< a pointer to the resulting field

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
!> returns whether a field with the given name exists
function soca_fields_has(self, name) result(res)
  class(soca_fields), intent(in) :: self
  character(len=*),   intent(in) :: name

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
!> reset all fields to zero
subroutine soca_fields_zeros(self)
  class(soca_fields), intent(inout) :: self
  integer :: i

  do i = 1, size(self%fields)
    self%fields(i)%val = 0.0_kind_real
  end do

  ! TODO generalize
  self%mld = 0.0_kind_real
  self%layer_depth = 0.0_kind_real
end subroutine soca_fields_zeros


! ------------------------------------------------------------------------------
! TODO, generalize by removing the hardcoded int=>field_name
subroutine soca_fields_dirac(self, f_conf)
  class(soca_fields),         intent(inout) :: self
  type(fckit_configuration), intent(in)    :: f_conf   !< Configuration

  integer :: isc, iec, jsc, jec
  integer :: ndir,n, z
  integer,allocatable :: ixdir(:),iydir(:),izdir(:),ifdir(:)
  type(fckit_mpi_comm) :: f_comm

  type(soca_field), pointer :: field

  ! Get MPI communicator
  f_comm = fckit_mpi_comm()

  ! Get Diracs size
  ndir = f_conf%get_size("ixdir")
  if (( f_conf%get_size("iydir") /= ndir ) .or. &
      ( f_conf%get_size("izdir") /= ndir ) .or. &
      ( f_conf%get_size("ifdir") /= ndir )) &
      call abor1_ftn('soca_fields_dirac: inconsistent sizes for ixdir, iydir, izdir, and ifdir')

  ! Allocation
  allocate(ixdir(ndir))
  allocate(iydir(ndir))
  allocate(izdir(ndir))
  allocate(ifdir(ndir))

  ! Get Diracs positions
  call f_conf%get_or_die("ixdir", ixdir)
  call f_conf%get_or_die("iydir", iydir)
  call f_conf%get_or_die("izdir", izdir)
  call f_conf%get_or_die("ifdir", ifdir)

  ! get PE domain bounds
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  ! Setup Diracs
  call self%zeros()
  do n=1,ndir
     ! skip this index if not in the bounds of this PE
     if (ixdir(n) > iec .or. ixdir(n) < isc) cycle
     if (iydir(n) > jec .or. iydir(n) < jsc) cycle

    field => null()
    select case(ifdir(n))
    case (1)
      call self%get("tocn", field)
    case (2)
      call self%get("socn", field)
    case (3)
      call self%get("ssh", field)
    case (4)
      call self%get("cicen", field)
    case (5)
      call self%get("hicen", field)
    case default
      ! TODO print error that out of range
    end select
    if (associated(field)) then
      z = 1
      if (field%nz > 1) z = izdir(n)
      field%val(ixdir(n),iydir(n),izdir(n)) = 1.0
    end if
  end do
end subroutine soca_fields_dirac


! ------------------------------------------------------------------------------
!> initialize fields with random normal distribution
subroutine soca_fields_random(self)
  class(soca_fields), intent(inout) :: self
  integer, parameter :: rseed = 1 ! constant for reproducability of tests
    ! NOTE: random seeds are not quite working the way expected,
    !  it is only set the first time normal_distribution() is called with a seed
  integer :: z, i

  type(soca_field), pointer :: field

  ! set random values
  do i = 1, size(self%fields)
    field => self%fields(i)
    ! NOTE: can't randomize "hocn", testIncrementInterpAD fails
    if (field%name == 'hocn') cycle
    call normal_distribution(field%val,  0.0_kind_real, 1.0_kind_real, rseed)
  end do

  ! mask out land, set to zero
  do i=1,size(self%fields)
    field => self%fields(i)
    if (.not. associated(field%mask) ) cycle
    do z=1,field%nz
      field%val(:,:,z) = field%val(:,:,z) * field%mask(:,:)
    end do
  end do

  ! update domains
  do i=1, size(self%fields)
    field => self%fields(i)
    call mpp_update_domains(field%val, self%geom%Domain%mpp_domain)
  end do
end subroutine soca_fields_random


! ------------------------------------------------------------------------------
!> add two sets of fields together
subroutine soca_fields_add(self,rhs)
  class(soca_fields), intent(inout) :: self
  type(soca_fields),     intent(in) :: rhs
  integer :: i

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! add
  do i=1,size(self%fields)
    self%fields(i)%val = self%fields(i)%val + rhs%fields(i)%val
  end do
end subroutine soca_fields_add


! ------------------------------------------------------------------------------
!> perform a shur product between two sets of fields
subroutine soca_fields_schur(self,rhs)
  class(soca_fields), intent(inout) :: self
  type(soca_fields),     intent(in) :: rhs
  integer :: i

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! schur product
  do i=1,size(self%fields)
    self%fields(i)%val = self%fields(i)%val * rhs%fields(i)%val
  end do
end subroutine soca_fields_schur


! ------------------------------------------------------------------------------
!> subtract two sets of fields
subroutine soca_fields_sub(self,rhs)
  class(soca_fields), intent(inout) :: self
  type(soca_fields),     intent(in) :: rhs
  integer :: i

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! subtract
  do i=1,size(self%fields)
    self%fields(i)%val = self%fields(i)%val - rhs%fields(i)%val
  end do
end subroutine soca_fields_sub


! ------------------------------------------------------------------------------
!> multiply a set of fields by a constant
subroutine soca_fields_mul(self,zz)
  class(soca_fields), intent(inout) :: self
  real(kind=kind_real),  intent(in) :: zz
  integer :: i

  do i=1,size(self%fields)
    self%fields(i)%val = zz * self%fields(i)%val
  end do
end subroutine soca_fields_mul


! ------------------------------------------------------------------------------
!> Add two fields (multiplying the rhs first)
subroutine soca_fields_axpy(self,zz,rhs)
  class(soca_fields), intent(inout) :: self
  real(kind=kind_real),  intent(in) :: zz
  type(soca_fields),     intent(in) :: rhs
  integer :: i

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  do i=1,size(self%fields)
    self%fields(i)%val = self%fields(i)%val + zz* rhs%fields(i)%val
  end do
end subroutine soca_fields_axpy

! ------------------------------------------------------------------------------
!> calculate the global dot product of two sets of fields
subroutine soca_fields_dotprod(fld1,fld2,zprod)
  class(soca_fields),     intent(in) :: fld1
  type(soca_fields),      intent(in) :: fld2
  real(kind=kind_real),  intent(out) :: zprod

  real(kind=kind_real) :: local_zprod
  integer :: ii, jj, kk, n
  type(fckit_mpi_comm) :: f_comm
  type(soca_field), pointer :: field1, field2

  ! make sure fields are same shape
  call fld1%check_congruent(fld2)

  ! loop over (almost) all fields
  local_zprod = 0.0_kind_real
  do n=1,size(fld1%fields)
    field1 => fld1%fields(n)
    field2 => fld2%fields(n)

    ! determine which fields to use
    ! TODO: use all fields (this will change answers in the ctests)
    select case(field1%name)
    case ("tocn","socn","ssh","uocn","vocn",&
          "hicen", "sw", "lhf", "shf", "lw", "us", "cicen")
      continue
    case default
      cycle
    end select

    ! add the given field to the dot product (only using the compute domain)
    do ii = fld1%geom%isc, fld1%geom%iec
      do jj = fld1%geom%jsc, fld1%geom%jec
        ! TODO masking is wrong, but this will change answers, should be:
        ! if (field1%masked .and. .not. fld1%geom%mask2d(ii,jj) == 1 ) cycle
        if (.not. fld1%geom%mask2d(ii,jj) == 1 ) cycle
        do kk=1,field1%nz
          local_zprod = local_zprod + field1%val(ii,jj,kk) * field2%val(ii,jj,kk)
        end do
      end do
    end do
  end do

  ! Get global dot product
  f_comm = fckit_mpi_comm()
  call f_comm%allreduce(local_zprod, zprod, fckit_mpi_sum())
end subroutine soca_fields_dotprod


! ------------------------------------------------------------------------------
!> add a set of increments to the set of fields
subroutine soca_fields_add_incr(self,rhs)
  class(soca_fields), intent(inout) :: self
  type(soca_fields),  intent(in)    :: rhs

  integer, save :: cnt_outer = 1
  character(len=800) :: filename, str_cnt
  type(soca_field), pointer :: fld, fld_r
  integer :: i, k

  real(kind=kind_real) :: min_ice = 1e-6_kind_real
  real(kind=kind_real) :: amin = 1e-6_kind_real
  real(kind=kind_real) :: amax = 10.0_kind_real
  real(kind=kind_real), allocatable :: alpha(:,:), aice_bkg(:,:), aice_ana(:,:)
  type(soca_geostrophy_type) :: geostrophy
  type(soca_fields) :: incr_geo
  type(soca_field), pointer :: t, s, u, v, h, dt, ds, du, dv

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! for each field
  do i=1,size(self%fields)
    fld => self%fields(i)
    fld_r => rhs%fields(i)

    select case (fld%name)
    case ('cicen')
      ! NOTE: sea ice concentration is special

      ! compute background rescaling
      allocate(alpha, mold=fld_r%val(:,:,1))
      allocate(aice_bkg, mold=alpha)
      allocate(aice_ana, mold=alpha)
      aice_bkg  = sum(fld%val(:,:,:), dim=3)
      aice_ana  = aice_bkg + sum(fld_r%val(:,:,:), dim=3)
      where (aice_ana < 0.0_kind_real) aice_ana = 0.0_kind_real
      where (aice_ana > 1.0_kind_real) aice_ana = 1.0_kind_real
      alpha = 1.0_kind_real
      where ( aice_bkg > min_ice) alpha = aice_ana / aice_bkg

      ! limit size of increment
      where ( alpha > amax ) alpha = amax
      where ( alpha < amin ) alpha = amin

      ! add fraction of increment
      do k=1,fld%nz
        fld%val(:,:,k) = alpha * fld%val(:,:,k)
      end do

    case default
      ! everyone else is normal
      fld%val = fld%val + fld_r%val
    end select
  end do

  ! Save increment for outer loop cnt_outer
  write(str_cnt,*) cnt_outer
  filename='incr.'//adjustl(trim(str_cnt))//'.nc'

  ! Compute geostrophic increment
  if (self%has('hocn').and.self%has('tocn').and.self%has('socn').and.&
       self%has('uocn').and.self%has('vocn')) then
     ! Get necessary background fields needed to compute geostrophic perturbations
     call self%get("tocn", t)
     call self%get("socn", s)
     call self%get("hocn", h)

     ! Make a copy of the increment and get the needed pointers
     call incr_geo%copy(rhs)
     call incr_geo%get("tocn", dt)
     call incr_geo%get("socn", ds)
     call incr_geo%get("uocn", du)
     call incr_geo%get("vocn", dv)

     ! Compute the geostrophic increment
     call geostrophy%setup(self%geom, h%val)
     call geostrophy%tl(h%val, t%val, s%val,&
          dt%val, ds%val, du%val, dv%val, self%geom)
     call geostrophy%delete()

     ! Add geostrophic increment to state
     call self%get("uocn", u)
     call self%get("vocn", v)
     where( (du%val**2 + dv%val**2)<0.25)
        u%val = u%val + du%val
        v%val = v%val + dv%val
     end where

     ! Save increment and clean memory
     call incr_geo%write_file(filename)
     call incr_geo%delete()
  else
     call rhs%write_file(filename)
  end if

  ! Update outer loop counter
  cnt_outer = cnt_outer + 1

end subroutine soca_fields_add_incr


! ------------------------------------------------------------------------------
!> subtract two sets of fields, saving the results separately
subroutine soca_fields_diff_incr(self,x1,x2)
  class(soca_fields), intent(inout) :: self
  type(soca_fields),  intent(in)    :: x1
  type(soca_fields),  intent(in)    :: x2
  integer :: i

  ! make sure fields are same shape
  call self%check_congruent(x1)
  call self%check_congruent(x2)

  ! subtract
  do i=1,size(self%fields)
    self%fields(i)%val = x1%fields(i)%val - x2%fields(i)%val
  end do
end subroutine soca_fields_diff_incr


! ------------------------------------------------------------------------------

subroutine soca_fields_read(fld, f_conf, vdate)
  class(soca_fields),         intent(inout) :: fld     !< Fields
  type(fckit_configuration), intent(in)    :: f_conf  !< Configuration
  type(datetime),            intent(inout) :: vdate   !< DateTime

  integer, parameter :: max_string_length=800
  character(len=max_string_length) :: ocn_filename, sfc_filename, filename
  character(len=:), allocatable :: basename, incr_filename
  real(kind=kind_real), allocatable :: cicen_val(:,:,:)
  integer :: iread = 0
  integer :: ii
  logical :: vert_remap=.false.
  character(len=max_string_length) :: remap_filename
  real(kind=kind_real), allocatable :: h_common(:,:,:)    !< layer thickness to remap to
  type(restart_file_type), target :: ocean_restart, sfc_restart, ice_restart
  type(restart_file_type) :: ocean_remap_restart
  type(restart_file_type), pointer :: restart
  integer :: idr
  integer :: isd, ied, jsd, jed
  integer :: isc, iec, jsc, jec
  integer :: i, j, nz, n
  type(remapping_CS)  :: remapCS
  character(len=:), allocatable :: str, seaice_model
  real(kind=kind_real), allocatable :: h_common_ij(:), hocn_ij(:), varocn_ij(:), varocn2_ij(:)
  logical :: read_sfc
  type(soca_field), pointer :: field, field2, hocn, cicen
  real(kind=kind_real) :: soca_rho_ice  = 905.0 !< [kg/m3]
  real(kind=kind_real) :: soca_rho_snow = 330.0 !< [kg/m3]

  if ( f_conf%has("read_from_file") ) &
      call f_conf%get_or_die("read_from_file", iread)

  call fld%get("hocn", hocn)

  ! Get Indices for data domain and allocate common layer depth array
  isd = fld%geom%isd ; ied = fld%geom%ied
  jsd = fld%geom%jsd ; jed = fld%geom%jed

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
          domain=fld%geom%Domain%mpp_domain)
     call restore_state(ocean_remap_restart, directory='')
     call free_restart_type(ocean_remap_restart)
     call fms_io_exit()
  end if

  ! iread = 0: Invent state
  if (iread==0) then
     call fld%zeros()
     call f_conf%get_or_die("date", str)
     call datetime_set(str, vdate)
  end if

  ! TODO redo this to be generic

  ! iread = 1 (state) or 3 (increment): Read restart file
  if ((iread==1).or.(iread==3)) then
    ! Read sea-ice
    seaice_model = ""
    if(f_conf%get("ice_filename", str)) then
      call f_conf%get_or_die("basename", basename)
      filename = trim(basename)//trim(str)
      if(.not. f_conf%get("seaice_model", seaice_model)) seaice_model = "sis2"
    end if

    select case(seaice_model)
    case ('sis2')
      call fms_io_init()
      do i=1,size(fld%fields)
        select case(fld%fields(i)%name)
        case ('cicen')
          allocate(cicen_val(isd:ied, jsd:jed, fld%fields(i)%nz + 1))
          cicen_val = 0.0_kind_real
          idr = register_restart_field(ice_restart, filename, 'part_size', &
                  cicen_val, &
                  domain=fld%geom%Domain%mpp_domain)
        case ('hicen')
          idr = register_restart_field(ice_restart, filename, 'h_ice', &
                  fld%fields(i)%val(:,:,:), &
                  domain=fld%geom%Domain%mpp_domain)
        case ('hsnon')
          idr = register_restart_field(ice_restart, filename, 'h_snow', &
                  fld%fields(i)%val(:,:,:), &
                  domain=fld%geom%Domain%mpp_domain)
        end select
      end do

      call restore_state(ice_restart, directory='')
      call free_restart_type(ice_restart)
      call fms_io_exit()
      do i=1,size(fld%fields)
        select case(fld%fields(i)%name)
        case ('cicen')
          fld%fields(i)%val = cicen_val(:,:,2:)
          deallocate(cicen_val)
        case ('hicen')
          fld%fields(i)%val = fld%fields(i)%val / soca_rho_ice
        case ('hsnon')
          fld%fields(i)%val = fld%fields(i)%val / soca_rho_snow
        end select
      end do

    case ('cice')
      call fms_io_init()
      do i=1,size(fld%fields)
        select case(fld%fields(i)%name)
        case ('cicen')
          idr = register_restart_field(ice_restart, filename, 'aicen', &
                  fld%fields(i)%val(:,:,:), &
                  domain=fld%geom%Domain%mpp_domain)
        case ('hicen')
          idr = register_restart_field(ice_restart, filename, 'vicen', &
                  fld%fields(i)%val(:,:,:), &
                  domain=fld%geom%Domain%mpp_domain)
        case ('hsnon')
          idr = register_restart_field(ice_restart, filename, 'vsnon', &
                  fld%fields(i)%val(:,:,:), &
                  domain=fld%geom%Domain%mpp_domain)
        end select
      end do
      call restore_state(ice_restart, directory='')
      call free_restart_type(ice_restart)
      call fms_io_exit()
      ! Convert to hicen and hsnon
      call fld%get("cicen", cicen)
      do i=1,size(fld%fields)
        select case(fld%fields(i)%name)
        case ('hicen','hsnon')
          where(cicen%val(:,:,:)>0.0_kind_real)
            fld%fields(i)%val  = fld%fields(i)%val / cicen%val(:,:,:)
          end where
        end select
      end do

    end select

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

    call fms_io_init()

    ! built-in variables
    do i=1,size(fld%fields)
      if(fld%fields(i)%io_name /= "") then
        ! which file are we reading from?
        select case(fld%fields(i)%io_file)
        case ('ocn')
          filename = ocn_filename
          restart => ocean_restart
        case ('sfc')
          if (sfc_filename == "") cycle ! we have sfc fields, but no file to read from
          filename = sfc_filename
          restart => sfc_restart
          read_sfc = .true.
        case default
          call abor1_ftn('read_file(): illegal io_file: '//fld%fields(i)%io_file)
        end select

      ! setup to read
        if (fld%fields(i)%nz == 1) then
          idr = register_restart_field(restart, filename, fld%fields(i)%io_name, &
              fld%fields(i)%val(:,:,1), domain=fld%geom%Domain%mpp_domain)
        else
          idr = register_restart_field(restart, filename, fld%fields(i)%io_name, &
              fld%fields(i)%val(:,:,:), domain=fld%geom%Domain%mpp_domain)
        end if
      end if
    end do

    call restore_state(ocean_restart, directory='')
    call free_restart_type(ocean_restart)
    if (read_sfc) then
      call restore_state(sfc_restart, directory='')
      call free_restart_type(sfc_restart)
    end if
    call fms_io_exit()

    ! Indices for compute domain
    isc = fld%geom%isc ; iec = fld%geom%iec
    jsc = fld%geom%jsc ; jec = fld%geom%jec

    ! Initialize mid-layer depth from layer thickness
    call fld%geom%thickness2depth(hocn%val, fld%layer_depth)

    ! Compute mixed layer depth TODO: Move somewhere else ...
    call fld%get("tocn", field)
    call fld%get("socn", field2)
    do i = isc, iec
      do j = jsc, jec
          fld%mld(i,j) = soca_mld(&
              &field2%val(i,j,:),&
              &field%val(i,j,:),&
              &fld%layer_depth(i,j,:),&
              &fld%geom%lon(i,j),&
              &fld%geom%lat(i,j))
      end do
    end do

    ! Remap layers if needed
    if (vert_remap) then
      allocate(h_common_ij(nz), hocn_ij(nz), varocn_ij(nz), varocn2_ij(nz))
      call initialize_remapping(remapCS,'PCM')
      do i = isc, iec
        do j = jsc, jec
          h_common_ij = h_common(i,j,:)
          hocn_ij = hocn%val(i,j,:)

          do n=1,size(fld%fields)
            field => fld%fields(n)
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
    do n=1,size(fld%fields)
      field => fld%fields(n)
      call mpp_update_domains(field%val, fld%geom%Domain%mpp_domain)
    end do

    ! Set vdate if reading state
    if (iread==1) then
      call f_conf%get_or_die("date", str)
      call datetime_set(str, vdate)
    end if

    return
  end if

  ! Read diagnostic file
  if (iread==2) then
    call f_conf%get_or_die("filename", str)
    incr_filename = str
    call fms_io_init()
    do ii = 1, size(fld%fields)
      field => fld%fields(ii)
      if (field%io_name == "" ) cycle
      if ( field%nz == 1) then
        call read_data(incr_filename, field%io_name, field%val(:,:,1), domain=fld%geom%Domain%mpp_domain)
      else
        call read_data(incr_filename, field%io_name, field%val(:,:,:), domain=fld%geom%Domain%mpp_domain)
      end if
    end do
    call fms_io_exit()
  endif
end subroutine soca_fields_read


! ------------------------------------------------------------------------------
!> calculate global statistics for each field (min, max, average)
subroutine soca_fields_gpnorm(fld, nf, pstat)
  class(soca_fields),      intent(in) :: fld
  integer,                 intent(in) :: nf
  real(kind=kind_real), intent(inout) :: pstat(3, nf) !> [min, max, average]

  logical :: mask(fld%geom%isc:fld%geom%iec, fld%geom%jsc:fld%geom%jec)
  real(kind=kind_real) :: ocn_count, local_ocn_count, tmp(3)
  integer :: jj, isc, iec, jsc, jec
  type(fckit_mpi_comm) :: f_comm
  type(soca_field), pointer :: field

  f_comm = fckit_mpi_comm()

  ! Indices for compute domain
  isc = fld%geom%isc ; iec = fld%geom%iec
  jsc = fld%geom%jsc ; jec = fld%geom%jec

  ! get the number of ocean grid cells
  local_ocn_count = sum(fld%geom%mask2d(isc:iec, jsc:jec))
  call f_comm%allreduce(local_ocn_count, ocn_count, fckit_mpi_sum())
  mask = fld%geom%mask2d(isc:iec,jsc:jec) > 0.0

  ! calculate global min, max, mean for each field
  do jj=1, size(fld%fields)

    ! get local min/max/sum of each variable
    ! TODO: use all fields (this will change answers in the ctests)
    select case(fld%fields(jj)%name)
    case("tocn", "socn", "ssh", "hocn", "uocn", "vocn", &
         "sw", "lw", "lhf", "shf", "us", "hicen", "hsnon", "cicen")
      continue
    case default
      cycle
    end select

    ! TODO only mask fields that should be masked (will change answers)
    call fld%get(fld%fields(jj)%name, field)
    call fldinfo(field%val(isc:iec,jsc:jec,:), mask, tmp)

    ! calculate global min/max/mean
    call f_comm%allreduce(tmp(1), pstat(1,jj), fckit_mpi_min())
    call f_comm%allreduce(tmp(2), pstat(2,jj), fckit_mpi_max())
    call f_comm%allreduce(tmp(3), pstat(3,jj), fckit_mpi_sum())
    pstat(3,jj) = pstat(3,jj)/ocn_count
  end do
end subroutine soca_fields_gpnorm


! ------------------------------------------------------------------------------
! TODO remove hardcoded number of variables
subroutine ug_size(self, ug)
  type(soca_fields),          intent(in) :: self
  type(unstructured_grid), intent(inout) :: ug

  integer :: isc, iec, jsc, jec
  integer :: igrid

  ! Indices for compute domain (no halo)
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  ! Set number of grids
  if (ug%colocated==1) then
     ! Colocated
     ug%ngrid = 1
  else
     ! Not colocated
     ug%ngrid = 2
  end if

  ! Allocate grid instances
  if (.not.allocated(ug%grid)) allocate(ug%grid(ug%ngrid))

  do igrid=1,ug%ngrid
     ! Set local number of points
     ug%grid(igrid)%nmga = (iec - isc + 1) * (jec - jsc + 1 )

     ! Set number of timeslots
     ug%grid(igrid)%nts = ug%nts
  end do

  if (ug%colocated==1) then
     ! Set number of levels
     ug%grid(1)%nl0 = self%geom%nzo

     ! Set number of variables
     ug%grid(1)%nv = self%geom%ncat*2 + 3
  else
     ! Set number of levels
     ug%grid(1)%nl0 = self%geom%nzo
     ug%grid(2)%nl0 = 1

     ! Set number of variables
     ug%grid(1)%nv = 2
     ug%grid(2)%nv = self%geom%ncat*2 + 1
  end if

end subroutine ug_size

! ------------------------------------------------------------------------------

subroutine soca_fields_ug_coord(self, ug)
  class(soca_fields), intent(in) :: self
  type(unstructured_grid), intent(inout) :: ug

  integer :: igrid
  integer :: isc, iec, jsc, jec, jz

  ! Indices for compute domain (no halo)
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  ! Define size
  call ug_size(self, ug)

  ! Allocate unstructured grid coordinates
  call allocate_unstructured_grid_coord(ug)

  ! Define coordinates for 3D grid
  igrid = 1
  ug%grid(igrid)%lon = &
       &reshape( self%geom%lon(isc:iec, jsc:jec), (/ug%grid(igrid)%nmga/) )
  ug%grid(igrid)%lat = &
       &reshape( self%geom%lat(isc:iec, jsc:jec), (/ug%grid(igrid)%nmga/) )
  ug%grid(igrid)%area = &
       &reshape( self%geom%cell_area(isc:iec, jsc:jec), (/ug%grid(igrid)%nmga/) )
  do jz = 1, ug%grid(igrid)%nl0
     ug%grid(igrid)%vunit(:,jz) = real(jz)
     ug%grid(igrid)%lmask(:,jz) = reshape( self%geom%mask2d(isc:iec, jsc:jec)==1, (/ug%grid(igrid)%nmga/) )
  end do

  if (ug%colocated==0) then
     ! Define coordinates for 2D grid
     igrid = 2
     ug%grid(igrid)%lon = &
          &reshape( self%geom%lon(isc:iec, jsc:jec), (/ug%grid(igrid)%nmga/) )
     ug%grid(igrid)%lat = &
          &reshape( self%geom%lat(isc:iec, jsc:jec), (/ug%grid(igrid)%nmga/) )
     ug%grid(igrid)%area = &
          &reshape( self%geom%cell_area(isc:iec, jsc:jec), (/ug%grid(igrid)%nmga/) )
     ug%grid(igrid)%vunit(:,1) = 0.0_kind_real
     ug%grid(igrid)%lmask(:,1) = reshape( self%geom%mask2d(isc:iec, jsc:jec)==1, (/ug%grid(igrid)%nmga/) )
  end if

end subroutine soca_fields_ug_coord

! ------------------------------------------------------------------------------
! TODO generalize to use all vars
subroutine soca_fields_to_ug(self, ug, its)
  class(soca_fields),         intent(in) :: self
  type(unstructured_grid), intent(inout) :: ug
  integer,                    intent(in) :: its

  integer :: isc, iec, jsc, jec, jk, igrid
  integer :: ni, nj, i, z

  type(soca_field), pointer :: field

  ! Indices for compute domain (no halo)
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  ni = iec - isc + 1
  nj = jec - jsc + 1

  ! Allocate unstructured grid field
  call ug_size(self, ug)
  call allocate_unstructured_grid_field(ug)

  ! Copy 3D field
  igrid = 1
  jk = 1
  ug%grid(igrid)%fld(:,:,:,its) = 0.0_kind_real
  do i=1,size(self%fields)
    field => self%fields(i)
    select case(field%name)
    case ('tocn','socn')
      do z = 1, field%nz
        ug%grid(igrid)%fld(1:ni*nj, z, jk, its) = &
          &reshape( field%val(isc:iec, jsc:jec,z), (/ug%grid(igrid)%nmga/) )
      end do
      jk = jk + 1
    end select
  end do

  if (ug%colocated==1) then
     ! 2D variables copied as 3D variables
     igrid = 1
  else
     ! 2D variables copied as 2D variables
     igrid = 2
     jk = 1
  end if

  ! 2d fields
  do i=1,size(self%fields)
    field => self%fields(i)
    select case(field%name)
    case ('hicen', 'cicen', 'ssh')
      do z = 1, field%nz
        ug%grid(igrid)%fld(1:ni*nj, 1, jk, its) = &
            &reshape( field%val(isc:iec, jsc:jec, z), (/ug%grid(igrid)%nmga/) )
        jk = jk + 1
      end do
    end select
  end do

end subroutine soca_fields_to_ug

! ------------------------------------------------------------------------------
! Generalize variable names used
subroutine soca_fields_from_ug(self, ug, its)
  class(soca_fields),   intent(inout) :: self
  type(unstructured_grid), intent(in) :: ug
  integer,                 intent(in) :: its

  integer :: isc, iec, jsc, jec, jk, igrid
  integer :: ni, nj, i, z
  type(soca_field), pointer :: field

  ! Indices for compute domain (no halo)
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  ni = iec - isc + 1
  nj = jec - jsc + 1

  igrid = 1
  jk = 1
  call self%zeros()

  ! 3d fields
  do i=1,size(self%fields)
    field => self%fields(i)
    select case(field%name)
    case ("tocn", "socn")
      do z = 1, field%nz
        field%val(isc:iec, jsc:jec,z) = &
          &reshape( ug%grid(igrid)%fld(1:ni*nj, z, jk, its), (/ni, nj/) )
      end do
      jk = jk + 1
    end select
  end do

  if (ug%colocated==1) then
     ! 2D variables copied as 3D variables
     igrid = 1
  else
     ! 2D variables copied as 2D variables
     igrid = 2
     jk = 1
  end if

  ! 2d fields
  do i=1,size(self%fields)
    field => self%fields(i)
    select case(field%name)
    case ('hicen', 'cicen', 'ssh')
      do z = 1, field%nz
        field%val(isc:iec, jsc:jec, z) = &
          &reshape( ug%grid(igrid)%fld(1:ni*nj, 1, jk, its), (/ni, nj/) )
        jk = jk + 1
      end do
    end select
  end do

end subroutine soca_fields_from_ug


! ------------------------------------------------------------------------------
!> make sure two sets of fields are the same shape
!> (same variables, same resolution)
!> TODO: make this more robust (allow for different number of fields?)
subroutine soca_fields_check_congruent(f1, f2)
  class(soca_fields), intent(in) :: f1, f2

  integer :: i, j

  ! number of fields should be the same
  if (size(f1%fields) /= size(f2%fields)) &
    call abor1_ftn("soca_fields: contains different number of fields")

  ! each field should match (name, size, shape)
  do i=1,size(f1%fields)
    if (f1%fields(i)%name /= f2%fields(i)%name) &
      call abor1_ftn("soca_fields: field have different names")
    do j = 1, size(shape(f1%fields(i)%val))
      if (size(f1%fields(i)%val, dim=j) /= size(f2%fields(i)%val, dim=j) ) then
        call abor1_ftn("soca_fields: field '"//f1%fields(i)%name//"' has different dimensions")
      end if
    end do
  end do
end subroutine soca_fields_check_congruent


! ------------------------------------------------------------------------------
!> Save soca fields to file using fms write_data
subroutine soca_fields_write_file(fld, filename)
  class(soca_fields),  intent(in) :: fld    !< Fields
  character(len=*),   intent(in) :: filename

  integer :: ii

  call fms_io_init()
  call set_domain( fld%geom%Domain%mpp_domain )

  ! write out all fields
  do ii = 1, size(fld%fields)
    call write_data( filename, fld%fields(ii)%name, fld%fields(ii)%val(:,:,:), fld%geom%Domain%mpp_domain)
  end do

  ! some other derived fields that should be written out
  call write_data( filename, "rossby_radius", fld%geom%rossby_radius, fld%geom%Domain%mpp_domain)

  call fms_io_exit()
end subroutine soca_fields_write_file

! ------------------------------------------------------------------------------
!> Save soca fields in a restart format
!> TODO this can be generalized even more
subroutine soca_fields_write_rst(fld, f_conf, vdate)
  class(soca_fields),         intent(inout) :: fld      !< Fields
  type(fckit_configuration), intent(in)    :: f_conf   !< Configuration
  type(datetime),            intent(inout) :: vdate    !< DateTime

  integer, parameter :: max_string_length=800
  character(len=:), allocatable :: seaice_model
  character(len=max_string_length) :: ocn_filename, sfc_filename, ice_filename, filename
  type(restart_file_type), target :: ocean_restart, sfc_restart, ice_restart
  type(restart_file_type), pointer :: restart
  integer :: idr, i
  type(soca_field), pointer :: field, cicen
  real(kind=kind_real), allocatable :: cicen_val(:,:,:), vicen(:,:,:), vsnon(:,:,:)
  logical :: write_sfc, write_ice

  write_ice = .false.
  write_sfc = .false.
  call fms_io_init()

  ! filenames
  ocn_filename = soca_genfilename(f_conf,max_string_length,vdate,"ocn")
  sfc_filename = soca_genfilename(f_conf,max_string_length,vdate,"sfc")

  ! Check what ice model we are writing a file to ('sis2' or 'cice')
  if ( .not. f_conf%get("seaice_model", seaice_model)) seaice_model = 'sis2'
  if ( seaice_model == 'sis2') then
    ice_filename = soca_genfilename(f_conf, max_string_length,vdate,"ice")
  else
    ice_filename = soca_genfilename(f_conf, max_string_length,vdate,"cice")
  end if

  ! built in variables
  do i=1,size(fld%fields)
    field => fld%fields(i)
    ! TODO move the ice calculations elsewhere
    ! these variable are handled specially (...annoying ice)
    if (field%name == "cicen" .and. seaice_model /= 'cice') then
      allocate(cicen_val(fld%geom%isd:fld%geom%ied, fld%geom%jsd:fld%geom%jed, field%nz+1))
      cicen_val(:,:,2:) = field%val
      cicen_val(:,:,1) = 1.0 - sum(cicen_val(:,:,2:), dim=3)
      idr = register_restart_field( ice_restart, ice_filename, "part_size", &
        cicen_val, domain=fld%geom%Domain%mpp_domain)

    else if (seaice_model == "cice" .and. field%name == "hicen") then
      allocate(vicen, mold=field%val)
      call fld%get("cicen", cicen)
      vicen = cicen%val * field%val
      idr = register_restart_field( ice_restart, ice_filename, "vicen", &
        vicen, domain=fld%geom%Domain%mpp_domain)

    else if (seaice_model == "cice" .and. field%name == "hsnon") then
      allocate(vsnon, mold=field%val)
      call fld%get("cicen", cicen)
      vsnon = cicen%val * field%val
      idr = register_restart_field( ice_restart, ice_filename, "vsnon", &
        vsnon, domain=fld%geom%Domain%mpp_domain)

    else if (field%io_file /= "") then
      ! which file are we writing to
      select case(field%io_file)
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
        ! TODO move io_name for ice variables into config file, since they
        ! depend on which model is used
        if ( seaice_model == 'sis2') then
          select case(field%name)
          case ('hicen')
            field%io_name = 'h_ice'
          case ('hsnon')
            field%io_name = 'h_snow'
          end select
        else if ( seaice_model == 'cice') then
          select case(field%name)
          case ('cicen')
            field%io_name = 'aicen'
          end select
        end if
      case default
        call abor1_ftn('soca_write_restart(): illegal io_file: '//field%io_file)
      end select

      ! write
      if (field%nz == 1) then
        idr = register_restart_field( restart, filename, field%io_name, &
          field%val(:,:,1), domain=fld%geom%Domain%mpp_domain)
      else
        idr = register_restart_field( restart, filename, field%io_name, &
        field%val(:,:,:), domain=fld%geom%Domain%mpp_domain)
      end if
    end if
  end do

  ! derived variables
  idr = register_restart_field(ocean_restart, ocn_filename, 'mld', fld%mld(:,:), &
       domain=fld%geom%Domain%mpp_domain)

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
  call fms_io_exit()

end subroutine soca_fields_write_rst

! ------------------------------------------------------------------------------

subroutine soca_fields_getpoint(self, geoiter, values)
  class(soca_fields),             intent(   in) :: self
  type(soca_geom_iter),           intent(   in) :: geoiter
  real(kind=kind_real),           intent(inout) :: values(:)

  integer :: ff, ii, nz
  type(soca_field), pointer :: field

  ! get values
  ! TODO generalize field names
  ii = 0
  do ff = 1, size(self%fields)
    field => self%fields(ff)
    select case(field%name)
    case("tocn", "socn", "ssh", "hocn", "cicen", "hicen","hsnon")
      nz = field%nz
      values(ii+1:ii+nz) = field%val(geoiter%iind, geoiter%jind,:)
      ii = ii + nz
    end select
  end do
end subroutine soca_fields_getpoint

! ------------------------------------------------------------------------------

subroutine soca_fields_setpoint(self, geoiter, values)
  class(soca_fields),             intent(inout) :: self
  type(soca_geom_iter),           intent(   in) :: geoiter
  real(kind=kind_real),           intent(   in) :: values(:)

  integer :: ff, ii, nz
  type(soca_field), pointer :: field

  ! Set values
  ! TODO generalize field names
  ii = 0
  do ff = 1, size(self%fields)
    field => self%fields(ff)
    select case(field%name)
    case("tocn", "socn", "ssh", "hocn", "cicen", "hicen","hsnon")
      nz = field%nz
      field%val(geoiter%iind, geoiter%jind,:) = values(ii+1:ii+nz)
      ii = ii + nz
    end select
  end do
end subroutine soca_fields_setpoint

end module soca_fields_mod
