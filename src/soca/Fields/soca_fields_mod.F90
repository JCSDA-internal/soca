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
use soca_ocnsfc_mod, only: soca_ocnsfc_type
use soca_seaice_mod, only: soca_seaice_type
use soca_utils, only: soca_mld

implicit none

private
public :: soca_fields, soca_field, &
          dirac, random,  &
          self_add, self_schur, self_sub, self_mul, axpy, &
          dot_prod, add_incr, diff_incr, &
          read_file, write_file, gpnorm, fldrms, soca_fld2file, &
          change_resol, check, &
          field_to_ug, field_from_ug, ug_coord, &
          soca_getpoint, soca_setpoint

! interface create
!    module procedure create_constructor, create_copy
! end interface create

! ------------------------------------------------------------------------------
!> Fortran derived type to hold fields
type :: soca_field
  character(len=:),     allocatable :: name
  !character(len=:),     allocatable :: io_name
  !character(len=:),     allocatable :: io_file
  !character(len=:),     allocatable :: vgrid
  !character(len=:),     allocatable :: hgrid
  integer                           :: nz
  logical :: masked = .false.
  real(kind=kind_real), allocatable :: val(:,:,:)  
contains 
  procedure :: delete => soca_field_delete
  procedure :: copy   => soca_field_copy
end type soca_field


type :: soca_fields
   type(soca_geom), pointer          :: geom           !< MOM6 Geometry
   integer                           :: nf             !< Number of fields
   !character(len=128)                :: gridfname      !< Grid file name
   !character(len=128)                :: cicefname      !< Fields file name for cice
   !character(len=128)                :: momfname       !< Fields file name for mom

   type(soca_field),     pointer :: fields(:) => null()

   ! Sea-ice state variables
   type(soca_seaice_type)            :: seaice         !< Sea-ice state

   real(kind=kind_real), allocatable :: hocn(:,:,:)    !< DA layer thickness (nx,ny,nzo)

   ! Ocean diagnostics
   real(kind=kind_real), allocatable :: mld(:,:)           !< Mixed layer depth (nx,ny)
   real(kind=kind_real), allocatable :: layer_depth(:,:,:) !< Mid-layer depth (nx,ny,nz0)

   ! Ocean surface fields
   type(soca_ocnsfc_type)            :: ocnsfc !< Surface fields needed for cool skin ufo

   character(len=5),     allocatable :: fldnames(:)    !< Variable identifiers             (nf)

contains 
  procedure :: create => soca_fields_create
  procedure :: copy   => soca_fields_copy  
  procedure :: delete => soca_fields_delete  
  procedure :: zeros  => soca_fields_zeros
  !procedure :: has    => soca_fields_has
  procedure :: get    => soca_fields_get
end type soca_fields

contains

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

subroutine soca_field_copy(self, rhs)
  ! TODO, check to see if allocated??
  class(soca_field), intent(inout) :: self
  type(soca_field),  intent(in)    :: rhs

  integer :: i

  if ( self%nz /= rhs%nz ) call abor1_ftn("soca_field::copy():  self%nz /= rhs%nz")

  ! make sure val array sizes are congruent
  if ( size(shape(self%val)) /= size(shape(rhs%val)) ) call abor1_ftn("soca_field::copy():  shape of self%val /= rhs%val")
  do i =1, size(shape(self%val))
    if (size(self%val, dim=i) /= size(rhs%val, dim=i)) &
      call abor1_ftn("soca_field::copy():  shape of self%val /= rhs%val")
  end do

  self%name = rhs%name
  self%nz = rhs%nz
  self%val = rhs%val

end subroutine

! ------------------------------------------------------------------------------

subroutine soca_field_delete(self)
  class(soca_field), intent(inout) :: self

  ! TODO, is this really necessary?
  deallocate(self%name)
  deallocate(self%val)
end subroutine

! ------------------------------------------------------------------------------

subroutine soca_fields_create(self, geom, vars)
  class(soca_fields),        intent(inout) :: self
  type(soca_geom),  pointer, intent(inout) :: geom
  type(oops_variables),         intent(in) :: vars

  integer :: i, nz

  self%geom => geom
  ! todo allocate extra for internal fields??
  allocate(self%fields(vars%nvars()))
  do i=1,vars%nvars()
    self%fields(i)%name = trim(vars%variable(i))
    
    ! determine number of levels, and allocate space
    select case(self%fields(i)%name)
    case ('tocn','socn','hocn')
      nz = geom%nzo
    case ('cicen','hicen')
      nz = geom%nzi
    case ('hsnon')
      nz = geom%nzs
    case ('ssh', 'sw', 'lhf', 'shf', 'lw', 'us')
      nz = 1
    case default
      call abor1_ftn('soca_fields::create(): unknown field '// self%fields(i)%name)
    end select
    self%fields(i)%masked = .true.
    self%fields(i)%nz = nz
    allocate(self%fields(i)%val(&
      geom%isd:geom%ied, &
      geom%jsd:geom%jed, &
      nz ))
  end do  

  ! TODO delete the following line eventually
  call create_constructor(self, geom, vars)

  call self%zeros()

end subroutine

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

subroutine soca_fields_delete(self)
  class(soca_fields), intent(inout) :: self
  
  integer :: i

  nullify(self%geom)
  do i = 1, size(self%fields)
    call self%fields(i)%delete()
  end do
  deallocate(self%fields)
  nullify(self%fields)

  ! TODO delete this line
  call delete(self)

end subroutine


! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

subroutine soca_fields_copy(self, rhs)
  ! Construct a field from an other field
  class(soca_fields), intent(inout) :: self
  type(soca_fields),  intent(in)    :: rhs
  
  integer :: i, nz


  ! TODO combine the following with create
  if (.not. associated(self%fields)) then
      self%geom => rhs%geom
    ! todo allocate extra for internal fields??
    allocate(self%fields(size(rhs%fields)))
    do i=1, size(rhs%fields)
      self%fields(i)%name = rhs%fields(i)%name
      
      ! determine number of levels, and allocate space
      select case(self%fields(i)%name)
      case ('tocn','socn','hocn')
        nz =  self%geom%nzo
      case ('cicen','hicen')
        nz =  self%geom%nzi
      case ('hsnon')
        nz =  self%geom%nzs
      case ('ssh', 'sw', 'lhf', 'shf', 'lw', 'us')
        nz = 1
      case default
        call abor1_ftn('soca_fields::create(): unknown field '// self%fields(i)%name)
      end select
      self%fields(i)%masked = .true.
      self%fields(i)%nz = nz
      allocate(self%fields(i)%val(&
        self%geom%isd: self%geom%ied, &
        self%geom%jsd: self%geom%jed, &
        nz ))
    end do    
  end if 

  ! copy values
  do i=1,size(self%fields)
    call self%fields(i)%copy(rhs%fields(i))
  end do

  ! TODO delete the following lines
  if (.not. allocated(self%fldnames)) then
    call create_copy(self, rhs)
  end if
  call copy(self, rhs)

end subroutine 

! ------------------------------------------------------------------------------

subroutine soca_fields_get(self, name, field)
  class(soca_fields),         intent(in) :: self
  character(len=*),           intent(in) :: name
  type(soca_field), pointer, intent(out) :: field
  
  integer :: i

  do i=1,size(self%fields)
    if (trim(name) == self%fields(i)%name) then
      field => self%fields(i)
      return
    end if
  end do

  call abor1_ftn("soca_fields::get():  cannot find field "//trim(name))

end subroutine

! ------------------------------------------------------------------------------
!> Create a field from geometry and variables
subroutine create_constructor(self, geom, vars)
  type(soca_fields),         intent(inout) :: self
  type(soca_geom),  pointer, intent(inout) :: geom
  type(oops_variables),         intent(in) :: vars

  integer :: i

  ! Allocate
  call soca_field_alloc(self, geom)

  ! Associate geometry
  self%geom => geom

  ! Set fields numbers and names
  self%nf   = vars%nvars()
  allocate(self%fldnames(self%nf))
  do i=1,self%nf
    self%fldnames(i)=vars%variable(i)
  end do 

  call check(self)

end subroutine create_constructor

! ------------------------------------------------------------------------------

subroutine create_copy(self, rhs_fld)
  ! Construct a field from an other field, lhs_fld=rhs_fld
  type(soca_fields), intent(inout) :: self
  type(soca_fields), intent(in)    :: rhs_fld

  ! Allocate and copy fields
  call soca_field_alloc(self, rhs_fld%geom)

  ! Associate geometry
  self%geom => rhs_fld%geom

  ! Set fields numbers and names
  self%nf   = rhs_fld%nf
  allocate(self%fldnames(self%nf))
  self%fldnames(:)=rhs_fld%fldnames(:)

  call check(self)

end subroutine create_copy

! ------------------------------------------------------------------------------

subroutine soca_field_alloc(self, geom)
  type (soca_fields),    intent(inout) :: self
  type(soca_geom), pointer, intent(in) :: geom

  integer :: isd, ied, jsd, jed, nzo

  ! Short cut to geometry
  isd = geom%isd ; ied = geom%ied
  jsd = geom%jsd ; jed = geom%jed
  nzo = geom%nzo

  ! Allocate ocean state
  allocate(self%hocn(isd:ied,jsd:jed,nzo))
  allocate(self%mld(isd:ied,jsd:jed))
  allocate(self%layer_depth(isd:ied,jsd:jed,nzo))

  ! Allocate sea-ice state
  call self%seaice%create(geom)

  ! Allocate surface fields for cool skin
  call self%ocnsfc%create(geom)

  call self%zeros()
end subroutine soca_field_alloc

! ------------------------------------------------------------------------------

subroutine delete(self)
  type (soca_fields), intent(inout) :: self

  ! Deallocate ocean state
  deallocate(self%hocn)
  deallocate(self%mld)
  deallocate(self%layer_depth)

  ! Deallocate sea-ice state
  call self%seaice%delete()


  ! Deallocate surface fields
  call self%ocnsfc%delete()

  ! Deassociate geometry
  nullify(self%geom)

end subroutine delete


! ------------------------------------------------------------------------------

subroutine soca_fields_zeros(self)
  class(soca_fields), intent(inout) :: self

  integer :: i
    
  do i = 1, size(self%fields)
    self%fields(i)%val = 0.0_kind_real
  end do
  
  ! TODO delete the following
  self%hocn = 0.0_kind_real
  self%mld = 0.0_kind_real

  call self%seaice%zeros()
  call self%ocnsfc%zeros()

end subroutine soca_fields_zeros

! ------------------------------------------------------------------------------

subroutine dirac(self, f_conf)
  type(soca_fields),         intent(inout) :: self
  type(fckit_configuration), intent(in)    :: f_conf   !< Configuration

  integer :: isc, iec, jsc, jec
  integer :: ndir,n, z
  integer,allocatable :: ixdir(:),iydir(:),izdir(:),ifdir(:)
  type(fckit_mpi_comm) :: f_comm

  type(soca_field), pointer :: field 

  call check(self)

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
    case default
      ! TODO print error that out range     
     if (ifdir(n)==4) self%seaice%cicen(ixdir(n),iydir(n),izdir(n)) = 1.0
     if (ifdir(n)==5) self%seaice%hicen(ixdir(n),iydir(n),izdir(n)) = 1.0
    end select
    if (associated(field)) then
      z = 1
      if (field%nz > 1) z = izdir(n)
      field%val(ixdir(n),iydir(n),izdir(n)) = 1.0
    end if

  end do

end subroutine dirac

! ------------------------------------------------------------------------------

subroutine random(self)
  type(soca_fields), intent(inout) :: self
  integer, parameter :: rseed = 1 ! constant for reproducability of tests
    ! NOTE: random seeds are not quite working the way expected,
    !  it is only set the first time normal_distribution() is called with a seed
  integer :: z, ff, i

  type(soca_field), pointer :: field

  call check(self)

  ! set random values
  do i = 1, size(self%fields)
    field => self%fields(i)
    select case(field%name)
    case("tocn", "socn", "ssh")
     call normal_distribution(field%val,  0.0_kind_real, 1.0_kind_real, rseed)
    !case("hocn")
    ! NOTE: can't randomize "hocn", testIncrementInterpAD fails
    end select
  end do

  ! mask out land, set to zero
  do i=1,size(self%fields)
    field => self%fields(i)
    if (.not. field%masked ) cycle
    do z=1,field%nz
      field%val(:,:,z) = field%val(:,:,z) * self%geom%mask2d
    end do    
  end do

  ! update domains
  do i=1, size(self%fields)
    field => self%fields(i)
    call mpp_update_domains(field%val, self%geom%Domain%mpp_domain)
  end do  

  ! do the same for the non-ocean fields
  call self%ocnsfc%random(self%fldnames)
  call self%seaice%random(self%fldnames)

end subroutine random

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)
  type(soca_fields), intent(inout) :: self
  type(soca_fields),    intent(in) :: rhs

  call check_resolution(self, rhs)

  !nf = common_vars(self, rhs)

  ! Associate geometry
  if (.not.associated(self%geom)) self%geom => rhs%geom

  ! Set fields numbers and names
  self%nf   = rhs%nf
  if (.not.allocated(self%fldnames)) allocate(self%fldnames(self%nf))
  self%fldnames(:)=rhs%fldnames(:)

  self%hocn  = rhs%hocn
  self%mld   = rhs%mld
  self%layer_depth   = rhs%layer_depth

  ! Sea-ice
  call self%seaice%copy(rhs%seaice)

  ! Ocean surface
  call self%ocnsfc%copy(rhs%ocnsfc)

  return
end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)
  type(soca_fields), intent(inout) :: self
  type(soca_fields),    intent(in) :: rhs

  integer :: i

  call check_resolution(self, rhs)

  self%hocn = self%hocn + rhs%hocn

  call self%seaice%add(rhs%seaice)
  call self%ocnsfc%add(rhs%ocnsfc)

  do i=1,size(self%fields)
    ! TODO, check to make sure size/name is the same
    self%fields(i)%val = self%fields(i)%val + rhs%fields(i)%val
  end do
end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)
  type(soca_fields), intent(inout) :: self
  type(soca_fields),    intent(in) :: rhs

  integer :: i

  call check_resolution(self, rhs)

  self%hocn=self%hocn*rhs%hocn

  call self%seaice%schur(rhs%seaice)
  call self%ocnsfc%schur(rhs%ocnsfc)

  do i=1,size(self%fields)
    ! TODO, check to make sure size/name is the same
    self%fields(i)%val = self%fields(i)%val * rhs%fields(i)%val
  end do
end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)
  type(soca_fields), intent(inout) :: self
  type(soca_fields),    intent(in) :: rhs

  integer :: i
  call check_resolution(self, rhs)

  self%hocn=self%hocn-rhs%hocn

  call self%seaice%sub(rhs%seaice)
  call self%ocnsfc%sub(rhs%ocnsfc)

  do i=1,size(self%fields)
    ! TODO, check to make sure size/name is the same
    self%fields(i)%val = self%fields(i)%val - rhs%fields(i)%val
  end do

end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)
  type(soca_fields), intent(inout) :: self
  real(kind=kind_real), intent(in) :: zz

  integer :: i

  call check(self)

  self%hocn = zz * self%hocn

  call self%seaice%mul(zz)
  call self%ocnsfc%mul(zz)

  do i=1,size(self%fields)
    ! TODO, check to make sure size/name is the same
    self%fields(i)%val = zz * self%fields(i)%val
  end do

end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)
  type(soca_fields), intent(inout) :: self
  real(kind=kind_real), intent(in) :: zz
  type(soca_fields),    intent(in) :: rhs

  integer :: i
  
  call check_resolution(self, rhs)

  self%hocn = self%hocn + zz * rhs%hocn

  call self%seaice%axpy(zz, rhs%seaice)
  call self%ocnsfc%axpy(zz, rhs%ocnsfc)

  do i=1,size(self%fields)
    ! TODO, check to make sure size/name is the same
    self%fields(i)%val = self%fields(i)%val + zz* rhs%fields(i)%val
  end do

end subroutine axpy

! ------------------------------------------------------------------------------

subroutine dot_prod(fld1,fld2,zprod)
  type(soca_fields),     intent(in) :: fld1
  type(soca_fields),     intent(in) :: fld2
  real(kind=kind_real), intent(out) :: zprod

  real(kind=kind_real) :: zprod_allpes
  integer :: ii, jj, kk, ff, n
  integer :: isc, iec, jsc, jec
  integer :: ncat, nzo, myrank
  type(fckit_mpi_comm) :: f_comm

  type(soca_field), pointer :: field1, field2


  ! Setup Communicator
  f_comm = fckit_mpi_comm()

  call check_resolution(fld1, fld2)
  if (fld1%nf /= fld2%nf .or. fld1%geom%nzo /= fld2%geom%nzo) then
     call abor1_ftn("soca_fields:field_prod error number of fields")
  endif

  ! Indices for compute domain (no halo)
  isc = fld1%geom%isc ; iec = fld1%geom%iec
  jsc = fld1%geom%jsc ; jec = fld1%geom%jec

  ! Get ice categories and ocean levels
  ncat = fld1%geom%ncat
  nzo = fld1%geom%nzo

  zprod = 0.0_kind_real

  !----- OCEAN
  do n=1,size(fld1%fields)
    field1 => fld1%fields(n)
    call fld2%get(field1%name, field2)
    select case(field1%name)
    case ("tocn","socn","ssh")
      do ii = isc, iec
        do jj = jsc, jec
          if (.not. fld1%geom%mask2d(ii,jj) == 1 ) cycle
          do kk=1,field1%nz
            zprod = zprod + field1%val(ii,jj,kk) * field2%val(ii,jj,kk)
          end do
         end do
      end do
    end select
  end do
  call f_comm%barrier() ! why is this barrier here???
  myrank = f_comm%rank()

  !----- SEA-ICE. TODO: Move to seaice module
  do ii = isc, iec
     do jj = jsc, jec
        do kk = 1, ncat
           if (fld1%geom%mask2d(ii,jj)==1) then
              zprod = zprod + fld1%seaice%cicen(ii,jj,kk+1)*fld2%seaice%cicen(ii,jj,kk+1) & !CICEN
                   + fld1%seaice%hicen(ii,jj,kk)*fld2%seaice%hicen(ii,jj,kk)       !HICEN
           end if
        end do
     end do
  end do

!!$    !----- OCEAN Surface. TODO: Move to ocnsfc module
  do ii = isc, iec
     do jj = jsc, jec
        if (fld1%geom%mask2d(ii,jj)==1) then
           zprod = zprod + fld1%ocnsfc%sw_rad(ii,jj)*fld2%ocnsfc%sw_rad(ii,jj) &
                         + fld1%ocnsfc%lw_rad(ii,jj)*fld2%ocnsfc%lw_rad(ii,jj) &
                         + fld1%ocnsfc%latent_heat(ii,jj)*fld2%ocnsfc%latent_heat(ii,jj) &
                         + fld1%ocnsfc%sens_heat(ii,jj)*fld2%ocnsfc%sens_heat(ii,jj) &
                         + fld1%ocnsfc%fric_vel(ii,jj)*fld2%ocnsfc%fric_vel(ii,jj)
        end if
     end do
  end do

  ! Get global dot product
  call f_comm%allreduce(zprod, zprod_allpes, fckit_mpi_sum())
  zprod = zprod_allpes

end subroutine dot_prod

! ------------------------------------------------------------------------------

subroutine add_incr(self,rhs)
  type(soca_fields), intent(inout) :: self
  type(soca_fields), intent(in)    :: rhs

  integer, save :: cnt_outer = 1
  character(len=800) :: filename, str_cnt

  integer :: i
  call check(self)
  call check(rhs)

  ! Add increment to field
  self%hocn = self%hocn + rhs%hocn

  call self%seaice%add_incr(rhs%seaice)
  call self%ocnsfc%add(rhs%ocnsfc)

  do i=1,size(self%fields)
    ! TODO, check to make sure size/name is the same
    self%fields(i)%val = self%fields(i)%val + rhs%fields(i)%val
  end do

  ! Save increment for outer loop cnt_outer
  write(str_cnt,*) cnt_outer
  filename='incr.'//adjustl(trim(str_cnt))//'.nc'
  call soca_fld2file(rhs, filename)

  ! Update outer loop counter
  cnt_outer = cnt_outer + 1

end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine diff_incr(self,x1,x2)
  type(soca_fields), intent(inout) :: self
  type(soca_fields), intent(in)    :: x1
  type(soca_fields), intent(in)    :: x2
  integer :: i

  call check(self)
  call check(x1)
  call check(x2)

  call self%zeros()

  self%hocn = x1%hocn - x2%hocn

  call self%seaice%diff_incr(x1%seaice, x2%seaice)
  call self%ocnsfc%diff_incr(x1%ocnsfc, x2%ocnsfc)

  do i=1,size(self%fields)
    ! TODO, check to make sure size/name is the same
    self%fields(i)%val = x1%fields(i)%val - x2%fields(i)%val
  end do
end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(fld,rhs)
  type(soca_fields), intent(inout) :: fld
  type(soca_fields), intent(in)    :: rhs

  call check(fld)
  call check(rhs)
  call fld%copy(rhs)
  call fld%ocnsfc%copy(rhs%ocnsfc)
  call fld%seaice%copy(rhs%seaice)

end subroutine change_resol

! ------------------------------------------------------------------------------

subroutine read_file(fld, f_conf, vdate)
  type(soca_fields),         intent(inout) :: fld     !< Fields
  type(fckit_configuration), intent(in)    :: f_conf  !< Configuration
  type(datetime),            intent(inout) :: vdate   !< DateTime

  integer, parameter :: max_string_length=800
  character(len=max_string_length) :: ocn_filename
  character(len=max_string_length) :: basename, incr_filename
  character(len=1024) :: buf
  integer :: iread = 0
  integer :: ii
  logical :: vert_remap=.false.
  character(len=max_string_length) :: remap_filename
  real(kind=kind_real), allocatable :: h_common(:,:,:)    !< layer thickness to remap to
  type(restart_file_type) :: ocean_restart
  type(restart_file_type) :: ocean_remap_restart
  integer :: idr_ocean
  integer :: isd, ied, jsd, jed
  integer :: isc, iec, jsc, jec
  integer :: i, j, nz, n
  type(remapping_CS)  :: remapCS
  character(len=:), allocatable :: str
  real(kind=kind_real), allocatable :: h_common_ij(:), hocn_ij(:), tsocn_ij(:), tsocn2_ij(:)

  type(soca_field), pointer :: field, field2

  if ( f_conf%has("read_from_file") ) &
      call f_conf%get_or_die("read_from_file", iread)

  ! Check if vertical remapping needs to be applied
  if ( f_conf%has("remap_filename") ) then
     vert_remap = .true.
     call f_conf%get_or_die("remap_filename", str)
     remap_filename = str

     ! Get Indices for data domain and allocate common layer depth array
     isd = fld%geom%isd ; ied = fld%geom%ied
     jsd = fld%geom%jsd ; jed = fld%geom%jed

     nz=size(fld%hocn, dim=3)
     allocate(h_common(isd:ied,jsd:jed,nz))

     ! Read common vertical coordinate from file
     call fms_io_init()
     idr_ocean = register_restart_field(ocean_remap_restart, remap_filename, 'h', fld%hocn(:,:,:), &
          domain=fld%geom%Domain%mpp_domain)
     call restore_state(ocean_remap_restart, directory='')
     call fms_io_exit()
     h_common = fld%hocn

  end if

  ! iread = 0: Invent state
  if (iread==0) then
     call fld%zeros()
     call f_conf%get_or_die("date", str)
     call datetime_set(str, vdate)
  end if

  ! iread = 1 (state) or 3 (increment): Read restart file
  if ((iread==1).or.(iread==3)) then
     ! Read ocean surface fields
     call fld%ocnsfc%read_restart(f_conf, fld%geom, fld%fldnames)

     ! Read sea-ice
     call fld%seaice%read_restart(f_conf, fld%geom, fld%fldnames)

     call f_conf%get_or_die("basename", str)
     basename = str
     call f_conf%get_or_die("ocn_filename", str)
     ocn_filename = trim(basename) // trim(str)

     call fms_io_init()
     do ii = 1, fld%nf
        select case(fld%fldnames(ii))
           ! Ocean
        case ('ssh')
          call fld%get("ssh", field)
           idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'ave_ssh', field%val(:,:,1), &
                domain=fld%geom%Domain%mpp_domain)
        case ('tocn')
          call fld%get("tocn", field)
           idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'Temp', field%val, &
                domain=fld%geom%Domain%mpp_domain)
        case ('socn')
          call fld%get("socn", field)
           idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'Salt', field%val, &
                domain=fld%geom%Domain%mpp_domain)
        case ('hocn')
           idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'h', fld%hocn(:,:,:), &
                domain=fld%geom%Domain%mpp_domain)
        case default
           call log%warning("soca_fields:read_file: Not reading var "//fld%fldnames(ii))
        end select
     end do
     call restore_state(ocean_restart, directory='')
     call fms_io_exit()

     ! Indices for compute domain
     isc = fld%geom%isc ; iec = fld%geom%iec
     jsc = fld%geom%jsc ; jec = fld%geom%jec

     ! Initialize mid-layer depth from layer thickness
     call fld%geom%thickness2depth(fld%hocn, fld%layer_depth)

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
        allocate(h_common_ij(nz), hocn_ij(nz), tsocn_ij(nz), tsocn2_ij(nz))
        call initialize_remapping(remapCS,'PCM')
        do i = isc, iec
           do j = jsc, jec
              if (fld%geom%mask2d(i,j).eq.1) then
                 h_common_ij = h_common(i,j,:)
                 hocn_ij = fld%hocn(i,j,:)

                 call fld%get("tocn", field)
                 tsocn_ij = field%val(i,j,:)
                 call remapping_core_h(remapCS, nz, h_common_ij, tsocn_ij,&
                      &nz, hocn_ij, tsocn2_ij)
                 field%val(i,j,:) = tsocn2_ij

                 call fld%get("socn", field)
                 tsocn_ij = field%val(i,j,:)
                 call remapping_core_h(remapCS, nz, h_common_ij, tsocn_ij,&
                      &nz, hocn_ij, tsocn2_ij)
                 field%val(i,j,:) = tsocn2_ij

              else
                call fld%get("tocn", field)
                field%val(i,j,:) = 0.0_kind_real
                call fld%get("socn", field)
                field%val(i,j,:) = 0.0_kind_real
              end if
           end do
        end do
        fld%hocn = h_common
        deallocate(h_common_ij, hocn_ij, tsocn_ij, tsocn2_ij)
     end if
     call end_remapping(remapCS)

     ! Update halo
     do n=1,size(fld%fields)
      field => fld%fields(n)
      select case(field%name)
      case ("tocn", "socn", "ssh")
        if (field%nz == 1) then
          call mpp_update_domains(field%val(:,:,1), fld%geom%Domain%mpp_domain)
        else
          call mpp_update_domains(field%val, fld%geom%Domain%mpp_domain)
        end if
      end select
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
     ! Read ocean surface fields
     call fld%ocnsfc%read_diag(f_conf, fld%geom, fld%fldnames)

     ! Read sea-ice fields
     call fld%ocnsfc%read_diag(f_conf, fld%geom, fld%fldnames)

     call f_conf%get_or_die("filename", str)
     incr_filename = str
     call fms_io_init()
     do ii = 1, fld%nf
        select case(fld%fldnames(ii))
           ! Ocean variables
        case ('ssh')
          call fld%get("ssh", field)          
           call read_data(incr_filename,"ssh",field%val(:,:,1),domain=fld%geom%Domain%mpp_domain)
        case ('tocn')
          call fld%get("tocn", field)
          call read_data(incr_filename,"temp",field%val(:,:,:),domain=fld%geom%Domain%mpp_domain)
        case ('socn')
          call fld%get("socn", field)
          call read_data(incr_filename,"salt",field%val(:,:,:),domain=fld%geom%Domain%mpp_domain)
        case ('hocn')
           call read_data(incr_filename,"h",fld%hocn(:,:,:),domain=fld%geom%Domain%mpp_domain)
        case default
           write(buf,*) 'soca_fields_mod::read_file::increment. Not reading '//fld%fldnames(ii)
           call log%info(buf,newl=.true.)
        end select
     end do
     call fms_io_exit()
  endif

  call check(fld)

  if (allocated(str)) deallocate(str)

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(fld, f_conf, vdate)
  type(soca_fields),         intent(inout) :: fld    !< Fields
  type(fckit_configuration), intent(in)    :: f_conf !< Configuration
  type(datetime),            intent(inout) :: vdate  !< DateTime

  integer, parameter :: max_string_length=800    ! Yuk!

  call check(fld)

  call soca_write_restart(fld, f_conf, vdate)

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine gpnorm(fld, nf, pstat)
  type(soca_fields),       intent(in) :: fld
  integer,                 intent(in) :: nf
  real(kind=kind_real), intent(inout) :: pstat(3, nf) !> [min, max, average]

  logical :: mask(fld%geom%isc:fld%geom%iec, fld%geom%jsc:fld%geom%jec)
  real(kind=kind_real) :: ocn_count, tmp(3)
  real(kind=kind_real) :: local_ocn_count
  integer :: jj
  integer :: isc, iec, jsc, jec
  type(fckit_mpi_comm) :: f_comm

  type(soca_field), pointer :: field  

  ! Setup Communicator
  f_comm = fckit_mpi_comm()

  call check(fld)

  ! Indices for compute domain
  isc = fld%geom%isc ; iec = fld%geom%iec
  jsc = fld%geom%jsc ; jec = fld%geom%jec

  ! get the number of ocean grid cells
  local_ocn_count = sum(fld%geom%mask2d(isc:iec, jsc:jec))
  call f_comm%allreduce(local_ocn_count, ocn_count, fckit_mpi_sum())
  mask = fld%geom%mask2d(isc:iec,jsc:jec) > 0.0


  ! calculate global min, max, mean for each field
  ! NOTE: "cicen" category 1 (no ice) is not included in the stats
  do jj=1, fld%nf
    tmp=0.0

    ! get local min/max/sum of each variable
    select case(fld%fldnames(jj))
    case("tocn", "socn", "ssh")
      call fld%get(fld%fldnames(jj), field)      
      call fldinfo(field%val(isc:iec,jsc:jec,:), mask, tmp)
    case("hocn")
      call fldinfo(fld%hocn(isc:iec,jsc:jec,:), mask, tmp)
    case("hicen")
       call fldinfo(fld%seaice%hicen(isc:iec,jsc:jec,:), mask, tmp)
    case("hsnon")
       call fldinfo(fld%seaice%hsnon(isc:iec,jsc:jec,:),  mask, tmp)
    case("cicen")
      call fldinfo(fld%seaice%cicen(isc:iec,jsc:jec,2:),  mask, tmp)
    case("sw")
      call fldinfo(fld%ocnsfc%sw_rad(isc:iec,jsc:jec),   mask, tmp)
    case("lw")
      call fldinfo(fld%ocnsfc%lw_rad(isc:iec,jsc:jec),   mask, tmp)
    case("lhf")
      call fldinfo(fld%ocnsfc%latent_heat(isc:iec,jsc:jec), mask, tmp)
    case("shf")
      call fldinfo(fld%ocnsfc%sens_heat(isc:iec,jsc:jec),   mask, tmp)
    case("us")
      call fldinfo(fld%ocnsfc%fric_vel(isc:iec,jsc:jec),    mask, tmp)
    end select

    ! calculate global min/max/mean
    call f_comm%allreduce(tmp(1), pstat(1,jj), fckit_mpi_min())
    call f_comm%allreduce(tmp(2), pstat(2,jj), fckit_mpi_max())
    call f_comm%allreduce(tmp(3), pstat(3,jj), fckit_mpi_sum())
    pstat(3,jj) = pstat(3,jj)/ocn_count

  end do
end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine fldrms(fld, prms)
  type(soca_fields),     intent(in) :: fld
  real(kind=kind_real), intent(out) :: prms

  call check(fld)

  call dot_prod(fld,fld,prms) ! Global value
  prms=sqrt(prms)

end subroutine fldrms

! ------------------------------------------------------------------------------

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

subroutine ug_coord(self, ug)
  type(soca_fields), intent(in) :: self
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

end subroutine ug_coord

! ------------------------------------------------------------------------------

subroutine field_to_ug(self, ug, its)
  type(soca_fields),          intent(in) :: self
  type(unstructured_grid), intent(inout) :: ug
  integer,                    intent(in) :: its

  integer :: isc, iec, jsc, jec, jk, incat, inzo, ncat, nzo, igrid
  integer :: ni, nj, i

  type(soca_field), pointer :: field

  ! Indices for compute domain (no halo)
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  ni = iec - isc + 1
  nj = jec - jsc + 1
  ncat = self%geom%ncat

  ! Define size
  call ug_size(self, ug)

  ! Allocate unstructured grid field
  call allocate_unstructured_grid_field(ug)
  ncat = self%geom%ncat
  nzo = ug%grid(1)%nl0

  ! Copy 3D field
  igrid = 1
  jk = 1
  ug%grid(igrid)%fld(:,:,:,its) = 0.0_kind_real

  ! tocn
  do i=1,size(self%fields)
    select case(self%fields(i)%name)
    case ('tocn','socn')
      call self%get(self%fields(i)%name, field)
      do inzo = 1, nzo
        ug%grid(igrid)%fld(1:ni*nj, inzo, jk, its) = &
          &reshape( field%val(isc:iec, jsc:jec,inzo), (/ug%grid(igrid)%nmga/) )
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

  ! cicen
  do incat = 1, ncat
     ug%grid(igrid)%fld(1:ni*nj, 1, jk, its) = &
          &reshape( self%seaice%cicen(isc:iec, jsc:jec, incat+1), (/ug%grid(igrid)%nmga/) )
     jk = jk + 1
  end do

  ! hicen
  do incat = 1, ncat
     ug%grid(igrid)%fld(1:ni*nj, 1, jk, its) = &
          &reshape( self%seaice%hicen(isc:iec, jsc:jec, incat), (/ug%grid(igrid)%nmga/) )
     jk = jk + 1
  end do

  ! ssh
  call self%get("ssh", field)
  ug%grid(igrid)%fld(1:ni*nj, 1, jk, its) = &
       &reshape( field%val(isc:iec, jsc:jec, 1), (/ug%grid(igrid)%nmga/) )
  jk = jk + 1

end subroutine field_to_ug

! ------------------------------------------------------------------------------

subroutine field_from_ug(self, ug, its)
  type(soca_fields),    intent(inout) :: self
  type(unstructured_grid), intent(in) :: ug
  integer,                 intent(in) :: its

  integer :: isc, iec, jsc, jec, jk, incat, inzo, ncat, nzo, igrid
  integer :: ni, nj, i
  type(soca_field), pointer :: field


  ! Indices for compute domain (no halo)
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  ni = iec - isc + 1
  nj = jec - jsc + 1
  ncat = self%geom%ncat
  nzo = self%geom%nzo

  ! Copy 3D field
  igrid = 1
  jk = 1
  call self%zeros()

  ! tocn
  do i=1,size(self%fields)
    select case(self%fields(i)%name)
    case ("tocn", "socn")
      call self%get(self%fields(i)%name, field)
      do inzo = 1, nzo
        field%val(isc:iec, jsc:jec,inzo) = &
          &reshape( ug%grid(igrid)%fld(1:ni*nj, inzo, jk, its), (/ni, nj/) )
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

  ! cicen
  do incat = 1, ncat
     self%seaice%cicen(isc:iec, jsc:jec, incat+1) = &
          &reshape( ug%grid(igrid)%fld(1:ni*nj, 1, jk, its), (/ni, nj/))
     jk = jk + 1
  end do

  ! hicen
  do incat = 1, ncat
     self%seaice%hicen(isc:iec, jsc:jec, incat) = &
          &reshape( ug%grid(igrid)%fld(1:ni*nj, 1, jk, its), (/ni, nj/) )
     jk = jk + 1
  end do

  ! ssh
  call self%get("ssh", field)
  field%val(isc:iec, jsc:jec, 1) = &
       &reshape( ug%grid(igrid)%fld(1:ni*nj, 1, jk, its), (/ni, nj/) )
  jk = jk + 1

end subroutine field_from_ug

! ------------------------------------------------------------------------------

function common_vars(x1, x2)
  type(soca_fields), intent(in) :: x1, x2

  integer :: common_vars
  integer :: jf

  ! We assume here that one set of fields is a subset of the other,
  ! that fields are always in the same order starting with x,
  ! and that the common fields are the first ones.

  common_vars = min(x1%nf, x2%nf)
  do jf = 1, common_vars
     if (x1%fldnames(jf)/=x2%fldnames(jf)) &
          & call abor1_ftn("common_vars: fields do not match")
  enddo
  if (x1%geom%nzo /= x2%geom%nzo) call abor1_ftn("common_vars: error number of levels")
  !common_vars = x1%geom%nzi * common_vars

end function common_vars

! ------------------------------------------------------------------------------

subroutine check_resolution(x1, x2)
  type(soca_fields), intent(in) :: x1, x2

  ! NEEDS WORK !!!
  if (x1%geom%nx /= x2%geom%nx .or.  &
       &x1%geom%ny /= x2%geom%ny ) then
     call abor1_ftn ("soca_fields: resolution error")
  endif
  call check(x1)
  call check(x2)

end subroutine check_resolution

! ------------------------------------------------------------------------------

subroutine check(self)
  type(soca_fields), intent(in) :: self

  logical :: bad

  ! Doesn't do any thing ...
  bad = .false.
  !bad = bad .or. (size(self%cicen, 1) /= self%geom%nx)

  ! add more test here ...

  if (bad) then
     write(0,*)'nx, ny, nf, nzi, nzo = ',self%geom%nx,self%geom%ny
     call abor1_ftn ("soca_fields: field not consistent")
  endif

end subroutine check

! ------------------------------------------------------------------------------
!> Save soca fields to file using fms write_data
subroutine soca_fld2file(fld, filename)
  type(soca_fields),  intent(in) :: fld    !< Fields
  character(len=800), intent(in) :: filename

  integer :: ii
  character(len=800) :: fname

  type(soca_field), pointer :: field

  fname = trim(filename)

  call check(fld)

  call fms_io_init()
  call set_domain( fld%geom%Domain%mpp_domain )
  do ii = 1, fld%nf
     select case(fld%fldnames(ii))

     case ('ssh')
      call fld%get("ssh", field)
        call write_data( fname, "ssh", field%val(:,:,1), fld%geom%Domain%mpp_domain)
        call write_data( fname, "rossby_radius", fld%geom%rossby_radius, fld%geom%Domain%mpp_domain)
     case ('tocn')
      call fld%get("tocn", field)
      call write_data( fname, "temp", field%val, fld%geom%Domain%mpp_domain)
     case ('socn')
      call fld%get("socn", field)
      call write_data( fname, "salt", field%val, fld%geom%Domain%mpp_domain)
     case ('hocn')
        call write_data( fname, "h", fld%hocn, fld%geom%Domain%mpp_domain)
     case ('hicen')
        call write_data( fname, "hicen", fld%seaice%hicen, fld%geom%Domain%mpp_domain)
     case ('cicen')
        call write_data(fname, "cicen", fld%seaice%cicen, fld%geom%Domain%mpp_domain)
     case ('sw')
        call write_data(fname, "sw", fld%ocnsfc%sw_rad, fld%geom%Domain%mpp_domain)
     case ('lw')
        call write_data(fname, "lw", fld%ocnsfc%lw_rad, fld%geom%Domain%mpp_domain)
     case ('lhf')
        call write_data(fname, "lhf", fld%ocnsfc%latent_heat, fld%geom%Domain%mpp_domain)
     case ('shf')
        call write_data(fname, "shf", fld%ocnsfc%sens_heat, fld%geom%Domain%mpp_domain)
     case ('us')
        call write_data(fname, "us", fld%ocnsfc%fric_vel, fld%geom%Domain%mpp_domain)

     case default

     end select

  end do
  call fms_io_exit()

end subroutine soca_fld2file

! ------------------------------------------------------------------------------
!> Save soca fields in a restart format
subroutine soca_write_restart(fld, f_conf, vdate)
  type(soca_fields),         intent(inout) :: fld      !< Fields
  type(fckit_configuration), intent(in)    :: f_conf   !< Configuration
  type(datetime),            intent(inout) :: vdate    !< DateTime

  integer, parameter :: max_string_length=800
  character(len=max_string_length) :: ocn_filename, ocnsfc_filename
  type(restart_file_type) :: ocean_restart
  type(restart_file_type) :: ocnsfc_restart
  integer :: idr, idr_ocean
  type(soca_field), pointer :: field

  
  ! Generate file names
  ocn_filename = soca_genfilename(f_conf,max_string_length,vdate,"ocn")
  ocnsfc_filename = soca_genfilename(f_conf,max_string_length,vdate,"sfc")

  call fms_io_init()
  ! Ocean State
  call fld%get("ssh", field)
  idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'ave_ssh', field%val(:,:,1), &
       domain=fld%geom%Domain%mpp_domain)
  call fld%get("tocn", field)       
  idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'Temp', field%val(:,:,:), &
       domain=fld%geom%Domain%mpp_domain)
  call fld%get("socn", field)       
  idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'Salt', field%val(:,:,:), &
       domain=fld%geom%Domain%mpp_domain)
  idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'h', fld%hocn(:,:,:), &
       domain=fld%geom%Domain%mpp_domain)
  idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'mld', fld%mld(:,:), &
       domain=fld%geom%Domain%mpp_domain)

  ! Ocean-surface
  idr = register_restart_field(ocnsfc_restart, ocnsfc_filename, &
                              'sw_rad', fld%ocnsfc%sw_rad, &
                              domain=fld%geom%Domain%mpp_domain)
  idr = register_restart_field(ocnsfc_restart, ocnsfc_filename, &
                              'lw_rad', fld%ocnsfc%lw_rad, &
                              domain=fld%geom%Domain%mpp_domain)
  idr = register_restart_field(ocnsfc_restart, ocnsfc_filename, &
                              'latent_heat', fld%ocnsfc%latent_heat, &
                              domain=fld%geom%Domain%mpp_domain)
  idr = register_restart_field(ocnsfc_restart, ocnsfc_filename, &
                              'sens_heat', fld%ocnsfc%sens_heat, &
                              domain=fld%geom%Domain%mpp_domain)
  idr = register_restart_field(ocnsfc_restart, ocnsfc_filename, &
                              'fric_vel', fld%ocnsfc%fric_vel, &
                              domain=fld%geom%Domain%mpp_domain)

  call save_restart(ocean_restart, directory='')
  call save_restart(ocnsfc_restart, directory='')
  call free_restart_type(ocean_restart)
  call free_restart_type(ocnsfc_restart)
  call fms_io_exit()

  ! Save sea-ice restart
  call fld%seaice%write_restart(f_conf, fld%geom, vdate)

  return

end subroutine soca_write_restart


! ------------------------------------------------------------------------------

subroutine soca_getpoint(self, geoiter, values)

  type(soca_fields),              intent(   in) :: self
  type(soca_geom_iter),           intent(   in) :: geoiter
  real(kind=kind_real),           intent(inout) :: values(:)
  integer :: ff, ii, nzo, ncat, nz

  type(soca_field), pointer :: field
  nzo = self%geom%nzo
  ncat = self%geom%ncat

  ! get values
  ii = 0 
  do ff = 1, self%nf
    select case(self%fldnames(ff))
    case("tocn", "socn", "ssh")
      call self%get(self%fldnames(ff), field)
      nz = field%nz
      values(ii+1:ii+nz) = field%val(geoiter%iind, geoiter%jind,:)
      ii = ii + nz
    case("hocn")
      values(ii+1:ii+nzo) = self%hocn(geoiter%iind, geoiter%jind,:)
      ii = ii + nzo
    case("cicen")
      values(ii+1:ii+ncat+1) = self%seaice%cicen(geoiter%iind, geoiter%jind,:)
      ii = ii + ncat + 1
    case("hicen")
      values(ii+1:ii+ncat) = self%seaice%hicen(geoiter%iind, geoiter%jind,:)
      ii = ii + ncat
    case("hsnon")
      values(ii+1:ii+ncat) = self%seaice%hsnon(geoiter%iind, geoiter%jind,:)
      ii = ii + ncat
    end select
  end do

end subroutine soca_getpoint

! ------------------------------------------------------------------------------

subroutine soca_setpoint(self, geoiter, values)

  ! Passed variables
  type(soca_fields),              intent(inout) :: self
  type(soca_geom_iter),           intent(   in) :: geoiter
  real(kind=kind_real),           intent(   in) :: values(:)
  integer :: ff, ii, nzo, ncat, nz

  type(soca_field), pointer :: field


  nzo = self%geom%nzo
  ncat = self%geom%ncat

  ! Set values
  ii = 0
  do ff = 1, self%nf
    select case(self%fldnames(ff))
    case("tocn", "socn", "ssh")
      call self%get(self%fldnames(ff), field)
      nz = field%nz
      field%val(geoiter%iind, geoiter%jind,:) = values(ii+1:ii+nz)
      ii = ii + nz
    case("hocn")
      self%hocn(geoiter%iind, geoiter%jind,:) = values(ii+1:ii+nzo)
      ii = ii + nzo
    case("cicen")
      self%seaice%cicen(geoiter%iind, geoiter%jind,:) = values(ii+1:ii+ncat+1)
      ii = ii + ncat + 1
    case("hicen")
      self%seaice%hicen(geoiter%iind, geoiter%jind,:) = values(ii+1:ii+ncat)
      ii = ii + ncat
    case("hsnon")
      self%seaice%hsnon(geoiter%iind, geoiter%jind,:) = values(ii+1:ii+ncat)
      ii = ii + ncat
    end select
  end do

end subroutine soca_setpoint

end module soca_fields_mod
