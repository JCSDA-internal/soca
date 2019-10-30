! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_seaice_mod

use fckit_configuration_module, only: fckit_configuration
use fms_io_mod, only: fms_io_init, fms_io_exit, &
                      register_restart_field, restart_file_type, &
                      restore_state, free_restart_type, save_restart
use fms_mod, only: read_data
use kinds, only: kind_real
use datetime_mod, only: datetime
use random_mod, only: normal_distribution
use soca_geom_mod, only: soca_geom
use soca_fieldsutils_mod, only: soca_genfilename

implicit none

private
public :: soca_seaice_type

type :: soca_seaice_type
   type(soca_geom), pointer          :: geom

   ! Sea-ice state  variables
   real(kind=kind_real), allocatable :: cicen(:,:,:) !< Ice Fraction
   real(kind=kind_real), allocatable :: hicen(:,:,:) !< Ice thickness
   real(kind=kind_real), allocatable :: hsnon(:,:,:) !< Snow thickness
   real(kind=kind_real), allocatable :: vicen(:,:,:) !< Ice volume
   real(kind=kind_real), allocatable :: vsnon(:,:,:) !< Snow volume

   integer :: isd, ied, jsd, jed                 !< Data domain indices
   integer :: ncat                               !< Number of categories

   ! TODO: Get densities of ice and snow from config
   real(kind=kind_real) :: soca_rho_ice  = 905.0 !< [kg/m3]
   real(kind=kind_real) :: soca_rho_snow = 330.0 !< [kg/m3]
   
 contains
   procedure :: create => soca_seaice_create
   procedure :: delete => soca_seaice_delete
   procedure :: zeros => soca_seaice_zeros
   procedure :: abs => soca_seaice_abs
   procedure :: random => soca_seaice_random
   procedure :: copy => soca_seaice_copy
   procedure :: add => soca_seaice_add
   procedure :: add_incr => soca_seaice_add_incr
   procedure :: schur => soca_seaice_schur
   procedure :: sub => soca_seaice_sub
   procedure :: mul => soca_seaice_mul
   procedure :: axpy => soca_seaice_axpy
   procedure :: diff_incr => soca_seaice_diff_incr
   procedure :: read_restart => soca_seaice_read_rst
   procedure :: read_diag => soca_seaice_read_diag
   procedure :: write_restart => soca_seaice_write_rst
end type soca_seaice_type

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
subroutine soca_seaice_create(self, geom)
  class(soca_seaice_type), intent(inout) :: self
  type(soca_geom), pointer,   intent(in) :: geom

  integer :: isd, ied, jsd, jed, ncat

  ! Associate geometry
  self%geom => geom

  ! Indices for data domain (with halo)
  isd = geom%isd ; self%isd = isd
  ied = geom%ied ; self%ied = ied
  jsd = geom%jsd ; self%jsd = jsd
  jed = geom%jed ; self%jed = jed

  ! Get number of categories
  ncat = geom%ncat ; self%ncat = ncat

  ! Allocate sea-ice state
  if (.not.allocated(self%cicen)) allocate(self%cicen(isd:ied,jsd:jed,geom%ice_column%ncat+1))
  if (.not.allocated(self%hicen)) allocate(self%hicen(isd:ied,jsd:jed,geom%ice_column%ncat))
  if (.not.allocated(self%hsnon)) allocate(self%hsnon(isd:ied,jsd:jed,geom%ice_column%ncat))

end subroutine soca_seaice_create

! ------------------------------------------------------------------------------
subroutine soca_seaice_delete(self)
  class(soca_seaice_type), intent(inout) :: self

  ! Deallocate default sea-ice field
  deallocate(self%cicen)
  deallocate(self%hicen)
  deallocate(self%hsnon)

  ! Deallocate volumes
  if (allocated(self%vicen)) deallocate(self%vicen)
  if (allocated(self%vsnon)) deallocate(self%vsnon)

end subroutine soca_seaice_delete

! ------------------------------------------------------------------------------
subroutine soca_seaice_zeros(self)
  class(soca_seaice_type), intent(inout) :: self

  self%cicen = 0.0_kind_real
  self%hicen = 0.0_kind_real
  self%hsnon = 0.0_kind_real

end subroutine soca_seaice_zeros

! ------------------------------------------------------------------------------
subroutine soca_seaice_abs(self)
  class(soca_seaice_type), intent(inout) :: self

  self%cicen = abs(self%cicen)
  self%hicen = abs(self%hicen)
  self%hsnon = abs(self%hsnon)

end subroutine soca_seaice_abs

! ------------------------------------------------------------------------------
subroutine soca_seaice_random(self)
  class(soca_seaice_type), intent(inout) :: self
  integer :: i
  integer, parameter :: rseed = 1

  ! set random values
  call normal_distribution(self%cicen, 0.0_kind_real, 1.0_kind_real, rseed)
  call normal_distribution(self%hicen, 0.0_kind_real, 1.0_kind_real, rseed)
  call normal_distribution(self%hsnon, 0.0_kind_real, 1.0_kind_real, rseed)

  ! mask out land, set to zero
  do i=1, self%geom%ice_column%ncat+1
     self%cicen(:,:,i) = self%cicen(:,:,i) * self%geom%mask2d
  end do
  do i=1, self%geom%ice_column%ncat
     self%hicen(:,:,i) = self%hicen(:,:,i) * self%geom%mask2d
     self%hsnon(:,:,i) = self%hsnon(:,:,i) * self%geom%mask2d
  end do

  ! update domains
  ! TODO: if/when ever needed
end subroutine soca_seaice_random

! ------------------------------------------------------------------------------
subroutine soca_seaice_copy(self, rhs)
  class(soca_seaice_type), intent(inout) :: self
  class(soca_seaice_type),    intent(in) :: rhs

  ! Associate geometry
  self%geom => rhs%geom

  ! Copy fields 
  self%cicen = rhs%cicen
  self%hicen = rhs%hicen
  self%hsnon = rhs%hsnon

end subroutine soca_seaice_copy

! ------------------------------------------------------------------------------
subroutine soca_seaice_add(self, other)
  class(soca_seaice_type), intent(inout) :: self
  class(soca_seaice_type),    intent(in) :: other

  self%cicen = self%cicen + other%cicen
  self%hicen = self%hicen + other%hicen
  self%hsnon = self%hsnon + other%hsnon

end subroutine soca_seaice_add

! ------------------------------------------------------------------------------
subroutine soca_seaice_add_incr(self, incr)
  class(soca_seaice_type), intent(inout) :: self  !< State
  class(soca_seaice_type),    intent(in) :: incr  !< Increment

  real(kind=kind_real), allocatable :: aice_bkg(:,:), aice_ana(:,:)
  real(kind=kind_real), allocatable :: aice_incr(:,:), alpha(:,:)
  real(kind=kind_real) :: amin = 1e-6_kind_real
  real(kind=kind_real) :: amax = 10.0_kind_real
  real(kind=kind_real) :: min_ice = 1e-6_kind_real
  integer :: k
  integer :: isd, ied, jsd, jed, ncat

  ! Short cut to domain indices
  isd = self%isd
  ied = self%ied
  jsd = self%jsd
  jed = self%jed
  ncat = self%ncat

  ! Allocate memory for temporary arrays
  allocate(aice_bkg(isd:ied,jsd:jed))
  allocate(aice_ana(isd:ied,jsd:jed))
  allocate(aice_incr(isd:ied,jsd:jed))
  allocate(alpha(isd:ied,jsd:jed))

  ! Allocate ice and snow volume
  if (.not.(allocated(self%vicen))) then
     allocate(self%vicen(isd:ied,jsd:jed,1:ncat))
     allocate(self%vsnon(isd:ied,jsd:jed,1:ncat))
  end if

  ! Compute ice and snow volume
  self%vicen = self%hicen * self%cicen(:,:,2:)
  self%vsnon = self%hsnon * self%cicen(:,:,2:)

  ! Initialize aggregate fields
  aice_bkg  = sum(self%cicen(:,:,2:), dim=3)
  aice_incr = sum(incr%cicen(:,:,2:), dim=3)
  aice_ana  = aice_bkg + aice_incr

  ! Fix out of bound values in aggregate ice fraction analysis
  where (aice_ana < 0.0_kind_real)
     aice_ana = 0.0_kind_real
  end where
  where (aice_ana > 1.0_kind_real)
     aice_ana = 1.0_kind_real
  end where

  ! Compute background rescaling
  alpha = 1.0_kind_real
  where (aice_bkg > min_ice)
     alpha = aice_ana / aice_bkg
  end where

  ! Limit size of increment
  where ( alpha > amax )
     alpha = amax
  end where
  where ( alpha < amin )
     alpha = amin
  end where

  ! Add increment for ice and snow thickness
  ! TODO: check bounds ...
  self%hicen = self%hicen + incr%hicen
  self%hsnon = self%hsnon + incr%hsnon

  ! Update ice and snow volume
  self%vicen = self%hicen * self%cicen(:,:,2:)
  self%vsnon = self%hsnon * self%cicen(:,:,2:)

  ! "Add" fraction increment and update volumes accordingly
  do k = 1, self%ncat
     self%cicen(:,:,k+1) = alpha * self%cicen(:,:,k+1)
     self%vicen(:,:,k) = alpha * self%vicen(:,:,k)
     self%vsnon(:,:,k) = alpha * self%vsnon(:,:,k)
  end do

  ! Clean-up memory
  deallocate(aice_bkg, aice_ana, aice_incr, alpha)

end subroutine soca_seaice_add_incr

! ------------------------------------------------------------------------------
subroutine soca_seaice_schur(self, other)
  class(soca_seaice_type), intent(inout) :: self
  class(soca_seaice_type),    intent(in) :: other

  self%cicen = self%cicen * other%cicen
  self%hicen = self%hicen * other%hicen
  self%hsnon = self%hsnon * other%hsnon

end subroutine soca_seaice_schur

! ------------------------------------------------------------------------------
subroutine soca_seaice_sub(self, other)
  class(soca_seaice_type), intent(inout) :: self
  class(soca_seaice_type),    intent(in) :: other

  self%cicen = self%cicen - other%cicen
  self%hicen = self%hicen - other%hicen
  self%hsnon = self%hsnon - other%hsnon

end subroutine soca_seaice_sub

! ------------------------------------------------------------------------------
subroutine soca_seaice_mul(self, zz)
  class(soca_seaice_type), intent(inout) :: self
  real(kind=kind_real),       intent(in) :: zz

  self%cicen = zz * self%cicen
  self%hicen = zz * self%hicen
  self%hsnon = zz * self%hsnon

end subroutine soca_seaice_mul

! ------------------------------------------------------------------------------
subroutine soca_seaice_axpy(self, zz, other)
  class(soca_seaice_type), intent(inout) :: self
  real(kind=kind_real),       intent(in) :: zz
  class(soca_seaice_type),    intent(in) :: other

  self%cicen = self%cicen + zz * other%cicen
  self%hicen = self%hicen + zz * other%hicen
  self%hsnon = self%hsnon + zz * other%hsnon

end subroutine soca_seaice_axpy

! ------------------------------------------------------------------------------
subroutine soca_seaice_diff_incr(self, x1, x2)
  class(soca_seaice_type), intent(inout) :: self
  class(soca_seaice_type),    intent(in) :: x1
  class(soca_seaice_type),    intent(in) :: x2

  self%cicen = x1%cicen - x2%cicen
  self%hicen = x1%hicen - x2%hicen
  self%hsnon = x1%hsnon - x2%hsnon

end subroutine soca_seaice_diff_incr

! ------------------------------------------------------------------------------
subroutine soca_seaice_read_rst(self, f_conf, geom, fldnames)
  class(soca_seaice_type), intent(inout) :: self
  type(fckit_configuration),  intent(in) :: f_conf
  type(soca_geom),            intent(in) :: geom
  character(len=5),           intent(in) :: fldnames(:)

  integer, parameter :: max_string_length=800
  integer :: idr, i
  character(len=max_string_length) :: filename, basename
  character(len=4) :: seaice_model
  type(restart_file_type) :: restart
  character(len=:), allocatable :: str

  ! Check what model we are reading a file from ('sis2' or 'cice')
  seaice_model = 'sis2' ! Default model is sis2

  if ( f_conf%has("seaice_model") ) then
      call f_conf%get_or_die("seaice_model", str)
      seaice_model = str
  endif

  if ( f_conf%has("ice_filename") ) then
      call f_conf%get_or_die("basename", str)
      basename = str
      call f_conf%get_or_die("ice_filename", str)
     filename = trim(basename)//trim(str)
  else
     ! Set seaice state to 0 if no file provided
     call self%zeros()
     return
  end if

  select case(seaice_model)
  case('sis2')
     call fms_io_init()
     do i = 1, size(fldnames)
        select case(fldnames(i))
        case('cicen')
           idr = register_restart_field(restart, filename, 'part_size', &
                self%cicen(:,:,:), &
                domain=geom%Domain%mpp_domain)
        case('hicen')
           idr = register_restart_field(restart, filename, 'h_ice', &
                self%hicen(:,:,:), &
                domain=geom%Domain%mpp_domain)
        case('hsnon')
           idr = register_restart_field(restart, filename, 'h_snow', &
                self%hsnon(:,:,:), &
                domain=geom%Domain%mpp_domain)
        end select
     end do
     call restore_state(restart, directory='')
     call free_restart_type(restart)
     call fms_io_exit()
     ! Convert hicen & hsnon from [kg/m2] to [m]
     self%hicen(:,:,:) = self%hicen(:,:,:)/self%soca_rho_ice
     self%hsnon(:,:,:) = self%hsnon(:,:,:)/self%soca_rho_snow

  case('cice')
     call fms_io_init()
     do i = 1, size(fldnames)
        select case(fldnames(i))
        case('cicen')
           idr = register_restart_field(restart, filename, 'aicen', &
                self%cicen(:,:,2:), &
                domain=geom%Domain%mpp_domain)
        case('hicen')
           idr = register_restart_field(restart, filename, 'vicen', &
                self%hicen(:,:,:), &
                domain=geom%Domain%mpp_domain)
        case('hsnon')
           idr = register_restart_field(restart, filename, 'vsnon', &
                self%hsnon(:,:,:), &
                domain=geom%Domain%mpp_domain)
        end select

     end do
     call restore_state(restart, directory='')
     call free_restart_type(restart)
     call fms_io_exit()

     ! Add ocean fraction
     self%cicen(:,:,1) = 1.0_kind_real - sum(self%cicen(:,:,2:), dim=3)

     ! Convert to hicen and hsnon
     where(self%cicen(:,:,2:)>0.0_kind_real)
        self%hicen(:,:,:) = self%hicen(:,:,:)/self%cicen(:,:,2:)
        self%hsnon(:,:,:) = self%hsnon(:,:,:)/self%cicen(:,:,2:)
     end where

  case default
     call abor1_ftn("soca_seaice_mod: Reading for seaice model "//trim(seaice_model)//" not implemented")
  end select

  if (allocated(str)) deallocate(str)

end subroutine soca_seaice_read_rst

! ------------------------------------------------------------------------------
subroutine soca_seaice_read_diag(self, f_conf, geom, fldnames)
  class(soca_seaice_type), intent(inout) :: self
  type(fckit_configuration),  intent(in) :: f_conf
  type(soca_geom),            intent(in) :: geom
  character(len=5),           intent(in) :: fldnames(:)

  integer, parameter :: max_string_length=800
  integer :: i
  character(len=max_string_length) :: filename
  character(len=:), allocatable :: str

  if ( f_conf%has("filename") ) then
      call f_conf%get_or_die("filename", str)
      filename = str
  else
     call self%zeros()
     return
  end if

  call fms_io_init()
  do i = 1, size(fldnames)
     select case(fldnames(i))
     case('cicen')
        call read_data(filename,"cicen", &
                       self%cicen(:,:,:), &
                       domain=geom%Domain%mpp_domain)
     case('hicen')
        call read_data(filename,"hicen", &
                       self%hicen(:,:,:), &
                       domain=geom%Domain%mpp_domain)
     case('hsnon')
        call read_data(filename,"hsnon", &
                       self%hsnon(:,:,:), &
                       domain=geom%Domain%mpp_domain)
     end select
  end do
  call fms_io_exit()

end subroutine soca_seaice_read_diag

! ------------------------------------------------------------------------------
subroutine soca_seaice_write_rst(self, f_conf, geom, vdate)
  class(soca_seaice_type), intent(inout) :: self
  type(fckit_configuration),  intent(in) :: f_conf
  type(soca_geom),            intent(in) :: geom
  type(datetime),          intent(inout) :: vdate

  integer, parameter :: max_string_length=800
  integer :: idr
  character(len=max_string_length) :: filename
  character(len=4) :: seaice_model
  type(restart_file_type) :: restart
  real(kind=kind_real), allocatable :: vicen(:,:,:), vsnon(:,:,:)
  real(kind=kind_real), allocatable :: aice(:,:), hice(:,:), hsno(:,:) ! Aggregates
  integer :: isd, ied, jsd, jed
  character(len=:), allocatable :: str

  ! Check what model we are reading a file from ('sis2' or 'cice')
  seaice_model = 'sis2' ! Default model is sis2
  if ( f_conf%has("seaice_model") ) then
      call f_conf%get_or_die("seaice_model", str)
      seaice_model = str
  endif

  ! Register sea-ice fields
  call fms_io_init()

  ! Allocate and compute aggregate variables
  isd = self%isd ; ied = self%ied
  jsd = self%jsd ; jed = self%jed
  allocate(aice(isd:ied,jsd:jed))
  allocate(hice(isd:ied,jsd:jed))
  allocate(hsno(isd:ied,jsd:jed))
  aice(:,:) = sum(self%cicen(:,:,2:), dim=3)
  hice(:,:) = sum(self%hicen(:,:,:), dim=3)
  hsno(:,:) = sum(self%hsnon(:,:,:), dim=3)

  select case(seaice_model)
  case('sis2')
     ! Generate file names
     filename = soca_genfilename(f_conf,max_string_length,vdate,"ice")

     ! Register sis2 variables
     idr = register_restart_field(restart, filename, 'part_size', self%cicen, &
          domain=geom%Domain%mpp_domain)
     ! TODO (Guillaume): Rescale with density
     idr = register_restart_field(restart, filename, 'h_ice', self%hicen, &
          domain=geom%Domain%mpp_domain)
     idr = register_restart_field(restart, filename, 'h_snow', self%hsnon, &
          domain=geom%Domain%mpp_domain)

  case('cice')
     ! Get ice volume
     allocate(vicen(isd:ied,jsd:jed,geom%ice_column%ncat))
     allocate(vsnon(isd:ied,jsd:jed,geom%ice_column%ncat))
     if (.not.allocated(self%vicen))then
        vicen(:,:,:) = self%cicen(:,:,2:)*self%hicen(:,:,:)
        vsnon(:,:,:) = self%cicen(:,:,2:)*self%hsnon(:,:,:)
     else
        vicen(:,:,:) = self%vicen(:,:,:)
        vsnon(:,:,:) = self%vsnon(:,:,:)
     end if

     ! Generate file names
     filename = soca_genfilename(f_conf,max_string_length,vdate,"cice")

     ! Register cice variables
     idr = register_restart_field(restart, filename, 'aicen', self%cicen(:,:,2:), &
          domain=geom%Domain%mpp_domain)
     idr = register_restart_field(restart, filename, 'vicen', vicen, &
          domain=geom%Domain%mpp_domain)
     idr = register_restart_field(restart, filename, 'vsnon', vsnon, &
          domain=geom%Domain%mpp_domain)
  end select

  ! Register aggregate variables
  idr = register_restart_field(restart, filename, 'aice', aice, &
       domain=geom%Domain%mpp_domain)
  idr = register_restart_field(restart, filename, 'hice', hice, &
       domain=geom%Domain%mpp_domain)
  idr = register_restart_field(restart, filename, 'hsno', hsno, &
       domain=geom%Domain%mpp_domain)

  ! Write restart to disk
  call save_restart(restart, directory='')
  call free_restart_type(restart)
  call fms_io_exit()

  ! Clean-up
  if (allocated(vicen)) deallocate(vicen)
  if (allocated(vsnon)) deallocate(vsnon)
  deallocate(aice, hice, hsno)

end subroutine soca_seaice_write_rst

end module soca_seaice_mod
