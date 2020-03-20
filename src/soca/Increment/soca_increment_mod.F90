! (C) Copyright 2020-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_increment_mod

use soca_fields_mod
use soca_geom_iter_mod, only : soca_geom_iter
use kinds, only: kind_real
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use random_mod, only: normal_distribution
use unstructured_grid_mod, only: unstructured_grid, &
                                 allocate_unstructured_grid_coord, &
                                 allocate_unstructured_grid_field

implicit none
private

type, public, extends(soca_fields) :: soca_increment

contains
  ! get/set a single point
  procedure :: getpoint  => soca_increment_getpoint
  procedure :: setpoint  => soca_increment_setpoint

  ! unstructured grid I/O
  procedure :: from_ug   => soca_increment_from_ug
  procedure :: to_ug     => soca_increment_to_ug
  procedure :: ug_coord  => soca_increment_ug_coord

  ! misc
  procedure :: dirac     => soca_increment_dirac
  procedure :: random    => soca_increment_random
  procedure :: schur     => soca_increment_schur
end type


contains

! ------------------------------------------------------------------------------
!> initialize fields with random normal distribution
subroutine soca_increment_random(self)
  class(soca_increment), intent(inout) :: self
  integer, parameter :: rseed = 1 ! constant for reproducability of tests
    ! NOTE: random seeds are not quite working the way expected,
    !  it is only set the first time normal_distribution() is called with a seed
  integer :: z, i

  type(soca_field), pointer :: field

  ! set random values
  do i = 1, size(self%fields)
    field => self%fields(i)
    ! TODO remove this once increment / state are fully separated
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
  call self%update_halos()
end subroutine soca_increment_random


! ------------------------------------------------------------------------------
!> perform a shur product between two sets of fields
subroutine soca_increment_schur(self,rhs)
  class(soca_increment), intent(inout) :: self
  class(soca_increment),    intent(in) :: rhs
  integer :: i

  ! make sure fields are same shape
  call self%check_congruent(rhs)

  ! schur product
  do i=1,size(self%fields)
    self%fields(i)%val = self%fields(i)%val * rhs%fields(i)%val
  end do
end subroutine soca_increment_schur

! ------------------------------------------------------------------------------

subroutine soca_increment_getpoint(self, geoiter, values)
  class(soca_increment), intent(   in) :: self
  type(soca_geom_iter),  intent(   in) :: geoiter
  real(kind=kind_real),  intent(inout) :: values(:)

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
end subroutine soca_increment_getpoint

! ------------------------------------------------------------------------------

subroutine soca_increment_setpoint(self, geoiter, values)
  class(soca_increment), intent(inout) :: self
  type(soca_geom_iter),  intent(   in) :: geoiter
  real(kind=kind_real),  intent(   in) :: values(:)

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
end subroutine soca_increment_setpoint


  ! ------------------------------------------------------------------------------
! TODO, generalize by removing the hardcoded int=>field_name
subroutine soca_increment_dirac(self, f_conf)
  class(soca_increment),        intent(inout) :: self
  type(fckit_configuration), value, intent(in):: f_conf   !< Configuration

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
end subroutine soca_increment_dirac


! ------------------------------------------------------------------------------
! TODO remove hardcoded number of variables
subroutine ug_size(self, ug)
  class(soca_increment),      intent(in) :: self
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

subroutine soca_increment_ug_coord(self, ug)
  class(soca_increment),      intent(in) :: self
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

end subroutine soca_increment_ug_coord


! ------------------------------------------------------------------------------
! TODO generalize to use all vars
subroutine soca_increment_to_ug(self, ug, its)
  class(soca_increment),      intent(in) :: self
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

end subroutine soca_increment_to_ug

! ------------------------------------------------------------------------------
! Generalize variable names used
subroutine soca_increment_from_ug(self, ug, its)
  class(soca_increment), intent(inout) :: self
  type(unstructured_grid),  intent(in) :: ug
  integer,                  intent(in) :: its

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

end subroutine soca_increment_from_ug


end module soca_increment_mod