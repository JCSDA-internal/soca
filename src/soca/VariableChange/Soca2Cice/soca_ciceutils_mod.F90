! (C) Copyright 2022-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_ciceutils_mod

use kinds, only: kind_real
use netcdf
use fckit_mpi_module, only: fckit_mpi_comm
use icepack_itd
use soca_utils, only: nc_check
use soca_geom_mod, only: soca_geom

implicit none
private

type, public :: agg_cice_state
   real(kind=kind_real),  allocatable :: aice(:)
   real(kind=kind_real),  allocatable :: lon(:)
   real(kind=kind_real),  allocatable :: lat(:)
   integer,  allocatable :: ij(:,:)
   integer :: n_src
end type agg_cice_state

type, public :: cice_state
   integer :: ncat, ni, nj, ice_lev=7, sno_lev=1
   integer :: isc, iec, jsc, jec
   character(len=:), allocatable :: rst_filename
   character(len=:), allocatable :: rst_out_filename
   real(kind=kind_real),  allocatable :: aicen(:,:,:)
   real(kind=kind_real),  allocatable :: vicen(:,:,:)
   real(kind=kind_real),  allocatable :: vsnon(:,:,:)
   real(kind=kind_real),  allocatable :: apnd(:,:,:)
   real(kind=kind_real),  allocatable :: hpnd(:,:,:)
   real(kind=kind_real),  allocatable :: ipnd(:,:,:)
   real(kind=kind_real),  allocatable :: tsfcn(:,:,:)
   real(kind=kind_real),  allocatable :: qice(:,:,:,:)
   real(kind=kind_real),  allocatable :: qsno(:,:,:,:)
   real(kind=kind_real),  allocatable :: sice(:,:,:,:)
   real(kind=kind_real),  allocatable :: iceumask(:,:)
   real(kind=kind_real),  allocatable :: aice(:,:)
   type(agg_cice_state):: agg
contains
  procedure :: init => soca_ciceutils_init
  procedure :: alloc => soca_ciceutils_alloc
  procedure :: read => soca_ciceutils_read
  procedure :: write => soca_ciceutils_write
  procedure :: broadcast => soca_ciceutils_broadcast
  procedure :: finialize => soca_ciceutils_finalize
end type cice_state

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
subroutine soca_ciceutils_init(self, geom, rst_filename, rst_out_filename)
  class(cice_state), intent(inout) :: self
  type(soca_geom), target, intent(in)  :: geom
  character(len=*),     intent(in) :: rst_filename
  character(len=*),     intent(in) :: rst_out_filename

  integer(kind=4) :: ncid, dimid, varid

  ! Save io info
  self%rst_filename = rst_filename
  self%rst_out_filename = rst_out_filename

  ! Get indices bounds for local arrays
  self%isc = geom%isc
  self%iec = geom%iec
  self%jsc = geom%jsc
  self%jec = geom%jec

  ! Get global seaice grid info
  call nc_check(nf90_open(self%rst_filename, nf90_nowrite, ncid))
  call nc_check(nf90_inq_dimid(ncid, 'ncat', dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len = self%ncat))
  call nc_check(nf90_inq_dimid(ncid, 'ni', dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len = self%ni))
  call nc_check(nf90_inq_dimid(ncid, 'nj', dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len = self%nj))
  call nc_check(nf90_close(ncid))

  ! Allocate seaice fields
  call self%alloc()

end subroutine soca_ciceutils_init

! ------------------------------------------------------------------------------
subroutine soca_ciceutils_broadcast(self, f_comm, root)
  class(cice_state), intent(inout) :: self
  type(fckit_mpi_comm), intent(in) :: f_comm
  integer,              intent(in) :: root

  ! TODO (G): Don't do that!
  ! broadcast size
  call f_comm%broadcast(self%ni, root)
  call f_comm%broadcast(self%nj, root)
  call f_comm%barrier()

  ! allocate cice fields
  !if (f_comm%rank().ne.root) call self%alloc()

  ! broadcast cice variables
  !call f_comm%broadcast(self%aice, root)
  !call f_comm%broadcast(self%aicen, root)
  !call f_comm%broadcast(self%vicen, root)
  !call f_comm%broadcast(self%vsnon, root)
  !call f_comm%barrier()

end subroutine soca_ciceutils_broadcast

! ------------------------------------------------------------------------------
subroutine soca_ciceutils_alloc(self)
  class(cice_state), intent(inout) :: self

  integer :: ni, nj, ncat, ice_lev, sno_lev
  integer :: isc, iec, jsc, jec

  ni = self%ni
  nj = self%nj
  ncat = self%ncat
  ice_lev = self%ice_lev
  sno_lev = self%sno_lev
  isc = self%isc
  iec = self%iec
  jsc = self%jsc
  jec = self%jec

  allocate(self%iceumask(isc:iec, jsc:jec));      self%iceumask = 0.0_kind_real

  allocate(self%aice(isc:iec, jsc:jec));          self%aice = 0.0_kind_real

  allocate(self%aicen(isc:iec, jsc:jec, ncat));   self%aicen = 0.0_kind_real
  allocate(self%vicen(isc:iec, jsc:jec, ncat));   self%vicen = 0.0_kind_real
  allocate(self%vsnon(isc:iec, jsc:jec, ncat));   self%vsnon = 0.0_kind_real

  allocate(self%apnd(isc:iec, jsc:jec, ncat));    self%apnd = 0.0_kind_real
  allocate(self%hpnd(isc:iec, jsc:jec, ncat));    self%hpnd = 0.0_kind_real
  allocate(self%ipnd(isc:iec, jsc:jec, ncat));    self%ipnd = 0.0_kind_real

  allocate(self%tsfcn(isc:iec, jsc:jec, ncat));           self%tsfcn = 0.0_kind_real
  allocate(self%qsno(isc:iec, jsc:jec, ncat, sno_lev));   self%qsno = 0.0_kind_real
  allocate(self%qice(isc:iec, jsc:jec, ncat, ice_lev));   self%qice = 0.0_kind_real
  allocate(self%sice(isc:iec, jsc:jec, ncat, ice_lev));   self%sice = 0.0_kind_real

end subroutine soca_ciceutils_alloc

! ------------------------------------------------------------------------------
subroutine soca_ciceutils_read(self, geom)
  class(cice_state), intent(inout) :: self
  type(soca_geom), target, intent(in)  :: geom

  integer(kind=4) :: ncid, dimid, varid, i, j, cnt, n_src
  real(kind=kind_real) :: aice0

  call nc_check(nf90_open(self%rst_filename, nf90_nowrite, ncid))

  ! ice mask
  call getvar2d('iceumask', self%iceumask, self%ni, self%nj, geom, ncid)

  ! dynamic variables
  call getvar3d('aicen', self%aicen, self%ni, self%nj, self%ncat, geom, ncid)
  call getvar3d('vicen', self%vicen, self%ni, self%nj, self%ncat, geom, ncid)
  call getvar3d('vsnon', self%vsnon, self%ni, self%nj, self%ncat, geom, ncid)

  ! pond stuff
  call getvar3d('apnd', self%apnd, self%ni, self%nj, self%ncat, geom, ncid)
  call getvar3d('hpnd', self%hpnd, self%ni, self%nj, self%ncat, geom, ncid)
  call getvar3d('ipnd', self%ipnd, self%ni, self%nj, self%ncat, geom, ncid)

  ! thermo variables
  call getvar3d('Tsfcn', self%tsfcn, self%ni, self%nj, self%ncat, geom, ncid)
  call getvar4d('qsno', self%qsno, self%ni, self%nj, self%ncat, self%sno_lev, geom, ncid)
  call getvar4d('qice', self%qice, self%ni, self%nj, self%ncat, self%ice_lev, geom, ncid)
  call getvar4d('sice', self%sice, self%ni, self%nj, self%ncat, self%ice_lev, geom, ncid)

  call nc_check(nf90_close(ncid))

  ! aggregate categories
  cnt = 1
  n_src = sum(self%iceumask)
  self%agg%n_src = n_src
  allocate(self%agg%lon(n_src), self%agg%lat(n_src), self%agg%aice(n_src))
  allocate(self%agg%ij(2,n_src))

  do j = geom%jsc, geom%jec
     do i = geom%isc, geom%iec
        if (self%iceumask(i,j).eq.1) then
           call aggregate_area (self%ncat, &
                                self%aicen(i,j,:), &
                                self%aice(i,j), aice0)
           self%agg%lon(cnt) = geom%lon(i,j)
           self%agg%lat(cnt) = geom%lat(i,j)
           self%agg%aice(cnt) = self%aice(i,j)
           self%agg%ij(1,cnt) = i
           self%agg%ij(2,cnt) = j
           cnt = cnt + 1
        end if
     end do
  end do

end subroutine soca_ciceutils_read


! ------------------------------------------------------------------------------
subroutine soca_ciceutils_write(self)
  class(cice_state), intent(inout) :: self

  integer(kind=4) :: ncid, dimid, ndims, varid, ni, nj, ncat, l
  character(len=12) :: cicevarname, varname
  character(len=1) :: strnum

  print *,self%rst_out_filename
  call nc_check( nf90_open(self%rst_out_filename, nf90_write, ncid) )
  call nc_check(nf90_inq_dimid(ncid, 'ncat', dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len = ncat))
  call nc_check(nf90_inq_dimid(ncid, 'ni', dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len = ni))
  call nc_check(nf90_inq_dimid(ncid, 'nj', dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len = nj))

  call nc_check(nf90_inq_varid(ncid,"aicen",varid))
  call nc_check(nf90_put_var(ncid,varid,self%aicen))

  call nc_check(nf90_inq_varid(ncid,"vicen",varid))
  call nc_check(nf90_put_var(ncid,varid,self%vicen))

  call nc_check(nf90_inq_varid(ncid,"vsnon",varid))
  call nc_check(nf90_put_var(ncid,varid,self%vsnon))

  call nc_check(nf90_inq_varid(ncid,"Tsfcn",varid))
  call nc_check(nf90_put_var(ncid,varid,self%tsfcn))

  call nc_check(nf90_inq_varid(ncid,"iceumask",varid))
  call nc_check(nf90_put_var(ncid,varid,self%iceumask))

  do l = 1, self%ice_lev
     write(strnum,'(i1)') l
     cicevarname = "qice00"//strnum
     call nc_check(nf90_inq_varid(ncid,cicevarname,varid))
     call nc_check(nf90_put_var(ncid,varid,self%qice(:,:,:,l)))

     write(strnum,'(i1)') l
     cicevarname = "sice00"//strnum
     call nc_check(nf90_inq_varid(ncid,cicevarname,varid))
     call nc_check(nf90_put_var(ncid,varid,self%sice(:,:,:,l)))
  end do

  do l = 1, self%sno_lev
     write(strnum,'(i1)') l
     cicevarname = "qsno00"//strnum
     call nc_check(nf90_inq_varid(ncid,cicevarname,varid))
     call nc_check(nf90_put_var(ncid,varid,self%qsno(:,:,:,l)))
  end do
  call nc_check(nf90_close(ncid))

end subroutine soca_ciceutils_write

! ------------------------------------------------------------------------------
subroutine soca_ciceutils_finalize(self)
  class(cice_state), intent(inout) :: self

  integer :: ni, nj, ncat, ice_lev, sno_lev

  ni = 0
  nj = 0
  ncat = 0
  ice_lev = 0
  sno_lev = 0

  deallocate(self%iceumask)

  deallocate(self%aicen)
  deallocate(self%vicen)
  deallocate(self%vsnon)

  deallocate(self%apnd) ! melt pond area fraction
  deallocate(self%hpnd) ! melt pond depth
  deallocate(self%ipnd) ! melt pond refrozen lid thickness

  deallocate(self%tsfcn)
  deallocate(self%qice)
  deallocate(self%sice)

end subroutine soca_ciceutils_finalize

! ------------------------------------------------------------------------------

subroutine getvar2d(varname, var, ni, nj, geom, ncid)
  character(len=*),                     intent(in) :: varname
  integer,                              intent(in) :: ni, nj, ncid
  type(soca_geom),             target, intent(in)  :: geom
  real(kind=kind_real), allocatable, intent(inout) :: var(:,:)

  integer(kind=4) :: varid
  real(kind=kind_real), allocatable :: tmpvar(:,:)

  if (.not. allocated(var)) allocate(var(geom%isc:geom%iec, geom%jsc:geom%jec))

  allocate(tmpvar(ni, nj))
  call nc_check(nf90_inq_varid(ncid,trim(varname),varid))
  call nc_check(nf90_get_var(ncid,varid,tmpvar))
  var(geom%isc:geom%iec, geom%jsc:geom%jec) = tmpvar(geom%isc:geom%iec, geom%jsc:geom%jec)

end subroutine getvar2d

subroutine getvar3d(varname, var, ni, nj, ncat, geom, ncid)
  character(len=*),                     intent(in) :: varname
  integer,                              intent(in) :: ni, nj, ncat, ncid
  type(soca_geom),             target, intent(in)  :: geom
  real(kind=kind_real), allocatable, intent(inout) :: var(:,:,:)

  integer(kind=4) :: varid
  real(kind=kind_real), allocatable :: tmpvar(:,:,:)

  if (.not. allocated(var)) allocate(var(geom%isc:geom%iec, geom%jsc:geom%jec, ncat))

  allocate(tmpvar(ni, nj, ncat))
  call nc_check(nf90_inq_varid(ncid,trim(varname),varid))
  call nc_check(nf90_get_var(ncid,varid,tmpvar))
  var(geom%isc:geom%iec, geom%jsc:geom%jec, :) = tmpvar(geom%isc:geom%iec, geom%jsc:geom%jec, :)

end subroutine getvar3d

subroutine getvar4d(varname, var, ni, nj, ncat, nl, geom, ncid)
  character(len=*),         intent(in) :: varname
  integer,                              intent(in) :: ni, nj, ncat, nl, ncid
  type(soca_geom),             target, intent(in)  :: geom
  real(kind=kind_real), allocatable, intent(inout) :: var(:,:,:,:)

  integer(kind=4) :: varid, l
  character(len=12) :: cicevarname
  character(len=1) :: strnum
  real(kind=kind_real), allocatable :: var3d(:,:,:)

  if (.not. allocated(var)) allocate(var(ni, nj, ncat, nl))
  do l = 1, nl
     write(strnum,'(i1)') l
     cicevarname = trim(varname)//"00"//strnum
     call getvar3d(cicevarname, var3d, ni, nj, ncat, geom, ncid)
     var(:,:,:,l) = var3d(:,:,:)
  end do

end subroutine getvar4d

end module soca_ciceutils_mod
