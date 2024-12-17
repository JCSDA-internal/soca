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
use mpp_domains_mod, only : mpp_update_domains
use mpp_mod, only : mpp_gather, mpp_root_pe, mpp_pe

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
   integer :: ncat, ni, nj, ice_lev, sno_lev
   integer :: isc, iec, jsc, jec   ! compute domain
   integer :: isd, ied, jsd, jed   ! data domain
   logical :: initialized = .false.
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
  procedure :: copydata => soca_ciceutils_copydata
  procedure :: read => soca_ciceutils_read
  procedure :: write => soca_ciceutils_write
  procedure :: gather => soca_ciceutils_gather
  procedure :: finialize => soca_ciceutils_finalize
end type cice_state

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
subroutine soca_ciceutils_init(self, geom, rst_filename, rst_out_filename, ice_lev, sno_lev, global)
  class(cice_state),    intent(inout) :: self
  type(soca_geom), target, intent(in) :: geom
  character(len=*),        intent(in) :: rst_filename
  character(len=*),        intent(in) :: rst_out_filename
  integer,                 intent(in) :: ice_lev, sno_lev
  logical,    optional, intent(inout) :: global
  integer(kind=4) :: ncid, dimid, varid
  logical :: isglobal

  ! Check if initializing global or local
  isglobal = .false.
  if (present(global)) isglobal = global

  ! Save io info
  self%rst_filename = rst_filename
  self%rst_out_filename = rst_out_filename

  ! Indices for the compute domain
  self%isc = geom%isc
  self%iec = geom%iec
  self%jsc = geom%jsc
  self%jec = geom%jec

  ! Indices for the data domain
  self%isd = geom%isd
  self%ied = geom%ied
  self%jsd = geom%jsd
  self%jed = geom%jed

  ! Get global seaice grid info
  call nc_check(nf90_open(self%rst_filename, nf90_nowrite, ncid))
  call nc_check(nf90_inq_dimid(ncid, 'ncat', dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len = self%ncat))
  call nc_check(nf90_inq_dimid(ncid, 'ni', dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len = self%ni))
  call nc_check(nf90_inq_dimid(ncid, 'nj', dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len = self%nj))
  call nc_check(nf90_close(ncid))

  self%ice_lev = ice_lev
  self%sno_lev = sno_lev

  ! Allocate seaice fields
  if (isglobal) then
     call self%alloc(1, self%ni, 1, self%nj)
  else
     call self%alloc(self%isd, self%ied, self%jsd, self%jed)
  end if

  self%initialized = .true.

end subroutine soca_ciceutils_init

! ------------------------------------------------------------------------------
subroutine soca_ciceutils_gather(self, glb, geom)
  class(cice_state),    intent(inout) :: self
  type(cice_state),    intent(inout) :: glb
  type(soca_geom), target, intent(in) :: geom

  integer :: n, i, pe, l
  integer, allocatable :: pelist(:)
  logical :: is_root_pe, isglobal = .true.
  real(kind=kind_real),  allocatable :: var2d(:,:)

  ! Prepare mpi info for mpp calls
  pe = geom%f_comm%rank()
  is_root_pe = .false.
  if (pe == mpp_root_pe()) is_root_pe = .true.
  allocate(pelist(geom%f_comm%size()))
  pelist = (/(i,i=0,geom%f_comm%size()-1)/)

  ! allocate the global cice_state data structure
  if (geom%f_comm%rank() == mpp_root_pe()) then
     call glb%init(geom, self%rst_filename, self%rst_out_filename, self%ice_lev, self%sno_lev, global=isglobal)
     allocate(var2d(self%ni, self%nj))
  end if

  ! gather 2D arrays
  call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
       self%iceumask(self%isc:self%iec,self%jsc:self%jec), &
       var2d, is_root_pe)
  if ( pe == mpp_root_pe()) glb%iceumask(1:self%ni,1:self%nj) = var2d
  call geom%f_comm%barrier()

  ! gather 2D arrays along categories, snow and ice levels
  do n = 1, self%ncat
     ! aicen
     call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
          self%aicen(self%isc:self%iec,self%jsc:self%jec,n), &
          var2d, is_root_pe)
     if ( pe == mpp_root_pe()) glb%aicen(1:self%ni,1:self%nj,n) = var2d*glb%iceumask

     ! vicen
     call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
          self%vicen(self%isc:self%iec,self%jsc:self%jec,n), &
          var2d, is_root_pe)
     if ( pe == mpp_root_pe()) glb%vicen(1:self%ni,1:self%nj,n) = var2d*glb%iceumask

     ! vsnon
     call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
          self%vsnon(self%isc:self%iec,self%jsc:self%jec,n), &
          var2d, is_root_pe)
     if ( pe == mpp_root_pe()) glb%vsnon(1:self%ni,1:self%nj,n) = var2d*glb%iceumask

     ! Tsfcn
     call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
          self%tsfcn(self%isc:self%iec,self%jsc:self%jec,n), &
          var2d, is_root_pe)
     if ( pe == mpp_root_pe()) glb%tsfcn(1:self%ni,1:self%nj,n) = var2d*glb%iceumask

     do l = 1, self%ice_lev
        ! Qice
        call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
             self%qice(self%isc:self%iec,self%jsc:self%jec,n,l), &
             var2d, is_root_pe)
        if ( pe == mpp_root_pe()) glb%qice(1:self%ni,1:self%nj,n,l) = var2d*glb%iceumask

        ! Sice
        call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
             self%sice(self%isc:self%iec,self%jsc:self%jec,n,l), &
             var2d, is_root_pe)
        if ( pe == mpp_root_pe()) glb%sice(1:self%ni,1:self%nj,n,l) = var2d*glb%iceumask
     end do

     do l = 1, self%sno_lev
        ! Qsno
        call mpp_gather(self%isc, self%iec, self%jsc, self%jec, pelist, &
             self%qsno(self%isc:self%iec,self%jsc:self%jec,n,l), &
             var2d, is_root_pe)
        if ( pe == mpp_root_pe()) glb%qsno(1:self%ni,1:self%nj,n,l) = var2d*glb%iceumask
     end do
  end do

end subroutine soca_ciceutils_gather

! ------------------------------------------------------------------------------
subroutine soca_ciceutils_alloc(self, isd, ied, jsd, jed)
  class(cice_state), intent(inout) :: self
  integer,              intent(in) :: isd, ied, jsd, jed
  integer :: ncat, ice_lev, sno_lev

  ncat = self%ncat
  ice_lev = self%ice_lev
  sno_lev = self%sno_lev

  allocate(self%iceumask(isd:ied, jsd:jed));      self%iceumask = 0.0_kind_real

  allocate(self%aice(isd:ied, jsd:jed));          self%aice = 0.0_kind_real

  allocate(self%aicen(isd:ied, jsd:jed, ncat));   self%aicen = 0.0_kind_real
  allocate(self%vicen(isd:ied, jsd:jed, ncat));   self%vicen = 0.0_kind_real
  allocate(self%vsnon(isd:ied, jsd:jed, ncat));   self%vsnon = 0.0_kind_real

  allocate(self%apnd(isd:ied, jsd:jed, ncat));    self%apnd = 0.0_kind_real
  allocate(self%hpnd(isd:ied, jsd:jed, ncat));    self%hpnd = 0.0_kind_real
  allocate(self%ipnd(isd:ied, jsd:jed, ncat));    self%ipnd = 0.0_kind_real

  allocate(self%tsfcn(isd:ied, jsd:jed, ncat));           self%tsfcn = 0.0_kind_real
  allocate(self%qsno(isd:ied, jsd:jed, ncat, sno_lev));   self%qsno = 0.0_kind_real
  allocate(self%qice(isd:ied, jsd:jed, ncat, ice_lev));   self%qice = 0.0_kind_real
  allocate(self%sice(isd:ied, jsd:jed, ncat, ice_lev));   self%sice = 0.0_kind_real

end subroutine soca_ciceutils_alloc

! ------------------------------------------------------------------------------
subroutine soca_ciceutils_copydata(self, other)
  class(cice_state), intent(inout) :: self
  class(cice_state),    intent(in) :: other
  integer :: isd, ied, jsd, jed
  integer :: ncat, ice_lev, sno_lev

  self%ncat = other%ncat; ncat = other%ncat
  self%ice_lev = other%ice_lev; ice_lev = other%ice_lev
  self%sno_lev = other%sno_lev; sno_lev = other%sno_lev
  self%ni = other%ni
  self%nj = other%nj
  self%isc = other%isc
  self%iec = other%iec
  self%jsc = other%jsc
  self%jec = other%jec
  self%isd = other%isd; isd = other%isd
  self%ied = other%ied; ied = other%ied
  self%jsd = other%jsd; jsd = other%jsd
  self%jed = other%jed; jed = other%jed

  allocate(self%iceumask(isd:ied, jsd:jed));      self%iceumask = other%iceumask

  allocate(self%aice(isd:ied, jsd:jed));          self%aice = other%aice

  allocate(self%aicen(isd:ied, jsd:jed, ncat));   self%aicen = other%aicen
  allocate(self%vicen(isd:ied, jsd:jed, ncat));   self%vicen = other%vicen
  allocate(self%vsnon(isd:ied, jsd:jed, ncat));   self%vsnon = other%vsnon

  allocate(self%apnd(isd:ied, jsd:jed, ncat));    self%apnd = other%apnd
  allocate(self%hpnd(isd:ied, jsd:jed, ncat));    self%hpnd = other%hpnd
  allocate(self%ipnd(isd:ied, jsd:jed, ncat));    self%ipnd = other%ipnd

  allocate(self%tsfcn(isd:ied, jsd:jed, ncat));           self%tsfcn = other%tsfcn
  allocate(self%qsno(isd:ied, jsd:jed, ncat, sno_lev));   self%qsno = other%qsno
  allocate(self%qice(isd:ied, jsd:jed, ncat, ice_lev));   self%qice = other%qice
  allocate(self%sice(isd:ied, jsd:jed, ncat, ice_lev));   self%sice = other%sice

end subroutine soca_ciceutils_copydata

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

  do j = geom%jsd, geom%jed
     do i = geom%isd, geom%ied
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
subroutine soca_ciceutils_write(self, geom)
  class(cice_state),     intent(inout) :: self
  type(soca_geom), target, intent(in)  :: geom

  integer(kind=4) :: ncid, dimid, ndims, varid, ni, nj, ncat, l
  character(len=12) :: cicevarname, varname
  character(len=1) :: strnum
  type(cice_state) :: glb

  ! MPI gather of all the sea-ice variables in the cice_state data-structure
  call self%gather(glb, geom)

  ! Over-write the CICE restart with what's in glb
  if (geom%f_comm%rank() == mpp_root_pe()) then
     call nc_check( nf90_open(self%rst_out_filename, nf90_write, ncid) )
     call nc_check(nf90_inq_dimid(ncid, 'ncat', dimid))
     call nc_check(nf90_inquire_dimension(ncid, dimid, len = ncat))
     call nc_check(nf90_inq_dimid(ncid, 'ni', dimid))
     call nc_check(nf90_inquire_dimension(ncid, dimid, len = ni))
     call nc_check(nf90_inq_dimid(ncid, 'nj', dimid))
     call nc_check(nf90_inquire_dimension(ncid, dimid, len = nj))

     call nc_check(nf90_inq_varid(ncid,"aicen",varid))
     call nc_check(nf90_put_var(ncid,varid,glb%aicen))

     call nc_check(nf90_inq_varid(ncid,"vicen",varid))
     call nc_check(nf90_put_var(ncid,varid,glb%vicen))

     call nc_check(nf90_inq_varid(ncid,"vsnon",varid))
     call nc_check(nf90_put_var(ncid,varid,glb%vsnon))

     call nc_check(nf90_inq_varid(ncid,"Tsfcn",varid))
     call nc_check(nf90_put_var(ncid,varid,glb%tsfcn))

     call nc_check(nf90_inq_varid(ncid,"iceumask",varid))
     call nc_check(nf90_put_var(ncid,varid,glb%iceumask))

     do l = 1, self%ice_lev
        write(strnum,'(i1)') l
        cicevarname = "qice00"//strnum
        call nc_check(nf90_inq_varid(ncid,cicevarname,varid))
        call nc_check(nf90_put_var(ncid,varid,glb%qice(:,:,:,l)))

        write(strnum,'(i1)') l
        cicevarname = "sice00"//strnum
        call nc_check(nf90_inq_varid(ncid,cicevarname,varid))
        call nc_check(nf90_put_var(ncid,varid,glb%sice(:,:,:,l)))
     end do

     do l = 1, self%sno_lev
        write(strnum,'(i1)') l
        cicevarname = "qsno00"//strnum
        call nc_check(nf90_inq_varid(ncid,cicevarname,varid))
        call nc_check(nf90_put_var(ncid,varid,glb%qsno(:,:,:,l)))
     end do
     call nc_check(nf90_close(ncid))
  end if
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

  if (.not. allocated(var)) allocate(var(geom%isd:geom%ied, geom%jsd:geom%jed))

  allocate(tmpvar(ni, nj))
  call nc_check(nf90_inq_varid(ncid,trim(varname),varid))
  call nc_check(nf90_get_var(ncid,varid,tmpvar))
  var(geom%isc:geom%iec, geom%jsc:geom%jec) = tmpvar(geom%isc:geom%iec, geom%jsc:geom%jec)
  call geom%f_comm%barrier()
  call mpp_update_domains(var, geom%Domain%mpp_domain)
end subroutine getvar2d

subroutine getvar3d(varname, var, ni, nj, ncat, geom, ncid)
  character(len=*),                     intent(in) :: varname
  integer,                              intent(in) :: ni, nj, ncat, ncid
  type(soca_geom),             target, intent(in)  :: geom
  real(kind=kind_real), allocatable, intent(inout) :: var(:,:,:)

  integer(kind=4) :: varid, k
  real(kind=kind_real), allocatable :: tmpvar(:,:,:)

  if (.not. allocated(var)) allocate(var(geom%isd:geom%ied, geom%jsd:geom%jed, ncat))

  allocate(tmpvar(ni, nj, ncat))
  call nc_check(nf90_inq_varid(ncid,trim(varname),varid))
  call nc_check(nf90_get_var(ncid,varid,tmpvar))
  var(geom%isc:geom%iec, geom%jsc:geom%jec, :) = tmpvar(geom%isc:geom%iec, geom%jsc:geom%jec, :)
  call geom%f_comm%barrier()
  do k = 1, ncat
     call mpp_update_domains(var(:,:,k), geom%Domain%mpp_domain)
  end do
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

  if (.not. allocated(var)) allocate(var(geom%isd:geom%ied, geom%jsd:geom%jed, ncat, nl))
  do l = 1, nl
     write(strnum,'(i1)') l
     cicevarname = trim(varname)//"00"//strnum
     call getvar3d(cicevarname, var3d, ni, nj, ncat, geom, ncid)
     var(geom%isd:geom%ied, geom%jsd:geom%jed,1:ncat,l) = var3d(geom%isd:geom%ied, geom%jsd:geom%jed,1:ncat)
  end do

end subroutine getvar4d

end module soca_ciceutils_mod
