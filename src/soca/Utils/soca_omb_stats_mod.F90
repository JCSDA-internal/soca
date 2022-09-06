! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> surface background error used by soca_bkgerrgodas_mod
module soca_omb_stats_mod

use fckit_log_module, only: fckit_log
use fckit_mpi_module, only: fckit_mpi_comm
use kinds, only: kind_real
use netcdf
use soca_utils, only: nc_check, soca_remap_idw

implicit none
private


!> domain indices used by soca_omb_stats
type, public :: soca_domain_indices ! TODO: Move elsewhere!
   integer :: is, ie, js, je     ! Compute domain indices
   integer :: isl, iel, jsl, jel ! Local compute domain indices
end type soca_domain_indices


!> interpolate surface background error file to grid
!!
!! Used by soca_bkgerrgodas_mod::soca_bkgerrgodas_config
type, public :: soca_omb_stats
   integer                            :: nlocs
   real(kind=kind_real),  allocatable :: lon(:)
   real(kind=kind_real),  allocatable :: lat(:)
   real(kind=kind_real),  allocatable :: bgerr(:)
   real(kind=kind_real),  allocatable :: bgerr_model(:,:)
   type(soca_domain_indices)          :: domain
 contains

   !> \copybrief soca_omb_stats_init \see soca_omb_stats_init
   procedure :: init => soca_omb_stats_init

   !> \copybrief soca_omb_stats_bin \see soca_omb_stats_bin
   procedure :: bin => soca_omb_stats_bin

   !> \copybrief soca_omb_stats_exit \see soca_omb_stats_exit
   procedure :: exit => soca_omb_stats_exit
end type soca_omb_stats


contains

! ------------------------------------------------------------------------------
!> constructor
!!
!! \relates soca_omb_stats_mod::soca_omb_stats
subroutine soca_omb_stats_init(self, domain, filename)
  class(soca_omb_stats),           intent(inout) :: self
  type(soca_domain_indices),       intent(in) :: domain
  character(len=:), allocatable,   intent(in) :: filename

  integer(kind=4) :: ncid
  integer(kind=4) :: dimid
  integer(kind=4) :: varid
  type(fckit_mpi_comm) :: f_comm
  integer :: myrank, root=0, ret

  ! Setup Communicator
  f_comm = fckit_mpi_comm()
  myrank = f_comm%rank()

  if (myrank.eq.root) then

     call fckit_log%info("Reading file "  // trim(filename))

     call nc_check(nf90_open(filename, nf90_nowrite, ncid))

     ! Get the size of the horizontal grid
     call nc_check(nf90_inq_dimid(ncid, 'nlocs', dimid))
     call nc_check(nf90_inquire_dimension(ncid, dimid, len = self%nlocs))

     allocate(self%lon(self%nlocs), self%lat(self%nlocs), self%bgerr(self%nlocs))

     ! Get longitude
     call nc_check(nf90_inq_varid(ncid,'longitude',varid))
     call nc_check(nf90_get_var(ncid,varid,self%lon))

     ! Get latitude
     call nc_check(nf90_inq_varid(ncid,'latitude',varid))
     call nc_check(nf90_get_var(ncid,varid,self%lat))

     ! Get omb stats
     call nc_check(nf90_inq_varid(ncid,'sst_bgerr',varid))
     call nc_check(nf90_get_var(ncid,varid,self%bgerr))

     ! Close netcdf file
     call nc_check(nf90_close(ncid))
  end if

  ! Broadcast to all workers
  call f_comm%broadcast(self%nlocs, root)
  if (myrank.ne.root) then
     allocate(self%lon(self%nlocs), self%lat(self%nlocs), self%bgerr(self%nlocs))
  end if
  call f_comm%broadcast(self%lon, root)
  call f_comm%broadcast(self%lat, root)
  call f_comm%broadcast(self%bgerr, root)
  call f_comm%barrier()

  ! Rotate longitude
  where (self%lon>180.0_kind_real)
     self%lon=self%lon-360.0_kind_real
  end where

  ! Compute domain info
  self%domain = domain

end subroutine soca_omb_stats_init


! ------------------------------------------------------------------------------
!> remap background error to grid
!!
!! \relates soca_omb_stats_mod::soca_omb_stats
subroutine soca_omb_stats_bin(self, lon, lat)
  class(soca_omb_stats), intent(inout) :: self
  real(kind=kind_real),     intent(in) :: lon(:,:)
  real(kind=kind_real),     intent(in) :: lat(:,:)

  integer :: is, ie, js, je
  integer :: isl, iel, jsl, jel

  ! Short cuts to global indices
  is = self%domain%is
  ie = self%domain%ie
  js = self%domain%js
  je = self%domain%je

  ! Short cuts to local indices
  isl = self%domain%isl
  iel = self%domain%iel
  jsl = self%domain%jsl
  jel = self%domain%jel

  allocate(self%bgerr_model(is:ie,js:je))
  self%bgerr_model = 0.0_kind_real
  call soca_remap_idw(self%lon, self%lat, self%bgerr,&
                      lon(isl:iel,jsl:jel), lat(isl:iel,jsl:jel), &
                      self%bgerr_model(is:ie,js:je))

end subroutine soca_omb_stats_bin


! ------------------------------------------------------------------------------
!> Destructor
!!
!! \relates soca_omb_stats_mod::soca_omb_stats
subroutine soca_omb_stats_exit(self)
  class(soca_omb_stats), intent(inout) :: self

  deallocate(self%lon, self%lat, self%bgerr, self%bgerr_model)

end subroutine soca_omb_stats_exit

end module soca_omb_stats_mod
