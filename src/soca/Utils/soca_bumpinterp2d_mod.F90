! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bumpinterp2d_mod

use fckit_mpi_module, only: fckit_mpi_comm
use fckit_log_module, only : fckit_log
use kinds, only: kind_real
use type_bump, only: bump_type
use type_obsop, only: obsop_type

implicit none

private
public :: soca_bumpinterp2d

type :: soca_bumpinterp2d
   type(obsop_type)                  :: obsop         !< bump interp object
   integer                           :: nobs          !< Number of values to interpolate
   real(kind=kind_real), allocatable :: lono(:)       !< Longitude of destination
   real(kind=kind_real), allocatable :: lato(:)       !< Latitude of destination
   logical                           :: initialized = .false. !< Initialization switch
 contains
   procedure :: initialize => interp_init
   procedure :: info => interp_info
   procedure :: apply => interp_apply
   procedure :: applyad => interpad_apply
   procedure :: finalize => interp_exit
end type soca_bumpinterp2d

! The bump grid here is stored separately from the interpolation weights
! (obsop inside soca_bumpinterp2d) because the unstructured grid only needs
! to be setup once
! TODO: these shouldn't be globals, put them somewhere correct
type(bump_type) :: bump
logical         :: bump_initialized = .false.

contains

!--------------------------------------------
subroutine interp_init(self, mod_lon, mod_lat, mod_mask, obs_lon, obs_lat, bumpid)
  ! Adapted from the fv3-jedi interface

  class(soca_bumpinterp2d), intent(inout) :: self
  real(kind=kind_real),      intent(in) :: mod_lon(:,:)
  real(kind=kind_real),      intent(in) :: mod_lat(:,:)
  real(kind=kind_real),      intent(in) :: mod_mask(:,:)
  real(kind=kind_real),      intent(in) :: obs_lon(:)
  real(kind=kind_real),      intent(in) :: obs_lat(:)
  integer,                   intent(in) :: bumpid
  !Locals
  integer :: ns, no, ni, nj
  real(kind=kind_real), allocatable :: area(:),vunit(:,:)
  real(kind=kind_real), allocatable :: tmp_lonmod(:), tmp_latmod(:)
  logical             , allocatable :: tmp_maskmod(:,:)

  character(len=5) :: cbumpcount
  character(len=23) :: bump_nam_prefix

  if (self%initialized) call interp_exit(self)

  ! Initialize the bump grid for the obs interpolation
  ! This only needs to be done once
  if(.not. bump_initialized) then
    ! Each bump%nam%prefix must be distinct, set bump id
    write(cbumpcount,"(I0.5)") bumpid
    bump_nam_prefix = 'soca_bump_interp_'//cbumpcount

    !Get the obs and state dimension
    ni = size(mod_lon, 1)
    nj = size(mod_lon, 2)
    ns = ni * nj

    ! Initialize bump parameters
    call bump%nam%init()
    bump%nam%prefix = bump_nam_prefix   ! Prefix for files output
    bump%nam%obsop_interp = 'bilin'     ! Interpolation type (bilinear)
    bump%nam%default_seed = .true.
    bump%nam%new_obsop = .true.
    bump%nam%verbosity = 'none'
    bump%nam%write_obsop = .false.
    bump%nam%sam_write = .false.
    bump%nam%strategy = 'specific_univariate'

    !Initialize geometry
    allocate(area(ns))
    allocate(vunit(ns,1))
    area = 1.0           ! Dummy area
    vunit = 1.0          ! Dummy vertical unit

    !Allocate temporary arrays
    allocate(tmp_lonmod(ns), tmp_latmod(ns), tmp_maskmod(ni, nj))
    tmp_lonmod(:) = reshape(mod_lon, (/ns/))
    tmp_latmod(:) = reshape(mod_lat, (/ns/))

    ! TODO: Fix interp with mask
  !!$   tmp_maskmod = .false.
  !!$    where(mod_mask==1.0)
  !!$       tmp_maskmod = .true.
  !!$    end where
    tmp_maskmod = reshape(tmp_maskmod, (/ns, 1/))
    tmp_maskmod = .true.

    ! Rotate longitudes
    where (tmp_lonmod < -180.0_kind_real)
       tmp_lonmod = tmp_lonmod + 360.0_kind_real
    end where

    !Initialize the bump geometry
    call bump%setup_online( ns, 1, 1, 1,&
         &tmp_lonmod, tmp_latmod, area, vunit, tmp_maskmod(:,1))

    !Release memory
    deallocate(area)
    deallocate(vunit)
    deallocate(tmp_lonmod, tmp_latmod, tmp_maskmod)

    bump_initialized = .true.
  end if

  ! initialize the interpolation weights
  no = size(obs_lon, 1)
  call self%obsop%from(no, obs_lon, obs_lat)
  call self%obsop%run_obsop(bump%mpl,bump%rng,bump%nam,bump%geom)

  self%initialized = .true.
  self%nobs = no

end subroutine interp_init

!--------------------------------------------
subroutine interp_apply(self, fld, obs)
  ! Forward interpolation: "fields to obs"
  ! obs = interp(fld)

  class(soca_bumpinterp2d), intent(inout) :: self
  real(kind=kind_real),     intent(in) :: fld(:,:)
  real(kind=kind_real),    intent(out) :: obs(:)

  ! Locals
  !real(kind=kind_real), allocatable :: tmp_fld(:,:)
  real(kind=kind_real), allocatable :: tmp_obs(:,:)
  integer :: ns

  !allocate(tmp_fld(size(fld,1),size(fld,2)))
  allocate(tmp_obs(size(obs,1),1))
  ns = size(fld,1)*size(fld,2)

  !tmp_fld = fld
  !tmp_fld = reshape(tmp_fld,(/ns, 1/))
  tmp_obs(:,1) = obs
  call self%obsop%apply(bump%mpl,bump%geom,fld,tmp_obs)
  obs = tmp_obs(:,1)

  deallocate(tmp_obs)

end subroutine interp_apply

!--------------------------------------------
subroutine interpad_apply(self, fld, obs)
  ! Backward interpolation: "obs to fields"
  ! fld = interpad(obs)

  class(soca_bumpinterp2d),    intent(inout) :: self
  real(kind=kind_real),     intent(inout) :: fld(:,:)
  real(kind=kind_real),        intent(in) :: obs(:)

  call self%obsop%apply_ad(bump%mpl,bump%geom,obs,fld)

end subroutine interpad_apply

!--------------------------------------------
subroutine interp_info(self)

  class(soca_bumpinterp2d), intent(in) :: self
  character(len=160) :: record

  write(record,*) "soca_bumpinterp2d%nobs: ",self%nobs
  call fckit_log%info(record)
  write(record,*) "soca_bumpinterp2d%initialized: ",self%initialized
  call fckit_log%info(record)

end subroutine interp_info

!--------------------------------------------
subroutine interp_exit(self)

  class(soca_bumpinterp2d), intent(out) :: self

  self%nobs = 0
  if (allocated(self%lono)) deallocate(self%lono)
  if (allocated(self%lato)) deallocate(self%lato)
  self%initialized = .false.

end subroutine interp_exit

end module soca_bumpinterp2d_mod
