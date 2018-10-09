!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bumpinterp2d_mod

  use kinds
  use type_bump, only: bump_type

  implicit none
  private

  type, public :: soca_bumpinterp2d
     type(bump_type)                   :: bump          !< bump interp object
     integer                           :: nobs          !< Number of values to interpolate
     real(kind=kind_real), allocatable :: lono(:)       !< Longitude of destination
     real(kind=kind_real), allocatable :: lato(:)       !< Latitude of destination 
     logical                           :: initialized   !< Initialization switch
   contains
     procedure :: initialize => interp_init
     procedure :: info => interp_info
     procedure :: apply => interp_apply     
     procedure :: applyad => interpad_apply
     procedure :: finalize => interp_exit
  end type soca_bumpinterp2d

contains

  !--------------------------------------------
  subroutine interp_init(self, mod_lon, mod_lat, mod_mask, obs_lon, obs_lat)
    ! Adapted from the fv3-jedi interface

    use fckit_mpi_module, only: fckit_mpi_comm
    use type_bump, only: bump_type
    use mpi
    
    implicit none

    class(soca_bumpinterp2d), intent(inout) :: self    
    real(kind=kind_real),      intent(in) :: mod_lon(:,:)
    real(kind=kind_real),      intent(in) :: mod_lat(:,:)
    real(kind=kind_real),      intent(in) :: mod_mask(:,:)    
    real(kind=kind_real),      intent(in) :: obs_lon(:)
    real(kind=kind_real),      intent(in) :: obs_lat(:)    
    
    !Locals
    integer :: ns, no, ni, nj
    real(kind=kind_real), allocatable :: area(:),vunit(:,:)
    real(kind=kind_real), allocatable :: tmp_lonmod(:), tmp_latmod(:)    
    logical             , allocatable :: tmp_maskmod(:,:)
    
    integer, save :: bumpcount = 0
    character(len=5) :: cbumpcount
    character(len=16) :: bump_nam_prefix

    type(fckit_mpi_comm) :: f_comm

    if (self%initialized) call interp_exit(self)
    
    f_comm = fckit_mpi_comm()

    ! Each bump%nam%prefix must be distinct
    ! -------------------------------------
    
    bumpcount = bumpcount + 1
    write(cbumpcount,"(I0.5)") bumpcount
    bump_nam_prefix = 'soca_bump_data_'//cbumpcount

    !Get the obs and state dimension
    !-------------------------------
    ni = size(mod_lon, 1)
    nj = size(mod_lon, 2)    
    ns = ni * nj
    no = size(obs_lon, 1)    
    
    !Calculate interpolation weight using BUMP
    !-----------------------------------------

    !Important namelist options
    call self%bump%nam%init()

    self%bump%nam%prefix = bump_nam_prefix   ! Prefix for files output
    self%bump%nam%obsop_interp = 'bilin'     ! Interpolation type (bilinear)
    self%bump%nam%default_seed = .true.
    self%bump%nam%new_obsop = .true.
    
    !Initialize geometry
    allocate(area(ns))
    allocate(vunit(ns,1))
    area = 1.0           ! Dummy area
    vunit = 1.0          ! Dummy vertical unit

    !Allocate temporary arrays
    allocate(tmp_lonmod(ns), tmp_latmod(ns), tmp_maskmod(ni, nj))
    tmp_lonmod(:) = reshape(mod_lon, (/ns/))
    tmp_latmod(:) = reshape(mod_lat, (/ns/))
    tmp_maskmod = .true.
    !where (mod_mask == 0.0)
    !   tmp_maskmod = .false.
    !end where
    tmp_maskmod = reshape(tmp_maskmod, (/ns, 1/))
    
    !Initialize BUMP
    call self%bump%setup_online( f_comm%communicator(), ns, 1, 1, 1,&
         &tmp_lonmod, tmp_latmod, area, vunit, tmp_maskmod(:,1),&
         &nobs=no, lonobs=obs_lon, latobs=obs_lat )

    self%initialized = .true.
    self%nobs = no
    
    !Release memory
    deallocate(area)
    deallocate(vunit)
    deallocate(tmp_lonmod, tmp_latmod, tmp_maskmod)

  end subroutine interp_init

  !--------------------------------------------  
  subroutine interp_apply(self, fld, obs)
    ! Forward interpolation: "fields to obs"
    ! obs = interp(fld)
    use kinds

    implicit none

    class(soca_bumpinterp2d), intent(in) :: self    
    real(kind=kind_real),     intent(in) :: fld(:,:)
    real(kind=kind_real),    intent(out) :: obs(:) 

    ! Locals
    real(kind=kind_real), allocatable :: tmp_fld(:,:)
    real(kind=kind_real), allocatable :: tmp_obs(:,:)    
    integer :: ns

    allocate(tmp_fld(size(fld,1),size(fld,2)))
    allocate(tmp_obs(size(obs,1),1))    
    ns = size(fld,1)*size(fld,2)

    tmp_fld = fld
    tmp_fld = reshape(tmp_fld,(/ns, 1/))
    tmp_obs(:,1) = obs
    call self%bump%apply_obsop(tmp_fld,tmp_obs)
    obs = tmp_obs(:,1)

    deallocate(tmp_fld, tmp_obs)

  end subroutine interp_apply

  !--------------------------------------------  
  subroutine interpad_apply(self, fld, obs)
    ! Backward interpolation: "obs to fields"
    ! fld = interpad(obs)    
    use kinds

    implicit none

    class(soca_bumpinterp2d),    intent(in) :: self    
    real(kind=kind_real),     intent(inout) :: fld(:,:)
    real(kind=kind_real),        intent(in) :: obs(:)

    ! Locals
    real(kind=kind_real), allocatable :: tmp_obs(:)
    integer :: nobs

    nobs = size(obs, 1)
    allocate(tmp_obs(nobs))
    tmp_obs = obs
    call self%bump%apply_obsop_ad(tmp_obs, fld)
    deallocate(tmp_obs)

  end subroutine interpad_apply

  !--------------------------------------------  
  subroutine interp_info(self)

    use fckit_log_module, only : fckit_log
    
    implicit none

    class(soca_bumpinterp2d), intent(in) :: self    
    character(len=160) :: record

    write(record,*) "soca_bumpinterp2d%nobs: ",self%nobs
    call fckit_log%info(record)
    write(record,*) "soca_bumpinterp2d%initialized: ",self%initialized
    call fckit_log%info(record)    

  end subroutine interp_info

  !--------------------------------------------  
  subroutine interp_exit(self)

    implicit none

    class(soca_bumpinterp2d), intent(out) :: self

    self%nobs = 0
    if (allocated(self%lono)) deallocate(self%lono)
    if (allocated(self%lato)) deallocate(self%lato)
    call self%bump%dealloc()
    self%initialized = .false.
    
  end subroutine interp_exit

end module soca_bumpinterp2d_mod

