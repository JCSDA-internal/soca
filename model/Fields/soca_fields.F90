
!> Handle fields for the  model

module soca_fields

  use config_mod
  use soca_geom_mod
  !use soca_goms_mod
  !use soca_locs_mod  
  use soca_vars_mod
  use type_linop
  !use tools_interp, only: interp_horiz
  !use type_randgen, only: rng,initialize_sampling,create_randgen
  !use module_namelist, only: namtype  
  use kinds
  use atmos_model_mod,         only: atmos_data_type
  use land_model_mod,          only: land_data_type    
  use ocean_model_mod,         only: ocean_public_type, ocean_state_type, ice_ocean_boundary_type
  use ice_model_mod,           only: ice_data_type
  use soca_mom6sis2, only : Coupled
  use MOM, only : MOM_control_struct
  implicit none
  private

  public :: soca_field, &
       & create, delete, zeros, dirac, random, copy, &
       & self_add, self_schur, self_sub, self_mul, axpy, &
       & dot_prod, add_incr, diff_incr, &
       & read_file, write_file, gpnorm, fldrms, &
       & change_resol, interp_tl, interp_ad, convert_to_ug, convert_from_ug
  public :: soca_field_registry

  ! ------------------------------------------------------------------------------

  !> Fortran derived type to hold fields
  type :: soca_field
     type (Coupled)                    :: AOGCM
     type(soca_geom), pointer          :: geom                  !< MOM5 & CICE5 Geometry
     integer                           :: nf                    !< Number of fields
     character(len=128)                :: gridfname             !< Grid file name
     character(len=128)                :: cicefname             !< Fields file name for cice
     character(len=128)                :: momfname              !< Fields file name for mom
     real(kind=kind_real), pointer     :: cicen(:,:,:) => NULL()  !< Sea-ice fraction                 (nx,ny,ncat)
     real(kind=kind_real), pointer     :: hicen(:,:,:) => NULL()          !< Sea-ice thickness                (nx,ny,ncat)
     real(kind=kind_real), pointer     :: vicen(:,:,:) => NULL()          !< Sea-ice volume                   (nx,ny,ncat)
     real(kind=kind_real), pointer     :: hsnon(:,:,:) => NULL()          !< Snow depth over sea-ice          (nx,ny,ncat)
     real(kind=kind_real), pointer     :: vsnon(:,:,:) => NULL()          !< Snow volume over sea-ice         (nx,ny,ncat) 
     real(kind=kind_real), pointer     :: tsfcn(:,:,:) => NULL()          !< Temperature over sea-ice or snow (nx,ny,ncat)
     real(kind=kind_real), pointer     :: qsnon(:,:,:,:) => NULL()        !< Enthalpy of snow                 (nx,ny,ncat)
     real(kind=kind_real), pointer     :: sicnk(:,:,:,:) => NULL()        !< Salin_wity of sea-ice            (nx,ny,ncat,nzi)
     real(kind=kind_real), pointer     :: socn(:,:,:) => NULL()           !< Ocean (surface) Salinity         (nx,ny,nzo)
     real(kind=kind_real), pointer     :: qicnk(:,:,:,:) => NULL()   !< Enthalpy of sea-ice              (nx,ny,ncat,nzi)
     real(kind=kind_real), pointer     :: tocn(:,:,:) => NULL()      !< Average temperature of grid cell (nx,ny,nzo)
     real(kind=kind_real), pointer     :: ssh(:,:) => NULL()         !< Sea-surface height (nx,ny,nzo)     
     character(len=5), allocatable     :: fldnames(:)           !< Variable identifiers             (nf)
     integer, allocatable              :: numfld_per_fldname(:) !< Number of 2d fields for each     (nf) 
                                                                !< element of fldnames

     type(linoptype)                   :: hinterp_op
     logical                           :: hinterp_initialized = .false. !True:  hinterp_op has been initialized
                                                                        !False: hinterp_op not initialized
  end type soca_field

#define LISTED_TYPE soca_field

  !> Linked list interface - defines registry_t type
#include "Utils/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_field_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "Utils/linkedList_c.f"

  ! ------------------------------------------------------------------------------

  subroutine create(self, geom, vars)

    use soca_mom6sis2, only: Coupled, soca_field_init, soca_geom_init
    use mpp_io_mod,              only: mpp_open, mpp_close
    use SIS_hor_grid, only: set_hor_grid, SIS_hor_grid_type
    !use SIS_get_input, only:directories, Get_SIS_Input
    use MOM_get_input,            only : directories    
    use MOM_file_parser, only : open_param_file, param_file_type, read_param
    use MOM, only : MOM_control_struct, initialize_MOM
    use MOM_time_manager,         only : time_type
    use ocean_model_mod,         only: update_ocean_model, ocean_model_init,  ocean_model_end
    use ice_model_mod,           only: ice_model_init, share_ice_domains, ice_model_end, ice_model_restart
    use ice_grid, only : ice_grid_type, set_ice_grid
    !use constants_mod,           only: constants_init
    !use atmos_model_mod,         only: atmos_data_type
    !use land_model_mod,          only: land_data_type    
    !use ocean_model_mod,         only: ocean_public_type, ocean_state_type, ice_ocean_boundary_type
    !use ice_model_mod,           only: ice_data_type
    use fms_io_mod,      only: fms_io_init, fms_io_exit
    !use MOM_time_manager,         only : time_type
    use time_manager_mod,        only: set_date, get_date, days_in_month, month_name, set_time, set_calendar_type, GREGORIAN
    
    implicit none
    type(soca_field), intent(inout)          :: self
    type(soca_geom),  pointer, intent(inout) :: geom
    type(soca_vars),  intent(in)             :: vars        

    !Stuff for sea-ice grid init
    !type(SIS_hor_grid_type) :: G        !< The horizontal grid type
    !type(param_file_type)    :: param_file !< Parameter file handle
    !type(hor_index_type) :: HI !< A hor_index_type for array extents
    !logical :: global_indexing !< If true use global index
                             !! values instead of having the data domain on each
                             !! processor start at 1.

    !type(time_type)   :: Time
    type(param_file_type) :: param_file
    type(directories) :: path
    !type(MOM_control_struct), pointer       :: CS
    
    !logical,               optional, intent(in)  :: check_params
    !character(len=*),      optional, intent(in)  :: component
    
    integer :: ivar, unit, nxny(2)
    real(kind=kind_real) :: kg_m2_to_H

    !integer :: date_init(6) = (/ 2016, 01, 01, 0, 0, 0 /) 
    
    !type(ice_grid_type) :: IG
    !type(ice_data_type) :: Ice            !< The ice data type that is being initialized.
    !type(time_type)       :: Time_Init      !< The starting time of the model integration
    !type(time_type)       :: Time           !< The current time
    !type(time_type)       :: Time_step_fast !< The time step for the ice_model_fast
    !type(time_type)      :: Time_step_slow !< The time step for the ice_model_slow    
    !call fms_io_init

    self%geom => geom
    self%nf   = vars%nv
    call soca_field_init(self%AOGCM, geom%ocean%G, geom%ocean%GV, geom%ocean%IG)

    ! Finish initializing ice grid

    ! Parse grid inputs
    !call Get_MOM_Input(param_file, dirs)    
    !call set_ice_grid(IG, param_file, NCat_dflt)
    !call set_calendar_type(GREGORIAN)
    !Time_init = set_date (date_init(1), date_init(2), date_init(3), &
    !     date_init(4), date_init(5), date_init(6))
    !Time = Time_init
    !Time_step_fast = set_time(1,0)
    !Time_step_slow = Time_step_fast
    !call ice_model_init(Ice, Time_Init, Time, Time_step_fast, Time_step_slow) 
    !call fms_io_exit
    
    !kg_m2_to_H = self%AOGCM%Ice%fCS%IG%kg_m2_to_H

    ! Assign convenience pointers
    !Ocean internal state    
    self%tocn => self%AOGCM%Ocn%T
    self%socn => self%AOGCM%Ocn%S
    self%ssh => self%AOGCM%Ocn%ssh

    !Sea-ice internal state
    self%cicen => self%AOGCM%Ice%part_size !(:,:,2:)
    self%hicen => self%AOGCM%Ice%h_ice
    self%hsnon => self%AOGCM%Ice%h_snow        
    self%tsfcn => self%AOGCM%Ice%T_skin
    self%sicnk => self%AOGCM%Ice%sal_ice
    self%qicnk => self%AOGCM%Ice%enth_ice
    self%qsnon => self%AOGCM%Ice%enth_snow

    call zeros(self)
    
    allocate(self%numfld_per_fldname(vars%nv))
    
    do ivar=1,vars%nv
       select case(vars%fldnames(ivar))
       case ('cicen','hicen','vicen','hsnon','vsnon','tsfcn')
          self%numfld_per_fldname(ivar)=geom%ocean%ncat
       case ('sicnk','qicnk')
          self%numfld_per_fldname(ivar)=geom%ocean%ncat*geom%ocean%nzi
       case ('qsnon')
          self%numfld_per_fldname(ivar)=geom%ocean%ncat*geom%ocean%nzs
       case ('socn','tocn')
          self%numfld_per_fldname(ivar)=geom%ocean%nzo
       case ('ssh')
         self%numfld_per_fldname(ivar)=1
       case default
          call abor1_ftn("c_soca_fields: undefined variables")
       end select
    end do

    if (self%nf>11) then
       call abor1_ftn ("soca_fields:create error number of fields")       
    endif
    allocate(self%fldnames(self%nf))
    self%fldnames(:)=vars%fldnames(:)

    call check(self)

  end subroutine create

  ! ------------------------------------------------------------------------------

  subroutine delete(self)
    use soca_mom6sis2    
    implicit none
    type(soca_field), intent(inout) :: self

    call soca_field_end(self%AOGCM)!, geom%ocean%G, geom%ocean%GV, geom%ocean%IG)
    
  end subroutine delete

  ! ------------------------------------------------------------------------------

  subroutine zeros(self)
    implicit none
    type(soca_field), intent(inout) :: self

    call check(self)

    self%cicen = 0.0_kind_real
    self%hicen = 0.0_kind_real
    self%hsnon = 0.0_kind_real
    self%tsfcn = 0.0_kind_real
    self%qsnon = 0.0_kind_real
    self%sicnk = 0.0_kind_real
    self%qicnk = 0.0_kind_real

    self%socn = 0.0_kind_real        
    self%tocn = 0.0_kind_real
    self%ssh = 0.0_kind_real
  end subroutine zeros

  ! ------------------------------------------------------------------------------

  subroutine dirac(self, c_conf)
    use iso_c_binding
    implicit none
    type(soca_field), intent(inout) :: self
    type(c_ptr), intent(in)       :: c_conf   !< Configuration
    integer :: ndir,idir,ildir,ifdir,ioff
    integer,allocatable :: ixdir(:),iydir(:)
    character(len=3) :: idirchar

    call check(self)

    ! Get Diracs positions
    ndir = config_get_int(c_conf,"ndir")
    allocate(ixdir(ndir))
    allocate(iydir(ndir))
    do idir=1,ndir
       write(idirchar,'(i3)') idir
       ixdir(idir) = config_get_int(c_conf,"ixdir("//trim(adjustl(idirchar))//")")
       iydir(idir) = config_get_int(c_conf,"iydir("//trim(adjustl(idirchar))//")")
    end do
    !ildir = config_get_int(c_conf,"ildir")
    !ifdir = config_get_int(c_conf,"ifdir")

    ! Check 
    !if (ndir<1) call abor1_ftn("qg_fields:dirac non-positive ndir")
    !if (any(ixdir<1).or.any(ixdir>self%geom%nx)) call abor1_ftn("qg_fields:dirac invalid ixdir")
    !if (any(iydir<1).or.any(iydir>self%geom%ny)) call abor1_ftn("qg_fields:dirac invalid iydir")
    !if ((ildir<1).or.(ildir>self%nl)) call abor1_ftn("qg_fields:dirac invalid ildir")
    !if ((ifdir<1).or.(ifdir>self%nf)) call abor1_ftn("qg_fields:dirac invalid ifdir")

    ! Setup Diracs
    call zeros(self)
    !ioff = (ifdir-1)*self%nl
    do idir=1,ndir
       !self%qicnk(ixdir(idir),iydir(idir),1,4) = 1.0 ! Surface temp incr for cat 1
       !self%tsfcn(ixdir(idir),iydir(idir),1) = 1.0 ! Surface temp incr for cat 1
       !self%tocn(ixdir(idir),iydir(idir)) = 1.0 ! Surface temp incr for cat 1
       !self%cicen(ixdir(idir),iydir(idir),3) = 1.0 ! Surface temp incr for cat 1
    end do

  end subroutine dirac

  ! ------------------------------------------------------------------------------

  ! ------------------------------------------------------------------------------

  subroutine random(self)

    implicit none
    type(soca_field), intent(inout) :: self

    call check(self)

    !call random_number(self%cicen); self%cicen=self%cicen-0.5_kind_real
    !call random_number(self%hicen); self%hicen=self%hicen-0.5_kind_real
    !call random_number(self%hsnon); self%hsnon=self%hsnon-0.5_kind_real        
    !call random_number(self%tsfcn); self%tsfcn=self%tsfcn-0.5_kind_real
    
    !call random_number(self%tocn); self%tocn=self%tocn-0.5_kind_real
    !call random_number(self%socn); self%tocn=self%socn-0.5_kind_real
    !self%ssh=0.1_kind_real
    call random_number(self%ssh); self%ssh=(self%ssh-0.5_kind_real)*self%geom%ocean%mask2d
    
  end subroutine random

  ! ------------------------------------------------------------------------------

  subroutine copy(self,rhs)
    implicit none
    type(soca_field), intent(inout) :: self
    type(soca_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen = rhs%cicen
    self%hicen = rhs%hicen
    self%hsnon = rhs%hsnon
    self%tsfcn = rhs%tsfcn
    self%qsnon = rhs%qsnon
    self%sicnk = rhs%sicnk
    self%qicnk = rhs%qicnk    

    self%socn  = rhs%socn
    self%tocn  = rhs%tocn
    self%ssh  = rhs%ssh    

    !call linop_copy(rhs%hinterp_op, self%hinterp_op)
    
    return
  end subroutine copy

  ! ------------------------------------------------------------------------------

  subroutine self_add(self,rhs)

    implicit none

    type(soca_field), intent(inout) :: self
    type(soca_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen = self%cicen + rhs%cicen
    self%hicen = self%hicen + rhs%hicen
    self%hsnon = self%hsnon + rhs%hsnon
    self%tsfcn = self%tsfcn + rhs%tsfcn
    self%qsnon = self%qsnon + rhs%qsnon
    self%sicnk = self%sicnk + rhs%sicnk
    self%qicnk = self%qicnk + rhs%qicnk

    self%tocn = self%tocn + rhs%tocn
    self%socn = self%socn + rhs%socn
    self%ssh = self%ssh + rhs%ssh

    !self%AOGCM%Ocn%ssh = self%AOGCM%Ocn%ssh + rhs%AOGCM%Ocn%ssh
    !return
  end subroutine self_add

  ! ------------------------------------------------------------------------------

  subroutine self_schur(self,rhs)
    implicit none
    type(soca_field), intent(inout) :: self
    type(soca_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen=self%cicen*rhs%cicen
    self%hicen=self%hicen*rhs%hicen
    self%hsnon=self%hsnon*rhs%hsnon
    self%tsfcn=self%tsfcn*rhs%tsfcn
    self%qsnon=self%qsnon*rhs%qsnon
    self%sicnk=self%sicnk*rhs%sicnk
    self%qicnk=self%qicnk*rhs%qicnk

    self%tocn=self%tocn*rhs%tocn
    self%socn=self%socn*rhs%socn
    self%ssh=self%ssh*rhs%ssh    
    
    return
  end subroutine self_schur

  ! ------------------------------------------------------------------------------

  subroutine self_sub(self,rhs)
    implicit none
    type(soca_field), intent(inout) :: self
    type(soca_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen=self%cicen-rhs%cicen
    self%hicen=self%hicen-rhs%hicen
    self%hsnon=self%hsnon-rhs%hsnon
    self%tsfcn=self%tsfcn-rhs%tsfcn
    self%qsnon=self%qsnon-rhs%qsnon
    self%sicnk=self%sicnk-rhs%sicnk
    self%qicnk=self%qicnk-rhs%qicnk
    
    self%socn=self%socn-rhs%socn
    self%tocn=self%tocn-rhs%tocn
    self%ssh=self%ssh-rhs%ssh    

  end subroutine self_sub

  ! ------------------------------------------------------------------------------

  subroutine self_mul(self,zz)
    implicit none
    type(soca_field), intent(inout) :: self
    real(kind=kind_real), intent(in) :: zz

    call check(self)

    self%cicen = zz * self%cicen
    self%hicen = zz * self%hicen
    self%hsnon = zz * self%hsnon
    self%tsfcn = zz * self%tsfcn
    self%qsnon = zz * self%qsnon
    self%sicnk = zz * self%sicnk
    self%qicnk = zz * self%qicnk
    
    self%tocn = zz * self%tocn
    self%socn = zz * self%socn
    self%ssh = zz * self%ssh
    
  end subroutine self_mul

  ! ------------------------------------------------------------------------------

  subroutine axpy(self,zz,rhs)
    implicit none
    type(soca_field), intent(inout) :: self
    real(kind=kind_real), intent(in) :: zz
    type(soca_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen=self%cicen + zz * rhs%cicen
    self%hicen=self%hicen + zz * rhs%hicen
    self%hsnon=self%hsnon + zz * rhs%hsnon
    self%tsfcn=self%tsfcn + zz * rhs%tsfcn
    self%qsnon=self%qsnon + zz * rhs%qsnon
    self%sicnk=self%sicnk + zz * rhs%sicnk
    self%qicnk=self%qicnk + zz * rhs%qicnk
    
    self%tocn=self%tocn + zz * rhs%tocn        
    self%socn=self%socn + zz * rhs%socn
    self%ssh=self%ssh + zz * rhs%ssh
    
  end subroutine axpy

  ! ------------------------------------------------------------------------------

  subroutine dot_prod(fld1,fld2,zprod)
    
    use mpp_mod,  only : mpp_pe, mpp_npes, mpp_root_pe, mpp_sync, mpp_sum, mpp_gather, mpp_broadcast
    
    implicit none
    type(soca_field), intent(in) :: fld1, fld2
    real(kind=kind_real), intent(out) :: zprod
    real(kind=kind_real),allocatable,dimension(:) :: zprod_allpes
    integer :: ii, jj, kk
    integer :: is, ie, js, je
    call check_resolution(fld1, fld2)
    if (fld1%nf /= fld2%nf .or. fld1%geom%ocean%nzo /= fld2%geom%ocean%nzo) then
       call abor1_ftn("soca_fields:field_prod error number of fields")
    endif

    ! Get compute domain
    is = fld1%geom%ocean%G%isc    
    ie = fld1%geom%ocean%G%iec
    js = fld1%geom%ocean%G%jsc    
    je = fld1%geom%ocean%G%jec    

    zprod = 0.0_kind_real
    do ii = is, ie
       do jj = js, je
          zprod = zprod + fld1%ssh(ii,jj)*fld2%ssh(ii,jj)*fld1%geom%ocean%mask2d(ii,jj)
       end do
    end do
    !zprod=sum(fld1%ssh(is:ie,js:je)*fld2%ssh(is:ie,js:je)*fld1%geom%ocean%mask2d(is:ie,js:je))
    allocate(zprod_allpes(mpp_npes()))

    call mpp_gather((/zprod/),zprod_allpes)
    call mpp_broadcast(zprod_allpes, mpp_npes(), mpp_root_pe())    
    zprod=sum(zprod_allpes)
    !print *,' ======================= DOT PROD GLOBAL= ',zprod
    
    deallocate(zprod_allpes)
    call mpp_sync()
    
  end subroutine dot_prod

  ! ------------------------------------------------------------------------------

  subroutine add_incr(self,rhs)
    implicit none
    type(soca_field), intent(inout) :: self
    type(soca_field), intent(in)    :: rhs

    call check(self)
    call check(rhs)

    call self_add(self,rhs)

    return
  end subroutine add_incr

  ! ------------------------------------------------------------------------------

  subroutine diff_incr(lhs,x1,x2)
    implicit none
    type(soca_field), intent(inout) :: lhs
    type(soca_field), intent(in)    :: x1
    type(soca_field), intent(in)    :: x2

    call check(lhs)
    call check(x1)
    call check(x2)

    call zeros(lhs)


    !if (x1%geom%nx==x2%geom%nx .and. x1%geom%ny==x2%geom%ny) then
    !   if (lhs%geom%nx==x1%geom%nx .and. lhs%geom%ny==x1%geom%ny) then
    lhs%cicen = x1%cicen - x2%cicen
    lhs%hicen = x1%hicen - x2%hicen
    lhs%hsnon = x1%hsnon - x2%hsnon
    lhs%tsfcn = x1%tsfcn - x2%tsfcn
    lhs%qsnon = x1%qsnon - x2%qsnon
    lhs%sicnk = x1%sicnk - x2%sicnk
    lhs%qicnk = x1%qicnk - x2%qicnk
    
    lhs%tocn = x1%tocn - x2%tocn
    lhs%socn = x1%socn - x2%socn
    lhs%ssh = x1%ssh - x2%ssh    
    !   else
    !      call abor1_ftn("soca_fields:diff_incr: not coded for low res increment yet")
    !   endif
    !else
    !   call abor1_ftn("soca_fields:diff_incr: states not at same resolution")
    !endif

  end subroutine diff_incr

  ! ------------------------------------------------------------------------------

  subroutine change_resol(fld,rhs)
    implicit none
    type(soca_field), intent(inout) :: fld
    type(soca_field), intent(in)    :: rhs
    real(kind=kind_real), allocatable :: ztmp(:,:)
    real(kind=kind_real) :: dy1, dy2, ya, yb, dx1, dx2, xa, xb
    integer :: jx, jy, jf, iy, ia, ib

    call check(fld)
    call check(rhs)
    call copy(fld,rhs)
    !call abor1_ftn("soca_fields:field_resol: untested code")

    return
  end subroutine change_resol

  ! ------------------------------------------------------------------------------

  subroutine read_file(fld, c_conf, vdate)

    use iso_c_binding
    use datetime_mod
    use fckit_log_module, only : log
    use netcdf
    use soca_thermo
    use soca_mom6sis2
    use fms_mod,                 only: read_data, write_data, set_domain
    use fms_io_mod,                only : fms_io_init, fms_io_exit
    use mpp_mod,  only : mpp_pe, mpp_npes, mpp_root_pe, mpp_sync, mpp_sum, mpp_gather, mpp_broadcast
    use fms_io_mod,       only : register_restart_field, restart_file_type
    use fms_io_mod,       only : restore_state, query_initialized

    use ufo_locs_mod
    use ufo_geovals_mod
    use ufo_vars_mod
    
    implicit none
    type(soca_field), intent(inout) :: fld      !< Fields
    type(c_ptr), intent(in)       :: c_conf   !< Configuration
    type(datetime), intent(inout) :: vdate    !< DateTime

    integer, parameter :: max_string_length=800 ! Yuk!
    character(len=max_string_length) :: ocn_sfc_filename, ocn_filename, ice_filename, basename
    character(len=20) :: sdate
    character(len=1024)  :: buf
    integer :: iread, ii

    type(restart_file_type) :: sis_restart    
    integer :: idr

    
    type(ufo_locs)    :: locs
    type(ufo_geovals)    :: gom
    type(ufo_vars)    :: vars
    integer            :: nobs, nval
    
    iread = 0
    if (config_element_exists(c_conf,"read_from_file")) then
       iread = config_get_int(c_conf,"read_from_file")
    endif
    if (iread==0) then
       call log%warning("soca_fields:read_file: Inventing State")
       call invent_state(fld,c_conf)
       sdate = config_get_string(c_conf,len(sdate),"date")
       WRITE(buf,*) 'validity date is: '//sdate
       call log%info(buf)
       call datetime_set(sdate, vdate)
    else
       print *,fld%fldnames
       basename = config_get_string(c_conf,len(basename),"basename")        
       ocn_sfc_filename = config_get_string(c_conf,len(ocn_filename),"ocn_sfc_filename")
       ocn_filename = config_get_string(c_conf,len(ocn_filename),"ocn_filename")       

       ocn_sfc_filename = trim(basename)//trim(ocn_sfc_filename)
       ocn_filename = trim(basename)//trim(ocn_filename)       
       ice_filename = config_get_string(c_conf,len(ice_filename),"ice_filename")       
       ice_filename = trim(basename)//trim(ice_filename)
       
       call fms_io_init()
       do ii = 1, fld%nf
          print *,fld%fldnames(ii)
          select case(fld%fldnames(ii))
          case ('ssh')
             call read_data(ocn_sfc_filename,"ave_ssh",fld%ssh(:,:),domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('tocn')
             call read_data(ocn_filename,"Temp",fld%tocn(:,:,:),domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('socn')             
             call read_data(ocn_filename,"Salt",fld%tocn(:,:,:),domain=fld%geom%ocean%G%Domain%mpp_domain)             
          case ('cicen')
             idr = register_restart_field(sis_restart, ice_filename, 'part_size', fld%AOGCM%Ice%part_size, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('hicen')
             idr = register_restart_field(sis_restart, ice_filename, 'h_ice', fld%AOGCM%Ice%h_ice, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('hsnon')
             idr = register_restart_field(sis_restart, ice_filename, 'h_snow', fld%AOGCM%Ice%h_snow, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)             
          case ('qicnk')             
             idr = register_restart_field(sis_restart, ice_filename, 'enth_ice', fld%AOGCM%Ice%enth_ice, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('qsnon')             
             idr = register_restart_field(sis_restart, ice_filename, 'enth_snow', fld%AOGCM%Ice%enth_snow, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('tsfcn')             
             idr = register_restart_field(sis_restart, ice_filename, 'T_skin', fld%AOGCM%Ice%T_skin, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('sicen')             
             idr = register_restart_field(sis_restart, ice_filename, 'sal_ice', fld%AOGCM%Ice%sal_ice, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)             
          case default
             print *,'Not reading var ',fld%fldnames(ii),' in file ',ocn_filename
             !call log%warning("soca_fields:read_file: Not reading var "//fld%fldnames(ii))
             !call abor1_ftn("soca_fields: undefined variables")             
          end select
       end do
       call restore_state(sis_restart, directory='')

       sdate = config_get_string(c_conf,len(sdate),"date")
       WRITE(buf,*) 'validity date is: '//sdate
       call log%info(buf)
       call datetime_set(sdate, vdate)       

    endif

    print *,'================ HACK TO TEST INTERP ================='

    nobs = 2
    print *,'Init locs'
    call ufo_locs_setup(locs, 2)
    locs%nlocs = nobs
    locs%lon(1) = 15.3_kind_real
    locs%lat(1) = 62.3_kind_real

    locs%lon(2) = 17.2_kind_real
    locs%lat(2) = 63.4_kind_real

    print *,'Init ufo vars'
    call ufo_vars_setup(vars, (/var_seaicefrac/))
    
    print *,'Init gom'
    call ufo_geovals_init(gom)
    call ufo_geovals_setup(gom, vars, nobs)
    nval = fld%geom%ocean%ncat
    gom%geovals(1)%nval = nval
    allocate(gom%geovals(1)%vals(nval,nobs))
    allocate(gom%geovals(1)%vals(nval,nobs))
    gom%linit = .true.    
    call ufo_geovals_zero(gom)
    print *,gom%lalloc
    !call ufo_geovals_print(gom, 1)
    !gom%nvar=1
    !gom%nobs=1
    !allocate(gom%variables(gom%nvar))
    !gom%variables(1)="cicen"
    print *,'================ INTERP ================='    
    call interp_tl(fld, locs, gom)
    
    call check(fld)
    !call mpp_sync()
    
  end subroutine read_file

  ! ------------------------------------------------------------------------------

  subroutine write_file(fld, c_conf, vdate)
    use iso_c_binding
    use datetime_mod
    use fckit_log_module, only : fckit_log
    use netcdf
    use fms_mod,                 only: read_data, write_data, set_domain
    use fms_io_mod,                only : fms_io_init, fms_io_exit
    
    implicit none
    type(soca_field), intent(inout) :: fld    !< Fields
    type(c_ptr), intent(in)    :: c_conf           !< Configuration
    type(datetime), intent(inout) :: vdate         !< DateTime
    integer, parameter :: max_string_length=800    ! Yuk!
    character(len=max_string_length) :: filename
    character(len=1024):: buf
    integer :: ii
    call check(fld)

    filename = genfilename(c_conf,max_string_length,vdate)
    WRITE(buf,*) 'field:write_file: writing '//filename
    call fckit_log%info(buf)
    
    call fms_io_init()
    call set_domain( fld%geom%ocean%G%Domain%mpp_domain )    
    do ii = 1, fld%nf
       print *,fld%fldnames(ii)
       select case(fld%fldnames(ii))
       case ('ssh')
          call write_data( filename, "ssh", fld%ssh, fld%geom%ocean%G%Domain%mpp_domain)              
       case ('tocn')
          call write_data( filename, "temp", fld%tocn, fld%geom%ocean%G%Domain%mpp_domain)          
       case ('socn')
          call write_data( filename, "salt", fld%socn, fld%geom%ocean%G%Domain%mpp_domain)          
       case ('hicen')
          call write_data( filename, "hicen", fld%hicen, fld%geom%ocean%G%Domain%mpp_domain)
       case ('hsnon')
          call write_data(filename, "hsnon", fld%hsnon, fld%geom%ocean%G%Domain%mpp_domain)
       case ('cicen')
          call write_data(filename, "cicen", fld%cicen, fld%geom%ocean%G%Domain%mpp_domain)
       case ('qicnk')
          call write_data(filename, "qicnk", fld%qicnk, fld%geom%ocean%G%Domain%mpp_domain)
       case ('sicnk')
          call write_data(filename, "sicnk", fld%sicnk, fld%geom%ocean%G%Domain%mpp_domain)
       case ('qsnon')
          call write_data(filename, "qsnon", fld%qsnon, fld%geom%ocean%G%Domain%mpp_domain)
       case ('tsfcn')
          call write_data(filename, "tsfcn", fld%tsfcn, fld%geom%ocean%G%Domain%mpp_domain)                    
       case default
          print *,'Not writing var ',fld%fldnames(ii),' in file ',filename
          !call log%warning("soca_fields:read_file: Not reading var "//fld%fldnames(ii))
          !call abor1_ftn("soca_fields: undefined variables")             
       end select
    end do
    call fms_io_exit()       
    
    print *,'===================== Done writting' 
    
  end subroutine write_file

  ! ------------------------------------------------------------------------------

  function genfilename (c_conf,length,vdate)
    use iso_c_binding
    use datetime_mod
    use duration_mod
    type(c_ptr), intent(in)    :: c_conf  !< Configuration
    integer, intent(in) :: length
    character(len=length) :: genfilename
    type(datetime), intent(in) :: vdate

    character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
         & prefix, mmb
    type(datetime) :: rdate
    type(duration) :: step
    integer lenfn

    ! here we should query the length and then allocate "string".
    ! But Fortran 90 does not allow variable-length allocatable strings.
    ! config_get_string checks the string length and aborts if too short.
    fdbdir = config_get_string(c_conf,len(fdbdir),"datadir")
    expver = config_get_string(c_conf,len(expver),"exp")
    typ    = config_get_string(c_conf,len(typ)   ,"type")

    if (typ=="ens") then
       mmb = config_get_string(c_conf, len(mmb), "member")
       lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
       prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ) // "." // TRIM(mmb)
    else
       lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
       prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ)
    endif

    if (typ=="fc" .or. typ=="ens") then
       referencedate = config_get_string(c_conf,len(referencedate),"date")
       call datetime_to_string(vdate, validitydate)
       call datetime_create(TRIM(referencedate), rdate)
       call datetime_diff(vdate, rdate, step)
       call duration_to_string(step, sstep)
       lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
       genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
    endif

    if (typ=="an") then
       call datetime_to_string(vdate, validitydate)
       lenfn = lenfn + 1 + LEN_TRIM(validitydate)
       genfilename = TRIM(prefix) // "." // TRIM(validitydate)
    endif

    if (lenfn>length) &
         & call abor1_ftn("fields:genfilename: filename too long")

  end function genfilename

  ! ------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------

  subroutine gpnorm(fld, nf, pstat) 
    implicit none
    type(soca_field), intent(in) :: fld
    integer, intent(in) :: nf
    real(kind=kind_real), intent(inout) :: pstat(3, nf) !> [average, min, max]
    real(kind=kind_real) :: zz
    integer :: jj,joff

    call check(fld)
    
    pstat(1,:) = minval(fld%ssh)
    pstat(2,:) = maxval(fld%ssh)
    !call fldrms(fld, zz)
    call dot_prod(fld, fld, zz)    
    pstat(3,:) = sqrt(zz)
    print *,'-------------------- GPNORM = ',zz
    !call abor1_ftn("soca_fields_gpnorm: error not implemented")

  end subroutine gpnorm

  ! ------------------------------------------------------------------------------

  subroutine fldrms(fld, prms)
    
    use mpp_mod,                   only : mpp_pe, mpp_npes, mpp_root_pe, mpp_sync
    
    implicit none
    
    type(soca_field), intent(in) :: fld
    real(kind=kind_real), intent(out) :: prms
    integer :: jf,jy,jx,ii
    real(kind=kind_real) :: zz, ns, n2dfld

    
    call check(fld)
    
    call dot_prod(fld,fld,prms) ! Global value 
    prms=sqrt(prms)
    print *,'PE# ',mpp_pe(),' ---------------- NORM = ',prms

  end subroutine fldrms

  ! ------------------------------------------------------------------------------

  subroutine interp_tl(fld, locs, gom)
    use ufo_locs_mod
    use ufo_geovals_mod
    implicit none
    type(soca_field), intent(inout)   :: fld
    type(ufo_locs), intent(in)    :: locs
    type(ufo_geovals), intent(inout) :: gom
    character(2)                        :: op_type='TL'

    call check(fld)

    call nicas_interph(fld, locs, gom, op_type)

  end subroutine interp_tl

  ! ------------------------------------------------------------------------------

  subroutine interp_ad(fld, locs, gom)
    use ufo_locs_mod
    use ufo_geovals_mod    
    implicit none
    type(soca_field), intent(inout) :: fld
    type(ufo_locs), intent(in)    :: locs
    type(ufo_geovals), intent(inout) :: gom    
    character(2)                        :: op_type='AD'

    call check(fld)
    call nicas_interph(fld, locs, gom, op_type)

  end subroutine interp_ad

  ! ------------------------------------------------------------------------------

  subroutine nicas_interph(fld, locs, gom, op_type)

    use type_linop
    !use tools_interp, only: compute_interp
    !use type_randgen, only: initialize_sampling, create_randgen!, randgentype
    !use type_nam, only: namtype
    !use tools_const, only : deg2rad
    !use horiz_interp_mod, only : horiz_interp_type, horiz_interp_new, horiz_interp
    !use horiz_interp_mod, only : horiz_interp_init, horiz_interp_end
    !use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_init, horiz_interp_bilinear
    use ufo_locs_mod_c  
    use ufo_locs_mod  
    use ufo_geovals_mod_c
    use ufo_geovals_mod
    
    type(soca_field), intent(inout)    :: fld
    type(ufo_locs), intent(in)     :: locs
    type(ufo_geovals), intent(inout)  :: gom
    character(2), intent(in)             :: op_type !('TL' or 'AD')

    integer :: Nc, No, var_type_index, Ncat
    integer :: ivar, gom_dim1, cnt_fld
    character(len=1024)  :: buf
    logical,allocatable :: mask(:), masko(:)               ! < mask (ncells, nlevels)
    real(kind=kind_real), allocatable :: lon(:), lat(:), lono(:), lato(:), fld_src(:), fld_dst(:)
    !type(namtype) :: nam !< Namelist variables
    !type(linoptype) :: hinterp_op
    

    Nc = fld%geom%ocean%nx*fld%geom%ocean%ny
    No = locs%nlocs   !< DOES NOT SEEM RIGHT, SHOULD BE TOTAL OBS IN da WINDOW
    Ncat = fld%geom%ocean%ncat
    if (No>0) then
       allocate(lon(Nc), lat(Nc), mask(Nc), fld_src(Nc))    ! <--- Not memory efficient ...
       allocate(masko(No), fld_dst(No), lono(No), lato(No)) ! <--- RECODE

       masko = .true.
       mask = .true.
       fld%hinterp_initialized = .false.
       if (.not.(fld%hinterp_initialized)) then
          print *,'INITIALIZE INTERP',gom%nobs, locs%nlocs
          print *,'LOCS=',locs%lat
          
          !lono = deg2rad*locs%lon(:)
          !lato = deg2rad*locs%lat(:)
          !lon = deg2rad*reshape(fld%geom%ocean%lon, (/Nc/))     ! Inline grid, structured to un-structured
          !lat = deg2rad*reshape(fld%geom%ocean%lat, (/Nc/))     ! and change to SI Units
          
          lono = locs%lon(:)
          lato = locs%lat(:)          
          lon = reshape(fld%geom%ocean%lon, (/Nc/))     ! Inline grid, structured to un-structured
          lat = reshape(fld%geom%ocean%lat, (/Nc/))     ! and change to SI Units


          ! Ben's interp initialization (BROKEN)
          !rng = create_randgen(nam)
          !call compute_interp(Nc, lon, lat, mask, No, lono, lato, masko, interp_type, fld%hinterp_op)          
          fld%hinterp_initialized = .true.

       end if

       !Finish Initializing gom (right place?)
!!$       if (.not.(gom%alloc)) then
!!$          
!!$          print *,'fld%numfld_per_fldname=',fld%numfld_per_fldname
!!$          print *,'gom%nvar=',gom%nvar
!!$          gom_dim1=sum(fld%numfld_per_fldname(1:gom%nvar)) ! WILL CREATE ISSUES:
!!$                                                           ! Assume the order of var type is preserved
!!$                                                           ! [cicen, hicen, ...]
!!$          print *,'gom_dim1=',gom_dim1
!!$          allocate(gom%geovals(gom_dim1,gom%nobs))
!!$       end if
!!$       if (.not.allocated(gom%numfld_per_fldname)) then
!!$          allocate(gom%numfld_per_fldname(gom%nvar))
!!$       end if
!!$       gom%numfld_per_fldname=fld%numfld_per_fldname ! Will be used in obs oper          
!!$       !end if !probably need to assert shape of gom%values==(gom_dim1,gom%nobs)

       select case (op_type)
       case ('TL')
!!$          if (.not.allocated(fld_src)) allocate(fld_src(Nc))
!!$          cnt_fld=0
!!$          !Loop through variable types
!!$          do var_type_index=1,gom%nvar
!!$             !Loop through variable fields
!!$             do ivar=1,fld%numfld_per_fldname(var_type_index)
!!$                cnt_fld=cnt_fld+1
!!$                print *,'Apply interp op to variable:',gom%variables(var_type_index),' field num:',cnt_fld
!!$                fld_src = reshape(fld%cicen(:,:,ivar), (/Nc/))
!!$                fld_dst = 0.0_kind_real
!!$                print *,'fld_dst=',fld_dst
!!$
!!$                print *,'shape src:',shape( fld%geom%ocean%lon )
!!$                print *,'shape field:', shape(fld%cicen(:,:,ivar))
!!$                print *,'shape dst:',shape( lono )                
!!$                
!!$                !call horiz_interp(fms_interp, fld%cicen(:,:,ivar), fld_dst)
!!$
!!$                print *,'out of interp:',fld_dst
!!$                !call apply_linop(fld%hinterp_op, fld_src, fld_dst)
!!$                gom%values(cnt_fld,gom%used:gom%used+No-1)=fld_dst(1:No)
!!$
!!$             end do
!!$          end do
!!$          deallocate(fld_src)
       case ('AD')
          call abor1_ftn("nicas_interph: Wrapper for adjoint not implemented yet")
          !call apply_linop_ad(hinterp_op,fld_dst,fld_src)
          !CHECK: WHERE DO WE TRIGGER THE ADJOINT TESTS? 
          !put fld_src
       end select
!!$       gom%used=gom%used+locs%nlocs
       print *,"========== OUT OF INTERP ====================="
    end if
  end subroutine nicas_interph

  ! ------------------------------------------------------------------------------

  subroutine convert_to_ug(self, ug)
    use unstructured_grid_mod
    use soca_thermo

    implicit none
    type(soca_field), intent(in) :: self
    type(unstructured_grid), intent(inout) :: ug
    real(kind=kind_real), allocatable :: zz(:)
    real(kind=kind_real), allocatable :: vv(:)    
    integer, allocatable :: cmask(:)
    integer :: jx,jy,jz,jk
    integer :: nz_total     ! Total number of levels in the 3D fields
    integer :: n_vars       ! Number of 3D variables 
    integer :: n_surf_vars  ! Number of surf vars (sould be 0 for ocean/ice)
    integer :: cat_num      ! !!!!!!!!! only doing 1 category for now !!!!!!!!!!!

!!$    cat_num = 1
!!$    nz_total = size(self%geom%level)
!!$    allocate(zz(nz_total))
!!$    allocate(vv(nz_total))
!!$    allocate(cmask(nz_total))
!!$    do jz = 1,nz_total
!!$       zz(jz) = real(self%geom%level(jz))
!!$    end do
!!$    call create_unstructured_grid(ug, nz_total, zz)

!!$    n_vars = 1      !!!!! START WITH ONLY ONE VAR !!!!!!!!! 
!!$    n_surf_vars = 0 !!!!! NO SURFACE VAR !!!!!!!!! 
!!$
!!$    do jy=1,self%geom%ny
!!$       do jx=1,self%geom%nx
!!$          jk = 1
!!$          cmask(jk) = int(self%geom%mask(jx,jy))       ! Surface T
!!$          vv(jk) = self%tsfcn(jx,jy,cat_num)
!!$          jk = jk + 1
!!$          do jz = 1,self%geom%nzs                              ! Snow T
!!$             cmask(jk) = int(self%geom%mask(jx,jy))    !
!!$             !vv(jk) = Ts_nl(self%qsnon(jx,jy,cat_num))
!!$             vv(jk) = self%qsnon(jx,jy,cat_num)             
!!$             jk = jk + 1
!!$          end do
!!$          do jz = 1,self%geom%nzi                              ! Ice T
!!$             cmask(jk) = int(self%geom%mask(jx,jy))    !
!!$             !vv(jk) = Ti_nl(self%qicnk(jx,jy,cat_num,jz),self%sicnk(jx,jy,cat_num,jz))
!!$             vv(jk) = self%qicnk(jx,jy,cat_num,jz)
!!$             jk = jk + 1
!!$          end do
!!$          cmask(jk) = int(self%geom%mask(jx,jy))       ! Ice/Ocean interface
!!$          !vv(jk) = Tm(self%socn(jx,jy))                 ! Tf = -mu * S
!!$          vv(jk) = self%socn(jx,jy)                      ! Tf = -mu * S          
!!$          jk = jk + 1          
!!$          do jz = 1,self%geom%nzo                              ! Ocean
!!$             cmask(jk) = int(self%geom%mask(jx,jy))    !             
!!$             !vv(jk) = self%tocn(jx,jy)                   !
!!$             jk = jk + 1
!!$          end do
!!$
!!$          !cmask(:) = int(self%geom%mask(jx,jy))           ! Some issues with the mask
!!$          !print *,'cmask=',cmask
!!$          !if (self%icemask(jx,jy)>0.0) then
!!$          !print *,vv(:)
!!$          !read(*,*)
!!$          !end if
!!$          !if (cmask(1)==1) read(*,*)
!!$
!!$          call add_column(ug, self%geom%lat(jx,jy), self%geom%lon(jx,jy), self%geom%cell_area(jx, jy), &
!!$               nz_total, &
!!$               n_vars, &
!!$               n_surf_vars, &
!!$               cmask, &
!!$               0)
!!$          ug%last%column%fld3d(:) = vv(:)
!!$       enddo
!!$    enddo
  end subroutine convert_to_ug

  ! ------------------------------------------------------------------------------

  subroutine convert_from_ug(self, ug)
    use unstructured_grid_mod
    implicit none
    type(soca_field), intent(inout) :: self
    type(unstructured_grid), intent(in) :: ug
    !type(column_element), pointer :: current
    real(kind=kind_real) :: dx, dy
    integer :: jx,jy,jz,jk
    integer :: n_vars       ! Number of 3D variables 
    integer :: n_surf_vars  ! Number of surf vars (sould be 0 for ocean/ice)
    integer :: cat_num      ! !!!!!!!!! only doing 1 category for now !!!!!!!!!!!

!!!!!!!!! code inverse of convert_to_ug !!!!!!!!!!!!!!!
!!$
!!$    current => ug%head
!!$    cat_num = 1 
!!$    n_vars = 1      !!!!! START WITH ONLY ONE VAR !!!!!!!!! 
!!$    n_surf_vars = 0 !!!!! NO SURFACE VAR !!!!!!!!! 
!!$
!!$    do jy=1,self%geom%ny
!!$       do jx=1,self%geom%nx
!!$          jk = 1
!!$          self%tsfcn(jx,jy,cat_num) = current%column%fld3d(jk)         ! Tsfcs
!!$          jk = jk + 1
!!$          do jz = 1,self%geom%nzs                                          ! Q Snow
!!$             self%qsnon(jx,jy,cat_num) = current%column%fld3d(jk)
!!$             jk = jk + 1
!!$          end do
!!$          do jz = 1,self%geom%nzi                                          ! Q Ice
!!$             self%qicnk(jx,jy,cat_num,jz) = current%column%fld3d(jk)
!!$             jk = jk + 1
!!$          end do
!!$          self%socn(jx,jy) = current%column%fld3d(jk)                 ! Ice/Ocean interface,
!!$          ! Tf = -mu * S          
!!$          jk = jk + 1          
!!$          do jz = 1,self%geom%nzo                                          ! Ocean SST
!!$             !self%tocn(jx,jy) = current%column%fld3d(jk)              !
!!$             jk = jk + 1
!!$          end do
!!$          current => current%next
!!$       enddo
!!$    enddo

  end subroutine convert_from_ug

  ! ------------------------------------------------------------------------------

  function common_vars(x1, x2)

    implicit none
    type(soca_field), intent(in) :: x1, x2
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
    if (x1%geom%ocean%nzo /= x2%geom%ocean%nzo) call abor1_ftn("common_vars: error number of levels")
    !common_vars = x1%geom%nzi * common_vars

  end function common_vars

  ! ------------------------------------------------------------------------------

  subroutine check_resolution(x1, x2)

    implicit none
    type(soca_field), intent(in) :: x1, x2

    ! NEEDS WORK !!!
    if (x1%geom%ocean%nx /= x2%geom%ocean%nx .or.  x1%geom%ocean%ny /= x2%geom%ocean%ny ) then
       call abor1_ftn ("soca_fields: resolution error")
    endif
    call check(x1)
    call check(x2)

  end subroutine check_resolution

  ! ------------------------------------------------------------------------------

  subroutine check(self)
    implicit none
    type(soca_field), intent(in) :: self
    logical :: bad

    bad = .false.
    !bad = bad .or. (size(self%cicen, 1) /= self%geom%nx)

    ! add more test here ...

    if (bad) then
       write(0,*)'nx, ny, nf, nzi, nzo = ',self%geom%ocean%nx,self%geom%ocean%ny
       call abor1_ftn ("soca_fields: field not consistent")
    endif

  end subroutine check

  ! ------------------------------------------------------------------------------

end module soca_fields
