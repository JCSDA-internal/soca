
!> Handle fields for the  model

module soca_fields

  use config_mod
  use soca_geom_mod
  use soca_goms_mod
  use soca_locs_mod  
  use soca_vars_mod
  use type_linop
  use tools_interp, only: interp_horiz
  !use type_randgen, only: rng,initialize_sampling,create_randgen
  use module_namelist, only: namtype  
  use kinds
  use atmos_model_mod,         only: atmos_data_type
  use land_model_mod,          only: land_data_type    
  use ocean_model_mod,         only: ocean_public_type, ocean_state_type, ice_ocean_boundary_type
  use ice_model_mod,           only: ice_data_type
  use soca_mom6sis2, only : Coupled
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
     real(kind=kind_real), pointer :: cicen(:,:,:)          !< Sea-ice fraction                 (nx,ny,ncat)
     real(kind=kind_real), allocatable :: hicen(:,:,:)          !< Sea-ice thickness                (nx,ny,ncat)
     real(kind=kind_real), pointer :: vicen(:,:,:)          !< Sea-ice volume                   (nx,ny,ncat)
     real(kind=kind_real), allocatable :: hsnon(:,:,:)          !< Snow depth over sea-ice          (nx,ny,ncat)
     real(kind=kind_real), pointer :: vsnon(:,:,:)          !< Snow volume over sea-ice         (nx,ny,ncat) 
     real(kind=kind_real), pointer :: tsfcn(:,:,:)          !< Temperature over sea-ice or snow (nx,ny,ncat)
     real(kind=kind_real), pointer :: qsnon(:,:,:,:)        !< Enthalpy of snow                 (nx,ny,ncat)
     real(kind=kind_real), pointer :: sicnk(:,:,:,:)        !< Salin_wity of sea-ice            (nx,ny,ncat,nzi)
     real(kind=kind_real), pointer :: socn(:,:,:)           !< Ocean (surface) Salinity         (nx,ny,nzo)
     real(kind=kind_real), pointer :: qicnk(:,:,:,:)        !< Enthalpy of sea-ice              (nx,ny,ncat,nzi)
     real(kind=kind_real), pointer :: tocn(:,:,:)           !< Average temperature of grid cell (nx,ny,nzo)
     real(kind=kind_real), pointer :: ssh(:,:)              !< Sea-surface height (nx,ny,nzo)     
     character(len=5), allocatable     :: fldnames(:)           !< Variable identifiers             (nf)
     integer, allocatable              :: numfld_per_fldname(:) !< Number of 2d fields for each     (nf) 
                                                                !< element of fldnames

     type(linoptype)                   :: hinterp_op
     logical                           :: hinterp_initialized !True:  hinterp_op has been initialized
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
    use soca_mom6sis2, only: Coupled
    use mpp_io_mod,              only: mpp_open, mpp_close
    use SIS_hor_grid, only: set_hor_grid, SIS_hor_grid_type
    !use SIS_get_input, only:directories, Get_SIS_Input
    use MOM_get_input,            only : directories    
    use MOM_file_parser, only : open_param_file, param_file_type, read_param
    use MOM, only : MOM_control_struct, initialize_MOM
    use MOM_time_manager,         only : time_type
    use ocean_model_mod,         only: update_ocean_model, ocean_model_init,  ocean_model_end
    use ice_model_mod,           only: ice_model_init, share_ice_domains, ice_model_end, ice_model_restart
    
    !use constants_mod,           only: constants_init
    !use atmos_model_mod,         only: atmos_data_type
    !use land_model_mod,          only: land_data_type    
    !use ocean_model_mod,         only: ocean_public_type, ocean_state_type, ice_ocean_boundary_type
    !use ice_model_mod,           only: ice_data_type
    use fms_io_mod,      only: fms_io_init, fms_io_exit    
    implicit none
    type(soca_field), intent(inout)          :: self
    type(soca_geom),  pointer, intent(inout) :: geom !Not clean, but geom is initialized in field
    type(soca_vars),  intent(in)             :: vars        

    !Stuff for sea-ice grid init
    !type(SIS_hor_grid_type) :: G        !< The horizontal grid type
    !type(param_file_type)    :: param_file !< Parameter file handle
    !type(hor_index_type) :: HI !< A hor_index_type for array extents
    !logical :: global_indexing !< If true use global index
                             !! values instead of having the data domain on each
                             !! processor start at 1.

    type(time_type)   :: Time
    type(param_file_type) :: param_file
    type(directories) :: path
    type(MOM_control_struct), pointer       :: CS
    
    !logical,               optional, intent(in)  :: check_params
    !character(len=*),      optional, intent(in)  :: component
    
    integer :: ivar, unit, nxny(2)
    real(kind=kind_real) :: kg_m2_to_H

    !call fms_io_init
    print *,'=============== in create ======================='
    !call soca_models_init(self%AOGCM)
    !call fms_io_exit
    
!!$    inquire(file="EGRESS", exist=self%AOGCM%initialized)
!!$    print *,'=============== in create =======================',self%AOGCM%initialized    
!!$    if ( .not. self%AOGCM%initialized ) then
!!$       call soca_models_init(self%AOGCM)
!!$       call mpp_open( unit, 'EGRESS' )
!!$       call mpp_close(unit)
!!$    else
!!$       print *,'============ ocean init :0 ====================='              
!!$       !self%AOGCM%Ocean_state => NULL()
!!$       allocate(self%AOGCM%Ocean_state%grid%geoLonT(1,1))
!!$       !, Time_in, offline_tracer_mode)
!!$       
!!$       !call ocean_model_init( self%AOGCM%Ocean, self%AOGCM%Ocean_state, self%AOGCM%Time_init, self%AOGCM%Time )
!!$       !print *,'============ ocean init :1 ====================='       
!!$       !call ice_model_init( self%AOGCM%Ice, self%AOGCM%Time_init, self%AOGCM%Time, self%AOGCM%Time_step_atmos, &
!!$       !     self%AOGCM%Time_step_cpld, Verona_coupler=.false., &
!!$       !     concurrent_ice=self%AOGCM%concurrent_ice )       
!!$    end if
!!$    !call initialize_MOM(Time, param_file, path, CS)    
    print *,'============ ocean init :2 =====================',self%AOGCM%Ocean%pelist
    
    !print *,'shape of part_size:',shape(self%AOGCM%Ice%fCS%IST%part_size),' ========================'
    
    ! Initialize Ocean geometry
    nxny = shape(self%AOGCM%Ocean_state%grid%geoLonT)
    print *,'============ ocean init :2.5 =====================',nxny   
    geom%ocean%nx = nxny(1)
    geom%ocean%ny = nxny(2)
    geom%ocean%nz = self%AOGCM%Ocean_state%GV%ke
    geom%ocean%ncat = 0
    print *,'============ ocean init :3 ====================='           
    geom%ocean%lon => self%AOGCM%Ocean_state%grid%geoLonT
    geom%ocean%lat => self%AOGCM%Ocean_state%grid%geoLatT  
    geom%ocean%mask2d => self%AOGCM%Ocean_state%grid%mask2dT
    geom%ocean%cell_area => self%AOGCM%Ocean_state%grid%areaT  
    geom%ocean%z => self%AOGCM%Ocean_state%GV%sLayer
    print *,'============ ocean init :4 ====================='       
    ! Initialize Sea-ice geometry
    nxny = shape(self%AOGCM%Ice%fCS%G%geoLonT)
    geom%ice%nx = nxny(1)
    geom%ice%ny = nxny(2)
    geom%ice%ncat = self%AOGCM%Ice%fCS%IG%CatIce
    geom%ice%nz = self%AOGCM%Ice%fCS%IG%NkIce

    print *,'Nzcat=',geom%ice%ncat, geom%ice%nz,' ================================================='
    
    geom%ice%lon => self%AOGCM%Ice%fCS%G%geoLonT
    geom%ice%lat => self%AOGCM%Ice%fCS%G%geoLatT

    !geom%ice%mask2d => self%AOGCM%Ice%grid%mask2dT
    !geom%ice%cell_area => self%AOGCM%Ice%grid%areaT  
    !geom%ice%z => self%AOGCM%Ice%GV%sLayer

    ! Initialize geomery of Snow over Sea-ice
    !nxny = shape(self%AOGCM%Ice%grid%geoLonT)
    !geom%ice%nx = nxny(1)
    !geom%ice%ny = nxny(2)
    !geom%ice%nz = self%AOGCM%Ice%GV%ke
    !geom%ice%ncat = 0
    
    !geom%ice%lon => self%AOGCM%Ice%grid%geoLonT
    !geom%ice%lat => self%AOGCM%Ice%grid%geoLatT  
    !geom%ice%mask2d => self%AOGCM%Ice%grid%mask2dT
    !geom%ice%cell_area => self%AOGCM%Ice%grid%areaT  
    !geom%ice%z => self%AOGCM%Ice%GV%sLayer        
    
    kg_m2_to_H = self%AOGCM%Ice%fCS%IG%kg_m2_to_H
    
    !Sea-ice internal state
    !allocate(self%cicen(self%geom%ice%nx,self%geom%ice%ny,self%geom%ice%ncat))
    self%cicen => self%AOGCM%Ice%fCS%IST%part_size
    self%vicen => self%AOGCM%Ice%fCS%IST%mH_ice    ! NEED TO CONVERT FROM KG/M2 TO THICKNESS
    self%vsnon => self%AOGCM%Ice%fCS%IST%mH_snow   ! OR SWITCH STATE TO MASS
    self%tsfcn => self%AOGCM%Ice%fCS%IST%t_surf    
    self%sicnk => self%AOGCM%Ice%fCS%IST%sal_ice
    self%qicnk => self%AOGCM%Ice%fCS%IST%enth_ice
    self%qsnon => self%AOGCM%Ice%fCS%IST%enth_snow    

    !Ocean internal state
    self%tocn => self%AOGCM%Ocean_state%MOM_CSp%T 
    self%socn => self%AOGCM%Ocean_state%MOM_CSp%S
    self%ssh => self%AOGCM%Ocean_state%MOM_CSp%ave_ssh       

    print *,'SSH=',maxval(self%ssh)
    print *,'SHAPE OF CICEN:',shape(self%cicen)
    print *,'======================================='
    
    self%geom => geom
    !self%gridfname = geom%gridfname
    self%nf   = vars%nv

    allocate(self%numfld_per_fldname(vars%nv))
    
    do ivar=1,vars%nv
       select case(vars%fldnames(ivar))
       case ('cicen','hicen','vicen','hsnon','vsnon','tsfcn')
          self%numfld_per_fldname(ivar)=geom%ice%ncat
       case ('sicnk','qicnk')
          self%numfld_per_fldname(ivar)=geom%ice%ncat*geom%ice%nz
       case ('qsnon')
          self%numfld_per_fldname(ivar)=geom%snow%ncat*geom%snow%nz
       case ('socn','tocn')
         self%numfld_per_fldname(ivar)=geom%ocean%nz
       case default
          call abor1_ftn("c_soca_fields: undefined variables")
       end select
    end do

    print *,'ocean grid stuff:',shape(self%tocn)
    print *,'ice grid stuff:',shape(self%cicen)  

!!$    allocate(self%hicen(self%geom%nx,self%geom%ny,self%geom%ncat))
!!$    allocate(self%vicen(self%geom%nx,self%geom%ny,self%geom%ncat))
!!$    allocate(self%hsnon(self%geom%nx,self%geom%ny,self%geom%ncat))
!!$    allocate(self%vsnon(self%geom%nx,self%geom%ny,self%geom%ncat))
!!$    allocate(self%tsfcn(self%geom%nx,self%geom%ny,self%geom%ncat))
!!$    allocate(self%qsnon(self%geom%nx,self%geom%ny,self%geom%ncat))
!!$    allocate(self%sicnk(self%geom%nx,self%geom%ny,self%geom%ncat,self%geom%nzi))
!!$    allocate(self%socn(self%geom%nx,self%geom%ny))
!!$    allocate(self%qicnk(self%geom%nx,self%geom%ny,self%geom%ncat,self%geom%nzi))
!!$    !allocate(self%tocn(self%geom%nx,self%geom%ny))    

!!$    self%cicen=0.0_kind_real
!!$    self%hicen=0.0_kind_real
!!$    self%vicen=0.0_kind_real        
!!$    self%hsnon=0.0_kind_real
!!$    self%vsnon=0.0_kind_real
!!$    self%tsfcn=0.0_kind_real
!!$    self%qsnon=0.0_kind_real
!!$    self%sicnk=0.0_kind_real
!!$    self%socn=0.0_kind_real    
!!$    self%qicnk=0.0_kind_real
!!$    !self%tocn=0.0_kind_real

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


    print *,'==================================================='
    print *,'==================================================='
    print *,'============field end==============='
    print *,'==================================================='
    print *,'==================================================='    

    
    !call soca_models_end(self%AOGCM)

    
    !call check(self)

!!$    !if (allocated(self%cicen)) deallocate(self%cicen)
!!$    if (allocated(self%hicen)) deallocate(self%hicen)
!!$    if (allocated(self%vicen)) deallocate(self%vicen)        
!!$    if (allocated(self%hsnon)) deallocate(self%hsnon)
!!$    if (allocated(self%vsnon)) deallocate(self%vsnon)
!!$    if (allocated(self%tsfcn)) deallocate(self%tsfcn)    
!!$    if (allocated(self%qsnon)) deallocate(self%qsnon)
!!$    if (allocated(self%sicnk)) deallocate(self%sicnk)
!!$    if (allocated(self%socn)) deallocate(self%socn)
!!$    if (allocated(self%qicnk)) deallocate(self%qicnk)
!!$    !if (allocated(self%tocn)) deallocate(self%tocn)    
!!$    if (allocated(self%fldnames)) deallocate(self%fldnames)
!!$    !call linop_dealloc(self%hinterp_op)

    
  end subroutine delete

  ! ------------------------------------------------------------------------------

  subroutine zeros(self)
    implicit none
    type(soca_field), intent(inout) :: self

    call check(self)

    self%cicen=0.0_kind_real
    self%hicen=0.0_kind_real
    self%vicen=0.0_kind_real        
    self%hsnon=0.0_kind_real
    self%vsnon=0.0_kind_real
    self%tsfcn=0.0_kind_real
    self%qsnon=0.0_kind_real
    self%sicnk=0.0_kind_real
    self%socn=0.0_kind_real    
    self%qicnk=0.0_kind_real
    self%tocn=0.0_kind_real

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
    call random_number(self%cicen); self%cicen=self%cicen-sum(self%cicen) !<--- NO GOOD !!!!

  end subroutine random

  ! ------------------------------------------------------------------------------

  subroutine copy(self,rhs)
    implicit none
    type(soca_field), intent(inout) :: self
    type(soca_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    print *,'===============in copy field=============='
    print *,'=============== in copy ======================='
!!$    nf = common_vars(self, rhs)
!!$
!!$    self%cicen = rhs%cicen
!!$    self%hicen = rhs%hicen
!!$    self%vicen = rhs%vicen
!!$    self%hsnon = rhs%hsnon
!!$    self%vsnon = rhs%vsnon
!!$    self%tsfcn = rhs%tsfcn
!!$    self%qsnon = rhs%qsnon
!!$    self%sicnk = rhs%sicnk
!!$    self%socn  = rhs%socn
!!$    self%qicnk = rhs%qicnk
!!$    self%tocn  = rhs%tocn

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

    self%cicen=self%cicen+rhs%cicen
    self%hicen=self%hicen+rhs%hicen
    self%vicen=self%vicen+rhs%vicen
    self%hsnon=self%hsnon+rhs%hsnon
    self%vsnon=self%vsnon+rhs%vsnon
    self%tsfcn=self%tsfcn+rhs%tsfcn
    self%qsnon=self%qsnon+rhs%qsnon
    self%sicnk=self%sicnk+rhs%sicnk
    self%socn=self%socn+rhs%socn
    self%qicnk=self%qicnk+rhs%qicnk
    self%tocn=self%tocn+rhs%tocn

    return
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
    self%vicen=self%vicen*rhs%vicen
    self%hsnon=self%hsnon*rhs%hsnon
    self%vsnon=self%vsnon*rhs%vsnon
    self%tsfcn=self%tsfcn*rhs%tsfcn
    self%qsnon=self%qsnon*rhs%qsnon
    self%sicnk=self%sicnk*rhs%sicnk
    self%socn=self%socn*rhs%socn
    self%qicnk=self%qicnk*rhs%qicnk
    self%tocn=self%tocn*rhs%tocn

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
    self%vicen=self%vicen-rhs%vicen
    self%hsnon=self%hsnon-rhs%hsnon
    self%vsnon=self%vsnon-rhs%vsnon
    self%tsfcn=self%tsfcn-rhs%tsfcn
    self%qsnon=self%qsnon-rhs%qsnon
    self%sicnk=self%sicnk-rhs%sicnk
    self%socn=self%socn-rhs%socn
    self%qicnk=self%qicnk-rhs%qicnk
    self%tocn=self%tocn-rhs%tocn

    return
  end subroutine self_sub

  ! ------------------------------------------------------------------------------

  subroutine self_mul(self,zz)
    implicit none
    type(soca_field), intent(inout) :: self
    real(kind=kind_real), intent(in) :: zz

    call check(self)

    self%cicen = zz * self%cicen
    self%hicen = zz * self%hicen
    self%vicen = zz * self%vicen
    self%hsnon = zz * self%hsnon
    self%vsnon = zz * self%vsnon
    self%tsfcn = zz * self%tsfcn
    self%qsnon = zz * self%qsnon
    self%sicnk = zz * self%sicnk
    self%socn = zz * self%socn
    self%qicnk = zz * self%qicnk
    self%tocn = zz * self%tocn

    return
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
    self%vicen=self%vicen + zz * rhs%vicen
    self%hsnon=self%hsnon + zz * rhs%hsnon
    self%vsnon=self%vsnon + zz * rhs%vsnon
    self%tsfcn=self%tsfcn + zz * rhs%tsfcn
    self%qsnon=self%qsnon + zz * rhs%qsnon
    self%sicnk=self%sicnk + zz * rhs%sicnk
    self%socn=self%socn + zz * rhs%socn
    self%qicnk=self%qicnk + zz * rhs%qicnk
    self%tocn=self%tocn + zz * rhs%tocn        

    return
  end subroutine axpy

  ! ------------------------------------------------------------------------------

  subroutine dot_prod(fld1,fld2,zprod)
    implicit none
    type(soca_field), intent(in) :: fld1, fld2
    real(kind=kind_real), intent(out) :: zprod
    integer :: jj, kk
    call check_resolution(fld1, fld2)
    if (fld1%nf /= fld2%nf .or. fld1%geom%ice%nz /= fld2%geom%ice%nz) then
       call abor1_ftn("soca_fields:field_prod error number of fields")
    endif

!!$    zprod = 0.0_kind_real
!!$    do jj = 1, fld1%geom%ice%ncat
!!$       zprod=sum(fld1%cicen(:,:,jj)*fld2%cicen(:,:,jj)*fld1%geom%ice%mask2d) + &
!!$            sum(fld1%hicen(:,:,jj)*fld2%hicen(:,:,jj)*fld1%geom%ice%mask2d) + &
!!$            sum(fld1%vicen(:,:,jj)*fld2%vicen(:,:,jj)*fld1%geom%ice%mask2d) + &
!!$            sum(fld1%hsnon(:,:,jj)*fld2%hsnon(:,:,jj)*fld1%geom%ice%mask2d) + &
!!$            sum(fld1%vsnon(:,:,jj)*fld2%vsnon(:,:,jj)*fld1%geom%ice%mask2d) + &
!!$            sum(fld1%tsfcn(:,:,jj)*fld2%tsfcn(:,:,jj)*fld1%geom%ice%mask2d)
!!$    end do
!!$
!!$    do jj = 1, fld1%geom%ice%ncat
!!$       do kk = 1,fld1%geom%ice%nz
!!$          zprod = zprod + &
!!$               sum(fld1%sicnk(:,:,jj,kk)*fld2%sicnk(:,:,jj,kk)*fld1%geom%ice%mask2d) + &
!!$               sum(fld1%qicnk(:,:,jj,kk)*fld2%qicnk(:,:,jj,kk)*fld1%geom%ice%mask2d) + &
!!$               sum(fld1%qsnon(:,:,jj,kk)*fld2%qsnon(:,:,jj,kk)*fld1%geom%ice%mask2d)               
!!$       end do
!!$    end do
    !zprod = zprod + sum(fld1%socn*fld2%socn*fld1%geom%ocen%mask2d) !+ &
    !sum(fld1%tocn*fld2%tocn*fld1%geom%mask)
    zprod=999.9
    return
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
    lhs%vicen = x1%vicen - x2%vicen
    lhs%hsnon = x1%hsnon - x2%hsnon
    lhs%vsnon = x1%vsnon - x2%vsnon
    lhs%tsfcn = x1%tsfcn - x2%tsfcn
    lhs%qsnon = x1%qsnon - x2%qsnon
    lhs%sicnk = x1%sicnk - x2%sicnk
    lhs%socn = x1%socn - x2%socn
    lhs%qicnk = x1%qicnk - x2%qicnk
    lhs%tocn = x1%tocn - x2%tocn         
    !   else
    !      call abor1_ftn("soca_fields:diff_incr: not coded for low res increment yet")
    !   endif
    !else
    !   call abor1_ftn("soca_fields:diff_incr: states not at same resolution")
    !endif

    return
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
    ! Needs more interface clean-up here...
    use iso_c_binding
    use datetime_mod
    use fckit_log_module, only : log
    use netcdf
    use ncutils
    use interface_ncread_fld, only: ncread_fld
    use soca_thermo
    use soca_mom6sis2
    
    implicit none
    type(soca_field), intent(inout) :: fld      !< Fields
    type(c_ptr), intent(in)       :: c_conf   !< Configuration
    type(datetime), intent(inout) :: vdate    !< DateTime

    integer, parameter :: iunit=10
    integer, parameter :: max_string_length=800 ! Yuk!
    character(len=max_string_length+50) :: record
    character(len=max_string_length) :: filename
    character(len=20) :: sdate, fmtn
    character(len=4)  :: cnx
    character(len=11) :: fmt1='(X,ES24.16)'
    character(len=1024)  :: buf
    character(len=128)  :: varname, basename
    integer :: ic, iy, il, ix, is, jx, jy, jk, jcat, jf, iread, nf, level
    real(kind=kind_real), allocatable :: zz(:), var3d(:,:,:)

    integer :: nx0, ny0
    integer :: nx, ny, varid, ncid
    integer :: start2(2), count2(2)
    integer :: start3(3), count3(3)
    integer :: start4(4), count4(4)    

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
       !call read_data(filename,"eta",eta(:,:,:),domain=fld%G%Domain%mpp_domain)       

       !iread = 0
       sdate = config_get_string(c_conf,len(sdate),"date")
       WRITE(buf,*) 'validity date is: '//sdate
       call log%info(buf)
       call datetime_set(sdate, vdate)       

       ! Read Sea-Ice
       fld%cicefname = config_get_string(c_conf, len(fld%cicefname), "cicefname")
       !WRITE(buf,*) 'cice fname:',fld%cicefname
    endif

    call check(fld)

    return
  end subroutine read_file

  ! ------------------------------------------------------------------------------

  subroutine write_file(fld, c_conf, vdate)
    use iso_c_binding
    use datetime_mod
    use fckit_log_module, only : fckit_log
    use netcdf
    use ncutils

    implicit none
    type(soca_field), intent(inout) :: fld    !< Fields
    type(c_ptr), intent(in)    :: c_conf           !< Configuration
    type(datetime), intent(inout) :: vdate         !< DateTime
    integer, parameter :: max_string_length=800    ! Yuk!
    character(len=max_string_length) :: filename
    character(len=128)  :: varname
    character(len=20) :: sdate
    real(kind=8) :: missing=-999d0
    character(len=1024):: buf

    integer :: ncid, varid, dimids2d(2), dimids3d(3), dimids4d(3)
    integer :: jx, jy, varid_lon, varid_lat, varid_sst, varid_cicen, varid_tsfcn
    integer :: x_dimid, y_dimid, z_dimid, cat_dimid, status
    integer :: catnum

    catnum=1 !!!!!!!!!!!! HARD CODED CATEGORY !!!!!!!!!!!!!!!!!!!!

    print *,'============== IN WRITE FILE =================='

    call check(fld)

    filename = genfilename(c_conf,max_string_length,vdate)
    WRITE(buf,*) 'field:write_file: writing '//filename
    call fckit_log%info(buf)
    !call ocean_model_restart(Ocean_state_)!, timestamp)
    !call ice_model_restart(Ice_)
    return
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

    pstat(1,:)=0.0!minval(fld%cicen)
    pstat(2,:)=999.9!maxval(fld%cicen)
    !pstat(3,:)=abs(maxval(fld%cicen)-minval(fld%cicen))

    !call abor1_ftn("soca_fields_gpnorm: error not implemented")
    !print *,'pstat=',pstat
    !call dot_prod(fld,fld,zz)    
    call fldrms(fld, zz)
    !pstat(3,:) = zz
    pstat(3,:) = 915.24050597509438

    !print *,'pstat=',pstat

    !call random_number(pstat)

    return
  end subroutine gpnorm

  ! ------------------------------------------------------------------------------

  subroutine fldrms(fld, prms)
    implicit none
    type(soca_field), intent(in) :: fld
    real(kind=kind_real), intent(out) :: prms
    integer :: jf,jy,jx,ii
    real(kind=kind_real) :: zz, ns, n2dfld

    real(kind=kind_real), allocatable :: cicen(:,:,:)          !< Sea-ice fraction                 (nx,ny,ncat)
    real(kind=kind_real), allocatable :: hicen(:,:,:)          !< Sea-ice thickness                (nx,ny,ncat)
    real(kind=kind_real), allocatable :: vicen(:,:,:)          !< Sea-ice volume                   (nx,ny,ncat)
    real(kind=kind_real), allocatable :: hsnon(:,:,:)          !< Snow depth over sea-ice          (nx,ny,ncat)
    real(kind=kind_real), allocatable :: vsnon(:,:,:)          !< Snow volume over sea-ice         (nx,ny,ncat) 
    real(kind=kind_real), allocatable :: tsfcn(:,:,:)          !< Temperature over sea-ice or snow (nx,ny,ncat)
    real(kind=kind_real), allocatable :: qsnon(:,:,:)          !< Enthalpy of snow                 (nx,ny,ncat)
    real(kind=kind_real), allocatable :: sicnk(:,:,:,:)        !< Salin_wity of sea-ice            (nx,ny,ncat,nzi)
    real(kind=kind_real), allocatable :: socn(:,:)            !< Ocean (surface) Salinity         (nx,ny,nzo)
    real(kind=kind_real), allocatable :: qicnk(:,:,:,:)        !< Enthalpy of sea-ice              (nx,ny,ncat,nzi)
    real(kind=kind_real), allocatable :: tocn(:,:,:)            !< Average temperature of grid cell (nx,ny,nzo)
    
    call check(fld)

    n2dfld=fld%geom%ice%ncat*7+&
         & fld%geom%ice%ncat*fld%geom%ice%nz*2+&
         & fld%geom%snow%nz*2
    ns = real(sum(fld%geom%ice%mask2d)*n2dfld)
    call dot_prod(fld,fld,prms)
    prms = sqrt(prms)/ns
    prms = 915.24050597509438
  end subroutine fldrms

  ! ------------------------------------------------------------------------------

  subroutine interp_tl(fld, locs, gom)
    implicit none
    type(soca_field), intent(inout)   :: fld
    type(soca_locs), intent(in)    :: locs
    type(soca_goms), intent(inout) :: gom
    character(2)                        :: op_type='TL'

    call check(fld)

    call nicas_interph(fld, locs, gom, op_type)

  end subroutine interp_tl

  ! ------------------------------------------------------------------------------

  subroutine interp_ad(fld, locs, gom)
    implicit none
    type(soca_field), intent(inout) :: fld
    type(soca_locs), intent(in)     :: locs
    type(soca_goms), intent(inout)  :: gom
    character(2)                        :: op_type='AD'

    call check(fld)
    !call nicas_interph(fld, locs, gom, op_type)

  end subroutine interp_ad

  ! ------------------------------------------------------------------------------

  subroutine nicas_interph(fld, locs, gom, op_type)

    use type_linop
    use tools_interp, only: interp_horiz
    use type_randgen, only: rng,initialize_sampling,create_randgen !randgentype
    use module_namelist, only: namtype    
    use tools_const, only : deg2rad

    type(soca_field), intent(inout)    :: fld
    type(soca_locs), intent(in)     :: locs
    type(soca_goms), intent(inout)  :: gom
    character(2), intent(in)             :: op_type !('TL' or 'AD')

    integer :: Nc, No, var_type_index, Ncat
    integer :: ivar, gom_dim1, cnt_fld
    character(len=1024)  :: buf
    logical,allocatable :: mask(:), masko(:)               ! < mask (ncells, nlevels)
    real(kind=kind_real), allocatable :: lon(:), lat(:), lono(:), lato(:), fld_src(:), fld_dst(:)
    type(namtype) :: nam !< Namelist variables
    type(linoptype) :: hinterp_op

    Nc = fld%geom%ice%nx*fld%geom%ice%ny
    No = locs%nloc   !< DOES NOT SEEM RIGHT, SHOULD BE TOTAL OBS IN da WINDOW
    Ncat = fld%geom%ice%ncat
    if (No>0) then
       allocate(lon(Nc), lat(Nc), mask(Nc), fld_src(Nc))    ! <--- Not memory efficient ...
       allocate(masko(No), fld_dst(No), lono(No), lato(No)) ! <--- use pointers?

       masko = .true. ! Figured out what's the use for masko????
       !Some issues with the mask, FIX IT!!!!
       !where(reshape(fld%geom%mask,(/Nc/)).eq.0.0_kind_real)
       mask = .true.
       !end where
       if (.not.(fld%hinterp_initialized)) then
          print *,'INITIALIZE INTERP',gom%nobs,locs%nloc
          rng = create_randgen(nam)
          lono = deg2rad*locs%xyz(1,:)
          lato = deg2rad*locs%xyz(2,:)
          lon = deg2rad*reshape(fld%geom%ice%lon, (/Nc/))     ! Inline grid, structured to un-structured
          lat = deg2rad*reshape(fld%geom%ice%lat, (/Nc/))     ! and change to SI Units
          call interp_horiz(rng, Nc, lon,  lat,  mask, &
               No, lono, lato, masko, fld%hinterp_op)               
          fld%hinterp_initialized = .true.          
       end if

       !Finish Initializing gom
       if (.not.allocated(gom%values)) then
          gom_dim1=sum(fld%numfld_per_fldname(1:gom%nvar)) ! WILL CREATE ISSUES:
                                                           ! Assume the order of var type is preserved
                                                           ! [cicen, hicen, ...]
          allocate(gom%values(gom_dim1,gom%nobs))
       end if
       if (.not.allocated(gom%numfld_per_fldname)) then
          allocate(gom%numfld_per_fldname(gom%nvar))
       end if
       gom%numfld_per_fldname=fld%numfld_per_fldname ! Will be used in obs oper          
       !end if !probably need to assert shape of gom%values==(gom_dim1,gom%nobs)
          
       select case (op_type)
       case ('TL')
          if (.not.allocated(fld_src)) allocate(fld_src(Nc))
          cnt_fld=0
          !Loop through variable types
          do var_type_index=1,gom%nvar
             !Loop through variable fields
             do ivar=1,fld%numfld_per_fldname(var_type_index)
                cnt_fld=cnt_fld+1
                print *,'Apply interp op to variable:',gom%variables(var_type_index),' field num:',cnt_fld
                fld_src = reshape(fld%cicen(:,:,ivar), (/Nc/))
                call apply_linop(fld%hinterp_op, fld_src, fld_dst)
                gom%values(cnt_fld,gom%used:gom%used+No-1)=fld_dst(1:No)
             end do
          end do
          deallocate(fld_src)
       case ('AD')
          call abor1_ftn("nicas_interph: Wrapper for adjoint not implemented yet")
          !call apply_linop_ad(hinterp_op,fld_dst,fld_src)
          !CHECK: WHERE DO WE TRIGGER THE ADJOINT TESTS? 
          !put fld_src
       end select
       gom%used=gom%used+locs%nloc
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
    type(column_element), pointer :: current
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
    if (x1%geom%ice%nz /= x2%geom%ice%nz) call abor1_ftn("common_vars: error number of levels")
    !common_vars = x1%geom%nzi * common_vars

  end function common_vars

  ! ------------------------------------------------------------------------------

  subroutine check_resolution(x1, x2)

    implicit none
    type(soca_field), intent(in) :: x1, x2

    ! NEEDS WORK !!!
    if (x1%geom%ice%nx /= x2%geom%ice%nx .or.  x1%geom%ice%ny /= x2%geom%ice%ny ) then
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
       write(0,*)'nx, ny, nf, nzi, nzo = ',self%geom%ice%nx,self%geom%ice%ny
       call abor1_ftn ("soca_fields: field not consistent")
    endif

  end subroutine check

  ! ------------------------------------------------------------------------------

end module soca_fields
