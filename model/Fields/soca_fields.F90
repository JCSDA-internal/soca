! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Handle fields for the  model

module soca_fields

  use config_mod
  use soca_geom_mod
  use soca_vars_mod
  use soca_interph_mod
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
       & create, delete, zeros, ones, dirac, random, copy, create_copy,&
       & self_add, self_schur, self_sub, self_mul, axpy, &
       & dot_prod, add_incr, diff_incr, &
       & read_file, write_file, gpnorm, fldrms, &
       & change_resol, interp_tl, interp_ad, convert_to_ug, convert_from_ug
  public :: soca_field_registry

  interface create
     module procedure create_constructor, create_copy
  end interface create
           
  ! ------------------------------------------------------------------------------
  !> Fortran derived type to hold fields
  type :: soca_field
     type (Coupled)                    :: AOGCM
     type(soca_geom), pointer          :: geom                  !< MOM6 & SIS2 Geometry
     integer                           :: nf                    !< Number of fields
     character(len=128)                :: gridfname             !< Grid file name
     character(len=128)                :: cicefname             !< Fields file name for cice
     character(len=128)                :: momfname              !< Fields file name for mom
     real(kind=kind_real), pointer     :: cicen(:,:,:) => NULL()   !< Sea-ice fraction                 (nx,ny,ncat+1)
     real(kind=kind_real), pointer     :: hicen(:,:,:) => NULL()   !< Sea-ice mass/m2                  (nx,ny,ncat) [kg/m2]
     real(kind=kind_real), pointer     :: vicen(:,:,:) => NULL()   !< Sea-ice volume                   (nx,ny,ncat)
     real(kind=kind_real), pointer     :: hsnon(:,:,:) => NULL()   !< Snow depth over sea-ice          (nx,ny,ncat)
     real(kind=kind_real), pointer     :: vsnon(:,:,:) => NULL()   !< Snow volume over sea-ice         (nx,ny,ncat) 
     real(kind=kind_real), pointer     :: tsfcn(:,:,:) => NULL()   !< Temperature over sea-ice or snow (nx,ny,ncat)
     real(kind=kind_real), pointer     :: qsnon(:,:,:,:) => NULL() !< Enthalpy of snow                 (nx,ny,ncat)
     real(kind=kind_real), pointer     :: sicnk(:,:,:,:) => NULL() !< Salin_wity of sea-ice            (nx,ny,ncat,nzi)
     real(kind=kind_real), pointer     :: socn(:,:,:) => NULL()    !< Ocean Practical Salinity         (nx,ny,nzo)
     real(kind=kind_real), pointer     :: qicnk(:,:,:,:) => NULL() !< Enthalpy of sea-ice              (nx,ny,ncat,nzi)
     real(kind=kind_real), pointer     :: tocn(:,:,:) => NULL()    !< Ocean Potential Temperature, ref to p=0      (nx,ny,nzo)
     real(kind=kind_real), pointer     :: ssh(:,:) => NULL()       !< Sea-surface height (nx,ny,nzo)
     real(kind=kind_real), pointer     :: hocn(:,:,:) => NULL()    !< Layer thickness (nx,ny,nzo)     
     character(len=5), allocatable     :: fldnames(:)              !< Variable identifiers             (nf)
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

  subroutine create_constructor(self, geom, vars)
    ! Construct a field from geom and vars
    use soca_mom6sis2, only: Coupled, soca_field_init, soca_geom_init
    use soca_vars_mod
    
    implicit none
    type(soca_field), intent(inout)          :: self
    type(soca_geom),  pointer, intent(inout) :: geom
    type(soca_vars),  intent(in)             :: vars        
    integer :: ivar!, unit, nxny(2)
    
    self%geom => geom
    self%nf   = vars%nv
    call soca_field_init(self%AOGCM, geom%ocean%G, geom%ocean%GV, geom%ocean%IG)

    ! Assign convenience pointers
    !Ocean internal state    
    self%tocn => self%AOGCM%Ocn%T
    self%socn => self%AOGCM%Ocn%S
    self%hocn => self%AOGCM%Ocn%H
    self%ssh => self%AOGCM%Ocn%ssh
    self%ssh = self%ssh*self%geom%ocean%mask2d

    !Sea-ice internal state
    self%cicen => self%AOGCM%Ice%part_size
    self%hicen => self%AOGCM%Ice%h_ice
    self%hsnon => self%AOGCM%Ice%h_snow        
    self%tsfcn => self%AOGCM%Ice%T_skin
    self%sicnk => self%AOGCM%Ice%sal_ice
    self%qicnk => self%AOGCM%Ice%enth_ice
    self%qsnon => self%AOGCM%Ice%enth_snow

    call zeros(self)

    if (self%nf>11) then
       print *,'Number of fields:',self%nf
       call abor1_ftn ("soca_fields:create error number of fields")       
    endif
    allocate(self%fldnames(self%nf))
    self%fldnames(:)=vars%fldnames(:)

    call check(self)

  end subroutine create_constructor

  ! ------------------------------------------------------------------------------
     
  subroutine create_copy(self, rhs_fld)
    ! Construct a field from an other field, lhs_fld=rhs_fld       
    use soca_mom6sis2, only: Coupled, soca_field_init, soca_geom_init
    use soca_vars_mod
    
    implicit none
    type(soca_field), intent(inout)          :: self
    type(soca_field), intent(inout)          :: rhs_fld
    integer :: ivar!, unit, nxny(2)
    
    self%geom => rhs_fld%geom
    self%nf   = rhs_fld%nf
    call soca_field_init(self%AOGCM, rhs_fld%geom%ocean%G, rhs_fld%geom%ocean%GV, rhs_fld%geom%ocean%IG)

    ! Assign convenience pointers
    !Ocean internal state    
    self%tocn => self%AOGCM%Ocn%T
    self%socn => self%AOGCM%Ocn%S
    self%ssh => self%AOGCM%Ocn%ssh
    self%hocn => self%AOGCM%Ocn%H
    
    !Sea-ice internal state
    self%cicen => self%AOGCM%Ice%part_size
    self%hicen => self%AOGCM%Ice%h_ice
    self%hsnon => self%AOGCM%Ice%h_snow        
    self%tsfcn => self%AOGCM%Ice%T_skin
    self%sicnk => self%AOGCM%Ice%sal_ice
    self%qicnk => self%AOGCM%Ice%enth_ice
    self%qsnon => self%AOGCM%Ice%enth_snow

    call zeros(self)

    if (self%nf>11) then
       call abor1_ftn ("soca_fields:create_copy error number of fields")       
    endif
    allocate(self%fldnames(self%nf))
    self%fldnames(:)=rhs_fld%fldnames(:)

    call check(self)

  end subroutine create_copy

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
    self%hocn = 0.0_kind_real    
  end subroutine zeros

  ! ------------------------------------------------------------------------------

  subroutine ones(self)
    implicit none
    type(soca_field), intent(inout) :: self

    call check(self)

    self%cicen = 1.0_kind_real
    self%hicen = 1.0_kind_real
    self%hsnon = 1.0_kind_real
    self%tsfcn = 1.0_kind_real
    self%qsnon = 1.0_kind_real
    self%sicnk = 1.0_kind_real
    self%qicnk = 1.0_kind_real

    self%socn = 1.0_kind_real        
    self%tocn = 1.0_kind_real
    self%ssh = 1.0_kind_real
    self%tocn = 1.0_kind_real
    
  end subroutine ones

  ! ------------------------------------------------------------------------------  

  subroutine dirac(self, c_conf)
    use iso_c_binding
    implicit none
    type(soca_field), intent(inout) :: self
    type(c_ptr), intent(in)       :: c_conf   !< Configuration
    integer :: ndir,idir,ildir,ifdir
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
    ildir = config_get_int(c_conf,"ildir")
    ifdir = config_get_int(c_conf,"ifdir")


    
    ! Check 
    if (ndir<1) call abor1_ftn("fields:dirac non-positive ndir")
    if (any(ixdir<1).or.any(ixdir>self%geom%ocean%nx)) call abor1_ftn("fields:dirac invalid ixdir")
    if (any(iydir<1).or.any(iydir>self%geom%ocean%ny)) call abor1_ftn("fields:dirac invalid iydir")

    ! Setup Diracs
    call zeros(self)
    !ioff = (ifdir-1)*self%nl
    do idir=1,ndir
       !self%qicnk(ixdir(idir),iydir(idir),1,4) = 1.0 ! Surface temp incr for cat 1
       !self%tsfcn(ixdir(idir),iydir(idir),1) = 1.0 ! Surface temp incr for cat 1
       !self%tocn(ixdir(idir),iydir(idir)) = 1.0 ! Surface temp incr for cat 1
       !self%cicen(ixdir(idir),iydir(idir),3) = 1.0 ! Surface temp incr for cat 1
       self%ssh(ixdir(idir),iydir(idir)) = 1.0 ! Surface temp incr for cat 1
    end do

  end subroutine dirac

  ! ------------------------------------------------------------------------------

  subroutine random(self)

    implicit none
    type(soca_field), intent(inout) :: self

    call check(self)

    call random_number(self%cicen); self%cicen=self%cicen-0.5_kind_real
    call random_number(self%hicen); self%hicen=self%hicen-0.5_kind_real
    call random_number(self%hsnon); self%hsnon=self%hsnon-0.5_kind_real        
    call random_number(self%tsfcn); self%tsfcn=self%tsfcn-0.5_kind_real
    
    call random_number(self%tocn); self%tocn=self%tocn-0.5_kind_real
    call random_number(self%socn); self%tocn=self%socn-0.5_kind_real
    call random_number(self%ssh); self%ssh=(self%ssh-0.5_kind_real)*self%geom%ocean%mask2d
    
  end subroutine random

  ! ------------------------------------------------------------------------------
  
  subroutine copy(self,rhs)
    
    implicit none
    type(soca_field), intent(inout) :: self
    type(soca_field), intent(in)    :: rhs
    integer :: nf

    call check_resolution(self, rhs)

    !nf = common_vars(self, rhs)

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
    self%hocn  = rhs%hocn
    
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
    integer :: is, ie, js, je, ncat, nzo
    call check_resolution(fld1, fld2)
    if (fld1%nf /= fld2%nf .or. fld1%geom%ocean%nzo /= fld2%geom%ocean%nzo) then
       print *,'STUFF:',fld1%nf,fld2%nf,fld1%geom%ocean%nzo, fld2%geom%ocean%nzo
       call abor1_ftn("soca_fields:field_prod error number of fields")
    endif

    ! Get compute domain
    is = fld1%geom%ocean%G%isc    
    ie = fld1%geom%ocean%G%iec
    js = fld1%geom%ocean%G%jsc    
    je = fld1%geom%ocean%G%jec    
    ncat = fld1%geom%ocean%ncat
    nzo = fld1%geom%ocean%nzo    
    
    zprod = 0.0_kind_real
    !----- OCEAN
    do ii = is, ie
       do jj = js, je
          zprod = zprod + fld1%ssh(ii,jj)*fld2%ssh(ii,jj)*fld1%geom%ocean%mask2d(ii,jj)       !SSH      
          do kk = 1, nzo !!! ADD 3D MASK!!!!!!
             zprod = zprod + fld1%tocn(ii,jj,kk)*fld2%tocn(ii,jj,kk)*fld1%geom%ocean%mask2d(ii,jj) &   !TOCN
                           + fld1%socn(ii,jj,kk)*fld2%socn(ii,jj,kk)*fld1%geom%ocean%mask2d(ii,jj)     !SOCN
          end do
       end do
    end do
    
    !----- SEA-ICE
    do ii = is, ie
       do jj = js, je
          do kk = 1, ncat
             zprod = zprod + fld1%cicen(ii,jj,kk)*fld2%cicen(ii,jj,kk)*fld1%geom%ocean%mask2d(ii,jj) & !CICEN
                           + fld1%hicen(ii,jj,kk)*fld2%hicen(ii,jj,kk)*fld1%geom%ocean%mask2d(ii,jj)   !HICEN          
          end do
       end do
    end do    
    
    allocate(zprod_allpes(mpp_npes()))

    call mpp_gather((/zprod/),zprod_allpes)
    call mpp_broadcast(zprod_allpes, mpp_npes(), mpp_root_pe())    
    zprod=sum(zprod_allpes)
    
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

    use ioda_locs_mod
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
    type(restart_file_type) :: ocean_restart    
    integer :: idr, idr_ocean

    
    type(ioda_locs)    :: locs
    type(ufo_geovals)    :: geovals
    !type(ufo_vars)    :: vars
    integer            :: nobs, nval
    
    iread = 0
    if (config_element_exists(c_conf,"read_from_file")) then
       iread = config_get_int(c_conf,"read_from_file")
    endif
    if (iread==0) then
       call log%warning("soca_fields:read_file: Inventing State")
       !call abor1_ftn("soca_fields: inventing state")                    
       !call invent_state(fld)
       call zeros(fld)
       call random_number(fld%cicen)

       sdate = config_get_string(c_conf,len(sdate),"date")
       WRITE(buf,*) 'validity date is: '//sdate
       call log%info(buf)
       call datetime_set(sdate, vdate)
    else
       basename = config_get_string(c_conf,len(basename),"basename")        
       !ocn_sfc_filename = config_get_string(c_conf,len(ocn_filename),"ocn_sfc_filename")
       ocn_filename = config_get_string(c_conf,len(ocn_filename),"ocn_filename")       

       !ocn_sfc_filename = trim(basename)//trim(ocn_sfc_filename)
       ocn_filename = trim(basename)//trim(ocn_filename)       
       ice_filename = config_get_string(c_conf,len(ice_filename),"ice_filename")       
       ice_filename = trim(basename)//trim(ice_filename)
       
       call fms_io_init()
       do ii = 1, fld%nf
          select case(fld%fldnames(ii))
          case ('ssh')
             idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'ave_ssh', fld%ssh(:,:), &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
             !call read_data(ocn_filename,"ave_ssh",fld%ssh(:,:),domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('tocn')
             call read_data(ocn_filename,"Temp",fld%tocn(:,:,:),domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('socn')             
             call read_data(ocn_filename,"Salt",fld%socn(:,:,:),domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('hocn')             
             call read_data(ocn_filename,"h",fld%hocn(:,:,:),domain=fld%geom%ocean%G%Domain%mpp_domain)             
          case ('cicen')
             !call read_data(ice_filename, 'part_size', fld%AOGCM%Ice%part_size, domain=fld%geom%ocean%G%Domain%mpp_domain)
             idr = register_restart_field(sis_restart, ice_filename, 'part_size', fld%AOGCM%Ice%part_size, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('hicen')
             !call read_data(ice_filename, 'h_ice', fld%AOGCM%Ice%h_ice, domain=fld%geom%ocean%G%Domain%mpp_domain)
             idr = register_restart_field(sis_restart, ice_filename, 'h_ice', fld%AOGCM%Ice%h_ice, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('hsnon')
             !call read_data(ice_filename, 'h_snow', fld%AOGCM%Ice%h_snow, domain=fld%geom%ocean%G%Domain%mpp_domain)             
             idr = register_restart_field(sis_restart, ice_filename, 'h_snow', fld%AOGCM%Ice%h_snow, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)             
          case ('qicnk')
             !call read_data(ice_filename, 'enth_ice', fld%AOGCM%Ice%enth_ice, domain=fld%geom%ocean%G%Domain%mpp_domain)
             idr = register_restart_field(sis_restart, ice_filename, 'enth_ice', fld%AOGCM%Ice%enth_ice, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('qsnon')
             !call read_data(ice_filename, 'enth_snow', fld%AOGCM%Ice%enth_snow, domain=fld%geom%ocean%G%Domain%mpp_domain)
             idr = register_restart_field(sis_restart, ice_filename, 'enth_snow', fld%AOGCM%Ice%enth_snow, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('tsfcn')             
             !call read_data(ice_filename, 'T_skin', fld%AOGCM%Ice%T_skin, domain=fld%geom%ocean%G%Domain%mpp_domain)
             idr = register_restart_field(sis_restart, ice_filename, 'T_skin', fld%AOGCM%Ice%T_skin, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)
          case ('sicnk')
             !call read_data(ice_filename, 'sal_ice', fld%AOGCM%Ice%sal_ice, domain=fld%geom%ocean%G%Domain%mpp_domain)             
             idr = register_restart_field(sis_restart, ice_filename, 'sal_ice', fld%AOGCM%Ice%sal_ice, &
                  domain=fld%geom%ocean%G%Domain%mpp_domain)             
          case default
             !print *,'Not reading var ',fld%fldnames(ii),' in file ',ocn_filename
             !call log%warning("soca_fields:read_file: Not reading var "//fld%fldnames(ii))
          end select
       end do
       call restore_state(sis_restart, directory='')
       call restore_state(ocean_restart, directory='')       

       sdate = config_get_string(c_conf,len(sdate),"date")
       WRITE(buf,*) 'validity date is: '//sdate
       call log%info(buf)
       call datetime_set(sdate, vdate)       

    endif

    call check(fld)
    
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

    !call geom_infotofile(fld%geom)
    
    filename = genfilename(c_conf,max_string_length,vdate)    
    WRITE(buf,*) 'field:write_file: writing '//filename
    call fckit_log%info(buf)

    call fms_io_init()
    call set_domain( fld%geom%ocean%G%Domain%mpp_domain )    
    do ii = 1, fld%nf
       select case(fld%fldnames(ii))
       case ('ssh')
          call write_data( filename, "ssh", fld%ssh, fld%geom%ocean%G%Domain%mpp_domain)
          call write_data( filename, "rossby_radius", fld%geom%ocean%rossby_radius, fld%geom%ocean%G%Domain%mpp_domain)          
       case ('tocn')
          call write_data( filename, "temp", fld%tocn, fld%geom%ocean%G%Domain%mpp_domain)          
       case ('socn')
          call write_data( filename, "salt", fld%socn, fld%geom%ocean%G%Domain%mpp_domain)
       case ('hocn')
          call write_data( filename, "h", fld%hocn, fld%geom%ocean%G%Domain%mpp_domain)
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
          !print *,'Not writing var ',fld%fldnames(ii),' in file ',filename
          !call log%warning("soca_fields:read_file: Not reading var "//fld%fldnames(ii))
          !call abor1_ftn("soca_fields: undefined variables")             
       end select

    end do
    call fms_io_exit()       
    
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
    integer :: jj, Nc2d

    call check(fld)

    Nc2d = sum(fld%geom%ocean%mask2d)
    
    pstat=0.0
    pstat(1,1) = minval(fld%cicen)
    pstat(2,1) = maxval(fld%cicen)
    pstat(3,1) = sqrt(sum(fld%cicen*fld%cicen)/real(Nc2d*fld%geom%ocean%ncat))

    pstat(1,2) = minval(fld%hicen)
    pstat(2,2) = maxval(fld%hicen)
    pstat(3,2) = sqrt(sum(fld%hicen*fld%hicen)/real(Nc2d*fld%geom%ocean%ncat))
    
    pstat(1,3) = minval(fld%hsnon)
    pstat(2,3) = maxval(fld%hsnon)
    pstat(2,3) = sqrt(sum(fld%hsnon*fld%hsnon)/real(Nc2d*fld%geom%ocean%ncat))

    pstat(1,4) = minval(fld%tsfcn)
    pstat(2,4) = maxval(fld%tsfcn)        

    pstat(1,5) = minval(fld%qsnon)
    pstat(2,5) = maxval(fld%qsnon)

    pstat(1,6) = minval(fld%sicnk)
    pstat(2,6) = maxval(fld%sicnk)        

    pstat(1,7) = minval(fld%qicnk)
    pstat(2,7) = maxval(fld%qicnk)

    pstat(1,8) = minval(fld%tocn)
    pstat(2,8) = maxval(fld%tocn)
    
    pstat(1,9) = minval(fld%socn)
    pstat(2,9) = maxval(fld%socn)    

    pstat(1,10) = minval(fld%ssh)
    pstat(2,10) = maxval(fld%ssh)    
    
    call dot_prod(fld, fld, zz)    
    !pstat(3,:) = sqrt(zz)

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

  end subroutine fldrms

  ! ------------------------------------------------------------------------------

  subroutine interp_tl(fld, locs, vars, geovals)
    use ioda_locs_mod
    use ufo_geovals_mod
    use ufo_vars_mod
    implicit none
    type(soca_field), intent(inout)   :: fld
    type(ioda_locs), intent(in)    :: locs
    type(ufo_vars),     intent(in)    :: vars    
    type(ufo_geovals), intent(inout) :: geovals

    call check(fld)

    call nicas_interph(fld, locs, vars, geovals)

  end subroutine interp_tl

  ! ------------------------------------------------------------------------------

  subroutine interp_ad(fld, locs, vars, geovals)
    use ioda_locs_mod
    use ufo_geovals_mod
    use ufo_vars_mod    
    implicit none
    type(soca_field), intent(inout) :: fld
    type(ioda_locs), intent(in)    :: locs
    type(ufo_vars),     intent(in)    :: vars        
    type(ufo_geovals), intent(inout) :: geovals    
    !character(2)                        :: op_type='AD'

    call check(fld)
    call nicas_interphad(fld, locs, vars, geovals)

  end subroutine interp_ad

  ! ------------------------------------------------------------------------------
  function get_obsop_index(horiz_interp, locs, interph_initialized)
    !> Returns index of interpolation object/obs operator
    use ioda_locs_mod  
    use soca_interph_mod
    implicit none

    type(soca_hinterp), dimension(10), intent(in) :: horiz_interp        !< HARD CODED ... HACK ...
    type(ioda_locs), intent(in)                    :: locs                !< HARD CODED ... HACK ...
    logical, dimension(10), intent(in)            :: interph_initialized !< HARD CODED ... HACK ...

    logical, dimension(10)            :: obs_type_test                   !< HARD CODED ... HACK ...
    integer :: nobs    
    integer :: cnt, cnt_obstype
    integer :: get_obsop_index

    ! Check nobs from interp object against nobs from locs 
    obs_type_test=(horiz_interp(:)%nobs.eq.locs%nlocs)

    ! Check for nobs matches between locs and horiz_interp
    cnt_obstype=count(obs_type_test)
    if (cnt_obstype.gt.1) then
       call abor1_ftn('Identification of obsop from nobs failed, needs further implementation')
    end if
    if (cnt_obstype.eq.0) then
       !New obs operator: add 1 to the last index for interp operator
       get_obsop_index=count(interph_initialized(:))+1
       return
    end if    

    ! Get index of interpolation object
    do cnt=1,size(obs_type_test)
       if (obs_type_test(cnt)) exit
    end do
    get_obsop_index=cnt

  end function get_obsop_index
  
  ! ------------------------------------------------------------------------------    

  subroutine initialize_interph(fld, locs, horiz_interp_p, interp_type)    
    use ioda_locs_mod  
    use soca_interph_mod

    implicit none

    type(soca_field), intent(in)             :: fld
    type(ioda_locs), intent(in)               :: locs
    type(soca_hinterp), pointer, intent(out) :: horiz_interp_p
    character(len=3), optional               :: interp_type     !< Forward: 'fwd' or adjoint: 'adj'
    
    logical, dimension(10), save                    :: interph_initialized = .false. !< HARD CODED ... HACK ... 
    type(soca_hinterp), dimension(10), save, target :: horiz_interp                  !< HARD CODED ... HACK ...
    logical, dimension(10)                          :: obs_type_test                 !< HARD CODED ... HACK ...    
    integer :: nobs
    integer :: cnt, cnt_obstype, obs_type_counter
    integer :: isc, iec, jsc, jec
    
    obs_type_counter = get_obsop_index(horiz_interp, locs, interph_initialized)

    ! Comming from adjoint, no need to initialize
    if (present(interp_type)) then   ! if interp_type is present, interp_type=adjoint
       horiz_interp_p => horiz_interp(obs_type_counter)
       return
    end if

    ! Compute interpolation weights if needed
    if (.not.(interph_initialized(obs_type_counter))) then       
       nobs = locs%nlocs
       ! Indices for compute domain (no halo)
       isc = fld%geom%ocean%G%isc
       iec = fld%geom%ocean%G%iec
       jsc = fld%geom%ocean%G%jsc
       jec = fld%geom%ocean%G%jec
       
       call horiz_interp(obs_type_counter)%interp_init(nobs)
       call horiz_interp(obs_type_counter)%interp_compute_weight(&
         &                       fld%geom%ocean%lon(isc:iec,jsc:jec),&
         &                       fld%geom%ocean%lat(isc:iec,jsc:jec),&
         &                       locs%lon,&
         &                       locs%lat)
       interph_initialized(obs_type_counter) = .true.
    end if
    horiz_interp_p => horiz_interp(obs_type_counter)

  end subroutine initialize_interph

  ! ------------------------------------------------------------------------------

  subroutine nicas_interph(fld, locs, ufovars, geovals)

    use ioda_locs_mod  
    use ufo_geovals_mod_c
    use ufo_geovals_mod
    use ufo_vars_mod
    use soca_constants, only : rho_i
    use soca_interph_mod
    use fckit_log_module, only : fckit_log
    
    implicit none

    type(soca_field), intent(inout)   :: fld
    type(ioda_locs), intent(in)       :: locs
    type(ufo_vars),     intent(in)    :: ufovars  
    type(ufo_geovals), intent(inout)  :: geovals

    integer :: ivar
    character(len=1024)  :: buf
    logical,allocatable :: mask(:), masko(:)               ! < mask (ncells, nlevels)
    real(kind=kind_real), allocatable :: lon(:), lat(:), lono(:), lato(:), fld_src(:), fld_dst(:)
    integer :: nobs
    integer :: nval

    ! interp stuff
    type(soca_hinterp), pointer :: horiz_interp_p

    integer, allocatable :: imask(:)
    real(kind=kind_real), allocatable :: area(:),vunit(:)
    real(kind=kind_real), allocatable :: obs_field(:,:), mod_field(:,:), obsout(:)

    integer :: icat, ilev
    character(len=160) :: record
    integer :: isc, iec, jsc, jec

!!$    !HACK TO TEST ADJOINT OF INTERPOLTION
!!$    type(ufo_locs)        :: locs_test
!!$    integer, parameter :: nobs_test=100000
!!$    real(kind_real) :: lats(nobs_test), lons(nobs_test)
!!$    real(kind=kind_real),allocatable :: dx1(:,:),dx2(:,:)
!!$    real(kind=kind_real) :: dy1(nobs_test),dy2(nobs_test)
!!$    
!!$    ! Test horizontal inerpolation
!!$    call random_number(lons)
!!$    lons=-30.0*lons
!!$    call random_number(lats)    
!!$    lats=75.0*lats
!!$    allocate(dx1(size(fld%cicen(:,:,1),1),size(fld%cicen(:,:,1),1)))
!!$    allocate(dx2(size(fld%cicen(:,:,1),1),size(fld%cicen(:,:,1),1)))
!!$    call ufo_locs_create(locs_test, nobs_test, lats, lons)
!!$    call random_number(dx1)
!!$    call random_number(dy1)
!!$    dx2=0.0
!!$    dy2=0.0
!!$    call initialize_interph(fld, locs_test, horiz_interp_p)
!!$    call horiz_interp_p%interp_apply(dx1, dy2)           ! dy2=I.dx1
!!$    call horiz_interp_p%interpad_apply(dx2, dy1)         ! dx2=I^T.dy1
!!$    print *,'dx1.dx2-dy1.dy2=',sum(dx2*dx1)-sum(dy2*dy1) ! dx1.dx2=dy1.dy2
!!$
!!$    call abor1_ftn("done testing interp")
!!$    !End test

    !call t_from_pt(pt_in,sp_in,p_in,lon_in,lat_in)
    
    nobs = locs%nlocs
    do ivar = 1, ufovars%nv
       write(record,*) "nicas_interph: ",ufovars%fldnames(ivar)
       call fckit_log%info(record)

       select case (trim(ufovars%fldnames(ivar)))
       case ("ice_concentration","ice_thickness")
          nval = fld%geom%ocean%ncat

       case ("steric_height","sea_surface_height_above_geoid","ocean_upper_level_temperature")
          nval = 1

       case ("ocean_potential_temperature","ocean_salinity")
          nval = fld%geom%ocean%nzo

       case default
          write(record,*) "nicas_interph: Doing nothing "
          call fckit_log%info(record)          

       end select
       
       ! Allocate GeoVaLs (fields at locations)
       if (nval.eq.0) call abor1_ftn("Wrong nval: nval = 0")
       geovals%geovals(ivar)%nval = nval
       if (allocated(geovals%geovals(ivar)%vals)) then
          deallocate(geovals%geovals(ivar)%vals)
       end if
       allocate(geovals%geovals(ivar)%vals(nval,nobs))
       geovals%geovals(ivar)%vals(:,:)=0.0_kind_real    
       geovals%lalloc = .true.       
       geovals%linit = .true.
       
       ! Initialize horizontal inerpolation
       call initialize_interph(fld, locs, horiz_interp_p)
       write(record,*) "nicas_interph: ",ufovars%fldnames(ivar)
       call fckit_log%info(record)
       ! Indices for compute domain (no halo)
       isc = fld%geom%ocean%G%isc
       iec = fld%geom%ocean%G%iec
       jsc = fld%geom%ocean%G%jsc
       jec = fld%geom%ocean%G%jec
              
       select case (trim(ufovars%fldnames(ivar)))
          
       case ("ice_concentration")
          do icat = 1,fld%geom%ocean%ncat
             call horiz_interp_p%interp_apply(fld%cicen(isc:iec,jsc:jec,icat+1)*fld%geom%ocean%mask2d(isc:iec,jsc:jec),&
                  &geovals%geovals(ivar)%vals(icat,:))
          end do
       case ("ice_thickness")
          do icat = 1,fld%geom%ocean%ncat
             call horiz_interp_p%interp_apply(fld%hicen(isc:iec,jsc:jec,icat)*fld%geom%ocean%mask2d(isc:iec,jsc:jec),&
                  &geovals%geovals(ivar)%vals(icat,:))
          end do
          
       case ("sea_surface_height_above_geoid","steric_height")
          call horiz_interp_p%interp_apply(fld%ssh(isc:iec,jsc:jec)*fld%geom%ocean%mask2d(isc:iec,jsc:jec),&
               &geovals%geovals(ivar)%vals(1,:))

       case ("ocean_potential_temperature")
          do ilev = 1, fld%geom%ocean%nzo
             call horiz_interp_p%interp_apply(fld%tocn(isc:iec,jsc:jec,ilev), geovals%geovals(ivar)%vals(ilev,:))
          end do
          
       case ("ocean_salinity")
          do ilev = 1, fld%geom%ocean%nzo
             call horiz_interp_p%interp_apply(fld%socn(isc:iec,jsc:jec,ilev), geovals%geovals(ivar)%vals(ilev,:))
          end do

       case ("ocean_layer_thickness")
          do ilev = 1, fld%geom%ocean%nzo
             call horiz_interp_p%interp_apply(fld%hocn(isc:iec,jsc:jec,ilev), geovals%geovals(ivar)%vals(ilev,:))
          end do

       case ("ocean_upper_level_temperature")          
          call horiz_interp_p%interp_apply(fld%tocn(isc:iec,jsc:jec,1), geovals%geovals(ivar)%vals(1,:))
          
       end select
    end do

  end subroutine nicas_interph

  ! ------------------------------------------------------------------------------
  
  subroutine nicas_interphad(fld, locs, ufovars, geovals)

    use ioda_locs_mod  
    use ufo_geovals_mod_c
    use ufo_geovals_mod
    use ufo_vars_mod
    use soca_constants, only : rho_i
    use soca_interph_mod
    use fckit_log_module, only : fckit_log

    implicit none
    
    type(soca_field), intent(inout)  :: fld
    type(ioda_locs), intent(in)       :: locs
    type(ufo_geovals), intent(in)    :: geovals    
    type(ufo_vars)    :: ufovars    

    integer :: ivar, geovals_dim1, cnt_fld
    character(len=1024)  :: buf
    logical,allocatable :: mask(:), masko(:)               ! < mask (ncells, nlevels)
    integer :: nobs, nval

    type(soca_hinterp), pointer :: horiz_interp_p
    integer :: icat, ilev
    character(len=160) :: record
    integer :: isc, iec, jsc, jec


    ! Indices for compute domain (no halo)
    isc = fld%geom%ocean%G%isc
    iec = fld%geom%ocean%G%iec
    jsc = fld%geom%ocean%G%jsc
    jec = fld%geom%ocean%G%jec

    call initialize_interph(fld, locs, horiz_interp_p, interp_type='adj')

    do ivar = 1, ufovars%nv
       write(record,*) "nicas_interphad: ",trim(ufovars%fldnames(ivar))
       call fckit_log%info(record)                           

       select case (trim(ufovars%fldnames(ivar)))
       case ("ice_concentration")
          do icat = 1,fld%geom%ocean%ncat
             call horiz_interp_p%interpad_apply(fld%cicen(isc:iec,jsc:jec,icat+1), geovals%geovals(ivar)%vals(icat,:))
          enddo

       case ("ice_thickness")
          do icat = 1,fld%geom%ocean%ncat
             call horiz_interp_p%interpad_apply(fld%hicen(isc:iec,jsc:jec,icat), geovals%geovals(ivar)%vals(icat,:))
          enddo

       case ("sea_surface_height_above_geoid","steric_height") !!!! steric height sould be  different case
          call horiz_interp_p%interpad_apply(fld%ssh(isc:iec,jsc:jec), geovals%geovals(ivar)%vals(1,:))
          
       case ("ocean_potential_temperature")
          do ilev = 1, fld%geom%ocean%nzo
             call horiz_interp_p%interpad_apply(fld%tocn(isc:iec,jsc:jec,ilev), geovals%geovals(ivar)%vals(ilev,:))
          end do
          
       case ("ocean_salinity")
          do ilev = 1, fld%geom%ocean%nzo
             call horiz_interp_p%interpad_apply(fld%socn(isc:iec,jsc:jec,ilev), geovals%geovals(ivar)%vals(ilev,:))
          end do

       case ("ocean_upper_level_temperature")
          call horiz_interp_p%interpad_apply(fld%tocn(isc:iec,jsc:jec,1), geovals%geovals(ivar)%vals(1,:))
          
       end select
    end do
  end subroutine nicas_interphad

  ! ------------------------------------------------------------------------------

  subroutine convert_to_ug(self, ug)
    use unstructured_grid_mod
    use soca_thermo
    use tools_const, only: deg2rad
    
    implicit none
    type(soca_field), intent(in) :: self
    type(unstructured_grid), intent(inout) :: ug
    integer :: jx,jy,jz,jk
    integer :: nz_total     ! Total number of levels in the 3D fields
    integer :: n_vars       ! Number of 3D variables 
    integer :: n_surf_vars  ! Number of surf vars (sould be 0 for ocean/ice)

    !nicas stuff
    integer :: nc0a, nl0, nv, nts
    real(kind=kind_real), allocatable :: lon(:), lat(:), area(:), vunit(:), rndnum(:)
    integer, allocatable :: imask(:,:)

    !Grid stuff
    integer :: isc, iec, jsc, jec, jjj

    !Get indices for compute domain (no halo)
    isc = self%geom%ocean%G%isc
    iec = self%geom%ocean%G%iec    
    jsc = self%geom%ocean%G%jsc
    jec = self%geom%ocean%G%jec
    
    nv = self%geom%ocean%ncat + 1
    nl0 = 1
    nts = 1
    nc0a = (iec - isc + 1) * (jec - jsc + 1 )

    allocate( lon(nc0a), lat(nc0a), area(nc0a) )
    allocate( vunit(nl0) )
    allocate( imask(nc0a, nl0) )    

    lon = deg2rad*reshape( self%geom%ocean%lon(isc:iec, jsc:jec), (/nc0a/) )
    lat = deg2rad*reshape( self%geom%ocean%lat(isc:iec, jsc:jec), (/nc0a/) ) 
    
    area = reshape( self%geom%ocean%cell_area(isc:iec, jsc:jec), (/nc0a/) )

    do jz = 1, nl0       
       vunit(jz) = real(jz)
       imask(1:nc0a,jz) = reshape( self%geom%ocean%mask2d(isc:iec, jsc:jec), (/nc0a/) )
    end do

    call create_unstructured_grid(ug, nc0a, nl0, nv, 1, lon, lat, area, vunit, imask)

    do jk = 1, nv - 1
       ug%fld(:, 1, jk, 1) = reshape( self%cicen(isc:iec, jsc:jec, jk+1), (/nc0a/) )
    end do
    jk = nv
    ug%fld(:, 1, jk, 1) = reshape( self%ssh(isc:iec, jsc:jec), (/nc0a/) )

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
    integer :: nv

    !Grid stuff
    integer :: isc, iec, jsc, jec
    
    !Get indices for compute domain (no halo)
    isc = self%geom%ocean%G%isc
    iec = self%geom%ocean%G%iec    
    jsc = self%geom%ocean%G%jsc
    jec = self%geom%ocean%G%jec
    
    nv = self%geom%ocean%ncat + 1    
    do jk = 1, nv - 1
       self%cicen(isc:iec, jsc:jec,jk+1) = reshape( ug%fld(:, 1, jk, 1), (/(iec - isc + 1 ), (jec - jsc + 1 )/) )
    end do
    jk = nv
    self%ssh(isc:iec, jsc:jec) = reshape( ug%fld(:, 1, jk, 1), (/(iec - isc + 1 ), (jec - jsc + 1 )/) )

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
