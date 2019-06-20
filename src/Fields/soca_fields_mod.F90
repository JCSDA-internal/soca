! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Handle fields for the  model

module soca_fields

  use config_mod
  use datetime_mod
  use duration_mod
  use fckit_log_module, only : log, fckit_log
  use fckit_mpi_module
  use fms_mod,    only: read_data, write_data, set_domain
  use fms_io_mod, only : fms_io_init, fms_io_exit,&
       &register_restart_field, restart_file_type,&
       &restore_state, query_initialized,&
       &free_restart_type, save_restart
  use iso_c_binding
  use kinds
  use MOM_remapping,       only : remapping_CS, initialize_remapping, remapping_core_h, end_remapping
  use mpp_domains_mod, only : mpp_update_domains
  use random_mod
  use soca_fieldsutils_mod
  use soca_geom_mod_c
  use soca_geom_mod, only : geom_get_domain_indices
  use soca_mom6
  use soca_utils
  use soca_bumpinterp2d_mod
  use soca_getvaltraj_mod
  use soca_ocnsfc_mod
  use ufo_vars_mod
  use unstructured_grid_mod
  use variables_mod

  implicit none

  private

  public :: soca_field, &
       & create, delete, zeros, dirac, random, copy, create_copy,&
       & self_add, self_schur, self_sub, self_mul, axpy, &
       & dot_prod, add_incr, diff_incr, &
       & read_file, write_file, gpnorm, fldrms, soca_fld2file, &
       & change_resol, check, &
       & field_to_ug, field_from_ug, ug_coord

  interface create
     module procedure create_constructor, create_copy
  end interface create

  ! ------------------------------------------------------------------------------
  !> Fortran derived type to hold fields
  type :: soca_field
     type(soca_geom), pointer          :: geom           !< MOM6 & SIS2 Geometry
     integer                           :: nf             !< Number of fields
     character(len=128)                :: gridfname      !< Grid file name
     character(len=128)                :: cicefname      !< Fields file name for cice
     character(len=128)                :: momfname       !< Fields file name for mom

     ! Sea-ice state variables
     real(kind=kind_real), allocatable :: cicen(:,:,:)   !< Sea-ice fraction                 (nx,ny,ncat+1)
     real(kind=kind_real), allocatable :: hicen(:,:,:)   !< Sea-ice mass/m2                  (nx,ny,ncat) [kg/m2]

     ! Ocean state variables
     real(kind=kind_real), allocatable :: socn(:,:,:)    !< Ocean Practical Salinity         (nx,ny,nzo)
     real(kind=kind_real), allocatable :: tocn(:,:,:)    !< Ocean Potential Temperature, ref to p=0      (nx,ny,nzo)
     real(kind=kind_real), allocatable :: ssh(:,:)       !< Sea-surface height (nx,ny)
     real(kind=kind_real), allocatable :: hocn(:,:,:)    !< DA layer thickness (nx,ny,nzo)

     ! Ocean diagnostics
     real(kind=kind_real), allocatable :: mld(:,:)           !< Mixed layer depth (nx,ny)
     real(kind=kind_real), allocatable :: layer_depth(:,:,:) !< Mid-layer depth (nx,ny,nz0)

     ! Ocean surface fields
     type(soca_ocnsfc_type) :: ocnsfc !< Surface fields needed for cool skin ufo

     character(len=5),     allocatable :: fldnames(:)    !< Variable identifiers             (nf)

  end type soca_field

contains

  ! ------------------------------------------------------------------------------
  !> Create a field from geometry and variables
  subroutine create_constructor(self, geom, vars)
    type(soca_field),          intent(inout) :: self
    type(soca_geom),  pointer, intent(inout) :: geom
    type(oops_vars),              intent(in) :: vars
    integer :: ivar

    ! Allocate
    call soca_field_alloc(self, geom)
    !call zeros(self)

    ! Associate geometry
    self%geom => geom

    ! Set fields numbers and names
    self%nf   = vars%nv
    allocate(self%fldnames(self%nf))
    self%fldnames(:)=vars%fldnames(:)

    call check(self)

  end subroutine create_constructor

  ! ------------------------------------------------------------------------------

  subroutine create_copy(self, rhs_fld)
    ! Construct a field from an other field, lhs_fld=rhs_fld
    type(soca_field), intent(inout) :: self
    type(soca_field), intent(in)    :: rhs_fld
    integer :: ivar!, unit, nxny(2)

    ! Allocate and copy fields
    call soca_field_alloc(self, rhs_fld%geom)
    call zeros(self)
    !call copy(self,rhs_fld)

    ! Associate geometry
    self%geom => rhs_fld%geom

    ! Set fields numbers and names
    self%nf   = rhs_fld%nf
    allocate(self%fldnames(self%nf))
    self%fldnames(:)=rhs_fld%fldnames(:)

    call check(self)

  end subroutine create_copy

  ! ------------------------------------------------------------------------------

  subroutine soca_field_alloc(self, geom)
    type (soca_field), intent(inout) :: self
    type(soca_geom),      intent(in) :: geom

    integer :: isd, ied, jsd, jed, nzo, nzi, nzs
    integer :: ncat, km
    character(7) :: domain_type

    ! Short cut to ice geometry
    ncat = geom%ice_column%ncat
    nzi = geom%ice_column%nzi
    nzs = geom%ice_column%nzs
    nzo = geom%nzo

    ! Indices for data domain (with halo)
    call geom_get_domain_indices(geom, "data   ", isd, ied, jsd, jed)

    ! Allocate ocean state
    allocate(self%tocn(isd:ied,jsd:jed,nzo))
    allocate(self%socn(isd:ied,jsd:jed,nzo))
    allocate(self%ssh(isd:ied,jsd:jed))
    allocate(self%hocn(isd:ied,jsd:jed,nzo))
    allocate(self%mld(isd:ied,jsd:jed))
    allocate(self%layer_depth(isd:ied,jsd:jed,nzo))

    ! Allocate sea-ice state
    km = ncat + 1
    allocate(self%cicen(isd:ied, jsd:jed, km))
    allocate(self%hicen(isd:ied, jsd:jed, ncat))

    ! Allocate surface fields for cool skin
    call self%ocnsfc%create(geom)

  end subroutine soca_field_alloc

  ! ------------------------------------------------------------------------------

  subroutine delete(self)
    type (soca_field), intent(inout) :: self

    ! Deallocate ocean state
    deallocate(self%tocn)
    deallocate(self%socn)
    deallocate(self%ssh)
    deallocate(self%hocn)
    deallocate(self%mld)
    deallocate(self%layer_depth)

    ! Deallocate sea-ice state
    deallocate(self%cicen)
    deallocate(self%hicen)

    ! Deallocate surface fields for cool skin
    call self%ocnsfc%delete()

    ! Deassociate geometry
    nullify(self%geom)

  end subroutine delete


  ! ------------------------------------------------------------------------------

  subroutine zeros(self)
    type(soca_field), intent(inout) :: self

    call check(self)

    self%cicen = 0.0_kind_real
    self%hicen = 0.0_kind_real

    self%socn = 0.0_kind_real
    self%tocn = 0.0_kind_real
    self%ssh = 0.0_kind_real
    self%hocn = 0.0_kind_real

    self%mld = 0.0_kind_real

    call self%ocnsfc%zeros()

  end subroutine zeros

  ! ------------------------------------------------------------------------------

  subroutine ones(self)
    type(soca_field), intent(inout) :: self

    call check(self)

    self%cicen = 1.0_kind_real
    self%hicen = 1.0_kind_real

    self%socn = 1.0_kind_real
    self%tocn = 1.0_kind_real
    self%ssh = 1.0_kind_real
    self%hocn = 1.0_kind_real

    self%mld = 1.0_kind_real
  end subroutine ones

  ! ------------------------------------------------------------------------------

  subroutine dirac(self, c_conf)
    type(soca_field), intent(inout) :: self
    type(c_ptr),         intent(in) :: c_conf   !< Configuration

    integer :: isc, iec, jsc, jec
    integer :: ndir,n,size,rank,info
    integer,allocatable :: ixdir(:),iydir(:),izdir(:),ifdir(:)
    character(len=3) :: idirchar
    type(fckit_mpi_comm) :: f_comm

    call check(self)

    ! Get MPI communicator
    f_comm = fckit_mpi_comm()

    ! Get Diracs size
    ndir = config_get_data_dimension(c_conf,'ixdir')
    if ((config_get_data_dimension(c_conf,'iydir')/=ndir) .or. &
         (config_get_data_dimension(c_conf,'izdir')/=ndir) .or. &
         (config_get_data_dimension(c_conf,'ifdir')/=ndir)) &
         & call abor1_ftn('qg_fields_dirac: inconsistent sizes for ixdir, iydir, izdir, ipdir and ifdir')

    ! Allocation
    allocate(ixdir(ndir))
    allocate(iydir(ndir))
    allocate(izdir(ndir))
    allocate(ifdir(ndir))

    ! Get Diracs positions
    call config_get_int_vector(c_conf,'ixdir',ixdir)
    call config_get_int_vector(c_conf,'iydir',iydir)
    call config_get_int_vector(c_conf,'izdir',izdir)
    call config_get_int_vector(c_conf,'ifdir',ifdir)

    ! get PE domain bounds
    call geom_get_domain_indices(self%geom, "compute", isc, iec, jsc, jec)

    ! Setup Diracs
    call zeros(self)
    do n=1,ndir
       ! skip this index if not in the bounds of this PE
       if (ixdir(n) > iec .or. ixdir(n) < isc) cycle
       if (iydir(n) > jec .or. iydir(n) < jsc) cycle

       if (ifdir(n)==1) self%tocn(ixdir(n),iydir(n),izdir(n)) = 1.0
       if (ifdir(n)==2) self%socn(ixdir(n),iydir(n),izdir(n)) = 1.0
       if (ifdir(n)==3) self%ssh(ixdir(n),iydir(n)) = 1.0
       if (ifdir(n)==4) self%cicen(ixdir(n),iydir(n),izdir(n)) = 1.0
       if (ifdir(n)==5) self%hicen(ixdir(n),iydir(n),izdir(n)) = 1.0
    end do

  end subroutine dirac

  ! ------------------------------------------------------------------------------

  subroutine random(self)
    type(soca_field), intent(inout) :: self
    integer, parameter :: rseed = 1 ! constant for reproducability of tests

    call check(self)

    call normal_distribution(self%cicen, 0.0_kind_real, 1.0_kind_real, rseed)
    call normal_distribution(self%hicen, 0.0_kind_real, 1.0_kind_real, rseed)
    call normal_distribution(self%tocn,  0.0_kind_real, 1.0_kind_real, rseed)
    call normal_distribution(self%socn,  0.0_kind_real, 1.0_kind_real, rseed)
    call normal_distribution(self%ssh,   0.0_kind_real, 1.0_kind_real, rseed)
    call self%ocnsfc%random()

  end subroutine random

  ! ------------------------------------------------------------------------------

  subroutine copy(self,rhs)
    type(soca_field), intent(inout) :: self
    type(soca_field),    intent(in) :: rhs

    integer :: nf

    call check_resolution(self, rhs)

    !nf = common_vars(self, rhs)

    ! Associate geometry
    if (.not.associated(self%geom)) self%geom => rhs%geom

    ! Set fields numbers and names
    self%nf   = rhs%nf
    if (.not.allocated(self%fldnames)) allocate(self%fldnames(self%nf))
    self%fldnames(:)=rhs%fldnames(:)

    self%cicen = rhs%cicen
    self%hicen = rhs%hicen

    self%socn  = rhs%socn
    self%tocn  = rhs%tocn
    self%ssh   = rhs%ssh
    self%hocn  = rhs%hocn
    self%mld   = rhs%mld
    self%layer_depth   = rhs%layer_depth

    ! Ocean surface
    call self%ocnsfc%copy(rhs%ocnsfc)

    return
  end subroutine copy

  ! ------------------------------------------------------------------------------

  subroutine self_add(self,rhs)
    type(soca_field), intent(inout) :: self
    type(soca_field),    intent(in) :: rhs

    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen = self%cicen + rhs%cicen
    self%hicen = self%hicen + rhs%hicen

    self%tocn = self%tocn + rhs%tocn
    self%socn = self%socn + rhs%socn
    self%ssh = self%ssh + rhs%ssh
    self%hocn = self%hocn + rhs%hocn

    ! Ocean surface
    call self%ocnsfc%add(rhs%ocnsfc)

  end subroutine self_add

  ! ------------------------------------------------------------------------------

  subroutine self_schur(self,rhs)
    type(soca_field), intent(inout) :: self
    type(soca_field),    intent(in) :: rhs

    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen=self%cicen*rhs%cicen
    self%hicen=self%hicen*rhs%hicen

    self%tocn=self%tocn*rhs%tocn
    self%socn=self%socn*rhs%socn
    self%ssh=self%ssh*rhs%ssh
    self%hocn=self%hocn*rhs%hocn

    ! Ocean surface
    call self%ocnsfc%schur(rhs%ocnsfc)

    return
  end subroutine self_schur

  ! ------------------------------------------------------------------------------

  subroutine self_sub(self,rhs)
    type(soca_field), intent(inout) :: self
    type(soca_field),    intent(in) :: rhs

    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen=self%cicen-rhs%cicen
    self%hicen=self%hicen-rhs%hicen

    self%socn=self%socn-rhs%socn
    self%tocn=self%tocn-rhs%tocn
    self%ssh=self%ssh-rhs%ssh
    self%hocn=self%hocn-rhs%hocn

    ! Ocean surface
    call self%ocnsfc%sub(rhs%ocnsfc)

  end subroutine self_sub

  ! ------------------------------------------------------------------------------

  subroutine self_mul(self,zz)
    type(soca_field),  intent(inout) :: self
    real(kind=kind_real), intent(in) :: zz

    call check(self)

    self%cicen = zz * self%cicen
    self%hicen = zz * self%hicen

    self%tocn = zz * self%tocn
    self%socn = zz * self%socn
    self%ssh = zz * self%ssh
    self%hocn = zz * self%hocn

    ! Ocean surface
    call self%ocnsfc%mul(zz)

  end subroutine self_mul

  ! ------------------------------------------------------------------------------

  subroutine axpy(self,zz,rhs)
    type(soca_field),  intent(inout) :: self
    real(kind=kind_real), intent(in) :: zz
    type(soca_field),     intent(in) :: rhs

    integer :: nf

    call check_resolution(self, rhs)

    nf = common_vars(self, rhs)

    self%cicen = self%cicen + zz * rhs%cicen
    self%hicen = self%hicen + zz * rhs%hicen

    self%tocn = self%tocn + zz * rhs%tocn
    self%socn = self%socn + zz * rhs%socn
    self%ssh = self%ssh + zz * rhs%ssh
    self%hocn = self%hocn + zz * rhs%hocn

    ! Ocean surface
    call self%ocnsfc%axpy(zz, rhs%ocnsfc)

  end subroutine axpy

  ! ------------------------------------------------------------------------------

  subroutine dot_prod(fld1,fld2,zprod)
    type(soca_field),      intent(in) :: fld1
    type(soca_field),      intent(in) :: fld2
    real(kind=kind_real), intent(out) :: zprod

    real(kind=kind_real) :: zprod_allpes
    integer :: ii, jj, kk
    integer :: is, ie, js, je, ncat, nzo, myrank
    type(fckit_mpi_comm) :: f_comm

    ! Setup Communicator
    f_comm = fckit_mpi_comm()

    call check_resolution(fld1, fld2)
    if (fld1%nf /= fld2%nf .or. fld1%geom%nzo /= fld2%geom%nzo) then
       call abor1_ftn("soca_fields:field_prod error number of fields")
    endif

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(fld1%geom, "compute", is, ie, js, je)

    ! Get ice categories and ocean levels
    ncat = fld1%geom%ncat
    nzo = fld1%geom%nzo

    zprod = 0.0_kind_real
    !----- OCEAN
    do ii = is, ie
       do jj = js, je
          if (fld1%geom%mask2d(ii,jj)==1) then
             zprod = zprod + fld1%ssh(ii,jj)*fld2%ssh(ii,jj)               ! SSH
             do kk = 1, nzo
                zprod = zprod + fld1%tocn(ii,jj,kk)*fld2%tocn(ii,jj,kk) &  ! TOCN
                     + fld1%socn(ii,jj,kk)*fld2%socn(ii,jj,kk)    ! SOCN
             end do
          end if
       end do
    end do
    call f_comm%barrier()
    myrank = f_comm%rank()

    !----- SEA-ICE
    do ii = is, ie
       do jj = js, je
          do kk = 1, ncat
             if (fld1%geom%mask2d(ii,jj)==1) then
                zprod = zprod + fld1%cicen(ii,jj,kk+1)*fld2%cicen(ii,jj,kk+1) & !CICEN
                     + fld1%hicen(ii,jj,kk)*fld2%hicen(ii,jj,kk)       !HICEN
             end if
          end do
       end do
    end do

!!$    !----- OCEAN Surface
    do ii = is, ie
       do jj = js, je
          if (fld1%geom%mask2d(ii,jj)==1) then
             zprod = zprod + fld1%ocnsfc%sw_rad(ii,jj)*fld2%ocnsfc%sw_rad(ii,jj) &
                           + fld1%ocnsfc%lw_rad(ii,jj)*fld2%ocnsfc%lw_rad(ii,jj) &
                           + fld1%ocnsfc%latent_heat(ii,jj)*fld2%ocnsfc%latent_heat(ii,jj) &
                           + fld1%ocnsfc%sens_heat(ii,jj)*fld2%ocnsfc%sens_heat(ii,jj) &
                           + fld1%ocnsfc%fric_vel(ii,jj)*fld2%ocnsfc%fric_vel(ii,jj)
          end if
       end do
    end do

    ! Get global dot product
    call f_comm%allreduce(zprod, zprod_allpes, fckit_mpi_sum())
    zprod = zprod_allpes

  end subroutine dot_prod

  ! ------------------------------------------------------------------------------

  subroutine add_incr(self,rhs)
    type(soca_field), intent(inout) :: self
    type(soca_field), intent(in)    :: rhs

    integer, save :: cnt_outer = 1
    real(kind=kind_real), allocatable :: incr(:,:,:)
    integer :: is, ie, js, je, i, j, nz
    character(len=800) :: filename, str_cnt

    call check(self)
    call check(rhs)

    ! Add increment to field
    call self_add(self,rhs)
    call self%ocnsfc%add(rhs%ocnsfc)

    ! Save increment for outer loop cnt_outer
    write(str_cnt,*) cnt_outer
    filename='incr.'//adjustl(trim(str_cnt))//'.nc'
    call soca_fld2file(rhs, filename)

    ! Update outer loop counter
    cnt_outer = cnt_outer + 1

    return
  end subroutine add_incr

  ! ------------------------------------------------------------------------------

  subroutine diff_incr(lhs,x1,x2)
    type(soca_field), intent(inout) :: lhs
    type(soca_field), intent(in)    :: x1
    type(soca_field), intent(in)    :: x2

    call check(lhs)
    call check(x1)
    call check(x2)

    call zeros(lhs)

    lhs%cicen = x1%cicen - x2%cicen
    lhs%hicen = x1%hicen - x2%hicen

    lhs%tocn = x1%tocn - x2%tocn
    lhs%socn = x1%socn - x2%socn
    lhs%ssh = x1%ssh - x2%ssh
    lhs%hocn = x1%hocn - x2%hocn

    call lhs%ocnsfc%diff_incr(x1%ocnsfc, x2%ocnsfc)

  end subroutine diff_incr

  ! ------------------------------------------------------------------------------

  subroutine change_resol(fld,rhs)
    type(soca_field), intent(inout) :: fld
    type(soca_field), intent(in)    :: rhs

    call check(fld)
    call check(rhs)
    call copy(fld,rhs)
    call fld%ocnsfc%copy(rhs%ocnsfc)

    return
  end subroutine change_resol

  ! ------------------------------------------------------------------------------

  subroutine read_file(fld, c_conf, vdate)
    type(soca_field), intent(inout) :: fld      !< Fields
    type(c_ptr),         intent(in) :: c_conf   !< Configuration
    type(datetime),   intent(inout) :: vdate    !< DateTime

    integer, parameter :: max_string_length=800
    character(len=max_string_length) :: ocn_filename
    character(len=max_string_length) :: ice_filename, basename, incr_filename
    character(len=20) :: sdate
    character(len=1024)  :: buf
    integer :: iread, ii
    logical :: vert_remap=.false.
    character(len=max_string_length) :: remap_filename
    real(kind=kind_real), allocatable :: h_common(:,:,:)    !< layer thickness to remap to
    type(restart_file_type) :: sis_restart
    type(restart_file_type) :: ocean_restart
    type(restart_file_type) :: ocean_remap_restart
    integer :: idr, idr_ocean
    integer            :: nobs, nval, pe, ierror
    integer :: is, ie, js, je, i, j, k, nl, nz
    type(remapping_CS)  :: remapCS

    ! Set default iread to 0
    iread = 0
    if (config_element_exists(c_conf,"read_from_file")) then
       iread = config_get_int(c_conf,"read_from_file")
    endif

    ! Check if vertical remapping needs to be applied
    if (config_element_exists(c_conf,"remap_filename")) then
       vert_remap = .true.
       remap_filename = config_get_string(c_conf,len(remap_filename),"remap_filename")
       remap_filename = trim(remap_filename)

       ! Get Indices for data domain and allocate common layer depth array
       call geom_get_domain_indices(fld%geom, "data   ", is, ie, js, je)
       nz=size(fld%hocn, dim=3)
       allocate(h_common(is:ie,js:je,nz))

       ! Read from file
       call fms_io_init()
       idr_ocean = register_restart_field(ocean_remap_restart, remap_filename, 'h', fld%hocn(:,:,:), &
            domain=fld%geom%G%Domain%mpp_domain)
       call restore_state(ocean_remap_restart, directory='')
       call fms_io_exit()
       h_common = fld%hocn

    end if

    ! iread = 0: Invent state
    if (iread==0) then
       call zeros(fld)
       sdate = config_get_string(c_conf,len(sdate),"date")
       call datetime_set(sdate, vdate)
    end if

    ! iread = 1 (state) or 3 (increment): Read restart file
    if ((iread==1).or.(iread==3)) then
       ! Read ocean surface fields
       call fld%ocnsfc%read_restart(c_conf, fld%geom, fld%fldnames)

       basename = config_get_string(c_conf,len(basename),"basename")
       ocn_filename = config_get_string(c_conf,len(ocn_filename),"ocn_filename")
       ocn_filename = trim(basename)//trim(ocn_filename)
       ice_filename = config_get_string(c_conf,len(ice_filename),"ice_filename")
       ice_filename = trim(basename)//trim(ice_filename)

       call fms_io_init()
       do ii = 1, fld%nf
          select case(fld%fldnames(ii))
             ! Ocean
          case ('ssh')
             idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'ave_ssh', fld%ssh(:,:), &
                  domain=fld%geom%G%Domain%mpp_domain)
          case ('tocn')
             idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'Temp', fld%tocn(:,:,:), &
                  domain=fld%geom%G%Domain%mpp_domain)
          case ('socn')
             idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'Salt', fld%socn(:,:,:), &
                  domain=fld%geom%G%Domain%mpp_domain)
          case ('hocn')
             idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'h', fld%hocn(:,:,:), &
                  domain=fld%geom%G%Domain%mpp_domain)
             ! Sea-ice
          case ('cicen')
             idr = register_restart_field(sis_restart, ice_filename, 'part_size', fld%cicen(:,:,:), &
                  domain=fld%geom%G%Domain%mpp_domain)
          case ('hicen')
             idr = register_restart_field(sis_restart, ice_filename, 'h_ice', fld%hicen, &
                  domain=fld%geom%G%Domain%mpp_domain)
          case default
             call log%warning("soca_fields:read_file: Not reading var "//fld%fldnames(ii))
          end select
       end do
       call restore_state(sis_restart, directory='')
       call restore_state(ocean_restart, directory='')
       call fms_io_exit()

       ! Indices for compute domain
       call geom_get_domain_indices(fld%geom, "compute", is, ie, js, je)

       ! Initialize mid-layer depth from layer thickness
       call fld%geom%thickness2depth(fld%hocn, fld%layer_depth)

       ! Compute mixed layer depth TODO: Move somewhere else ...
       do i = is, ie
          do j = js, je
             fld%mld(i,j) = soca_mld(fld%socn(i,j,:),&
                  &fld%tocn(i,j,:),&
                  &fld%layer_depth(i,j,:),&
                  &fld%geom%lon(i,j),&
                  &fld%geom%lat(i,j))
          end do
       end do

       ! Remap layers if needed
       if (vert_remap) then
          call initialize_remapping(remapCS,'PCM')
          do i = is, ie
             do j = js, je
                if (fld%geom%mask2d(i,j).eq.1) then
                   call remapping_core_h(remapCS, nz, h_common(i,j,:), fld%tocn(i,j,:),&
                        &nz, fld%hocn(i,j,:), fld%tocn(i,j,:))
                   call remapping_core_h(remapCS, nz, h_common(i,j,:), fld%socn(i,j,:),&
                        &nz, fld%hocn(i,j,:), fld%socn(i,j,:))

                else
                   fld%tocn(i,j,:) = 0.0_kind_real
                   fld%socn(i,j,:) = 0.0_kind_real
                end if
             end do
          end do
          fld%hocn = h_common
       end if
       call end_remapping(remapCS)

       ! Update halo
       call mpp_update_domains(fld%tocn, fld%geom%G%Domain%mpp_domain)
       call mpp_update_domains(fld%socn, fld%geom%G%Domain%mpp_domain)
       call mpp_update_domains(fld%ssh, fld%geom%G%Domain%mpp_domain)

       ! Set vdate if reading state
       if (iread==1) then
          sdate = config_get_string(c_conf,len(sdate),"date")
          call datetime_set(sdate, vdate)
       end if

       return
    end if

    ! Read diagnostic file
    if (iread==2) then
       ! Read ocean surface fields
       call fld%ocnsfc%read_diag(c_conf, fld%geom, fld%fldnames)

       incr_filename = config_get_string(c_conf,len(incr_filename),"filename")
       call fms_io_init()
       do ii = 1, fld%nf
          select case(fld%fldnames(ii))
             ! Ocean variables
          case ('ssh')
             call read_data(incr_filename,"ssh",fld%ssh(:,:),domain=fld%geom%G%Domain%mpp_domain)
          case ('tocn')
             call read_data(incr_filename,"temp",fld%tocn(:,:,:),domain=fld%geom%G%Domain%mpp_domain)
          case ('socn')
             call read_data(incr_filename,"salt",fld%socn(:,:,:),domain=fld%geom%G%Domain%mpp_domain)
          case ('hocn')
             call read_data(incr_filename,"h",fld%hocn(:,:,:),domain=fld%geom%G%Domain%mpp_domain)

             ! Sea-ice variables
          case ('cicen')
             call read_data(incr_filename, 'cicen', fld%cicen, domain=fld%geom%G%Domain%mpp_domain)
          case ('hicen')
             call read_data(incr_filename, 'hicen', fld%hicen, domain=fld%geom%G%Domain%mpp_domain)
          case default
             write(buf,*) 'soca_fields_mod::read_file::increment. Not reading '//fld%fldnames(ii)
             call log%info(buf,newl=.true.)

          end select
       end do
       call fms_io_exit()
    endif

    call check(fld)

  end subroutine read_file

  ! ------------------------------------------------------------------------------

  subroutine write_file(fld, c_conf, vdate)
    type(soca_field), intent(inout) :: fld    !< Fields
    type(c_ptr),         intent(in) :: c_conf !< Configuration
    type(datetime),   intent(inout) :: vdate  !< DateTime

    integer, parameter :: max_string_length=800    ! Yuk!
    character(len=max_string_length) :: filename
    character(len=1024):: buf
    integer :: ii

    call check(fld)

    !call geom_infotofile(fld%geom)
    call soca_write_restart(fld, c_conf, vdate)

!!$    filename = soca_genfilename(c_conf,max_string_length,vdate)
!!$    WRITE(buf,*) 'field:write_file: writing '//filename
!!$    call fckit_log%info(buf)
!!$
!!$    call soca_fld2file(fld, filename)

  end subroutine write_file

  ! ------------------------------------------------------------------------------

  subroutine gpnorm(fld, nf, pstat)
    type(soca_field),        intent(in) :: fld
    integer,                 intent(in) :: nf
    real(kind=kind_real), intent(inout) :: pstat(3, nf) !> [min, max, average]

    real(kind=kind_real) :: ocn_count, tmp(3)
    integer :: jj,  myrank, is, ie, js, je
    type(fckit_mpi_comm) :: f_comm

    ! Setup Communicator
    f_comm = fckit_mpi_comm()

    call check(fld)

    ! get domain bounds
    call geom_get_domain_indices(fld%geom, "compute", is, ie, js, je )

    ! get the total number of ocean grid cells
    tmp(1) = sum(fld%geom%mask2d(is:ie, js:je))
    call f_comm%allreduce(tmp(1), ocn_count, fckit_mpi_sum())

    ! calculate global min, max, mean for each field
    ! Note: The following code makes object oriented programmers cry a little.
    !  Most of the functions in this module should be rewritten to be
    !  agnostic to the actual number/names of variables, sigh.
    do jj=1, fld%nf
      tmp=0.0

      ! get local min/max/sum of each variable
      select case(fld%fldnames(jj))
      case("tocn")
        call fldinfo(fld%tocn(is:ie,js:je,:), tmp)
      case("socn")
        call fldinfo(fld%socn(is:ie,js:je,:), tmp)
      case("hocn")
        call fldinfo(fld%hocn(is:ie,js:je,:), tmp)
      case("ssh")
        call fldinfo(fld%ssh(is:ie,js:je), tmp)
      case("hicen")
        call fldinfo(fld%hicen(is:ie,js:je,:), tmp)
      case("cicen")
        call fldinfo(fld%cicen(is:ie,js:je,:), tmp)
      case("sw")
        call fldinfo(fld%ocnsfc%sw_rad(is:ie,js:je), tmp)
      case("lw")
        call fldinfo(fld%ocnsfc%lw_rad(is:ie,js:je), tmp)
      case("lhf")
        call fldinfo(fld%ocnsfc%latent_heat(is:ie,js:je), tmp)
      case("shf")
        call fldinfo(fld%ocnsfc%sens_heat(is:ie,js:je), tmp)
      case("us")
        call fldinfo(fld%ocnsfc%fric_vel(is:ie,js:je), tmp)
      end select

      ! calculate global min/max/mean
      call f_comm%allreduce(tmp(1), pstat(1,jj), fckit_mpi_min())
      call f_comm%allreduce(tmp(2), pstat(2,jj), fckit_mpi_max())
      call f_comm%allreduce(tmp(3), pstat(3,jj), fckit_mpi_sum())
      pstat(3,jj) = pstat(3,jj)/ocn_count

    end do
  end subroutine gpnorm

  ! ------------------------------------------------------------------------------

  subroutine fldrms(fld, prms)
    type(soca_field),      intent(in) :: fld
    real(kind=kind_real), intent(out) :: prms

    integer :: jf,jy,jx,ii
    real(kind=kind_real) :: zz, ns, n2dfld

    call check(fld)

    call dot_prod(fld,fld,prms) ! Global value
    prms=sqrt(prms)

  end subroutine fldrms

  ! ------------------------------------------------------------------------------

  subroutine ug_size(self, ug)
    type(soca_field),           intent(in) :: self
    type(unstructured_grid), intent(inout) :: ug

    integer :: isc, iec, jsc, jec
    integer :: igrid

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(self%geom, "compute", isc, iec, jsc, jec)

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

  subroutine ug_coord(self, ug)
    type(soca_field), intent(in) :: self
    type(unstructured_grid), intent(inout) :: ug

    integer :: igrid
    integer :: isc, iec, jsc, jec, jz

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(self%geom, "compute", isc, iec, jsc, jec)

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

  end subroutine ug_coord

  ! ------------------------------------------------------------------------------

  subroutine field_to_ug(self, ug, its)
    type(soca_field),           intent(in) :: self
    type(unstructured_grid), intent(inout) :: ug
    integer,                    intent(in) :: its

    integer :: isc, iec, jsc, jec, jk, incat, inzo, ncat, nzo, igrid
    integer :: ni, nj

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(self%geom, "compute", isc, iec, jsc, jec)

    ni = iec - isc + 1
    nj = jec - jsc + 1
    ncat = self%geom%ncat

    ! Define size
    call ug_size(self, ug)

    ! Allocate unstructured grid field
    call allocate_unstructured_grid_field(ug)
    ncat = self%geom%ncat
    nzo = ug%grid(1)%nl0

    ! Copy 3D field
    igrid = 1
    jk = 1
    ug%grid(igrid)%fld(:,:,:,its) = 0.0_kind_real

    ! tocn
    do inzo = 1, nzo
       ug%grid(igrid)%fld(1:ni*nj, inzo, jk, its) = &
            &reshape( self%tocn(isc:iec, jsc:jec,inzo), (/ug%grid(igrid)%nmga/) )
    end do
    jk = jk + 1

    ! socn
    do inzo = 1, nzo
       ug%grid(igrid)%fld(1:ni*nj, inzo, jk, its) = &
            &reshape( self%socn(isc:iec, jsc:jec,inzo), (/ug%grid(igrid)%nmga/) )
    end do
    jk = jk + 1

    if (ug%colocated==1) then
       ! 2D variables copied as 3D variables
       igrid = 1
    else
       ! 2D variables copied as 2D variables
       igrid = 2
       jk = 1
    end if

    ! cicen
    do incat = 1, ncat
       ug%grid(igrid)%fld(1:ni*nj, 1, jk, its) = &
            &reshape( self%cicen(isc:iec, jsc:jec, incat+1), (/ug%grid(igrid)%nmga/) )
       jk = jk + 1
    end do

    ! hicen
    do incat = 1, ncat
       ug%grid(igrid)%fld(1:ni*nj, 1, jk, its) = &
            &reshape( self%hicen(isc:iec, jsc:jec, incat), (/ug%grid(igrid)%nmga/) )
       jk = jk + 1
    end do

    ! ssh
    ug%grid(igrid)%fld(1:ni*nj, 1, jk, its) = &
         &reshape( self%ssh(isc:iec, jsc:jec), (/ug%grid(igrid)%nmga/) )
    jk = jk + 1

  end subroutine field_to_ug

  ! ------------------------------------------------------------------------------

  subroutine field_from_ug(self, ug, its)
    type(soca_field),     intent(inout) :: self
    type(unstructured_grid), intent(in) :: ug
    integer,                 intent(in) :: its

    integer :: isc, iec, jsc, jec, jk, incat, inzo, ncat, nzo, igrid
    integer :: ni, nj

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(self%geom, "compute", isc, iec, jsc, jec)

    ni = iec - isc + 1
    nj = jec - jsc + 1
    ncat = self%geom%ncat
    nzo = self%geom%nzo

    ! Copy 3D field
    igrid = 1
    jk = 1
    call zeros(self)

    ! tocn
    do inzo = 1, nzo
       self%tocn(isc:iec, jsc:jec,inzo) = &
            &reshape( ug%grid(igrid)%fld(1:ni*nj, inzo, jk, its), (/ni, nj/) )
    end do
    jk = jk + 1

    ! socn
    do inzo = 1, nzo
       self%socn(isc:iec, jsc:jec,inzo) = &
            &reshape( ug%grid(igrid)%fld(1:ni*nj, inzo, jk, its), (/ni, nj/) )
    end do
    jk = jk + 1

    if (ug%colocated==1) then
       ! 2D variables copied as 3D variables
       igrid = 1
    else
       ! 2D variables copied as 2D variables
       igrid = 2
       jk = 1
    end if

    ! cicen
    do incat = 1, ncat
       self%cicen(isc:iec, jsc:jec, incat+1) = &
            &reshape( ug%grid(igrid)%fld(1:ni*nj, 1, jk, its), (/ni, nj/))
       jk = jk + 1
    end do

    ! hicen
    do incat = 1, ncat
       self%hicen(isc:iec, jsc:jec, incat) = &
            &reshape( ug%grid(igrid)%fld(1:ni*nj, 1, jk, its), (/ni, nj/) )
       jk = jk + 1
    end do

    ! ssh
    self%ssh(isc:iec, jsc:jec) = &
         &reshape( ug%grid(igrid)%fld(1:ni*nj, 1, jk, its), (/ni, nj/) )
    jk = jk + 1

  end subroutine field_from_ug

  ! ------------------------------------------------------------------------------

  function common_vars(x1, x2)
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
    if (x1%geom%nzo /= x2%geom%nzo) call abor1_ftn("common_vars: error number of levels")
    !common_vars = x1%geom%nzi * common_vars

  end function common_vars

  ! ------------------------------------------------------------------------------

  subroutine check_resolution(x1, x2)
    type(soca_field), intent(in) :: x1, x2

    ! NEEDS WORK !!!
    if (x1%geom%nx /= x2%geom%nx .or.  &
         &x1%geom%ny /= x2%geom%ny ) then
       call abor1_ftn ("soca_fields: resolution error")
    endif
    call check(x1)
    call check(x2)

  end subroutine check_resolution

  ! ------------------------------------------------------------------------------

  subroutine check(self)
    type(soca_field), intent(in) :: self

    logical :: bad

    ! Doesn't do any thing ...
    bad = .false.
    !bad = bad .or. (size(self%cicen, 1) /= self%geom%nx)

    ! add more test here ...

    if (bad) then
       write(0,*)'nx, ny, nf, nzi, nzo = ',self%geom%nx,self%geom%ny
       call abor1_ftn ("soca_fields: field not consistent")
    endif

  end subroutine check

  ! ------------------------------------------------------------------------------
  !> Save soca fields to file using fms write_data
  subroutine soca_fld2file(fld, filename)
    type(soca_field),   intent(in) :: fld    !< Fields
    character(len=800), intent(in) :: filename

    integer :: ii
    character(len=1024):: buf
    character(len=800) :: fname

    fname = trim(filename)

    call check(fld)

    call fms_io_init()
    call set_domain( fld%geom%G%Domain%mpp_domain )
    do ii = 1, fld%nf
       select case(fld%fldnames(ii))

       case ('ssh')
          call write_data( fname, "ssh", fld%ssh, fld%geom%G%Domain%mpp_domain)
          call write_data( fname, "rossby_radius", fld%geom%rossby_radius, fld%geom%G%Domain%mpp_domain)
       case ('tocn')
          call write_data( fname, "temp", fld%tocn, fld%geom%G%Domain%mpp_domain)
       case ('socn')
          call write_data( fname, "salt", fld%socn, fld%geom%G%Domain%mpp_domain)
       case ('hocn')
          call write_data( fname, "h", fld%hocn, fld%geom%G%Domain%mpp_domain)
       case ('hicen')
          call write_data( fname, "hicen", fld%hicen, fld%geom%G%Domain%mpp_domain)
       case ('cicen')
          call write_data(fname, "cicen", fld%cicen, fld%geom%G%Domain%mpp_domain)
       case ('sw')
          call write_data(fname, "sw", fld%ocnsfc%sw_rad, fld%geom%G%Domain%mpp_domain)
       case ('lw')
          call write_data(fname, "lw", fld%ocnsfc%lw_rad, fld%geom%G%Domain%mpp_domain)
       case ('lhf')
          call write_data(fname, "lhf", fld%ocnsfc%latent_heat, fld%geom%G%Domain%mpp_domain)
       case ('shf')
          call write_data(fname, "shf", fld%ocnsfc%sens_heat, fld%geom%G%Domain%mpp_domain)
       case ('us')
          call write_data(fname, "us", fld%ocnsfc%fric_vel, fld%geom%G%Domain%mpp_domain)

       case default

       end select

    end do
    call fms_io_exit()

  end subroutine soca_fld2file

  ! ------------------------------------------------------------------------------
  !> Save soca fields in a restart format
  subroutine soca_write_restart(fld, c_conf, vdate)
    type(soca_field), intent(inout) :: fld      !< Fields
    type(c_ptr),         intent(in) :: c_conf   !< Configuration
    type(datetime),   intent(inout) :: vdate    !< DateTime

    integer, parameter :: max_string_length=800
    character(len=max_string_length) :: ocn_filename, ocnsfc_filename
    character(len=max_string_length) :: ice_filename, basename, incr_filename
    character(len=20) :: sdate
    character(len=1024)  :: buf
    integer :: iread, ii

    type(restart_file_type) :: ice_restart
    type(restart_file_type) :: ocean_restart
    type(restart_file_type) :: ocnsfc_restart
    integer :: idr, idr_ocean

    integer            :: nobs, nval, pe, ierror

    ! Generate file names
    ocn_filename = soca_genfilename(c_conf,max_string_length,vdate,"ocn")
    ice_filename = soca_genfilename(c_conf,max_string_length,vdate,"ice")
    ocnsfc_filename = soca_genfilename(c_conf,max_string_length,vdate,"sfc")

    call fms_io_init()
    ! Ocean State
    idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'ave_ssh', fld%ssh(:,:), &
         domain=fld%geom%G%Domain%mpp_domain)
    idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'Temp', fld%tocn(:,:,:), &
         domain=fld%geom%G%Domain%mpp_domain)
    idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'Salt', fld%socn(:,:,:), &
         domain=fld%geom%G%Domain%mpp_domain)
    idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'h', fld%hocn(:,:,:), &
         domain=fld%geom%G%Domain%mpp_domain)
    idr_ocean = register_restart_field(ocean_restart, ocn_filename, 'mld', fld%mld(:,:), &
         domain=fld%geom%G%Domain%mpp_domain)

    ! Sea-Ice
    idr = register_restart_field(ice_restart, ice_filename, 'part_size', fld%cicen, &
         domain=fld%geom%G%Domain%mpp_domain)
    idr = register_restart_field(ice_restart, ice_filename, 'h_ice', fld%hicen, &
         domain=fld%geom%G%Domain%mpp_domain)

    ! Ocean-surface
    idr = register_restart_field(ocnsfc_restart, ocnsfc_filename, &
                                'sw_rad', fld%ocnsfc%sw_rad, &
                                domain=fld%geom%G%Domain%mpp_domain)
    idr = register_restart_field(ocnsfc_restart, ocnsfc_filename, &
                                'lw_rad', fld%ocnsfc%lw_rad, &
                                domain=fld%geom%G%Domain%mpp_domain)
    idr = register_restart_field(ocnsfc_restart, ocnsfc_filename, &
                                'latent_heat', fld%ocnsfc%latent_heat, &
                                domain=fld%geom%G%Domain%mpp_domain)
    idr = register_restart_field(ocnsfc_restart, ocnsfc_filename, &
                                'sens_heat', fld%ocnsfc%sens_heat, &
                                domain=fld%geom%G%Domain%mpp_domain)
    idr = register_restart_field(ocnsfc_restart, ocnsfc_filename, &
                                'fric_vel', fld%ocnsfc%fric_vel, &
                                domain=fld%geom%G%Domain%mpp_domain)

    call save_restart(ocean_restart, directory='')
    call save_restart(ice_restart, directory='')
    call save_restart(ocnsfc_restart, directory='')
    call free_restart_type(ice_restart)
    call free_restart_type(ocean_restart)
    call free_restart_type(ocnsfc_restart)
    call fms_io_exit()

    return

  end subroutine soca_write_restart

  ! ------------------------------------------------------------------------------
  !> Generate filename (based on oops/qg)
  function soca_genfilename (c_conf,length,vdate, domain_type)
    type(c_ptr),                intent(in) :: c_conf
    integer,                    intent(in) :: length
    type(datetime),             intent(in) :: vdate
    character(len=3), optional, intent(in) :: domain_type

    character(len=length)                  :: soca_genfilename
    character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
         & prefix, mmb
    type(datetime) :: rdate
    type(duration) :: step
    integer lenfn

    fdbdir = config_get_string(c_conf,len(fdbdir),"datadir")
    expver = config_get_string(c_conf,len(expver),"exp")
    typ    = config_get_string(c_conf,len(typ)   ,"type")

    if (present(domain_type)) then
       expver = trim(domain_type)//"."//expver
    else
       expver = "ocn.ice."//expver
    end if
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
       soca_genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
    endif

    if (typ=="an") then
       call datetime_to_string(vdate, validitydate)
       lenfn = lenfn + 1 + LEN_TRIM(validitydate)
       soca_genfilename = TRIM(prefix) // "." // TRIM(validitydate)
    endif

    if (typ=="incr") then
       soca_genfilename = 'test-incr.nc'
    endif

    if (lenfn>length) &
         & call abor1_ftn("fields:genfilename: filename too long")

  end function soca_genfilename

  ! ------------------------------------------------------------------------------

end module soca_fields
