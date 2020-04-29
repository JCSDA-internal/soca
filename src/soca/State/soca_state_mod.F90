! (C) Copyright 2020-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_state_mod

use soca_geom_mod
use soca_fields_mod
use soca_increment_mod
use oops_variables_mod
use soca_geostrophy_mod
use kinds, only: kind_real
use fckit_log_module, only: log, fckit_log
use fms_mod,    only: read_data, set_domain
use fms_io_mod, only: fms_io_init, fms_io_exit
use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h 
use MOM_coms, only : max_across_PEs
use MOM_domains, only : pass_var, sum_across_PEs, root_PE
use mpp_mod, only     : mpp_broadcast, mpp_sync, mpp_sync_self
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use mpp_domains_mod, only  : mpp_global_field, mpp_get_global_domain
use horiz_interp_mod, only : horiz_interp_new, horiz_interp,horiz_interp_type
use horiz_interp_mod, only : horiz_interp_init, horiz_interp_del

implicit none
private

#include <MOM_memory.h>

type, public, extends(soca_fields) :: soca_state

contains

  ! constructors / destructors
  procedure :: create => soca_state_create

  ! interactions with increment
  procedure :: diff_incr=> soca_state_diff_incr
  procedure :: add_incr => soca_state_add_incr

  ! misc
  procedure :: rotate => soca_state_rotate
  procedure :: convert => soca_state_convert

end type

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine soca_state_create(self, geom, vars)
  class(soca_state),         intent(inout) :: self
  type(soca_geom),  pointer, intent(inout) :: geom
  type(oops_variables),      intent(inout) :: vars

  ! additional internal fields need to be created
  call vars%push_back("mld")
  call vars%push_back("layer_depth")

  ! continue with normal fields initialization
  call self%soca_fields%create(geom, vars)
end subroutine soca_state_create


! ------------------------------------------------------------------------------
!> Rotate horizontal vector
subroutine soca_state_rotate(self, coordinate, uvars, vvars)
  class(soca_state),  intent(inout) :: self
  character(len=*),      intent(in) :: coordinate ! "north" or "grid"
  type(oops_variables),  intent(in) :: uvars
  type(oops_variables),  intent(in) :: vvars

  integer :: z, i
  type(soca_field), pointer :: uocn, vocn
  real(kind=kind_real), allocatable :: un(:,:,:), vn(:,:,:)
  character(len=64) :: u_names, v_names

  do i=1, uvars%nvars()
    ! get (u, v) pair and make a copy
    u_names = trim(uvars%variable(i))
    v_names = trim(vvars%variable(i))
    if (self%has(u_names).and.self%has(v_names)) then
      call fckit_log%info("rotating "//trim(u_names)//" "//trim(v_names))
      call self%get(u_names, uocn)
      call self%get(v_names, vocn)
    else
      ! Skip if no pair found.
      call fckit_log%info("not rotating "//trim(u_names)//" "//trim(v_names))
      cycle
    end if
    allocate(un(size(uocn%val,1),size(uocn%val,2),size(uocn%val,3)))
    allocate(vn(size(uocn%val,1),size(uocn%val,2),size(uocn%val,3)))
    un = uocn%val
    vn = vocn%val

    select case(trim(coordinate))
    case("north")   ! rotate (uocn, vocn) to geo north
      do z=1,uocn%nz
        uocn%val(:,:,z) = &
        (self%geom%cos_rot(:,:)*un(:,:,z) + self%geom%sin_rot(:,:)*vn(:,:,z)) * uocn%mask(:,:)
        vocn%val(:,:,z) = &
        (- self%geom%sin_rot(:,:)*un(:,:,z) + self%geom%cos_rot(:,:)*vn(:,:,z)) * vocn%mask(:,:)
      end do
    case("grid")
      do z=1,uocn%nz
        uocn%val(:,:,z) = &
        (self%geom%cos_rot(:,:)*un(:,:,z) - self%geom%sin_rot(:,:)*vn(:,:,z)) * uocn%mask(:,:)
        vocn%val(:,:,z) = &
        (self%geom%sin_rot(:,:)*un(:,:,z) + self%geom%cos_rot(:,:)*vn(:,:,z)) * vocn%mask(:,:)
      end do
    end select
    deallocate(un, vn)

    ! update halos
    call uocn%update_halo(self%geom)
    call vocn%update_halo(self%geom)
  end do
end subroutine soca_state_rotate


! ------------------------------------------------------------------------------
!> add a set of increments to the set of fields
subroutine soca_state_add_incr(self, rhs)
  class(soca_state),  intent(inout) :: self
  class(soca_increment), intent(in) :: rhs

  integer, save :: cnt_outer = 1
  character(len=800) :: filename, str_cnt
  type(soca_field), pointer :: fld, fld_r
  integer :: i, k

  real(kind=kind_real) :: min_ice = 1e-6_kind_real
  real(kind=kind_real) :: amin = 1e-6_kind_real
  real(kind=kind_real) :: amax = 10.0_kind_real
  real(kind=kind_real), allocatable :: alpha(:,:), aice_bkg(:,:), aice_ana(:,:)
  type(soca_geostrophy_type) :: geostrophy
  type(soca_fields) :: incr
  type(soca_field), pointer :: t, s, u, v, h, dt, ds, du, dv


  ! make sure rhs is a subset of self
  call rhs%check_subset(self)

  ! Make a copy of the increment
  call incr%copy(rhs)

  ! Compute geostrophic increment
  ! TODO Move inside of the balance operator.
  !      Will need to be removed when assimilating ocean current or when using
  !      ensemble derived increments (ensemble cross-covariances that include currents)
  if (self%has('hocn').and.self%has('tocn').and.self%has('socn').and.&
      self%has('uocn').and.self%has('vocn')) then
     ! Get necessary background fields needed to compute geostrophic perturbations
     call self%get("tocn", t)
     call self%get("socn", s)
     call self%get("hocn", h)

     ! Make a copy of the increment and get the needed pointers
     !call incr_geo%copy(rhs)
     call incr%get("tocn", dt)
     call incr%get("socn", ds)
     call incr%get("uocn", du)
     call incr%get("vocn", dv)

     ! Compute the geostrophic increment
     call geostrophy%setup(self%geom, h%val)
     call geostrophy%tl(h%val, t%val, s%val,&
          dt%val, ds%val, du%val, dv%val, self%geom)
     call geostrophy%delete()
  end if

  ! Colocate increment fields with h-grid
  call incr%colocate('h')

  ! for each field that exists in incr, add to self
  do i=1,size(incr%fields)
    fld_r => incr%fields(i)
    call self%get(fld_r%name, fld)

    select case (fld%name)
    case ('cicen')
      ! NOTE: sea ice concentration is special

      ! compute background rescaling
      allocate(alpha, mold=fld_r%val(:,:,1))
      allocate(aice_bkg, mold=alpha)
      allocate(aice_ana, mold=alpha)
      aice_bkg  = sum(fld%val(:,:,:), dim=3)
      aice_ana  = aice_bkg + sum(fld_r%val(:,:,:), dim=3)
      where (aice_ana < 0.0_kind_real) aice_ana = 0.0_kind_real
      where (aice_ana > 1.0_kind_real) aice_ana = 1.0_kind_real
      alpha = 1.0_kind_real
      where ( aice_bkg > min_ice) alpha = aice_ana / aice_bkg

      ! limit size of increment
      where ( alpha > amax ) alpha = amax
      where ( alpha < amin ) alpha = amin

      ! add fraction of increment
      do k=1,fld%nz
        fld%val(:,:,k) = alpha * fld%val(:,:,k)
      end do

    case default
      ! everyone else is normal
      fld%val = fld%val + fld_r%val
    end select
  end do

  ! Save increment for outer loop cnt_outer
  write(str_cnt,*) cnt_outer
  filename='incr.'//adjustl(trim(str_cnt))//'.nc'

  ! Save increment and clean memory
  call incr%write_file(filename)
  call incr%delete()

  ! Update outer loop counter
  cnt_outer = cnt_outer + 1

end subroutine soca_state_add_incr


! ------------------------------------------------------------------------------
!> subtract two sets of fields, saving the results separately
subroutine soca_state_diff_incr(x1, x2, inc)
  class(soca_state),      intent(in)    :: x1
  class(soca_state),      intent(in)    :: x2
  class(soca_increment), intent(inout)  :: inc

  integer :: i
  type(soca_field), pointer :: f1, f2

  ! make sure fields correct shapes
  call inc%check_subset(x2)
  call x2%check_subset(x1)

  ! subtract
  do i=1,size(inc%fields)
    call x1%get(inc%fields(i)%name, f1)
    call x2%get(inc%fields(i)%name, f2)
    inc%fields(i)%val = f1%val - f2%val
  end do
end subroutine soca_state_diff_incr

!module for change of resolution
subroutine soca_state_convert(self, rhs)
  class(soca_state),  intent(inout) :: self
  class(soca_state), intent(in) :: rhs
  !
  integer :: isc1, iec1, jsc1, jec1
  integer :: isd1, ied1, jsd1, jed1
  integer :: isc2, iec2, jsc2, jec2
  integer :: i, j, k, n, nz, nz2
  integer :: iverbose
  real(kind=kind_real) :: missing
  real(kind=kind_real) :: missing_value = -1.e20
  real(kind=kind_real) :: roundoff = 1.e-3
  real(kind=kind_real), allocatable :: lon_in(:,:), lat_in(:,:)
  real(kind=kind_real), allocatable :: lon_out(:,:), lat_out(:,:)
  real(kind=kind_real), allocatable :: lonu_in(:,:), latu_in(:,:)
  real(kind=kind_real), allocatable :: lonu_out(:,:), latu_out(:,:) 
  real(kind=kind_real), allocatable :: lonv_in(:,:), latv_in(:,:)
  real(kind=kind_real), allocatable :: lonv_out(:,:), latv_out(:,:)
  real(kind=kind_real), allocatable :: mask2d(:,:,:), mask_ice(:,:,:)
  type(soca_field), pointer :: field1, field2, hocn, hocn2
  !
  character(len=256) :: remap_filename
  real(kind=kind_real), allocatable :: z_bot(:), h_common(:), tmp(:,:,:), bathy(:,:)
  type(remapping_CS)  :: remapCS
  integer :: tmp_nz, tmp_nz2, woa09_nz
  real(kind=kind_real), allocatable :: gdata(:,:,:)
  integer :: isg, ieg, jsg, jeg
  ! Thicknesses [m] that give level centers corresponding to table 2 of WOA09
  real(kind=kind_real), dimension(40) :: woa09_dz = (/ 5.,  10.,  10.,  15.,  22.5, 25., 25.,  25.,  &
                                                      37.5, 50.,  50.,  75., 100., 100., 100., 100., &
                                                     100., 100., 100., 100., 100., 100., 100., 175., &
                                                     250., 375., 500., 500., 500., 500., 500., 500., &
                                                     500., 500., 500., 500., 500., 500., 500., 500. /)

  ! Indices for compute domain for model1 (rhs)
  isc1 = rhs%geom%isc ; iec1 = rhs%geom%iec
  jsc1 = rhs%geom%jsc ; jec1 = rhs%geom%jec
  isd1 = rhs%geom%isd ; ied1 = rhs%geom%ied
  jsd1 = rhs%geom%jsc ; jed1 = rhs%geom%jed

  ! Indices for compute domain for model2 (self)
  isc2 = self%geom%isc ; iec2 = self%geom%iec
  jsc2 = self%geom%jsc ; jec2 = self%geom%jec

  ! allocate gloc for model1 domain (rhs)
  call mpp_get_global_domain(rhs%geom%Domain%mpp_domain, isg, ieg, jsg, jeg)
  allocate(lon_in(isg:ieg,jsg:jeg))
  allocate(lat_in(isg:ieg,jsg:jeg))
  call mpp_global_field (rhs%geom%Domain%mpp_domain, rhs%geom%lon(:,:), lon_in(:,:) )
  call mpp_global_field (rhs%geom%Domain%mpp_domain, rhs%geom%lat(:,:), lat_in(:,:) )

  !u Gloc for model1 (rhs)
  allocate(lonu_in(isg:ieg,jsg:jeg))
  allocate(latu_in(isg:ieg,jsg:jeg))
  call mpp_global_field (rhs%geom%Domain%mpp_domain, rhs%geom%lonu(:,:), lonu_in(:,:) )
  call mpp_global_field (rhs%geom%Domain%mpp_domain, rhs%geom%latu(:,:), latu_in(:,:) )

  !v Gloc for model1 (rhs)
  allocate(lonv_in(isg:ieg,jsg:jeg))
  allocate(latv_in(isg:ieg,jsg:jeg))
  call mpp_global_field (rhs%geom%Domain%mpp_domain, rhs%geom%lonv(:,:), lonv_in(:,:) )
  call mpp_global_field (rhs%geom%Domain%mpp_domain, rhs%geom%latv(:,:), latv_in(:,:) )

  !allocate gloc for model2 domain (self)
  allocate(lon_out(SZI_(self%geom),SZJ_(self%geom)), lat_out(SZI_(self%geom),SZJ_(self%geom)))
  lon_out(:,:) = self%geom%lon(:,:)
  lat_out(:,:) = self%geom%lat(:,:)

  !u Gloc for model2 (self)
  allocate(lonu_out(SZI_(self%geom),SZJ_(self%geom)), latu_out(SZI_(self%geom),SZJ_(self%geom)))
  lonu_out(:,:) = self%geom%lonu(:,:)
  latu_out(:,:) = self%geom%latu(:,:)

  !v Gloc for model2 (self)
  allocate(lonv_out(SZI_(self%geom),SZJ_(self%geom)), latv_out(SZI_(self%geom),SZJ_(self%geom)))
  lonv_out(:,:) = self%geom%lonv(:,:)
  latv_out(:,:) = self%geom%latv(:,:)

  ! Get vertical coordinate layer for model1 (rhs)
  call rhs%get("hocn", hocn)
  nz = hocn%nz
  ! Get vertical coordinate layer for model2 (self)
  call self%get("hocn", hocn2)
  nz2 = hocn2%nz

  if (rhs%change_resol > 0 .and. allocated(rhs%vertical_grid_fixz)) then !fix-z case

    allocate(h_common(1:nz2))     
    ! Get user-defined vertical cell thickness 
    if (trim(rhs%vertical_grid_fixz)=="woa09") then
      woa09_nz = size(woa09_dz,1)
      if (woa09_nz /= nz2) call MOM_error(FATAL, "woa09: number of vertical layers in model2 needs to be 40!")
      h_common(:) = woa09_dz(:)
    else
      remap_filename = trim(rhs%vertical_grid_fixz)
      call fms_io_init()
      call read_data(remap_filename, 'dz', h_common,domain=self%geom%Domain%mpp_domain)
      call fms_io_exit()
    end if  ! rhs%vertical_grid_fixz     

    ! Read common vertical coordinate from file
    allocate(z_bot(1:nz2))
    allocate(bathy(SZI_(self%geom),SZJ_(self%geom)))
    allocate(mask2d(SZI_(self%geom),SZJ_(self%geom),1:nz2))

    !
    z_bot(1) = h_common(1)
    do k=2,nz2
      z_bot(k) = z_bot(k-1)+h_common(k)
    end do

    !read bathy
    if (allocated(rhs%bathy_file)) then
      remap_filename = trim(rhs%bathy_file)
      call fms_io_init()
      call read_data(remap_filename, 'depth', bathy, domain=self%geom%Domain%mpp_domain)
      call fms_io_exit()
    else
      call MOM_error(FATAL, "Need to define bathy_file in yaml file if change_resol > 0")  
    end if

    ! constructure new h for model2 
    mask2d = 0.d0
    do i = isc2, iec2
      do j = jsc2, jec2
        if (hocn2%mask(i,j) < 1.d0 ) cycle
        do k = 1,nz2
          if (k==1) then
            hocn2%val(i,j,k) = h_common(k) 
            mask2d(i,j,k) = 1.d0 
          elseif (k>1 .and. z_bot(k) <= bathy(i,j)) then
            hocn2%val(i,j,k) = h_common(k)
            mask2d(i,j,k) = 1.d0
          end if
        end do ! k
      end do !j
    end do !i

    ! Initialize vertical remapping 
    call initialize_remapping(remapCS,'PPM_IH4')

  end if  ! change_resol>0 

  !
  allocate(gdata(isg:ieg,jsg:jeg,1:nz2))
  allocate(tmp(SZI_(rhs%geom),SZJ_(rhs%geom),1:nz2))

  ! loop through fields
  do n=1,size(rhs%fields)
    field1 => rhs%fields(n)
    select case(field1%name)

    !
    case ('ssh')
    call self%get(trim(field1%name),field2)
    tmp = 0.d0
    tmp(:,:,1:field1%nz) = field1%val(:,:,1:field1%nz) !2D only case
    ! get global data field for model1 
    call mpp_global_field (rhs%geom%Domain%mpp_domain, tmp(:,:,1:field1%nz), gdata(:,:,1:field1%nz) )
    if (rhs%change_resol==1) field2%val(:,:,1:field2%nz) = tmp(:,:,1:field2%nz)*mask2d(:,:,1:field2%nz)
    ! define missing value and do horizontal remapping
    missing = 0.0d0
    if (rhs%change_resol==2) &
            call soca_hinterp(self,field2,gdata,mask2d,field1%nz,missing,lon_in,lat_in,lon_out,lat_out)

    !
    case ('tocn','socn')
    call self%get(trim(field1%name),field2)
    tmp = 0.d0
    if (rhs%change_resol > 0) then
      do i = isc1, iec1
        do j = jsc1, jec1
          if (field1%mask(i,j) .eq. 0.d0) cycle
          tmp_nz = count(hocn%val(i,j,:) > 0.001d0)
          call remapping_core_h(remapCS, tmp_nz, hocn%val(i,j,1:tmp_nz), field1%val(i,j,1:tmp_nz),&
                     nz2, h_common, tmp(i,j,:))
          forall (k=1:nz2, sum(h_common(1:k)) .gt. sum(hocn%val(i,j,1:tmp_nz))) tmp(i,j,k) = 0.d0
        end do !i
      end do !j
    end if !vert_remap

    ! get global field for model1
    call mpp_global_field (rhs%geom%Domain%mpp_domain, tmp(:,:,1:nz2), gdata(:,:,1:nz2) )
    if (rhs%change_resol==1) field2%val(:,:,1:field2%nz) = tmp(:,:,1:field2%nz)*mask2d(:,:,1:field2%nz)
    ! define missing value and do horizontal remapping
    missing = 0.0d0
    if (rhs%change_resol==2) &
            call soca_hinterp(self,field2,gdata,mask2d,field2%nz,missing,lon_in,lat_in,lon_out,lat_out)

     !
    case ('uocn')
    call self%get(trim(field1%name),field2)
    tmp = 0.d0
    if (rhs%change_resol > 0) then
      do i = isc1, iec1
        do j = jsc1, jec1
          if (field1%mask(i,j) .eq. 0.d0) cycle
          tmp_nz = count(hocn%val(i,j,:) > 0.001d0)
          call remapping_core_h(remapCS, tmp_nz, hocn%val(i,j,1:tmp_nz), field1%val(i,j,1:tmp_nz),&
                     nz2, h_common, tmp(i,j,:))
          forall (k=1:nz2, sum(h_common(1:k)) .gt. sum(hocn%val(i,j,1:tmp_nz))) tmp(i,j,k) = 0.d0
        end do !i
      end do !j
    end if !vert_remap

    ! get global field for model1
    call mpp_global_field (rhs%geom%Domain%mpp_domain, tmp(:,:,1:nz2), gdata(:,:,1:nz2) )
    if (rhs%change_resol==1) field2%val(:,:,1:field2%nz) = tmp(:,:,1:field2%nz)*mask2d(:,:,1:field2%nz)
    ! define missing value and do horizontal remapping
    missing = 0.0d0
    if (rhs%change_resol==2) &
            call soca_hinterp(self,field2,gdata,mask2d,field2%nz,missing,lonu_in,latu_in,lonu_out,latu_out)
   
    !
    case ('vocn')
    call self%get(trim(field1%name),field2)
    tmp = 0.d0
    if (rhs%change_resol > 0) then
      do i = isc1, iec1
        do j = jsc1, jec1
          if (field1%mask(i,j) .eq. 0.d0) cycle
          tmp_nz = count(hocn%val(i,j,:) > 0.001d0)
          call remapping_core_h(remapCS, tmp_nz, hocn%val(i,j,1:tmp_nz), field1%val(i,j,1:tmp_nz),&
                     nz2, h_common, tmp(i,j,:))
          forall (k=1:nz2, sum(h_common(1:k)) .gt. sum(hocn%val(i,j,1:tmp_nz))) tmp(i,j,k) = 0.d0
        end do !i
      end do !j
    end if !vert_remap

    ! get global field for model1
    call mpp_global_field (rhs%geom%Domain%mpp_domain, tmp(:,:,1:nz2), gdata(:,:,1:nz2) )
    if (rhs%change_resol==1) field2%val(:,:,1:field2%nz) = tmp(:,:,1:field2%nz)*mask2d(:,:,1:field2%nz)
    ! defie missing value and do horizontal remapping
    missing = 0.0d0
    if (rhs%change_resol==2) &
            call soca_hinterp(self,field2,gdata,mask2d,field2%nz,missing,lonv_in,latv_in,lonv_out,latv_out)

    case ('cicen')
    call self%get(trim(field1%name),field2)
    tmp = 0.d0
    tmp(:,:,1:field1%nz) = field1%val(:,:,1:field1%nz) 
    call mpp_global_field (rhs%geom%Domain%mpp_domain, tmp(:,:,1:field1%nz), gdata(:,:,1:field1%nz) )
    if (rhs%change_resol==1) field2%val(:,:,1:field2%nz) = tmp(:,:,1:field2%nz)*mask2d(:,:,1:field2%nz)
    
    !ice_remapping
    !TODO more sophiscated mask for ice may be needed.
    allocate(mask_ice(SZI_(self%geom),SZJ_(self%geom),1:field2%nz))
    mask_ice = 1.0d0 !tmp fix
    !
    missing = missing_value
    if (rhs%change_resol==2) &
            call soca_hinterp(self,field2,gdata,mask_ice,field2%nz,missing,lonu_in,latv_in,lonu_out,latv_out) 

    end select
  end do ! n

end subroutine soca_state_convert 

subroutine fill_miss_2d(aout, good, fill, prev, G, smooth, num_pass, relc, crit, debug, answers_2018)
  use MOM_coms, only : sum_across_PEs

  class(soca_state),         intent(inout) :: G
  real(kind=kind_real), dimension(SZI_(G%geom),SZJ_(G%geom)), intent(inout) :: aout
  real(kind=kind_real), dimension(SZI_(G%geom),SZJ_(G%geom)), intent(in) :: good,fill
  real(kind=kind_real), dimension(SZI_(G%geom),SZJ_(G%geom)), optional, intent(in) :: prev
  logical,     optional, intent(in)    :: smooth ! If present and true, apply a number of
                                                 ! Laplan iterations to the interpolated data
  integer,     optional, intent(in)    :: num_pass ! The maximum number of iterations
  real(kind=kind_real),        optional, intent(in)    :: relc ! A relaxation coefficient for Laplacian (ND)
  real(kind=kind_real),        optional, intent(in)    :: crit ! A minimal value for deltas between iterations.
  logical,     optional, intent(in)    :: debug ! If true, write verbose debugging messages.
  logical,     optional, intent(in)    :: answers_2018 ! If true, use expressions that give the same
                                                ! answers as the code did in late 2018.  Otherwise
                                                ! add parentheses for rotational symmetry.

  real(kind=kind_real), dimension(SZI_(G%geom),SZJ_(G%geom)) :: b,r
  real(kind=kind_real), dimension(SZI_(G%geom),SZJ_(G%geom)) :: fill_pts, good_, good_new

  character(len=256) :: mesg  ! The text of an error message
  integer :: i,j,k
  real(kind=kind_real)    :: east,west,north,south,sor
  real(kind=kind_real)    :: ge,gw,gn,gs,ngood
  logical :: do_smooth,siena_bug
  real(kind=kind_real)    :: nfill, nfill_prev
  integer, parameter :: num_pass_default = 10000
  real(kind=kind_real), parameter :: relc_default = 0.25, crit_default = 1.e-3

  integer :: npass
  integer :: is, ie, js, je
  real(kind=kind_real)    :: relax_coeff, acrit, ares
  logical :: debug_it, ans_2018

  debug_it=.false.
  if (PRESENT(debug)) debug_it=debug

  is = G%geom%isc ; ie = G%geom%iec ; js = G%geom%jsc ; je = G%geom%jec

  npass = num_pass_default
  if (PRESENT(num_pass)) npass = num_pass

  relax_coeff = relc_default
  if (PRESENT(relc)) relax_coeff = relc

  acrit = crit_default
  if (PRESENT(crit)) acrit = crit

  do_smooth=.false.
  if (PRESENT(smooth)) do_smooth=smooth

  ans_2018 = .true. ; if (PRESENT(answers_2018)) ans_2018 = answers_2018

  fill_pts(:,:) = fill(:,:)

  nfill = sum(fill(is:ie,js:je))
  call sum_across_PEs(nfill)

  nfill_prev = nfill
  good_(:,:) = good(:,:)
  r(:,:) = 0.0

  do while (nfill > 0.0)

    call pass_var(good_, G%geom%Domain)  
    call pass_var(aout, G%geom%Domain)   

    b(:,:)=aout(:,:)
    good_new(:,:)=good_(:,:)

    do j=js,je ; do i=is,ie

      if (good_(i,j) == 1.0 .or. fill(i,j) == 0.) cycle

      ge=good_(i+1,j) ; gw=good_(i-1,j)
      gn=good_(i,j+1) ; gs=good_(i,j-1)
      east=0.0 ; west=0.0 ; north=0.0 ; south=0.0
      if (ge == 1.0) east = aout(i+1,j)*ge
      if (gw == 1.0) west = aout(i-1,j)*gw
      if (gn == 1.0) north = aout(i,j+1)*gn
      if (gs == 1.0) south = aout(i,j-1)*gs

      if (ans_2018) then
        ngood = ge+gw+gn+gs
      else
        ngood = (ge+gw) + (gn+gs)
      endif
      if (ngood > 0.) then
        if (ans_2018) then
          b(i,j)=(east+west+north+south)/ngood
        else
          b(i,j) = ((east+west) + (north+south))/ngood
        endif
        fill_pts(i,j) = 0.0
        good_new(i,j) = 1.0
      endif
    enddo ; enddo

    aout(is:ie,js:je) = b(is:ie,js:je)
    good_(is:ie,js:je) = good_new(is:ie,js:je)
    nfill_prev = nfill
    nfill = sum(fill_pts(is:ie,js:je))
    call sum_across_PEs(nfill)

    if (nfill == nfill_prev .and. PRESENT(prev)) then
      do j=js,je ; do i=is,ie ; if (fill_pts(i,j) == 1.0) then
        aout(i,j) = prev(i,j)
        fill_pts(i,j) = 0.0
      endif ; enddo ; enddo
    elseif (nfill == nfill_prev) then
      call MOM_error(WARNING, &
           'Unable to fill missing points using either data at the same vertical level from a connected basin'//&
           'or using a point from a previous vertical level.  Make sure that the original data has some valid'//&
           'data in all basins.', .true.)
      write(mesg,*) 'nfill=',nfill
      call MOM_error(WARNING, mesg, .true.)
    endif

    nfill = sum(fill_pts(is:ie,js:je))
    call sum_across_PEs(nfill)  

  enddo

  if (do_smooth) then ; do k=1,npass
    call pass_var(aout,G%geom%Domain)
    do j=js,je ; do i=is,ie
      if (fill(i,j) == 1) then
        east = max(good(i+1,j),fill(i+1,j)) ; west = max(good(i-1,j),fill(i-1,j))
        north = max(good(i,j+1),fill(i,j+1)) ; south = max(good(i,j-1),fill(i,j-1))
        if (ans_2018) then
          r(i,j) = relax_coeff*(south*aout(i,j-1)+north*aout(i,j+1) + &
                                west*aout(i-1,j)+east*aout(i+1,j) - &
                               (south+north+west+east)*aout(i,j))
        else
          r(i,j) = relax_coeff*( ((south*aout(i,j-1) + north*aout(i,j+1)) + &
                                  (west*aout(i-1,j)+east*aout(i+1,j))) - &
                                 ((south+north)+(west+east))*aout(i,j) )
        endif
      else
        r(i,j) = 0.
      endif
    enddo ; enddo
    ares = 0.0
    do j=js,je ; do i=is,ie
      aout(i,j) = r(i,j) + aout(i,j)
      ares = max(ares, abs(r(i,j)))
    enddo ; enddo
    call max_across_PEs(ares) 
    if (ares <= acrit) exit
  enddo ; endif

  do j=js,je ; do i=is,ie
    if (good_(i,j) == 0.0 .and. fill_pts(i,j) == 1.0) then
      write(mesg,*) 'In fill_miss, fill, good,i,j=',fill_pts(i,j),good_(i,j),i,j
      call MOM_error(WARNING, mesg, .true.)
      call MOM_error(FATAL, "MOM_initialize: "// &
           "fill is true and good is false after fill_miss, how did this happen?")
    endif
 enddo ; enddo

end subroutine fill_miss_2d
  
subroutine soca_hinterp(self,field2,gdata,mask2d,nz,missing,lon_in,lat_in,lon_out,lat_out)
  class(soca_state),  intent(inout) :: self
  type(soca_field), pointer, intent(inout) :: field2
  real(kind=kind_real), dimension(:,:,:), intent(in) :: gdata
  real(kind=kind_real), dimension(SZI_(self%geom),SZJ_(self%geom),1:field2%nz), intent(in) :: mask2d
  integer, intent(in) :: nz
  real(kind=kind_real), intent(in) :: missing
  real(kind=kind_real), dimension(:,:), intent(in) :: lon_in, lat_in 
  real(kind=kind_real), dimension(SZI_(self%geom),SZJ_(self%geom)), intent(in) :: lon_out, lat_out 

  !local variables
  integer :: i, j, k, isg, ieg, jsg, jeg
  integer :: isc2, iec2, jsc2, jec2
  real(kind=kind_real) :: missing_value = -1.e20
  real(kind=kind_real) :: roundoff = 1.e-3
  real(kind=kind_real) :: PI_180
  type(horiz_interp_type) :: Interp
  real(kind_real), dimension(:,:), allocatable :: mask_in_, tr_inp
  real(kind_real), dimension(:,:), allocatable :: tr_out, fill, good, prev, mask_out_, tr_outf
  !
  integer  :: jeg1
  real(kind=kind_real) :: max_lat,pole,npole
  real(kind=kind_real) :: min_lat
  real(kind=kind_real), dimension(:), allocatable :: last_row, first_row 
  real(kind=kind_real), dimension(:,:), allocatable :: lon_inp, lat_inp, wild, tr_in 
  logical :: add_np, add_sp
  !
  PI_180=atan(1.0d0)/45.0d0 

  !
  allocate(tr_out(SZI_(self%geom),SZJ_(self%geom)))
  allocate(tr_outf(SZI_(self%geom),SZJ_(self%geom)))
  allocate(fill(SZI_(self%geom),SZJ_(self%geom)))
  allocate(good(SZI_(self%geom),SZJ_(self%geom)))
  allocate(prev(SZI_(self%geom),SZJ_(self%geom)))
  allocate(mask_out_(SZI_(self%geom),SZJ_(self%geom)))
 
  !
  isg = 1; jsg = 1;
  ieg = size(gdata,1); jeg = size(gdata,2)

  ! Indices for compute domain for regional model
  isc2 = self%geom%isc ; iec2 = self%geom%iec
  jsc2 = self%geom%jsc ; jec2 = self%geom%jec

  ! extrapolate the input data to the north pole using the northerm-most latitude
  max_lat = maxval(lat_in(:,jeg))
  add_np=.false.
  if (max_lat < 90.d0) then
    add_np=.true.
    jeg1=jeg+1
    allocate(lat_inp(isg:ieg,jsg:jeg1))
    allocate(lon_inp(isg:ieg,jsg:jeg1))
    lat_inp(isg:ieg,jsg:jeg)=lat_in(isg:ieg,jsg:jeg)
    lat_inp(isg:ieg,jeg1)=90.d0
    lon_inp(isg:ieg,jsg:jeg)=lon_in(isg:ieg,jsg:jeg)
    lon_inp(isg:ieg,jeg1)=lon_in(isg:ieg,jeg)
  else
    jeg1=jeg
    allocate(lat_inp(isg:ieg,jsg:jeg))
    allocate(lon_inp(isg:ieg,jsg:jeg))
    lat_inp(isg:ieg,jsg:jeg)=lat_in(isg:ieg,jsg:jeg)    
    lon_inp(isg:ieg,jsg:jeg)=lon_in(isg:ieg,jsg:jeg) 
  endif

  !TODO do same thing for the south pole
  !Right now simply repeat first row if min_lat < -90
  min_lat = minval(lat_in(:,jsg))
  add_sp=.false.
  if (min_lat > -90.d0) then
    add_sp=.true.
    jeg1=jeg1+1
    allocate(wild(isg:ieg,jsg:jeg1-1))
    wild(:,:) = lat_inp(:,:)
    if (allocated(lat_inp)) deallocate(lat_inp)
    allocate(lat_inp(isg:ieg,jsg:jeg1))
    lat_inp(isg:ieg,jsg+1:jeg1)=wild(isg:ieg,jsg:jeg1-1)
    lat_inp(isg:ieg,jsg) = -90.d0 
    !
    wild(:,:) = lon_inp
    if (allocated(lon_inp)) deallocate(lon_inp)
    allocate(lon_inp(isg:ieg,jsg:jeg1))
    lon_inp(isg:ieg,jsg+1:jeg1)=lon_in(isg:ieg,jsg:jeg1-1)
    lon_inp(isg:ieg,jsg)=lon_in(isg:ieg,jsg)
  endif

  !
  allocate(tr_in(isg:ieg,jsg:jeg))
  allocate(tr_inp(isg:ieg,jsg:jeg1))
  allocate(mask_in_(isg:ieg,jsg:jeg1))
  allocate(last_row(isg:ieg))

  do k = 1, nz
    ! extrapolate the input data to the north pole using the northerm-most latitude 
    if (is_root_pe()) then
      tr_in(isg:ieg,jsg:jeg) = gdata(isg:ieg,jsg:jeg,k)

      if (add_np) then
        last_row(:)=tr_in(:,jeg); pole=0.d0; npole=0.d0
        do i=isg,ieg
          if (abs(tr_in(i,jeg)-missing) > abs(roundoff)) then
               pole = pole+last_row(i)
               npole = npole+1.d0
          endif
        enddo
        if (npole > 0) then
            pole=pole/npole
        else
            pole=missing_value
        endif
        if (add_sp) then
          tr_inp(:,jsg) = tr_in(:,jsg)      
          tr_inp(:,jsg+1:jeg1-1) = tr_in(:,:) 
          tr_inp(:,jeg1) = pole       
        else
          tr_inp(:,jsg:jeg) = tr_in(:,:)
          tr_inp(:,jeg1) = pole
        endif !add_sp
      else
         tr_inp(isg:ieg,jsg:jeg) = tr_in(isg:ieg,jsg:jeg)     
      endif !add_np    
    
    end if !root_pe 

    call mpp_sync()
    call mpp_broadcast(tr_inp, ieg*jeg1, root_PE())
    call mpp_sync_self()


    mask_in_=0.d0

    do j=jsg,jeg1 ; do i=isg,ieg
      if (abs(tr_inp(i,j)-missing) > abs(roundoff)) then
        mask_in_(i,j)=1.d0
      else
        tr_inp(i,j) = missing_value
      endif
    enddo ; enddo

    tr_out(:,:) = 0.d0

    ! initialize horizontal remapping 
    if (k==1) call horiz_interp_new(Interp, lon_inp(:,:)*PI_180, lat_inp(:,:)*PI_180, lon_out(isc2:iec2,jsc2:jec2)*PI_180, &
       lat_out(isc2:iec2,jsc2:jec2)*PI_180, interp_method='bilinear', src_modulo=.true.)

    call horiz_interp(Interp,tr_inp,tr_out(isc2:iec2,jsc2:jec2),&
         missing_value=missing_value, new_missing_handle=.true.)

    mask_out_=0.d0
    do j=jsc2,jec2
      do i=isc2,iec2
        if (abs(tr_out(i,j)-missing_value) > abs(roundoff) .and. mask2d(i,j,k) == 1.d0) &
          mask_out_(i,j)=1.d0
      enddo
    enddo

    fill(:,:) = 0.d0
    good(:,:) = 0.d0
    do j=jsc2,jec2
      do i=isc2,iec2
        if (mask_out_(i,j) < 1.0d0) then
          tr_out(i,j) = missing_value
        else
          good(i,j) = 1.0d0
        endif
          if (mask2d(i,j,k) == 1.d0 .and. mask_out_(i,j) < 1.0d0) fill(i,j)=1.d0
      end do !i
    end do !j
    call pass_var(fill, self%geom%Domain)
    call pass_var(good, self%geom%Domain)

    tr_outf(:,:) = tr_out(:,:)
    if (k==1) prev(:,:) = tr_outf(:,:)

    call fill_miss_2d(tr_outf, good, fill, prev=prev, G=self, &
         smooth=.true.,answers_2018=.true.)

    field2%val(:,:,k) = tr_outf(:,:)*mask2d(:,:,k)
    prev(:,:) = field2%val(:,:,k)

  end do

end subroutine soca_hinterp


end module
