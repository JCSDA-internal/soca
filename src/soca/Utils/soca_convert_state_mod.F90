! (C) Copyright 2020-2020 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! ------------------------------------------------------------------------------

module soca_convert_state_mod

  use soca_geom_mod
  use soca_fields_mod
  use soca_utils, only: soca_remap_idw
  use kinds, only: kind_real
  use fms_io_mod, only: read_data, write_data, fms_io_init, fms_io_exit
  use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h
  use MOM_domains, only : pass_var, root_PE, sum_across_pes
  use mpp_mod, only     : mpp_broadcast, mpp_sync, mpp_sync_self
  use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
  use mpp_domains_mod, only  : mpp_global_field, mpp_update_domains
  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_type
  use MOM_horizontal_regridding, only : meshgrid, fill_miss_2d
  use MOM_grid, only : ocean_grid_type
  use fckit_exception_module, only: fckit_exception

  implicit none
  private

  type, public :: soca_convertstate_type
     real(kind=kind_real), allocatable, dimension(:,:,:) :: hocn_src, hocn_des

   contains
     procedure :: setup => soca_convertstate_setup
     procedure :: change_resol => soca_convertstate_change_resol
     procedure :: change_resol2d => soca_convertstate_change_resol2d
     procedure :: clean => soca_convertstate_delete
  end type soca_convertstate_type

! ------------------------------------------------------------------------------
contains

subroutine soca_convertstate_setup(self, src, des, hocn, hocn2)
  class(soca_convertstate_type), intent(inout) :: self
  type(soca_geom),              intent(inout) :: src, des
  type(soca_field),             intent(inout) :: hocn, hocn2

  !local
  integer :: tmp(1)

  call fms_io_init()

  call read_data(trim(src%geom_grid_file), 'nzo_zstar', tmp(1), domain=src%Domain%mpp_domain)
  src%nzo_zstar = tmp(1)

  call read_data(trim(des%geom_grid_file), 'nzo_zstar', tmp(1), domain=des%Domain%mpp_domain)
  des%nzo_zstar = tmp(1)

  if (des%nzo_zstar /= src%nzo_zstar) call fckit_exception%abort(&
     "target nzo_zstar /= source nzo_zstar! Reset maximum depth in target grid MOM_input file and re-run soca gridgen")


  if (allocated(src%h_zstar)) deallocate(src%h_zstar)
  allocate(src%h_zstar(src%isd:src%ied,src%jsd:src%jed,1:src%nzo_zstar))
  call read_data(trim(src%geom_grid_file), 'h_zstar', src%h_zstar, domain=src%Domain%mpp_domain)

  if (allocated(des%h_zstar)) deallocate(des%h_zstar)
  allocate(des%h_zstar(des%isd:des%ied,des%jsd:des%jed,1:des%nzo_zstar))
  call read_data(trim(des%geom_grid_file), 'h_zstar', des%h_zstar, domain=des%Domain%mpp_domain)

  call fms_io_exit()

  allocate(self%hocn_src(src%isd:src%ied,src%jsd:src%jed,1:src%nzo))
  allocate(self%hocn_des(des%isd:des%ied,des%jsd:des%jed,1:des%nzo))

  ! set hocn for target grid
  hocn2%val = des%h
  self%hocn_src = hocn%val
  self%hocn_des = hocn2%val

end subroutine soca_convertstate_setup

! ------------------------------------------------------------------------------
!> Cleanup
  subroutine soca_convertstate_delete(self)
    class(soca_convertstate_type), intent(inout) :: self

    deallocate(self%hocn_src)
    deallocate(self%hocn_des)

  end subroutine soca_convertstate_delete

! ------------------------------------------------------------------------------
subroutine soca_convertstate_change_resol2d(self, field_src, field_des, geom_src, geom_des)
  class(soca_convertstate_type),  intent(inout) :: self
  type(soca_field), pointer,      intent(inout) :: field_src, field_des
  type(soca_geom),                intent(inout) :: geom_src, geom_des

  !local
  integer :: i, j, k, tmp_nz, nz_
  integer :: isc1, iec1, jsc1, jec1, isd1, ied1, jsd1, jed1, isg, ieg, jsg, jeg
  integer :: isc2, iec2, jsc2, jec2, isd2, ied2, jsd2, jed2
  type(remapping_CS)  :: remapCS2
  type(horiz_interp_type) :: Interp
  real(kind=kind_real) :: missing = 0.d0
  real(kind=kind_real) :: z_tot
  real(kind=kind_real), dimension(geom_src%isg:geom_src%ieg) :: lon_in
  real(kind=kind_real), dimension(geom_src%jsg:geom_src%jeg) :: lat_in
  real(kind=kind_real), dimension(geom_des%isd:geom_des%ied,geom_des%jsd:geom_des%jed) :: lon_out, lat_out
  real(kind=kind_real), dimension(geom_des%isd:geom_des%ied,geom_des%jsd:geom_des%jed) :: mask_
  real(kind=kind_real), allocatable :: tmp(:,:,:), tmp2(:,:,:), gdata(:,:,:)

  ! Indices for compute, data, and global domain for source
  isc1 = geom_src%isc ; iec1 = geom_src%iec ; jsc1 = geom_src%jsc ; jec1 = geom_src%jec
  isd1 = geom_src%isd ; ied1 = geom_src%ied ; jsd1 = geom_src%jsd ; jed1 = geom_src%jed
  isg = geom_src%isg ; ieg = geom_src%ieg ; jsg = geom_src%jsg ; jeg = geom_src%jeg

  ! Indices for compute and data domain for des
  isc2 = geom_des%isc ; iec2 = geom_des%iec ; jsc2 = geom_des%jsc ; jec2 = geom_des%jec
  isd2 = geom_des%isd ; ied2 = geom_des%ied ; jsd2 = geom_des%jsd ; jed2 = geom_des%jed

  lon_in = geom_src%lonh ; lat_in = geom_src%lath
  if (field_src%name == "uocn" .and. field_des%name == "uocn") lon_in = geom_src%lonq
  if (field_src%name == "vocn" .and. field_des%name == "vocn") lat_in = geom_src%latq

  ! Initialize work arrays
  nz_ = field_src%nz
  allocate(tmp(isd1:ied1,jsd1:jed1,1:nz_),gdata(isg:ieg,jsg:jeg,1:nz_),tmp2(isd2:ied2,jsd2:jed2,1:nz_))
  tmp = 0.d0 ; gdata = 0.d0 ; tmp2 = 0.d0;
  tmp(:,:,1:nz_) = field_src%val(:,:,1:nz_)

  ! Reconstruct global input field
  call mpp_update_domains(tmp, geom_src%Domain%mpp_domain)
  mask_ = field_des%mask
  call mpp_global_field (geom_src%Domain%mpp_domain, tmp(:,:,1:nz_), gdata(:,:,1:nz_) )

  ! Interpolate to destination geometry
  call soca_hinterp(geom_des,field_des%val,gdata,mask_(:,:),nz_,missing,lon_in,lat_in,field_des%lon,field_des%lat)

  ! Update halos
  call mpp_update_domains(field_des%val, geom_des%Domain%mpp_domain)

end subroutine soca_convertstate_change_resol2d

! ------------------------------------------------------------------------------
subroutine soca_convertstate_change_resol(self, field_src, field_des, geom_src, geom_des)
  class(soca_convertstate_type),  intent(inout) :: self
  type(soca_field), pointer,      intent(inout) :: field_src, field_des
  type(soca_geom),                intent(inout) :: geom_src, geom_des

  !local
  integer :: i, j, k, tmp_nz, nz_
  integer :: isc1, iec1, jsc1, jec1, isd1, ied1, jsd1, jed1, isg, ieg, jsg, jeg
  integer :: isc2, iec2, jsc2, jec2, isd2, ied2, jsd2, jed2
  type(remapping_CS)  :: remapCS2
  type(horiz_interp_type) :: Interp
  real(kind=kind_real) :: missing = 0.d0
  real(kind=kind_real) :: PI_180, z_tot
  real(kind=kind_real), dimension(geom_src%isg:geom_src%ieg) :: lon_in
  real(kind=kind_real), dimension(geom_src%jsg:geom_src%jeg) :: lat_in
  real(kind=kind_real), dimension(geom_des%isd:geom_des%ied,geom_des%jsd:geom_des%jed) :: lon_out, lat_out
  real(kind=kind_real), dimension(geom_des%isd:geom_des%ied,geom_des%jsd:geom_des%jed) :: mask_
  real(kind=kind_real), allocatable :: tmp(:,:,:), tmp2(:,:,:), gdata(:,:,:)
  real(kind=kind_real), allocatable :: h1(:), h2(:)
  real(kind=kind_real), dimension(geom_src%isd:geom_src%ied,geom_src%jsd:geom_src%jed,1:geom_src%nzo_zstar) :: h_new1
  real(kind=kind_real), dimension(geom_des%isd:geom_des%ied,geom_des%jsd:geom_des%jed,1:geom_des%nzo_zstar) :: h_new2

  PI_180=atan(1.0d0)/45.0d0

  ! Indices for compute, data, and global domain for source
  isc1 = geom_src%isc ; iec1 = geom_src%iec ; jsc1 = geom_src%jsc ; jec1 = geom_src%jec
  isd1 = geom_src%isd ; ied1 = geom_src%ied ; jsd1 = geom_src%jsd ; jed1 = geom_src%jed
  isg = geom_src%isg ; ieg = geom_src%ieg ; jsg = geom_src%jsg ; jeg = geom_src%jeg

  ! Indices for compute and data domain for des
  isc2 = geom_des%isc ; iec2 = geom_des%iec ; jsc2 = geom_des%jsc ; jec2 = geom_des%jec
  isd2 = geom_des%isd ; ied2 = geom_des%ied ; jsd2 = geom_des%jsd ; jed2 = geom_des%jed

  ! Initialize vertical remapping
  call initialize_remapping(remapCS2,'PPM_IH4')

  ! Set grid thickness based on zstar level for src & target grid
  if (field_des%metadata%io_file=="ocn".or.field_des%metadata%io_file=='ice') then
    mask_ = field_des%mask
    h_new1(isc1:iec1,jsc1:jec1,1:geom_src%nzo_zstar) = geom_src%h_zstar(isc1:iec1,jsc1:jec1,1:geom_src%nzo_zstar)
    h_new2(isc2:iec2,jsc2:jec2,1:geom_des%nzo_zstar) = geom_des%h_zstar(isc2:iec2,jsc2:jec2,1:geom_des%nzo_zstar)
    call mpp_update_domains(mask_, geom_des%Domain%mpp_domain)
    call mpp_update_domains(h_new1, geom_src%Domain%mpp_domain)
    call mpp_update_domains(h_new2, geom_des%Domain%mpp_domain)
  else
    mask_ = 1.d0
  end if

  ! target hocn has been set in setup
  if (field_des%name == "hocn" ) then
    return
  end if
  lon_in = geom_src%lonh ; lat_in = geom_src%lath
  if (field_src%name == "uocn" .and. field_des%name == "uocn") lon_in = geom_src%lonq
  if (field_src%name == "vocn" .and. field_des%name == "vocn") lat_in = geom_src%latq
!  call meshgrid(geom_des%lonh(isd2:ied2),geom_des%lath(jsd2:jed2),lon_out,lat_out)
!  if (field_des%name == "uocn") call meshgrid(geom_des%lonq(isd2:ied2),geom_des%lath(jsd2:jed2),lon_out,lat_out)
!  if (field_des%name == "vocn") call meshgrid(geom_des%lonh(isd2:ied2),geom_des%latq(jsd2:jed2),lon_out,lat_out)

  ! Converts src grid to zstar coordinate
  nz_ = geom_src%nzo_zstar
  if (field_src%nz == 1 .or. field_src%metadata%io_file=="ice") nz_ = field_src%nz
  allocate(tmp(isd1:ied1,jsd1:jed1,1:nz_),gdata(isg:ieg,jsg:jeg,1:nz_),tmp2(isd2:ied2,jsd2:jed2,1:nz_))
  allocate(h1(field_src%nz),h2(nz_))
  tmp = 0.d0 ; gdata = 0.d0 ; tmp2 = 0.d0;
  if ( field_src%nz > 1 .and. field_src%metadata%io_file/="ice") then
    do j = jsc1, jec1
      do i = isc1, iec1
        tmp_nz = field_src%nz
        if(field_src%name =="uocn") then
          if (field_src%mask(i,j)>0.) then
            h1(1:tmp_nz) = 0.5 * ( self%hocn_src(i,j,1:tmp_nz) + self%hocn_src(i+1,j,1:tmp_nz) )
            h2(1:nz_) = 0.5 * ( h_new1(i,j,1:nz_) + h_new1(i+1,j,1:nz_) )
            call remapping_core_h(remapCS2, tmp_nz, h1(1:tmp_nz), field_src%val(i,j,1:tmp_nz), &
                                  nz_, h2(1:nz_), tmp(i,j,1:nz_))
          endif
        else if(field_src%name =="vocn") then
          if (field_src%mask(i,j)>0.) then
            h1(1:tmp_nz) = 0.5 * ( self%hocn_src(i,j,1:tmp_nz) + self%hocn_src(i,j+1,1:tmp_nz) )
            h2(1:nz_) = 0.5 * ( h_new1(i,j,1:nz_) + h_new1(i,j+1,1:nz_) )
            call remapping_core_h(remapCS2, tmp_nz, h1(1:tmp_nz), field_src%val(i,j,1:tmp_nz), &
                                  nz_, h2(1:nz_), tmp(i,j,1:nz_))
          endif
        else
          if (field_src%mask(i,j) > 0.d0) then
            call remapping_core_h(remapCS2, tmp_nz, self%hocn_src(i,j,1:tmp_nz), field_src%val(i,j,1:tmp_nz), &
                                  nz_, h_new1(i,j,1:nz_), tmp(i,j,1:nz_))
          endif
        end if
      end do !i
    end do !j
  else
    if (field_src%metadata%io_file=="ocn") tmp(:,:,1) = field_src%val(:,:,1) !*field_src%mask(:,:) !2D
    if (field_src%metadata%io_file=="sfc") tmp(:,:,1) = field_src%val(:,:,1) !2D no mask
    if (field_src%metadata%io_file=="ice") tmp(:,:,1:nz_) = field_src%val(:,:,1:nz_)
  end if ! field_src%nz > 1
  call mpp_update_domains(tmp, geom_src%Domain%mpp_domain)

  ! Convert src field to target field at zstar coord
  call mpp_global_field (geom_src%Domain%mpp_domain, tmp(:,:,1:nz_), gdata(:,:,1:nz_) )
  call soca_hinterp(geom_des,tmp2(:,:,1:nz_),gdata,mask_(:,:),nz_,missing,lon_in,lat_in,field_des%lon,field_des%lat)

  call mpp_update_domains(tmp2, geom_des%Domain%mpp_domain)

  ! Final step: vertical remapping to desired vertical coordinate
  if(allocated(h1)) deallocate(h1)
  if(allocated(h2)) deallocate(h2)
  allocate(h1(nz_),h2(field_des%nz))
  if ( field_des%nz > 1 .and. field_des%metadata%io_file/="ice") then
    do j = jsc2, jec2
      do i = isc2, iec2
        tmp_nz = nz_ !assume geom_src%nzo_zstar == geom%des%nzo_zstar
        if(field_des%name =="uocn") then
          if (field_des%mask(i,j)>0.) then
            h1(1:tmp_nz) = 0.5 * ( h_new2(i,j,1:tmp_nz) + h_new2(i+1,j,1:tmp_nz) )
            h2(1:field_des%nz) = 0.5 * ( self%hocn_des(i,j,1:field_des%nz) + self%hocn_des(i+1,j,1:field_des%nz) )
            call remapping_core_h(remapCS2, tmp_nz, h1(1:tmp_nz), tmp2(i,j,1:tmp_nz), &
                                  field_des%nz, h2(1:field_des%nz), field_des%val(i,j,1:field_des%nz))
          end if
        else if (field_des%name =="vocn") then
           if (field_des%mask(i,j)>0.) then
             h1(1:tmp_nz) = 0.5 * ( h_new2(i,j,1:tmp_nz) + h_new2(i,j+1,1:tmp_nz) )
             h2(1:field_des%nz) = 0.5 * ( self%hocn_des(i,j,1:field_des%nz) + self%hocn_des(i,j+1,1:field_des%nz) )
             call remapping_core_h(remapCS2, tmp_nz, h1(1:tmp_nz), tmp2(i,j,1:tmp_nz), &
                                   field_des%nz, h2(1:field_des%nz), field_des%val(i,j,1:field_des%nz))
           end if
        else
          if (field_des%mask(i,j)>0.) then
            call remapping_core_h(remapCS2, tmp_nz, h_new2(i,j,1:tmp_nz), tmp2(i,j,1:tmp_nz), &
                                  field_des%nz, self%hocn_des(i,j,1:field_des%nz), field_des%val(i,j,1:field_des%nz))
          end if
        end if
      end do !j
    end do !i
  else
   if (field_des%metadata%io_file=="ocn") field_des%val(:,:,1) = tmp2(:,:,1)*field_des%mask(:,:) ! 2D
   if (field_des%metadata%io_file=="sfc") field_des%val(:,:,1) = tmp2(:,:,1) ! 2D no mask
   if (field_des%metadata%io_file=="ice") field_des%val(:,:,1:field_des%nz) = tmp2(:,:,1:field_des%nz)
  end if ! nz > 1

  call mpp_update_domains(field_des%val, geom_des%Domain%mpp_domain)

end subroutine soca_convertstate_change_resol


! ------------------------------------------------------------------------------
subroutine soca_hinterp(self,field2,gdata,mask2,nz,missing,lon_in,lat_in,lon_out,lat_out)
  class(soca_geom),  intent(inout) :: self
  real(kind=kind_real), dimension(self%isd:self%ied,self%jsd:self%jed,1:nz), intent(inout) :: field2
  real(kind=kind_real), dimension(:,:,:), intent(in) :: gdata
  real(kind=kind_real), dimension(self%isd:self%ied,self%jsd:self%jed), intent(in) :: mask2
  integer, intent(in) :: nz
  real(kind=kind_real), intent(in) :: missing
  real(kind=kind_real), dimension(:), intent(in) :: lon_in, lat_in
  real(kind=kind_real), dimension(self%isd:self%ied,self%jsd:self%jed), intent(in) :: lon_out, lat_out

  !local variables
  integer :: i, j, k, isg, ieg, jsg, jeg, jeg1
  integer :: isc2, iec2, jsc2, jec2, npoints
  real(kind=kind_real) :: roundoff = 1.e-5
  real(kind=kind_real) :: PI_180
  type(horiz_interp_type) :: Interp
  type(ocean_grid_type) :: grid
  real(kind_real), dimension(:), allocatable :: lath_inp
  real(kind_real), dimension(:,:), allocatable :: lon_inp, lat_inp, tr_inp, mask_in_
  real(kind_real), dimension(self%isd:self%ied,self%jsd:self%jed) :: tr_out, fill, good, prev, mask_out_
  real(kind=kind_real) :: max_lat,min_lat, pole, npole, varavg
  real(kind=kind_real), dimension(:), allocatable :: last_row, lonh, lath
  logical :: add_np, add_sp

  PI_180=atan(1.0d0)/45.0d0

  isg = 1; jsg = 1;
  ieg = size(gdata,1); jeg = size(gdata,2)

  ! Indices for compute domain for regional model
  isc2 = self%isc ; iec2 = self%iec ; jsc2 = self%jsc ; jec2 = self%jec

  grid%isc = self%isc ; grid%iec = self%iec ; grid%jsc = self%jsc ; grid%jec = self%jec
  grid%isd = self%isd ; grid%ied = self%ied ; grid%jsd = self%jsd ; grid%jed = self%jed
  grid%Domain => self%Domain

  jeg1=jeg
  max_lat = maxval(lat_in)
  add_np=.false.
  if (max_lat < 90.0) then
    add_np=.true.
    jeg1=jeg1+1
    allocate(lath(jsg:jeg1))
    lath(jsg:jeg)=lat_in(:)
    lath(jeg1)=90.d0
  else
    allocate(lath(jsg:jeg1))
    lath(:) = lat_in
  endif
  min_lat = minval(lat_in)
  add_sp=.false.
  if (min_lat > -90.0) then
    add_sp=.true.
    jeg1=jeg1+1
    if (allocated(lath_inp)) deallocate(lath_inp)
    allocate(lath_inp(jeg1))
    lath_inp(jsg+1:jeg1)=lath(:)
    lath_inp(jsg)=-90.d0
    if (allocated (lath)) deallocate(lath)
    allocate(lath(jsg:jeg1))
    lath(:)=lath_inp(:)
  endif

  allocate(lonh(isg:ieg))
  lonh(:) = lon_in(:)

  allocate(lon_inp(isg:ieg,jsg:jeg1))
  allocate(lat_inp(isg:ieg,jsg:jeg1))
  call meshgrid(lonh,lath,lon_inp,lat_inp)

  allocate(mask_in_(isg:ieg,jsg:jeg1))
  allocate(tr_inp(isg:ieg,jsg:jeg1)) ; allocate(last_row(isg:ieg))
  do k = 1, nz
    ! extrapolate the input data to the north pole using the northerm-most latitude
    if (is_root_pe()) then
      if (add_np) then
        last_row(:)=gdata(:,jeg,k); pole=0.d0; npole=0.d0
        do i=isg,ieg
          if (abs(gdata(i,jeg,k)-missing) > abs(roundoff)) then
            pole = pole+last_row(i)
            npole = npole+1.d0
          endif
        enddo
        if (npole > 0) then
          pole=pole/npole
        else
          pole=missing
        endif

        if (add_sp) then
          tr_inp(:,jsg) = gdata(:,jsg,k)
          tr_inp(:,jsg+1:jeg1-1) = gdata(:,:,k)
          tr_inp(:,jeg1) = pole
        else
          tr_inp(:,jsg:jeg) = gdata(:,:,k)
          tr_inp(:,jeg1) = pole
        endif !add_sp

      else
        if (add_sp) then
          tr_inp(:,jsg) = gdata(:,jsg,k)
          tr_inp(:,jsg+1:jeg1) = gdata(:,:,k)
        else
          tr_inp(isg:ieg,jsg:jeg) = gdata(isg:ieg,jsg:jeg,k)
        endif !add_sp
      endif !add_np
    end if !root_pe

    call mpp_sync()
    call mpp_broadcast(tr_inp, ieg*jeg1, root_PE())
    call mpp_sync_self()

    mask_in_ = 1.d0
    do j=jsg,jeg1 ; do i=isg,ieg
      if (abs(tr_inp(i,j)-missing) <= abs(roundoff)) then
        tr_inp(i,j) = missing
        mask_in_(i,j) = 0.d0;
      end if
    enddo ; enddo

    tr_out(:,:) = 0.d0
    ! initialize horizontal remapping
    if (k==1) call horiz_interp_new(Interp, lon_inp(:,:)*PI_180, lat_inp(:,:)*PI_180, lon_out(isc2:iec2,jsc2:jec2)*PI_180, &
       lat_out(isc2:iec2,jsc2:jec2)*PI_180, interp_method='bilinear', src_modulo=.true., mask_in=mask_in_)

    call horiz_interp(Interp, tr_inp, tr_out(isc2:iec2,jsc2:jec2), mask_in=mask_in_, missing_value=missing, missing_permit=3)

    mask_out_ = 1.d0 ; fill = 0.d0 ; good = 0.d0
    npoints = 0 ; varavg = 0.d0
    do j=jsc2,jec2
      do i=isc2,iec2
        if (abs(tr_out(i,j)-missing) < abs(roundoff)) mask_out_(i,j)=0.d0
        if (mask_out_(i,j) < 1.0d0) then
          tr_out(i,j) = missing
        else
          good(i,j) = 1.0d0
          npoints = npoints + 1
          varavg = varavg + tr_out(i,j)
        endif
        if (mask2(i,j) == 1.d0 .and. mask_out_(i,j) < 1.0d0) fill(i,j) = 1.d0
      end do !i
    end do !j
    call pass_var(fill, self%Domain) ; call pass_var(good, self%Domain)
    call sum_across_pes(npoints) ; call sum_across_pes(varavg)
    if (npoints > 0) then
      varavg = varavg/real(npoints)
    end if

    if (k==1) prev(:,:) = tr_out(:,:)
    call fill_miss_2d(tr_out, good, fill, prev=prev, G=grid, smooth=.true.)

    !TODO: In case fill_miss_2d failed at surface (k=1), use IDW to fill data pt that is located in ocean mask
    !Problem: IDW is compiler-dependent

    field2(:,:,k) = tr_out(:,:)*mask2(:,:)
    prev(:,:) = field2(:,:,k)

  end do

end subroutine soca_hinterp

end module soca_convert_state_mod
