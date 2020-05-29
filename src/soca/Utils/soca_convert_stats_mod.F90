! (C) Copyright 2020-2020 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! ------------------------------------------------------------------------------

module soca_convert_statse_mod

  use soca_geom_mod
  use soca_fields_mod
  use kinds, only: kind_real
  use MOM_coms, only : max_across_PEs ! will remove in the future
  use fms_io_mod, only: read_data, write_data, fms_io_init, fms_io_exit 
  use MOM_remapping, only : remapping_CS, initialize_remapping, remapping_core_h
  use MOM_domains, only : pass_var, sum_across_PEs, root_PE
  use mpp_mod, only     : mpp_broadcast, mpp_sync, mpp_sync_self
  use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
  use mpp_domains_mod, only  : mpp_global_field, mpp_update_domains
  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_type
!  use MOM_horizontal_regridding, only : meshgrid, fill_miss_2d 
  use MOM_grid, only : ocean_grid_type 
!
  implicit none
  private

  type, public :: soca_convertstate_type
     real(kind=kind_real), allocatable, dimension(:,:,:) :: hocn_src, hocn_des 

   contains
     procedure :: setup => soca_convertstate_setup
     procedure :: change_resol => soca_convertstate_change_resol
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

  !
  call fms_io_init()
  
  call read_data(trim(src%geom_grid_file), 'nzo_zstar', tmp(1), domain=src%Domain%mpp_domain)
  src%nzo_zstar = tmp(1)

  call read_data(trim(des%geom_grid_file), 'nzo_zstar', tmp(1), domain=des%Domain%mpp_domain)
  des%nzo_zstar = tmp(1)
  
  if (des%nzo_zstar /= src%nzo_zstar) & 
    call MOM_error(FATAL, "target nzo_zstar /= source nzo_zstar! Reset maximum depth in target grid MOM_input file and re-run soca gridgen")

  if (allocated(src%h_zstar)) deallocate(src%h_zstar)
  allocate(src%h_zstar(src%isd:src%ied,src%jsd:src%jed,1:src%nzo_zstar))
  call read_data(trim(src%geom_grid_file), 'h_zstar', src%h_zstar, domain=src%Domain%mpp_domain)

  if (allocated(des%h_zstar)) deallocate(des%h_zstar)
  allocate(des%h_zstar(des%isd:des%ied,des%jsd:des%jed,1:des%nzo_zstar))
  call read_data(trim(des%geom_grid_file), 'h_zstar', des%h_zstar, domain=des%Domain%mpp_domain)

  call fms_io_exit()  

  !
  allocate(self%hocn_src(src%isd:src%ied,src%jsd:src%jed,1:src%nzo))
  allocate(self%hocn_des(des%isd:des%ied,des%jsd:des%jed,1:des%nzo))

  ! set hocn for target grid
  hocn2%val = des%h
  !
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
subroutine soca_convertstate_change_resol(self, field_src, field_des, geom_src, geom_des)
  class(soca_convertstate_type),  intent(inout) :: self
  type(soca_field), pointer,      intent(inout) :: field_src, field_des 
  type(soca_geom),                intent(inout) :: geom_src, geom_des 

  !local
  integer :: i, j, k, n, tmp_nz, nz_
  integer :: isc1, iec1, jsc1, jec1, isd1, ied1, jsd1, jed1, isg, ieg, jsg, jeg 
  integer :: isc2, iec2, jsc2, jec2, isd2, ied2, jsd2, jed2 
  type(remapping_CS)  :: remapCS2
  type(horiz_interp_type) :: Interp
  real(kind=kind_real) :: missing = 0.d0
  real(kind=kind_real) :: PI_180, z_tot 
  real(kind=kind_real), dimension(geom_src%isg:geom_src%ieg) :: lon_in
  real(kind=kind_real), dimension(geom_src%jsg:geom_src%jeg) :: lat_in 
  real(kind=kind_real), dimension(geom_des%isd:geom_des%ied,geom_des%jsd:geom_des%jed) :: mask_
  real(kind=kind_real), allocatable :: tmp(:,:,:), tmp2(:,:,:), gdata(:,:,:)
  real(kind=kind_real), allocatable :: h1(:), h2(:)
  real(kind=kind_real), dimension(geom_src%isd:geom_src%ied,geom_src%jsd:geom_src%jed,1:geom_src%nzo_zstar) :: h_new1
  real(kind=kind_real), dimension(geom_des%isd:geom_des%ied,geom_des%jsd:geom_des%jed,1:geom_des%nzo_zstar) :: h_new2

  !
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
  if (field_des%io_file=="ocn") then
    mask_ = field_des%mask 
    h_new1 = geom_src%h_zstar 
    h_new2 = geom_des%h_zstar  
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

  ! Converts src grid to zstar coordinate 
  nz_ = geom_src%nzo_zstar 
  if (field_src%nz == 1 .or. field_src%io_file=="ice") nz_ = field_src%nz
  allocate(tmp(isd1:ied1,jsd1:jed1,1:nz_),gdata(isg:ieg,jsg:jeg,1:nz_),tmp2(isd2:ied2,jsd2:jed2,1:nz_))
  allocate(h1(field_src%nz),h2(nz_))
  tmp = 0.d0 ; gdata = 0.d0 ; tmp2 = 0.d0;
  if ( field_src%nz > 1 .and. field_src%io_file/="ice") then
    do i = isc1-1, iec1+1
      do j = jsc1-1, jec1+1
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
      end do !j
    end do !i
  else
    if (field_src%io_file=="ocn") tmp(:,:,1) = field_src%val(:,:,1)*field_src%mask(:,:) !2D 
    if (field_src%io_file=="sfc") tmp(:,:,1) = field_src%val(:,:,1) !2D no mask
    if (field_src%io_file=="ice") tmp(:,:,1:nz_) = field_src%val(:,:,1:nz_)
  end if ! field_src%nz > 1
  call mpp_update_domains(tmp, geom_src%Domain%mpp_domain)

  ! horizontal interp: convert src field to target field at zstar coord 
  call mpp_global_field (geom_src%Domain%mpp_domain, tmp(:,:,1:nz_), gdata(:,:,1:nz_) )
  call soca_hinterp(geom_des,tmp2(:,:,1:nz_),gdata,mask_(:,:),nz_,missing,lon_in,lat_in,field_des%lon,field_des%lat) 

  call mpp_update_domains(tmp2, geom_des%Domain%mpp_domain)

  ! Final step: vertical remapping to desired vertical coordinate
  if(allocated(h1)) deallocate(h1)
  if(allocated(h2)) deallocate(h2) 
  allocate(h1(nz_),h2(field_des%nz))
  if ( field_des%nz > 1 .and. field_des%io_file/="ice") then
    do i = isc2-1, iec2+1
      do j = jsc2-1, jec2+1
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
      end do !i
    end do !j
  else
   if (field_des%io_file=="ocn") field_des%val(:,:,1) = tmp2(:,:,1)*field_des%mask(:,:) ! 2D
   if (field_des%io_file=="sfc") field_des%val(:,:,1) = tmp2(:,:,1) ! 2D no mask
   if (field_des%io_file=="ice") field_des%val(:,:,1:field_des%nz) = tmp2(:,:,1:field_des%nz)
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
  integer :: isc2, iec2, jsc2, jec2
  integer :: npoints
  real(kind=kind_real) :: roundoff = 1.e-9
  real(kind=kind_real) :: PI_180
  type(horiz_interp_type) :: Interp
  type(ocean_grid_type) :: grid 
  real(kind_real), dimension(:), allocatable :: lath_inp 
  real(kind_real), dimension(:,:), allocatable :: lon_inp, lat_inp, tr_inp
  real(kind_real), dimension(self%isd:self%ied,self%jsd:self%jed) :: tr_out, fill, good, prev, mask_out_ 
  real(kind=kind_real) :: max_lat,min_lat,pole,npole, varavg
  real(kind=kind_real), dimension(:), allocatable :: last_row, lonh, lath 
  logical :: add_np, add_sp 
  !
  PI_180=atan(1.0d0)/45.0d0 
!  roundoff=3.0*epsilon(missing)
  !
  isg = 1; jsg = 1;
  ieg = size(gdata,1); jeg = size(gdata,2)

  ! Indices for compute domain for regional model
  isc2 = self%isc ; iec2 = self%iec ; jsc2 = self%jsc ; jec2 = self%jec

  !
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
  call soca_meshgrid(lonh,lath,lon_inp,lat_inp)
  
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

    do j=jsg,jeg1 ; do i=isg,ieg
      if (abs(tr_inp(i,j)-missing) <= abs(roundoff)) tr_inp(i,j) = missing 
    enddo ; enddo

    tr_out(:,:) = 0.d0
    ! initialize horizontal remapping 
    if (k==1) call horiz_interp_new(Interp, lon_inp(:,:)*PI_180, lat_inp(:,:)*PI_180, lon_out(isc2:iec2,jsc2:jec2)*PI_180, &
       lat_out(isc2:iec2,jsc2:jec2)*PI_180, interp_method='bilinear', src_modulo=.true.)

    call horiz_interp(Interp,tr_inp,tr_out(isc2:iec2,jsc2:jec2),&
         missing_value=missing, new_missing_handle=.true.)

    mask_out_ = 1.d0 ; fill = 0.d0 ; good = 0.d0
    do j=jsc2,jec2
      do i=isc2,iec2
        if (abs(tr_out(i,j)-missing) < abs(roundoff)) mask_out_(i,j)=0.d0
      end do
    end do 
    npoints = 0 ; varavg = 0.d0 
    do j=jsc2,jec2
      do i=isc2,iec2
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
    call soca_fill_miss_2d(tr_out, good, fill, prev=prev, G=grid, smooth=.true.)

    field2(:,:,k) = tr_out(:,:)*mask2(:,:) 
    prev(:,:) = field2(:,:,k)

  end do

end subroutine soca_hinterp

! ------------------------------------------------------------------------------
!TODO: Remove the below 2 subroutines: meshgrid, fill_miss_2d once they become public in MOM6
! ------------------------------------------------------------------------------
!> Create a 2d-mesh of grid coordinates from 1-d arrays.
subroutine soca_meshgrid(x, y, x_T, y_T)
  real(kind=kind_real), dimension(:),                   intent(in)    :: x  !< input 1-dimensional vector
  real(kind=kind_real), dimension(:),                   intent(in)    :: y  !< input 1-dimensional vector
  real(kind=kind_real), dimension(size(x,1),size(y,1)), intent(inout) :: x_T !< output 2-dimensional array
  real(kind=kind_real), dimension(size(x,1),size(y,1)), intent(inout) :: y_T !< output 2-dimensional array

  integer :: ni,nj,i,j

  ni=size(x,1) ; nj=size(y,1)

  do j=1,nj ; do i=1,ni
    x_T(i,j) = x(i)
  enddo ; enddo

  do j=1,nj ; do i=1,ni
    y_T(i,j) = y(j)
  enddo ; enddo

end subroutine soca_meshgrid

! ------------------------------------------------------------------------------
subroutine soca_fill_miss_2d(aout, good, fill, prev, G, smooth, num_pass, relc, crit )

  type(ocean_grid_type), intent(inout) :: G
  real(kind=kind_real), dimension(G%isd:G%ied,G%jsd:G%jed), intent(inout) :: aout
  real(kind=kind_real), dimension(G%isd:G%ied,G%jsd:G%jed), intent(in) :: good,fill
  real(kind=kind_real), dimension(G%isd:G%ied,G%jsd:G%jed), optional, intent(in) :: prev
  logical,     optional, intent(in)    :: smooth ! If present and true, apply a number of
                                                 ! Laplan iterations to the interpolated data
  integer,     optional, intent(in)    :: num_pass ! The maximum number of iterations
  real(kind=kind_real),        optional, intent(in)    :: relc ! A relaxation coefficient for Laplacian (ND)
  real(kind=kind_real),        optional, intent(in)    :: crit ! A minimal value for deltas between iterations.

  real(kind=kind_real), dimension(G%isd:G%ied,G%jsd:G%jed) :: b,r
  real(kind=kind_real), dimension(G%isd:G%ied,G%jsd:G%jed) :: fill_pts, good_, good_new

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

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  npass = num_pass_default
  if (PRESENT(num_pass)) npass = num_pass

  relax_coeff = relc_default
  if (PRESENT(relc)) relax_coeff = relc

  acrit = crit_default
  if (PRESENT(crit)) acrit = crit

  do_smooth=.false.
  if (PRESENT(smooth)) do_smooth=smooth

  fill_pts(:,:) = fill(:,:)

  nfill = sum(fill(is:ie,js:je))
  call sum_across_PEs(nfill)

  nfill_prev = nfill
  good_(:,:) = good(:,:)
  r(:,:) = 0.0

  do while (nfill > 0.0)

    call pass_var(good_, G%Domain)
    call pass_var(aout, G%Domain)

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

      ngood = ge+gw+gn+gs

      if (ngood > 0.) then
        b(i,j)=(east+west+north+south)/ngood
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
    call pass_var(aout,G%Domain)
    do j=js,je ; do i=is,ie
      if (fill(i,j) == 1) then
        east = max(good(i+1,j),fill(i+1,j)) ; west = max(good(i-1,j),fill(i-1,j))
        north = max(good(i,j+1),fill(i,j+1)) ; south = max(good(i,j-1),fill(i,j-1))
        r(i,j) = relax_coeff*(south*aout(i,j-1)+north*aout(i,j+1) + &
                              west*aout(i-1,j)+east*aout(i+1,j) - &
                             (south+north+west+east)*aout(i,j))
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

end subroutine soca_fill_miss_2d

end module soca_convert_statse_mod  
