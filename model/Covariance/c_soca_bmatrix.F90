! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! ------------------------------------------------------------------------------

!> Setup for the SOCA model's background error covariance matrix

subroutine c_soca_b_setup(c_key_self, c_conf, c_key_geom) &
     & bind (c,name='soca_b_setup_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_geom_mod
  use tools_const, only: deg2rad
  use hdiag_nicas_mod
  use type_nam, only: namtype,namcheck
  use model_oops, only: model_oops_coord
  use obsop_parameters, only: compute_parameters
  use type_nam, only: namtype, namncwrite
  use type_geom, only: geomtype,compute_grid_mesh
  use type_odata, only: odatatype
  use type_randgen, only: create_randgen  
  use type_bpar
  use type_bdata
  use driver_nicas, only: run_nicas
  use nicas_apply_nicas
  use nicas_parameters, only: nicas_compute_parameters => compute_parameters
  use driver_hdiag, only: run_hdiag
  use netcdf
  use tools_nc, only: ncfloat,ncerr
  use nicas_parameters_convol
  use type_mpl, only: mpl,mpl_start
  
  implicit none

  integer(c_int), intent(inout) :: c_key_self   !< The background covariance structure
  type(c_ptr), intent(in)    :: c_conf          !< The configuration
  integer(c_int), intent(in) :: c_key_geom      !< Geometry
  type(soca_3d_covar_config), pointer :: self
  type(soca_geom),  pointer :: geom

  !nicas stuff
  integer :: nc0a, nl0, nv, nts
  real(kind=kind_real), allocatable :: lon(:), lat(:), area(:), vunit(:), rndnum(:)
  integer, allocatable :: imask(:,:)
  integer :: ens1_ne = 4
  type(namtype) :: nam
  integer :: ncid

  !Grid stuff
  integer :: isc, iec, jsc, jec, jjj, jz, il, ib
  character(len=1024) :: subr = 'model_write'

  call mpl_start
  
  call soca_geom_registry%get(c_key_geom, geom)
  call soca_3d_cov_registry%init()
  call soca_3d_cov_registry%add(c_key_self)
  call soca_3d_cov_registry%get(c_key_self, self)

  !--- Initialize geometry to be passed to NICAS
  ! Indices for compute domain (no halo)
  isc = geom%ocean%G%isc
  iec = geom%ocean%G%iec
  jsc = geom%ocean%G%jsc
  jec = geom%ocean%G%jec

  nv = geom%ocean%ncat + 1                   !< Number of variables
  nl0 = 1                                    !< Number of independent levels
  nts = 1                                    !< Number of time slots
  nc0a = (iec - isc + 1) * (jec - jsc + 1 )  !< Total number of grid cells in the compute domain

  allocate( lon(nc0a), lat(nc0a), area(nc0a) )
  allocate( vunit(nl0) )
  allocate( imask(nc0a, nl0) )    
  lon = deg2rad*reshape( geom%ocean%lon(isc:iec, jsc:jec), (/nc0a/) )
  lat = deg2rad*reshape( geom%ocean%lat(isc:iec, jsc:jec), (/nc0a/) ) 
  area = reshape( geom%ocean%cell_area(isc:iec, jsc:jec), (/nc0a/) )
  do jz = 1, nl0       
     vunit(jz) = real(jz)
     imask(1:nc0a,jz) = reshape( geom%ocean%mask2d(isc:iec, jsc:jec), (/nc0a/) )
  end do
  vunit = 1.0                      !< Dummy vertical unit
  self%nicasB%geom%nc0a = nc0a     !< Number of grid points (local)
  self%nicasB%geom%nl0 = 1         !< Number of levels: only one level here (same interpolation for all levels)
  self%nicasB%geom%nlev = self%nicasB%geom%nl0 ! Copy
  call model_oops_coord(self%nicasB%geom,lon,lat,area,vunit,imask)
  
  !--- Initialize namelist (parse yml file/hardcode)
  call hdiag_nicas_read_conf(c_conf,self%nicasB%nam)
  self%nicasB%nam%model='oops'  
  self%nicasB%nam%datadir='./Data'
  self%nicasB%nam%nl = nl0
  self%nicasB%nam%levs = nl0
  self%nicasB%nam%nv = nv
  self%nicasB%nam%ens1_ne_offset = 0  
  self%nicasB%nam%ens1_nsub = 1
  self%nicasB%nam%ens1_ne = ens1_ne
  do il=1,self%nicasB%nam%nl
     self%nicasB%nam%levs(il) = il
  end do
  self%nicasB%nam%nts = nts

  !--- Check parameters
  call namcheck(self%nicasB%nam)   
  
  !--- Initialize random number generator
  call create_randgen(nam)

  !--- Write parameters
  call ncerr(subr,nf90_create('nam.nc',or(nf90_clobber,nf90_64bit_offset),ncid))
  call namncwrite(self%nicasB%nam,ncid)
  call ncerr(subr,nf90_close(ncid))
  
  !--- Setup nicas geometry
  call compute_grid_mesh(self%nicasB%nam,self%nicasB%geom)  

  deallocate(lon,lat,area)
  deallocate(vunit)
  deallocate(imask)

  !--- Initialize block parameters: bpar
  call bpar_alloc(self%nicasB%nam,self%nicasB%geom,self%nicasB%bpar)

  !--- Initialize bdata
  self%nicasB%bpar%diag_block = .true.  ! Needed to allocate bdata
  call bdata_alloc(self%nicasB%nam,self%nicasB%geom,self%nicasB%bpar,self%nicasB%bdata)  
  do ib = 1,self%nicasB%bpar%nb + 1
     self%nicasB%bdata(ib)%rh0=self%nicasB%nam%dc
     self%nicasB%bdata(ib)%rv0=1000.0!self%nicasB%nam%dc
     self%nicasB%bdata(ib)%rh0s=self%nicasB%nam%dc
     self%nicasB%bdata(ib)%rv0s=1000.0!self%nicasB%nam%dc
  end do  
  call run_nicas(self%nicasB%nam,self%nicasB%geom,self%nicasB%bpar,self%nicasB%bdata,self%nicasB%ndata)
  
  print *,'nsb=',self%nicasB%ndata(1)%nsb
  print *,'convolution: ',self%nicasB%ndata(1)%c%S
  print *,'----------- OUT OF COVAR SETUP --------'
  !call abor1_ftn("c_soca_bmatrix: Done covar_sqrt_mult")
end subroutine c_soca_b_setup

! ------------------------------------------------------------------------------
!> Delete for the SOCA model's background error covariance matrix

subroutine c_soca_b_delete(c_key_self) bind (c,name='soca_b_delete_f90')

  use iso_c_binding
  use soca_covariance_mod

  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure
  type(soca_3d_covar_config), pointer :: self

  call soca_3d_cov_registry%get(c_key_self,self)
  call soca_3d_covar_delete(c_key_self)
  call soca_3d_cov_registry%remove(c_key_self)

end subroutine c_soca_b_delete

! ------------------------------------------------------------------------------

!> Multiply by inverse of covariance

subroutine c_soca_b_inv_mult(c_key_conf, c_key_in, c_key_out) bind(c,name='soca_b_invmult_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_fields
  use kinds

  implicit none
  integer(c_int), intent(in) :: c_key_conf  !< covar config structure
  integer(c_int), intent(in) :: c_key_in    !< Streamfunction: psi
  integer(c_int), intent(in) :: c_key_out   !< Streamfunction: psi
  type(soca_3d_covar_config), pointer :: conf
  type(soca_field), pointer :: xin
  type(soca_field), pointer :: xout

  call soca_3d_cov_registry%get(c_key_conf,conf)
  call soca_field_registry%get(c_key_in,xin)
  call soca_field_registry%get(c_key_out,xout)

  print *,"[[[[[[[[[[[[[[[[[[[[[[[[ IN B INV MULT: NOT IMPLEMENTED ]]]]]]]]]]]]]]]]]]]]]]]]"

  call ones(xout)
  call self_schur(xout, xin)

end subroutine c_soca_b_inv_mult

! ------------------------------------------------------------------------------

!> Multiply by covariance

subroutine c_soca_b_mult(c_key_conf, c_key_in, c_key_out) bind(c,name='soca_b_mult_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_fields
  use kinds
  use soca_Butils
  use nicas_apply_nicas
  use nicas_apply_localization
  use type_cv
  
  implicit none
  integer(c_int), intent(in) :: c_key_conf  !< covar config structure
  integer(c_int), intent(in) :: c_key_in    !< Streamfunction: psi
  integer(c_int), intent(in) :: c_key_out   !< Streamfunction: psi
  type(soca_3d_covar_config), pointer :: B
  type(soca_field), pointer :: xin
  type(soca_field), pointer :: xout
  integer :: ncat, k
  !real(kind=kind_real) :: Lx=5.0, Ly=1.0, sig_sic=0.01, sig_sit=0.5
  real(kind=kind_real) :: Lx=1.0, Ly=.5, sig_sic=0.05, sig_sit=150.0

  real(kind_real), allocatable :: fld(:,:,:)!(geom%nc0a,geom%nl0,nam%nv,nam%nts) !< Field
  
  !real(kind_real), allocatable :: fld(:,:)
  type(cvtype), allocatable :: cv(:) !control vector 
  real(kind_real), allocatable :: alpha(:)
  integer :: ll  
  call soca_3d_cov_registry%get(c_key_conf,B)
  call soca_field_registry%get(c_key_in,xin)
  call soca_field_registry%get(c_key_out,xout)

  call zeros(xout)
  print *,"[[[[[[[[[[[[[[[[[[[[[[[[ IN B MULT ]]]]]]]]]]]]]]]]]]]]]]]]"

  !allocate(cv(B%nicasB%bpar%nb+1))
  !call cv_alloc(B%nicasB%bpar,B%nicasB%ndata,cv)
  
  !allocate(fld(B%nicasB%geom%nc0a,B%nicasB%geom%nl0))
  allocate(fld(B%nicasB%geom%nc0a,B%nicasB%geom%nl0,B%nicasB%nam%nv))!,B%nicasB%nam%nts)) !< Field
  fld=0.0
  fld(1:100,1,1)=1.0
  print *,'shape ndata:',shape(B%nicasB%ndata)
  print *,'nc0a,nl0:',B%nicasB%geom%nc0a,B%nicasB%geom%nl0
  print *,'norm:',B%nicasB%ndata(1)%norm
  !call apply_nicas_from_sqrt(B%nicasB%geom,B%nicasB%ndata(1),fld)
  !call apply_localization(B%nicasB%nam,B%nicasB%geom,B%nicasB%bpar,B%nicasB%ndata,fld)
    
  !call apply_localization(B%nicasB%nam,B%nicasB%geom,B%nicasB%bpar,B%nicasB%ndata,fld)
  !allocate(alpha(B%nicasB%ndata(1)%nsa))
  
  !call apply_nicas_sqrt_ad(B%nicasB%geom,B%nicasB%ndata(1),fld,alpha)

  print *,'nsb=',B%nicasB%ndata(1)%nsb
  print *,'convolution: ',B%nicasB%ndata(1)%c%S

  call abor1_ftn("c_soca_bmatrix: Done covar_sqrt_mult")  
  do ll = 1,B%nicasB%geom%nc0a
     write(800,*) fld(ll,1,1)
  end do
  
  !call soca_3d_covar_sqrt_mult(xout, xin, conf)

  ncat = xin%geom%ocean%ncat

  !cicen, !hicen
  do k=2, ncat
     print *,'category:',k
     call gauss(xin%cicen(:,:,k), xout%cicen(:,:,k), xin%geom%ocean%lon, xin%geom%ocean%lat, lx, ly)
     call gauss(xin%hicen(:,:,k), xout%hicen(:,:,k), xin%geom%ocean%lon, xin%geom%ocean%lat, lx, ly)   
  end do
  xout%cicen = sig_sic**2 * xout%cicen
  xout%hicen = sig_sit**2 * xout%hicen

  !ssh
  print *,'ssh'
  call gauss(xin%ssh, xout%ssh, xin%geom%ocean%lon, xin%geom%ocean%lat, lx, ly)

end subroutine c_soca_b_mult

! ------------------------------------------------------------------------------

!> Generate randomized increment

subroutine c_soca_b_randomize(c_key_conf, c_key_out) bind(c,name='soca_b_randomize_f90')

  use iso_c_binding
  use soca_covariance_mod
  use soca_fields
  use random_vectors_mod
  use kinds

  implicit none
  integer(c_int), intent(in) :: c_key_conf  !< covar config structure
  integer(c_int), intent(in) :: c_key_out   !< Randomized increment
  type(soca_3d_covar_config), pointer :: conf
  type(soca_field), pointer :: xout
  real(kind=kind_real) :: prms ! Control vector
  !real(kind=kind_real), allocatable :: xctl(:,:,:) ! Control vector

  call soca_3d_cov_registry%get(c_key_conf,conf)
  call soca_field_registry%get(c_key_out,xout)

  !allocate(xctl(conf%nx, conf%ny, 2))

  !call random_vector(xctl(:,:,:))
  !call zeros(xout)
  call ones(xout)
  call random(xout)
  call fldrms(xout, prms)

  !call soca_3d_covar_sqrt_mult(conf%nx,conf%ny,xout,xctl,conf)

  !deallocate(xctl)

end subroutine c_soca_b_randomize

! ------------------------------------------------------------------------------
