! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> various utility functions
module soca_utils

use atlas_module, only: atlas_geometry, atlas_indexkdtree
use fckit_exception_module, only: fckit_exception
use gsw_mod_toolbox, only : gsw_rho, gsw_sa_from_sp, gsw_ct_from_pt, gsw_mlp
use kinds, only: kind_real
use netcdf

implicit none

private
public :: write2pe, soca_str2int, soca_adjust, &
          soca_rho, soca_diff, soca_mld, nc_check, &
          soca_stencil_interp, soca_stencil_neighbors

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> calculate density from temp/salinity profile
!!
!! \param sp: practical salinity profile
!! \param pt: potential temperature profile
!! \param p: pressure (depth) profile
!! \param lon: longitude
!! \param lat: latitude
elemental function soca_rho(sp, pt, p, lon, lat)
  real(kind=kind_real), intent(in)  :: pt, sp, p, lon, lat
  real(kind=kind_real) :: sa, ct, lon_rot, soca_rho

  !Rotate longitude if necessary
  lon_rot = lon
  if (lon<-180.0) lon_rot=lon+360.0
  if (lon>180.0) lon_rot=lon-360.0

  ! Convert practical salinity to absolute salinity
  sa = gsw_sa_from_sp (sp, p, lon_rot, lat)

  ! Convert potential temperature to concervative temperature
  ct = gsw_ct_from_pt (sa, pt)

  ! Insitu density
  soca_rho = gsw_rho(sa,ct,p)

  return
end function soca_rho


! ------------------------------------------------------------------------------
!> calculate mixed layer depth from temp/salinity profile
!!
!! \param pt: potential temperature profile
!! \param sp: practical salinity profile
!! \param p: pressure (depth) profile
!! \param lon: longitude
!! \param lat: latitude
function soca_mld(sp, pt, p, lon, lat)
  real(kind=kind_real), intent(in)  :: pt(:), sp(:), p(:), lon, lat

  real(kind=kind_real) :: lon_rot, soca_mld
  real(kind=kind_real), allocatable :: sa(:), ct(:)

  !Rotate longitude if necessary
  lon_rot = lon
  if (lon<-180.0) lon_rot=lon+360.0
  if (lon>180.0) lon_rot=lon-360.0

  ! Allocate memory
  allocate(sa(size(sp,1)),ct(size(sp,1)))

  ! Convert practical salinity to absolute salinity
  sa = gsw_sa_from_sp (sp, p, lon_rot, lat)

  ! Convert potential temperature to concervative temperature
  ct = gsw_ct_from_pt (sa, pt)

  ! Mixed layer depth
  soca_mld = gsw_mlp(sa,ct,p)
  if (soca_mld>9999.9_kind_real) soca_mld = p(1)
  soca_mld = max(soca_mld, p(1))

  deallocate(sa, ct)

  return
end function soca_mld


! ------------------------------------------------------------------------------
subroutine soca_diff(dvdz,v,h)

  real(kind=kind_real), intent(in)  :: v(:), h(:)
  real(kind=kind_real), intent(out) :: dvdz(:)

  integer :: k, ik

  k = size(v,1)

  do ik = 2, k-1
     dvdz(ik) = (v(ik+1)-v(ik-1))/(h(ik)+0.5*(h(ik+1)+h(ik-1)))
  end do
  dvdz(1) = dvdz(2)
  dvdz(k) = dvdz(k-1)

end subroutine soca_diff


! ------------------------------------------------------------------------------
!> per-PE file output, used by soca_geom to write per-PE grid
subroutine write2pe(vec,varname,filename,append)

  real(kind=kind_real), intent(in) :: vec(:)
  character(len=*),     intent(in) :: varname
  character(len=256),   intent(in) :: filename
  logical,              intent(in) :: append

  !netcdf stuff
  integer(kind=4) :: iNcid
  integer(kind=4) :: iDim_ID
  integer(kind=4) :: iVar_ID
  integer         :: ndims=1, ns

  ns=size(vec)
  if (append) then  ! If file exists, append to it
     call nc_check( nf90_open(filename, NF90_WRITE, iNcid) )
     call nc_check( nf90_inquire(iNcid, nDimensions = ndims) )
     call nc_check( nf90_inq_dimid(iNcid, "ns", iDim_ID) )
     call nc_check( nf90_redef(iNcid) )
  else
     call nc_check(nf90_create(filename, NF90_CLOBBER, iNcid))
     call nc_check(nf90_def_dim(iNcid, "ns", ns, iDim_ID))
  end if

  ! Define of variables.
  call nc_check( nf90_def_var(iNcid, trim(varname), NF90_DOUBLE,  (/ iDim_ID /), iVar_ID) )

  ! End define mode.
  call nc_check(nf90_enddef(iNcid))

  ! Writing
  call nc_check(nf90_put_var(iNcid, iVar_ID , vec))

  ! Close file.
  call nc_check(nf90_close(iNcid))

end subroutine write2pe


! ------------------------------------------------------------------------------
!> wrapper for netcdf calls
subroutine nc_check(status)

  integer(4), intent ( in) :: status
  if(status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop "Stopped"
  end if
end subroutine nc_check


! ------------------------------------------------------------------------------
!> Apply bounds
elemental function soca_adjust(std, minstd, maxstd)
  real(kind=kind_real), intent(in)  :: std, minstd, maxstd
  real(kind=kind_real) :: soca_adjust

  soca_adjust = min( max(std, minstd), maxstd)

end function soca_adjust


! ------------------------------------------------------------------------------
subroutine soca_str2int(str, int)
  character(len=*),intent(in) :: str
  integer,intent(out)         :: int

  read(str,*)  int
end subroutine soca_str2int


! ------------------------------------------------------------------------------
!> find the 6 stencil neighbors
subroutine soca_stencil_neighbors(fromto, i, j, ij)
  character(len=4), intent(in) :: fromto !< "u2h", "v2h"
  integer,          intent(in) :: i, j
  integer,          intent(out):: ij(2,6)

  integer :: sti

  select case(fromto)
  case("vtoh")
     ! return the 6 v-grid neighbors of h-grid(i,j)
     ij(1,1) = i;   ij(2,1) = j
     ij(1,2) = i;   ij(2,2) = j-1
     ij(1,3) = i+1; ij(2,3) = j
     ij(1,4) = i-1; ij(2,4) = j
     ij(1,5) = i-1; ij(2,5) = j-1
     ij(1,6) = i+1; ij(2,6) = j-1
  case("utoh")
     ! return the 6 u-grid neighbors of h-grid(i,j)
     ij(1,1) = i;   ij(2,1) = j
     ij(1,2) = i;   ij(2,2) = j-1
     ij(1,3) = i;   ij(2,3) = j+1
     ij(1,4) = i-1; ij(2,4) = j+1
     ij(1,5) = i-1; ij(2,5) = j
     ij(1,6) = i-1; ij(2,6) = j-1
  case default
     call fckit_exception%abort(fromto, "option is not implemented")
  end select

end subroutine soca_stencil_neighbors

! ------------------------------------------------------------------------------
!> idw interpolation given known neighbors
subroutine soca_stencil_interp(lon_src, lat_src, lon_dst, lat_dst, data, data_out, nn)
  real(kind_real),  intent(in) :: lon_src(:)
  real(kind_real),  intent(in) :: lat_src(:)
  real(kind_real),  intent(in) :: lon_dst
  real(kind_real),  intent(in) :: lat_dst
  real(kind_real), intent(in) :: data(:,:)
  real(kind_real), intent(inout) :: data_out(:)
  integer, intent(in) :: nn

  type(atlas_geometry) :: ageometry
  integer :: i, n, nz, k
  real(kind_real) :: val
  real(kind_real), allocatable :: w(:)

  ! Initialize atlas geometry on the sphere
  ageometry = atlas_geometry("UnitSphere")

  ! nn cannot be larger than 6
  if (nn > 6 ) call fckit_exception%abort( "Using more than 6 neighbors is not allowed")

  ! compute idw weights
  allocate(w(nn))
  do i = 1, nn
     w(i) = 1_kind_real/ageometry%distance(lon_src(i), lat_src(i), lon_dst, lat_dst)
  end do

  nz = size(data, dim=2)
  data_out = 0_kind_real

  do k = 1, nz
     val = 0_kind_real
     do i = 1, nn
        val = val + w(i)*data(i,k)
     end do
     data_out(k) = val/sum(w(1:nn))
  end do

end subroutine soca_stencil_interp

! ------------------------------------------------------------------------------
end module soca_utils
