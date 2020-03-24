! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_utils

use netcdf
use kinds, only: kind_real
use gsw_mod_toolbox, only : gsw_rho, gsw_sa_from_sp, gsw_ct_from_pt, gsw_mlp
use fckit_kdtree_module, only: kdtree, kdtree_create, kdtree_destroy, &
                               kdtree_k_nearest_neighbors
use fckit_geometry_module, only: sphere_distance
use fckit_exception_module, only: fckit_exception

implicit none

private
public :: write2pe, soca_str2int, soca_adjust, &
          soca_rho, soca_diff, soca_mld, nc_check, soca_remap_idw

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

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
     dvdz(ik) = (v(ik+1)-v(ik-1))/(h(ik)+0.5*h(ik+1)+h(ik-1))
  end do
  dvdz(1) = dvdz(2)
  dvdz(k) = dvdz(k-1)

end subroutine soca_diff

! ------------------------------------------------------------------------------

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
! inverse distance weighted remaping (modified Shepard's method)
subroutine soca_remap_idw(lon_src, lat_src, data_src, lon_dst, lat_dst, data_dst)
  real(kind_real), intent(in) :: lon_src(:)
  real(kind_real), intent(in) :: lat_src(:)
  real(kind_real), intent(in) :: data_src(:)
  real(kind_real), intent(in) :: lon_dst(:,:)
  real(kind_real), intent(in) :: lat_dst(:,:)
  real(kind_real), intent(inout) :: data_dst(:,:)

  integer, parameter :: nn_max = 10
  real(kind_real), parameter :: idw_pow = 2.0

  integer :: idx(nn_max)
  integer :: n_src, i, j, n, nn
  real(kind_real) :: dmax, r, w(nn_max),  dist(nn_max)
  type(kdtree) :: kd

  ! create kd tree
  n_src = size(lon_src)
  kd = kdtree_create(n_src, lon_src, lat_src)

  ! remap
  do i = 1, size(data_dst, dim=1)
    do j = 1, size(data_dst, dim=2)

      ! get nn_max nearest neighbors
      call kdtree_k_nearest_neighbors(kd, lon_dst(i,j), lat_dst(i,j), nn_max, idx)

      ! get distances. Add a small offset so there is never any 0 values
      do n=1,nn_max
        dist(n) = sphere_distance(lon_dst(i,j), lat_dst(i,j), &
                                  lon_src(idx(n)), lat_src(idx(n)))
      end do
      dist = dist + 1e-6

      ! truncate the list if the last points are the same distance.
      ! This is needed to ensure reproducibility across machines.
      ! The last point is always removed (becuase we don't know if it would
      ! have been identical to the one after it)
      nn=nn_max-1
      do n=nn_max-1, 1, -1
        if (dist(n) /= dist(nn_max)) exit
        nn = n-1
      end do
      if (nn <= 0 ) call fckit_exception%abort( &
        "No valid points found in IDW remapping, uh oh.")

      ! calculate weights based on inverse distance
      dmax = maxval(dist(1:nn))
      w = 0.0
      do n=1,nn
        w(n) = ((dmax-dist(n)) / (dmax*dist(n))) ** idw_pow
      end do
      w = w / sum(w)

      ! calculate final value
      r = 0.0
      do n=1,nn
        r = r + data_src(idx(n))*w(n)
      end do
      data_dst(i,j) = r

    end do
  end do

  ! done, cleanup
  call kdtree_destroy(kd)
end subroutine soca_remap_idw

! ------------------------------------------------------------------------------
end module soca_utils
