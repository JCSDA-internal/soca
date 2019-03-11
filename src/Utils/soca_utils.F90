! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_utils
  use kinds
  
  implicit none

  private
  public :: write2pe, htoz, soca_str2int
  public :: soca_rho, soca_diff, soca_mld, nc_check

contains

  ! ------------------------------------------------------------------------------

  elemental function soca_rho(sp, pt, p, lon, lat)
    use kinds
    use gsw_mod_toolbox, only : gsw_rho, gsw_sa_from_sp, gsw_ct_from_pt
    real(kind=kind_real), intent(in)  :: pt, sp, p, lon, lat
    real(kind=kind_real) :: sa, ct, lon_rot, soca_rho, soca_mld

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
    use kinds
    use gsw_mod_toolbox, only : gsw_rho, gsw_sa_from_sp, gsw_ct_from_pt, gsw_mlp
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
    use kinds

    implicit none

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

  subroutine htoz(h, z)
    real(kind=kind_real),  intent(in) :: h(:) ! Layer thickness
    real(kind=kind_real), intent(out) :: z(:) ! Mid-layer depth

    integer :: nlev, ilev

    nlev = size(h,1)
    z(1)=0.5_kind_real*h(1)
    do ilev = 2, nlev
       z(ilev)=sum(h(1:ilev-1))+0.5_kind_real*h(ilev)
    end do
  end subroutine htoz
    
  ! ------------------------------------------------------------------------------

  subroutine write2pe(vec,varname,filename,append)

    use netcdf
    use kinds
    
    implicit none

    real(kind=kind_real), intent(in) :: vec(:)
    character(len=256),   intent(in) :: varname
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

    use netcdf
    IMPLICIT NONE
    integer(4), intent ( in) :: status
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine nc_check

  ! ------------------------------------------------------------------------------  
  
  subroutine soca_str2int(str, int)
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int

    read(str,*)  int
  end subroutine soca_str2int  

end module soca_utils
