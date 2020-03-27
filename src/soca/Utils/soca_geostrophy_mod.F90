! (C) Copyright 2017-2020 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


module soca_geostrophy_mod
  use kinds
  use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad
  use soca_geom_mod, only : soca_geom
  use tools_const, only: deg2rad, req
  use soca_utils, only : soca_rho
  use soca_ksshts_mod

  implicit none

  private

  !> Coordinate transform (logical <-> spherical) for tracer grid points only
  !> Fortran derived type for differentiation operators
  type, public :: soca_geostrophy_type
     real(kind=kind_real), allocatable, dimension(:,:,:) :: Jac
     real(kind=kind_real), allocatable, dimension(:,:,:) :: dzdk

     real(kind=kind_real), allocatable, dimension(:,:) :: didlon
     real(kind=kind_real), allocatable, dimension(:,:) :: didlat
     real(kind=kind_real), allocatable, dimension(:,:) :: djdlon
     real(kind=kind_real), allocatable, dimension(:,:) :: djdlat
     real(kind=kind_real), allocatable, dimension(:,:,:) :: dkdlon
     real(kind=kind_real), allocatable, dimension(:,:,:) :: dkdlat
     real(kind=kind_real), allocatable, dimension(:,:) :: dkdz

     real(kind=kind_real), allocatable, dimension(:,:) :: icorio
     real(kind=kind_real), allocatable, dimension(:,:) :: dxwgt
     real(kind=kind_real), allocatable, dimension(:,:) :: dywgt

   contains
     procedure :: setup => soca_geostrophy_setup
     procedure :: delete => soca_geostrophy_delete
     procedure :: grad => soca_grad
     procedure :: tl => soca_geostrophy_tl
  end type soca_geostrophy_type

  ! ------------------------------------------------------------------------------
contains

  ! ------------------------------------------------------------------------------
  !> Setup transformation: logical grid <=> physical curvilinear grid
  subroutine soca_geostrophy_setup(self, geom, h)
    class(soca_geostrophy_type), intent(inout) :: self
    type(soca_geom),            intent(in) :: geom
    real(kind=kind_real), dimension(:,:,:) :: h  ! Layer thicknesses

    real(kind=kind_real), allocatable, dimension(:,:) :: dlondi
    real(kind=kind_real), allocatable, dimension(:,:) :: dlondj
    real(kind=kind_real), allocatable, dimension(:,:) :: dlatdi
    real(kind=kind_real), allocatable, dimension(:,:) :: dlatdj
    real(kind=kind_real), allocatable, dimension(:,:,:) :: dzdi
    real(kind=kind_real), allocatable, dimension(:,:,:) :: dzdj
    real(kind=kind_real), allocatable, dimension(:,:,:) :: z
    real(kind=kind_real), allocatable, dimension(:,:) :: wb
    real(kind=kind_real) :: Lb
    integer :: i, j, k

    ! Allocate inverse transformation metrics
    allocate(dlondi(geom%isd:geom%ied,geom%jsd:geom%jed))
    allocate(dlondj(geom%isd:geom%ied,geom%jsd:geom%jed))
    allocate(dlatdi(geom%isd:geom%ied,geom%jsd:geom%jed))
    allocate(dlatdj(geom%isd:geom%ied,geom%jsd:geom%jed))
    allocate(dzdi(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%nzo))
    allocate(dzdj(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%nzo))

    ! Allocate transformation metrics
    allocate(self%Jac(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%nzo))
    allocate(self%dzdk(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%nzo))

    allocate(self%didlon(geom%isd:geom%ied,geom%jsd:geom%jed))
    allocate(self%didlat(geom%isd:geom%ied,geom%jsd:geom%jed))
    allocate(self%djdlon(geom%isd:geom%ied,geom%jsd:geom%jed))
    allocate(self%djdlat(geom%isd:geom%ied,geom%jsd:geom%jed))
    allocate(self%dkdlon(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%nzo))
    allocate(self%dkdlat(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%nzo))
    allocate(self%dkdz(geom%isd:geom%ied,geom%jsd:geom%jed))

    ! Allocate weights for zonal and meridional derivatives
    allocate(self%dxwgt(geom%isd:geom%ied,geom%jsd:geom%jed))
    allocate(self%dywgt(geom%isd:geom%ied,geom%jsd:geom%jed))

    ! Allocate wf/f^-1
    allocate(self%icorio(geom%isd:geom%ied,geom%jsd:geom%jed))

    ! Compute horizontal inverse transformation metrics
    do j = geom%jsc, geom%jec
       do i = geom%isc, geom%iec
          dlondi(i,j) = 0.5*deg2rad*(geom%lon(i+1,j)-geom%lon(i-1,j))
          dlondj(i,j) = 0.5*deg2rad*(geom%lon(i,j+1)-geom%lon(i,j-1))
          dlatdi(i,j) = 0.5*deg2rad*(geom%lat(i+1,j)-geom%lat(i-1,j))
          dlatdj(i,j) = 0.5*deg2rad*(geom%lat(i,j+1)-geom%lat(i,j-1))
       end do
    end do

    ! Allocate and get mid-layer depth
    allocate(z(geom%isd:geom%ied,geom%jsd:geom%jed,1:geom%nzo))
    call geom%thickness2depth(h, z)

    ! Compute horizontal inverse transformation metrics
    do k = 1, geom%nzo
       do j = geom%jsc, geom%jec
          do i = geom%isc, geom%iec
             dzdi(i,j,k) = 0.5*(z(i+1,j,k)-z(i-1,j,k))
             dzdj(i,j,k) = 0.5*(z(i,j+1,k)-z(i,j-1,k))
          end do
       end do
    end do
    self%dzdk = h

    ! Compute transformation metrics
    self%Jac = 0.0
    do k = 1, geom%nzo
       self%Jac(:,:,k) = dlondi(:,:)*dlatdj(:,:)*self%dzdk(:,:,k) - &
                         dlondj(:,:)*dlatdi(:,:)*self%dzdk(:,:,k)
    end do

    ! Needs to be * dzdk(:,:,k)/Jac(:,:,k)
    self%didlon =   dlatdj
    self%didlat = - dlondj

    self%djdlon = - dlatdi
    self%djdlat =   dlondi

    ! Needs to be * 1/Jac(:,:,k)
    do k = 1, geom%nzo
       self%dkdlon(:,:,k) =  dlatdi*dzdj(:,:,k)-dlatdj*dzdi(:,:,k)
       self%dkdlat(:,:,k) = -dlondi*dzdj(:,:,k)+dlondj*dzdi(:,:,k)
    end do
    self%dkdz = dlondi*dlatdj-dlondj*dlatdi

    ! Zero out weights close to equator (no beta-plane proxi yet)
    Lb = 1.55
    allocate(wb(geom%isd:geom%ied,geom%jsd:geom%jed))
    wb  =exp(-geom%lat**2/(2*Lb**2))
    self%icorio = 0.0
    where (geom%lat.ne.0.0)
       self%icorio = (1.0-wb)/(2.0*1035.0*sin(deg2rad*geom%lat)*7.2921e-5)
    end where

    ! Compute weights for zonal and meridional derivatives
    self%dxwgt = 0.0
    where (geom%lat.ne.0.0)
       self%dxwgt = 1.0/(req*cos(deg2rad*geom%lat))
    end where
    self%dywgt = 1.0/req

    ! Fill halo of grid metric parameters
    call mpp_update_domains(self%Jac, geom%Domain%mpp_domain)
    call mpp_update_domains(self%didlon, geom%Domain%mpp_domain)
    call mpp_update_domains(self%didlat, geom%Domain%mpp_domain)
    call mpp_update_domains(self%djdlon, geom%Domain%mpp_domain)
    call mpp_update_domains(self%djdlat, geom%Domain%mpp_domain)
    call mpp_update_domains(self%dkdlon, geom%Domain%mpp_domain)
    call mpp_update_domains(self%dkdlat, geom%Domain%mpp_domain)
    call mpp_update_domains(self%dkdz, geom%Domain%mpp_domain)

  end subroutine soca_geostrophy_setup

  ! ------------------------------------------------------------------------------
  !> Cleanup
  subroutine soca_geostrophy_delete(self)
    class(soca_geostrophy_type), intent(inout) :: self

    deallocate(self%Jac)
    deallocate(self%dzdk)

    deallocate(self%didlon)
    deallocate(self%didlat)
    deallocate(self%djdlon)
    deallocate(self%djdlat)
    deallocate(self%dkdlon)
    deallocate(self%dkdlat)
    deallocate(self%dkdz)

    deallocate(self%icorio)

    deallocate(self%dxwgt)
    deallocate(self%dywgt)

  end subroutine soca_geostrophy_delete

  ! ------------------------------------------------------------------------------
  !> Linear geostrophy  (Loosely based on Weaver, 2005)
  subroutine soca_geostrophy_tl(self, h, t, s, dt, ds, dug, dvg, geom)
    class(soca_geostrophy_type),              intent(in) :: self
    real(kind=kind_real),    allocatable, intent(in) :: h(:,:,:)
    real(kind=kind_real),    allocatable, intent(in) :: dt(:,:,:)
    real(kind=kind_real),    allocatable, intent(in) :: ds(:,:,:)
    real(kind=kind_real),    allocatable, intent(in) :: t(:,:,:)
    real(kind=kind_real),    allocatable, intent(in) :: s(:,:,:)
    real(kind=kind_real), allocatable, intent(inout) :: dug(:,:,:)
    real(kind=kind_real), allocatable, intent(inout) :: dvg(:,:,:)
    type(soca_geom),                      intent(in) :: geom

    integer :: i, j, k, nl
    real(kind=kind_real), allocatable :: dpressure(:,:,:)
    real(kind=kind_real), allocatable :: drhoh(:,:,:)
    real(kind=kind_real), allocatable :: z(:,:,:)
    real(kind=kind_real), allocatable :: mask(:,:), rho0(:,:)
    real(kind=kind_real) :: jac(2)
    integer :: isd, ied, jsd, jed

    ! Indices for compute/data domain
    isd = geom%isd ;  ied = geom%ied ; jsd = geom%jsd; jed = geom%jed
    allocate(dpressure(isd:ied,jsd:jed,1:geom%nzo))
    allocate(drhoh(isd:ied,jsd:jed,1:geom%nzo))
    allocate(z(isd:ied,jsd:jed,1:geom%nzo))
    allocate(mask(isd:ied,jsd:jed))
    allocate(rho0(isd:ied,jsd:jed))

    nl = geom%nzo

    ! Mid-layer depth
    call geom%thickness2depth(h, z)

    ! Density increment scaled by layer thickness at level k
    drhoh = 0.0
    do k = 1, nl
       rho0 = soca_rho(s(:,:,k), t(:,:,k), z(:,:,k), geom%lon, geom%lat)
       do j = geom%jsc, geom%jec
          do i = geom%isc, geom%iec
             call soca_steric_jacobian (jac,&
                                  t(i,j,k), s(i,j,k),z(i,j,k), h(i,j,k),&
                                  geom%lon(i,j),geom%lat(i,j))
             drhoh(i,j,k) = rho0(i,j)*(jac(1)*dt(i,j,k) + jac(2)*ds(i,j,k))
          end do
       end do
    end do

    ! Pressure increment at level k
    dpressure(:,:,nl) = 9.81*drhoh(:,:,nl)
    do k = nl-1, 1, -1
       dpressure(:,:,k) = dpressure(:,:,k+1) + 9.81*drhoh(:,:,k)
    end do

    ! Pressure gradient
    dug = 0.0_kind_real
    dvg = 0.0_kind_real
    call soca_grad(self, dpressure, dug, geom, 'lat')
    call soca_grad(self, dpressure, dvg, geom, 'lon')

    ! Baroclinic currents
    do k = 1,nl
       mask = 1.0
       where( h(:,:,k)<0.1 )
          mask = 0.0
       end where
       dug(:,:,k) = -mask*dug(:,:,k)*self%icorio
       dvg(:,:,k) =  mask*dvg(:,:,k)*self%icorio
    end do

  end subroutine soca_geostrophy_tl

  ! ------------------------------------------------------------------------------
  !> Gradient operator for 3D fields
  subroutine soca_grad(self, fld, dflddxy, geom, direction)
    class(soca_geostrophy_type),              intent(in) :: self
    real(kind=kind_real),    allocatable, intent(in) :: fld(:,:,:)
    real(kind=kind_real), allocatable, intent(inout) :: dflddxy(:,:,:)
    type(soca_geom),                      intent(in) :: geom
    character(len=3),                     intent(in) :: direction ! lon or lat

    integer :: i, j, k, levels
    integer :: isc, iec, jsc, jec
    integer :: isd, ied, jsd, jed
    real(kind=kind_real), allocatable :: fldtmp(:,:,:)

    ! Indices for compute/data domain
    isd = geom%isd ;  ied = geom%ied ; jsd = geom%jsd; jed = geom%jed
    isc = geom%isc ;  iec = geom%iec ; jsc = geom%jsc; jec = geom%jec

    ! Make a temporary copy of fld
    levels = size(fld,3)
    allocate(fldtmp(isd:ied,jsd:jed,0:levels+1))
    fldtmp(:,:,1:levels) = fld(:,:,1:levels)
    fldtmp(:,:,0) = fld(:,:,1)
    fldtmp(:,:,levels+1) = fld(:,:,levels)

    ! Update halo points of input array
    call mpp_update_domains(fldtmp, geom%Domain%mpp_domain, complete=.true.)

    select case (direction)
       case("lon")
          ! Zonal derivative
          dflddxy = 0.0_kind_real
          do k = 1, levels
             do j = jsc, jec
                do i = isc, iec
                   dflddxy(i,j,k) = &
                        self%dzdk(i,j,k)*&
                        ( self%didlon(i,j)*(fldtmp(i+1,j,k)-fldtmp(i-1,j,k)) + &
                        self%djdlon(i,j)*(fldtmp(i,j+1,k)-fldtmp(i,j-1,k)) ) + &
                        self%dkdlon(i,j,k)*(fldtmp(i,j,k+1)-fldtmp(i,j,k-1))
                end do
             end do
             where (self%Jac(:,:,k).ne.0.0)
                dflddxy(:,:,k) = 0.5*geom%mask2d*self%dxwgt*dflddxy(:,:,k)/self%Jac(:,:,k)
             end where
          end do

       case("lat")
          ! Meridional derivative
          dflddxy = 0.0_kind_real
          do k = 1, levels
             do j = jsc, jec
                do i = isc, iec
                   dflddxy(i,j,k) = &
                        self%dzdk(i,j,k)*&
                        ( self%didlat(i,j)*(fldtmp(i+1,j,k)-fldtmp(i-1,j,k)) + &
                        self%djdlat(i,j)*(fldtmp(i,j+1,k)-fldtmp(i,j-1,k)) ) + &
                        self%dkdlat(i,j,k)*(fldtmp(i,j,k+1)-fldtmp(i,j,k-1))
                end do
             end do
             where (self%Jac(:,:,k).ne.0.0)
                dflddxy(:,:,k) = 0.5*geom%mask2d*self%dywgt*dflddxy(:,:,k)/self%Jac(:,:,k)
             end where
          end do
       end select

    ! Update halo
    call mpp_update_domains(dflddxy, geom%Domain%mpp_domain, complete=.true.)

    deallocate(fldtmp)

  end subroutine soca_grad

end module soca_geostrophy_mod
