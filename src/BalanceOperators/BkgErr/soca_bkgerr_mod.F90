!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bkgerr_mod

  use kinds
  implicit none

  !> Fortran derived type to hold configuration D
  type :: soca_bkgerr_config
     type(soca_field),         pointer :: bkg
     real(kind=kind_real), allocatable :: z(:,:,:)
     real(kind=kind_real), allocatable :: sig_temp(:,:,:)
     real(kind=kind_real), allocatable :: sig_salt(:,:,:)
     real(kind=kind_real), allocatable :: sig_ssh(:,:,:)
     integer              :: isc, iec, jsc, jec       !> Compute domain 
  end type soca_bkgerr_config

#define LISTED_TYPE soca_bkgerr_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_bkgerr_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------
  subroutine soca_bkgerr_setup(c_conf, self, bkg)

    use kinds
    use iso_c_binding
    use config_mod
    use soca_fields
    use soca_model_geom_type, only : geom_get_domain_indices

    implicit none

    type(soca_vertconv),   intent(inout) :: self
    type(soca_field), target, intent(in) :: bkg
    type(c_ptr),              intent(in) :: c_conf

    integer :: isc, iec, jsc, jec, i, j, k, nl
    real(kind=kind_real), allocatable :: dvdz(:), v(:), h(:)
    real(kind=kind_real) :: dt, ds, t0, s0, p, h, lon, lat
    real(kind=kind_real) :: detas
    
    nl = size(bkg%hocn,3)

    ! Get configuration for vertical convolution
    ! No config for now!

    ! Store background
    self%bkg => bkg

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(bkg%geom%ocean, "compute", isc, iec, jsc, jec)
    self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec

    allocate(self%sig_temp(isc:iec, jsc:jec, nl))
    allocate(self%sig_salt(isc:iec, jsc:jec, nl))    
    allocate(self%sig_ssh(isc:iec, jsc:jec))

    ! Initialize local ocean depth from layer thickness  
    allocate(self%z(isc:iec, jsc:jec, nl))
    do i = isc, iec
       do j = jsc, jec
          if (bkg%geom%ocean%mask2d(i,j).eq.1) then
             do k = 1, nl
                self%z(i,j,k) = sum(bkg%hocn(i,j,k:))
             end do
          end if
       end do
    end do

    ! Estimate std. dev. of background error for temp and salt
    ! from vertical gradients
    do i = isc, iec
       do j = jsc, jec
          if (bkg%geom%ocean%mask2d(i,j).eq.1) then
             v(:) = bkg%tocn(i,j,:)
             h(:) = bkg%hocn(i,j,:)           
             call soca_diff(dvdz,v,h)
             do k = 1, nl              
                self%sig_temp(i,j,k) =  dvdz(k)
             end do
             v(:)=bkg%socn(i,j,:)
             call soca_diff(dvdz,v,h)
             do k = 1, nl
                self%sig_salt(i,j,k) =  dvdz(k)              
             end do
          end if
       end do
    end do

    ! Derive sig_ssh from dynamic height assumption
    do i = isc, iec
       do j = jsc, jec
          h=traj%hocn(i,j,k)
          p=z
          call soca_steric_tl(deta, dt, ds, tb, sb, p, h,&
               &traj%geom%ocean%lon(i,j), traj%geom%ocean%lat(i,j))
          dxm%ssh(i,j)=dxm%ssh(i,j)+deta
          
          call soca_steric_tl (detas, dt, ds, t0, s0, p, h, lon, lat)
    
  end subroutine soca_bkgerr_setup
  ! ------------------------------------------------------------------------------

end module soca_bkgerr_mod
