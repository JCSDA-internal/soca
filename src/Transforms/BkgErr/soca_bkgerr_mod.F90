!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bkgerr_mod

  use kinds
  use soca_fields
  
  implicit none

  !> Fortran derived type to hold configuration D
  type :: soca_bkgerr_config
     type(soca_field),         pointer :: bkg
     type(soca_field)                  :: std_bkgerr     
     real(kind=kind_real), allocatable :: z(:,:,:)
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
    use soca_kst_mod
    use datetime_mod
    use mpi
    
    implicit none

    type(soca_bkgerr_config), intent(inout) :: self
    type(soca_field),    target, intent(in) :: bkg
    type(c_ptr),                 intent(in) :: c_conf

    integer :: isc, iec, jsc, jec, i, j, k, nl
    real(kind=kind_real), allocatable :: dvdz(:), v(:), h(:)
    real(kind=kind_real) :: dt, ds, t0, s0, p, lon, lat
    real(kind=kind_real) :: detas
    type(datetime) :: vdate
    
    nl = size(bkg%hocn,3)

    ! Read background error
    call create_copy(self%std_bkgerr, bkg)
    call read_file(self%std_bkgerr, c_conf, vdate)
    
    ! Store background
    self%bkg => bkg

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(bkg%geom%ocean, "compute", isc, iec, jsc, jec)
    self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec

    ! Initialize local ocean depth from layer thickness  
    call bkg%geom%ocean%thickness2depth(bkg%hocn, self%z)
    
    ! Limit background error
    do i = isc, iec
       do j = jsc, jec
          !if (self%bkg%geom%ocean%mask2d(i,j).eq.1) then
             ! ocean
             self%std_bkgerr%ssh(i,j) = adjusted_std(abs(self%std_bkgerr%ssh(i,j)), 0.1d0, 10.0d0)
             do k = 1, nl
                self%std_bkgerr%tocn(i,j,k) = 1.0*exp(-self%z(i,j,k)/300d0)
                    !adjusted_std(abs(self%std_bkgerr%tocn(i,j,k)), 0.d0, 2.0d0)
                self%std_bkgerr%socn(i,j,k) = 0.2*exp(-self%z(i,j,k)/300d0)
                    !adjusted_std(abs(self%std_bkgerr%socn(i,j,k)), 0.0d0, 0.1d0)
             end do
             ! sea-ice
             self%std_bkgerr%cicen(i,j,:) = adjusted_std(abs(self%std_bkgerr%cicen(i,j,:)), 0.01d0, 0.5d0)
             self%std_bkgerr%hicen(i,j,:) = adjusted_std(abs(self%std_bkgerr%hicen(i,j,:)), 10d0, 100.0d0)             
          !end if
       end do
    end do
    
!!$    ! Setup arrays of backgroud error std dev
!!$    allocate(self%sig_temp(isc:iec, jsc:jec, nl))
!!$    allocate(self%sig_salt(isc:iec, jsc:jec, nl))    
!!$    allocate(self%sig_ssh(isc:iec, jsc:jec))
!!$
!!$    ! Initialize local ocean depth from layer thickness  
!!$    call bkg%geom%ocean%thickness2depth(bkg%hocn, self%z)
!!$
!!$    ! Estimate std. dev. of background error for temp and salt
!!$    ! from vertical gradients
!!$    allocate(h(nl),v(nl),dvdz(nl))
!!$    do i = isc, iec
!!$       do j = jsc, jec
!!$          if (self%bkg%geom%ocean%mask2d(i,j).eq.1) then
!!$             v(:) = self%bkg%tocn(i,j,:)
!!$             h(:) = self%bkg%hocn(i,j,:)           
!!$             call soca_diff(dvdz,v,h)
!!$             do k = 1, nl              
!!$                self%sig_temp(i,j,k) =  1.0d0 !min(abs(dvdz(k)),5.0d0)
!!$             end do
!!$             v(:)=bkg%socn(i,j,:)
!!$             call soca_diff(dvdz,v,h)
!!$             do k = 1, nl
!!$                self%sig_salt(i,j,k) =  0.1d0 !min(abs(dvdz(k)),5.0d0)              
!!$             end do
!!$          end if
!!$       end do
!!$    end do
!!$
!!$    deallocate(h,v,dvdz)
!!$
!!$    ! Derive sig_ssh from dynamic height assumption
!!$    ! TODO !!!!!!!!!!!!
!!$    !self%sig_ssh(isc:iec, jsc:jec) = 0.1d0
!!$    self%sig_ssh = 0.1d0    
    
  end subroutine soca_bkgerr_setup

  ! ------------------------------------------------------------------------------

  subroutine soca_bkgerr_mult(self, dxa, dxm)

    use kinds
    use soca_model_geom_type, only : geom_get_domain_indices

    implicit none

    type(soca_bkgerr_config),    intent(in) :: self    
    type(soca_field),            intent(in) :: dxa
    type(soca_field),         intent(inout) :: dxm

    integer :: isc, iec, jsc, jec, i, j, k
    
    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(self%bkg%geom%ocean, "compute", isc, iec, jsc, jec)

    do i = isc, iec
       do j = jsc, jec
          if (self%bkg%geom%ocean%mask2d(i,j).eq.1) then
             dxm%ssh(i,j) = self%std_bkgerr%ssh(i,j) * dxa%ssh(i,j)
             dxm%tocn(i,j,:) = self%std_bkgerr%tocn(i,j,:) * dxa%tocn(i,j,:)
             dxm%socn(i,j,:) = self%std_bkgerr%socn(i,j,:)  * dxa%socn(i,j,:)

             dxm%cicen(i,j,:) =  self%std_bkgerr%cicen(i,j,:) * dxa%cicen(i,j,:)
             dxm%hicen(i,j,:) =  self%std_bkgerr%hicen(i,j,:) * dxa%hicen(i,j,:)
          end if
       end do
    end do

  end subroutine soca_bkgerr_mult

  elemental function adjusted_std(std, minstd, maxstd)
    
    implicit none
    
    real(kind=kind_real), intent(in)  :: std, minstd, maxstd
    real(kind=kind_real) :: adjusted_std
    
    adjusted_std = min( max(std, minstd), maxstd)
    
  end function adjusted_std
  
end module soca_bkgerr_mod


