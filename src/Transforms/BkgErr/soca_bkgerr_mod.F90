!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bkgerr_mod

  use kinds
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices
  
  implicit none

  type :: soca_bkgerror_bounds
     real(kind=kind_real) :: t_min, t_max
     real(kind=kind_real) :: s_min, s_max
     real(kind=kind_real) :: ssh_min, ssh_max
  end type soca_bkgerror_bounds
  
  !> Fortran derived type to hold configuration D
  type :: soca_bkgerr_config
     type(soca_field),         pointer :: bkg
     type(soca_field)                  :: std_bkgerr
     type(soca_bkgerror_bounds)        :: bounds
     integer                           :: isc, iec, jsc, jec
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
    character(len=800) :: fname = 'soca_bkgerror.nc'
    logical :: read_from_file = .false.
    
    nl = size(bkg%hocn,3)

    ! Read background error
    call create_copy(self%std_bkgerr, bkg)
    if (config_element_exists(c_conf,"ocn_filename")) then
       read_from_file = .true.
       call read_file(self%std_bkgerr, c_conf, vdate)
    else
       read_from_file = .false.       
    end if       

    ! Get bounds from configuration
    self%bounds%t_min   = config_get_real(c_conf,"t_min")
    self%bounds%t_max   = config_get_real(c_conf,"t_max")
    self%bounds%s_min   = config_get_real(c_conf,"s_min")
    self%bounds%s_max   = config_get_real(c_conf,"s_max")    
    self%bounds%ssh_min = config_get_real(c_conf,"ssh_min")
    self%bounds%ssh_max = config_get_real(c_conf,"ssh_max")
    
    ! Store background
    self%bkg => bkg

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(bkg%geom%ocean, "compute", isc, iec, jsc, jec)
    self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec

    ! Std of bkg error for temperature based on dT/dz
    if (.not.read_from_file) call soca_bkgerr_tocn(self)
    
    ! Apply config bounds to background error
    do i = isc, iec
       do j = jsc, jec
          ! Ocean
          self%std_bkgerr%ssh(i,j) = adjusted_std(abs(self%std_bkgerr%ssh(i,j)), &
               &self%bounds%ssh_min,&
               &self%bounds%ssh_max)
          do k = 1, nl
             self%std_bkgerr%tocn(i,j,k) = adjusted_std(abs(self%std_bkgerr%tocn(i,j,k)),&
               &self%bounds%t_min,&
               &self%bounds%t_max)
             self%std_bkgerr%socn(i,j,k) = adjusted_std(abs(self%std_bkgerr%socn(i,j,k)),&
               &self%bounds%s_min,&
               &self%bounds%s_max)
          end do

          ! Sea-ice
          self%std_bkgerr%cicen(i,j,:) = adjusted_std(abs(self%std_bkgerr%cicen(i,j,:)), 0.01d0, 0.5d0)
          self%std_bkgerr%hicen(i,j,:) = adjusted_std(abs(self%std_bkgerr%hicen(i,j,:)), 10d0, 100.0d0)
       end do
    end do

    ! Save filtered background error
    call soca_fld2file(self%std_bkgerr, fname)
    
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

  ! ------------------------------------------------------------------------------
  
  subroutine soca_bkgerr_tocn(self)
    use kinds
    use soca_utils
    
    implicit none
    type(soca_bkgerr_config), intent(inout) :: self

    real(kind=kind_real), allocatable :: temp(:), vmask(:)
    real(kind=kind_real) :: jac(2)
    real(kind=kind_real) :: delta_z = 10.0_kind_real ![m] TODO: Move to config
    integer :: is, ie, js, je, i, j, k

    ! Set all fields to zero
    call zeros(self%std_bkgerr)
    !call ones(self%std_bkgerr)    
    !call self_mul(self%std_bkgerr, 1d-8)
    
    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(self%bkg%geom%ocean, "compute", is, ie, js, je)

    ! Compute temperature gradient
    allocate(temp(self%bkg%geom%ocean%nzo))
    allocate(vmask(self%bkg%geom%ocean%nzo))

    do i = is, ie
       do j = js, je
          if (self%bkg%geom%ocean%mask2d(i,j).eq.1) then

             ! Make sure Temp values in thin layers are realistic
             temp(:) = self%bkg%tocn(i,j,:)
             call soca_clean_vertical(self%bkg%hocn(i,j,:), temp(:))
             
             ! T Bkg error from dT/dz 
             call soca_diff(self%std_bkgerr%tocn(i,j,:),&
                           &self%bkg%tocn(i,j,:),&
                           &self%bkg%hocn(i,j,:))

             ! Apply vertical mask
             vmask = 1.0_kind_real
             where (self%bkg%hocn(i,j,:)<1.0d-6)
                vmask = 0.0_kind_real
             end where

             ! Scale background error
             self%std_bkgerr%tocn(i,j,:) = abs(delta_z * &
                  &(self%std_bkgerr%tocn(i,j,:) * vmask))
             
          end if
       end do
    end do

    ! Make up background error for salt and ssh
    self%std_bkgerr%socn = 0.1_kind_real * self%std_bkgerr%tocn
    self%std_bkgerr%ssh = 0.1_kind_real    
    
    ! Release memory
    deallocate(temp, vmask)
    
  end subroutine soca_bkgerr_tocn

  ! ------------------------------------------------------------------------------
  
end module soca_bkgerr_mod


