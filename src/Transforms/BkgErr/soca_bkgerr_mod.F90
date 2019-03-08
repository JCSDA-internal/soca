!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bkgerr_mod
  use config_mod
  use datetime_mod  
  use iso_c_binding
  use kinds
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices
  use soca_utils
  use soca_omb_stats_mod
  use fckit_mpi_module
  
  implicit none

  type :: soca_bkgerror_bounds
     real(kind=kind_real) :: t_min, t_max, t_ml
     real(kind=kind_real) :: s_min, s_max
     real(kind=kind_real) :: ssh_min, ssh_max
  end type soca_bkgerror_bounds
  
  !> Fortran derived type to hold configuration D
  type :: soca_bkgerr_config
     type(soca_field),         pointer :: bkg
     type(soca_field)                  :: std_bkgerr
     type(soca_bkgerror_bounds)        :: bounds
     real(kind=kind_real)              :: delta_z     
     real(kind=kind_real)              :: efold_z
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
  !> Setup the static background error
  subroutine soca_bkgerr_setup(c_conf, self, bkg)
    type(soca_bkgerr_config), intent(inout) :: self
    type(soca_field),    target, intent(in) :: bkg
    type(c_ptr),                 intent(in) :: c_conf

    integer :: isc, iec, jsc, jec, i, j, k, nl
    real(kind=kind_real), allocatable :: dvdz(:), v(:), h(:)
    real(kind=kind_real) :: dt, ds, t0, s0, p, lon, lat
    real(kind=kind_real) :: detas, efold
    type(datetime) :: vdate
    character(len=800) :: fname = 'soca_bkgerror.nc'
    logical :: read_from_file = .false.

    ! Get number of ocean levels
    nl = size(bkg%hocn,3)

    ! Allocate memory for bkgerror
    call create_copy(self%std_bkgerr, bkg)
    if (config_element_exists(c_conf,"ocn_filename")) then
       read_from_file = .true.
       
       ! Read variance       
       call read_file(self%std_bkgerr, c_conf, vdate)

       ! Convert to standard deviation       
       self%std_bkgerr%tocn = sqrt(self%std_bkgerr%tocn)
       self%std_bkgerr%socn = sqrt(self%std_bkgerr%socn)       
       self%std_bkgerr%ssh = sqrt(self%std_bkgerr%ssh)
    else
       read_from_file = .false.
       self%delta_z   = config_get_real(c_conf,"delta_z")
    end if

    ! Vertical e-folding scale
    self%efold_z   = config_get_real(c_conf,"efold_z")

    ! Get bounds from configuration
    self%bounds%t_min   = config_get_real(c_conf,"t_min")
    self%bounds%t_max   = config_get_real(c_conf,"t_max")
    if (config_element_exists(c_conf,"t_ml")) then
       self%bounds%t_ml = config_get_real(c_conf,"t_ml")
    else
       self%bounds%t_ml = 0.5_kind_real
    endif    
    self%bounds%s_min   = config_get_real(c_conf,"s_min")
    self%bounds%s_max   = config_get_real(c_conf,"s_max")
    self%bounds%ssh_min = config_get_real(c_conf,"ssh_min")
    self%bounds%ssh_max = config_get_real(c_conf,"ssh_max")

    ! Associate background
    self%bkg => bkg

    ! Indices for compute domain (no halo)
    call geom_get_domain_indices(bkg%geom%ocean, "compute", isc, iec, jsc, jec)
    self%isc=isc; self%iec=iec; self%jsc=jsc; self%jec=jec

    ! Std of bkg error for temperature based on dT/dz
    if (.not.read_from_file) then
       call soca_bkgerr_tocn(self)
    end if
    
    ! Apply config bounds to background error
    do i = isc, iec
       do j = jsc, jec
          ! Ocean
          self%std_bkgerr%ssh(i,j) = adjusted_std(abs(self%std_bkgerr%ssh(i,j)), &
               &self%bounds%ssh_min,&
               &self%bounds%ssh_max)
          do k = 1, nl
             efold = exp(-bkg%layer_depth(i,j,k)/self%efold_z)
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
  !> Apply background error: dxm = D dxa
  subroutine soca_bkgerr_mult(self, dxa, dxm)
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

  ! ------------------------------------------------------------------------------
  !> Apply bounds
  elemental function adjusted_std(std, minstd, maxstd)
    
    implicit none
    
    real(kind=kind_real), intent(in)  :: std, minstd, maxstd
    real(kind=kind_real) :: adjusted_std
    
    adjusted_std = min( max(std, minstd), maxstd)
    
  end function adjusted_std

  ! ------------------------------------------------------------------------------
  !> Derive background error from vertial gradient of temperature 
  subroutine soca_bkgerr_tocn(self)
    type(soca_bkgerr_config), intent(inout) :: self

    real(kind=kind_real), allocatable :: temp(:), sig1(:), sig2(:)
    type(soca_domain_indices), target :: domain    
    integer :: is, ie, js, je, i, j, k
    integer :: ins, ns = 1, iter, niter = 1
    type(soca_omb_stats) :: sst

    ! Set all fields to zero
    call zeros(self%std_bkgerr)

    ! Get compute domain indices
    call geom_get_domain_indices(self%bkg%geom%ocean, "compute", &
         &domain%is, domain%ie, domain%js, domain%je)
    
    ! Get local compute domain indices
    call geom_get_domain_indices(self%bkg%geom%ocean, "compute", &
         &domain%isl, domain%iel, domain%jsl, domain%jel, local=.true.)

       
    ! Allocate temporary arrays
    allocate(temp(self%bkg%geom%ocean%nzo))
    allocate(sig1(self%bkg%geom%ocean%nzo), sig2(self%bkg%geom%ocean%nzo))

       
    ! Initialize sst background error to previously computed std of omb's
    ! Currently hard-coded to read GODAS file
    call sst%init(domain)
    call sst%bin(self%bkg%geom%ocean%lon, self%bkg%geom%ocean%lat)

    ! Loop over compute domain
    do i = domain%is, domain%ie
       do j = domain%js, domain%je
          if (self%bkg%geom%ocean%mask2d(i,j).eq.1) then

             ! Make sure Temp values in thin layers are realistic
             temp(:) = self%bkg%tocn(i,j,:)
             !call soca_clean_vertical(self%bkg%hocn(i,j,:), temp(:))

             ! Step 1: sigb from dT/dz 
             call soca_diff(sig1(:), temp(:), self%bkg%hocn(i,j,:))
             sig1(:) = self%delta_z * abs(sig1) ! Rescale background error

             ! Step 2: sigb based on efolding scale
             sig2(:) = sst%bgerr_model(i,j)*&
                  &exp(-self%bkg%layer_depth(i,j,:)/self%efold_z)

             ! Step 3: sigb = max(sig1, sig2)
             do k = 1, self%bkg%geom%ocean%nzo
                self%std_bkgerr%tocn(i,j,k) = max(sig1(k), sig2(k))
             end do
             
             ! Step 4: Vertical smoothing
             do iter = 1, niter
                do k = 2, self%bkg%geom%ocean%nzo-1
                   self%std_bkgerr%tocn(i,j,k) = &
                        &( self%std_bkgerr%tocn(i,j,k-1)*self%bkg%hocn(i,j,k-1) +&
                        &  self%std_bkgerr%tocn(i,j,k)*self%bkg%hocn(i,j,k) +&
                        &  self%std_bkgerr%tocn(i,j,k+1)*self%bkg%hocn(i,j,k+1) )/&
                        & (sum(self%bkg%hocn(i,j,k-1:k+1)))
                end do
             end do

             ! Apply bounds from yaml config
             do k = 1, self%bkg%geom%ocean%nzo
                if (self%bkg%layer_depth(i,j,k)<self%bkg%mld(i,j)) then
                   self%std_bkgerr%tocn(i,j,k) = max(self%std_bkgerr%tocn(i,j,k),&
                        &self%bounds%t_ml)
                else
                   self%std_bkgerr%tocn(i,j,k) = max(self%std_bkgerr%tocn(i,j,k),&
                        &self%bounds%t_min)                   
                end if
             end do
             
          end if
       end do
    end do

    ! Release memory
    call sst%exit()
    deallocate(temp, sig1, sig2)

  end subroutine soca_bkgerr_tocn

  ! ------------------------------------------------------------------------------
  
end module soca_bkgerr_mod


