!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bkgerrgodas_mod
  use config_mod
  use datetime_mod  
  use iso_c_binding
  use kinds
  use soca_bkgerrutil_mod  
  use soca_fields
  use soca_model_geom_type, only : geom_get_domain_indices
  use soca_utils
  use soca_omb_stats_mod
  use fckit_mpi_module
  use tools_const, only : pi

  implicit none

  !> Fortran derived type to hold configuration D
  type :: soca_bkgerrgodas_config
     type(soca_field),         pointer :: bkg
     type(soca_field)                  :: std_bkgerr
     type(soca_bkgerr_bounds_type)     :: bounds         ! Bounds for bkgerrgodas
     real(kind=kind_real)              :: t_dz           ! For rescaling of the vertical gradient
     real(kind=kind_real)              :: t_efold        ! E-folding scale for surf based T min
     real(kind=kind_real)              :: ssh_phi_ex
     integer                           :: isc, iec, jsc, jec
  end type soca_bkgerrgodas_config

#define LISTED_TYPE soca_bkgerrgodas_config

  !> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: soca_bkgerrgodas_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "oops/util/linkedList_c.f"
  ! ------------------------------------------------------------------------------
  !> Setup the static background error
  subroutine soca_bkgerrgodas_setup(c_conf, self, bkg)
    type(soca_bkgerrgodas_config), intent(inout) :: self
    type(soca_field),         target, intent(in) :: bkg
    type(c_ptr),                      intent(in) :: c_conf

    type(datetime) :: vdate
    character(len=800) :: fname = 'soca_bkgerrgodas.nc'

    ! Allocate memory for bkgerrgodasor
    call create_copy(self%std_bkgerr, bkg)

    ! Get bounds from configuration
    call self%bounds%read(c_conf)

    ! get parameters not already included in self%bounds
    self%t_dz   = config_get_real(c_conf,"t_dz")
    self%t_efold   = config_get_real(c_conf,"t_efold")
    self%ssh_phi_ex = config_get_real(c_conf, "ssh_phi_ex")

    ! Associate background
    self%bkg => bkg

    ! Set all fields to zero
    call zeros(self%std_bkgerr)

    ! Std of bkg error for T/S/SSH based on background.
    ! S and SSH error are only for the unbalanced portion of S/SSH
    call soca_bkgerrgodas_tocn(self)
    call soca_bkgerrgodas_socn(self)
    call soca_bkgerrgodas_ssh(self)

    ! Apply config bounds to background error
    call self%bounds%apply(self%std_bkgerr)

    ! Save
    call soca_fld2file(self%std_bkgerr, fname)
    
  end subroutine soca_bkgerrgodas_setup

  ! ------------------------------------------------------------------------------
  !> Apply background error: dxm = D dxa
  subroutine soca_bkgerrgodas_mult(self, dxa, dxm)
    type(soca_bkgerrgodas_config),    intent(in) :: self    
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

  end subroutine soca_bkgerrgodas_mult

  ! ------------------------------------------------------------------------------
  !> Derive T background error from vertial gradient of temperature
  subroutine soca_bkgerrgodas_tocn(self)
    type(soca_bkgerrgodas_config),     intent(inout) :: self

    real(kind=kind_real), allocatable :: sig1(:), sig2(:)
    type(soca_domain_indices), pointer :: domain
    integer :: is, ie, js, je, i, j, k
    integer :: ins, ns = 1, iter, niter = 1
    type(soca_omb_stats) :: sst

    ! Get compute domain indices
    call geom_get_domain_indices(self%bkg%geom%ocean, "compute", &
         &domain%is, domain%ie, domain%js, domain%je)

    ! Get local compute domain indices
    call geom_get_domain_indices(self%bkg%geom%ocean, "compute", &
         &domain%isl, domain%iel, domain%jsl, domain%jel, local=.true.)

    ! Allocate temporary arrays
    allocate(sig1(self%bkg%geom%ocean%nzo), sig2(self%bkg%geom%ocean%nzo))
       
    ! Initialize sst background error to previously computed std of omb's
    ! Currently hard-coded to read GODAS file
    call sst%init(domain)
    call sst%bin(self%bkg%geom%ocean%lon, self%bkg%geom%ocean%lat)
    
    ! Loop over compute domain
    do i = domain%is, domain%ie
       do j = domain%js, domain%je
          if (self%bkg%geom%ocean%mask2d(i,j).eq.1) then

             ! Step 1: sigb from dT/dz
             call soca_diff(sig1(:), self%bkg%tocn(i,j,:), self%bkg%hocn(i,j,:))
             sig1(:) = self%t_dz * abs(sig1) ! Rescale background error

             ! Step 2: sigb based on efolding scale
             sig2(:) = self%bounds%t_min + (sst%bgerr_model(i,j)-self%bounds%t_min)*&
                  &exp((self%bkg%layer_depth(i,j,1)-self%bkg%layer_depth(i,j,:))&
                    &/self%t_efold)

             ! Step 3: sigb = max(sig1, sig2)
             do k = 1, self%bkg%geom%ocean%nzo
                self%std_bkgerr%tocn(i,j,k) = min( max(sig1(k), sig2(k)), &
                  & self%bounds%t_max)
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

          end if
       end do
    end do

    ! Release memory
    call sst%exit()
    deallocate(sig1, sig2)

  end subroutine soca_bkgerrgodas_tocn

  ! ------------------------------------------------------------------------------
  !> Derive unbalanced SSH background error, based on latitude
  subroutine soca_bkgerrgodas_ssh(self)
    type(soca_bkgerrgodas_config),     intent(inout) :: self
    type(soca_domain_indices), target :: domain
    integer :: i, j, k

    ! Get compute domain indices
    call geom_get_domain_indices(self%bkg%geom%ocean, "compute", &
         &domain%is, domain%ie, domain%js, domain%je)

    ! Loop over compute domain
    do i = domain%is, domain%ie
       do j = domain%js, domain%je
            if (self%bkg%geom%ocean%mask2d(i,j) .ne. 1) cycle

            if ( abs(self%bkg%geom%ocean%lat(i,j)) >= self%ssh_phi_ex) then
              ! if in extratropics, set to max value
              self%std_bkgerr%ssh(i,j) = self%bounds%ssh_max
            else
              ! otherwise, taper to min value (0.0) toward equator
              self%std_bkgerr%ssh(i,j) = self%bounds%ssh_min + 0.5 * &
                (self%bounds%ssh_max - self%bounds%ssh_min) * &
                (1 - cos(pi * self%bkg%geom%ocean%lat(i,j) / self%ssh_phi_ex))
            end if
       end do
    end do
  end subroutine soca_bkgerrgodas_ssh


  ! ------------------------------------------------------------------------------
  !> Derive unbalanced S background error, based on MLD
  subroutine soca_bkgerrgodas_socn(self)
    type(soca_bkgerrgodas_config),     intent(inout) :: self
    !
    real(kind=kind_real), allocatable :: dsdz(:), dtdz(:)
    type(soca_domain_indices), target :: domain
    integer :: is, ie, js, je, i, j, k
    real(kind=kind_real) :: r


    ! Get compute domain indices
    call geom_get_domain_indices(self%bkg%geom%ocean, "compute", &
         &domain%is, domain%ie, domain%js, domain%je)
    !
    ! Get local compute domain indices
    call geom_get_domain_indices(self%bkg%geom%ocean, "compute", &
         &domain%isl, domain%iel, domain%jsl, domain%jel, local=.true.)

    ! TODO read in a precomputed surface S background error

    ! Loop over compute domain
    do i = domain%is, domain%ie
      do j = domain%js, domain%je
        if (self%bkg%geom%ocean%mask2d(i,j) /= 1)  cycle

        do k = 1, self%bkg%geom%ocean%nzo
          if ( self%bkg%layer_depth(i,j,k) <= self%bkg%mld(i,j)) then
            ! if in the mixed layer, set to the maximum value
            self%std_bkgerr%socn(i,j,k) = self%bounds%s_max
          else
            ! otherwise, taper to the minium value below MLD
            r = 0.1 + 0.45 * (1-tanh( 2 * log( &
              & self%bkg%layer_depth(i,j,k) / self%bkg%mld(i,j) )))
            self%std_bkgerr%socn(i,j,k) = max(self%bounds%s_min, r*self%bounds%s_max)
          end if
        end do
      end do
    end do
  end subroutine soca_bkgerrgodas_socn

end module soca_bkgerrgodas_mod
