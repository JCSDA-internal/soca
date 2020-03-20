! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerrgodas_mod

use fckit_configuration_module, only: fckit_configuration
use tools_const, only : pi
use datetime_mod, only: datetime
use kinds, only: kind_real
use soca_fields_mod
use soca_state_mod
use soca_increment_mod
use soca_utils, only: soca_diff
use soca_bkgerrutil_mod, only: soca_bkgerr_bounds_type
use soca_omb_stats_mod, only: soca_omb_stats, soca_domain_indices

implicit none

private
public :: soca_bkgerrgodas_config, &
          soca_bkgerrgodas_setup, soca_bkgerrgodas_mult, &
          soca_bkgerrgodas_tocn, soca_bkgerrgodas_socn, &
          soca_bkgerrgodas_ssh

!> Fortran derived type to hold configuration D
type :: soca_bkgerrgodas_config
   type(soca_state),         pointer :: bkg
   type(soca_fields)                 :: std_bkgerr
   type(soca_bkgerr_bounds_type)     :: bounds         ! Bounds for bkgerrgodas
   real(kind=kind_real)              :: t_dz           ! For rescaling of the vertical gradient
   real(kind=kind_real)              :: t_efold        ! E-folding scale for surf based T min
   real(kind=kind_real)              :: ssh_phi_ex
   integer                           :: isc, iec, jsc, jec
end type soca_bkgerrgodas_config

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup the static background error
subroutine soca_bkgerrgodas_setup(f_conf, self, bkg)
  type(fckit_configuration),        intent(in) :: f_conf
  type(soca_bkgerrgodas_config), intent(inout) :: self
  type(soca_state),         target, intent(in) :: bkg

  type(soca_field), pointer :: field, field_bkg
  integer :: i
  character(len=800) :: fname = 'soca_bkgerrgodas.nc'

  ! Allocate memory for bkgerrgodasor
  call self%std_bkgerr%copy(bkg)
  !call create_copy(self%std_bkgerr, bkg)

  ! Get bounds from configuration
  call self%bounds%read(f_conf)

  ! get parameters not already included in self%bounds
  call f_conf%get_or_die("t_dz", self%t_dz)
  call f_conf%get_or_die("t_efold", self%t_efold)
  call f_conf%get_or_die("ssh_phi_ex", self%ssh_phi_ex)

  ! Associate background
  self%bkg => bkg

  ! Set all fields to zero
  call self%std_bkgerr%zeros()

  ! Std of bkg error for T/S/SSH based on background.
  ! S and SSH error are only for the unbalanced portion of S/SSH
  call soca_bkgerrgodas_tocn(self)
  call soca_bkgerrgodas_socn(self)
  call soca_bkgerrgodas_ssh(self)

  ! Invent background error for ocnsfc fields: set it
  ! to 10% of the background for now ...
  do i=1,size(self%std_bkgerr%fields)
    field => self%std_bkgerr%fields(i)
    select case(field%name)
    case ('sw','lw','lhf','shf','us')
      call bkg%get(field%name, field_bkg)
      field%val = abs(field_bkg%val)
      field%val = 0.1_kind_real * field%val
    end select
  end do

  ! Apply config bounds to background error
  call self%bounds%apply(self%std_bkgerr)

  ! Save
  call self%std_bkgerr%write_file(fname)

end subroutine soca_bkgerrgodas_setup

! ------------------------------------------------------------------------------
!> Apply background error: dxm = D dxa
subroutine soca_bkgerrgodas_mult(self, dxa, dxm)
  type(soca_bkgerrgodas_config),  intent(in) :: self
  type(soca_increment),           intent(in) :: dxa
  type(soca_increment),        intent(inout) :: dxm

  type(soca_field), pointer :: field_m, field_e, field_a
  integer :: isc, iec, jsc, jec, i, j, n

  ! make sure fields are the right shape
  call dxa%check_congruent(dxm)
  call dxa%check_subset(self%std_bkgerr)

  ! Indices for compute domain (no halo)
  isc = self%bkg%geom%isc ; iec = self%bkg%geom%iec
  jsc = self%bkg%geom%jsc ; jec = self%bkg%geom%jec

  do n=1,size(dxa%fields)
    field_a => dxa%fields(n)
    call self%std_bkgerr%get(field_a%name, field_e)
    call dxm%get(field_a%name, field_m)
    do i = isc, iec
      do j = jsc, jec
        if (self%bkg%geom%mask2d(i,j).eq.1) then
          field_m%val(i,j,:) = field_e%val(i,j,:) * field_a%val(i,j,:)
        end if
      end do
    end do
  end do
end subroutine soca_bkgerrgodas_mult

! ------------------------------------------------------------------------------
!> Derive T background error from vertial gradient of temperature
subroutine soca_bkgerrgodas_tocn(self)
  type(soca_bkgerrgodas_config),     intent(inout) :: self

  real(kind=kind_real), allocatable :: sig1(:), sig2(:)
  type(soca_domain_indices) :: domain
  integer :: i, j, k
  integer :: iter, niter = 1
  type(soca_omb_stats) :: sst
  type(soca_field), pointer :: tocn_b, tocn_e, hocn, layer_depth

  ! Get compute domain indices
  domain%is = self%bkg%geom%isc ; domain%ie = self%bkg%geom%iec
  domain%js = self%bkg%geom%jsc ; domain%je = self%bkg%geom%jec

  ! Get local compute domain indices
  domain%isl = self%bkg%geom%iscl ; domain%iel = self%bkg%geom%iecl
  domain%jsl = self%bkg%geom%jscl ; domain%jel = self%bkg%geom%jecl

  ! Allocate temporary arrays
  allocate(sig1(self%bkg%geom%nzo), sig2(self%bkg%geom%nzo))

  ! Initialize sst background error to previously computed std of omb's
  ! Currently hard-coded to read GODAS file
  call sst%init(domain)
  call sst%bin(self%bkg%geom%lon, self%bkg%geom%lat)

  call self%bkg%get("tocn", tocn_b)
  call self%std_bkgerr%get("tocn", tocn_e)
  call self%bkg%get("hocn", hocn)
  call self%bkg%get("layer_depth",layer_depth)

  ! Loop over compute domain
  do i = domain%is, domain%ie
     do j = domain%js, domain%je
        if (self%bkg%geom%mask2d(i,j).eq.1) then

           ! Step 1: sigb from dT/dz
           call soca_diff(sig1(:), tocn_b%val(i,j,:), hocn%val(i,j,:))
           sig1(:) = self%t_dz * abs(sig1) ! Rescale background error

           ! Step 2: sigb based on efolding scale
           sig2(:) = self%bounds%t_min + (sst%bgerr_model(i,j)-self%bounds%t_min)*&
                &exp((layer_depth%val(i,j,1)-layer_depth%val(i,j,:))&
                  &/self%t_efold)

           ! Step 3: sigb = max(sig1, sig2)
           do k = 1, self%bkg%geom%nzo
              tocn_e%val(i,j,k) = min( max(sig1(k), sig2(k)), &
                & self%bounds%t_max)
           end do

           ! Step 4: Vertical smoothing
           do iter = 1, niter
              do k = 2, self%bkg%geom%nzo-1
                 tocn_e%val(i,j,k) = &
                      &( tocn_e%val(i,j,k-1)*hocn%val(i,j,k-1) +&
                      &  tocn_e%val(i,j,k)*hocn%val(i,j,k) +&
                      &  tocn_e%val(i,j,k+1)*hocn%val(i,j,k+1) )/&
                      & (sum(hocn%val(i,j,k-1:k+1)))
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
  integer :: i, j
  type(soca_field), pointer :: ssh

  ! Get compute domain indices
  domain%is = self%bkg%geom%isc ; domain%ie = self%bkg%geom%iec
  domain%js = self%bkg%geom%jsc ; domain%je = self%bkg%geom%jec

  call self%std_bkgerr%get("ssh", ssh)

  ! Loop over compute domain
  do i = domain%is, domain%ie
     do j = domain%js, domain%je
          if (self%bkg%geom%mask2d(i,j) .ne. 1) cycle

          if ( abs(self%bkg%geom%lat(i,j)) >= self%ssh_phi_ex) then
            ! if in extratropics, set to max value
            ssh%val(i,j,:) = self%bounds%ssh_max
          else
            ! otherwise, taper to min value (0.0) toward equator
            ssh%val(i,j,:) = self%bounds%ssh_min + 0.5 * &
              (self%bounds%ssh_max - self%bounds%ssh_min) * &
              (1 - cos(pi * self%bkg%geom%lat(i,j) / self%ssh_phi_ex))
          end if
     end do
  end do
end subroutine soca_bkgerrgodas_ssh


! ------------------------------------------------------------------------------
!> Derive unbalanced S background error, based on MLD
subroutine soca_bkgerrgodas_socn(self)
  type(soca_bkgerrgodas_config),     intent(inout) :: self
  !
  type(soca_domain_indices), target :: domain
  type(soca_field), pointer :: field, mld, layer_depth
  integer :: i, j, k
  real(kind=kind_real) :: r


  ! Get compute domain indices
  domain%is = self%bkg%geom%isc ; domain%ie = self%bkg%geom%iec
  domain%js = self%bkg%geom%jsc ; domain%je = self%bkg%geom%jec
  !
  ! Get local compute domain indices
  domain%isl = self%bkg%geom%iscl ; domain%iel = self%bkg%geom%iecl
  domain%jsl = self%bkg%geom%jscl ; domain%jel = self%bkg%geom%jecl

  ! TODO read in a precomputed surface S background error

  ! Loop over compute domain
  call self%std_bkgerr%get("socn", field)
  call self%bkg%get("mld", mld)
  call self%bkg%get("layer_depth", layer_depth)

  do i = domain%is, domain%ie
    do j = domain%js, domain%je
      if (self%bkg%geom%mask2d(i,j) /= 1)  cycle

      do k = 1, self%bkg%geom%nzo
        if ( layer_depth%val(i,j,k) <= mld%val(i,j,1)) then
          ! if in the mixed layer, set to the maximum value
          field%val(i,j,k) = self%bounds%s_max
        else
          ! otherwise, taper to the minium value below MLD
          r = 0.1 + 0.45 * (1-tanh( 2 * log( &
            & layer_depth%val(i,j,k) / mld%val(i,j,1) )))
          field%val(i,j,k) = max(self%bounds%s_min, r*self%bounds%s_max)
        end if
      end do
    end do
  end do
end subroutine soca_bkgerrgodas_socn

end module soca_bkgerrgodas_mod
