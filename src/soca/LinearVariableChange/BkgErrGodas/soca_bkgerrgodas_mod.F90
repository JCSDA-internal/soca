! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> variable transform: background error
module soca_bkgerrgodas_mod

use fckit_configuration_module, only: fckit_configuration
use tools_const, only : pi
use kinds, only: kind_real

! soca modules
use soca_bkgerrutil_mod, only: soca_bkgerr_bounds_type
use soca_fields_mod, only: soca_field, soca_fields
use soca_geom_mod, only: soca_geom
use soca_increment_mod, only: soca_increment
use soca_omb_stats_mod, only: soca_omb_stats, soca_domain_indices
use soca_state_mod, only: soca_state
use soca_utils, only: soca_diff

implicit none
private


!> Variable transform for background error (D), GODAS version
type, public :: soca_bkgerrgodas
  type(soca_state),         pointer :: bkg
  type(soca_fields)                 :: std_bkgerr

  ! private members
  type(soca_geom), pointer, private :: geom
  type(soca_bkgerr_bounds_type)     :: bounds         !< Bounds for bkgerrgodas
  real(kind=kind_real), private     :: t_dz           !< For rescaling of the vertical gradient
  real(kind=kind_real), private     :: t_efold        !< E-folding scale for surf based T min
  real(kind=kind_real), private     :: ssh_phi_ex     !< latitude scale of ssh error

contains

  !> \copybrief soca_bkgerrgodas_setup \see soca_bkgerrgodas_setup
  procedure :: setup => soca_bkgerrgodas_setup

  !> \copybrief soca_bkgerrgodas_mult \see soca_bkgerrgodas_mult
  procedure :: mult => soca_bkgerrgodas_mult

end type soca_bkgerrgodas


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Setup the static background error
!!
!! \relates soca_bkgerrgodas_mod::soca_bkgerrgodas
subroutine soca_bkgerrgodas_setup(self, f_conf, bkg, geom)
  class(soca_bkgerrgodas), target,  intent(inout) :: self
  type(fckit_configuration),        intent(in)    :: f_conf !< configuration
  type(soca_state),         target, intent(in)    :: bkg !< background state
  type(soca_geom),          target, intent(in)    :: geom !< model geometry

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

  ! Associate background and geometry
  self%bkg => bkg
  self%geom => geom

  ! Set all fields to zero
  call self%std_bkgerr%zeros()

  ! Std of bkg error for T/S/SSH based on background.
  ! S and SSH error are only for the unbalanced portion of S/SSH
  call soca_bkgerrgodas_tocn(self)
  call soca_bkgerrgodas_socn(self)
  call soca_bkgerrgodas_ssh(self)

  ! Invent background error for ocnsfc, wav and ocn_bgc fields: set
  ! it to 10% or 20% of the background for now ...
  do i=1,size(self%std_bkgerr%fields)
    field => self%std_bkgerr%fields(i)
    select case(field%name)
    case ('sw','lw','lhf','shf','us','swh')
      call bkg%get(field%name, field_bkg)
      field%val = abs(field_bkg%val)
      field%val = 0.1_kind_real * field%val
    case ('chl','biop')
      call bkg%get(field%name, field_bkg)
      field%val = abs(field_bkg%val) * 0.2_kind_real
    end select
  end do

  ! Apply config bounds to background error
  call self%bounds%apply(self%std_bkgerr)

  ! Save
  call self%std_bkgerr%write_file(fname)

end subroutine soca_bkgerrgodas_setup


! ------------------------------------------------------------------------------
!> Apply background error: dxm = D dxa
!!
!! \relates soca_bkgerrgodas_mod::soca_bkgerrgodas
subroutine soca_bkgerrgodas_mult(self, dxa, dxm)
  class(soca_bkgerrgodas),      intent(in)    :: self
  type(soca_increment), target, intent(in)    :: dxa !< input increment
  type(soca_increment), target, intent(inout) :: dxm !< output increment

  type(soca_field), pointer :: field_m, field_e, field_a
  integer :: isc, iec, jsc, jec, i, j, n

  ! make sure fields are the right shape
  call dxa%check_congruent(dxm)
  call dxa%check_subset(self%std_bkgerr)

  ! Indices for compute domain (no halo)
  isc = self%geom%isc ; iec = self%geom%iec
  jsc = self%geom%jsc ; jec = self%geom%jec

  do n=1,size(dxa%fields)
    field_a => dxa%fields(n)
    call self%std_bkgerr%get(field_a%name, field_e)
    call dxm%get(field_a%name, field_m)
    do i = isc, iec
      do j = jsc, jec
        if (self%geom%mask2d(i,j).eq.1) then
          field_m%val(i,j,:) = field_e%val(i,j,:) * field_a%val(i,j,:)
        end if
      end do
    end do
  end do
end subroutine soca_bkgerrgodas_mult


! ------------------------------------------------------------------------------
!> Derive T background error from vertial gradient of temperature
!!
!! \relates soca_bkgerrgodas_mod::soca_bkgerrgodas
subroutine soca_bkgerrgodas_tocn(self)
  class(soca_bkgerrgodas),     intent(inout) :: self

  real(kind=kind_real), allocatable :: sig1(:), sig2(:)
  type(soca_domain_indices) :: domain
  integer :: i, j, k
  integer :: iter, niter = 1
  type(soca_omb_stats) :: sst
  type(soca_field), pointer :: tocn_b, tocn_e, hocn, layer_depth

  ! Get compute domain indices
  domain%is = self%geom%isc ; domain%ie = self%geom%iec
  domain%js = self%geom%jsc ; domain%je = self%geom%jec

  ! Get local compute domain indices
  domain%isl = self%geom%iscl ; domain%iel = self%geom%iecl
  domain%jsl = self%geom%jscl ; domain%jel = self%geom%jecl

  ! Allocate temporary arrays
  allocate(sig1(self%geom%nzo), sig2(self%geom%nzo))

  ! Initialize sst background error to previously computed std of omb's
  ! Currently hard-coded to read GODAS file
  call sst%init(domain)
  call sst%bin(self%geom%lon, self%geom%lat)

  call self%bkg%get("tocn", tocn_b)
  call self%std_bkgerr%get("tocn", tocn_e)
  call self%bkg%get("hocn", hocn)
  call self%bkg%get("layer_depth",layer_depth)

  ! Loop over compute domain
  do i = domain%is, domain%ie
     do j = domain%js, domain%je
        if (self%geom%mask2d(i,j).eq.1) then

           ! Step 1: sigb from dT/dz
           call soca_diff(sig1(:), tocn_b%val(i,j,:), hocn%val(i,j,:))
           sig1(:) = self%t_dz * abs(sig1) ! Rescale background error

           ! Step 2: sigb based on efolding scale
           sig2(:) = self%bounds%t_min + (sst%bgerr_model(i,j)-self%bounds%t_min)*&
                &exp((layer_depth%val(i,j,1)-layer_depth%val(i,j,:))&
                  &/self%t_efold)

           ! Step 3: sigb = max(sig1, sig2)
           do k = 1, self%geom%nzo
              tocn_e%val(i,j,k) = min( max(sig1(k), sig2(k)), &
                & self%bounds%t_max)
           end do

           ! Step 4: Vertical smoothing
           do iter = 1, niter
              do k = 2, self%geom%nzo-1
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
!!
!! \relates soca_bkgerrgodas_mod::soca_bkgerrgodas
subroutine soca_bkgerrgodas_ssh(self)
  class(soca_bkgerrgodas),     intent(inout) :: self
  type(soca_domain_indices), target :: domain
  integer :: i, j
  type(soca_field), pointer :: ssh

  ! Get compute domain indices
  domain%is = self%geom%isc ; domain%ie = self%geom%iec
  domain%js = self%geom%jsc ; domain%je = self%geom%jec

  call self%std_bkgerr%get("ssh", ssh)

  ! Loop over compute domain
  do i = domain%is, domain%ie
     do j = domain%js, domain%je
          if (self%geom%mask2d(i,j) .ne. 1) cycle

          if ( abs(self%geom%lat(i,j)) >= self%ssh_phi_ex) then
            ! if in extratropics, set to max value
            ssh%val(i,j,:) = self%bounds%ssh_max
          else
            ! otherwise, taper to min value (0.0) toward equator
            ssh%val(i,j,:) = self%bounds%ssh_min + 0.5 * &
              (self%bounds%ssh_max - self%bounds%ssh_min) * &
              (1 - cos(pi * self%geom%lat(i,j) / self%ssh_phi_ex))
          end if
     end do
  end do
end subroutine soca_bkgerrgodas_ssh


! ------------------------------------------------------------------------------
!> Derive unbalanced S background error, based on MLD
!!
!! \relates soca_bkgerrgodas_mod::soca_bkgerrgodas
subroutine soca_bkgerrgodas_socn(self)
  class(soca_bkgerrgodas),     intent(inout) :: self
  !
  type(soca_domain_indices), target :: domain
  type(soca_field), pointer :: field, mld, layer_depth
  integer :: i, j, k
  real(kind=kind_real) :: r


  ! Get compute domain indices
  domain%is = self%geom%isc ; domain%ie = self%geom%iec
  domain%js = self%geom%jsc ; domain%je = self%geom%jec
  !
  ! Get local compute domain indices
  domain%isl = self%geom%iscl ; domain%iel = self%geom%iecl
  domain%jsl = self%geom%jscl ; domain%jel = self%geom%jecl

  ! TODO read in a precomputed surface S background error

  ! Loop over compute domain
  call self%std_bkgerr%get("socn", field)
  call self%bkg%get("mld", mld)
  call self%bkg%get("layer_depth", layer_depth)

  do i = domain%is, domain%ie
    do j = domain%js, domain%je
      if (self%geom%mask2d(i,j) /= 1)  cycle

      do k = 1, self%geom%nzo
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
