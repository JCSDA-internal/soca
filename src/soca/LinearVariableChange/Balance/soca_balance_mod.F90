! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_balance_mod

use fckit_configuration_module, only: fckit_configuration
use fms_io_mod, only: fms_io_init, fms_io_exit
use fms_mod, only: read_data
use kinds, only: kind_real

! soca modules
use soca_fields_mod, only: soca_field
use soca_geom_mod, only: soca_geom
use soca_increment_mod, only: soca_increment
use soca_ksshts_mod, only: soca_ksshts, soca_steric_jacobian
use soca_kst_mod, only: soca_kst, soca_soft_jacobian
use soca_state_mod, only: soca_state

implicit none
private


!> Variable transform for the balance operators (K)
!!
!! The core of the balance transformations are provided by
!! soca_ksshts_mod::soca_ksshts and soca_kst_mod::soca_kst
type, public :: soca_balance
  ! private members
  type(soca_kst), private             :: kst                 !< T/S balance
  type(soca_ksshts), private          :: ksshts              !< SSH/T/S balance
  real(kind=kind_real), private, allocatable :: kct(:,:)     !< C/T Jacobian
  type(soca_geom),  pointer, private       :: geom !< geometry

contains
  !> \copybrief soca_balance_setup \see soca_balance_setup
  procedure :: setup => soca_balance_setup

  !> \copybrief soca_balance_delete \see soca_balance_delete
  procedure :: delete => soca_balance_delete

  !> \copybrief soca_balance_mult \see soca_balance_mult
  procedure :: mult => soca_balance_mult

  !> \copybrief soca_balance_multad \see soca_balance_multad
  procedure :: multad => soca_balance_multad

  !> \copybrief soca_balance_multinv \see soca_balance_multinv
  procedure :: multinv => soca_balance_multinv

  !> \copybrief soca_balance_multinvad \see soca_balance_multinvad
  procedure :: multinvad => soca_balance_multinvad

end type soca_balance


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


function soca_tanh_filt(l, l0) result (coef)
  real(kind=kind_real), intent(in) :: l
  real(kind=kind_real), intent(in) :: l0

  real(kind=kind_real) :: coef

  coef = 0.5_kind_real*(tanh(l-l0)+1.0_kind_real)

end function soca_tanh_filt

! ------------------------------------------------------------------------------
!> Initialization of the balance operator and its trajectory.
!!
!! - balances always used: T,S,SSH
!! - optional balances depending on input fields: cice
!! \relates soca_balance_mod::soca_balance
subroutine soca_balance_setup(self, f_conf, traj, geom)
  class(soca_balance),       intent(inout) :: self
  type(fckit_configuration),   intent(in)  :: f_conf !< configuration
  type(soca_state),    target, intent(in)  :: traj !< trajectory
  type(soca_geom),     target, intent(in)  :: geom !< geometry

  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: i, j, k, nl
  real(kind=kind_real), allocatable :: jac(:), coef_mld, coef_layers
  type(soca_field), pointer :: tocn, socn, hocn, cice, mld, layer_depth

  ! declarations related to the dynamic height Jacobians
  character(len=:), allocatable :: filename
  real(kind=kind_real) :: threshold

  ! declarations related to the sea-ice Jacobian
  character(len=:), allocatable :: kct_name
  real(kind=kind_real), allocatable :: kct(:,:) !> dc/dT

  self%geom => geom

  ! Indices for compute domain
  isc=geom%isc; iec=geom%iec
  jsc=geom%jsc; jec=geom%jec
  isd=geom%isd; ied=geom%ied
  jsd=geom%jsd; jed=geom%jed

  ! Get required fields
  call traj%get("sea_water_potential_temperature", tocn)
  call traj%get("sea_water_salinity", socn)
  call traj%get("sea_water_cell_thickness", hocn)
  call traj%get("ocean_mixed_layer_thickness", mld)
  call traj%get("sea_water_depth", layer_depth)
  if (traj%has("sea_ice_area_fraction"))  call traj%get("sea_ice_area_fraction", cice)

  ! allocate space
  nl = hocn%nz
  allocate(self%kst%jacobian(isc:iec,jsc:jec,geom%nzo))
  self%kst%jacobian=0.0

  ! Setup Kst if in the configuration
  if ( f_conf%has("kst") ) then
     allocate(jac(nl))
     call f_conf%get_or_die("kst.dsdtmax", self%kst%dsdtmax)
     call f_conf%get_or_die("kst.dsdzmin", self%kst%dsdzmin)
     call f_conf%get_or_die("kst.dtdzmin", self%kst%dtdzmin)
     call f_conf%get_or_die("kst.nlayers", self%kst%nlayers)

     ! Compute and store Jacobian of Kst
     do i = isc, iec
        do j = jsc, jec
           ! do nothing if on land
           if ( geom%mask2d(i, j) == 0 ) cycle

           ! compute dS(T)/dT
           call soca_soft_jacobian(jac,&
                &tocn%val(i,j,:),&
                &socn%val(i,j,:),&
                &hocn%val(i,j,:),&
                &self%kst%dsdtmax, self%kst%dsdzmin, self%kst%dtdzmin)

           ! filter out the Jacobian as specified in the configuration
           do k=1,nl
              coef_mld = soca_tanh_filt(layer_depth%val(i,j,k),mld%val(i,j,1))
              coef_layers = soca_tanh_filt(real(k, kind=kind_real), real(self%kst%nlayers, kind=kind_real))
              self%kst%jacobian(i,j,k) = jac(k)*coef_mld*coef_layers
           end do
        end do
     end do
     deallocate(jac)
  end if

  ! Get configuration for Ksshts
  self%ksshts%nlayers = -999   ! No filtering by default
  if ( f_conf%has("ksshts") ) call f_conf%get_or_die("ksshts.nlayers", self%ksshts%nlayers)

  ! Compute Jacobian of Ksshts
  allocate(self%ksshts%kssht, mold=self%kst%jacobian)
  allocate(self%ksshts%ksshs, mold=self%kst%jacobian)
  allocate(jac(2))
  self%ksshts%kssht=0.0_kind_real
  self%ksshts%ksshs=0.0_kind_real
  do i = isc, iec
    do j = jsc, jec
      if (geom%mask2d(i,j) == 0.0) cycle
      do k = 1, nl
        call soca_steric_jacobian (jac, &
        tocn%val(i,j,k), &
        socn%val(i,j,k), &
        &layer_depth%val(i,j,k),&
        &hocn%val(i,j,k),&
        &geom%lon(i,j),&
        &geom%lat(i,j))
        coef_layers = soca_tanh_filt(real(k, kind=kind_real), real(self%ksshts%nlayers, kind=kind_real))
        self%ksshts%kssht(i,j,k) = jac(1)*coef_layers
        self%ksshts%ksshs(i,j,k) = jac(2)*coef_layers
     end do
    end do
  end do
  deallocate(jac)

  ! Compute Kct
  if (traj%has("sea_ice_area_fraction")) then
    ! Setup dc/dT
    allocate(kct(isd:ied,jsd:jed))
    kct = 0.0_kind_real
    if ( f_conf%has("dcdt") ) then
      call f_conf%get_or_die("dcdt.filename", filename)
      call f_conf%get_or_die("dcdt.name", kct_name)
      call fms_io_init()
      call read_data(filename, kct_name, kct, domain=geom%Domain%mpp_domain)
      call fms_io_exit()
    end if
    allocate(self%kct(isc:iec,jsc:jec))
    self%kct = 0.0_kind_real
    do i = isc, iec
      do j = jsc, jec
          if (cice%val(i,j,1) > 1.0e-3_kind_real) then
            self%kct = kct(i,j)
          end if
      end do
    end do
  end if

end subroutine soca_balance_setup


! ------------------------------------------------------------------------------
!> Destructor for the balance oprator
!!
!! \relates soca_balance_mod::soca_balance
subroutine soca_balance_delete(self)
  class(soca_balance), intent(inout) :: self

  ! the following always exist
  deallocate(self%kst%jacobian)
  deallocate(self%ksshts%kssht)
  deallocate(self%ksshts%ksshs)

  ! only exists if cice was given
  if (allocated(self%kct)) deallocate(self%kct)
end subroutine soca_balance_delete


! ------------------------------------------------------------------------------
!> Apply forward balance operator
!!
!! \relates soca_balance_mod::soca_balance
subroutine soca_balance_mult(self, dxa, dxm)
  class(soca_balance),          intent(in)    :: self
  type(soca_increment), target, intent(in)    :: dxa !< input increment
  type(soca_increment), target, intent(inout) :: dxm !< output increment

  type(soca_field), pointer :: fld_m, fld_a
  type(soca_field), pointer :: tocn_a, socn_a

  integer :: i, j, k, n

  !>    [ I       0   0  0 ]
  !>    [ Kst     I   0  0 ]
  !> K= [ Ketat Ketas I  0 ]
  !>    [ Kct     0   0  I ]

  call dxa%get("sea_water_potential_temperature",tocn_a)
  call dxa%get("sea_water_salinity",socn_a)

  do n=1, size(dxm%fields)
    fld_m => dxm%fields(n)
    fld_a => dxa%fields(n)

    do i = self%geom%isc, self%geom%iec
      do j = self%geom%jsc, self%geom%jec
        select case(fld_m%name)
        case default
          fld_m%val(i,j,:) = fld_a%val(i,j,:)

        case("sea_water_salinity") ! Salinity
          fld_m%val(i,j,:) = fld_a%val(i,j,:) + &
            & self%kst%jacobian(i,j,:) * tocn_a%val(i,j,:)

        case ("sea_surface_height_above_geoid") ! SSH
          fld_m%val(i,j,:) = fld_a%val(i,j,:)
          do k = 1, tocn_a%nz
            fld_m%val(i,j,:) = fld_m%val(i,j,:) + &
              & self%ksshts%kssht(i,j,k) * tocn_a%val(i,j,k) + &
              & self%ksshts%ksshs(i,j,k) * socn_a%val(i,j,k)
          end do

        case ("sea_ice_area_fraction") ! Ice fraction
            fld_m%val(i,j,1) = fld_a%val(i,j,1) + &
              & self%kct(i,j) * tocn_a%val(i,j,1)

        end select
      end do
    end do
  end do
end subroutine soca_balance_mult


! ------------------------------------------------------------------------------
!> Apply backward balance operator
!!
!! \relates soca_balance_mod::soca_balance
subroutine soca_balance_multad(self, dxa, dxm)
  class(soca_balance),          intent(in)    :: self
  type(soca_increment), target, intent(in)    :: dxm !< input increment
  type(soca_increment), target, intent(inout) :: dxa !< output increment

  type(soca_field), pointer :: fld_a, fld_m
  type(soca_field), pointer :: socn_m, ssh_m, cice_m
  integer :: i, j, n

  cice_m => null()

  call dxm%get("sea_water_salinity", socn_m)
  call dxm%get("sea_surface_height_above_geoid",  ssh_m)
  if (dxm%has("sea_ice_area_fraction")) call dxm%get("sea_ice_area_fraction",cice_m)

  do n = 1, size(dxa%fields)
    fld_a => dxa%fields(n)
    fld_m => dxm%fields(n)

    do i = self%geom%isc, self%geom%iec
      do j = self%geom%jsc, self%geom%jec
        select case(fld_a%name)
        case default
          fld_a%val(i,j,:) = fld_m%val(i,j,:)

        case ("sea_water_potential_temperature") ! Temperature
          fld_a%val(i,j,:) = fld_m%val(i,j,:) + &
            & self%kst%jacobian(i,j,:) * socn_m%val(i,j,:) + &
            & self%ksshts%kssht(i,j,:) * ssh_m%val(i,j,1)

          if (associated(cice_m)) then ! use cice only if present
            fld_a%val(i,j,1) = fld_a%val(i,j,1) + &
              & self%kct(i,j) * cice_m%val(i,j,1)
          end if

        case ("sea_water_salinity") ! Salinity
          fld_a%val(i,j,:) = fld_m%val(i,j,:) + &
            & self%ksshts%ksshs(i,j,:) * ssh_m%val(i,j, 1)

        end select
      end do
    end do
  end do
end subroutine soca_balance_multad


! ------------------------------------------------------------------------------
!> Apply inverse of the forward balance operator
!!
!! \relates soca_balance_mod::soca_balance
subroutine soca_balance_multinv(self, dxa, dxm)
  class(soca_balance),          intent(in)    :: self
  type(soca_increment), target, intent(in)    :: dxm !< input increment
  type(soca_increment), target, intent(inout) :: dxa !< output increment

  integer :: i, j, k, n
  type(soca_field), pointer :: fld_m, fld_a
  type(soca_field), pointer :: tocn_m, socn_m

  call dxm%get("sea_water_potential_temperature", tocn_m)
  call dxm%get("sea_water_salinity", socn_m)

  do n = 1, size(dxa%fields)
    fld_a => dxa%fields(n)
    fld_m => dxm%fields(n)

    do i = self%geom%isc, self%geom%iec
      do j = self%geom%jsc, self%geom%jec
        select case(fld_a%name)
        case default
          fld_a%val(i,j,:) = fld_m%val(i,j,:)

        case ('sea_water_salinity') ! Salinity
          fld_a%val(i,j,:) = fld_m%val(i,j,:) - &
            & self%kst%jacobian(i,j,:) * tocn_m%val(i,j,:)

        case ('sea_surface_height_above_geoid') ! SSH
          fld_a%val(i,j, :) = fld_m%val(i,j, :)
          do k = 1, tocn_m%nz
            fld_a%val(i,j,:) = fld_a%val(i,j,:) + &
              & ( self%ksshts%ksshs(i,j,k) * self%kst%jacobian(i,j,k) - &
              & self%ksshts%kssht(i,j,k) ) *  tocn_m%val(i,j,k) - &
              & self%ksshts%ksshs(i,j,k) * socn_m%val(i,j,k)
          end do

        case ('sea_ice_area_fraction') ! Ice fraction
          fld_a%val(i,j,:) =  fld_m%val(i,j,:)
            fld_a%val(i,j,1) = fld_a%val(i,j,1) - &
              & self%kct(i,j) * tocn_m%val(i,j,1)

        end select
      end do
    end do
  end do
end subroutine soca_balance_multinv


! ------------------------------------------------------------------------------
!> Apply inverse of the backward balance operator
!!
!! \relates soca_balance_mod::soca_balance
subroutine soca_balance_multinvad(self, dxa, dxm)
  class(soca_balance),          intent(in)    :: self
  type(soca_increment), target, intent(inout) :: dxm !< output increment
  type(soca_increment), target, intent(in)    :: dxa !< input increment

  integer :: i, j, n
  type(soca_field), pointer :: fld_a, fld_m
  type(soca_field), pointer :: socn_a, ssh_a, cice_a

  cice_a => null()

  call dxa%get("sea_water_salinity", socn_a)
  call dxa%get("sea_surface_height_above_geoid",  ssh_a)
  if (dxa%has("sea_ice_area_fraction")) call dxa%get("sea_ice_area_fraction",cice_a)

  do n = 1, size(dxm%fields)
    fld_m => dxm%fields(n)
    fld_a => dxa%fields(n)

    do i = self%geom%isc, self%geom%iec
      do j = self%geom%jsc, self%geom%jec
        select case (fld_m%name)
        case default
          fld_m%val(i,j,:) = fld_a%val(i,j,:)

        case ('sea_water_potential_temperature') ! Temperature
          fld_m%val(i,j,:) = fld_a%val(i,j,:) &
            & - self%kst%jacobian(i,j,:) * socn_a%val(i,j,:) &
            & + ( self%ksshts%ksshs(i,j,:) * self%kst%jacobian(i,j,:) &
            &     - self%ksshts%kssht(i,j,:) ) * ssh_a%val(i,j,1)

          if (associated(cice_a)) then ! use cice only if present
            fld_m%val(i,j,1) = fld_m%val(i,j,1) &
              & - self%kct(i,j) * cice_a%val(i,j,1)
          end if

        case ('sea_water_salinity') ! Salinity
          fld_m%val(i,j,:) = fld_a%val(i,j,:) - &
            & self%ksshts%ksshs(i,j,:) * ssh_a%val(i,j,1)

        end select
      end do
    end do
  end do
end subroutine soca_balance_multinvad

end module soca_balance_mod
