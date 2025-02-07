! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_balance_mod

use fckit_configuration_module, only: fckit_configuration
use fms_mod, only: read_data
use kinds, only: kind_real
use atlas_module, only: atlas_field

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
  ! TODO the jacobians should really be stored in atlas fields, but
  !  I didn't feel like dealing with all that refactoring
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
  integer :: i, j, k, nl, idx
  real(kind=kind_real), allocatable :: jac(:), coef_mld, coef_layers

  type(atlas_field) :: tocn, socn, hocn, cice, mld, layer_depth
  real(kind=kind_real), pointer :: data_tocn(:,:), data_socn(:,:), data_hocn(:,:)
  real(kind=kind_real), pointer :: data_cice(:,:) => null(), data_mld(:,:), data_layer_depth(:,:)
  real(kind=kind_real), allocatable :: col_tocn(:), col_socn(:), col_hocn(:)

  ! declarations related to the dynamic height Jacobians
  character(len=:), allocatable :: filename

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
  tocn = traj%afieldset%field("sea_water_potential_temperature")
  socn = traj%afieldset%field("sea_water_salinity")
  hocn = traj%afieldset%field("sea_water_cell_thickness")
  mld = traj%afieldset%field("ocean_mixed_layer_thickness")
  layer_depth = traj%afieldset%field("sea_water_depth")
  call tocn%data(data_tocn)
  call socn%data(data_socn)
  call hocn%data(data_hocn)
  call mld%data(data_mld)
  call layer_depth%data(data_layer_depth)
  if (traj%has("sea_ice_area_fraction")) then
    cice = traj%afieldset%field("sea_ice_area_fraction")
    call cice%data(data_cice)
  end if

  ! allocate space
  nl = hocn%shape(1)
  allocate(self%kst%jacobian(isc:iec,jsc:jec,nl))
  self%kst%jacobian=0.0
  allocate(col_tocn(nl), col_socn(nl), col_hocn(nl))

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
          idx = geom%atlas_ij2idx(i,j)

          ! do nothing if on land
          if ( geom%mask2d(i, j) == 0 ) cycle

          ! compute dS(T)/dT
          do k=1,nl
             col_tocn(k) = data_tocn(k, idx)
             col_socn(k) = data_socn(k, idx)
             col_hocn(k) = data_hocn(k, idx)
          end do
          call soca_soft_jacobian(jac, col_tocn, col_socn, col_hocn, &
            self%kst%dsdtmax, self%kst%dsdzmin, self%kst%dtdzmin)

          ! filter out the Jacobian as specified in the configuration
          do k=1,nl
            coef_mld = soca_tanh_filt(data_layer_depth(k, idx), data_mld(1, idx))
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
      idx = geom%atlas_ij2idx(i,j)
      do k = 1, nl
        call soca_steric_jacobian (jac, &
          data_tocn(k, idx), data_socn(k, idx), data_layer_depth(k, idx), &
          data_hocn(k,idx), geom%lon(i,j), geom%lat(i,j))
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
      call read_data(filename, kct_name, kct, domain=geom%Domain%mpp_domain)
    end if
    allocate(self%kct(isc:iec,jsc:jec))
    self%kct = 0.0_kind_real
    do i = isc, iec
      do j = jsc, jec
        idx = geom%atlas_ij2idx(i,j)
        if (data_cice(1, idx) > 1.0e-3_kind_real) then
          self%kct = kct(i,j)
        end if
      end do
    end do
  end if

  ! Finalize fields
  call tocn%final()
  call socn%final()
  call hocn%final()
  call mld%final()
  call layer_depth%final()
  if (associated(data_cice)) call cice%final()

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

  type(atlas_field) :: fld_m, fld_a, tocn_a, socn_a
  real(kind=kind_real), pointer :: data_m(:,:), data_a(:,:), data_tocn(:,:), data_socn(:,:)
  integer :: i, j, k, n, idx

  !>    [ I       0   0  0 ]
  !>    [ Kst     I   0  0 ]
  !> K= [ Ketat Ketas I  0 ]
  !>    [ Kct     0   0  I ]

  tocn_a = dxa%afieldset%field("sea_water_potential_temperature")
  socn_a = dxa%afieldset%field("sea_water_salinity")
  call tocn_a%data(data_tocn)
  call socn_a%data(data_socn)

  do n=1, dxm%afieldset%size()
    fld_m = dxm%afieldset%field(n)
    fld_a = dxa%afieldset%field(n)
    call fld_m%data(data_m)
    call fld_a%data(data_a)

    select case(fld_m%name())
    case default
      do i = 1, fld_m%shape(2)
        do k = 1, fld_m%shape(1)
          data_m(k, i) = data_a(k, i)
        end do
      end do

    case ("sea_water_salinity")
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          do k = 1, fld_m%shape(1)
            data_m(k, idx) = data_a(k, idx) + &
              & self%kst%jacobian(i,j,k) * data_tocn(k, idx)
          end do
        end do
      end do
      call fld_m%set_dirty()

    case ("sea_surface_height_above_geoid")
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          data_m(1, idx) = data_a(1, idx)
          do k = 1, tocn_a%shape(1)
            data_m(1, idx) = data_m(1, idx) + &
              self%ksshts%kssht(i,j,k) * data_tocn(k, idx) +&
              self%ksshts%ksshs(i,j,k) * data_socn(k, idx)
          end do
        end do
      end do
      call fld_m%set_dirty()

    case ("sea_ice_area_fraction")
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          data_m(1, idx) = data_a(1, idx) + &
            self%kct(i,j) * data_tocn(1, idx)
        end do
      end do
      call fld_m%set_dirty()

    end select

    call fld_m%final()
    call fld_a%final()
  end do

  call tocn_a%final()
  call socn_a%final()

end subroutine soca_balance_mult


! ------------------------------------------------------------------------------
!> Apply backward balance operator
!!
!! \relates soca_balance_mod::soca_balance
subroutine soca_balance_multad(self, dxa, dxm)
  class(soca_balance),          intent(in)    :: self
  type(soca_increment), target, intent(in)    :: dxm !< input increment
  type(soca_increment), target, intent(inout) :: dxa !< output increment

  type(atlas_field) :: fld_a, fld_m, socn_m, ssh_m, cice_m
  real(kind=kind_real), pointer :: data_a(:,:), data_m(:,:), data_socn(:,:)
  real(kind=kind_real), pointer :: data_ssh(:,:), data_cice(:,:) => null()
  integer :: i, j, n, k, idx


  socn_m = dxm%afieldset%field("sea_water_salinity")
  ssh_m = dxm%afieldset%field("sea_surface_height_above_geoid")
  call socn_m%data(data_socn)
  call ssh_m%data(data_ssh)
  if (dxm%afieldset%has("sea_ice_area_fraction")) then
    cice_m = dxm%afieldset%field("sea_ice_area_fraction")
    call cice_m%data(data_cice)
  end if

  do n = 1, dxa%afieldset%size()
    fld_a = dxa%afieldset%field(n)
    fld_m = dxm%afieldset%field(n)
    call fld_a%data(data_a)
    call fld_m%data(data_m)

    select case(fld_a%name())
    case default
      do i = 1, fld_a%shape(2)
        do k = 1, fld_a%shape(1)
          data_a(k, i) = data_m(k, i)
        end do
      end do

    case ("sea_water_salinity")
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          do k = 1, fld_a%shape(1)
            data_a(k, idx) = data_m(k, idx) + &
              self%ksshts%ksshs(i,j,k) * data_ssh(1, idx)
          end do
        end do
      end do
      call fld_a%set_dirty()

    case ("sea_water_potential_temperature")
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          do k = 1, fld_a%shape(1)
            data_a(k, idx) = data_m(k, idx) + &
              self%kst%jacobian(i,j,k) * data_socn(k, idx) + &
              self%ksshts%kssht(i,j,k) * data_ssh(1, idx)
          end do
          if (associated(data_cice)) then
            data_a(1, idx) = data_a(1, idx) + &
              self%kct(i,j) * data_cice(1, idx)
          end if
        end do
      end do
      call fld_a%set_dirty()
    end select

    call fld_a%final()
    call fld_m%final()
  end do

  call socn_m%final()
  call ssh_m%final()
  if ( associated(data_cice) ) call cice_m%final()

end subroutine soca_balance_multad


! ------------------------------------------------------------------------------
!> Apply inverse of the forward balance operator
!!
!! \relates soca_balance_mod::soca_balance
subroutine soca_balance_multinv(self, dxa, dxm)
  class(soca_balance),          intent(in)    :: self
  type(soca_increment), target, intent(in)    :: dxm !< input increment
  type(soca_increment), target, intent(inout) :: dxa !< output increment

  integer :: i, j, k, n, idx

  type(atlas_Field) :: fld_m, fld_a, tocn_m, socn_m
  real(kind=kind_real), pointer :: data_m(:,:), data_a(:,:), data_tocn(:,:), data_socn(:,:)

  tocn_m = dxm%afieldset%field("sea_water_potential_temperature")
  socn_m = dxm%afieldset%field("sea_water_salinity")
  call tocn_m%data(data_tocn)
  call socn_m%data(data_socn)

  do n = 1, dxa%afieldset%size()
    fld_m = dxm%afieldset%field(n)
    fld_a = dxa%afieldset%field(n)
    call fld_m%data(data_m)
    call fld_a%data(data_a)

    select case(fld_a%name())
    case default
      do i = 1, fld_a%shape(2)
        do k = 1, fld_a%shape(1)
          data_a(k, i) = data_m(k, i)
        end do
      end do

    case ("sea_water_salinity")
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          do k = 1, fld_a%shape(1)
            data_a(k, idx) = data_m(k, idx) - &
              self%kst%jacobian(i,j,k) * data_tocn(k, idx)
          end do
        end do
      end do
      call fld_a%set_dirty()

    case ("sea_surface_height_above_geoid")
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          data_a(1, idx) = data_m(1, idx)
          do k = 1, tocn_m%shape(1)
            data_a(1, idx) = data_a(1, idx) + &
              ( self%ksshts%ksshs(i,j,k) * self%kst%jacobian(i,j,k) - &
              self%ksshts%kssht(i,j,k) ) *  data_tocn(k, idx) - &
              self%ksshts%ksshs(i,j,k) * data_socn(k, idx)
          end do
        end do
      end do
      call fld_a%set_dirty()

    case ("sea_ice_area_fraction")
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          data_a(1, idx) = data_m(1, idx) - &
            self%kct(i,j) * data_tocn(1, idx)
        end do
      end do
      call fld_a%set_dirty()

    end select
    call fld_m%final()
    call fld_a%final()
  end do
  call tocn_m%final()
  call socn_m%final()

end subroutine soca_balance_multinv


! ------------------------------------------------------------------------------
!> Apply inverse of the backward balance operator
!!
!! \relates soca_balance_mod::soca_balance
subroutine soca_balance_multinvad(self, dxa, dxm)
  class(soca_balance),          intent(in)    :: self
  type(soca_increment), target, intent(inout) :: dxm !< output increment
  type(soca_increment), target, intent(in)    :: dxa !< input increment

  integer :: i, j, k, n, idx

  type(atlas_field) :: fld_a, fld_m, socn_a, ssh_a, cice_a
  real(kind=kind_real), pointer :: data_a(:,:), data_m(:,:), data_socn(:,:)
  real(kind=kind_real), pointer :: data_ssh(:,:), data_cice(:,:) => null()

  socn_a = dxa%afieldset%field("sea_water_salinity")
  ssh_a = dxa%afieldset%field("sea_surface_height_above_geoid")
  call socn_a%data(data_socn)
  call ssh_a%data(data_ssh)
  if (dxa%afieldset%has("sea_ice_area_fraction")) then
    cice_a = dxa%afieldset%field("sea_ice_area_fraction")
    call cice_a%data(data_cice)
  end if

  do n = 1, dxm%afieldset%size()
    fld_m = dxm%afieldset%field(n)
    fld_a = dxa%afieldset%field(n)
    call fld_m%data(data_m)
    call fld_a%data(data_a)

    select case(fld_m%name())
    case default
      do i = 1, fld_m%shape(2)
        do k = 1, fld_m%shape(1)
          data_m(k, i) = data_a(k, i)
        end do
      end do

    case ("sea_water_potential_temperature")
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          do k = 1, fld_m%shape(1)
            data_m(k, idx) = data_a(k, idx) - &
              self%kst%jacobian(i,j,k) * data_socn(k, idx) + &
              ( self%ksshts%ksshs(i,j,k) * self%kst%jacobian(i,j,k) - &
              self%ksshts%kssht(i,j,k) ) * data_ssh(1, idx)
          end do
          if (associated(data_cice)) then
            data_m(1, idx) = data_m(1, idx) - &
              self%kct(i,j) * data_cice(1, idx)
          end if
        end do
      end do
      call fld_m%set_dirty()

    case ("sea_water_salinity")
      do j = self%geom%jsc, self%geom%jec
        do i = self%geom%isc, self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          do k = 1, fld_m%shape(1)
            data_m(k, idx) = data_a(k, idx) - &
              self%ksshts%ksshs(i,j,k) * data_ssh(1, idx)
          end do
        end do
      end do
      call fld_m%set_dirty()

    end select
    call fld_m%final()
    call fld_a%final()
  end do

  call socn_a%final()
  call ssh_a%final()
  if (associated(data_cice)) call cice_a%final()

end subroutine soca_balance_multinvad

end module soca_balance_mod
