! (C) Copyright 2017-2019 UCAR.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


module soca_horizfilt_mod
  use config_mod
  use fckit_configuration_module, only: fckit_configuration
  use fckit_geometry_module, only: sphere_distance
  use iso_c_binding
  use kinds
  use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad
  use soca_fields_mod
  use soca_geom_mod_c
  use soca_geom_mod, only : soca_geom
  use soca_utils
  use tools_func, only: gc99
  use type_mpl, only: mpl_type
  use random_mod
  use variables_mod

  implicit none

  private
  public :: soca_horizfilt_setup, soca_horizfilt_delete
  public :: soca_horizfilt_mult, soca_horizfilt_multad

  !> Fortran derived type to hold configuration data for horizfilt
  type, public :: soca_horizfilt_type
     type(soca_field),         pointer :: bkg            !< Background field (or first guess)
     type(oops_vars)                   :: vars           !< Apply filtering to vars
     real(kind=kind_real), allocatable :: wgh(:,:,:,:)   !< Filtering weight
     character(len=:),     allocatable :: filter_type    !< Two choices:
                                                         !<   9 point weighted average: "wavg"
                                                         !<   flow dependent: "flow"
     real(kind=kind_real) :: l_ssh  !< Used with "flow" filter, sea surface height decorrelation scale
     integer :: isc, iec, jsc, jec  !< indices of compute domain
     integer :: isd, ied, jsd, jed  !< indices of data domain
   contains
     procedure :: setup => soca_horizfilt_setup
     procedure :: delete => soca_horizfilt_delete
     procedure :: mult => soca_horizfilt_mult
     procedure :: multad => soca_horizfilt_multad
  end type soca_horizfilt_type

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Setup for the horizfilt operator
  subroutine soca_horizfilt_setup(self, f_conf, geom, traj, vars)
    class(soca_horizfilt_type), intent(inout) :: self   !< The horizfilt structure
    type(fckit_configuration),     intent(in) :: f_conf !< The configuration
    type(soca_geom),               intent(in) :: geom   !< Geometry
    type(soca_field),              intent(in) :: traj   !< Trajectory
    type(oops_vars),               intent(in) :: vars   !< List of variables

    character(len=3)  :: domain
    integer :: isc, iec, jsc, jec, i, j, ivar, ii, jj
    logical :: init_seaice, init_ocean
    real(kind=kind_real) :: sum_dist, dist(-1:1,-1:1), sum_w, dist_tmp
    type(mpl_type) :: mpl
    character(len=4) :: filter_type

    ! Setup list of variables to apply filtering on
    self%vars = vars

    ! Get filter type from configuration
    call f_conf%get_or_die("filter_type", self%filter_type )

    if (self%filter_type=="flow") then
       call f_conf%get_or_die("l_ssh", self%l_ssh)
    endif

    ! Indices for compute/data domain
    self%isd = geom%isd ;  self%ied = geom%ied ; self%jsd = geom%jsd; self%jed = geom%jed
    self%isc = geom%isc ;  self%iec = geom%iec ; self%jsc = geom%jsc; self%jec = geom%jec

    ! Allocate and compute filtering weights
    allocate(self%wgh(self%isd:self%ied,self%jsd:self%jed,-1:1,-1:1))

    ! Compute distance based weights
    self%wgh = 0.0_kind_real
    do i = self%isc, self%iec
       do j = self%jsc, self%jec
          ! Great circle distance
          do ii = -1,1
             do jj = -1,1
                select case (trim(self%filter_type))
                case ("wavg")
                   ! Great circle distance
                   dist_tmp = sphere_distance(geom%lon(i,j), geom%lat(i,j), geom%lon(i+ii,j+jj), geom%lat(i+ii,j+jj) )
                case ("flow")
                   ! Following sea surface height
                   dist_tmp = gc99(mpl, abs((traj%ssh(i,j) - traj%ssh(i+ii,j+jj))/(self%l_ssh)))
                end select
                dist(ii,jj) = geom%mask2d(i+ii,j+jj) * dist_tmp
             end do
          end do

          ! Weighted to distance
          sum_dist=sum(dist)
          do ii = -1,1
             do jj = -1,1
                self%wgh(i,j,ii,jj) = geom%mask2d(i+ii,j+jj) * ( sum_dist - dist(ii,jj) )
             end do
          end do

          ! Normalize
          sum_w = sum(self%wgh(i,j,:,:))
          if (sum_w>0.0_kind_real) then
             do ii = -1,1
                do jj = -1,1
                   self%wgh(i,j,ii,jj) = self%wgh(i,j,ii,jj) / sum_w
                end do
             end do
          else
             self%wgh(i,j,ii,jj) = 0.0_kind_real
          end if
       end do
    end do

  end subroutine soca_horizfilt_setup

  ! ------------------------------------------------------------------------------
  !> Delete horizfilt
  subroutine soca_horizfilt_delete(self)
    class(soca_horizfilt_type), intent(inout) :: self       !< The horizfilt structure

    deallocate(self%wgh)
    nullify(self%bkg)

  end subroutine soca_horizfilt_delete

  ! ------------------------------------------------------------------------------
  !> Forward filtering
  subroutine soca_horizfilt_mult(self, dxin, dxout, geom)
    class(soca_horizfilt_type), intent(inout) :: self  !< The horizfilt structure
    type(soca_field),              intent(in) :: dxin  !< Input: Increment
    type(soca_field),           intent(inout) :: dxout !< Output: filtered Increment
    type(soca_geom),               intent(in) :: geom

    integer :: k, ivar, iter
    real(kind=kind_real), allocatable, dimension(:,:) :: dxi, dxo

    allocate(dxi(self%isd:self%ied,self%jsd:self%jed))
    allocate(dxo(self%isd:self%ied,self%jsd:self%jed))

    do ivar = 1, self%vars%nv
       select case (trim(self%vars%fldnames(ivar)))

       case ("ssh")
          dxi = dxin%ssh(:,:)
          call soca_filt2d(self, dxi, dxo, geom)
          dxout%ssh(:,:) = dxo

       case ("tocn")
          do k = 1, geom%nzo
             dxi = dxin%tocn(:,:,k)
             call soca_filt2d(self, dxi, dxo, geom)
             dxout%tocn(:,:,k) = dxo
          end do

       case ("socn")
          do k = 1, geom%nzo
             dxi = dxin%socn(:,:,k)
             call soca_filt2d(self, dxi, dxo, geom)
             dxout%socn(:,:,k) = dxo
          end do

       case ("cicen")
          do k = 1, geom%ncat
             dxi = dxin%seaice%cicen(:,:,k)
             call soca_filt2d(self, dxi, dxo, geom)
             dxout%seaice%cicen(:,:,k) = dxo
          end do

       case ("hicen")
          do k = 1, geom%ncat
             dxi = dxin%seaice%hicen(:,:,k)
             call soca_filt2d(self, dxi, dxo, geom)
             dxout%seaice%hicen(:,:,k) = dxo
          end do

       end select
    end do
    deallocate(dxi, dxo)

  end subroutine soca_horizfilt_mult

  ! ------------------------------------------------------------------------------
  !> Backward filtering
  subroutine soca_horizfilt_multad(self, dxin, dxout, geom)
    class(soca_horizfilt_type), intent(inout) :: self  !< The horizfilt structure
    type(soca_field),              intent(in) :: dxin  !< Input:
    type(soca_field),           intent(inout) :: dxout !< Output:
    type(soca_geom),               intent(in) :: geom

    integer :: k, ivar, iter
    real(kind=kind_real), allocatable, dimension(:,:) :: dxi, dxo

    allocate(dxi(self%isd:self%ied,self%jsd:self%jed))
    allocate(dxo(self%isd:self%ied,self%jsd:self%jed))

    do ivar = 1, self%vars%nv
       select case (trim(self%vars%fldnames(ivar)))

       case ("ssh")
          dxi = dxin%ssh(:,:)
          call soca_filt2d_ad(self, dxi, dxo, geom)
          dxout%ssh(:,:) = dxo

       case ("tocn")
          do k = 1, geom%nzo
             dxi = dxin%tocn(:,:,k)
             call soca_filt2d_ad(self, dxi, dxo, geom)
             dxout%tocn(:,:,k) = dxo
          end do

       case ("socn")
          do k = 1, geom%nzo
             dxi = dxin%socn(:,:,k)
             call soca_filt2d_ad(self, dxi, dxo, geom)
             dxout%socn(:,:,k) = dxo
          end do

       case ("cicen")
          do k = 1, geom%ncat
             dxi = dxin%seaice%cicen(:,:,k)
             call soca_filt2d_ad(self, dxi, dxo, geom)
             dxout%seaice%cicen(:,:,k) = dxo
          end do

       case ("hicen")
          do k = 1, geom%ncat
             dxi = dxin%seaice%hicen(:,:,k)
             call soca_filt2d_ad(self, dxi, dxo, geom)
             dxout%seaice%hicen(:,:,k) = dxo
          end do

       end select
    end do
    deallocate(dxi, dxo)

  end subroutine soca_horizfilt_multad

  ! ------------------------------------------------------------------------------
  !> Forward filtering for 2D array
  subroutine soca_filt2d(self, dxin, dxout, geom)
    class(soca_horizfilt_type),           intent(in) :: self
    real(kind=kind_real),    allocatable, intent(in) :: dxin(:,:)
    real(kind=kind_real), allocatable, intent(inout) :: dxout(:,:)
    type(soca_geom),                      intent(in) :: geom

    integer :: i, j
    real(kind=kind_real), allocatable :: dxtmp(:,:)

    ! Make a temporary copy of dxin
    allocate(dxtmp(self%isd:self%ied,self%jsd:self%jed))
    dxtmp = 0.0_kind_real
    dxtmp(self%isc:self%iec,self%jsc:self%jec) = dxin(self%isc:self%iec,self%jsc:self%jec)

    ! Update halo points of input array
    call mpp_update_domains(dxtmp, geom%Domain%mpp_domain, complete=.true.)

    ! 9-point distance weighted average
    dxout = 0.0_kind_real
    do i = self%isc, self%iec
       do j = self%jsc, self%jec
          if (geom%mask2d(i,j)==1) then
             dxout(i,j) = &
               self%wgh(i,j,-1,1)*dxtmp(i-1,j+1) + &
               self%wgh(i,j,0,1)*dxtmp(i,j+1) + &
               self%wgh(i,j,1,1)*dxtmp(i+1,j+1) + &
               self%wgh(i,j,-1,0)*dxtmp(i-1,j) + &
               self%wgh(i,j,0,0)*dxtmp(i,j) + &
               self%wgh(i,j,1,0)*dxtmp(i+1,j) + &
               self%wgh(i,j,-1,-1)*dxtmp(i-1,j-1) + &
               self%wgh(i,j,0,-1)*dxtmp(i,j-1) + &
               self%wgh(i,j,1,-1)*dxtmp(i+1,j-1)
          end if
       end do
    end do

    ! Update halo
    call mpp_update_domains(dxout, geom%Domain%mpp_domain, complete=.true.)

    deallocate(dxtmp)

  end subroutine soca_filt2d

  ! ------------------------------------------------------------------------------
  !> Backward filtering for 2D array
  subroutine soca_filt2d_ad(self, dxin, dxout, geom)
    class(soca_horizfilt_type),           intent(in) :: self
    real(kind=kind_real),    allocatable, intent(in) :: dxin(:,:)
    real(kind=kind_real), allocatable, intent(inout) :: dxout(:,:)
    type(soca_geom),                      intent(in) :: geom

    integer :: i, j
    real(kind=kind_real), allocatable :: dxtmp(:,:)

    ! Make a temporary copy of dxin
    allocate(dxtmp(self%isd:self%ied,self%jsd:self%jed))
    dxtmp = 0.0_kind_real
    dxtmp(self%isc:self%iec,self%jsc:self%jec) = dxin(self%isc:self%iec,self%jsc:self%jec)

    dxout = 0.0_kind_real
    ! Adjoint of 9-point weighted average
    do i =  self%iec, self%isc, -1
       do j = self%jec, self%jsc, -1
          if (geom%mask2d(i,j)==1) then
             dxout(i-1,j+1) = dxout(i-1,j+1) + self%wgh(i,j,-1, 1)*dxtmp(i,j)
             dxout(i,j+1)   = dxout(i,j+1)   + self%wgh(i,j, 0, 1)*dxtmp(i,j)
             dxout(i+1,j+1) = dxout(i+1,j+1) + self%wgh(i,j, 1, 1)*dxtmp(i,j)
             dxout(i-1,j)   = dxout(i-1,j)   + self%wgh(i,j,-1, 0)*dxtmp(i,j)
             dxout(i,j)     = dxout(i,j)     + self%wgh(i,j, 0, 0)*dxtmp(i,j)
             dxout(i+1,j)   = dxout(i+1,j)   + self%wgh(i,j, 1, 0)*dxtmp(i,j)
             dxout(i-1,j-1) = dxout(i-1,j-1) + self%wgh(i,j,-1,-1)*dxtmp(i,j)
             dxout(i,j-1)   = dxout(i,j-1)   + self%wgh(i,j, 0,-1)*dxtmp(i,j)
             dxout(i+1,j-1) = dxout(i+1,j-1) + self%wgh(i,j, 1,-1)*dxtmp(i,j)
          end if
       end do
    end do

    ! Adjoint of halo update
    call mpp_update_domains_ad(dxout, geom%Domain%mpp_domain, complete=.true.)

    deallocate(dxtmp)

  end subroutine soca_filt2d_ad

end module soca_horizfilt_mod
