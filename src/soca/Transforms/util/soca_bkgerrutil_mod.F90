! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerrutil_mod
  use fckit_configuration_module, only: fckit_configuration
  use kinds, only: kind_real
  use soca_fields_mod, only: soca_field
  use soca_utils

  implicit none
  private

  type, public :: soca_bkgerr_bounds_type
     real(kind=kind_real) :: t_min, t_max
     real(kind=kind_real) :: s_min, s_max
     real(kind=kind_real) :: ssh_min, ssh_max
     real(kind=kind_real) :: cicen_min, cicen_max
     real(kind=kind_real) :: hicen_min, hicen_max
   contains
     procedure :: read => soca_bkgerr_readbounds
     procedure :: apply => soca_bkgerr_applybounds
  end type soca_bkgerr_bounds_type

contains

  ! ------------------------------------------------------------------------------
  !> Read bounds from config
  subroutine soca_bkgerr_readbounds(self, f_conf)
    class(soca_bkgerr_bounds_type), intent(inout) :: self
    type(fckit_configuration),      intent(in)    :: f_conf

    ! Get bounds from configuration
    call f_conf%get_or_die("t_min", self%t_min)
    call f_conf%get_or_die("t_max", self%t_max)
    call f_conf%get_or_die("s_min", self%s_min)
    call f_conf%get_or_die("s_max", self%s_max)
    call f_conf%get_or_die("ssh_min", self%ssh_min)
    call f_conf%get_or_die("ssh_max", self%ssh_max)
    call f_conf%get_or_die("cicen_min", self%cicen_min)
    call f_conf%get_or_die("cicen_max", self%cicen_max)
    call f_conf%get_or_die("hicen_min", self%hicen_min)
    call f_conf%get_or_die("hicen_max", self%hicen_max)

  end subroutine soca_bkgerr_readbounds

  ! ------------------------------------------------------------------------------
  !> Setup the static background error
  subroutine soca_bkgerr_applybounds(self, fld)
    class(soca_bkgerr_bounds_type), intent(inout) :: self
    type(soca_field),               intent(inout) :: fld

    integer :: isc, iec, jsc, jec, i, j

    ! Apply config bounds to background error
    isc = fld%geom%isc ; iec = fld%geom%iec
    jsc = fld%geom%jsc ; jec = fld%geom%jec

    do i = isc, iec
       do j = jsc, jec
          ! Apply bounds
          fld%ssh(i,j) = soca_adjust(fld%ssh(i,j), &
                                     &self%ssh_min,&
                                     &self%ssh_max)
          fld%tocn(i,j,:) = soca_adjust(fld%tocn(i,j,:),&
                                        &self%t_min,&
                                        &self%t_max)
          fld%socn(i,j,:) = soca_adjust(fld%socn(i,j,:),&
                                        &self%s_min,&
                                        &self%s_max)
          fld%seaice%cicen(i,j,:) = soca_adjust(fld%seaice%cicen(i,j,:),&
                                        &self%cicen_min,&
                                        &self%cicen_max)
          fld%seaice%hicen(i,j,:) = soca_adjust(fld%seaice%hicen(i,j,:),&
                                        &self%hicen_min,&
                                        &self%hicen_max)
       end do
    end do


  end subroutine soca_bkgerr_applybounds

end module soca_bkgerrutil_mod
