! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_bkgerrutil_mod
  use config_mod
  use iso_c_binding  
  use kinds
  use soca_fields
  use soca_utils
  use soca_model_geom_type, only : geom_get_domain_indices

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
  subroutine soca_bkgerr_readbounds(self, c_conf)
    class(soca_bkgerr_bounds_type), intent(inout) :: self
    type(c_ptr),                       intent(in) :: c_conf
    
    ! Get bounds from configuration
    self%t_min   = config_get_real(c_conf,"t_min")
    self%t_max   = config_get_real(c_conf,"t_max")
    self%s_min   = config_get_real(c_conf,"s_min")
    self%s_max   = config_get_real(c_conf,"s_max")
    self%ssh_min = config_get_real(c_conf,"ssh_min")
    self%ssh_max = config_get_real(c_conf,"ssh_max")
    self%cicen_min = config_get_real(c_conf,"cicen_min")
    self%cicen_max = config_get_real(c_conf,"cicen_max")
    self%hicen_min = config_get_real(c_conf,"hicen_min")
    self%hicen_max = config_get_real(c_conf,"hicen_max")
    
  end subroutine soca_bkgerr_readbounds

  ! ------------------------------------------------------------------------------
  !> Setup the static background error
  subroutine soca_bkgerr_applybounds(self, fld)
    class(soca_bkgerr_bounds_type), intent(inout) :: self    
    type(soca_field),               intent(inout) :: fld

    integer :: isc, iec, jsc, jec, i, j
    
    ! Apply config bounds to background error
    call geom_get_domain_indices(fld%geom%ocean, "compute", isc, iec, jsc, jec)    
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
