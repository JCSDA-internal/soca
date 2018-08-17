!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_bkgerr_mod

  use kinds
  implicit none

  !> Fortran derived type to hold configuration D
  type :: soca_bkgerr_config

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

  ! ------------------------------------------------------------------------------

  subroutine soca_bkgerr_setup(c_conf, config, bkg)
    use iso_c_binding
    use config_mod
    use kinds
    use soca_fields
    use fms_io_mod,      only : fms_io_init, fms_io_exit

    implicit none

    type(c_ptr),                 intent(in) :: c_conf   !< The configuration
    type(soca_bkgerr_config), intent(inout) :: config   !< Config parameters for D

    type(soca_field) :: bkg
    type(soca_field) :: lensca    
    integer :: inzo, incat
    character(len=800) :: filename

    call create_copy(lensca, bkg)

    ! Setup horizontal length scale
    ! Ocean
    call bkg%geom%ocean%get_rossby_radius()    

    do inzo = 1,bkg%geom%ocean%nzo
       lensca%tocn(:,:,inzo) = max(200e3, 5.0*bkg%geom%ocean%rossby_radius)
       lensca%socn(:,:,inzo) = lensca%tocn(:,:,inzo)
       lensca%hocn(:,:,inzo) = lensca%tocn(:,:,inzo)
    end do
    lensca%ssh = lensca%tocn(:,:,1)

    ! Sea-ice
    do incat = 1,bkg%geom%ocean%ncat
       lensca%cicen(:,:,incat+1) = 100e3
       lensca%hicen(:,:,incat) = 100e3
    end do

    filename="rh.nc"
    call fld2file(lensca, filename)

    ! Setup vertical length scale
    ! Ocean
    do inzo = 1,bkg%geom%ocean%nzo
       lensca%tocn(:,:,inzo) = 20
       lensca%socn(:,:,inzo) = 20
       lensca%hocn(:,:,inzo) = 20
    end do
    lensca%ssh = 50

    ! Sea-ice
    do incat = 1,bkg%geom%ocean%ncat
       lensca%cicen(:,:,incat+1) = 1
       lensca%hicen(:,:,incat) = 1
    end do

    filename="rv.nc"
    call fld2file(lensca, filename)    
    
  end subroutine soca_bkgerr_setup
end module soca_bkgerr_mod
