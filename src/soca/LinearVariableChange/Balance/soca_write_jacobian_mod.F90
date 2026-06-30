! (C) Copyright 2017-2026 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Utility module to write Jacobian fields to NetCDF files
module soca_write_jacobian_mod

use kinds, only: kind_real
use soca_geom_mod, only: soca_geom
use soca_io_mod, only: soca_io_writer

implicit none
private
public :: write_jacobian_to_netcdf

contains

!> Write one or two Jacobian fields to a NetCDF file
!!
!! @param[in] jacobian  First Jacobian field to write (3D array)
!! @param[in] geom      SOCA geometry object
!! @param[in] filename  Output NetCDF filename
!! @param[in] varname   Name of the first variable in the NetCDF file
!! @param[in] jacobian2 Second optional Jacobian field to write (3D array)
!! @param[in] varname2  Optional name of the second variable in the NetCDF file
subroutine write_jacobian_to_netcdf(jacobian, geom, filename, varname, jacobian2, varname2)
    real(kind=kind_real), intent(in) :: jacobian(:,:,:)
    type(soca_geom), intent(in) :: geom
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname
    real(kind=kind_real), intent(in), optional :: jacobian2(:,:,:)
    character(len=*), intent(in), optional :: varname2

    type(soca_io_writer) :: writer
    real(kind=kind_real), target, allocatable :: buf1(:,:,:), buf2(:,:,:)

    ! Allocate and fill data-domain bounds buffer for variable 1
    allocate(buf1(geom%isd:geom%ied, geom%jsd:geom%jed, size(jacobian, 3)))
    buf1 = 0.0_kind_real
    buf1(geom%isc:geom%iec, geom%jsc:geom%jec, :) = jacobian

    ! Initialize writer with geometry domain and file
    call writer%init(geom%Domain%mpp_domain, filename)

    ! Enqueue the first jacobian
    call writer%enqueue(varname, buf1)

    ! Enqueue the second jacobian if provided
    if (present(jacobian2) .and. present(varname2)) then
        allocate(buf2(geom%isd:geom%ied, geom%jsd:geom%jed, size(jacobian2, 3)))
        buf2 = 0.0_kind_real
        buf2(geom%isc:geom%iec, geom%jsc:geom%jec, :) = jacobian2
        call writer%enqueue(varname2, buf2)
    end if

    ! Commit the write (this does the mpp gathers and netcdf output)
    call writer%commit()

    ! Deallocate local buffers
    deallocate(buf1)
    if (allocated(buf2)) deallocate(buf2)

end subroutine write_jacobian_to_netcdf

end module soca_write_jacobian_mod
