! (C) Copyright 2017-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_write_jacobian_mod

use kinds, only: kind_real
use netcdf
use soca_geom_mod, only: soca_geom
use fms_mod, only: write_data, set_domain
use fms_io_mod, only: fms_io_init, fms_io_exit

implicit none
private
public :: write_jacobian_to_netcdf

contains

subroutine write_jacobian_to_netcdf(jacobian, geom, filename, varname)
    real(kind=kind_real), intent(in) :: jacobian(:,:,:)
    type(soca_geom), intent(in) :: geom
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: varname


    call fms_io_init()
    call set_domain(geom%Domain%mpp_domain)

    ! Write the jacobian data from all PEs
    call write_data(filename, varname, jacobian, geom%Domain%mpp_domain)

    call fms_io_exit()
end subroutine write_jacobian_to_netcdf

end module soca_write_jacobian_mod
