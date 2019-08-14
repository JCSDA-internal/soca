! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module SOCALOCALIZATION
  use iso_c_binding
  use soca_constants
  use soca_geom_mod, only: soca_geom
  use soca_geom_interface_mod, only: soca_geom_registry
  use soca_fields_mod, only: soca_field
  use soca_fields_interface_mod, only: soca_fields_registry
  use soca_covariance_mod

  implicit none

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------

  subroutine soca_localization_mult(c_key_conf, c_key_xincr) bind(c,name='soca_localization_mult_f90')
    integer(c_int), intent(in) :: c_key_conf
    integer(c_int), intent(in) :: c_key_xincr

    type(soca_cov),   pointer :: conf   !< Config structure
    type(soca_field), pointer :: xincr

    call soca_cov_registry%get(c_key_conf,conf)
    call soca_field_registry%get(c_key_xincr,xincr)
    call abor1_ftn("localization: not implemented")

  end subroutine soca_localization_mult

  ! ------------------------------------------------------------------------------

  subroutine soca_localization_setup(c_key_conf, c_model, c_key_geom) bind(c,name='soca_localization_setup_f90')
    integer(c_int), intent(inout) :: c_key_conf
    type(c_ptr), intent(in)    :: c_model  !< The configuration
    integer(c_int), intent(in) :: c_key_geom !< Geometry
    type(soca_cov), pointer :: conf !< covar structure
    type(soca_geom), pointer :: geom     !< Geometry

    call abor1_ftn("localization: not implemented")

    call soca_cov_registry%init()
    call soca_cov_registry%add(c_key_conf)
    call soca_cov_registry%get(c_key_conf, conf)
    call soca_geom_registry%get(c_key_geom, geom)

    return
  end subroutine soca_localization_setup

  ! ------------------------------------------------------------------------------

  subroutine soca_localization_delete(c_key_self) bind(c,name='soca_localization_delete_f90')
    integer(c_int), intent(inout) :: c_key_self
    type(soca_cov), pointer :: self

    call abor1_ftn("localization: not implemented")
    call soca_cov_registry%get(c_key_self, self)

  end subroutine soca_localization_delete

  ! ------------------------------------------------------------------------------

end module SOCALOCALIZATION
