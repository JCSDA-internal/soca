! (C) Copyright 2020-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------
!> C++ interfaces for soca_getvalues_mod::soca_getvalues
module soca_getvalues_mod_c

use datetime_mod, only: datetime, c_f_datetime
use iso_c_binding
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod, only: ufo_geovals
use ufo_locations_mod, only: ufo_locations

! soca modules
use soca_geom_mod, only: soca_geom
use soca_geom_mod_c, only: soca_geom_registry
use soca_getvalues_mod, only: soca_getvalues
use soca_getvalues_reg, only: soca_getvalues_registry
use soca_state_mod, only: soca_state
use soca_state_reg, only: soca_state_registry
use soca_increment_mod, only: soca_increment
use soca_increment_reg, only: soca_increment_registry

implicit none
private


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> C++ interface for soca_getvalues_mod::soca_getvalues::create()
subroutine soca_getvalues_create_c(c_key_self, c_key_geom, c_locs) &
           bind (c, name='soca_getvalues_create_f90')
integer(c_int),     intent(inout) :: c_key_self      !< Key to self
integer(c_int),     intent(in)    :: c_key_geom      !< Key to geometry
type(c_ptr), value, intent(in)    :: c_locs

type(soca_getvalues), pointer :: self
type(soca_geom),      pointer :: geom
type(ufo_locations)           :: locs

! Create object
call soca_getvalues_registry%init()
call soca_getvalues_registry%add(c_key_self)
call soca_getvalues_registry%get(c_key_self, self)

! Others
call soca_geom_registry%get(c_key_geom, geom)
locs = ufo_locations(c_locs)

! Call method
call self%create(geom, locs)

end subroutine soca_getvalues_create_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for soca_getvalues_mod::soca_getvalues::delete()
subroutine soca_getvalues_delete_c(c_key_self) bind (c, name='soca_getvalues_delete_f90')
integer(c_int), intent(inout) :: c_key_self !< Key to self

type(soca_getvalues), pointer :: self

! Get object
call soca_getvalues_registry%get(c_key_self, self)

! Call method
call self%delete()

! Remove object
call soca_getvalues_registry%remove(c_key_self)

end subroutine soca_getvalues_delete_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for soca_getvalues_mod::soca_getvalues::fill_geovals()
subroutine soca_getvalues_fill_geovals_c(c_key_self, c_key_geom, c_key_state, c_t1, c_t2, &
                                         c_locs, c_key_geovals) &
           bind (c, name='soca_getvalues_fill_geovals_f90')

integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geom
integer(c_int),     intent(in) :: c_key_state
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
type(c_ptr), value, intent(in) :: c_locs
integer(c_int),     intent(in) :: c_key_geovals

type(soca_getvalues), pointer :: self
type(soca_geom),      pointer :: geom
type(soca_state),     pointer :: state
type(datetime)                :: t1
type(datetime)                :: t2
type(ufo_locations)           :: locs
type(ufo_geovals),    pointer :: geovals

! Get objects
call soca_getvalues_registry%get(c_key_self, self)
call soca_geom_registry%get(c_key_geom, geom)
call soca_state_registry%get(c_key_state, state)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
locs = ufo_locations(c_locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals(geom, state, t1, t2, locs, geovals)

end subroutine soca_getvalues_fill_geovals_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for soca_getvalues_mod::soca_getvalues::fill_geovals()
subroutine soca_getvalues_fill_geovals_tl_c(c_key_self, c_key_geom, c_key_incr, c_t1, c_t2, &
                                         c_locs, c_key_geovals) &
           bind (c, name='soca_getvalues_fill_geovals_tl_f90')

integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geom
integer(c_int),     intent(in) :: c_key_incr
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
type(c_ptr), value, intent(in) :: c_locs
integer(c_int),     intent(in) :: c_key_geovals

type(soca_getvalues), pointer :: self
type(soca_geom),      pointer :: geom
type(soca_increment), pointer :: incr
type(datetime)                :: t1
type(datetime)                :: t2
type(ufo_locations)           :: locs
type(ufo_geovals),    pointer :: geovals

! Get objects
call soca_getvalues_registry%get(c_key_self, self)
call soca_geom_registry%get(c_key_geom, geom)
call soca_increment_registry%get(c_key_incr, incr)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
locs = ufo_locations(c_locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals(geom, incr, t1, t2, locs, geovals)

end subroutine soca_getvalues_fill_geovals_tl_c


! --------------------------------------------------------------------------------------------------
!> C++ interface for soca_getvalues_mod::soca_getvalues::fill_geovals_ad()
subroutine soca_getvalues_fill_geovals_ad_c(c_key_self, c_key_geom, c_key_incr, c_t1, c_t2, &
                                            c_locs, c_key_geovals) &
           bind (c, name='soca_getvalues_fill_geovals_ad_f90')

integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_geom
integer(c_int),     intent(in) :: c_key_incr
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
type(c_ptr), value, intent(in) :: c_locs
integer(c_int),     intent(in) :: c_key_geovals

type(soca_getvalues), pointer :: self
type(soca_geom),      pointer :: geom
type(soca_increment), pointer :: incr
type(datetime)                :: t1
type(datetime)                :: t2
type(ufo_locations)           :: locs
type(ufo_geovals),    pointer :: geovals

! Get objects
call soca_getvalues_registry%get(c_key_self, self)
call soca_geom_registry%get(c_key_geom, geom)
call soca_increment_registry%get(c_key_incr, incr)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
locs = ufo_locations(c_locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals_ad(geom, incr, t1, t2, locs, geovals)

end subroutine soca_getvalues_fill_geovals_ad_c

end module soca_getvalues_mod_c
