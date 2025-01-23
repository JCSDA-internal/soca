! (C) Copyright 2020-2024 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> State fields
module soca_state_mod

use logger_mod
use kinds, only: kind_real
use oops_variables_mod
use atlas_module, only: atlas_field

! soca modules
use soca_geom_mod
use soca_fields_mod
use soca_increment_mod

implicit none
private


!-------------------------------------------------------------------------------
!> State fields.
!!
!! Any procedures that are shared with soca_increment are implemented
!! in the soca_fields base class
type, public, extends(soca_fields) :: soca_state

contains


  !> \name misc
  !! \{

  !! TODO(travis) These remaning subroutines should probably be removed, and instead
  !! live on as a non-linear variable change or saber outer block.... someday

  !> \copybrief soca_state_rotate \see soca_state_rotate
  procedure :: rotate => soca_state_rotate

  !> \copybrief soca_state_logexpon \see soca_state_logexpon
  procedure :: logexpon => soca_state_logexpon

  !> \}

end type


!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Rotate horizontal vector
!!
!! One or more sets of vectors, represented by corresponding u and v variables
!! in the \p uvars and \p vvars lists are rotated to north (if \p coordinate == "north")
!! or rotated back to the grid (if \p coordinate == "grid")
!! \relates soca_state_mod::soca_state
subroutine soca_state_rotate(self, coordinate, uvars, vvars)
  class(soca_state),  intent(inout) :: self
  character(len=*),      intent(in) :: coordinate !< "north" or "grid"
  type(oops_variables),  intent(in) :: uvars !< list of one or more U variables
  type(oops_variables),  intent(in) :: vvars !< list of one or more V variables

  integer :: i, j, k, n, idx

  type(atlas_field) :: ufield, vfield
  real(kind=kind_real), pointer :: udata(:,:), vdata(:,:)
  real(kind=kind_real) :: u, v

  character(len=64) :: u_names, v_names

  do n=1, uvars%nvars()
    ! get (u, v) pair
    u_names = trim(uvars%variable(n))
    v_names = trim(vvars%variable(n))

    if (.not. (self%has(u_names).and.self%has(v_names))) then
      ! Skip if no pair found.
      call oops_log%info("not rotating "//trim(u_names)//" "//trim(v_names))
      cycle
    end if

    call oops_log%info("rotating "//trim(u_names)//" "//trim(v_names))
    ufield = self%afieldset%field(u_names)
    vfield = self%afieldset%field(v_names)
    call ufield%data(udata)
    call vfield%data(vdata)

    select case(trim(coordinate))
    case("north")
      do j=self%geom%jsc,self%geom%jec
        do i=self%geom%isc,self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          do k=1,ufield%shape(1)
            u = udata(k, idx)
            v = vdata(k, idx)
            udata(k, idx) = self%geom%cos_rot(i,j)*u + self%geom%sin_rot(i,j) * v
            vdata(k, idx) = -self%geom%sin_rot(i,j)*u + self%geom%cos_rot(i,j) * v
            ! TODO should apply mask, not sure if I care enough to do it
          end do
        end do
      end do
    case("grid")
      do j=self%geom%jsc,self%geom%jec
        do i=self%geom%isc,self%geom%iec
          idx = self%geom%atlas_ij2idx(i,j)
          do k=1,ufield%shape(1)
            print *, "DBG grid "
            u = udata(k, idx)
            v = vdata(k, idx)
            udata(k, idx) = self%geom%cos_rot(i,j)*u - self%geom%sin_rot(i,j) * v
            vdata(k, idx) = self%geom%sin_rot(i,j)*u + self%geom%cos_rot(i,j) * v
            ! TODO should apply mask, not sure if I care enough to do it
          end do
        end do
      end do
    end select

    call ufield%set_dirty()
    call vfield%set_dirty()
    call ufield%final()
    call vfield%final()
  end do
end subroutine soca_state_rotate


! ------------------------------------------------------------------------------
!> Apply logarithmic and exponential transformations
!!
!! \relates soca_state_mod::soca_state
subroutine soca_state_logexpon(self, transfunc, trvars)
  class(soca_state),  intent(inout) :: self
  character(len=*),      intent(in) :: transfunc !< "log" or "expon"
  type(oops_variables),  intent(in) :: trvars !< list of variables to transform

  integer :: z, i
  type(soca_field), pointer :: trocn
  real(kind=kind_real), allocatable :: trn(:,:,:)
  real(kind=kind_real) :: min_val = 1e-6_kind_real
  character(len=64) :: tr_names

  do i=1, trvars%nvars()
    ! get a list variables to be transformed and make a copy
    tr_names = trim(trvars%variable(i))
    if (self%has(tr_names)) then
      call oops_log%info("transforming "//trim(tr_names))
      call self%get(tr_names, trocn)
    else
      ! Skip if no variable found.
      call oops_log%info("not transforming "//trim(tr_names))
      cycle
    end if
    allocate(trn(size(trocn%val,1),size(trocn%val,2),size(trocn%val,3)))
    trn = trocn%val

    select case(trim(transfunc))
    case("log")   ! apply logarithmic transformation
      trocn%val = log(trn + min_val)
    case("expon") ! Apply exponential transformation
      trocn%val = exp(trn) - min_val
    end select

    ! update halos
    call trocn%update_halo(self%geom)

    ! deallocate trn for next variable
    deallocate(trn)
  end do
end subroutine soca_state_logexpon
! ------------------------------------------------------------------------------


end module
