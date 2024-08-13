! (C) Copyright 2020-2024 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> State fields
module soca_state_mod

use logger_mod
use kinds, only: kind_real
use oops_variables_mod

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

  !> \copybrief soca_state_convert \see soca_state_convert
  procedure :: convert => soca_state_convert

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

  integer :: z, i
  type(soca_field), pointer :: uocn, vocn
  real(kind=kind_real), allocatable :: un(:,:,:), vn(:,:,:)
  character(len=64) :: u_names, v_names

  do i=1, uvars%nvars()
    ! get (u, v) pair and make a copy
    u_names = trim(uvars%variable(i))
    v_names = trim(vvars%variable(i))
    if (self%has(u_names).and.self%has(v_names)) then
      call oops_log%info("rotating "//trim(u_names)//" "//trim(v_names))
      call self%get(u_names, uocn)
      call self%get(v_names, vocn)
    else
      ! Skip if no pair found.
      call oops_log%info("not rotating "//trim(u_names)//" "//trim(v_names))
      cycle
    end if
    allocate(un(size(uocn%val,1),size(uocn%val,2),size(uocn%val,3)))
    allocate(vn(size(uocn%val,1),size(uocn%val,2),size(uocn%val,3)))
    un = uocn%val
    vn = vocn%val

    select case(trim(coordinate))
    case("north")   ! rotate (uocn, vocn) to geo north
      do z=1,uocn%nz
        uocn%val(:,:,z) = &
        (self%geom%cos_rot(:,:)*un(:,:,z) + self%geom%sin_rot(:,:)*vn(:,:,z)) * uocn%mask(:,:)
        vocn%val(:,:,z) = &
        (- self%geom%sin_rot(:,:)*un(:,:,z) + self%geom%cos_rot(:,:)*vn(:,:,z)) * vocn%mask(:,:)
      end do
    case("grid")
      do z=1,uocn%nz
        uocn%val(:,:,z) = &
        (self%geom%cos_rot(:,:)*un(:,:,z) - self%geom%sin_rot(:,:)*vn(:,:,z)) * uocn%mask(:,:)
        vocn%val(:,:,z) = &
        (self%geom%sin_rot(:,:)*un(:,:,z) + self%geom%cos_rot(:,:)*vn(:,:,z)) * vocn%mask(:,:)
      end do
    end select
    deallocate(un, vn)

    ! update halos
    call uocn%update_halo(self%geom)
    call vocn%update_halo(self%geom)
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
