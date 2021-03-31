! (C) Copyright 2020-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_state_mod

use soca_geom_mod
use soca_fields_mod
use soca_increment_mod
use soca_convert_state_mod
use oops_variables_mod
use kinds, only: kind_real
use fckit_log_module, only: fckit_log

implicit none
private

type, public, extends(soca_fields) :: soca_state

contains

  ! constructors / destructors
  procedure :: create => soca_state_create

  ! interactions with increment
  procedure :: diff_incr=> soca_state_diff_incr
  procedure :: add_incr => soca_state_add_incr

  ! misc
  procedure :: rotate => soca_state_rotate
  procedure :: convert => soca_state_convert
  procedure :: logexpon => soca_state_logexpon

end type

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine soca_state_create(self, geom, vars)
  class(soca_state),         intent(inout) :: self
  type(soca_geom),  pointer, intent(inout) :: geom
  type(oops_variables),      intent(inout) :: vars

  ! continue with fields initialization by base class
  call self%soca_fields%create(geom, vars)
end subroutine soca_state_create


! ------------------------------------------------------------------------------
!> Rotate horizontal vector
subroutine soca_state_rotate(self, coordinate, uvars, vvars)
  class(soca_state),  intent(inout) :: self
  character(len=*),      intent(in) :: coordinate ! "north" or "grid"
  type(oops_variables),  intent(in) :: uvars
  type(oops_variables),  intent(in) :: vvars

  integer :: z, i
  type(soca_field), pointer :: uocn, vocn
  real(kind=kind_real), allocatable :: un(:,:,:), vn(:,:,:)
  character(len=64) :: u_names, v_names

  do i=1, uvars%nvars()
    ! get (u, v) pair and make a copy
    u_names = trim(uvars%variable(i))
    v_names = trim(vvars%variable(i))
    if (self%has(u_names).and.self%has(v_names)) then
      call fckit_log%info("rotating "//trim(u_names)//" "//trim(v_names))
      call self%get(u_names, uocn)
      call self%get(v_names, vocn)
    else
      ! Skip if no pair found.
      call fckit_log%info("not rotating "//trim(u_names)//" "//trim(v_names))
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
!> add a set of increments to the set of fields
subroutine soca_state_add_incr(self, rhs)
  class(soca_state),  intent(inout) :: self
  class(soca_increment), intent(in) :: rhs

  type(soca_field), pointer :: fld, fld_r
  integer :: i, k

  real(kind=kind_real) :: min_ice = 1e-6_kind_real
  real(kind=kind_real) :: amin = 1e-6_kind_real
  real(kind=kind_real) :: amax = 10.0_kind_real
  real(kind=kind_real), allocatable :: alpha(:,:), aice_bkg(:,:), aice_ana(:,:)
  type(soca_fields) :: incr

  ! make sure rhs is a subset of self
  call rhs%check_subset(self)

  ! Make a copy of the increment
  call incr%copy(rhs)


  ! for each field that exists in incr, add to self
  do i=1,size(incr%fields)
    fld_r => incr%fields(i)
    call self%get(fld_r%name, fld)
    fld%val = fld%val + fld_r%val
  end do

end subroutine soca_state_add_incr


! ------------------------------------------------------------------------------
!> subtract two sets of fields, saving the results separately
subroutine soca_state_diff_incr(x1, x2, inc)
  class(soca_state),      intent(in)    :: x1
  class(soca_state),      intent(in)    :: x2
  class(soca_increment), intent(inout)  :: inc

  integer :: i
  type(soca_field), pointer :: f1, f2

  ! make sure fields correct shapes
  call inc%check_subset(x2)
  call x2%check_subset(x1)

  ! subtract
  do i=1,size(inc%fields)
    call x1%get(inc%fields(i)%name, f1)
    call x2%get(inc%fields(i)%name, f2)
    inc%fields(i)%val = f1%val - f2%val
  end do
end subroutine soca_state_diff_incr

! ------------------------------------------------------------------------------
!> ConvertState app
subroutine soca_state_convert(self, rhs)
  class(soca_state), intent(inout) :: self ! target
  class(soca_state), intent(in)   :: rhs   ! source
  integer :: n
  type(soca_convertstate_type) :: convert_state
  type(soca_field), pointer :: field1, field2, hocn1, hocn2

  call rhs%get("hocn", hocn1)
  call self%get("hocn", hocn2)
  call convert_state%setup(rhs%geom, self%geom, hocn1, hocn2)
  do n = 1, size(rhs%fields)
    field1 => rhs%fields(n)
    call self%get(trim(field1%name),field2)
    if (field1%fieldspec%io_file=="ocn" .or. field1%fieldspec%io_file=="sfc" .or. field1%fieldspec%io_file=="ice")  &
    call convert_state%change_resol(field1, field2, rhs%geom, self%geom)
  end do !n
  call convert_state%clean()
end subroutine soca_state_convert

! ------------------------------------------------------------------------------
!> Apply logarithmic and exponential transformations
subroutine soca_state_logexpon(self, transfunc, trvars)
  class(soca_state),  intent(inout) :: self
  character(len=*),      intent(in) :: transfunc ! "log" or "expon"
  type(oops_variables),  intent(in) :: trvars

  integer :: z, i
  type(soca_field), pointer :: trocn
  real(kind=kind_real), allocatable :: trn(:,:,:)
  real(kind=kind_real) :: min_val = 1e-6_kind_real
  character(len=64) :: tr_names

  do i=1, trvars%nvars()
    ! get a list variables to be transformed and make a copy
    tr_names = trim(trvars%variable(i))
    if (self%has(tr_names)) then
      call fckit_log%info("transforming "//trim(tr_names))
      call self%get(tr_names, trocn)
    else
      ! Skip if no variable found.
      call fckit_log%info("not transforming "//trim(tr_names))
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
