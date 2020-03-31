! (C) Copyright 2020-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_state_mod

use soca_geom_mod
use soca_fields_mod
use soca_increment_mod
use oops_variables_mod
use soca_geostrophy_mod
use kinds, only: kind_real
use fckit_log_module, only: log, fckit_log

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

end type

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
subroutine soca_state_create(self, geom, vars)
  class(soca_state),         intent(inout) :: self
  type(soca_geom),  pointer, intent(inout) :: geom
  type(oops_variables),      intent(inout) :: vars

  ! additional internal fields need to be created
  call vars%push_back("mld")
  call vars%push_back("layer_depth")

  ! continue with normal fields initialization
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

  integer, save :: cnt_outer = 1
  character(len=800) :: filename, str_cnt
  type(soca_field), pointer :: fld, fld_r
  integer :: i, k

  real(kind=kind_real) :: min_ice = 1e-6_kind_real
  real(kind=kind_real) :: amin = 1e-6_kind_real
  real(kind=kind_real) :: amax = 10.0_kind_real
  real(kind=kind_real), allocatable :: alpha(:,:), aice_bkg(:,:), aice_ana(:,:)
  type(soca_geostrophy_type) :: geostrophy
  type(soca_fields) :: incr
  type(soca_field), pointer :: t, s, u, v, h, dt, ds, du, dv


  ! make sure rhs is a subset of self
  call rhs%check_subset(self)

  ! Make a copy of the increment
  call incr%copy(rhs)

  ! Compute geostrophic increment
  ! TODO Move inside of the balance operator.
  !      Will need to be removed when assimilating ocean current or when using
  !      ensemble derived increments (ensemble cross-covariances that include currents)
  if (self%has('hocn').and.self%has('tocn').and.self%has('socn').and.&
      self%has('uocn').and.self%has('vocn')) then
     ! Get necessary background fields needed to compute geostrophic perturbations
     call self%get("tocn", t)
     call self%get("socn", s)
     call self%get("hocn", h)

     ! Make a copy of the increment and get the needed pointers
     !call incr_geo%copy(rhs)
     call incr%get("tocn", dt)
     call incr%get("socn", ds)
     call incr%get("uocn", du)
     call incr%get("vocn", dv)

     ! Compute the geostrophic increment
     call geostrophy%setup(self%geom, h%val)
     call geostrophy%tl(h%val, t%val, s%val,&
          dt%val, ds%val, du%val, dv%val, self%geom)
     call geostrophy%delete()
  end if

  ! Colocate increment fields with h-grid
  call incr%colocate('h')

  ! for each field that exists in incr, add to self
  do i=1,size(incr%fields)
    fld_r => incr%fields(i)
    call self%get(fld_r%name, fld)

    select case (fld%name)
    case ('cicen')
      ! NOTE: sea ice concentration is special

      ! compute background rescaling
      allocate(alpha, mold=fld_r%val(:,:,1))
      allocate(aice_bkg, mold=alpha)
      allocate(aice_ana, mold=alpha)
      aice_bkg  = sum(fld%val(:,:,:), dim=3)
      aice_ana  = aice_bkg + sum(fld_r%val(:,:,:), dim=3)
      where (aice_ana < 0.0_kind_real) aice_ana = 0.0_kind_real
      where (aice_ana > 1.0_kind_real) aice_ana = 1.0_kind_real
      alpha = 1.0_kind_real
      where ( aice_bkg > min_ice) alpha = aice_ana / aice_bkg

      ! limit size of increment
      where ( alpha > amax ) alpha = amax
      where ( alpha < amin ) alpha = amin

      ! add fraction of increment
      do k=1,fld%nz
        fld%val(:,:,k) = alpha * fld%val(:,:,k)
      end do

    case default
      ! everyone else is normal
      fld%val = fld%val + fld_r%val
    end select
  end do

  ! Save increment for outer loop cnt_outer
  write(str_cnt,*) cnt_outer
  filename='incr.'//adjustl(trim(str_cnt))//'.nc'

  ! Save increment and clean memory
  call incr%write_file(filename)
  call incr%delete()

  ! Update outer loop counter
  cnt_outer = cnt_outer + 1

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


end module
