!
! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_ocnsfc_mod

  use kinds
  use MOM_forcing_type,    only : forcing  
  use soca_geom_mod_c
  use random_mod
  use soca_model_geom_type, only : geom_get_domain_indices

  implicit none
  private

  type, public :: soca_ocnsfc_type
     real(kind=kind_real), allocatable :: sw_rad(:,:)
     real(kind=kind_real), allocatable :: lw_rad(:,:)     
     real(kind=kind_real), allocatable :: latent_heat(:,:)       
     real(kind=kind_real), allocatable :: sens_heat(:,:)
     real(kind=kind_real), allocatable :: fric_vel(:,:)
     real(kind=kind_real), allocatable :: mask(:,:)     
   contains
     procedure :: create => soca_ocnsfc_create
     procedure :: delete => soca_ocnsfc_delete
     procedure :: zeros => soca_ocnsfc_zeros
     procedure :: ones => soca_ocnsfc_ones
     procedure :: random => soca_ocnsfc_random
     procedure :: copy => soca_ocnsfc_copy
     procedure :: add => soca_ocnsfc_add
     procedure :: schur => soca_ocnsfc_schur
     procedure :: sub => soca_ocnsfc_sub
     procedure :: mul => soca_ocnsfc_mul
     procedure :: axpy => soca_ocnsfc_axpy
     procedure :: diff_incr => soca_ocnsfc_diff_incr
     procedure :: read_file => soca_ocnsfc_read_file
     procedure :: getforcing => soca_ocnsfc_getforcing
     procedure :: pushforcing => soca_ocnsfc_pushforcing     
     procedure :: applymask => soca_ocnsfc_applymask     
  end type soca_ocnsfc_type

contains

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_create(self, geom)
    class(soca_ocnsfc_type), intent(inout) :: self
    type(soca_geom),            intent(in) :: geom

    integer :: isd, ied, jsd, jed

    ! Indices for data domain (with halo)
    call geom_get_domain_indices(geom%ocean, "data   ", isd, ied, jsd, jed)    

    ! Allocate ocean state
    if (.not.allocated(self%sw_rad)) allocate(self%sw_rad(isd:ied,jsd:jed))
    if (.not.allocated(self%lw_rad)) allocate(self%lw_rad(isd:ied,jsd:jed))
    if (.not.allocated(self%latent_heat)) allocate(self%latent_heat(isd:ied,jsd:jed))
    if (.not.allocated(self%sens_heat)) allocate(self%sens_heat(isd:ied,jsd:jed))    
    if (.not.allocated(self%fric_vel)) allocate(self%fric_vel(isd:ied,jsd:jed))
    if (.not.allocated(self%mask)) allocate(self%mask(isd:ied,jsd:jed))
    self%mask = geom%ocean%mask2d
    
  end subroutine soca_ocnsfc_create

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_delete(self)
    class(soca_ocnsfc_type), intent(inout) :: self

    ! Deallocate all ocean surface state
    deallocate(self%sw_rad)
    deallocate(self%lw_rad)
    deallocate(self%latent_heat)
    deallocate(self%sens_heat)
    deallocate(self%fric_vel)

  end subroutine soca_ocnsfc_delete

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_zeros(self)
    class(soca_ocnsfc_type), intent(inout) :: self

    self%sw_rad      = 0.0_kind_real
    self%lw_rad      = 0.0_kind_real
    self%latent_heat = 0.0_kind_real
    self%sens_heat   = 0.0_kind_real
    self%fric_vel    = 0.0_kind_real

  end subroutine soca_ocnsfc_zeros

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_ones(self)
    class(soca_ocnsfc_type), intent(inout) :: self

    self%sw_rad      = 1.0_kind_real
    self%lw_rad      = 1.0_kind_real
    self%latent_heat = 1.0_kind_real
    self%sens_heat   = 1.0_kind_real
    self%fric_vel    = 1.0_kind_real

  end subroutine soca_ocnsfc_ones

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_random(self)
    class(soca_ocnsfc_type), intent(inout) :: self

    integer :: rseed = 1
    
    ! Deallocate all ocean surface state
    call normal_distribution(self%sw_rad, 0.0_kind_real, 1.0_kind_real, rseed)
    call normal_distribution(self%lw_rad, 0.0_kind_real, 1.0_kind_real, rseed)
    call normal_distribution(self%latent_heat, 0.0_kind_real, 1.0_kind_real, rseed)
    call normal_distribution(self%sens_heat, 0.0_kind_real, 1.0_kind_real, rseed)
    call normal_distribution(self%fric_vel, 0.0_kind_real, 1.0_kind_real, rseed)
    
  end subroutine soca_ocnsfc_random

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_copy(self, rhs)
    class(soca_ocnsfc_type), intent(inout) :: self
    class(soca_ocnsfc_type),    intent(in) :: rhs    

    self%sw_rad      = rhs%sw_rad
    self%lw_rad      = rhs%lw_rad
    self%latent_heat = rhs%latent_heat
    self%sens_heat   = rhs%sens_heat
    self%fric_vel    = rhs%fric_vel
    
  end subroutine soca_ocnsfc_copy

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_add(self, other)
    class(soca_ocnsfc_type), intent(inout) :: self
    class(soca_ocnsfc_type),    intent(in) :: other    

    self%sw_rad      = self%sw_rad      + other%sw_rad
    self%lw_rad      = self%lw_rad      + other%lw_rad
    self%latent_heat = self%latent_heat + other%latent_heat
    self%sens_heat   = self%sens_heat   + other%sens_heat
    self%fric_vel    = self%fric_vel    + other%fric_vel
    
  end subroutine soca_ocnsfc_add

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_schur(self, other)
    class(soca_ocnsfc_type), intent(inout) :: self
    class(soca_ocnsfc_type),    intent(in) :: other    

    self%sw_rad      = self%sw_rad      * other%sw_rad
    self%lw_rad      = self%lw_rad      * other%lw_rad
    self%latent_heat = self%latent_heat * other%latent_heat
    self%sens_heat   = self%sens_heat   * other%sens_heat
    self%fric_vel    = self%fric_vel    * other%fric_vel
    
  end subroutine soca_ocnsfc_schur
  
  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_sub(self, other)
    class(soca_ocnsfc_type), intent(inout) :: self
    class(soca_ocnsfc_type),    intent(in) :: other    

    self%sw_rad      = self%sw_rad      - other%sw_rad
    self%lw_rad      = self%lw_rad      - other%lw_rad
    self%latent_heat = self%latent_heat - other%latent_heat
    self%sens_heat   = self%sens_heat   - other%sens_heat
    self%fric_vel    = self%fric_vel    - other%fric_vel
    
  end subroutine soca_ocnsfc_sub
  
  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_mul(self, zz)
    class(soca_ocnsfc_type), intent(inout) :: self
    real(kind=kind_real),       intent(in) :: zz    

    self%sw_rad      = zz * self%sw_rad
    self%lw_rad      = zz * self%lw_rad
    self%latent_heat = zz * self%latent_heat
    self%sens_heat   = zz * self%sens_heat
    self%fric_vel    = zz * self%fric_vel
    
  end subroutine soca_ocnsfc_mul

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_axpy(self, zz, other)
    class(soca_ocnsfc_type), intent(inout) :: self
    real(kind=kind_real),       intent(in) :: zz
    class(soca_ocnsfc_type),    intent(in) :: other

    self%sw_rad      = self%sw_rad      + zz * other%sw_rad
    self%lw_rad      = self%lw_rad      + zz * other%lw_rad
    self%latent_heat = self%latent_heat + zz * other%latent_heat
    self%sens_heat   = self%sens_heat   + zz * other%sens_heat
    self%fric_vel    = self%fric_vel    + zz * other%fric_vel
    
  end subroutine soca_ocnsfc_axpy

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_diff_incr(self, x1, x2)
    class(soca_ocnsfc_type), intent(inout) :: self
    class(soca_ocnsfc_type),    intent(in) :: x1
    class(soca_ocnsfc_type),    intent(in) :: x2    

    self%sw_rad      = x1%sw_rad      - x2%sw_rad
    self%lw_rad      = x1%lw_rad      - x2%lw_rad
    self%latent_heat = x1%latent_heat - x2%latent_heat
    self%sens_heat   = x1%sens_heat   - x2%sens_heat
    self%fric_vel    = x1%fric_vel    - x2%fric_vel
    
  end subroutine soca_ocnsfc_diff_incr

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_getforcing(self, fluxes)
    class(soca_ocnsfc_type), intent(inout) :: self
    type(forcing),              intent(in) :: fluxes !< Thermodynamic forcing
    
    ! Get ocnsfc from mom6 forcing
    self%sw_rad      = - real(fluxes%sw, kind=kind_real)
    self%lw_rad      = - real(fluxes%lw, kind=kind_real)
    self%latent_heat = - real(fluxes%latent, kind=kind_real)
    self%sens_heat   = - real(fluxes%sens, kind=kind_real)
    self%fric_vel    = real(fluxes%ustar, kind=kind_real)
    call soca_ocnsfc_applymask(self)
    
  end subroutine soca_ocnsfc_getforcing

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_pushforcing(self, fluxes)
    class(soca_ocnsfc_type), intent(in) :: self
    type(forcing),        intent(inout) :: fluxes !< Thermodynamic forcing
    
    ! Push ocnsfc into mom6 forcing
    fluxes%sw     = - real(self%sw_rad, kind=8)
    fluxes%lw     = - real(self%lw_rad, kind=8)
    fluxes%latent = - real(self%latent_heat, kind=8)
    fluxes%sens   = - real(self%sens_heat, kind=8)
    fluxes%ustar  = real(self%fric_vel, kind=8)
    
  end subroutine soca_ocnsfc_pushforcing

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_applymask(self)
    class(soca_ocnsfc_type), intent(inout) :: self
    
    self%sw_rad      = self%sw_rad      * self%mask 
    self%lw_rad      = self%lw_rad      * self%mask 
    self%latent_heat = self%latent_heat * self%mask 
    self%sens_heat   = self%sens_heat   * self%mask 
    self%fric_vel    = self%fric_vel    * self%mask 
    
  end subroutine soca_ocnsfc_applymask

  ! ------------------------------------------------------------------------------  
  subroutine soca_ocnsfc_read_file(self)
    ! HACK, TODO: Do something, like read a file!
    class(soca_ocnsfc_type), intent(inout) :: self

    self%sw_rad      = -57.1_kind_real
    self%lw_rad      = 600.0_kind_real
    self%latent_heat = 103.0_kind_real
    self%sens_heat   = 7.7_kind_real
    self%fric_vel    = 0.08_kind_real
    
  end subroutine soca_ocnsfc_read_file
  
end module soca_ocnsfc_mod
