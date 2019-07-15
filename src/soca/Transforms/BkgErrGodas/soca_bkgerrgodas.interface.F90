
!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

! ------------------------------------------------------------------------------
!> Constructor for D (standard deviation of background error)
subroutine c_soca_bkgerrgodas_setup(c_key_self, c_conf, c_key_bkg) &
     &bind(c,name='soca_bkgerrgodas_setup_f90')
  use iso_c_binding
  use soca_bkgerrgodas_mod
  use soca_fields_mod_c

  integer(c_int), intent(inout) :: c_key_self   !< The D structure
  type(c_ptr),       intent(in) :: c_conf       !< The configuration
  integer(c_int), intent(in)    :: c_key_bkg    !< Background field

  type(soca_field), pointer :: bkg
  type(soca_bkgerrgodas_config), pointer :: self

  call soca_bkgerrgodas_registry%init()
  call soca_bkgerrgodas_registry%add(c_key_self)
  call soca_bkgerrgodas_registry%get(c_key_self, self)
  call soca_field_registry%get(c_key_bkg, bkg)

  call soca_bkgerrgodas_setup(c_conf, self, bkg)
  
end subroutine c_soca_bkgerrgodas_setup

! ------------------------------------------------------------------------------
!> Destructor for D
subroutine c_soca_bkgerrgodas_delete(c_key_self) bind(c,name='soca_bkgerrgodas_delete_f90')
  use iso_c_binding
  use soca_bkgerrgodas_mod
  use soca_fields_mod_c

  implicit none
  integer(c_int), intent(inout) :: c_key_self
  type(soca_bkgerrgodas_config), pointer :: self

  call soca_bkgerrgodas_registry%get(c_key_self, self)
  if (associated(self%bkg)) nullify(self%bkg)
  call delete(self%std_bkgerr)
  
  call soca_bkgerrgodas_registry%remove(c_key_self)

end subroutine c_soca_bkgerrgodas_delete

! ------------------------------------------------------------------------------
!> Multiplication forward and adjoint
subroutine c_soca_bkgerrgodas_mult_f90(c_key_self, c_key_a, c_key_m)&
     &bind(c,name='soca_bkgerrgodas_mult_f90')
  use iso_c_binding
  use soca_bkgerrgodas_mod
  use soca_fields_mod_c
  use kinds
  use soca_kst_mod
  
  implicit none
  integer(c_int), intent(in) :: c_key_a     !<    "   to Increment in
  integer(c_int), intent(in) :: c_key_m     !<    "   to Increment out 
  integer(c_int), intent(in) :: c_key_self 

  type(soca_field), pointer :: dxa
  type(soca_field), pointer :: dxm
  type(soca_bkgerrgodas_config), pointer :: self
  
  call soca_field_registry%get(c_key_a,dxa)
  call soca_field_registry%get(c_key_m,dxm)
  call soca_bkgerrgodas_registry%get(c_key_self,self)  

  !< Computes dxm = D dxa
  call copy(dxm, dxa)
  call soca_bkgerrgodas_mult(self, dxa, dxm)
  
end subroutine c_soca_bkgerrgodas_mult_f90

