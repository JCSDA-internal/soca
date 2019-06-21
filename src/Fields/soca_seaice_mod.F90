!
! (C) Copyright 2017-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_seaice_mod

  use config_mod
  use datetime_mod  
  use fms_io_mod, only : fms_io_init, fms_io_exit,&
       &register_restart_field, restart_file_type,&
       &restore_state, free_restart_type, save_restart
  use fms_mod,    only: read_data  
  use iso_c_binding  
  use kinds
  use MOM_forcing_type,    only : forcing  
  use random_mod
  use soca_fieldsutils_mod
  use soca_geom_mod_c  
  use soca_geom_mod

  implicit none
  private

  type, public :: soca_seaice_type
     real(kind=kind_real), allocatable :: cicen(:,:,:)
     real(kind=kind_real), allocatable :: hicen(:,:,:)
   contains
     procedure :: create => soca_seaice_create
     procedure :: delete => soca_seaice_delete
     procedure :: zeros => soca_seaice_zeros
     procedure :: ones => soca_seaice_ones
     procedure :: abs => soca_seaice_abs
     procedure :: random => soca_seaice_random
     procedure :: copy => soca_seaice_copy
     procedure :: add => soca_seaice_add
     procedure :: schur => soca_seaice_schur
     procedure :: sub => soca_seaice_sub
     procedure :: mul => soca_seaice_mul
     procedure :: axpy => soca_seaice_axpy
     procedure :: diff_incr => soca_seaice_diff_incr
     procedure :: read_restart => soca_seaice_read_rst
     procedure :: read_diag => soca_seaice_read_diag
     procedure :: write_restart => soca_seaice_write_rst     
  end type soca_seaice_type

  real(kind=kind_real), parameter :: soca_rho_ice = 905.0 ! [kg/m3]
  
contains

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_create(self, geom)
    class(soca_seaice_type), intent(inout) :: self
    type(soca_geom),            intent(in) :: geom

    integer :: isd, ied, jsd, jed

    ! Indices for data domain (with halo)
    call geom_get_domain_indices(geom, "data   ", isd, ied, jsd, jed)    

    ! Allocate sea-ice state
    if (.not.allocated(self%cicen)) allocate(self%cicen(isd:ied,jsd:jed,geom%ice_column%ncat+1))
    if (.not.allocated(self%hicen)) allocate(self%hicen(isd:ied,jsd:jed,geom%ice_column%ncat))
    
  end subroutine soca_seaice_create

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_delete(self)
    class(soca_seaice_type), intent(inout) :: self

    ! Deallocate all ocean surface state
    deallocate(self%cicen)
    deallocate(self%hicen)

  end subroutine soca_seaice_delete

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_zeros(self)
    class(soca_seaice_type), intent(inout) :: self

    self%cicen      = 0.0_kind_real
    self%hicen      = 0.0_kind_real

  end subroutine soca_seaice_zeros

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_ones(self)
    class(soca_seaice_type), intent(inout) :: self

    self%cicen      = 1.0_kind_real
    self%hicen      = 1.0_kind_real

  end subroutine soca_seaice_ones

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_abs(self)
    class(soca_seaice_type), intent(inout) :: self

    self%cicen      = abs(self%cicen)
    self%hicen      = abs(self%hicen)

  end subroutine soca_seaice_abs

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_random(self)
    class(soca_seaice_type), intent(inout) :: self

    integer :: rseed = 1
    
    call normal_distribution(self%cicen, 0.0_kind_real, 1.0_kind_real, rseed)
    call normal_distribution(self%hicen, 0.0_kind_real, 1.0_kind_real, rseed)
    
  end subroutine soca_seaice_random

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_copy(self, rhs)
    class(soca_seaice_type), intent(inout) :: self
    class(soca_seaice_type),    intent(in) :: rhs    

    self%cicen      = rhs%cicen
    self%hicen      = rhs%hicen
    
  end subroutine soca_seaice_copy

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_add(self, other)
    class(soca_seaice_type), intent(inout) :: self
    class(soca_seaice_type),    intent(in) :: other    

    self%cicen      = self%cicen      + other%cicen
    self%hicen      = self%hicen      + other%hicen
    
  end subroutine soca_seaice_add

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_schur(self, other)
    class(soca_seaice_type), intent(inout) :: self
    class(soca_seaice_type),    intent(in) :: other    

    self%cicen      = self%cicen      * other%cicen
    self%hicen      = self%hicen      * other%hicen
    
  end subroutine soca_seaice_schur
  
  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_sub(self, other)
    class(soca_seaice_type), intent(inout) :: self
    class(soca_seaice_type),    intent(in) :: other    

    self%cicen      = self%cicen      - other%cicen
    self%hicen      = self%hicen      - other%hicen
    
  end subroutine soca_seaice_sub
  
  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_mul(self, zz)
    class(soca_seaice_type), intent(inout) :: self
    real(kind=kind_real),       intent(in) :: zz    

    self%cicen      = zz * self%cicen
    self%hicen      = zz * self%hicen
    
  end subroutine soca_seaice_mul

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_axpy(self, zz, other)
    class(soca_seaice_type), intent(inout) :: self
    real(kind=kind_real),       intent(in) :: zz
    class(soca_seaice_type),    intent(in) :: other

    self%cicen      = self%cicen      + zz * other%cicen
    self%hicen      = self%hicen      + zz * other%hicen
    
  end subroutine soca_seaice_axpy

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_diff_incr(self, x1, x2)
    class(soca_seaice_type), intent(inout) :: self
    class(soca_seaice_type),    intent(in) :: x1
    class(soca_seaice_type),    intent(in) :: x2    

    self%cicen      = x1%cicen      - x2%cicen
    self%hicen      = x1%hicen      - x2%hicen
    
  end subroutine soca_seaice_diff_incr

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_read_rst(self, c_conf, geom, fldnames)
    class(soca_seaice_type), intent(inout) :: self
    type(c_ptr),                intent(in) :: c_conf
    type(soca_geom),            intent(in) :: geom    
    character(len=5),           intent(in) :: fldnames(:)

    integer, parameter :: max_string_length=800
    integer :: idr, i
    character(len=max_string_length) :: filename, basename
    character(len=4) :: seaice_model
    type(restart_file_type) :: restart

    ! Check what model we are reading a file from ('sis2' or 'cice')
    seaice_model = 'sis2' ! Default model is sis2    
    if (config_element_exists(c_conf,"seaice_model")) then
       seaice_model = config_get_string(c_conf,len(seaice_model),"seaice_model")
    end if

    if (config_element_exists(c_conf,"ice_filename")) then
       basename = config_get_string(c_conf,len(basename),"basename")
       filename = config_get_string(c_conf,len(filename),"ice_filename")
       filename = trim(basename)//trim(filename)       
    else
       ! Set seaice state to 0 if no file provided
       call self%zeros()
       return
    end if

    select case(seaice_model)
    case('sis2')
       call fms_io_init()
       do i = 1, size(fldnames)
          select case(fldnames(i))
          case('cicen')
             idr = register_restart_field(restart, filename, 'part_size', &
                  self%cicen(:,:,:), &
                  domain=geom%G%Domain%mpp_domain)
          case('hicen')
             idr = register_restart_field(restart, filename, 'h_ice', &
                  self%hicen(:,:,:), &
                  domain=geom%G%Domain%mpp_domain)
          end select
       end do
       call restore_state(restart, directory='')
       call free_restart_type(restart)    
       call fms_io_exit()
       ! Convert hicen from [kg/m2] to [m]
       self%hicen(:,:,:) = self%hicen(:,:,:)/soca_rho_ice

    case('cice')
       call fms_io_init()
       do i = 1, size(fldnames)
          select case(fldnames(i))
          case('cicen')
             idr = register_restart_field(restart, filename, 'aicen', &
                  self%cicen(:,:,2:), &
                  domain=geom%G%Domain%mpp_domain)
          case('hicen')
             idr = register_restart_field(restart, filename, 'vicen', &
                  self%hicen(:,:,:), &
                  domain=geom%G%Domain%mpp_domain)
          end select
          ! Add ocean fraction
          self%cicen(:,:,1) = 1.0_kind_real - sum(self%cicen(:,:,2:), dim=3)
          ! Convert to hicen
          where(self%cicen(:,:,2:)>0.0_kind_real)
             self%hicen(:,:,:) = self%hicen(:,:,:)/self%cicen(:,:,2:)
          end where
       end do
       call restore_state(restart, directory='')
       call free_restart_type(restart)    
       call fms_io_exit()
       
    case default
       call abor1_ftn("soca_seaice_mod: Reading for seaice model "//trim(seaice_model)//" not implemented")
    end select
  end subroutine soca_seaice_read_rst

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_read_diag(self, c_conf, geom, fldnames)
    class(soca_seaice_type), intent(inout) :: self
    type(c_ptr),                intent(in) :: c_conf
    type(soca_geom),            intent(in) :: geom    
    character(len=5),           intent(in) :: fldnames(:)

    integer, parameter :: max_string_length=800
    integer :: i
    character(len=max_string_length) :: filename

    if (config_element_exists(c_conf,"filename")) then
       filename = config_get_string(c_conf,len(filename),"filename")
    else
       call self%zeros()
       return
    end if

    call fms_io_init()
    do i = 1, size(fldnames)
       select case(fldnames(i))
       case('cicen')
          call read_data(filename,"cicen", &
                         self%cicen(:,:,:), &
                         domain=geom%G%Domain%mpp_domain)
       case('hicen')
          call read_data(filename,"hicen", &
                         self%hicen(:,:,:), &
                         domain=geom%G%Domain%mpp_domain)
       end select
    end do
    call fms_io_exit()

  end subroutine soca_seaice_read_diag

  ! ------------------------------------------------------------------------------  
  subroutine soca_seaice_write_rst(self, c_conf, geom, vdate)
    class(soca_seaice_type), intent(inout) :: self
    type(c_ptr),                intent(in) :: c_conf
    type(soca_geom),            intent(in) :: geom
    type(datetime),          intent(inout) :: vdate

    integer, parameter :: max_string_length=800
    integer :: idr
    character(len=max_string_length) :: filename
    character(len=4) :: seaice_model
    type(restart_file_type) :: restart
    real(kind=kind_real), allocatable :: vicen(:,:,:)
    real(kind=kind_real), allocatable :: aice(:,:), hice(:,:) ! Aggregates
    integer :: isd, ied, jsd, jed

    ! Check what model we are reading a file from ('sis2' or 'cice')
    seaice_model = 'sis2' ! Default model is sis2    
    if (config_element_exists(c_conf,"seaice_model")) then
       seaice_model = config_get_string(c_conf,len(seaice_model),"seaice_model")
    end if

    ! Register sea-ice fields
    call fms_io_init()

    ! Allocate and compute aggregate variables
    call geom_get_domain_indices(geom, "data   ", isd, ied, jsd, jed)    
    allocate(aice(isd:ied,jsd:jed))
    allocate(hice(isd:ied,jsd:jed))
    aice(:,:) = sum(self%cicen(:,:,2:), dim=3)
    hice(:,:) = sum(self%hicen(:,:,:), dim=3)
    
    select case(seaice_model)
    case('sis2')    
       ! Generate file names
       filename = soca_genfilename(c_conf,max_string_length,vdate,"ice")    

       ! Register sis2 variables
       idr = register_restart_field(restart, filename, 'part_size', self%cicen, &
            domain=geom%G%Domain%mpp_domain)
       idr = register_restart_field(restart, filename, 'h_ice', self%hicen, &
            domain=geom%G%Domain%mpp_domain)

    case('cice')
       ! Get ice volume
       call geom_get_domain_indices(geom, "data   ", isd, ied, jsd, jed)    
       allocate(vicen(isd:ied,jsd:jed,geom%ice_column%ncat))
       vicen(:,:,:) = self%cicen(:,:,2:)*self%hicen(:,:,:)
       
       ! Generate file names
       filename = soca_genfilename(c_conf,max_string_length,vdate,"cice")

       ! Register cice variables
       idr = register_restart_field(restart, filename, 'aicen', self%cicen(:,:,2:), &
            domain=geom%G%Domain%mpp_domain)

       idr = register_restart_field(restart, filename, 'vicen', vicen, &
            domain=geom%G%Domain%mpp_domain)

       deallocate(vicen)
    end select

    ! Register aggregate variables
    idr = register_restart_field(restart, filename, 'aice', aice, &
         domain=geom%G%Domain%mpp_domain)
    idr = register_restart_field(restart, filename, 'hice', hice, &
         domain=geom%G%Domain%mpp_domain)

    ! Write restart to disk
    call save_restart(restart, directory='')
    call free_restart_type(restart)    
    call fms_io_exit()

    
  end subroutine soca_seaice_write_rst
  
end module soca_seaice_mod
