! (C) Copyright 2023-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module soca_diffusion_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_min, fckit_mpi_max, fckit_mpi_sum
use kinds, only: kind_real
use logger_mod
use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad
use random_mod
use fms_io_mod
use oops_variables_mod

use soca_increment_mod
use soca_geom_mod, only : soca_geom

implicit none
private

! ------------------------------------------------------------------------------
! convenience macros, because I'm TIRED of the unreadability of endless isc, jsc .....
#define DOMAIN                  self%geom%isc:self%geom%iec,self%geom%jsc:self%geom%jec
#define DOMAIN_WITH_HALO        self%geom%isd:self%geom%ied,self%geom%jsd:self%geom%jed
#define LOOP_DOMAIN_I           self%geom%isc, self%geom%iec
#define LOOP_DOMAIN_J           self%geom%jsc, self%geom%jec
#define LOOP_DOMAIN_WITH_HALO_I self%geom%isd, self%geom%ied
#define LOOP_DOMAIN_WITH_HALO_J self%geom%jsd, self%geom%jed

! ------------------------------------------------------------------------------
! The diffusion parameters needed by a group of variables.
! More than one model variable can be part of a group, they will all share the
! same correlation lengths.
type :: soca_diffusion_group_params
 character(len=:),allocatable :: name 
 real(kind_real), allocatable :: KhDt(:,:)               !< horizontal diffusion coefficient
 real(kind_real), allocatable :: KvDt(:,:,:)             !< vertical diffusion
 real(kind_real), allocatable :: normalization_hz(:,:)   !< horizontal normalization constant
 real(kind_real), allocatable :: normalization_vt(:,:,:) !< vertical normalization constant
 integer                      :: niter_hz = -1           !< number of iterations for horizontal diffusion,
 integer                      :: niter_vt = -1           !< number of iterations for vertical diffusion,
                                                         ! (set to -1 to indicate has not been initialized) 
end type soca_diffusion_group_params

! ------------------------------------------------------------------------------
! do the mapping between a group name and one or more variables names
type :: soca_diffusion_group_mapping
character(len=:), allocatable :: group_name
type(oops_variables)          :: variables
end type

! ------------------------------------------------------------------------------
! The main EXPLICIT_DIFFUSION class.
type, public :: soca_diffusion
 private
  ! grid metrics
  real(kind_real), allocatable :: inv_sqrt_area(:,:) !< 1/sqrt(area)
  real(kind_real), allocatable :: dx(:,:)     !< cell spacing, x-direction (m)
  real(kind_real), allocatable :: dy(:,:)     !< cell spacing, y-direction (m)
  real(kind_real), allocatable :: pmon_u(:,:) !< pm/pn at u points (pm = 1/dx, pn = 1/dy)
  real(kind_real), allocatable :: pnom_v(:,:) !< pn/pm at v points
  real(kind_real), allocatable :: mask(:,:)   !< 1.0 where water, 0.0 where land
  
  ! parameters calculated during calibration() / read()
  ! "group" contains the weights / normalizations for one or more variables
  type(soca_diffusion_group_params), allocatable :: group(:)
  type(soca_diffusion_group_mapping), allocatable :: group_mapping(:)

  class(soca_geom), pointer :: geom    

contains
  procedure :: init => soca_diffusion_init
  procedure :: calibrate => soca_diffusion_calibrate
  procedure :: multiply => soca_diffusion_multiply
  procedure :: write_params => soca_diffusion_write_params
  procedure :: read_params => soca_diffusion_read_params

  ! calibration of horizontal parameters
  procedure, private :: calibrate_hz => soca_diffusion_calibrate_hz
  procedure, private :: calibrate_norm_hz_bruteforce => soca_diffusion_calibrate_norm_hz_bruteforce
  procedure, private :: calibrate_norm_hz_randomization => soca_diffusion_calibrate_norm_hz_randomization
  
  ! calibration of vertical parameters
  procedure, private :: calibrate_vt => soca_diffusion_calibrate_vt
  procedure, private :: calibrate_norm_vt => soca_diffusion_calibrate_norm_vt

  ! diffusion operators (hz and vt) (tl and ad)
  procedure, private :: diffusion_hz_tl => soca_diffusion_hz_tl
  procedure, private :: diffusion_hz_ad => soca_diffusion_hz_ad
  procedure, private :: diffusion_vt_tl => soca_diffusion_vt_tl
  procedure, private :: diffusion_vt_ad => soca_diffusion_vt_ad
  
  ! helper function to calculate stats (min/max/mean) of a field
  generic, private :: calc_stats => &
    soca_diffusion_calc_stats_2D, &
    soca_diffusion_calc_stats_3D
  procedure, private :: soca_diffusion_calc_stats_2D
  procedure, private :: soca_diffusion_calc_stats_3D

end type soca_diffusion

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
! calculate the masked global min/max/mean values for a given field
! TODO: is there any chance we need to also operate on unmasked fields?!?
subroutine soca_diffusion_calc_stats_3D(self, field, stats)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable, intent(in) :: field(:,:,:)
  real(kind=kind_real),             intent(out) :: stats(3)

  real(kind=kind_real), allocatable :: field_2d(:,:)
  integer :: z
  real(kind=kind_real) :: stats_tmp(3)

  stats(1) =huge(stats(1))
  stats(2) = -huge(stats(2))
  stats(3) = 0.0
  allocate(field_2d(size(field, dim=1), size(field,dim=2)))
  do z =1, size(field, dim=3)
    field_2d = field(:,:,z)
    call self%calc_stats(field_2d, stats_tmp)
    stats(1) = min(stats(1), stats_tmp(1))
    stats(2) = max(stats(2), stats_tmp(2))
    stats(3) = (stats(3)*(z-1) + stats_tmp(3) ) / z ! eh, close enough
  end do  
end subroutine

subroutine soca_diffusion_calc_stats_2D(self, field, stats)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable, intent(in) :: field(:,:)
  real(kind=kind_real),             intent(out) :: stats(3)

  real(kind=kind_real) :: l_min, l_max, l_sum, l_count, g_count

  l_min =  minval(field(DOMAIN), &
                  mask=self%mask(DOMAIN)==1.0)
  l_max =  maxval(field(DOMAIN), &
                  mask=self%mask(DOMAIN)==1.0)
  l_sum =     sum(field(DOMAIN), &
                  mask=self%mask(DOMAIN)==1.0)
  l_count = count(self%mask(DOMAIN)==1.0)

  call self%geom%f_comm%allreduce(l_min, stats(1), fckit_mpi_min())
  call self%geom%f_comm%allreduce(l_max, stats(2), fckit_mpi_max())
  call self%geom%f_comm%allreduce(l_sum, stats(3), fckit_mpi_sum())
  call self%geom%f_comm%allreduce(l_count, g_count, fckit_mpi_sum())
  stats(3) = stats(3) / g_count
end subroutine

! ------------------------------------------------------------------------------
! Initialize the grid, allocate memory, for the diffusion operator.
! A call to either calibrate() or read() should be done before
! this class is ready to use.
subroutine soca_diffusion_init(self, geom, f_conf)
  class(soca_diffusion), intent(inout) :: self
  class(soca_geom), target, intent(in) :: geom
  type(fckit_configuration), intent(in) :: f_conf

  type(fckit_configuration), allocatable :: f_conf_list(:)

  real(kind=kind_real) :: stats(3) ! min, max, mean
  character(len=1024) :: str
  character(len=:), allocatable :: str_list(:)
  integer :: i, j

  call oops_log%trace("soca_diffusion::init() starting", flush=.true.)
  call oops_log%info("")
  self%geom => geom

  ! read variable -> group_name mapping
  if (f_conf%get("group mapping", f_conf_list)) then
    call oops_log%info("ExplicitDiffusion: variable groups")    
    allocate(self%group_mapping(size(f_conf_list)))
    do i=1,size(f_conf_list)
      call f_conf_list(i)%get_or_die("name", self%group_mapping(i)%group_name)
      call oops_log%info(" group name: "//self%group_mapping(i)%group_name)
      call oops_log%info("   variables:")
      call f_conf_list(i)%get_or_die("variables", str_list)
      self%group_mapping(i)%variables = oops_variables()
      call self%group_mapping(i)%variables%push_back(str_list)
      do j=1,self%group_mapping(i)%variables%nvars()
        call oops_log%info("     "//self%group_mapping(i)%variables%variable(j))
      end do
    end do
  end if

  ! grid and derived grid parameters
  !---------------------------------------------------------------------------
  call oops_log%info("ExplicitDiffusion: Initializing grid")
  
  allocate(self%mask(DOMAIN_WITH_HALO))
  allocate(self%dx(DOMAIN_WITH_HALO))
  allocate(self%dy(DOMAIN_WITH_HALO))
  ! NOTE MOM6 likes to leave unused halo regions undefined, 
  ! so we explicitly set safe default values here (should probably do this in geom instead??)
  self%mask = 0.0
  self%mask(DOMAIN) = self%geom%mask2d(DOMAIN)
  self%dx = 1.0e-5   
  self%dx(DOMAIN) = geom%dx(DOMAIN)
  self%dy = 1.0e-5
  self%dy(DOMAIN) = geom%dy(DOMAIN)
  call mpp_update_domains(self%mask, self%geom%Domain%mpp_domain, complete=.true.)
  call mpp_update_domains(self%dx, self%geom%Domain%mpp_domain, complete=.true.)
  call mpp_update_domains(self%dy, self%geom%Domain%mpp_domain, complete=.true.)

  call self%calc_stats(self%dx,stats)
    write (str, *) "   dx:   min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)
  call self%calc_stats(self%dy, stats)
    write (str, *) "   dy:   min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)
  call self%calc_stats(self%geom%cell_area, stats)
    write (str, *) "   area: min=", stats(1), "max=", stats(2), "mean=", stats(3)
    call oops_log%info(str)

  ! calculate derived parameters
  allocate(self%pmon_u(DOMAIN_WITH_HALO))
  allocate(self%pnom_v(DOMAIN_WITH_HALO))
  allocate(self%inv_sqrt_area(DOMAIN_WITH_HALO))
  self%pmon_u = 0.0
  self%pnom_v = 0.0
  self%inv_sqrt_area = 0.0

  do j = LOOP_DOMAIN_J
    do i = LOOP_DOMAIN_I
      self%inv_sqrt_area(i,j) = 1.0 / sqrt(self%dx(i,j)*self%dy(i,j))
      self%pmon_u(i,j) = (1.0/self%dx(i-1,j) + 1.0/self%dx(i,j)) / (1.0/self%dy(i-1,j) + 1.0/self%dy(i,j))
      self%pnom_v(i,j) = (1.0/self%dy(i,j-1) + 1.0/self%dy(i,j)) / (1.0/self%dx(i,j-1) + 1.0/self%dx(i,j))
    end do
  end do

  call mpp_update_domains(self%inv_sqrt_area, self%geom%Domain%mpp_domain, complete=.true.)
  call mpp_update_domains(self%pmon_u, self%geom%Domain%mpp_domain, complete=.true.)
  call mpp_update_domains(self%pnom_v, self%geom%Domain%mpp_domain, complete=.true.)   

  call oops_log%trace("soca_diffusion::init() done", flush=.true.)
end subroutine

! ------------------------------------------------------------------------------
! Calibration for the explicit diffuion operator.
! The following are calculated and saved to a file:
! 1) number of iterations required
! 2) the diffusion constants
! 3) the normalization constants
subroutine soca_diffusion_calibrate(self, f_conf)
  class(soca_diffusion),  intent(inout) :: self
  type(fckit_configuration), intent(in) :: f_conf
  
  type(fckit_configuration) :: params_conf, norm_conf
  type(fckit_configuration), allocatable :: group_conf(:)
  character(len=1024) :: str
  integer :: ngroup, grp

  call oops_log%trace("soca_diffusion::calibrate() starting", flush=.true.)
  call oops_log%info("ExplicitDiffusion: running calibration")
  
  ! initialize the groups
  ! TODO move this to init() so we don't have to do the same thing when reading?    
  call f_conf%get_or_die('normalization', norm_conf)
  call f_conf%get_or_die('groups', group_conf)
  ngroup = size(group_conf)
  allocate(self%group(ngroup))
  do grp=1,ngroup
    ! get group name
    call group_conf(grp)%get_or_die("name", self%group(grp)%name)
    write (str, '(A,I2,A,I2)') " group ", grp, " of ", size(self%group)
    call oops_log%info(str)
    write (str, *) " name: ", self%group(grp)%name
    call oops_log%info(str)  
  
    ! horizontal calibration
    if(group_conf(grp)%get("horizontal", params_conf)) then
      ! allocate space and initialize with safe values    
      allocate(self%group(grp)%KhDt(DOMAIN_WITH_HALO))
      allocate(self%group(grp)%normalization_hz(DOMAIN_WITH_HALO))
      self%group(grp)%KhDt = 0.0
      self%group(grp)%normalization_hz = 1.0
      call self%calibrate_hz(self%group(grp), params_conf, norm_conf)
    end if

    ! vertical calibration
    if(group_conf(grp)%get("vertical", params_conf)) then
      allocate(self%group(grp)%KvDt(DOMAIN_WITH_HALO, self%geom%nzo))    
      allocate(self%group(grp)%normalization_vt(DOMAIN_WITH_HALO, self%geom%nzo))
      self%group(grp)%KvDt = 0.0
      self%group(grp)%normalization_vt = 1.0
      call self%calibrate_vt(self%group(grp), params_conf)
    end if
  end do
  call oops_log%trace("soca_diffusion::calibrate() done", flush=.true.)
end subroutine

! ------------------------------------------------------------------------------
subroutine soca_diffusion_calibrate_vt(self, params, vt_conf)
  class(soca_diffusion),              intent(inout) :: self
  class(soca_diffusion_group_params), intent(inout) :: params
  type(fckit_configuration),          intent(in)    :: vt_conf

  real(kind=kind_real) :: stats(3) ! min, max, mean  
  real(kind=kind_real), allocatable :: vt_scales(:,:,:), r_tmp(:,:,:)
  real(kind=kind_real) :: fixed_scale
  integer :: i, j
  character(len=1024) :: str  
  logical :: b

  call oops_log%info("  vertical calibration:")

  allocate(vt_scales(DOMAIN_WITH_HALO, self%geom%nzo))
  allocate(r_tmp(DOMAIN_WITH_HALO, self%geom%nzo))

  ! Get input lengthscales. Either from:
  !  1) a fixed length scale used globally
  !  2) read in from a file
  ! the result is hz_scales containing the length scales (defined as 1 sigma of a guassian)
  vt_scales = 1e10  
  if (.not. vt_conf%has("fixed value") .neqv. vt_conf%has("from file")) then
    ! that was an XOR opperation above, if you were curious
    call abor1_ftn("ERROR: calibration.scales[] must define 1 of 'fixed value' or 'from file'")
  end if
  if ( vt_conf%has("fixed value")) then
    ! used a single fixed value globally
    call oops_log%info("    Using fixed length scales")
    call vt_conf%get_or_die("fixed value", fixed_scale)
    vt_scales = fixed_scale
  else
    call abor1_ftn("ERROR: reading vt scales from file not yet supported")
  end if
  if (.not. vt_conf%get("as gaussian", b)) b = .false.
  write (str,*) "   input values as gaussian (vs GC half width): ", b
  call oops_log%info(str)
  if (.not. b) then
    ! by default, a gaspari cohn half width is expected in the config.
    ! (but the rest of this code asssumes gaussian 1 sigma)
    ! Do the conversion if needed
    vt_scales = vt_scales / 3.57_kind_real
  end if    

  ! TODO make sure we are handling bottom mask
  
  ! print some stats
  call self%calc_stats(vt_scales, stats)
  write (str, '(4X,A,EN10.1,A,EN10.1,A,EN10.1)') &
      "L_vt: min=", stats(1), "  max=", stats(2), "  mean=", stats(3)
  call oops_log%info(str)

  ! calculate the minimum number of iterations needed, rounding up to the
  ! nearest even number.
  !  M >= 2.0 *(L^2)
  r_tmp = 0.0
  do j = LOOP_DOMAIN_J
    do i = LOOP_DOMAIN_I
      if (self%mask(i,j) == 0.0) then
        vt_scales(i,j,:) = 0.0
        cycle
      end if
      r_tmp(i,j,:) = 2.0 * vt_scales(i,j,:)**2
    end do
  end do
  call self%calc_stats(r_tmp, stats)
  i = ceiling(stats(2))
  if (mod(i,2) == 1) i = i + 1
  params%niter_vt = i
  write (str, '(4X,A,I5)') "minimum iterations: ", i
  call oops_log%info(str)    

  ! calculate KvDt based on scales and number of iterations
  params%KvDt = vt_scales**2 / (2.0 * params%niter_vt)

  ! calculate normalization
  call oops_log%info("    Calculating vertical normalization...")
  params%normalization_vt = 1.0
  call self%calibrate_norm_vt(params)
end subroutine

! ------------------------------------------------------------------------------
subroutine soca_diffusion_calibrate_hz(self, params, hz_conf, norm_conf)
  class(soca_diffusion),              intent(inout) :: self
  class(soca_diffusion_group_params), intent(inout) :: params
  type(fckit_configuration),          intent(in)    :: hz_conf, norm_conf

  real(kind=kind_real) :: stats(3) ! min, max, mean  
  character(len=1024) :: str  
  character(len=:), allocatable :: str2, str3
  integer :: i, j, idr
  type(restart_file_type) :: restart_file
  logical :: b

  real(kind=kind_real) :: fixed_scale
  real(kind=kind_real), allocatable :: hz_scales(:,:), r_tmp(:,:)

  call oops_log%info("  horizontal calibration:")
  
  ! allocate things for use later in the loops
  allocate(hz_scales(DOMAIN_WITH_HALO))
  allocate(r_tmp(DOMAIN_WITH_HALO))
    
  ! Get input lengthscales. Either from:
  !  1) a fixed length scale used globally
  !  2) read in from a file
  ! the result is hz_scales containing the length scales (defined as 1 sigma of a guassian)
  hz_scales = 1e10  
  if (.not. hz_conf%has("fixed value") .neqv. hz_conf%has("from file")) then
    ! that was an XOR opperation above, if you were curious
    call abor1_ftn("ERROR: calibration.scales[] must define 1 of 'fixed value' or 'from file'")
  end if
  if ( hz_conf%has("fixed value")) then
    ! used a single fixed value globally
    call oops_log%info("    Using fixed length scales")
    call hz_conf%get_or_die("fixed value", fixed_scale)
    hz_scales = fixed_scale
  else
    ! read lengths from a file. a 2d field is expected
    call oops_log%info("    Reading length scales from file")
    call hz_conf%get_or_die("from file.filename", str2)
    call hz_conf%get_or_die("from file.variable name", str3)
    call fms_io_init()
    idr = register_restart_field(restart_file, str2, str3, &
    hz_scales, domain=self%geom%Domain%mpp_domain)
    call restore_state(restart_file, directory='')
    call free_restart_type(restart_file)
    call fms_io_exit()
  end if
  if(.not. hz_conf%get("as gaussian", b)) b = .false.  
  write (str,*) "   input values as gaussian (vs GC half width): ", b
  call oops_log%info(str)
  if (.not. b) then
    ! by default, a gaspari cohn half width is expected in the config.
    ! (but the rest of this code asssumes gaussian 1 sigma)
    ! Do the conversion if needed
    hz_scales = hz_scales / 3.57_kind_real
  end if

  ! make sure halos are up to date
  call mpp_update_domains(hz_scales, self%geom%Domain%mpp_domain, complete=.true.)

  ! print some stats
  call self%calc_stats(hz_scales, stats)
  write (str, '(4X,A,EN10.1,A,EN10.1,A,EN10.1)') &
      "L_hz: min=", stats(1), "  max=", stats(2), "  mean=", stats(3)
  call oops_log%info(str)

  ! calculate the minimum number of iterations needed, rounding up to the
  ! nearest even number.
  !  M >= 2.0 *(L^2 / (1/dx^2 + 1/dy^2))
  r_tmp = 0.0
  do j = LOOP_DOMAIN_J
    do i = LOOP_DOMAIN_I
      if (self%mask(i,j) == 0.0 ) cycle
      r_tmp(i,j) = 2.0 * hz_scales(i,j)**2 * (1.0/(self%dx(i,j)**2) + 1.0/(self%dy(i,j)**2))
    end do
  end do
  call self%calc_stats(r_tmp, stats)
  i = ceiling(stats(2))
  if (mod(i,2) == 1) i = i + 1
  params%niter_hz = i
  write (str, '(4X,A,I5)') "minimum iterations: ", i
  call oops_log%info(str)

  ! calculate KhDt based on scales and number of iterations    
  params%KhDt = hz_scales**2 / (2.0 * params%niter_hz)

  ! calculate normalization
  call oops_log%info("    Calculating horizontal normalization:")
  call norm_conf%get_or_die("method", str2)
  call oops_log%info("      method: "//str2)
  params%normalization_hz = 1.0
  if (str2 == "brute force") then
    call self%calibrate_norm_hz_bruteforce(params)
  else if (str2 == "randomization") then
    call norm_conf%get_or_die("iterations", i)
    call self%calibrate_norm_hz_randomization(i, params)
  else
    call abor1_ftn("ERROR: normalization.method must be 'brute force' or 'randomization'")
  end if
end subroutine

! ------------------------------------------------------------------------------
! Perform horizontal diffusion on each level
subroutine soca_diffusion_multiply(self, dx)
  class(soca_diffusion), intent(inout) :: self
  type(soca_increment),  intent(inout) :: dx

  real(kind=kind_real), allocatable :: tmp2d(:,:)
  character(len=1024) :: str
  integer :: f, z, grp, g, g2

  call oops_log%trace("soca_diffusion::multiply() starting", flush=.true.)

  if (.not. allocated(self%group) .or. .not. allocated(self%group_mapping)) then
    ! uninitialized, calibrate or read should have been called before now
    call abor1_ftn("ERROR: soca_diffusion has not been initialized.")
  end if
  
  allocate(tmp2d(DOMAIN_WITH_HALO))
  do f=1, size(dx%fields)
    ! find the group associated with this variable
    ! TODO, simplify this
    call oops_log%debug("running multiply on "//dx%fields(f)%name)
    grp = -1
    do g=1, size(self%group_mapping)
      if (self%group_mapping(g)%variables%has(dx%fields(f)%name)) then
        do g2=1, size(self%group)
          if (self%group(g2)%name == self%group_mapping(g)%group_name) then
            grp = g2
          end if
        end do
      end if
    end do
    if (grp == -1) then
      call abor1_ftn("ERROR: could not find a valid group for the variable "//dx%fields(f)%name)
    end if

    ! normalization (horizontal + vertical)
    do z = 1, dx%fields(f)%nz
      dx%fields(f)%val(DOMAIN,z)  = dx%fields(f)%val(DOMAIN,z) * self%group(grp)%normalization_hz(DOMAIN)
    end do
    if (self%group(grp)%niter_vt > 0) then
      dx%fields(f)%val(DOMAIN,:)  = dx%fields(f)%val(DOMAIN,:) * self%group(grp)%normalization_vt(DOMAIN,:)
    end if  
     
    ! vertical diffusion AD
    if (self%group(grp)%niter_vt > 0) then
      call self%diffusion_vt_ad(dx%fields(f)%val, self%group(grp))
    end if   

    do z = 1, dx%fields(f)%nz
      tmp2d = dx%fields(f)%val(:,:,z)
      
      ! horizontal diffusion AD      
      call self%diffusion_hz_ad(tmp2d, self%group(grp))

      ! TODO grid metric
      ! tmp2d = tmp2d * self%inv_sqrt_area

      ! horizontal diffusion TL
      call self%diffusion_hz_tl(tmp2d, self%group(grp))

      dx%fields(f)%val(DOMAIN,z) = tmp2d(DOMAIN)
    end do
  
    ! vertical diffusion TL    
    if (self%group(grp)%niter_vt > 0) then
      call self%diffusion_vt_tl(dx%fields(f)%val, self%group(grp))
    end if
    
    ! normalization (horizontal + vertical)
    if (self%group(grp)%niter_vt > 0) then    
      dx%fields(f)%val(DOMAIN,:)  = dx%fields(f)%val(DOMAIN,:) * self%group(grp)%normalization_vt(DOMAIN,:)    
    end if
    do z = 1, dx%fields(f)%nz
      dx%fields(f)%val(DOMAIN,z)  = dx%fields(f)%val(DOMAIN,z) * self%group(grp)%normalization_hz(DOMAIN)
    end do

  end do

  call oops_log%trace("soca_diffusion::multiply() done", flush=.true.)
end subroutine

! ------------------------------------------------------------------------------
! Apply half the required iterations of diffusion
subroutine soca_diffusion_hz_tl(self, field, params)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable, intent(inout) :: field(:,:)
  type(soca_diffusion_group_params), intent(in) :: params

  real(kind=kind_real), allocatable :: flux_x(:,:), flux_y(:,:), hfac(:,:)
  integer :: i, j, iter, niter
 
  ! note, number of iterations is half of what is required
  ! (the other half come from application of adjoint)
  niter = params%niter_hz / 2

  ! NOTE: flux_x(i,j) is the flux through the western edge of the grid cell.
  !  this is opposite of MOM6 conventions where u(i,j) points are east of t(i,j) points.
  !  (dosen't really matter, just something to note)
  allocate(flux_x(DOMAIN_WITH_HALO))
  allocate(flux_y(DOMAIN_WITH_HALO))
  flux_x = 0.0
  flux_y = 0.0

  ! calculate some needed constants  (TODO, move to initialization?)
  allocate(hfac(DOMAIN_WITH_HALO))
  hfac = 0.0
  do j=LOOP_DOMAIN_J
    do i=LOOP_DOMAIN_I
      hfac(i,j) = (1.0/self%dy(i,j)/self%dx(i,j))
    end do
  end do
  call mpp_update_domains(field, self%geom%Domain%mpp_domain, complete=.true.)

  do iter=1,niter
      ! calculate diffusive flux on each edge of a grid box. masking out where there is land
    do j=LOOP_DOMAIN_J 
      do i=LOOP_DOMAIN_I+1 ! assume halo size is >= 1, and skip doing a halo update
        flux_x(i,j) = self%pmon_u(i,j) * 0.5 * (params%KhDt(i,j) + params%KhDt(i-1,j)) * (field(i,j) - field(i-1,j))
        flux_x(i,j) = flux_x(i,j) * self%mask(i,j) * self%mask(i-1,j)
      end do
    end do
    do j=LOOP_DOMAIN_J+1 ! assume halo size is >= 1, and skip doing a halo update
      do i=LOOP_DOMAIN_I
        flux_y(i,j) = self%pnom_v(i,j) * 0.5 * (params%KhDt(i,j) + params%KhDt(i,j-1)) * (field(i,j) - field(i,j-1))
        flux_y(i,j) = flux_y(i,j) * self%mask(i,j) * self%mask(i,j-1)
      end do
    end do

    ! time-step hz diffusion terms
    do j=LOOP_DOMAIN_J
      do i=LOOP_DOMAIN_I
        field(i,j) = field(i,j) + hfac(i,j) * &
          (flux_x(i+1, j) - flux_x(i,j) + flux_y(i, j+1) - flux_y(i,j))
      end do
    end do
    call mpp_update_domains(field, self%geom%Domain%mpp_domain, complete=.true.)
  end do
end subroutine


! ------------------------------------------------------------------------------
subroutine soca_diffusion_hz_ad(self, field, params)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable, intent(inout) :: field(:,:) 
  type(soca_diffusion_group_params), intent(in) :: params

  real(kind=kind_real), allocatable :: wrk_old(:,:), wrk_new(:,:), tmp(:,:)
  real(kind=kind_real), allocatable :: flux_x(:,:), flux_y(:,:), hfac(:,:)
  real(kind=kind_real) :: adfac
  integer :: i, j, iter, niter
  
  ! note, number of iterations is half of what is required
  ! (the other half come from application of tl)
  niter = params%niter_hz / 2

  ! NOTE: flux_x(i,j) is the flux through the western edge of the grid cell.
  !  this is opposite of MOM6 conventions where u(i,j) points are east of t(i,j) points.
  !  (dosen't really matter, just something to note)
  allocate(wrk_new(DOMAIN_WITH_HALO))
  allocate(wrk_old(DOMAIN_WITH_HALO))
  allocate(tmp(DOMAIN_WITH_HALO))
  allocate(flux_x(DOMAIN_WITH_HALO))
  allocate(flux_y(DOMAIN_WITH_HALO))
  wrk_old = 0.0
  wrk_new = 0.0
  flux_x = 0.0
  flux_y = 0.0

  ! calculate some needed constants (TODO, move to initialization?)
  allocate(hfac(DOMAIN_WITH_HALO))
  hfac = 0.0
  do j=LOOP_DOMAIN_J
    do i=LOOP_DOMAIN_I
      hfac(i,j) = (1.0/self%dy(i,j)/self%dx(i,j))
    end do
  end do  

  ! adjoint of convoled solution  
  do j=LOOP_DOMAIN_J
    do i=LOOP_DOMAIN_I
      wrk_old(i,j) = wrk_old(i,j) + field(i,j)
      field(i,j) = 0.0
    end do
  end do

  ! integrate adjoint hz diffusion terms
  do iter=1,niter
    tmp = wrk_new
    wrk_new = wrk_old
    wrk_old = tmp    
    flux_x = 0.0
    flux_y = 0.0

    ! time-step adjoint hz diffusion terms
    do j=LOOP_DOMAIN_J
      do i=LOOP_DOMAIN_I
        adfac=hfac(i,j)*wrk_new(i,j)
        flux_y(i,j  ) = flux_y(i,j  ) - adfac
        flux_y(i,j+1) = flux_y(i,j+1) + adfac
        flux_x(i  ,j) = flux_x(i  ,j) - adfac
        flux_x(i+1,j) = flux_x(i+1,j) + adfac
        wrk_old(i,j) = wrk_new(i,j)
      end do
    end do
    wrk_new = 0.0

    ! compute adjoint diffusive flux
    do j=LOOP_DOMAIN_J+1
      do i=LOOP_DOMAIN_I
        flux_y(i,j) = flux_y(i,j) * self%mask(i,j) * self%mask(i,j-1)
        adfac = self%pnom_v(i,j) * 0.5*(params%KhDt(i,j-1)+params%KhDt(i,j)) * flux_y(i,j)
        wrk_old(i,j-1) = wrk_old(i,j-1) - adfac
        wrk_old(i,j  ) = wrk_old(i,j  ) + adfac
        flux_y(i,j) = 0.0
      end do
    end do
    do j=LOOP_DOMAIN_J
      do i=LOOP_DOMAIN_I+1
        flux_x(i,j) = flux_x(i,j) * self%mask(i,j) * self%mask(i-1,j)
        adfac = self%pmon_u(i,j) * 0.5*(params%KhDt(i-1,j)+params%KhDt(i,j)) * flux_x(i,j)
        wrk_old(i-1,j) = wrk_old(i-1,j) - adfac
        wrk_old(i,  j) = wrk_old(i,  j) + adfac
        flux_x(i,j) = 0.0
      end do
    end do
    call mpp_update_domains_ad(wrk_old, self%geom%Domain%mpp_domain, complete=.true.)
  end do

  ! set adjoint initial conditions
  field = field + wrk_old
  wrk_old = 0.0

end subroutine

! ------------------------------------------------------------------------------
subroutine soca_diffusion_vt_tl(self, field, params)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable, intent(inout) :: field(:,:,:) 
  type(soca_diffusion_group_params), intent(in) :: params

  real(kind=kind_real), allocatable :: flux(:)
  real(kind=kind_real) :: vfac
  integer :: i, j, k, iter, niter, nz

  ! TODO grid metric?

  ! apply M/2 iterations of diffusion
  niter = params%niter_vt / 2
  nz = size(field, dim=3)
  vfac = 1.0
  allocate(flux(nz+1))

  do j=LOOP_DOMAIN_J
    do i=LOOP_DOMAIN_I
      if(self%mask(i,j) == 0.0) then
        cycle
      end if

      do iter=1, niter        
        ! calculate diffusive flux at the edgescall mpp_update_domains(field, self%geom%Domain%mpp_domain, complete=.true.)
        flux = 0.0
        do k=2,nz
          flux(k) = 0.5 * (params%KvDt(i,j,k) + params%KvDt(i,j,k-1)) * (field(i,j,k)-field(i,j,k-1))
        end do

        ! time-step vt diffusion terms  
        do k=1,nz
          field(i,j,k) = field(i,j,k) + vfac * (flux(k+1) - flux(k))
        end do
      end do
    end do
  end do
end subroutine

! ------------------------------------------------------------------------------
subroutine soca_diffusion_vt_ad(self, field, params)
  class(soca_diffusion), intent(inout) :: self
  real(kind=kind_real), allocatable, intent(inout) :: field(:,:,:) 
  type(soca_diffusion_group_params), intent(in) :: params

  real(kind=kind_real), allocatable :: flux(:), wrk_old(:), wrk_new(:), tmp(:)
  real(kind=kind_real) :: adfac, vfac
  integer :: i, j, k, iter, niter, nz


  vfac = 1.0 ! TODO check this
  niter = params%niter_vt / 2
  nz = size(field, dim=3)
  allocate(flux(nz+1))
  allocate(wrk_old(nz), wrk_new(nz), tmp(nz))

  do j=LOOP_DOMAIN_J
    do i=LOOP_DOMAIN_I
      if(self%mask(i,j) == 0.0) cycle
      wrk_old = 0.0
      wrk_new = 0.0
      flux = 0.0
      
      wrk_old(:) = field(i,j,:)
      field(i,j,:) = 0.0

      ! integrate adjoint vt diffusion terms
      do iter=1,niter
        tmp(:) = wrk_new(:)
        wrk_new(:) = wrk_old(:)
        wrk_old(:) = tmp(:)
        flux(:) = 0.0
        
        ! timestep adjoint vt diffusion term
        do k=1,nz
          adfac=vfac*wrk_new(k)
          flux(k)   = flux(k)   - adfac
          flux(k+1) = flux(k+1) + adfac
          wrk_old(k) = wrk_new(k)
        end do
        wrk_new = 0.0

        ! adjoint diffusive flux
        flux(1) = 0.0
        do k=2,nz
          adfac = 0.5*(params%KvDt(i,j,k-1) + params%KvDt(i,j,k)) * flux(k)
          wrk_old(k-1) = wrk_old(k-1) - adfac
          wrk_old(k)   = wrk_old(k)   + adfac
          flux(k) = 0.0
        end do
      end do
      
      field(i,j,:) = field(i,j,:) + wrk_old(:)
    end do
  end do
end subroutine

! ------------------------------------------------------------------------------
! Calculate the exact normalization weights using the brute force method
! (creating a dirac at every SINGLE point).
! You probably don't want to use this, except for testing. 
! Use randomization instead.
! ------------------------------------------------------------------------------
subroutine soca_diffusion_calibrate_norm_hz_bruteforce(self, params)
  class(soca_diffusion), intent(inout) :: self
  type(soca_diffusion_group_params), intent(inout) :: params

  integer :: i, j, n, n10pct
  character(len=1024) :: str  
  logical :: local
  real(kind=kind_real), allocatable :: r_tmp(:,:), norm(:,:)

  call oops_log%info("WARNING: make sure you really want to be using bruteforce!")

  allocate(r_tmp(DOMAIN_WITH_HALO))
  allocate(norm(DOMAIN_WITH_HALO))

  n=1
  n10pct = (self%geom%jeg-self%geom%jsg+1) * (self%geom%ieg-self%geom%isg+1) / 10
  params%normalization_hz = 1.0
  call oops_log%info("      normalization progress: ")
  do j=self%geom%jsg, self%geom%jeg
    do i=self%geom%isg, self%geom%ieg

      if (mod(n, n10pct) == 0) then
        write (str, '(8X,I3,A)') 10*n/n10pct, "%"
        ! hmm, odd, it doesn't flush if I set newl=.false.
        call oops_log%info(str)
      end if
      n = n + 1

      r_tmp = 0.0
      local = i >= self%geom%isc .and. i <= self%geom%iec .and. &
              j >= self%geom%jsc .and. j <= self%geom%jec
      if (local) r_tmp(i,j) = 1.0
      
      call self%diffusion_hz_ad(r_tmp, params)
      ! TODO grid metric
      call self%diffusion_hz_tl(r_tmp, params)
      
      if (local) then
        if(self%mask(i,j) == 0.0) cycle
        norm(i,j) = 1.0 / sqrt(r_tmp(i,j))
      end if
    end do
  end do
  call mpp_update_domains(norm, self%geom%Domain%mpp_domain, complete=.true.)
  params%normalization_hz = norm
end subroutine

! ------------------------------------------------------------------------------
subroutine soca_diffusion_calibrate_norm_vt(self, params)
  class(soca_diffusion), intent(inout) :: self
  type(soca_diffusion_group_params), intent(inout) :: params

  integer :: k, nz
  real(kind=kind_real), allocatable :: r_tmp(:,:,:), norm(:,:,:)

  nz = self%geom%nzo
  allocate(r_tmp(DOMAIN_WITH_HALO, nz))
  allocate(norm(DOMAIN_WITH_HALO, nz))

  norm = 1.0
  do k=1, nz
    r_tmp = 0.0
    r_tmp(:,:,k) = 1.0
    call self%diffusion_vt_ad(r_tmp, params)
    call self%diffusion_vt_tl(r_tmp, params)
    norm(:,:,k) = 1.0 / sqrt(r_tmp(:,:, k))
  end do
  
  params%normalization_vt = norm
end subroutine

! ------------------------------------------------------------------------------
! Estimate the normalization weights by creating random vectors (normally distributed)
! applying the diffusion TL, and keeping a running statistic of the variance of
! those results.
!
! Typically a good number of iterations is around 10,000
! ------------------------------------------------------------------------------
subroutine soca_diffusion_calibrate_norm_hz_randomization(self, iter, params)
  class(soca_diffusion), intent(inout) :: self
  integer, intent(in) :: iter
  type(soca_diffusion_group_params), intent(inout) :: params

  real(kind=kind_real), allocatable :: field(:,:)
  real(kind=kind_real), allocatable :: s(:,:)
  real(kind=kind_real), allocatable :: m(:,:), new_m(:,:)

  integer :: n, n10pct, rnd
  character(len=1024) :: str  
  
  allocate(field(DOMAIN_WITH_HALO))
  allocate(s(DOMAIN_WITH_HALO))
  allocate(m(DOMAIN_WITH_HALO))
  allocate(new_m(DOMAIN_WITH_HALO))

  s = 0.0
  m = 0.0
  n10pct = iter/10 !< ouput info to screen every 10% 

  call oops_log%info("      normalization progress: ")
  do n=1,iter
    if (mod(n, n10pct) == 0) then
      write (str, '(8X,I3,A)') 10*n/n10pct, "%"
      ! hmm, odd, it doesn't flush if I set newl=.false.
      call oops_log%info(str)
    end if
  
    ! create a random vector
    ! Ensure random number are different on each PE.
    ! TODO: when this is all refactored into saber, it would 
    !   be nice to generate the random numbers on 1 PE then scatter,
    !   this will ensure answers don't change when PE counts change
    rnd=n*self%geom%f_comm%size() + self%geom%f_comm%rank()
    call normal_distribution(field, 0.0_kind_real, 1.0_kind_real, rnd, .true.) 

    ! grid metric
    ! TODO

    ! apply the diffusion TL
    call self%diffusion_hz_tl(field, params)

    ! keep track of the stats needed for a running variance calculation
    ! (Welford 1962 algorithm)
    new_m = m + (field-m)/n
    s = s + (field - m)*(field - new_m)
    m = new_m
  end do
  
  ! calculate final variance
  field = (s/(iter-1)) 

  ! normalization (where ocean) is 1/sqrt(variance)
  where (self%mask == 1.0)  params%normalization_hz = 1.0 / sqrt(field)
  
  call mpp_update_domains(params%normalization_hz, self%geom%Domain%mpp_domain, complete=.true.)
end subroutine

! ------------------------------------------------------------------------------
! write out the parameters to a restart file
subroutine soca_diffusion_write_params(self, f_conf)
  class(soca_diffusion), intent(inout) :: self
  type(fckit_configuration), intent(in) :: f_conf  
  
  type(fckit_configuration), allocatable :: f_conf_list(:)
  character(len=:), allocatable   :: filename, group_name
  type(restart_file_type) :: restart_file
  character(len=1024) :: str
  integer :: idr, grp

  call f_conf%get_or_die("groups", f_conf_list)

  ! write to file
  call fms_io_init()
  do grp=1,size(f_conf_list)
    call f_conf_list(grp)%get_or_die('name', group_name)
    call f_conf_list(grp)%get_or_die('write.filename', filename)

    write (str, *) "Writing group: " // group_name // " to '" // filename // "'"
    call oops_log%info(str)
   
    ! write horizontal parameters
    if (self%group(grp)%niter_hz > 0) then
      str = "iterations_hz"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%niter_hz, domain=self%geom%Domain%mpp_domain)

      str = "khdt"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%KhDt, domain=self%geom%Domain%mpp_domain)

      str = "normalization_hz"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%normalization_hz, domain=self%geom%Domain%mpp_domain)  
    end if
    
    ! write vertical parameters
    if (self%group(grp)%niter_vt > 0) then
      str = "iterations_vt"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%niter_vt, domain=self%geom%Domain%mpp_domain)

      str = "kvdt"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%KvDt, domain=self%geom%Domain%mpp_domain)  

      str = "normalization_vt"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%normalization_vt, domain=self%geom%Domain%mpp_domain)        
    end if  

    call save_restart(restart_file, directory='')
    call free_restart_type(restart_file)
  end do 
  call fms_io_exit()  
end subroutine

! ------------------------------------------------------------------------------
! read in the parameters from a restart file
subroutine soca_diffusion_read_params(self, f_conf)
  class(soca_diffusion), intent(inout) :: self
  type(fckit_configuration), intent(in) :: f_conf

  type(fckit_configuration), allocatable :: f_conf_list(:)
  type(fckit_configuration) :: params_conf
  character(len=:), allocatable   :: filename, group_name
  type(restart_file_type) :: restart_file
  integer :: idr, grp
  character(len=1024) :: str

  call f_conf%get_or_die("groups", f_conf_list)

  ! make sure we havent read in parameters already
  if ( allocated(self%group)) then
    call abor1_ftn("ERROR: soca_diffusion has already been initialized.")
  end if
  allocate(self%group(size(f_conf_list)))
  
  ! read from file
  call fms_io_init()
  do grp=1,size(self%group)
    call f_conf_list(grp)%get_or_die('name', group_name)
    self%group(grp)%name = trim(group_name)    
    write (str, *) "Reading group: " // group_name
    call oops_log%info("")
    call oops_log%info(str)

    ! read horizontal
    if (f_conf_list(grp)%get('horizontal', params_conf)) then
      call params_conf%get_or_die('filename', filename)
      call oops_log%info("  horizontal parameters:")
      call oops_log%info("    filename: "//filename)
    
      allocate(self%group(grp)%KhDt(DOMAIN_WITH_HALO))
      allocate(self%group(grp)%normalization_hz(DOMAIN_WITH_HALO))

      str = "iterations_hz"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%niter_hz, domain=self%geom%Domain%mpp_domain)

      str = "khdt"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%KhDt, domain=self%geom%Domain%mpp_domain)

      str = "normalization_hz"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%normalization_hz, domain=self%geom%Domain%mpp_domain)
    
      call restore_state(restart_file, directory='')
      call free_restart_type(restart_file)

      write (str, '(4X,A,I5)') "minimum iterations: ", self%group(grp)%niter_hz
      call oops_log%info(str)

      call mpp_update_domains(self%group(grp)%normalization_hz, self%geom%Domain%mpp_domain, complete=.true.)
      call mpp_update_domains(self%group(grp)%KhDt, self%geom%Domain%mpp_domain, complete=.true.)      
    end if

    ! read vertical
    if (f_conf_list(grp)%get('vertical', params_conf)) then
      call params_conf%get_or_die('filename', filename)
      call oops_log%info("  vertical parameters:")
      call oops_log%info("    filename: "//filename)

      allocate(self%group(grp)%KvDt(DOMAIN_WITH_HALO, self%geom%nzo))          
      allocate(self%group(grp)%normalization_vt(DOMAIN_WITH_HALO, self%geom%nzo))

      str = "iterations_vt"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%niter_vt, domain=self%geom%Domain%mpp_domain)

      str = "kvdt"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%KvDt, domain=self%geom%Domain%mpp_domain)    

      str = "normalization_vt"
      idr = register_restart_field(restart_file, filename, str, &
        self%group(grp)%normalization_vt, domain=self%geom%Domain%mpp_domain)          

      call restore_state(restart_file, directory='')
      call free_restart_type(restart_file)  

      write (str, '(4X,A,I5)') "minimum iterations: ", self%group(grp)%niter_vt
      call oops_log%info(str)
    end if
  end do
  call fms_io_exit()
end subroutine

! ------------------------------------------------------------------------------

end module soca_diffusion_mod
