
!> Structure holding configuration variables for the  model

module mom5cice5_configs

  use kinds
  implicit none
  private
  public :: mom5cice5_config
  public :: mom5cice5_config_registry

  !> Fortran derived type to hold configuration data for the  model
  type :: mom5cice5_config
     integer :: nx     !< Zonal grid dimension
     integer :: ny     !< Meridional grid dimension
     ! dimensional parameters
     real(kind=kind_real) :: dt0       !< dimensional time (seconds)
  end type mom5cice5_config

#define LISTED_TYPE mom5cice5_config

  !> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

  !> Global registry
  type(registry_t) :: mom5cice5_config_registry

  ! ------------------------------------------------------------------------------
contains
  ! ------------------------------------------------------------------------------
  !> Linked list implementation
#include "util/linkedList_c.f"

end module mom5cice5_configs
