
!> Constants for the MOM5CICE5 model

module mom5cice5_constants

  use kinds
  implicit none

  real(kind=kind_real),parameter :: rho_s=330.0_kind_real       !< Snow density [kg/m3]
  real(kind=kind_real),parameter :: rho_i=917.0_kind_real       !< Sea-Ice density [kg/m3]
  real(kind=kind_real),parameter :: c0 = 2106.0_kind_real       !< Specific heat of fresh ice at 0C [J/kg/deg]
  real(kind=kind_real),parameter :: L0 = 3.34e5_kind_real       !< Latent heat of fusion of fresh ice at 0C [J/kg]
  real(kind=kind_real),parameter :: cw = 4218.0_kind_real       !< Specific heat of sea-water ???? <--- Check
  real(kind=kind_real),parameter :: mu = 0.054_kind_real        !< Liquidus ratio between the freezing temperature
                                                                !< and salinity of brine [deg/psu]
  real(kind=kind_real),parameter :: K0 = 2.03_kind_real         !< Conductivity of fresh ice [W/m/deg]
  real(kind=kind_real),parameter :: beta = 0.13_kind_real       !< Empirical constant [W/m/psu]
  real(kind=kind_real),parameter :: Ks = 0.30_kind_real         !< Conductivity of snow [W/m/deg]
  real(kind=kind_real),parameter :: hs_min = 1.0e-4_kind_real   !< Minimum snow depth [m]
  real(kind=kind_real),parameter :: hi_min = 0.01_kind_real     !< Minimum sea-ice thickness [m]

end module mom5cice5_constants
