module mom5cice5_thermo

  use mom5cice5_constants
  use mom5cice5_fields
  !use mpi
  implicit none

  
contains

  function Tm(S)
    implicit none
    real(kind=kind_real):: S
    real(kind=kind_real):: Tm

    Tm = -mu * S

    return
  end function Tm

  ! ------------------------------------------------------------------------------

  function Qi_nl(Ti,Si)
    ! Enthalpy of sea-ice given sea-ice temperature and salinity
    !
    !  Args:
    !      Ti   (kind_real): Sea-Ice temperature
    !      Si   (kind_real): Sea-Ice salinity
    !
    !  Return: 
    !         (kind_real): Enthalpy of Sea-Ice     
    !
    real(kind=kind_real) :: Si
    real(kind=kind_real) :: Ti
    real(kind=kind_real) :: Qi_nl
    
    real(kind=kind_real) :: dum1, dum2, dum3
    
    dum1 = c0 * (Tm(Si)-Ti)
    dum2 = L0 * (1-Tm(Si)/Ti)
    dum3 = cw * Tm(Si)

    Qi_nl = -rho_i*(dum1+dum2-dum3)

  end function Qi_nl

  ! ------------------------------------------------------------------------------

  function dQi_tl(dTi, dSi, Ti,Si)
    ! Sea-ice enthalpy increment given sea-ice temperature and salinity increment
    ! and background
    !dQ    μ⋅ρᵢ⋅(-L₀ + 2⋅T⋅c₀)
    !── =  ───────────────────
    !dS            T         
    !
    !  Args:
    !      Ti   (kind_real): Sea-Ice temperature
    !      Si   (kind_real): Sea-Ice salinity
    !
    !  Return: 
    !         (kind_real): Enthalpy of Sea-Ice     
    !
    real(kind=kind_real) :: Si
    real(kind=kind_real) :: Ti
    real(kind=kind_real) :: dSi
    real(kind=kind_real) :: dTi
    real(kind=kind_real) :: dQi_tl

    real(kind=kind_real) :: dqdt, dqds

    dqdt = L0*mu*rho_i*(Si/Ti**2)
    dqds = mu*rho_i*(2.0_kind_Real*c0*Ti-L0)/Ti

    dQi_tl = dqdt * dTi + dqds * dSi

  end function dQi_tl

  ! ------------------------------------------------------------------------------

  function Ti_nl(Qi,Si)
    ! Temperature of sea-ice given sea-ice enthalpy and salinity ( Qi_nl^-1 )
    !
    !  Args:
    !      Qi   (kind_real): Sea-Ice enthalpy
    !      Si   (kind_real): Sea-Ice salinity
    !
    !  Return: 
    !         (kind_real): Temperature of Sea-Ice     
    !
    real(kind=kind_real) :: Si
    real(kind=kind_real) :: Qi
    real(kind=kind_real) :: Ti_nl

    real(kind=kind_real):: a, b, c
    
    a = c0
    b = (cw-c0) * Tm(Si) - (Qi/rho_i) - L0
    c = L0 * Tm(Si)
    
    Ti_nl = -(b+(b**2.0_kind_Real-4.0_kind_real*a*c)**0.5_kind_real)/(2.0_kind_real*a)

  end function Ti_nl

  ! ------------------------------------------------------------------------------

  function dTi_tl(dQi, dSi, Qi,Si)
    ! Sea-ice temperature increment given sea-ice enthalpy and salinity increment 
    ! and background
    !  Args:
    !      dQi   (kind_real): Sea-Ice enthalpy increment
    !      dSi   (kind_real): Sea-Ice salinity increment
    !      Qi    (kind_real): Background Sea-Ice enthalpy
    !      Si    (kind_real): Background Sea-Ice salinity
    !
    !  Return: 
    !         (kind_real): Temperature of Sea-Ice     
    !
    real(kind=kind_real) :: Si
    real(kind=kind_real) :: Qi
    real(kind=kind_real) :: dSi
    real(kind=kind_real) :: dQi
    real(kind=kind_real) :: dTi_tl

    real(kind=kind_real):: a, b, c
    
    a = c0
    b = (cw-c0) * Tm(Si) - (qi/rho_i) - L0
    c = L0 * Tm(Si)
    
    dTi_tl = -(b+(b**2.0_kind_Real-4.0_kind_real*a*c)**0.5_kind_real)/(2.0_kind_real*a)

  end function dTi_tl

  ! ------------------------------------------------------------------------------

  function Ki_nl(Ti,Si)
    ! Thermal conductivity of sea-ice given sea-ice temperature and salinity (conduct='MU71')
    !
    !  Args:
    !      Ti   (kind_real): Sea-Ice temperature
    !      Si   (kind_real): Sea-Ice salinity
    !
    !  Return: 
    !         (kind_real): Thermal conductivity of Sea-Ice     
    !
    real(kind=kind_real):: Si
    real(kind=kind_real):: Ti
    real(kind=kind_real) :: Ki_nl

    Ki_nl = K0 + beta * (Si/Ti)

  end function Ki_nl

end module mom5cice5_thermo
