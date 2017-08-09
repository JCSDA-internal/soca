module mom5cice5_thermo

  use mom5cice5_constants
  use mom5cice5_fields

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

  function Ti_tl(dQi, dSi, Qi,Si)
    ! Temperature of sea-ice given sea-ice enthalpy and salinity ( Qi_nl^-1 )
    !
    !  Args:
    !      dQi   (kind_real): Background Sea-Ice enthalpy
    !      dSi   (kind_real): Background Sea-Ice salinity
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
    real(kind=kind_real) :: Ti_tl

    real(kind=kind_real):: a, b, c
    
    a = c0
    b = (cw-c0) * Tm(Si) - (qi/rho_i) - L0
    c = L0 * Tm(Si)
    
    Ti_tl = -(b+(b**2.0_kind_Real-4.0_kind_real*a*c)**0.5_kind_real)/(2.0_kind_real*a)

  end function Ti_tl

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

  ! ------------------------------------------------------------------------------

  subroutine Kop_inv(dx, xb, dxi)
    implicit none
    type(mom5cice5_field), intent(in) :: dx
    type(mom5cice5_field), intent(in) :: xb
    type(mom5cice5_field), intent(out) :: dxi

    real(kind=kind_real), allocatable :: aice(:,:)
    real(kind=kind_real) :: A, B, C
    
    allocate( aice(xb%nx,xb%ny) )
    
    aice = sum(xb%cicen,3)

    A = (rho_s*c0)**(-1.0_kind_real)
    B = 0.0_kind_real !-mu*rho_i(L0/xb%
    C = 0.0_kind_real

    call check_resolution(dxi, xb)
    
    !A = (rho_s*c0)**(-1.0_kind_real)
    !B = 0.0_kind_real !-mu*rho_i(L0/xb%
    !C = 0.0_kind_real
    !D = -mu*sum(xb%cicen,3)
    !E = -mu*xb%sssoc
    !F = 1.0_kind_real-sum(xb%cicen,3)
    
    call zeros(dxi)
    dxi%cicen = dx%cicen
    dxi%hicen = dx%hicen
    dxi%vicen = xb%cicen * dx%hicen + xb%hicen * dx%cicen
    dxi%hsnon = dx%hsnon
    dxi%vsnon = xb%cicen * dx%hsnon + xb%hsnon * dx%cicen
    dxi%tsfcn = dx%tsfcn
    dxi%qsnon = A*dx%qsnon    ! snow temperature from enthalpy    
    dxi%sicnk = dx%sicnk
    dxi%qicnk = (dx%qicnk - B * dx%sicnk)/C
    dxi%sssoc = dx%sssoc
    dxi%tlioc = dx%tlioc
    dxi%sstoc = -mu*aice * dx%sssoc -mu*xb%sssoc * sum(dx%cicen,3) + (1.0_kind_real-aice) * dx%tlioc  

    return
  end subroutine Kop_inv

  ! ------------------------------------------------------------------------------

  subroutine Kop_inv_ad(dx_ad, xb, dxi_ad)
    implicit none
    type(mom5cice5_field), intent(in) :: dx_ad
    type(mom5cice5_field), intent(in) :: xb
    type(mom5cice5_field), intent(out) :: dxi_ad

    real(kind=kind_real), allocatable :: aice(:,:)
    real(kind=kind_real) :: A, B, C
    
    allocate( aice(xb%nx,xb%ny) )
    
    aice = sum(xb%cicen,3)

    A = (rho_s*c0)**(-1.0_kind_real)
    B = 0.0_kind_real !-mu*rho_i(L0/xb%
    C = 0.0_kind_real

    call zeros(dxi_ad)
    dxi_ad%cicen = dx_ad%cicen + xb%hicen * dx_ad%vicen! + xb%hsnon * dx_ad%vsnon - mu*xb%sssoc * dx_ad%sstoc
    dxi_ad%hicen = dx_ad%hicen + xb%cicen * dx_ad%vicen
    dxi_ad%vicen = 0.0_kind_real
    dxi_ad%hsnon = dx_ad%hsnon + xb%cicen * dx_ad%vsnon
    dxi_ad%vsnon = 0.0_kind_real
    dxi_ad%tsfcn = dx_ad%tsfcn      !<---- Wrong below that line ... CHECK
    dxi_ad%qsnon = 0.0_kind_real 
    dxi_ad%sicnk = 0.0_kind_real 
    dxi_ad%qicnk = 0.0_kind_real 
    dxi_ad%sssoc = dx_ad%sssoc - mu*aice * dx_ad%sstoc
    dxi_ad%tlioc = dx_ad%tlioc + (1.0_kind_real-aice) * dx_ad%sstoc
    dxi_ad%sstoc = 0.0_kind_real

    deallocate(aice)

    return
  end subroutine Kop_inv_ad


end module mom5cice5_thermo
