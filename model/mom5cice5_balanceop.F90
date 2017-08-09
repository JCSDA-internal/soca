module mom5cice5_balanceop

  use mom5cice5_constants
  use mom5cice5_fields

  implicit none
  
contains

  ! ------------------------------------------------------------------------------

  subroutine Kop_inv(dx, xb, dxi)

    ! dxi = Kop_inv(xb) * dx 
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
!!$
!!$  subroutine Kop_inv_ad(dx_ad, xb, dxi_ad)
!!$    ! dxi_ad = Kop_inv(xb)^T * dx_ad 
!!$    implicit none    
!!$    type(mom5cice5_field), intent(in) :: dx_ad
!!$    type(mom5cice5_field), intent(in) :: xb
!!$    type(mom5cice5_field), intent(out) :: dxi_ad
!!$
!!$    real(kind=kind_real), allocatable :: aice(:,:)
!!$    real(kind=kind_real) :: A, B, C
!!$    
!!$    allocate( aice(xb%nx,xb%ny) )
!!$    
!!$    aice = sum(xb%cicen,3)
!!$
!!$    A = (rho_s*c0)**(-1.0_kind_real)
!!$    B = 0.0_kind_real !-mu*rho_i(L0/xb%
!!$    C = 0.0_kind_real
!!$
!!$    call zeros(dxi_ad)
!!$    dxi_ad%cicen = dx_ad%cicen + xb%hicen * dx_ad%vicen! + xb%hsnon * dx_ad%vsnon - mu*xb%sssoc * dx_ad%sstoc
!!$    dxi_ad%hicen = dx_ad%hicen + xb%cicen * dx_ad%vicen
!!$    dxi_ad%vicen = 0.0_kind_real
!!$    dxi_ad%hsnon = dx_ad%hsnon + xb%cicen * dx_ad%vsnon
!!$    dxi_ad%vsnon = 0.0_kind_real
!!$    dxi_ad%tsfcn = dx_ad%tsfcn      !<---- Wrong below that line ... CHECK
!!$    dxi_ad%qsnon = 0.0_kind_real 
!!$    dxi_ad%sicnk = 0.0_kind_real 
!!$    dxi_ad%qicnk = 0.0_kind_real 
!!$    dxi_ad%sssoc = dx_ad%sssoc - mu*aice * dx_ad%sstoc
!!$    dxi_ad%tlioc = dx_ad%tlioc + (1.0_kind_real-aice) * dx_ad%sstoc
!!$    dxi_ad%sstoc = 0.0_kind_real
!!$
!!$    deallocate(aice)
!!$
!!$    return
!!$  end subroutine Kop_inv_ad


end module mom5cice5_balanceop
