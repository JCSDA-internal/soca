!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_Butils
  implicit none

contains

    ! ------------------------------------------------------------------------------  
  
  subroutine soca_k_mult(dx, traj)
    use kinds
    use soca_fields
    use soca_balanceop
    use soca_seaice_balanceop
    
    implicit none
    type(soca_field), intent(in) :: dx       !< Increment            
    type(soca_field), intent(in) :: traj     !< trajectory    

    !Grid stuff
    integer :: isc, iec, jsc, jec, i, j, k
    real(kind=kind_real) :: tb, sb, dt, ds,deta, z, p, h
    real(kind=kind_real), allocatable :: dcn(:), cnb(:), dtv(:), dsv(:)

    ! Indices for compute domain (no halo)
    isc = traj%geom%ocean%G%isc
    iec = traj%geom%ocean%G%iec
    jsc = traj%geom%ocean%G%jsc
    jec = traj%geom%ocean%G%jec

    ! Steric height/density balance
    do i = isc, iec
       do j = jsc, jec
          do k = 1, traj%geom%ocean%nzo
             tb=traj%tocn(i,j,k)
             sb=traj%socn(i,j,k)
             dt=dx%tocn(i,j,k)
             ds=dx%socn(i,j,k)
             if (k.eq.1) then
                z=traj%hocn(i,j,k)
             else
                z=sum(traj%hocn(i,j,1:k-1))+0.5_kind_real*traj%hocn(i,j,k)
             end if
             h=traj%hocn(i,j,k)
             p=z
             call soca_steric_tl(deta, dt, ds, tb, sb, p, h,&
                  &traj%geom%ocean%lon(i,j), traj%geom%ocean%lat(i,j))
             dx%ssh(i,j)=dx%ssh(i,j)+deta    
          end do
       end do
    end do

    ! T-S balance
    allocate(dsv(traj%geom%ocean%nzo),dtv(traj%geom%ocean%nzo))
    dsv=0.0
    dtv=0.0    
    do i = isc, iec
       do j = jsc, jec
          dtv = dx%tocn(i,j,:)
          call soca_soft_tl (dsv,dtv,&
                            &traj%tocn(i,j,:),&
                            &traj%socn(i,j,:),&
                            &traj%hocn(i,j,:))
          dx%socn(i,j,:) = dx%socn(i,j,:) + dsv
       end do
    end do
    
!!$    ! T/C balance
!!$    allocate(dcn(traj%geom%ocean%ncat),cnb(traj%geom%ocean%ncat))
!!$    do i = isc, iec
!!$       do j = jsc, jec
!!$          !tb=traj%tocn(i,j,1)
!!$          !sb=traj%socn(i,j,1)
!!$          !cnb=traj%cicen(i,j,2:)
!!$          dcn=dx%cicen(i,j,2:)          
!!$          !call tofc_tl (dt, dcn, tb, sb, cnb)
!!$          dx%tocn(i,j,1)=-sum(dcn)!+dt
!!$       end do
!!$    end do
!!$
!!$    deallocate(dcn)    
  end subroutine soca_k_mult

  ! ------------------------------------------------------------------------------  
  !> Adjoint of balance operators
  
  subroutine soca_k_mult_ad(dy, traj)
    use kinds
    use soca_fields
    use soca_balanceop
    use soca_seaice_balanceop
    
    implicit none
    !type(soca_field), intent(inout) :: KTdy     !< K^T dy
    type(soca_field), intent(inout) :: dy       !< Input:Increment
                                                !< Input:K^T dy
    type(soca_field), intent(in)    :: traj     !< trajectory    

    ! Grid stuff
    integer :: isc, iec, jsc, jec, i, j, k

    ! Convenience variables
    ! Adjoint of steric height: Inputs
    real(kind=kind_real) :: tb   !< Background potential temperature [C]
    real(kind=kind_real) :: sb   !< Background practical salinity [psu]
    real(kind=kind_real) :: z    !< Mid-layer depth [m]
    real(kind=kind_real) :: p    !< Pressure at mid-layer depth [dbar]
    real(kind=kind_real) :: h    !< Layer thickness [m]
    real(kind=kind_real) :: deta !< Sea surface height increment [m]
    
    ! Adjoint of steric height: Outputs 
    real(kind=kind_real) :: dt   !< Potential temperature increment [C] 
    real(kind=kind_real) :: ds   !< Practical salinity increment [psu]

    real(kind=kind_real), allocatable :: dcn(:), cnb(:), dtv(:), dsv(:)

    ! Indices for compute domain (no halo)
    isc = traj%geom%ocean%G%isc
    iec = traj%geom%ocean%G%iec
    jsc = traj%geom%ocean%G%jsc
    jec = traj%geom%ocean%G%jec

!!$    ! T/C balance
!!$    allocate(dcn(traj%geom%ocean%ncat),cnb(traj%geom%ocean%ncat))
!!$    do i = isc, iec
!!$       do j = jsc, jec
!!$          !tb = traj%tocn(i,j,1)
!!$          !sb = traj%socn(i,j,1)
!!$          !cnb = traj%cicen(i,j,2:)
!!$          dt = dy%tocn(i,j,1)
!!$          !call tofc_ad (dt, dcn, tb, sb, cnb)
!!$          do k = 1, traj%geom%ocean%ncat
!!$             dy%cicen(i,j,k+1) = dy%cicen(i,j,k+1) - dt! + dcn(k)
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    deallocate(dcn, cnb)

    ! T-S balance
    allocate(dsv(traj%geom%ocean%nzo),dtv(traj%geom%ocean%nzo))
    dsv=0.0
    dtv=0.0
    do i = isc, iec
       do j = jsc, jec
          dsv = dy%socn(i,j,:)          
          call soca_soft_ad (dsv,dtv,&
                            &traj%tocn(i,j,:),&
                            &traj%socn(i,j,:),&
                            &traj%hocn(i,j,:))
          dy%tocn(i,j,:) = dy%tocn(i,j,:) + dtv
       end do
    end do
    deallocate(dtv,dsv)
    
    ! Steric height/density balance
    do i = isc, iec
       do j = jsc, jec
          do k = traj%geom%ocean%nzo, 1, -1
             tb=traj%tocn(i,j,k)
             sb=traj%socn(i,j,k)
             if (k.eq.1) then
                z=traj%hocn(i,j,k)
             else
                z=sum(traj%hocn(i,j,1:k-1))+0.5_kind_real*traj%hocn(i,j,k)
             end if
             h=traj%hocn(i,j,k)
             p=z
             deta=dy%ssh(i,j)
             call soca_steric_ad(deta, dt, ds, tb, sb, p, h,&
                  &traj%geom%ocean%lon(i,j), traj%geom%ocean%lat(i,j))
             dy%tocn(i,j,k)=dy%tocn(i,j,k)+dt
             dy%socn(i,j,k)=dy%socn(i,j,k)+ds
          end do
       end do
    end do

  end subroutine soca_k_mult_ad

  
end module soca_Butils
