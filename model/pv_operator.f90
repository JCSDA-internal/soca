
!> Potential vorticity operator

!> Applies the linear operator part of the PV calculation:
!! \f{eqnarray*}{
!! q_1 &=& \nabla^2 \psi_1 - F_1 (\psi_1 -\psi_2 ) \\\\
!! q_2 &=& \nabla^2 \psi_2 - F_2 (\psi_2 -\psi_1 )
!! \f}
!!
!! (Note: The full potential vorticity calculation is done in calc_pv, and
!! includes additional beta and orography terms.)

subroutine pv_operator (x,pv,nx,ny,F1,F2,deltax,deltay)

!--- The part of the pv calculation that acts on internal streamfunction.

use kinds

implicit none
integer, intent(in) :: nx        !< Zonal grid dimension
integer, intent(in) :: ny        !< Meridional grid dimension
real(kind=kind_real), intent(in)  :: x(nx,ny,2)  !< Streamfunction
real(kind=kind_real), intent(out) :: pv(nx,ny,2) !< Result of applying the operator to x
real(kind=kind_real), intent(in) :: F1           !< Parameter in the PV operator
real(kind=kind_real), intent(in) :: F2           !< Parameter in the PV operator
real(kind=kind_real), intent(in) :: deltax       !< Zonal grid spacing (non-dimensional)
real(kind=kind_real), intent(in) :: deltay       !< Meridional grid spacing (non-dimensional)

!--- del-squared of the streamfunction

call laplacian_2d (x(:,:,1),pv(:,:,1),nx,ny,deltax,deltay)
call laplacian_2d (x(:,:,2),pv(:,:,2),nx,ny,deltax,deltay)

!--- vertical differences:

pv(:,:,1) = pv(:,:,1) -F1*(x(:,:,1)-x(:,:,2))
pv(:,:,2) = pv(:,:,2) -F2*(x(:,:,2)-x(:,:,1))

end subroutine pv_operator
