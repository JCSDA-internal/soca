
!> Invent an initial state for the QG model.

!> This routine invent an initial state for the QG model. It is used to
!! initialise the "truth run". The initial state consists of a horizontally
!! uniform wind in each layer, with a vertical shear sufficient to produce
!! baroclinic instability. Povided the orography is non-zero and is not
!! symmetrically place in the domain, this is sufficient to generate a
!! non-trivial flow after a few days of integration.
!!
!! Two slightly different initial states may be created (according to whether
!! or not ctype is set to 'f').

subroutine invent_state(flds,config)

use mom5cice5_fields
use iso_c_binding
use config_mod
use fckit_log_module, only : log
use mom5cice5_constants, only: u1,u2,bet,worog, domain_zonal, domain_meridional, &
                      & dlogtheta,f0,g,horog,scale_length,rossby_number
use kinds

implicit none

type(mom5cice5_field), intent(inout) :: flds    !< Model fields
type(c_ptr), intent(in)       :: config  !< Configuration structure

integer :: jx,jy,icentre,jcentre,ii,jj,ipert
real(kind=kind_real) :: distx,disty,d1,d2,f1,f2,rsmax,deltax0,deltay0,deltax,deltay,zz
real(kind=kind_real), allocatable :: pv(:,:,:), rs(:,:)
character(len=160) :: record

! ------------------------------------------------------------------------------

d1  = config_get_real(config,"top_layer_depth")
d2  = config_get_real(config,"bottom_layer_depth")

f1 = f0*f0*scale_length*scale_length/(g*dlogtheta*d1)
f2 = f0*f0*scale_length*scale_length/(g*dlogtheta*d2)
rsmax = horog/(rossby_number*d2)
deltax0 = domain_zonal/real(flds%nx,kind_real)
deltay0 = domain_meridional/real(flds%ny+1,kind_real)
deltax = deltax0/scale_length
deltay = deltay0/scale_length

allocate(pv(flds%nx,flds%ny,2))
allocate(rs(flds%nx,flds%ny))

!--- Uniform wind in each layer. A vertical shear outside the range
!--- -bet/F1 to bet/F2 should produce baroclinic instability.

flds%x_south(1) = 0.0_kind_real
flds%x_south(2) = 0.0_kind_real
flds%x_north(1) = -real(flds%ny+1,kind_real)*deltay*u1
flds%x_north(2) = -real(flds%ny+1,kind_real)*deltay*u2

do jy=1,flds%ny
do jx=1,flds%nx
  flds%x(jx,jy,1) = -real(jy,kind_real)*deltay*u1
  flds%x(jx,jy,2) = -real(jy,kind_real)*deltay*u2
enddo
enddo

ipert = config_get_int(config,"perturb")
if (ipert/=0) then
  write(record,*)"mom5cice5_invent_state_f90: Perturbing invented state by ",ipert,"%."
  call log%info(record)
  zz=real(ipert,kind_real)/100.0_kind_real
  do jy=1,flds%ny
  do jx=1,flds%nx
    flds%x(jx,jy,1) = (1.0_kind_real+zz)*flds%x(jx,jy,1)
    flds%x(jx,jy,2) = (1.0_kind_real-zz)*flds%x(jx,jy,2)
  enddo
  enddo
endif

icentre=flds%nx/4
jcentre=3*flds%ny/4
do jj=1,flds%ny
  do ii=1,flds%nx
    distx = real(min(icentre-ii,flds%nx-(icentre-ii)),kind_real) * deltax0
    disty = real(abs(jj-jcentre),kind_real) * deltay0
    rs(ii,jj) = rsmax*exp(-(distx*distx+disty*disty)/(worog*worog))
  enddo
enddo

call calc_pv(flds%nx,flds%ny,pv,flds%x,flds%x_north,flds%x_south, &
           & f1,f2,deltax,deltay,bet,rs)

do jx=1,flds%nx
  flds%q_south(jx,1) = 2.0_kind_real*pv(jx,1,1)-pv(jx,2,1)
  flds%q_south(jx,2) = 2.0_kind_real*pv(jx,1,2)-pv(jx,2,2)
enddo
do jx=1,flds%nx
  flds%q_north(jx,1) = 2.0_kind_real*pv(jx,flds%ny,1)-pv(jx,flds%ny-1,1)
  flds%q_north(jx,2) = 2.0_kind_real*pv(jx,flds%ny,2)-pv(jx,flds%ny-1,2)
enddo

deallocate(pv,rs)
! ------------------------------------------------------------------------------

return
end subroutine invent_state
