
subroutine generate_locations(c_conf,nlocs,ntimes,bgn,step,times,obsloc)

use iso_c_binding
use config_mod
use datetime_mod
use duration_mod
use mom5cice5_obs_vectors
use kinds

implicit none
type(c_ptr), intent(in) :: c_conf
integer, intent(in)     :: nlocs
integer, intent(in)     :: ntimes
type(datetime), intent(in) :: bgn
type(duration), intent(in) :: step
type(datetime), intent(inout) :: times(nlocs*ntimes)
type(obs_vect), intent(inout) :: obsloc

integer :: ijk, jobs, intinv, iobs, jstep, ii, jj, kk, ij, nx, ny
real(kind=kind_real) :: goldinv, dx, dy
real(kind=kind_real), allocatable :: xxyyzz(:,:)
type(datetime) :: now

nx = config_get_int(c_conf, "nx");
ny = config_get_int(c_conf, "ny");
dx=1.0_kind_real/real(nx,kind_real)
dy=1.0_kind_real/real(ny,kind_real)

goldinv = 0.5_kind_real*(sqrt(5.0_kind_real)-1.0_kind_real) ! 1/(golden ratio)
intinv = nint(goldinv*real(nx*ny*2,kind_real))

allocate(xxyyzz(3,nlocs))
ijk=0
do jobs=1,nlocs
  ijk=ijk+intinv
  ijk=modulo(ijk,2*nx*(ny-2))
  kk=ijk/(nx*(ny-2))
  ij=ijk-kk*nx*(ny-2)
  jj=ij/nx
  ii=ij-jj*nx
  jj=jj+1
  if (ii<0 .or. ii>nx-1) call abor1_ftn ("generate_locations: error ii")
  if (jj<1 .or. jj>ny-2) call abor1_ftn ("generate_locations: error jj")
  if (kk<0 .or. kk>1)  call abor1_ftn ("generate_locations: error kk")
  xxyyzz(1,jobs)=real(ii,kind_real)*dx
  xxyyzz(2,jobs)=real(jj,kind_real)*dy
  xxyyzz(3,jobs)=real(kk+1,kind_real)
enddo

call obsvec_setup(obsloc,3,nlocs*ntimes)

now = bgn
iobs=0
do jstep=1,ntimes
  do jobs=1,nlocs
    iobs=iobs+1
    times(iobs)=now
    obsloc%values(:,iobs)=xxyyzz(:,jobs)
  enddo
  call datetime_update(now,step)
enddo

call datetime_delete(now)
deallocate(xxyyzz)

end subroutine generate_locations
