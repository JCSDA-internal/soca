
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

  integer :: ijk, jobs, intinv, iobs, jstep, ii, jj, kk, ij
  real(kind=kind_real), allocatable :: xxyyzz(:,:)
  type(datetime) :: now

  !!!!!!!!!!!! NOT DONE YET !!!!!!!!!!!!
  allocate(xxyyzz(3,nlocs))
  print *,'-----GENERATE ',nlocs,' OBS -----'
  ijk=0
  call random_number(xxyyzz)

  do jobs=1,nlocs
     xxyyzz(1,jobs)=360.0_kind_real*xxyyzz(1,jobs)
     xxyyzz(2,jobs)=180.0_kind_real*xxyyzz(2,jobs)-90.0_kind_real
     xxyyzz(3,jobs)=0.0_kind_real
     print *,'generate location index:',jobs,xxyyzz(:,jobs)
     
  enddo

  call obsvec_setup(obsloc,3,nlocs*ntimes)

  now = bgn
  iobs=0
  do jstep=1,ntimes
     do jobs=1,nlocs
        iobs=iobs+1
        times(iobs)=now
        obsloc%values(1:3,iobs)=xxyyzz(1:3,jobs)
     enddo
     call datetime_update(now,step)
  enddo

  call datetime_delete(now)
  deallocate(xxyyzz)

end subroutine generate_locations
