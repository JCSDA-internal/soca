!
! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

module soca_interph_mod

  use kinds
  
  implicit none
  private

  type, public :: soca_hinterp
     integer, allocatable              :: index(:,:,:)  !< Indices of nearest neighbors
                                                        !< (nobs x 2 x number of neighbors)
     real(kind=kind_real), allocatable :: wgh(:,:)      !< Interp weight
                                                        !< (nobs x number of neighbors)     
     integer                           :: nn            !< Number of neighbors
     real(kind=kind_real)              :: lx, ly        !< Decorrelation length scales
     character(len=3)                  :: wgt_type      !< 'avg': wgt=1 (NOT IMPLEMENTED ANYMORE) 
                                                        !< 'bar': bilinear-interp weight 
     integer                           :: nobs          !< Number of values to interpolate
     real(kind=kind_real), allocatable :: lono(:)       !< Longitude of destination
     real(kind=kind_real), allocatable :: lato(:)       !< Latitude of destination 
     logical                           :: initialized   !< Initialization switch
     logical                           :: alloc         !< ... Not sure if this is still used ...
   contains
     procedure :: interp_init
     procedure :: interp_compute_weight
     procedure :: interp_apply     
     procedure :: interpad_apply
     procedure :: interp_exit
  end type soca_hinterp

contains

  !--------------------------------------------
  subroutine interp_init(self, nobs, nn, lx, ly, wgt_type)

    implicit none

    integer, intent(in)              :: nobs     !< Number of obs
    integer, optional                :: nn       !< Number of neighbors
    real(kind=kind_real), optional   :: lx       !< Decorrelation scales
    real(kind=kind_real), optional   :: ly       !< for dist wghted interp
    character(len=3), optional       :: wgt_type !< 'avg' or 'bar'     
    class(soca_hinterp), intent(out) :: self
        
    self%nn=5;          if (present(nn)) self%nn = nn
    self%wgt_type='bar'; if (present(nn)) self%wgt_type = wgt_type    
    self%lx=5e-2;        if (present(lx)) self%lx = lx
    self%ly=1e-2;        if (present(ly)) self%ly = ly

    self%nobs = nobs
    allocate(self%lono(nobs),self%lato(nobs))
    allocate(self%index(nobs,2,self%nn))
    allocate(self%wgh(nobs,self%nn))
    self%initialized = .false.
    self%alloc = .true.
    
  end subroutine interp_init

  !--------------------------------------------  
  subroutine interp_compute_weight(self, lon, lat, lono, lato)

    use kinds
    use type_ctree, only: ctree_type!,ctree_create,delete_ctree,find_nearest_neighbors
    use tools_const, only: pi,req,deg2rad,rad2deg
    use iso_fortran_env
    use mpi
    
    implicit none
    
    class(soca_hinterp), intent(inout)               :: self    
    real(kind=kind_real), dimension(:,:), intent(in) :: lon, lat
    real(kind=kind_real), dimension(:), intent(in)   :: lono, lato

    integer :: nobs, ni, nj, k, l, ij(2), cnt
    integer :: n, nn
    logical, allocatable :: mask(:)
    type(ctree_type) :: cover_tree
    real(kind=kind_real), allocatable :: nn_dist(:,:), tmplon(:), tmplat(:)
    real(kind=kind_real), allocatable :: tmplono(:), tmplato(:)
    integer, allocatable :: nn_index(:,:)              ! nobsxnn
    real(kind=kind_real) :: offset, dist
    integer :: ii,jj
    integer :: my_rank, ierr
    ni = size(lon,1)
    nj = size(lon,2)

    !--- Save destination locations
    self%lono = lono
    self%lato = lato
    
    !--- Initialize kd-tree
    n=ni*nj
    allocate(mask(n))
    mask=.true.
    allocate(tmplon(n),tmplat(n))
    tmplon=deg2rad*reshape(lon,(/n/))
    tmplat=deg2rad*reshape(lat,(/n/))
    call cover_tree%create(n,tmplon,tmplat,mask)

    !--- Find nn nearest neighbors
    nn = self%nn
    allocate(nn_dist(self%nobs,nn),nn_index(self%nobs,nn))
    allocate(tmplono(self%nobs),tmplato(self%nobs))
    tmplono=deg2rad*lono
    tmplato=deg2rad*lato

    call MPI_Comm_rank ( MPI_COMM_WORLD, my_rank, ierr )
    cnt = 0
    do k = 1, self%nobs
       if (my_rank.eq.0) then
          write(*,FMT="(A1,A,t21,F6.2,A,A)",ADVANCE="NO") achar(13), &
               & " Percent Complete: ", (real(k)/real(self%nobs))*100.0, "% for ",self%wgt_type
       end if
       call cover_tree%find_nearest_neighbors(tmplono(k),tmplato(k),nn,nn_index(k,:),nn_dist(k,:))       
       !nn_dist(k,:)=exp(-(nn_dist(k,:)/self%lx)**2)
       dist=sum(nn_dist(k,:))
       do l = 1, nn
          self%index(k,1,l)=max(1,mod(nn_index(k,l),ni))
          self%index(k,2,l)=max(1,nn_index(k,l)/ni+1)
          self%wgh(k,l)=(dist-nn_dist(k,l))
       end do
       self%wgh(k,:)=self%wgh(k,:)/sum(self%wgh(k,:))       
    end do

    call mpi_barrier(MPI_COMM_WORLD,ierr)
    self%initialized = .false.

    deallocate(mask,tmplon,tmplat,tmplono,tmplato)

  end subroutine interp_compute_weight

  !--------------------------------------------  
  subroutine interp_apply(self, fld, obs)
    ! Forward interpolation: "fields to obs"
    ! obs = interp(fld)
    use kinds

    implicit none

    class(soca_hinterp), intent(in) :: self    
    real(kind=kind_real), dimension(:,:), intent(in) :: fld
    real(kind=kind_real), dimension(:), intent(out) :: obs    
    integer :: k,l,i
    obs = 0.0
    do k = 1, self%nobs
       do l = 1, self%nn
          obs(k) = obs(k) + self%wgh(k,l)*fld(self%index(k,1,l),self%index(k,2,l))          
       end do
    end do
  end subroutine interp_apply

  !--------------------------------------------  
  subroutine interpad_apply(self, fld, obs)
    ! Backward interpolation: "obs to fields"
    ! fld = interpad(obs)    
    use kinds

    implicit none

    class(soca_hinterp), intent(in) :: self    
    real(kind=kind_real), dimension(:,:), intent(inout) :: fld
    real(kind=kind_real), dimension(:), intent(in) :: obs    
    integer :: k,l

    do k = 1, self%nobs
       do l = 1, self%nn
          fld(self%index(k,1,l),self%index(k,2,l))=fld(self%index(k,1,l),self%index(k,2,l))+&
               &self%wgh(k,l)*obs(k)
       end do
       !obs(k) = 0.0
    end do
    
  end subroutine interpad_apply
  
  !--------------------------------------------  
  subroutine interp_exit(self)

    implicit none
    
    class(soca_hinterp), intent(out) :: self

    deallocate(self%index,self%lono,self%lato,self%wgh)
    self%initialized = .false.
    self%alloc = .true.
    self%nobs = 0
    !More cleaning/deallocating needs to be done
    
  end subroutine interp_exit

end module soca_interph_mod

