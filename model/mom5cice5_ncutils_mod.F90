module ncutils

  use netcdf

  implicit none

contains
  subroutine nccheck(status)

    integer, intent (in) :: status

    if (status /= nf90_noerr) then
       !print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if

  end subroutine nccheck

  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i3.3)') i
    res = trim(tmp)
  end function itoa

  subroutine fld_name_int2str(basename, level, varname)
    integer, intent(in) :: level
    character(len=128), intent(in) :: basename
    character(len=8), intent(inout) :: varname
    varname=trim(basename)//itoa(level)
  end subroutine fld_name_int2str

end module ncutils

module interface_ncread_fld
  
  use kinds
  use netcdf
  use ncutils

  implicit none

  interface ncread_fld
     module procedure read2d
     module procedure read3d
     module procedure read2d_from4d
  end interface ncread_fld
contains

  subroutine read2d(filename, varname, VAR, n1, n2)
    ! Read 2d field from 2d nc field
    implicit none

    character(len=128), intent(in)                                    :: filename
    character(len=128), intent(in)                                         :: varname
    real(kind=kind_real), allocatable, dimension(:,:), intent(inout)  :: VAR 
    integer                                                           :: varid, fid_in, n1, n2

    call nccheck(nf90_open(filename, nf90_nowrite, fid_in))
    call nccheck(nf90_inq_varid(fid_in, varname, varid))
    call nccheck(nf90_get_var(fid_in, varid, VAR))
    call nccheck(nf90_close(fid_in))

  end subroutine read2d

  subroutine read3d(filename, varname, VAR, n1, n2, n3)
    ! Read 3d field from 3d nc field
    implicit none

    character(len=128), intent(in)                                          :: filename
    character(len=128), intent(in)                                          :: varname
    real(kind=kind_real), allocatable, dimension(:,:,:), intent(inout)      :: VAR 
    integer                                                                 :: varid, fid_in, n1, n2, n3

    print *,'Reading ',varname
    call nccheck(nf90_open(filename, nf90_nowrite, fid_in) )
    call nccheck(nf90_inq_varid(fid_in, varname, varid))
    call nccheck(nf90_get_var(fid_in, varid, VAR))
    call nccheck(nf90_close(fid_in))

  end subroutine read3d

  subroutine read2d_from4d(filename, varname, VAR, n1, n2, start, count)
    ! Read 2d field from 4d nc field
    implicit none

    character(len=128), intent(in)                                     :: filename
    character(len=128), intent(in)                                     :: varname
    real(kind=kind_real), allocatable, dimension(:, :), intent(inout)  :: VAR 
    integer, intent(in)                                                :: n1, n2, start(4), count(4)
    integer                                                            :: varid, fid_in

    call nccheck(nf90_open(filename, nf90_nowrite, fid_in) )
    call nccheck(nf90_inq_varid(fid_in, varname, varid))
    call nccheck(nf90_get_var(fid_in, varid, VAR))
    call nccheck(nf90_close(fid_in))

  end subroutine read2d_from4d

end module interface_ncread_fld
