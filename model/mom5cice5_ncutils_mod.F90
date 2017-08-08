module ncutils

  !use netcdf
  implicit none

contains
  subroutine nccheck(status)

    !use netcdf

    integer, intent (in) :: status

    !if (status /= nf90_noerr) then
       !print *, trim(nf90_strerror(status))
    !   stop "Stopped"
    !end if

  end subroutine nccheck

end module ncutils
