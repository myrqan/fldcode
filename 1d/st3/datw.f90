module datw
  implicit none
contains
  subroutine dataclean()
    use,intrinsic :: iso_fortran_env,only:error_unit
    ! clean data
    call execute_command_line("rm dat/*")
  end subroutine dataclean

  subroutine w0d(qq,fname)
    ! write array(ar) in unformatted style
    real(8),intent(in) :: qq
    character(len=*),intent(in) :: fname
    integer :: fnum

    open(newunit=fnum,file=fname,position='append',&
      & form='unformatted')
    write(fnum) qq
    close(fnum)
  end subroutine w0d

  subroutine w1d(ix,ar,fname)
    ! write array(ar) in unformatted style
    integer,intent(in) :: ix
    real(8),intent(in) :: ar(ix)
    character(len=*),intent(in) :: fname
    integer :: fnum

    open(newunit=fnum,file=fname,position='append',&
      & form='unformatted')
    write(fnum) ar
    close(fnum)
  end subroutine w1d
end module datw
