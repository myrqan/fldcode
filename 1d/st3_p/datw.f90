module datw
  use vars, only : nprocs,my_rank,ierr
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

    !open(newunit=fnum,file=fname,position='append',&
    !  & form='unformatted')

    open(newunit=fnum,file=fname,position='append',&
      & form='unformatted',access='stream')
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

  !subroutine w1d_mpi(ifile,qq,fname,lix,mg)
  !  integer,intent(in) :: ifile,lix,mg
  !  character(len=*),intent(in)  :: fname
  !  real(8),intent(in) :: qq(lix)
  !end subroutine w1d_mpi

end module datw
