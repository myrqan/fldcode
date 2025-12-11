module datw
  use vars, only : nprocs,my_rank,ierr
  use mpi
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

  subroutine w1d_b_mpi(fname,ifqq,lix,mg)
    character(len=*),intent(in) :: fname
    integer,intent(in) :: ifqq,lix,mg

    integer :: write_cnt
    integer(kind=MPI_OFFSET_KIND) :: disp,offset,file_end_pos
    write_cnt = lix-2*mg
    disp=int(write_cnt,kind=MPI_OFFSET_KIND)*8_MPI_OFFSET_KIND
    offset=int(my_rank,kind=MPI_OFFSET_KIND)*disp

  
  end subroutine w1d_b_mpi


end module datw
