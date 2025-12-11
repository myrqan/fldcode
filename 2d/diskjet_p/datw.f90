module datw
  use vars,only: nprocs,my_rank,ierr
  use mpi
  implicit none
contains
  subroutine dataclean
    ! remove dat/*.dat
    use,intrinsic:: iso_fortran_env,only:error_unit
    call EXECUTE_COMMAND_LINE("rm dat/*.dat")
    call EXECUTE_COMMAND_LINE("rm param.txt")
  end subroutine dataclean
end module datw
