module vars
  use mpi
  implicit none

  integer,save :: nprocs,my_rank,ierr
  integer::mplx,mplz
end module vars
