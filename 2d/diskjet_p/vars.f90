module vars
  use mpi
  implicit none

  integer,save :: nprocs,my_rank,ierr
  integer::mplx,mplz
  integer,parameter:: mpx=2
  integer,parameter:: mpz=2
  integer:: mpall=mpx*mpz
end module vars
