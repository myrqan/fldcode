module vars
  use mpi
  implicit none

  integer,save :: nprocs,my_rank,ierr
  integer::mplx,mplz
  integer,parameter:: mpx=1
  integer,parameter:: mpz=1
  integer:: mpall=mpx*mpz
end module vars
