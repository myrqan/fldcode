program main
  use mpi
  implicit none
  integer :: nprocs,my_rank,ierr

  integer,parameter :: nn=10000
  real(8),allocatable :: vec1(:), vec2(:)
  integer :: n_local
  real(8) :: global_ipr,local_ipr
  integer :: chunk,istart,iend

  integer :: ii
  !========================================================!

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

  !========================================================!

  chunk = nn / nprocs
  istart = my_rank * chunk + 1
  if (my_rank == nprocs-1) then
    iend = nn
  else
    iend = istart + chunk - 1
  endif

  n_local = iend - istart + 1
  allocate(vec1(n_local), vec2(n_local))

  do ii = 1,n_local
    vec1(ii) = dsin((istart+ii-1)*0.1d0)
    vec2(ii) = dcos((istart+ii-1)*0.1d0)
  enddo

  local_ipr = 0.d0
  !global_ipr = 0.d0

  do ii = 1,n_local
    local_ipr = local_ipr + vec1(ii)*vec2(ii)
  enddo

  call MPI_ALLREDUCE(local_ipr,global_ipr,1,&
    & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  if(my_rank == 0) then
    write(*,*) "inner prod. = ", global_ipr
  endif


  !========================================================!

  call MPI_FINALIZE(ierr)
  stop

end program main
