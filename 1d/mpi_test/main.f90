program main
  use mpi
  implicit none
  integer :: nprocs,my_rank,ierr
  integer,PARAMETER :: gx = 16
  integer,PARAMETER :: mpix = 4
  integer,PARAMETER :: mg = 1
  integer,PARAMETER :: lix = 2*mg+gx/mpix
  integer,PARAMETER :: gix = 2*mg+gx
  integer :: lx,ii
  integer :: ifro,ifx
  integer :: write_cnt

  integer(kind=MPI_OFFSET_KIND) :: disp,offset,file_end_pos
  integer :: status(MPI_STATUS_SIZE)

  real(8) :: xx(lix),xm(lix),ro(lix)
  real(8) :: xmin,xmax,dx,g_idx_ofs



  !========================================================!

  call MPI_Init(ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD,nprocs,ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)

  !========================================================!

  !! clean data
  if(my_rank == 0) then
    call EXECUTE_COMMAND_LINE("rm dat/*")
  end if


  !========================================================!

  if(nprocs /= mpix) then
    write(*,*) "number of processes is not consistent"
    write(*,*) "change mpix in main.f90 or MPINUM in makefile"
    write(*,*) "nprocs= ", nprocs
    write(*,*) "mpix= ", mpix
    stop
  end if

  !========================================================!

  !gx = gix-2*margin
  lx = lix-2*mg
  xmin = 0.d0
  xmax = 1.d0
  dx = (xmax-xmin)/dble(gx-1)

  do ii = 1, lix
    g_idx_ofs = dble(my_rank*lx+ii-mg-1)
    xx(ii) = xmin+g_idx_ofs*dx
  end do

  do ii = 1, lix-1
    xm(ii) = 0.5d0 * (xx(ii) + xx(ii+1))
  end do
  xm(lix) = xm(lix-1)+dx


  do ii = 1,lix
    if(xx(ii) <= 0.5d0) then
      ro(ii) = 1.d0
    else
      ro(ii) = 0.125d0
    end if
  end do

  !========================================================!
  !output~
  ! do ii = 0,nprocs-1
  !   call MPI_Barrier(MPI_Comm_World,ierr)
  !   if(ii == my_rank) then
  !     write(*,*) "rank", my_rank
  !     write(*,*) ro
  !     call flush(6)
  !     !call EXECUTE_COMMAND_LINE("sleep 1")
  !   end if
  ! end do
  !========================================================!
  ! file output test~
  write_cnt = lix-2*mg
  disp = int(write_cnt,kind=MPI_OFFSET_KIND)* 8_MPI_OFFSET_KIND
  offset=int(my_rank,kind=MPI_OFFSET_KIND)*disp

  call MPI_File_Open(MPI_Comm_World,'dat/ro.dat',&
    MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifro,ierr)
  call MPI_File_Open(MPI_Comm_World,'dat/x.dat',&
    MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifx,ierr)

  call MPI_FILE_WRITE_AT_ALL(ifro,offset,ro(1+mg),write_cnt,&
    MPI_DOUBLE_PRECISION, status, ierr)
  call MPI_FILE_WRITE_AT_ALL(ifx,offset,xx(1+mg),write_cnt,&
    MPI_DOUBLE_PRECISION, status, ierr)

  call MPI_File_Close(ifro,ierr)
  call MPI_File_Close(ifx,ierr)
  !========================================================!


  do ii = 1,lix
    ro(ii) = ro(ii) * 1.5d0
  end do

  call MPI_File_Open(MPI_Comm_World,'dat/ro.dat',&
    MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifro,ierr)
  call MPI_File_Get_Size(ifro,file_end_pos,ierr)
  offset = file_end_pos + (int(my_rank,kind=MPI_OFFSET_KIND)*disp)
  call MPI_FILE_WRITE_AT_ALL(ifro,offset,ro(1+mg),write_cnt,&
    MPI_DOUBLE_PRECISION, status, ierr)

  call MPI_File_Close(ifro,ierr)

  call MPI_FINALIZE(ierr)
  stop

end program main
