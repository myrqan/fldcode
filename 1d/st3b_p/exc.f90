module exc
  use vars, only : nprocs,my_rank,ierr
  use mpi
  implicit none
contains
  subroutine exchange_1d(qq,lix,mg)
    integer,intent(in) :: lix,mg
    real(8),intent(inout) :: qq(lix)
    integer :: nrank,prank
    integer :: status(MPI_Status_Size)
    integer :: tag = 0

    nrank = my_rank+1
    prank = my_rank-1

    if(prank<0) then
      prank = MPI_Proc_Null
    end if
    if(nrank>=nprocs) then
      nrank = MPI_Proc_Null
    end if
    
    call MPI_Sendrecv(&
      qq(lix-2*mg+1),mg,MPI_Double_Precision,nrank,tag,&
      qq(1)         ,mg,MPI_Double_Precision,prank,tag,&
      MPI_Comm_World,status,ierr)
    call MPI_Sendrecv(&
      qq(mg+1)    ,mg,MPI_Double_Precision,prank,tag,&
      qq(lix-mg+1),mg,MPI_Double_Precision,nrank,tag,&
      MPI_Comm_World,status,ierr)
    !sendrecvは2回呼ぶ必要がある
  end subroutine exchange_1d

end module
