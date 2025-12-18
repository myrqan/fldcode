module datw
  use vars
  use mpi
  implicit none
contains
  subroutine dataclean()
    ! remove dat/*.dat
    use,intrinsic:: iso_fortran_env,only: error_unit
    call EXECUTE_COMMAND_LINE("rm dat/*.dat")
    call EXECUTE_COMMAND_LINE("rm param.txt")
  end subroutine dataclean

  subroutine put_param_int(char,num)
    character(*),intent(in):: char
    integer,intent(in):: num
    integer:: unum
    open(newunit=unum,file='param.txt',form='formatted',&
      position='append')
    write(unum,*) char, num
    close(unum)
  end subroutine put_param_int

  subroutine put_param_real(char,num)
    character(*),intent(in):: char
    real(8),intent(in):: num
    integer:: unum
    open(newunit=unum,file='param.txt',form='formatted',&
      position='append')
    write(unum,*) char, num
    close(unum)
  end subroutine put_param_real

  subroutine put_grid_data(lix,liz,mg,xx,zz)
    integer,intent(in)::lix,liz,mg
    real(8),intent(in)::xx(lix),zz(liz)
    integer:: write_cntx,write_cntz,write_cntxz
    integer(kind=MPI_OFFSET_KIND):: dispx,dispz,dispxz,&
      offsetx,offsetz
    integer:: ifxx,ifzz
    integer:: status(MPI_STATUS_SIZE)
   
    ! x-dir grid
    if(mplz==0) then
      write_cntx=(lix-2*mg) ! output data size (number)
    else
      write_cntx=0
    end if
    ! output type : real(8)
    dispx = int(lix-2*mg,kind=MPI_OFFSET_KIND)*8_MPI_OFFSET_KIND
    offsetx=int(mplx,kind=MPI_OFFSET_KIND)*dispx
    call MPI_File_Open(MPI_Comm_World,'dat/x.dat',&
      MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifxx,ierr)
    call MPI_File_Write_at_All(ifxx,offsetx,xx(1+mg),write_cntx,&
      MPI_Double_Precision,status,ierr)
    call MPI_File_Close(ifxx,ierr)

    ! z-dir grid
    if(mplx==0) then
      write_cntz=(liz-2*mg) ! output data size (number)
    else
      write_cntz=0
    end if
    ! output type: real(8)
    dispz = int(liz-2*mg,kind=MPI_OFFSET_KIND)*8_MPI_OFFSET_KIND
    offsetz=int(mplz,kind=MPI_OFFSET_KIND)*dispz
    call MPI_File_Open(MPI_Comm_World,'dat/z.dat',&
      MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifzz,ierr)
    call MPI_File_Write_at_All(ifzz,offsetz,zz(1+mg),write_cntz,&
      MPI_Double_Precision,status,ierr)
    call MPI_File_Close(ifzz,ierr)
  end subroutine put_grid_data

  subroutine put_time_data(fname,time)
    character(*),intent(in):: fname
    real(8),intent(in):: time
    integer:: fnum
    open(newunit=fnum,file=fname,position='append',&
      form='unformatted',access='stream')
    write(fnum) time
    close(fnum)
  end subroutine put_time_data

  subroutine put_2d_data_each_rank(fname,lix,liz,mg,qq)
    character(*),intent(in):: fname
    integer,intent(in):: lix,liz,mg
    real(8),intent(in):: qq(lix,liz)
    integer:: memtype,ifh
    integer:: lx,lz
    integer:: msize(2),lsize(2),start(2)
    character(80):: rank_fname
    
    lx=lix-2*mg
    lz=liz-2*mg
    msize=(/lix,liz/)
    lsize=(/lx,lz/)
    start=(/mg,mg/)

    write(rank_fname,'(A,"_",i0,"_",i0,".dat")')&
      trim(fname),mplx,mplz

    call MPI_Type_Create_Subarray(2,msize,lsize,start,&
      MPI_Order_fortran,MPI_Double_Precision,memtype,ierr)
    call MPI_Type_Commit(memtype,ierr)
    call MPI_File_Open(MPI_Comm_SELF,trim(rank_fname),&
      MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifh,ierr)
    call MPI_File_Write(ifh,qq,1,memtype,MPI_Status_Ignore,ierr)
    call MPI_File_Close(ifh,ierr)
    call MPI_Type_Free(memtype,ierr)
  end subroutine put_2d_data_each_rank
end module datw
