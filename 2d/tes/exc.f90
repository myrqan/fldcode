module exc
  use mpi
  use vars
  implicit none
contains
  subroutine exc_2d_13var( &
      q1,q2,q3,q4,q5,q6,q7, &
      q8,q9,q10,q11,q12,q13, &
      lix,liz,mg)

    !---------------------------------------
    ! halo exchange based on exc_9
    !---------------------------------------
    integer,intent(in) :: lix,liz,mg
    real(8),intent(inout) :: &
      q1(lix,liz), q2(lix,liz), q3(lix,liz), &
      q4(lix,liz), q5(lix,liz), q6(lix,liz), &
      q7(lix,liz), q8(lix,liz), q9(lix,liz), &
      q10(lix,liz),q11(lix,liz),q12(lix,liz),q13(lix,liz)

    integer :: rank_xm, rank_xp
    integer :: rank_zm, rank_zp
    integer :: i,j,idx
    integer :: cntx, cntz
    integer :: status(MPI_STATUS_SIZE)

    real(8),allocatable :: sndx(:),rcvx(:)
    real(8),allocatable :: sndz(:),rcvz(:)

    !---------------------------------------
    ! neighbor ranks
    !---------------------------------------
    rank_xm = my_rank - 1
    rank_xp = my_rank + 1
    rank_zm = my_rank - mpx
    rank_zp = my_rank + mpx

    if (mplx == 0)       rank_xm = MPI_PROC_NULL
    if (mplx == mpx-1)   rank_xp = MPI_PROC_NULL
    if (mplz == 0)       rank_zm = MPI_PROC_NULL
    if (mplz == mpz-1)   rank_zp = MPI_PROC_NULL

    !---------------------------------------
    ! buffer sizes
    !---------------------------------------
    cntx = mg * liz * 13
    cntz = mg * lix * 13

    allocate(sndx(cntx), rcvx(cntx))
    allocate(sndz(cntz), rcvz(cntz))

    !==================================================================
    ! X-direction exchange
    !==================================================================

    ! send -> x-, recv <- x+
    idx = 0
    do j=1,liz
    do i=mg+1,2*mg
      idx=idx+1; sndx(idx)=q1(i,j)
      idx=idx+1; sndx(idx)=q2(i,j)
      idx=idx+1; sndx(idx)=q3(i,j)
      idx=idx+1; sndx(idx)=q4(i,j)
      idx=idx+1; sndx(idx)=q5(i,j)
      idx=idx+1; sndx(idx)=q6(i,j)
      idx=idx+1; sndx(idx)=q7(i,j)
      idx=idx+1; sndx(idx)=q8(i,j)
      idx=idx+1; sndx(idx)=q9(i,j)
      idx=idx+1; sndx(idx)=q10(i,j)
      idx=idx+1; sndx(idx)=q11(i,j)
      idx=idx+1; sndx(idx)=q12(i,j)
      idx=idx+1; sndx(idx)=q13(i,j)
    end do
    end do

    call MPI_Sendrecv( &
      sndx,cntx,MPI_DOUBLE_PRECISION,rank_xm,0, &
      rcvx,cntx,MPI_DOUBLE_PRECISION,rank_xp,0, &
      MPI_COMM_WORLD,status,ierr )

    if (rank_xp /= MPI_PROC_NULL) then
      idx=0
      do j=1,liz
      do i=lix-mg+1,lix
        idx=idx+1; q1(i,j)=rcvx(idx)
        idx=idx+1; q2(i,j)=rcvx(idx)
        idx=idx+1; q3(i,j)=rcvx(idx)
        idx=idx+1; q4(i,j)=rcvx(idx)
        idx=idx+1; q5(i,j)=rcvx(idx)
        idx=idx+1; q6(i,j)=rcvx(idx)
        idx=idx+1; q7(i,j)=rcvx(idx)
        idx=idx+1; q8(i,j)=rcvx(idx)
        idx=idx+1; q9(i,j)=rcvx(idx)
        idx=idx+1; q10(i,j)=rcvx(idx)
        idx=idx+1; q11(i,j)=rcvx(idx)
        idx=idx+1; q12(i,j)=rcvx(idx)
        idx=idx+1; q13(i,j)=rcvx(idx)
      end do
      end do
    end if

    ! send -> x+, recv <- x-
    idx = 0
    do j=1,liz
    do i=lix-2*mg+1,lix-mg
      idx=idx+1; sndx(idx)=q1(i,j)
      idx=idx+1; sndx(idx)=q2(i,j)
      idx=idx+1; sndx(idx)=q3(i,j)
      idx=idx+1; sndx(idx)=q4(i,j)
      idx=idx+1; sndx(idx)=q5(i,j)
      idx=idx+1; sndx(idx)=q6(i,j)
      idx=idx+1; sndx(idx)=q7(i,j)
      idx=idx+1; sndx(idx)=q8(i,j)
      idx=idx+1; sndx(idx)=q9(i,j)
      idx=idx+1; sndx(idx)=q10(i,j)
      idx=idx+1; sndx(idx)=q11(i,j)
      idx=idx+1; sndx(idx)=q12(i,j)
      idx=idx+1; sndx(idx)=q13(i,j)
    end do
    end do

    call MPI_Sendrecv( &
      sndx,cntx,MPI_DOUBLE_PRECISION,rank_xp,1, &
      rcvx,cntx,MPI_DOUBLE_PRECISION,rank_xm,1, &
      MPI_COMM_WORLD,status,ierr )

    if (rank_xm /= MPI_PROC_NULL) then
      idx=0
      do j=1,liz
      do i=1,mg
        idx=idx+1; q1(i,j)=rcvx(idx)
        idx=idx+1; q2(i,j)=rcvx(idx)
        idx=idx+1; q3(i,j)=rcvx(idx)
        idx=idx+1; q4(i,j)=rcvx(idx)
        idx=idx+1; q5(i,j)=rcvx(idx)
        idx=idx+1; q6(i,j)=rcvx(idx)
        idx=idx+1; q7(i,j)=rcvx(idx)
        idx=idx+1; q8(i,j)=rcvx(idx)
        idx=idx+1; q9(i,j)=rcvx(idx)
        idx=idx+1; q10(i,j)=rcvx(idx)
        idx=idx+1; q11(i,j)=rcvx(idx)
        idx=idx+1; q12(i,j)=rcvx(idx)
        idx=idx+1; q13(i,j)=rcvx(idx)
      end do
      end do
    end if

    !==================================================================
    ! Z-direction exchange (x-ghost INCLUDED)
    !==================================================================

    ! send -> z-, recv <- z+
    idx=0
    do j=mg+1,2*mg
    do i=1,lix
      idx=idx+1; sndz(idx)=q1(i,j)
      idx=idx+1; sndz(idx)=q2(i,j)
      idx=idx+1; sndz(idx)=q3(i,j)
      idx=idx+1; sndz(idx)=q4(i,j)
      idx=idx+1; sndz(idx)=q5(i,j)
      idx=idx+1; sndz(idx)=q6(i,j)
      idx=idx+1; sndz(idx)=q7(i,j)
      idx=idx+1; sndz(idx)=q8(i,j)
      idx=idx+1; sndz(idx)=q9(i,j)
      idx=idx+1; sndz(idx)=q10(i,j)
      idx=idx+1; sndz(idx)=q11(i,j)
      idx=idx+1; sndz(idx)=q12(i,j)
      idx=idx+1; sndz(idx)=q13(i,j)
    end do
    end do

    call MPI_Sendrecv( &
      sndz,cntz,MPI_DOUBLE_PRECISION,rank_zm,10, &
      rcvz,cntz,MPI_DOUBLE_PRECISION,rank_zp,10, &
      MPI_COMM_WORLD,status,ierr )

    if (rank_zp /= MPI_PROC_NULL) then
      idx=0
      do j=liz-mg+1,liz
      do i=1,lix
        idx=idx+1; q1(i,j)=rcvz(idx)
        idx=idx+1; q2(i,j)=rcvz(idx)
        idx=idx+1; q3(i,j)=rcvz(idx)
        idx=idx+1; q4(i,j)=rcvz(idx)
        idx=idx+1; q5(i,j)=rcvz(idx)
        idx=idx+1; q6(i,j)=rcvz(idx)
        idx=idx+1; q7(i,j)=rcvz(idx)
        idx=idx+1; q8(i,j)=rcvz(idx)
        idx=idx+1; q9(i,j)=rcvz(idx)
        idx=idx+1; q10(i,j)=rcvz(idx)
        idx=idx+1; q11(i,j)=rcvz(idx)
        idx=idx+1; q12(i,j)=rcvz(idx)
        idx=idx+1; q13(i,j)=rcvz(idx)
      end do
      end do
    end if

    ! send -> z+, recv <- z-
    idx=0
    do j=liz-2*mg+1,liz-mg
    do i=1,lix
      idx=idx+1; sndz(idx)=q1(i,j)
      idx=idx+1; sndz(idx)=q2(i,j)
      idx=idx+1; sndz(idx)=q3(i,j)
      idx=idx+1; sndz(idx)=q4(i,j)
      idx=idx+1; sndz(idx)=q5(i,j)
      idx=idx+1; sndz(idx)=q6(i,j)
      idx=idx+1; sndz(idx)=q7(i,j)
      idx=idx+1; sndz(idx)=q8(i,j)
      idx=idx+1; sndz(idx)=q9(i,j)
      idx=idx+1; sndz(idx)=q10(i,j)
      idx=idx+1; sndz(idx)=q11(i,j)
      idx=idx+1; sndz(idx)=q12(i,j)
      idx=idx+1; sndz(idx)=q13(i,j)
    end do
    end do

    call MPI_Sendrecv( &
      sndz,cntz,MPI_DOUBLE_PRECISION,rank_zp,11, &
      rcvz,cntz,MPI_DOUBLE_PRECISION,rank_zm,11, &
      MPI_COMM_WORLD,status,ierr )

    if (rank_zm /= MPI_PROC_NULL) then
      idx=0
      do j=1,mg
      do i=1,lix
        idx=idx+1; q1(i,j)=rcvz(idx)
        idx=idx+1; q2(i,j)=rcvz(idx)
        idx=idx+1; q3(i,j)=rcvz(idx)
        idx=idx+1; q4(i,j)=rcvz(idx)
        idx=idx+1; q5(i,j)=rcvz(idx)
        idx=idx+1; q6(i,j)=rcvz(idx)
        idx=idx+1; q7(i,j)=rcvz(idx)
        idx=idx+1; q8(i,j)=rcvz(idx)
        idx=idx+1; q9(i,j)=rcvz(idx)
        idx=idx+1; q10(i,j)=rcvz(idx)
        idx=idx+1; q11(i,j)=rcvz(idx)
        idx=idx+1; q12(i,j)=rcvz(idx)
        idx=idx+1; q13(i,j)=rcvz(idx)
      end do
      end do
    end if

    deallocate(sndx,rcvx,sndz,rcvz)

  end subroutine exc_2d_13var
end module exc
