module exc
  use mpi
  use vars  ! my_rank, mpx, mpz, mplx, mplz, ierr が定義されていること
  implicit none

contains

  subroutine exc_2d_13var( &
      qq1, qq2, qq3, qq4, qq5, qq6, qq7, &
      qq8, qq9, qq10, qq11, qq12, qq13, &
      lix, liz, mg)

    integer, intent(in) :: lix, liz, mg
    real(8), intent(inout) :: &
      qq1(lix,liz),  qq2(lix,liz),  qq3(lix,liz),  &
      qq4(lix,liz),  qq5(lix,liz),  qq6(lix,liz),  &
      qq7(lix,liz),  qq8(lix,liz),  qq9(lix,liz),  &
      qq10(lix,liz), qq11(lix,liz), qq12(lix,liz), qq13(lix,liz)

    integer :: rank_xm, rank_xp, rank_zm, rank_zp
    integer :: i, j, idx, cntx, cntz
    integer :: status(MPI_STATUS_SIZE)
    real(8), allocatable :: sndx(:), rcvx(:), sndz(:), rcvz(:)

    !---------------------------------------------------------
    ! 1. 隣接プロセスの特定 (MPI_PROC_NULL は通信を自動スキップする)
    !---------------------------------------------------------
    rank_xm = my_rank - 1;   if (mplx == 0)       rank_xm = MPI_PROC_NULL
    rank_xp = my_rank + 1;   if (mplx == mpx - 1) rank_xp = MPI_PROC_NULL
    rank_zm = my_rank - mpx; if (mplz == 0)       rank_zm = MPI_PROC_NULL
    rank_zp = my_rank + mpx; if (mplz == mpz - 1) rank_zp = MPI_PROC_NULL

    ! バッファサイズの計算 (13変数分)
    cntx = mg * liz * 13
    cntz = mg * lix * 13

    allocate(sndx(cntx), rcvx(cntx))
    allocate(sndz(cntz), rcvz(cntz))

    !==================================================================
    ! 2. X-direction exchange (垂直境界)
    !==================================================================

    ! [A] Send to Left (XM), Receive from Right (XP)
    idx = 0
    do j = 1, liz
      do i = mg + 1, 2 * mg
        idx = idx + 1; sndx(idx) = qq1(i,j);  idx = idx + 1; sndx(idx) = qq2(i,j)
        idx = idx + 1; sndx(idx) = qq3(i,j);  idx = idx + 1; sndx(idx) = qq4(i,j)
        idx = idx + 1; sndx(idx) = qq5(i,j);  idx = idx + 1; sndx(idx) = qq6(i,j)
        idx = idx + 1; sndx(idx) = qq7(i,j);  idx = idx + 1; sndx(idx) = qq8(i,j)
        idx = idx + 1; sndx(idx) = qq9(i,j);  idx = idx + 1; sndx(idx) = qq10(i,j)
        idx = idx + 1; sndx(idx) = qq11(i,j); idx = idx + 1; sndx(idx) = qq12(i,j)
        idx = idx + 1; sndx(idx) = qq13(i,j)
      end do
    end do

    call MPI_Sendrecv(sndx, cntx, MPI_DOUBLE_PRECISION, rank_xm, 0, &
                      rcvx, cntx, MPI_DOUBLE_PRECISION, rank_xp, 0, &
                      MPI_COMM_WORLD, status, ierr)

    if (rank_xp /= MPI_PROC_NULL) then
      idx = 0
      do j = 1, liz
        do i = lix - mg + 1, lix
          idx = idx + 1; qq1(i,j) = rcvx(idx);  idx = idx + 1; qq2(i,j) = rcvx(idx)
          idx = idx + 1; qq3(i,j) = rcvx(idx);  idx = idx + 1; qq4(i,j) = rcvx(idx)
          idx = idx + 1; qq5(i,j) = rcvx(idx);  idx = idx + 1; qq6(i,j) = rcvx(idx)
          idx = idx + 1; qq7(i,j) = rcvx(idx);  idx = idx + 1; qq8(i,j) = rcvx(idx)
          idx = idx + 1; qq9(i,j) = rcvx(idx);  idx = idx + 1; qq10(i,j) = rcvx(idx)
          idx = idx + 1; qq11(i,j) = rcvx(idx); idx = idx + 1; qq12(i,j) = rcvx(idx)
          idx = idx + 1; qq13(i,j) = rcvx(idx)
        end do
      end do
    end if

    ! [B] Send to Right (XP), Receive from Left (XM)
    idx = 0
    do j = 1, liz
      do i = lix - 2 * mg + 1, lix - mg
        idx = idx + 1; sndx(idx) = qq1(i,j);  idx = idx + 1; sndx(idx) = qq2(i,j)
        idx = idx + 1; sndx(idx) = qq3(i,j);  idx = idx + 1; sndx(idx) = qq4(i,j)
        idx = idx + 1; sndx(idx) = qq5(i,j);  idx = idx + 1; sndx(idx) = qq6(i,j)
        idx = idx + 1; sndx(idx) = qq7(i,j);  idx = idx + 1; sndx(idx) = qq8(i,j)
        idx = idx + 1; sndx(idx) = qq9(i,j);  idx = idx + 1; sndx(idx) = qq10(i,j)
        idx = idx + 1; sndx(idx) = qq11(i,j); idx = idx + 1; sndx(idx) = qq12(i,j)
        idx = idx + 1; sndx(idx) = qq13(i,j)
      end do
    end do

    call MPI_Sendrecv(sndx, cntx, MPI_DOUBLE_PRECISION, rank_xp, 1, &
                      rcvx, cntx, MPI_DOUBLE_PRECISION, rank_xm, 1, &
                      MPI_COMM_WORLD, status, ierr)

    if (rank_xm /= MPI_PROC_NULL) then
      idx = 0
      do j = 1, liz
        do i = 1, mg
          idx = idx + 1; qq1(i,j) = rcvx(idx);  idx = idx + 1; qq2(i,j) = rcvx(idx)
          idx = idx + 1; qq3(i,j) = rcvx(idx);  idx = idx + 1; qq4(i,j) = rcvx(idx)
          idx = idx + 1; qq5(i,j) = rcvx(idx);  idx = idx + 1; qq6(i,j) = rcvx(idx)
          idx = idx + 1; qq7(i,j) = rcvx(idx);  idx = idx + 1; qq8(i,j) = rcvx(idx)
          idx = idx + 1; qq9(i,j) = rcvx(idx);  idx = idx + 1; qq10(i,j) = rcvx(idx)
          idx = idx + 1; qq11(i,j) = rcvx(idx); idx = idx + 1; qq12(i,j) = rcvx(idx)
          idx = idx + 1; qq13(i,j) = rcvx(idx)
        end do
      end do
    end if

    !==================================================================
    ! 3. Z-direction exchange (水平境界 + Xのゴースト分を含む)
    !==================================================================
    

    ! [A] Send to Down (ZM), Receive from Up (ZP)
    idx = 0
    do j = mg + 1, 2 * mg
      do i = 1, lix
        idx = idx + 1; sndz(idx) = qq1(i,j);  idx = idx + 1; sndz(idx) = qq2(i,j)
        idx = idx + 1; sndz(idx) = qq3(i,j);  idx = idx + 1; sndz(idx) = qq4(i,j)
        idx = idx + 1; sndz(idx) = qq5(i,j);  idx = idx + 1; sndz(idx) = qq6(i,j)
        idx = idx + 1; sndz(idx) = qq7(i,j);  idx = idx + 1; sndz(idx) = qq8(i,j)
        idx = idx + 1; sndz(idx) = qq9(i,j);  idx = idx + 1; sndz(idx) = qq10(i,j)
        idx = idx + 1; sndz(idx) = qq11(i,j); idx = idx + 1; sndz(idx) = qq12(i,j)
        idx = idx + 1; sndz(idx) = qq13(i,j)
      end do
    end do

    call MPI_Sendrecv(sndz, cntz, MPI_DOUBLE_PRECISION, rank_zm, 10, &
                      rcvz, cntz, MPI_DOUBLE_PRECISION, rank_zp, 10, &
                      MPI_COMM_WORLD, status, ierr)

    if (rank_zp /= MPI_PROC_NULL) then
      idx = 0
      do j = liz - mg + 1, liz
        do i = 1, lix
          idx = idx + 1; qq1(i,j) = rcvz(idx);  idx = idx + 1; qq2(i,j) = rcvz(idx)
          idx = idx + 1; qq3(i,j) = rcvz(idx);  idx = idx + 1; qq4(i,j) = rcvz(idx)
          idx = idx + 1; qq5(i,j) = rcvz(idx);  idx = idx + 1; qq6(i,j) = rcvz(idx)
          idx = idx + 1; qq7(i,j) = rcvz(idx);  idx = idx + 1; qq8(i,j) = rcvz(idx)
          idx = idx + 1; qq9(i,j) = rcvz(idx);  idx = idx + 1; qq10(i,j) = rcvz(idx)
          idx = idx + 1; qq11(i,j) = rcvz(idx); idx = idx + 1; qq12(i,j) = rcvz(idx)
          idx = idx + 1; qq13(i,j) = rcvz(idx)
        end do
      end do
    end if

    ! [B] Send to Up (ZP), Receive from Down (ZM)
    idx = 0
    do j = liz - 2 * mg + 1, liz - mg
      do i = 1, lix
        idx = idx + 1; sndz(idx) = qq1(i,j);  idx = idx + 1; sndz(idx) = qq2(i,j)
        idx = idx + 1; sndz(idx) = qq3(i,j);  idx = idx + 1; sndz(idx) = qq4(i,j)
        idx = idx + 1; sndz(idx) = qq5(i,j);  idx = idx + 1; sndz(idx) = qq6(i,j)
        idx = idx + 1; sndz(idx) = qq7(i,j);  idx = idx + 1; sndz(idx) = qq8(i,j)
        idx = idx + 1; sndz(idx) = qq9(i,j);  idx = idx + 1; sndz(idx) = qq10(i,j)
        idx = idx + 1; sndz(idx) = qq11(i,j); idx = idx + 1; sndz(idx) = qq12(i,j)
        idx = idx + 1; sndz(idx) = qq13(i,j)
      end do
    end do

    call MPI_Sendrecv(sndz, cntz, MPI_DOUBLE_PRECISION, rank_zp, 11, &
                      rcvz, cntz, MPI_DOUBLE_PRECISION, rank_zm, 11, &
                      MPI_COMM_WORLD, status, ierr)

    if (rank_zm /= MPI_PROC_NULL) then
      idx = 0
      do j = 1, mg
        do i = 1, lix
          idx = idx + 1; qq1(i,j) = rcvz(idx);  idx = idx + 1; qq2(i,j) = rcvz(idx)
          idx = idx + 1; qq3(i,j) = rcvz(idx);  idx = idx + 1; qq4(i,j) = rcvz(idx)
          idx = idx + 1; qq5(i,j) = rcvz(idx);  idx = idx + 1; qq6(i,j) = rcvz(idx)
          idx = idx + 1; qq7(i,j) = rcvz(idx);  idx = idx + 1; qq8(i,j) = rcvz(idx)
          idx = idx + 1; qq9(i,j) = rcvz(idx);  idx = idx + 1; qq10(i,j) = rcvz(idx)
          idx = idx + 1; qq11(i,j) = rcvz(idx); idx = idx + 1; qq12(i,j) = rcvz(idx)
          idx = idx + 1; qq13(i,j) = rcvz(idx)
        end do
      end do
    end if

    deallocate(sndx, rcvx, sndz, rcvz)

  end subroutine exc_2d_13var
end module exc
