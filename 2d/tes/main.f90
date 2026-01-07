program test_halo_exchange
    use mpi
    use vars  ! my_rank, mpx, mpz, mplx, mplz 等が定義されていると想定
    use exc
    implicit none

    ! 格子サイズとゴースト幅の設定
    integer, parameter :: lix = 12, liz = 12, mg = 2
    real(8) :: qq1(lix,liz),  qq2(lix,liz),  qq3(lix,liz),  &
               qq4(lix,liz),  qq5(lix,liz),  qq6(lix,liz),  &
               qq7(lix,liz),  qq8(lix,liz),  qq9(lix,liz),  &
               qq10(lix,liz), qq11(lix,liz), qq12(lix,liz), qq13(lix,liz)
    
    integer :: i, j, k, gi, gj
    logical :: success, global_success

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    ! 2x2 分割のテスト（最小4プロセス必要）
    !mpx = 2
    !mpz = nprocs / mpx
    mplx = mod(my_rank, mpx)
    mplz = my_rank / mpx

    !---------------------------------------------------------
    ! 1. データの初期化
    !---------------------------------------------------------
    ! 実計算領域のみに「座標依存の一意な値」を代入し、ゴーストは0で初期化
    call initialize_data(qq1, 1);  call initialize_data(qq2, 2)
    call initialize_data(qq3, 3);  call initialize_data(qq4, 4)
    call initialize_data(qq5, 5);  call initialize_data(qq6, 6)
    call initialize_data(qq7, 7);  call initialize_data(qq8, 8)
    call initialize_data(qq9, 9);  call initialize_data(qq10, 10)
    call initialize_data(qq11, 11); call initialize_data(qq12, 12)
    call initialize_data(qq13, 13)

    !---------------------------------------------------------
    ! 2. 通信ルーチンの実行
    !---------------------------------------------------------
    if (my_rank == 0) print *, "Starting halo exchange for 13 variables..."
    
    call exc_2d_13var( &
        qq1, qq2, qq3, qq4, qq5, qq6, qq7, &
        qq8, qq9, qq10, qq11, qq12, qq13, &
        lix, liz, mg)

    !---------------------------------------------------------
    ! 3. 検証（一例として qq1 の境界をチェック）
    !---------------------------------------------------------
    success = .true.
    ! 右側のランク(mplx+1)からデータが正しく届いているか確認
    if (mplx < mpx - 1) then
        ! 受信した右側ゴースト(i = lix-mg+1)の値をチェック
        i = lix - mg + 1
        j = liz / 2
        ! 期待値: 隣のランク(mplx+1)の左端の実計算領域(local i = mg+1)の値
        gi = (mplx + 1) * (lix - 2*mg) + 1
        gj = mplz * (liz - 2*mg) + (j - mg)
        if (abs(qq1(i, j) - (dble(gi * 1000 + gj) + 0.01d0)) > 1d-8) success = .false.
    end if

    ! 全ランクの成功判定を集計
    call MPI_ALLREDUCE(success, global_success, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)

    if (my_rank == 0) then
        if (global_success) then
            print *, "RESULT: SUCCESS! All variables exchanged correctly."
        else
            print *, "RESULT: FAILURE! Data mismatch in ghost cells."
        end if
    end if

    call MPI_FINALIZE(ierr)

contains

    subroutine initialize_data(q, k_idx)
        real(8), intent(out) :: q(lix, liz)
        integer, intent(in) :: k_idx
        integer :: i, j, gi, gj
        
        q = 0.0d0
        do j = mg + 1, liz - mg
            do i = mg + 1, lix - mg
                ! グローバルインデックスの計算
                gi = mplx * (lix - 2*mg) + (i - mg)
                gj = mplz * (liz - 2*mg) + (j - mg)
                ! 値 = (1000 * x座標 + z座標) + 変数番号/100
                q(i, j) = dble(gi * 1000 + gj) + dble(k_idx) / 100.0d0
            end do
        end do
    end subroutine initialize_data

end program test_halo_exchange
!program main
!  !--------------------------------------------------
!  ! include modules
!  !--------------------------------------------------
!  use vars ! for mpi
!  use datw
!  use const
!  ! use init
!  use model
!  use calc
!  use cfl
!  use mlw
!  use bnd2d
!  use arvis
!  use exc
!  use mpi
!  use,intrinsic::ieee_arithmetic,only:ieee_is_nan
!
!  !--------------------------------------------------
!  ! variables definitions
!  !--------------------------------------------------
!  implicit none
!  ! message
!  character(len=100)::foutmsg,fstpmsg
!
!  ! for mpi
!  !integer,parameter:: mpx=3
!  !integer,parameter:: mpz=2
!  !integer:: mpall=mpx*mpz
!  !integer:: mplx,mplz
!
!  ! for programm (time stepping et al.)
!  integer:: ns,nd
!  real(8):: t,dt,tend,tout,dtout
!  integer:: ii,jj
!
!  ! numerical variables (parameters)
!  integer,parameter:: mg=3
!  integer,parameter:: gx=144
!  integer,parameter:: gz=100
!  ! grid # for one mpi cells
!  integer,parameter:: lix=2*mg+gx/mpx
!  integer,parameter:: liz=2*mg+gz/mpz
!  ! grid # for global cells
!  integer,parameter:: gix=2*mg+gx
!  integer,parameter:: giz=2*mg+gz
!  ! specific heat
!  real(8),parameter:: gm=5.d0/3.d0
!
!  ! array definitions
!  real(8):: xx(lix),xm(lix),dx
!  real(8):: zz(liz),zm(liz),dz
!  ! physical variables
!  real(8):: ro(lix,liz),&
!    vx(lix,liz),vy(lix,liz),vz(lix,liz),&
!    bx(lix,liz),by(lix,liz),bz(lix,liz),&
!    ex(lix,liz),ey(lix,liz),ez(lix,liz),&
!    pr(lix,liz),etot(lix,liz),ay(lix,liz)
!  real(8):: ron(lix,liz),&
!    vxn(lix,liz),vyn(lix,liz),vzn(lix,liz),&
!    bxn(lix,liz),byn(lix,liz),bzn(lix,liz),&
!    exn(lix,liz),eyn(lix,liz),ezn(lix,liz),&
!    prn(lix,liz),etotn(lix,liz)
!  real(8):: rvx(lix,liz),rvy(lix,liz),rvz(lix,liz),&
!    rvxn(lix,liz),rvyn(lix,liz),rvzn(lix,liz)
!  real(8):: dro(lix,liz),detot(lix,liz),&
!    drvx(lix,liz),drvy(lix,liz),drvz(lix,liz),&
!    dbx(lix,liz),dby(lix,liz),dbz(lix,liz)
!  real(8):: grx(lix,liz),grz(lix,liz),&
!    grxm(lix,liz),grzm(lix,liz)
!  real(8):: v2(lix,liz),v2n(lix,liz),b2(lix,liz),b2n(lix,liz)
!
!  ! for mlw array
!  real(8):: qfx(lix,liz),qfz(lix,liz),qr(lix,liz)
!  !real(8):: qfxm(lix,liz),qfzm(lix,liz)! <- not needed?
!
!  ! for artificial viscosity
!  real(8),parameter:: qv=10.d0
!  real(8):: kx(lix,liz),kz(lix,liz)
!
!  ! for flooring
!  ! real(8),parameter::ro_floor=1.d-6
!  ! real(8),parameter::pr_floor=1.d-7
!
!  ! for file output
!  character(80):: dir
!
!  !--------------------------------------------------
!  ! setup for mpi (message passing interface)
!  !--------------------------------------------------
!  call MPI_INIT(ierr)
!  call MPI_COMM_SIZE(MPI_Comm_World,nprocs,ierr)
!  call MPI_COMM_RANK(MPI_Comm_World,my_rank,ierr)
!  ! check mpi number is consistent
!  if(my_rank==0) then
!    ! check number of mpi
!    if(nprocs/=mpx*mpz) then
!      write(*,*) "number of processes is not consistent"
!      write(*,*) "change mpx and mpz in vars.f90 or MPINUM in makefile"
!      write(*,*) "nprocs= ", nprocs
!      write(*,*) "(mpx,mpz)=", mpx,mpz
!      call MPI_FINALIZE(ierr)
!      stop
!    end if
!    ! check if grid size can be divided by mpi number
!    if(mod(gx,mpx)+mod(gz,mpz)/=0) then
!      write(*,*) "gx/mpix or gz/mpiz is not integer"
!      write(*,*) "gx=", gx
!      write(*,*) "mpx= ", mpx
!      write(*,*) "gz=", gz
!      write(*,*) "mpx= ", mpz
!      call MPI_FINALIZE(ierr)
!      stop
!    end if
!  end if
!
!  !--------------------------------------------------
!  ! setup and initialize
!  !--------------------------------------------------
!  ! make terminal output format
!  !--------------------------------------------------
!  foutmsg='(1x," WRITE  ","STEP=  ",i7," time=",e10.3," nd=",i4)'
!  fstpmsg='(1x," STOP  ","STEP=  ",i7," time=",e10.3)'
!  !--------------------------------------------------
!  ! initialize variables
!  !--------------------------------------------------
!  ro=0.d0;vx=0.d0;vy=0.d0;vz=0.d0;bx=0.d0;by=0.d0;bz=0.d0
!  ex=0.d0;ey=0.d0;ez=0.d0;pr=0.d0;etot=0.d0;ay=0.d0
!  ron=0.d0;vxn=0.d0;vyn=0.d0;vzn=0.d0;bxn=0.d0;byn=0.d0;bzn=0.d0
!  exn=0.d0;eyn=0.d0;ezn=0.d0;prn=0.d0;etotn=0.d0
!  rvx=0.d0;rvy=0.d0;rvz=0.d0;rvxn=0.d0;rvyn=0.d0;rvzn=0.d0
!  dro=0.d0;detot=0.d0;drvx=0.d0;drvy=0.d0;drvz=0.d0
!  dbx=0.d0;dby=0.d0;dbz=0.d0
!  v2=0.d0;v2n=0.d0;b2=0.d0;b2n=0.d0
!  kx=0.d0;kz=0.d0
!  !--------------------------------------------------
!  ! clean past data and initialize parameter files
!  ! (dat/*.dat and param.txt)
!  !--------------------------------------------------
!  if(my_rank==0) then
!    call dataclean()
!    call put_param_int("mgn:",mg)
!    call put_param_int("gix:",gix)
!    call put_param_int("giz:",giz)
!    call put_param_int("mpx:",mpx)
!    call put_param_int("mpz:",mpz)
!    call put_param_real("specific heat:",gm)
!  end if
!  ! call MPI_BARRIER(MPI_Comm_World,ierr)
!  !--------------------------------------------------
!  ! time paramters
!  ! ns = # of steps; nd = # of written data
!  !--------------------------------------------------
!  tend=30.d0;dtout=1.d0
!  t=0.d0;tout=0.d0
!  ns=0;nd=0
!  !--------------------------------------------------
!  ! setup model
!  !--------------------------------------------------
!  call make_grid(xx,xm,dx,zz,zm,dz,gix,lix,giz,liz,mg)
!  !call model_diskjet(ro,vx,vy,vz,bx,by,bz,pr,&
!    !gm,grx,grxm,grz,grzm,xx,zz,lix,liz,dx,dz,mg)
!  !call calc_ay(ay,bz,xx,zz,lix,liz)
!  !call calc_efield(bx,by,bz,vx,vy,vz,ex,ey,ez)
!  !call calc_etot(ro,pr,vx,vy,vz,bx,by,bz,etot,gm)
!  !--------------------------------------------------
!  ! output (todo module (subroutine)) 1st time
!  !--------------------------------------------------
!  write(dir,'(A,i3.3,"_")') 'dat/', nd
!  if(my_rank==0) then
!    call put_time_data('dat/t.dat',t)
!  end if
!  !call put_2d_data_each_rank(&
!  !  trim(dir)//'ro',lix,liz,mg,ro)
!  !call put_2d_data_each_rank(&
!  !  trim(dir)//'vx',lix,liz,mg,vx)
!  !call put_2d_data_each_rank(&
!  !  trim(dir)//'vy',lix,liz,mg,vy)
!  !call put_2d_data_each_rank(&
!  !  trim(dir)//'vz',lix,liz,mg,vz)
!  !call put_2d_data_each_rank(&
!  !  trim(dir)//'bx',lix,liz,mg,bx)
!  !call put_2d_data_each_rank(&
!  !  trim(dir)//'by',lix,liz,mg,by)
!  !call put_2d_data_each_rank(&
!  !  trim(dir)//'bz',lix,liz,mg,bz)
!  !call put_2d_data_each_rank(&
!  !  trim(dir)//'pr',lix,liz,mg,pr)
!  !call put_2d_data_each_rank(&
!  !  trim(dir)//'ay',lix,liz,mg,ay)
!  !if(my_rank==0) then
!  !  write(*,foutmsg) ns,t,nd
!  !end if
!  !nd=nd+1
!
!
!  !--------------------------------------------------
!  ! initialize
!  !--------------------------------------------------
!  ro=0.d0;vx=0.d0;vy=0.d0;vz=0.d0;bx=0.d0;by=0.d0;bz=0.d0
!  ex=0.d0;ey=0.d0;ez=0.d0;pr=0.d0;etot=0.d0;ay=0.d0
!  ron=0.d0;vxn=0.d0;vyn=0.d0;vzn=0.d0;bxn=0.d0;byn=0.d0;bzn=0.d0
!  exn=0.d0;eyn=0.d0;ezn=0.d0;prn=0.d0;etotn=0.d0
!  rvx=0.d0;rvy=0.d0;rvz=0.d0;rvxn=0.d0;rvyn=0.d0;rvzn=0.d0
!  dro=0.d0;detot=0.d0;drvx=0.d0;drvy=0.d0;drvz=0.d0
!  dbx=0.d0;dby=0.d0;dbz=0.d0
!  v2=0.d0;v2n=0.d0;b2=0.d0;b2n=0.d0
!  kx=0.d0;kz=0.d0
!
!  
!
!  ro(:,:) = -10.d0
!
!  do jj = mg+1, liz-mg
!  do ii = mg+1, lix-mg
!    ro(ii,jj) = dble(my_rank)
!  end do
!  end do
!
!  call exc_2d_13var( &
!    ro,ro,ro,ro,ro,ro,ro, &
!    ro,ro,ro,ro,ro,ro, &
!   lix,liz,mg)
!
!  call put_2d_data_each_rank('dat/ro',lix,liz,mg,ro)
!  !--------------------------------------------------
!  ! ending message
!  !--------------------------------------------------
!  call MPI_FINALIZE(ierr)
!
!
!
!end program main
