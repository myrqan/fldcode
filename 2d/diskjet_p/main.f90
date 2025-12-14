program main
  !--------------------------------------------------
  ! include modules
  !--------------------------------------------------
  use vars ! for mpi
  use datw
  ! use init
  use model
  ! use datw
  ! use cfl
  ! use mlw
  ! use bnd2d
  ! use arvis
  ! use exc
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_is_nan

  !--------------------------------------------------
  ! variables definitions
  !--------------------------------------------------
  implicit none
  ! message
  character(len=100)::foutmsg,fstpmsg

  ! for mpi
  integer,parameter:: mpx=3
  integer,parameter:: mpz=2
  integer:: mpall=mpx*mpz
  !integer:: mplx,mplz

  ! for programm (time stepping et al.)
  integer:: ns,nd
  real(8):: t,dt,tend,tout,dtout
  integer:: ii,jj

  ! numerical variables (parameters)
  integer,parameter:: mg=2
  integer,parameter:: gx=300
  integer,parameter:: gz=200
    ! grid # for one mpi cells
  integer,parameter:: lix=2*mg+gx/mpx
  integer,parameter:: liz=2*mg+gz/mpz
    ! grid # for global cells
  integer,parameter:: gix=2*mg+gx
  integer,parameter:: giz=2*mg+gz
    ! specific heat
  real(8),parameter:: gm=5.d0/3.d0

  ! array definitions
  real(8):: xx(lix),xm(lix),dx
  real(8):: zz(liz),zm(liz),dz
    ! physical variables
  real(8):: ro(lix,liz),&
    vx(lix,liz),vy(lix,liz),vz(lix,liz),&
    bx(lix,liz),by(lix,liz),bz(lix,liz),&
    ex(lix,liz),ey(lix,liz),ez(lix,liz),&
    pr(lix,liz),etot(lix,liz),ay(lix,liz)
  real(8):: ron(lix,liz),&
    vxn(lix,liz),vyn(lix,liz),vzn(lix,liz),&
    bxn(lix,liz),byn(lix,liz),bzn(lix,liz),&
    exn(lix,liz),eyn(lix,liz),ezn(lix,liz),&
    prn(lix,liz),etotn(lix,liz)
  real(8):: rvx(lix,liz),rvy(lix,liz),rvz(lix,liz),&
    rvxn(lix,liz),rvyn(lix,liz),rvzn(lix,liz)
  real(8):: dro(lix,liz),detot(lix,liz),&
    drvx(lix,liz),drvy(lix,liz),drvz(lix,liz),&
    dbx(lix,liz),dby(lix,liz),dbz(lix,liz)
  real(8):: v2(lix,liz),v2n(lix,liz),b2(lix,liz),b2n(lix,liz)

    ! for mlw array
  real(8):: qu(lix,liz),qfx(lix,liz),qfz(lix,liz),qr(lix,liz)
  real(8):: qfxm(lix,liz),qfzm(lix,liz)! <- not needed?

    ! for artificial viscosity
  real(8),parameter:: qv=10.d0
  real(8):: kx(lix,liz),kz(lix,liz)

  ! mathematical constants
  real(8),parameter:: pi=4.d0*datan(1.d0)
  
  ! for flooring
  ! real(8),parameter::ro_floor=1.d-6
  ! real(8),parameter::pr_floor=1.d-7

  ! for file output
  integer:: ifxx,ifzz,ifro,ifvx,ifvy,ifvz,ifbx,ifby,ifbz,ifpr
  integer:: write_cntx,write_cntz,write_cntxz
  integer(kind=MPI_OFFSET_KIND):: dispx,dispz,dispxz,&
    offsetx,offsetz,offsetxz,file_end_pos
  integer:: status(MPI_STATUS_SIZE)

  !--------------------------------------------------
  ! setup for mpi (message passing interface)
  !--------------------------------------------------
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_Comm_World,nprocs,ierr)
  call MPI_COMM_RANK(MPI_Comm_World,my_rank,ierr)
    ! check mpi number is consistent
  if(my_rank==0) then
    ! check number of mpi
    if(nprocs/=mpx*mpz) then
      write(*,*) "number of processes is not consistent"
      write(*,*) "change mpx and mpz in main.f90 or MPINUM in makefile"
      write(*,*) "nprocs= ", nprocs
      write(*,*) "(mpx,mpz)=", mpx,mpz
      call MPI_FINALIZE(ierr)
      stop
    end if
    ! check if grid size can be divided by mpi number
    if(mod(gx,mpx)+mod(gz,mpz)/=0) then
      write(*,*) "gx/mpix or gz/mpiz is not integer"
      write(*,*) "gx=", gx
      write(*,*) "mpx= ", mpx
      write(*,*) "gz=", gz
      write(*,*) "mpx= ", mpz
      call MPI_FINALIZE(ierr)
      stop
    end if
  end if

  ! to determine file position when outputting data
  !write_cnt=(lix-2*mg)*(liz-2*mg) ! output data size (number)
  !disp = int(write_cnt,kind=MPI_OFFSET_KIND)*8_MPI_OFFSET_KIND
  !offset=int(my_rank,kind=MPI_OFFSET_KIND)*disp
  !--------------------------------------------------
  ! setup and initialize
  !--------------------------------------------------
  ! make terminal output format
  !--------------------------------------------------
  foutmsg='(1x," WRITE  ","STEP=  ",i7," time=",e10.3," nd=",i4)'
  fstpmsg='(1x," STOP  ","STEP=  ",i7," time=",e10.3)'
  !--------------------------------------------------
  ! initialize variables
  !--------------------------------------------------
  ro=0.d0;vx=0.d0;vy=0.d0;vz=0.d0;bx=0.d0;by=0.d0;bz=0.d0
  ex=0.d0;ey=0.d0;ez=0.d0;pr=0.d0;etot=0.d0;ay=0.d0
  ron=0.d0;vxn=0.d0;vyn=0.d0;vzn=0.d0;bxn=0.d0;byn=0.d0;bzn=0.d0
  exn=0.d0;eyn=0.d0;ezn=0.d0;prn=0.d0;etotn=0.d0
  rvx=0.d0;rvy=0.d0;rvz=0.d0;rvxn=0.d0;rvyn=0.d0;rvzn=0.d0
  dro=0.d0;detot=0.d0;drvx=0.d0;drvy=0.d0;drvz=0.d0
  dbx=0.d0;dby=0.d0;dbz=0.d0
  v2=0.d0;v2n=0.d0;b2=0.d0;b2n=0.d0
  kx=0.d0;kz=0.d0
  !--------------------------------------------------
  ! clean past data and initialize parameter files
  ! (dat/*.dat and param.txt)
  !--------------------------------------------------
  if(my_rank==0) then
    call dataclean()
    call put_param_int("mgn:",mg)
    call put_param_int("gix:",gix)
    call put_param_int("giz:",giz)
    call put_param_int("mpx:",mpx)
    call put_param_int("mpz:",mpz)
    call put_param_real("specific heat:",gm)
  end if
  call MPI_BARRIER(MPI_Comm_World,ierr)
  !--------------------------------------------------
  ! time paramters
  ! ns = # of steps; nd = # of written data
  !--------------------------------------------------
  tend=30.d0;dtout=1.d0
  t=0.d0;tout=0.d0
  ns=0;nd=0
  !--------------------------------------------------
  ! setup model
  !--------------------------------------------------
  call make_grid(xx,xm,dx,zz,zm,dz,gix,lix,giz,liz,mg)
  ! todo
  !call model_diskjet
  !call calculate_ay

  !call apply_bnd

  !call data_output
  if(my_rank==0) then
    write(*,foutmsg) ns,t,nd
  end if




  !--------------------------------------------------
  ! ending message
  !--------------------------------------------------
  call MPI_FINALIZE(ierr)



end program main
