program main
  !--------------------------------------------------
  ! include modules
  !--------------------------------------------------
  use vars ! for mpi
  use datw
  use const
  ! use init
  use model
  use calc
  use cfl
  ! use mlw
  use bnd2d
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
  !integer,parameter:: mpx=3
  !integer,parameter:: mpz=2
  !integer:: mpall=mpx*mpz
  !integer:: mplx,mplz

  ! for programm (time stepping et al.)
  integer:: ns,nd
  real(8):: t,dt,tend,tout,dtout
  integer:: ii,jj

  ! numerical variables (parameters)
  integer,parameter:: mg=2
  integer,parameter:: gx=600
  integer,parameter:: gz=1200
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
  real(8):: grx(lix,liz),grz(lix,liz),&
    grxm(lix,liz),grzm(lix,liz)
  real(8):: v2(lix,liz),v2n(lix,liz),b2(lix,liz),b2n(lix,liz)

    ! for mlw array
  real(8):: qu(lix,liz),qfx(lix,liz),qfz(lix,liz),qr(lix,liz)
  real(8):: qfxm(lix,liz),qfzm(lix,liz)! <- not needed?

    ! for artificial viscosity
  real(8),parameter:: qv=10.d0
  real(8):: kx(lix,liz),kz(lix,liz)

  ! for flooring
  ! real(8),parameter::ro_floor=1.d-6
  ! real(8),parameter::pr_floor=1.d-7

  ! for file output
  character(80):: dir

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
      write(*,*) "change mpx and mpz in vars.f90 or MPINUM in makefile"
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
  call model_diskjet(ro,vx,vy,vz,bx,by,bz,pr,&
    gm,grx,grxm,grz,grzm,xx,zz,lix,liz,dx,dz,mg)
  call calc_ay(ay,bz,xx,zz,lix,liz)
  call calc_efield(bx,by,bz,vx,vy,vz,ex,ey,ez)
  call calc_etot(ro,pr,vx,vy,vz,bx,by,bz,etot,gm)
  !--------------------------------------------------
  ! output (todo module (subroutine)) 1st time
  !--------------------------------------------------
  write(dir,'(A,i3.3,"_")') 'dat/', nd
  if(my_rank==0) then
    call put_time_data('dat/t.dat',t)
  end if
  call put_2d_data_each_rank(&
    trim(dir)//'ro',lix,liz,mg,ro)
  call put_2d_data_each_rank(&
    trim(dir)//'vx',lix,liz,mg,vx)
  call put_2d_data_each_rank(&
    trim(dir)//'vy',lix,liz,mg,vy)
  call put_2d_data_each_rank(&
    trim(dir)//'vz',lix,liz,mg,vz)
  call put_2d_data_each_rank(&
    trim(dir)//'bx',lix,liz,mg,bx)
  call put_2d_data_each_rank(&
    trim(dir)//'by',lix,liz,mg,by)
  call put_2d_data_each_rank(&
    trim(dir)//'bz',lix,liz,mg,bz)
  call put_2d_data_each_rank(&
    trim(dir)//'pr',lix,liz,mg,pr)
  call put_2d_data_each_rank(&
    trim(dir)//'ay',lix,liz,mg,ay)
  if(my_rank==0) then
    write(*,foutmsg) ns,t,nd
  end if
  nd=nd+1
  tout = tout+dtout

  !--------------------------------------------------
  ! main loop (time integration)
  !--------------------------------------------------
  main_loop:&
    do while(t<tend)
  !--------------------------------------------------
  ! determine time interval using CFL cond.
  !--------------------------------------------------
  call calc_dt(ro,vx,vy,vz,pr,bx,by,bz,gm,lix,liz,dt,dx,dz,mg)

  !--------------------------------------------------
  ! for debug
  !--------------------------------------------------
  call MPI_FINALIZE(ierr)
  stop

  end do main_loop





  !--------------------------------------------------
  ! ending message
  !--------------------------------------------------
  call MPI_FINALIZE(ierr)



end program main
