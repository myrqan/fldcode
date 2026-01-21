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
  use mlw
  use bnd2d
  use arvis
  use exc
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
  integer,parameter:: mg=4
  integer,parameter:: gx=500
  integer,parameter:: gz=500
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
  real(8):: qfx(lix,liz),qfz(lix,liz),qr(lix,liz)
  !real(8):: qfxm(lix,liz),qfzm(lix,liz)! <- not needed?

    ! for artificial viscosity
  real(8),parameter:: qv=10.d0
  real(8):: kx(lix,liz),kz(lix,liz)

  ! for flooring
  real(8),parameter::ro_floor=1.d-6
  real(8),parameter::pr_floor=1.d-7

  ! for file output
  character(80):: dir

  !--------------------------------------------------
  ! setup for mpi (message passing interface)
  !--------------------------------------------------
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_Comm_World,nprocs,ierr)
  call MPI_COMM_RANK(MPI_Comm_World,my_rank,ierr)
    ! check mpi number is consistent
    ! check number of mpi
    if(nprocs/=mpx*mpz) then
      if(my_rank==0) then
        write(*,*) "number of processes is not consistent"
        write(*,*) "change mpx and mpz in vars.f90 or MPINUM in makefile"
        write(*,*) "nprocs= ", nprocs
        write(*,*) "(mpx,mpz)=", mpx,mpz
      end if
      call MPI_FINALIZE(ierr)
      stop
    end if
    ! check if grid size can be divided by mpi number
    if(mod(gx,mpx)+mod(gz,mpz)/=0) then
      if(my_rank==0) then
        write(*,*) "gx/mpix or gz/mpiz is not integer"
        write(*,*) "gx=", gx
        write(*,*) "mpx= ", mpx
        write(*,*) "gz=", gz
        write(*,*) "mpx= ", mpz
      end if
      call MPI_FINALIZE(ierr)
      stop
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
  ! call MPI_BARRIER(MPI_Comm_World,ierr)
  !--------------------------------------------------
  ! time paramters
  ! ns = # of steps; nd = # of written data
  !--------------------------------------------------
  tend=30.d0;
  dtout = 1.d0
  !dtout=0.d0 ! todo changed
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
  ! solve with modified two-step Lax-Wendroff scheme
  !--------------------------------------------------
  main_loop:&
    do while(t<tend)
  !--------------------------------------------------
  ! determine time interval using CFL cond.
  !--------------------------------------------------
  call calc_dt(ro,vx,vy,vz,pr,bx,by,bz,gm,lix,liz,dt,dx,dz,mg)



  !--------------------------------------------------
  ! initialize variables
  !--------------------------------------------------
  ron=0.d0;rvxn=0.d0;rvyn=0.d0;rvzn=0.d0
  bxn=0.d0;byn=0.d0;bzn=0.d0
  exn=0.d0;eyn=0.d0;ezn=0.d0;etotn=0.d0

  dro=0.d0;detot=0.d0;drvx=0.d0;drvy=0.d0;drvz=0.d0
  dbx=0.d0;dby=0.d0;dbz=0.d0

  rvx(:,:) = ro(:,:)*vx(:,:)
  rvy(:,:) = ro(:,:)*vy(:,:)
  rvz(:,:) = ro(:,:)*vz(:,:)
  b2(:,:) = bx(:,:)**2+by(:,:)**2+bz(:,:)**2
  v2(:,:) = vx(:,:)**2+vy(:,:)**2+vz(:,:)**2

  !--------------------------------------------------
  ! 1st step
  !--------------------------------------------------
  ! mass conservation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = rvx(ii,jj)
    qfz(ii,jj) = rvz(ii,jj)
    qr(ii,jj) = -qfx(ii,jj)/xx(ii)
  end do
  end do
  call mlw2d1st(ro,qfx,qfz,ron,dro,lix,liz,dt,dx,dz)
  call mlw2dsrc1st(qr,ron,dro,lix,liz,dt,dx,dz)

  ! x-momentum conservation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = ro(ii,jj)*vx(ii,jj)**2+pr(ii,jj)&
    +0.5d0*(-bx(ii,jj)**2+by(ii,jj)**2+bz(ii,jj)**2)
    qfz(ii,jj) = ro(ii,jj)*vz(ii,jj)*vx(ii,jj)-bz(ii,jj)*bx(ii,jj)
    qr(ii,jj) = -(ro(ii,jj)*(vx(ii,jj)**2-vy(ii,jj)**2) &
      -(bx(ii,jj)**2-by(ii,jj)**2))/xx(ii) + ro(ii,jj)*grx(ii,jj)
  end do
  end do
  call mlw2d1st(rvx,qfx,qfz,rvxn,drvx,lix,liz,dt,dx,dz)
  call mlw2dsrc1st(qr,rvxn,drvx,lix,liz,dt,dx,dz)

  ! y-momentum conservation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = ro(ii,jj)*vx(ii,jj)*vy(ii,jj)-bx(ii,jj)*by(ii,jj)
    qfz(ii,jj) = ro(ii,jj)*vz(ii,jj)*vy(ii,jj)-bz(ii,jj)*by(ii,jj)
    qr(ii,jj) = -2.d0*qfx(ii,jj)/xx(ii)
  end do
  end do
  call mlw2d1st(rvy,qfx,qfz,rvyn,drvy,lix,liz,dt,dx,dz)
  call mlw2dsrc1st(qr,rvyn,drvy,lix,liz,dt,dx,dz)

  ! z-momentum conservation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = ro(ii,jj)*vx(ii,jj)*vz(ii,jj)-bx(ii,jj)*bz(ii,jj)
    qfz(ii,jj) = ro(ii,jj)*vz(ii,jj)**2+pr(ii,jj)&
      +0.5d0*(bx(ii,jj)**2+by(ii,jj)**2-bz(ii,jj)**2)
    qr(ii,jj) = - qfx(ii,jj)/xx(ii) + ro(ii,jj)*grz(ii,jj)
  end do
  end do
  call mlw2d1st(rvz,qfx,qfz,rvzn,drvz,lix,liz,dt,dx,dz)
  call mlw2dsrc1st(qr,rvzn,drvz,lix,liz,dt,dx,dz)

  ! energy conservation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = (gm/(gm-1.d0)*pr(ii,jj) + 0.5d0*ro(ii,jj)*v2(ii,jj))*vx(ii,jj)& 
      + (ey(ii,jj)*bz(ii,jj)-ez(ii,jj)*by(ii,jj))
    qfz(ii,jj) = (gm/(gm-1.d0)*pr(ii,jj) + 0.5d0*ro(ii,jj)*v2(ii,jj))*vz(ii,jj)&
      + (ex(ii,jj)*by(ii,jj)-ey(ii,jj)*bx(ii,jj))
    qr(ii,jj) = -qfx(ii,jj)/xx(ii) &
      + ro(ii,jj)*grx(ii,jj)*vx(ii,jj)+ro(ii,jj)*grz(ii,jj)*vz(ii,jj)
  end do
  end do
  call mlw2d1st(etot,qfx,qfz,etotn,detot,lix,liz,dt,dx,dz)
  call mlw2dsrc1st(qr,etotn,detot,lix,liz,dt,dx,dz)

  ! x-induction equation
  qfx(:,:)=0.d0
  qfz(:,:)=-ey(:,:)
  qr(:,:)=0.d0
  call mlw2d1st(bx,qfx,qfz,bxn,dbx,lix,liz,dt,dx,dz)
  !call mlw2dsrc1st()

  ! y-induction equation
  qfx(:,:) = -ez(:,:)
  qfz(:,:) = ex(:,:)
  qr(:,:) = 0.d0
  call mlw2d1st(by,qfx,qfz,byn,dby,lix,liz,dt,dx,dz)
  !call mlw2dsrc1st()

  ! z-induction equation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = ey(ii,jj)
    qfz(ii,jj) = 0.d0
    qr(ii,jj) = -qfx(ii,jj)/xx(ii)
  end do
  end do
  call mlw2d1st(bz,qfx,qfz,bzn,dbz,lix,liz,dt,dx,dz)
  call mlw2dsrc1st(qr,bzn,dbz,lix,liz,dt,dx,dz)

  ! get half step value
  vxn(:,:) = rvxn(:,:)/ron(:,:)
  vyn(:,:) = rvyn(:,:)/ron(:,:)
  vzn(:,:) = rvzn(:,:)/ron(:,:)
  v2n(:,:) = vxn(:,:)**2+vyn(:,:)**2+vzn(:,:)**2
  b2n(:,:) = bxn(:,:)**2+byn(:,:)**2+bzn(:,:)**2
  call calc_pr(ron,prn,vxn,vyn,vzn,bxn,byn,bzn,etotn,gm,lix,liz)
  call calc_efield(bxn,byn,bzn,vxn,vyn,vzn,exn,eyn,ezn)

  !! check if negative pressure / density
  ! re-calculate total energy (after flooring pressure and density)A
  ! calc_etot()


  !--------------------------------------------------
  ! 2nd step
  !--------------------------------------------------
  ! mass conservation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = rvxn(ii,jj)
    qfz(ii,jj) = rvzn(ii,jj)
    qr(ii,jj) = -qfx(ii,jj)/xm(ii)
  end do
  end do
  call mlw2d2nd(qfx,qfz,dro,lix,liz,dt,dx,dz)
  call mlw2dsrc2nd(dro,qr,lix,liz,dt,dx,dz)

  ! x-momentum conservation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = ron(ii,jj)*vxn(ii,jj)**2+prn(ii,jj)&
    +0.5d0*(-bxn(ii,jj)**2+byn(ii,jj)**2+bzn(ii,jj)**2)
    qfz(ii,jj) = ron(ii,jj)*vzn(ii,jj)*vxn(ii,jj)-bzn(ii,jj)*bxn(ii,jj)
    qr(ii,jj) = -(ron(ii,jj)*(vxn(ii,jj)**2-vyn(ii,jj)**2) &
      -(bxn(ii,jj)**2-byn(ii,jj)**2))/xm(ii) + ron(ii,jj)*grxm(ii,jj)
  end do
  end do
  call mlw2d2nd(qfx,qfz,drvx,lix,liz,dt,dx,dz)
  call mlw2dsrc2nd(drvx,qr,lix,liz,dt,dx,dz)

  ! y-momentum conservation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = ron(ii,jj)*vxn(ii,jj)*vyn(ii,jj)-bxn(ii,jj)*byn(ii,jj)
    qfz(ii,jj) = ron(ii,jj)*vzn(ii,jj)*vyn(ii,jj)-bzn(ii,jj)*byn(ii,jj)
    qr(ii,jj) = -2.d0*qfx(ii,jj)/xm(ii)
  end do
  end do
  call mlw2d2nd(qfx,qfz,drvy,lix,liz,dt,dx,dz)
  call mlw2dsrc2nd(drvy,qr,lix,liz,dt,dx,dz)

  ! z-momentum conservation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = ron(ii,jj)*vxn(ii,jj)*vzn(ii,jj)-bxn(ii,jj)*bzn(ii,jj)
    qfz(ii,jj) = ron(ii,jj)*vzn(ii,jj)**2+prn(ii,jj)&
      +0.5d0*(bxn(ii,jj)**2+byn(ii,jj)**2-bzn(ii,jj)**2)
    qr(ii,jj) = - qfx(ii,jj)/xm(ii) + ron(ii,jj)*grzm(ii,jj)
  end do
  end do
  call mlw2d2nd(qfx,qfz,drvz,lix,liz,dt,dx,dz)
  call mlw2dsrc2nd(drvz,qr,lix,liz,dt,dx,dz)

  ! energy conservation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = (gm/(gm-1.d0)*prn(ii,jj) + 0.5d0*ron(ii,jj)*v2n(ii,jj))*vxn(ii,jj)& 
      + (eyn(ii,jj)*bzn(ii,jj)-ezn(ii,jj)*byn(ii,jj))
    qfz(ii,jj) = (gm/(gm-1.d0)*prn(ii,jj) + 0.5d0*ron(ii,jj)*v2n(ii,jj))*vzn(ii,jj)&
      + (exn(ii,jj)*byn(ii,jj)-eyn(ii,jj)*bxn(ii,jj))
    qr(ii,jj) = -qfx(ii,jj)/xm(ii) &
      + ron(ii,jj)*grxm(ii,jj)*vxn(ii,jj)+ron(ii,jj)*grzm(ii,jj)*vzn(ii,jj)
  end do
  end do
  call mlw2d2nd(qfx,qfz,detot,lix,liz,dt,dx,dz)
  call mlw2dsrc2nd(detot,qr,lix,liz,dt,dx,dz)

  ! x-induction equation
  qfx(:,:)=0.d0
  qfz(:,:)=-eyn(:,:)
  qr(:,:)=0.d0
  call mlw2d2nd(qfx,qfz,dbx,lix,liz,dt,dx,dz)
  !call mlw2dsrc2nd()

  ! y-induction equation
  qfx(:,:) = -ezn(:,:)
  qfz(:,:) = exn(:,:)
  qr(:,:) = 0.d0
  call mlw2d2nd(qfx,qfz,dby,lix,liz,dt,dx,dz)
  !call mlw2dsrc2nd()

  ! z-induction equation
  do jj=1,liz
  do ii=1,lix
    qfx(ii,jj) = eyn(ii,jj)
    qfz(ii,jj) = 0.d0
    qr(ii,jj) = -qfx(ii,jj)/xm(ii)
  end do
  end do
  call mlw2d2nd(qfx,qfz,dbz,lix,liz,dt,dx,dz)
  call mlw2dsrc2nd(dbz,qr,lix,liz,dt,dx,dz)

  !! check if negative pressure / density
  ! re-calculate total energy (after flooring pressure and density)
  !call calc_etot()

  !--------------------------------------------------
  ! apply atrificial viscosity
  !--------------------------------------------------
  if(ns==0) then
    call put_param_real("qv:",qv)
  end if

  ! calculate coefficients kx, kz
  call calc_coef_k(qv,kx,kz,vx,vz,lix,liz,dx,dz)

  ! apply artificial viscosity
  call arvis2d(kx,kz,ro,dro,lix,liz,dt,dx,dz)
  call arvis2d(kx,kz,rvx,drvx,lix,liz,dt,dx,dz)
  call arvis2d(kx,kz,rvy,drvy,lix,liz,dt,dx,dz)
  call arvis2d(kx,kz,rvz,drvz,lix,liz,dt,dx,dz)
  call arvis2d(kx,kz,etot,detot,lix,liz,dt,dx,dz)
  call arvis2d(kx,kz,bx,dbx,lix,liz,dt,dx,dz)
  call arvis2d(kx,kz,by,dby,lix,liz,dt,dx,dz)
  call arvis2d(kx,kz,bz,dbz,lix,liz,dt,dx,dz)


  !--------------------------------------------------
  ! update values
  !--------------------------------------------------
  ro(:,:)=ro(:,:)+dro(:,:)
  rvx(:,:)=rvx(:,:)+drvx(:,:)
  rvy(:,:)=rvy(:,:)+drvy(:,:)
  rvz(:,:)=rvz(:,:)+drvz(:,:)
  vx(:,:)=rvx(:,:)/ro(:,:)
  vy(:,:)=rvy(:,:)/ro(:,:)
  vz(:,:)=rvz(:,:)/ro(:,:)
  etot(:,:)=etot(:,:)+detot(:,:)
  bx(:,:)=bx(:,:)+dbx(:,:)
  by(:,:)=by(:,:)+dby(:,:)
  bz(:,:)=bz(:,:)+dbz(:,:)


  !--------------------------------------------------
  ! exchange values each mpi cells
  !--------------------------------------------------
  ! todo
  if  (nprocs > 1) then
    call exc_2d_13var(&
        ro,vx,vy,vz,bx,by,bz,&
        etot,pr,ex,ey,ez,ay,&
        lix,liz,mg)
  end if

  !--------------------------------------------------
  ! apply boundary condition (todo : module)
  !--------------------------------------------------
  call bnd_2d_symm(ro,0,mg,lix,liz)
  call bnd_2d_free(ro,1,mg,lix,liz)
  call bnd_2d_free(ro,2,mg,lix,liz)
  call bnd_2d_symm(ro,3,mg,lix,liz)

  call bnd_2d_symm(vx,0,mg,lix,liz)
  call bnd_2d_free(vx,1,mg,lix,liz)
  call bnd_2d_free(vx,2,mg,lix,liz)
  call bnd_2d_asym(vx,3,mg,lix,liz)

  call bnd_2d_symm(vy,0,mg,lix,liz)
  call bnd_2d_free(vy,1,mg,lix,liz)
  call bnd_2d_free(vy,2,mg,lix,liz)
  call bnd_2d_asym(vy,3,mg,lix,liz)

  call bnd_2d_asym(vz,0,mg,lix,liz)
  call bnd_2d_free(vz,1,mg,lix,liz)
  call bnd_2d_free(vz,2,mg,lix,liz)
  call bnd_2d_symm(vz,3,mg,lix,liz)

  call bnd_2d_asym(bx,0,mg,lix,liz)
  call bnd_2d_free(bx,1,mg,lix,liz)
  call bnd_2d_free(bx,2,mg,lix,liz)
  call bnd_2d_asym(bx,3,mg,lix,liz)

  call bnd_2d_asym(by,0,mg,lix,liz)
  call bnd_2d_free(by,1,mg,lix,liz)
  call bnd_2d_free(by,2,mg,lix,liz)
  call bnd_2d_asym(by,3,mg,lix,liz)

  call bnd_2d_symm(bz,0,mg,lix,liz)
  call bnd_2d_free(bz,1,mg,lix,liz)
  call bnd_2d_free(bz,2,mg,lix,liz)
  call bnd_2d_symm(bz,3,mg,lix,liz)

  call bnd_2d_symm(pr,0,mg,lix,liz)
  call bnd_2d_free(pr,1,mg,lix,liz)
  call bnd_2d_free(pr,2,mg,lix,liz)
  call bnd_2d_symm(pr,3,mg,lix,liz)

  call bnd_2d_symm(etot,0,mg,lix,liz)
  call bnd_2d_free(etot,1,mg,lix,liz)
  call bnd_2d_free(etot,2,mg,lix,liz)
  call bnd_2d_symm(etot,3,mg,lix,liz)
  

  call calc_pr(ro,pr,vx,vy,vz,bx,by,bz,etot,gm,lix,liz)
  call calc_efield(bx,by,bz,vx,vy,vz,ex,ey,ez)
  ! magnetic field line
  ay(:,:) = ay(:,:)-ey(:,:)*dt

  !! check if negative pressure / density
  ! re-calculate total energy (after flooring pressure and density)
  do jj = 1, liz
  do ii = 1, lix
    if(pr(ii,jj) <= pr_floor) then
      pr(ii,jj) = pr_floor
    end if
    if(ro(ii,jj) <= ro_floor) then
      ro(ii,jj) = ro_floor
    end if
  end do
  end do
  call calc_etot(ro,pr,vx,vy,vz,bx,by,bz,etot,gm)

  ns= ns+1
  t = t+dt
  !--------------------------------------------------
  ! time integration is finished.
  !--------------------------------------------------
  !--------------------------------------------------
  ! output
  !--------------------------------------------------
  if(t > tout) then
  !if(mod(ns,10)==0) then
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
  end if

  ! terminal output
  if(mod(ns,100)==0) then
  !if(mod(ns,1)==0) then
    if(my_rank==0) then
      write(*,*) 'step:', ns, t
    end if
  end if
  
  !--------------------------------------------------
  ! for debug
  !--------------------------------------------------
  !if(nd>=10) then 
    !call MPI_FINALIZE(ierr)
    !stop
  !end if
  !--------------------------------------------------
  !--------------------------------------------------
  ! for debug
  !--------------------------------------------------
  ! call MPI_FINALIZE(ierr)
  ! stop
  !--------------------------------------------------
  ! debug part end
  !--------------------------------------------------

  end do main_loop

  if(my_rank==0) then
    write(*,fstpmsg) ns,t
    write(*,*) 'Normal Termination'
  end if


  !--------------------------------------------------
  ! ending message
  !--------------------------------------------------
  call MPI_FINALIZE(ierr)



end program main
