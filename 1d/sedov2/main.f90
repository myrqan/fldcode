program main
  !======================================================================|
  !     use MODULES
  !======================================================================|
  use file_output
  use initial_condition
  use mlw
  use boundary
  use arvis
  use cfl
  !======================================================================|
  !     array definitions
  !======================================================================|
  IMPLICIT NONE
  character(len=100) :: formatoutdata,formatoutmessage,formatstopmessage
  integer :: mfile
  integer,parameter :: margin = 1
  integer,parameter :: grid_x = 1000
  integer,parameter :: ix=2*margin+grid_x
  double precision :: x(ix),xm(ix),dx
  double precision :: u(ix),um(ix),un(ix)
  double precision :: r(ix),rm(ix),rn(ix)
  double precision :: f(ix),f0(ix),fm(ix)
  double precision :: t,dt,tend
  integer :: ns, nd
  double precision :: tout,dtout
  double precision, parameter :: pi = 4 * atan(1.d0)
  double precision :: rho(ix),rhom(ix),rhon(ix)
  double precision :: vx(ix),vxm(ix),vxn(ix)
  double precision :: p(ix),pm(ix),pn(ix)
  double precision :: eps(ix),epsm(ix),epsn(ix)
  double precision :: s(ix), sm(ix), ds(ix), dsm(ix)
  double precision :: rs(ix),rsm(ix),drs(ix)
  double precision :: rvxs(ix),rvxsm(ix),drvxs(ix)
  double precision :: es(ix),esm(ix),des(ix)
  double precision :: kappa(ix),qv
  double precision, parameter :: gamma = 5.d0/3.d0

  integer :: i
  !======================================================================|
  !     prologue
  !======================================================================|
  !----------------------------------------------------------------------|
  !   output format
  formatoutdata= '(1x,e13.5e3,1x,e13.5e3,1x,e13.5e3,1x,e13.5e3,1x,e13.5e3)'
  formatoutmessage= '(1x," write    ","step=",i8," t=",e10.3," nd =",i3)'
  formatstopmessage= '(1x," stop     ","step=",i8," t=",e10.3)'
  !----------------------------------------------------------------------|
  !   initialize

  call init1d(x); call init1d(xm)
  call init1d(u); call init1d(um); call init1d(un)
  call init1d(r); call init1d(rm); call init1d(rn)
  call init1d(f); call init1d(f0)
  call init1d(rho); call init1d(rhom); call init1d(rhon)
  call init1d(vx); call init1d(vxm); call init1d(rhon)
  call init1d(p); call init1d(pm); call init1d(pn)
  call init1d(eps); call init1d(epsm); call init1d(epsn)
  call init1d(s); call init1d(sm); call init1d(ds); call init1d(dsm)
  call init1d(rs); call init1d(rsm); call init1d(drs)
  call init1d(rvxs); call init1d(rvxsm); call init1d(drvxs)
  call init1d(es); call init1d(esm); call init1d(des)
  call init1d(kappa)

  !----------------------------------------------------------------------|
  !   time control parameters
  tend=6.d0  ! time for end of calculation
  dtout=0.5d0! time spacing for data output
  !----------------------------------------------------------------------|
  !  initialize counters
  t=0.d0
  ns=0      ! number of tme steps
  tout=0.d0  ! time for next output
  nd=1        ! number of output data sets
  !----------------------------------------------------------------------|

  ! initialize output files
  call reset(ix)

  !   setup numerical model (grid, initial conditions, etc.)
  call gridx(ix,margin,dx,x,xm)
  !call init_shocktube(ix,x,rho,p)
  call init_sedov1d(ix,x,rho,vx,p)

  call cross_section(ix,s,sm,ds,dsm,x,xm,dx)

 ! do i = 1, ix
 ! s(i) = x(i)**2
 ! if (i < ix) then
 !   sm(i) = xm(i)**2
 ! endif
 ! enddo


  !! boundary condition
  call bnd1d_f_lr(rho,ix)
  call bnd1d_fix_l(vx,ix)
  call bnd1d_f_r(vx,ix)
  call bnd1d_f_lr(p,ix)

  call bnd1d_f_lr(s,ix)
  call bnd1d_fix_l(ds,ix)
  call bnd1d_f_r(ds,ix)


  !----------------------------------------------------------------------|
  !! BINARY OUTPUT
  ! calc eps
  eps(:) = p(:)/(gamma-1.d0) + 0.5d0*rho(:)*vx(:)**2

  call put1dreal(11,'x.dac',x)
  call put0dreal(10,'t.dac',t)
  call put1dreal(15,'rho.dac',rho)
  call put1dreal(16,'vx.dac',vx)
  call put1dreal(17,'p.dac',p)
  call put1dreal(18,'eps.dac',eps)
  write(*,formatoutmessage) ns,t,nd
  tout = tout + dtout
  nd = nd + 1

  !======================================================================|
  !     time integration 
  !======================================================================|
  do while (t < tend)
  ns=ns+1
  !----------------------------------------------------------------------|
  !     time spacing

  !dt = 0.5d-4
  CALL calc_dt(rho,vx,p,gamma,ix,dt,dx)
  !! temporary
!  dtout =  10 * dt
!  tend = 10 * dtout
  !---------------------------------------------------------------------|
  ! FIRST STEP
  !----------------------------------------------------------------------|
  !eps(:) = p(:)/(gamma-1.d0) + 0.5d0*rho(:)*vx(:)**2
  ! initialize d** 
  CALL init1d(drs); CALL init1d(drvxs); CALL init1d(des)


  rs(:) = rho(:)*s(:)
  f0(:) = rho(:)*vx(:)*s(:)
  !r(:) = 0.d0
  CALL mlw1d1st(rs,f0,rsm,drs,ix,dt,dx)

  rvxs(:) = rho(:)*vx(:)*s(:)
  f0(:) = (rho(:)*vx(:)**2+p(:))*s(:)
  r(:) = p(:)*ds(:)
  CALL mlw1d1st(rvxs,f0,rvxsm,drvxs,ix,dt,dx)
  CALL mlw1dsrc1st(rvxsm,drvxs,r,ix,dt)

  es(:) = eps(:)*s(:)
  f0(:) = (eps(:)+p(:))*vx(:)*s(:)
  !r(:) = 0.d0
  call mlw1d1st(es,f0,esm,des,ix,dt,dx)

  rhom(:) = rsm(:)/sm(:)
  vxm(:) = rvxsm(:)/rhom(:)/sm(:)
  epsm(:) = esm(:)/sm(:)

  pm(:) = (gamma-1.d0) * (epsm(:)-0.5d0*rhom(:)*vxm(:)**2)


  !----------------------------------------------------------------------|
  ! SECOND STEP
  !----------------------------------------------------------------------|

  fm(:) = rhom(:)*vxm(:)*sm(:)
  CALL mlw1d2nd(fm,drs,ix,dt,dx)

  fm(:) = (rhom(:)*vxm(:)**2+pm(:))*sm(:)
  rm(:) = pm(:)*ds(:)
  CALL mlw1d2nd(fm,drvxs,ix,dt,dx)
  CALL mlw1dsrc2nd(drvxs,rm,ix,dt)

  fm(:) = (epsm(:)+pm(:))*vxm(:)*sm(:)
  CALL mlw1d2nd(fm,des,ix,dt,dx)

  rs(:) = rs(:)+drs(:)
  rvxs(:) = rvxs(:) + drvxs(:)
  es(:) = es(:) + des(:)

  rhon(:) = rs(:)/s(:)
  vxn(:) = rvxs(:)/rhon(:)/s(:)
  epsn(:) = es(:)/s(:)

  pn(:) = (gamma - 1.d0) * (epsn(:) - 0.5d0*rhon(:)*vxn(:)**2)

  !----------------------------------------------------------------------|
  !     update
  rho(:) = rhon(:); vx(:) = vxn(:); eps(:) = epsn(:); p(:) = pn(:)

  !----------------------------------------------------------------------|
  !     boundary condition
  call bnd1d_f_lr(rho,ix)
  !call bnd1d_f_lr(vx,ix)
  call bnd1d_fix_l(vx,ix)
  call bnd1d_f_r(vx,ix)
  call bnd1d_f_lr(eps,ix)
  call bnd1d_f_lr(p,ix)

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! artificial viscosity
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  qv=3.d0

  if(nd == 1) then
    call qv_param(qv)
  endif

  call calc_coef_k(qv,kappa,vx,ix,dx)

  u=rho(:)*s(:)
  call arvis1d(kappa,u,um,ix,dt,dx)
  rhom(:)=um(:)/s(:)

  u(:) = rho(:)*vx(:)*s(:)
  call arvis1d(kappa,u,um,ix,dt,dx)
  vxm(:) = um(:)/rhom(:)/s(:)

  u(:) = eps(:)*s(:)
  call arvis1d(kappa,u,um,ix,dt,dx)
  epsm(:) = um(:)/s(:)

  ! update
  rho(:)=rhom(:)
  vx(:)=vxm(:)
  eps(:) = epsm(:)

  !----------------------------------------------------------------------|
  !     boundary condition
  call bnd1d_f_lr(rho,ix)
  call bnd1d_fix_l(vx,ix)
  call bnd1d_f_r(vx,ix)
  !call bnd1d_f_lr(vx,ix)
  call bnd1d_f_lr(eps,ix)

  !----------------------------------------------------------------------|

  t=t+dt
  !     data output 
  if (t >= tout) then


    call put0dreal(10,'t.dac',t)
    call put1dreal(15,'rho.dac',rho)
    call put1dreal(16,'vx.dac',vx)
    call put1dreal(17,'p.dac',p)
    call put1dreal(18,'eps.dac',eps)
    write(*,formatoutmessage) ns,t,nd

    tout=tout+dtout
    nd=nd+1
  endif

  enddo ! do while (t < tend)

  !======================================================================|
  !     ending msg
  !======================================================================|
  write(*,formatstopmessage) ns,t
  write(*,*) '  ### normal stop ###'

  stop
end program main
