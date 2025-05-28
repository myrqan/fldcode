program main
  !======================================================================|
  !     use MODULES
  !======================================================================|
  use file_output
  use initial_condition
  use mlw
  use boundary2d
  use arvis
  use cfl
  !======================================================================|
  !     array definitions
  !======================================================================|
  IMPLICIT NONE
  character(len=100) :: formatoutdata,formatoutmessage,formatstopmessage
  integer :: mfile
  integer,parameter :: margin = 1
  integer,parameter :: grid_x = 100
  integer,parameter :: grid_z = 100
  integer,parameter :: ix=2*margin+grid_x
  integer,parameter :: jx=2*margin+grid_z
  double precision :: x(ix,jx),xm(ix,jx),dx
  double precision :: z(ix,jx),zm(ix,jx),dz
  double precision :: u(ix,jx),um(ix,jx),un(ix,jx)
  double precision :: r(ix,jx),rm(ix,jx),rn(ix,jx)
  double precision :: f(ix,jx),f0(ix,jx),fm(ix,jx)
  double precision :: t,dt,tend
  integer :: ns, nd
  double precision :: tout,dtout
  double precision, parameter :: pi = 4 * atan(1.d0)
  double precision :: rho(ix,jx),rhom(ix,jx),rhon(ix,jx)
  double precision :: vx(ix,jx),vxm(ix,jx),vxn(ix,jx)
  double precision :: vz(ix,jx),vzm(ix,jx),vzn(ix,jx)
  double precision :: p(ix,jx),pm(ix,jx),pn(ix,jx)
  double precision :: eps(ix,jx),epsm(ix,jx),epsn(ix,jx)
  !double precision :: s(ix,jx), sm(ix,jx), ds(ix,jx), dsm(ix,jx)
  !double precision :: r(ix,jx),rm(ix,jx),dr(ix,jx)
  double precision :: dro(ix,jx),drvx(ix,jx),drvz(ix,jx),de(ix,jx)
  !double precision :: rvx(ix,jx),rvxm(ix,jx),drvx(ix,jx)
  !double precision :: rvz(ix,jx),rvzm(ix,jx),drvz(ix,jx)
  !double precision :: e(ix,jx),em(ix,jx),des(ix,jx)
  double precision :: kappa(ix,jx),qv
  double precision, parameter :: gamma = 5.d0/3.d0

  integer :: i,j
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

  call init2d(x); call init2d(xm)
  call init2d(u); call init2d(um); call init2d(un)
  call init2d(r); call init2d(rm); call init2d(rn)
  call init2d(f); call init2d(f0)
  call init2d(rho); call init2d(rhom); call init2d(rhon)
  call init2d(vx); call init2d(vxm); call init2d(vxn)
  call init2d(vz); call init2d(vzm); call init2d(vzn)
  call init2d(p); call init2d(pm); call init2d(pn)
  call init2d(eps); call init2d(epsm); call init2d(epsn)
  
  call init2d(dro); call init2d(drvx); call init2d(drvz); call init2d(de)
  !call init2d(s); call init2d(sm); call init2d(ds); call init2d(dsm)
  !call init2d(rs); call init2d(rsm); call init2d(drs)
  !call init2d(rvxs); call init2d(rvxsm); call init2d(drvxs)
  !call init2d(es); call init2d(esm); call init2d(des)
  call init2d(kappa)

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
  call reset(ix,jx)

  ! setup a numerical model (grid, initial conditions, etc.)
  call grid2d(ix,jx,margin,dx,dz,x,z,xm,zm)

  call sedov2d(ix,jx,x,z,rho,vx,vz,p)

  !! apply boundary conditions
  call bc_free_b(rho,ix,jx); call bc_free_u(rho,ix,jx)
  call bc_free_l(rho,ix,jx); call bc_free_r(rho,ix,jx)

  call bc_free_b(vx,ix,jx); call bc_free_u(vx,ix,jx)
  call bc_fix_l(vx,ix,jx); call bc_free_r(vx,ix,jx)
  
  call bc_fix_b(vz,ix,jx); call bc_free_u(vz,ix,jx)
  call bc_free_l(vz,ix,jx); call bc_free_r(vz,ix,jx)

  call bc_free_b(p,ix,jx); call bc_free_u(p,ix,jx)
  call bc_free_l(p,ix,jx); call bc_free_r(p,ix,jx)


  !----------------------------------------------------------------------|
  !! BINARY OUTPUT
  ! calc eps
  eps(:,:) = p(:,:)/(gamma-1.d0) + 0.5d0*rho(:,:)*(vx(:,:)**2+vz(:,:)**2)

  call put1dreal(11,'x.dac',x)
  call put0dreal(10,'t.dac',t)
  call put2dreal(16,'vx.dac',vx)
  call put2dreal(17,'vz.dac',vx)
  call put2dreal(20,'rho.dac',rho)
  call put2dreal(21,'p.dac',p)
  call put2dreal(22,'eps.dac',eps)
  write(*,formatoutmessage) ns,t,nd
  tout = tout + dtout
  nd = nd + 1



  !!!!!!!!!!!!!!!!!!!!ここまでやった
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
    call put2dreal(16,'vx.dac',vx)
    call put2dreal(17,'vz.dac',vz)
    call put2dreal(20,'rho.dac',rho)
    call put2dreal(21,'p.dac',p)
    call put2dreal(22,'eps.dac',eps)
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
