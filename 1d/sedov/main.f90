program main
  !======================================================================|
  !     use MODULES
  !======================================================================|
  use file_output
  use initial_condition
  use mlw
  use boundary
  use arvis
  !======================================================================|
  !     array definitions
  !======================================================================|
  IMPLICIT NONE
  character(len=100) :: formatoutdata,formatoutmessage,formatstopmessage
  integer :: mfile
  integer,parameter :: margin = 1
  integer,parameter :: grid_x = 1024
  integer,parameter :: ix=2*margin+grid_x
  double precision :: x(ix),dx
  double precision :: u(ix),um(ix),un(ix)
  double precision :: r(ix),rm(ix),rn(ix)
  double precision :: f(ix),f0(ix)
  double precision :: t,dt,tend
  integer :: ns, nd
  double precision :: tout,dtout
  double precision, parameter :: pi = 4 * atan(1.d0)
  double precision :: rho(ix),rhom(ix),rhon(ix)
  double precision :: vx(ix),vxm(ix),vxn(ix)
  double precision :: p(ix),pm(ix),pn(ix)
  double precision :: eps(ix),epsm(ix),epsn(ix)
  double precision :: s(ix), sm(ix)
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

  call init1d(x)
  call init1d(u); call init1d(um); call init1d(un)
  call init1d(r); call init1d(rm); call init1d(rn)
  call init1d(f); call init1d(f0)
  call init1d(rho); call init1d(rhom); call init1d(rhon)
  call init1d(vx); call init1d(vxm); call init1d(rhon)
  call init1d(p); call init1d(pm); call init1d(pn)
  call init1d(eps); call init1d(epsm); call init1d(epsn)
  call init1d(s); call init1d(sm)
  call init1d(kappa)

  !----------------------------------------------------------------------|
  !   time control parameters
  tend=0.6d1  ! time for end of calculation
  dtout=0.1d1 ! time spacing for data output
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
  call gridx(ix,margin,dx,x)
  call init_shocktube(ix,x,rho,p)
  call init_sedov1d(ix,x,rho,p)

  do i = 1, ix
  s(i) = x(i)**2
  if (i < ix) then
    sm(i) = (0.5d0*(x(i) + x(i+1)))**2
  endif
  enddo

  !----------------------------------------------------------------------|

  call put1dreal(11,'x.dac',x)
  call put0dreal(10,'t.dac',t)
  call put1dreal(15,'rho.dac',rho)
  call put1dreal(16,'vx.dac',vx)
  call put1dreal(17,'p.dac',p)
  call put1dreal(18,'eps.dac',eps)

  !======================================================================|
  !     time integration 
  !======================================================================|
  do while (t < tend)
  ns=ns+1
  !----------------------------------------------------------------------|
  !     time spacing

  dt = 2.d-4
  !----------------------------------------------------------------------|
  !     solve difference equations
  !     by two-step Lax-Wendroff scheme

  !----------------------------------------------------------------------|
  ! FIRST STEP
  !----------------------------------------------------------------------|

  eps = p / (gamma - 1.d0) + 0.5d0 * rho * vx**2

  f0 = rho*vx*s
  u = rho*s
  call mlw1d1st(u,um,f0,ix,dt,dx)
  rhom = um / sm

  f0 = (eps + p) * vx * s
  u = eps * s
  call mlw1d1st(u,um,f0,ix,dt,dx)
  epsm = um / sm

  f0 = (rho * vx**2 + p) * s
  u = rho * vx * s
  r = p * 2.d0 * x
  call mlw1d1st(u,um,f0,ix,dt,dx)

  ! CORRECTION
  do i = 1, ix
  r(i) = p(i) * 2.d0 * x(i)
  rm(i) = pm(i) * 2.d0 * x(i)
  enddo
  call mlw1dsrc1st(um,r,rm,ix,dt,dx)
  vxm=um/rhom/sm

  pm = (gamma-1.d0) * (epsm - 0.5d0 * rhom * vxm**2)



  !----------------------------------------------------------------------|
  ! SECOND STEP
  !----------------------------------------------------------------------|

  !pm = (gamma-1.d0) * (epsm - 0.5d0 * rhom * vxm**2)

  ! MASS CONSERVATION
  f = rhom*vxm*sm
  u = rho*s
  call mlw1d2nd(u,un,f,ix,dt,dx)
  rhon = un / s

  ! ENERGY CONSERVATION 
  f = (epsm + pm) * vxm *sm
  u = eps*s
  call mlw1d2nd(u,un,f,ix,dt,dx)
  epsn = un/s

  ! MOMENTUM CONSERVATION
  f = (rhom*vxm**2+pm)*sm
  u = rho*vx*s
  call mlw1d2nd(u,un,f,ix,dt,dx)
  vxn=un/rhon/s
  pn = (gamma - 1.d0) * (epsn - 0.5d0 * rhon * vxn**2)

  ! CORRECTION
  do i = 1, ix
  ! SOURCE TERM
  r(i) = p(i) * 2.d0 * x(i)
  rn(i) = pn(i) * 2.d0 * x(i)
  enddo
  call mlw1dsrc2nd(un,r,rn,ix,dt,dx)
  vxn = un/rhon/s

  pn = (gamma - 1.d0) * (epsn - 0.5d0 * rhon * vxn**2)

  !----------------------------------------------------------------------|
  !     update
  rho(:) = rhon(:)
  vx(:) = vxn(:)
  eps(:) = epsn(:)
  p(:) = pn(:)

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

  qv=3.0d0

  if(nd == 1) then
    call qv_param(qv)
  endif

  call calc_coef_k(qv,kappa,vx,ix,dx)

  u=rho*s
  call arvis1d(kappa,u,um,ix,dt,dx)
  rhom=um/s

  u = rho*vx*s
  call arvis1d(kappa,u,um,ix,dt,dx)
  vxm = um/rhom/s

  u = eps*s
  call arvis1d(kappa,u,um,ix,dt,dx)
  epsm = um/s

  ! update
  rho(:)=rhom(:)
  vx(:)=vxm(:)
  eps(:) = epsm(:)

  !----------------------------------------------------------------------|
  !     boundary condition
  call bnd1d_f_lr(rho,ix)
  call bnd1d_fix_l(vx,ix)
  call bnd1d_f_r(vx,ix)
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
