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
  character(len=80) :: filename
  character(len=30) :: filenum
  integer :: mfile
  integer,parameter :: ix=1000
  double precision :: x(0:ix),dx
  double precision :: u(0:ix),um(0:ix),un(0:ix)
  double precision :: f(0:ix),f0(0:ix)
  double precision :: t,dt,tend
  integer :: ns, nd
  double precision :: tout,dtout
  double precision, parameter :: pi = 4 * atan(1.d0)
  double precision :: rho(0:ix),rhom(0:ix),rhon(0:ix)
  double precision :: vx(0:ix),vxm(0:ix),vxn(0:ix)
  double precision :: p(0:ix),pm(0:ix),pn(0:ix)
  double precision :: eps(0:ix),epsm(0:ix),epsn(0:ix)
  double precision :: kappa(0:ix),qv
  double precision, parameter :: gamma = 1.4d0
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
  call init1d(f); call init1d(f0)
  call init1d(rho); call init1d(rhom); call init1d(rhon)
  call init1d(vx); call init1d(vxm); call init1d(rhon)
  call init1d(p); call init1d(pm); call init1d(pn)
  call init1d(eps); call init1d(epsm); call init1d(epsn)
  call init1d(kappa)

  !----------------------------------------------------------------------|
  !   time control parameters
  tend=0.141d0  ! time for end of calculation
  dtout=0.02d0 ! time spacing for data output
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
  call gridx(ix,dx,x)
  call init_shocktube(ix,x,rho,p)

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
  dt=1.d-4
  !----------------------------------------------------------------------|
  !     solve difference equations
  !     by two-step Lax-Wendroff scheme

  !----------------------------------------------------------------------|
  ! FIRST STEP
  !----------------------------------------------------------------------|

  eps = p / (gamma - 1.d0) + 0.5d0 * rho * vx**2

  f0 = rho*vx
  u = rho
  call mlw1d1st(u,um,f0,ix,dt,dx)
  rhom = um

  f0 = rho * vx**2 + p
  u = rho * vx
  call mlw1d1st(u,um,f0,ix,dt,dx)
  vxm=um/rhom

  f0 = (eps + p) * vx
  u = eps
  call mlw1d1st(u,um,f0,ix,dt,dx)
  epsm = um

  !----------------------------------------------------------------------|
  ! SECOND STEP
  !----------------------------------------------------------------------|

  pm = (gamma-1.d0) * (epsm - 0.5d0 * rhom * vxm**2)

  f = rhom*vxm
  u = rho
  call mlw1d2nd(u,un,f,ix,dt,dx)
  rhon = un

  f = rhom*vxm**2 + pm
  u = rho*vx
  call mlw1d2nd(u,un,f,ix,dt,dx)
  vxn = un/rhon

  f = (epsm + pm) * vxm
  u = eps
  call mlw1d2nd(u,un,f,ix,dt,dx)
  epsn = un

  pn = (gamma - 1.d0) * (epsn - 0.5d0 * rhon * vxn**2)

  !----------------------------------------------------------------------|
  !     update
  rho(:) = rhon(:)
  vx(:) = vxn(:)
  eps(:) = epsn(:)
  p(:) = pn(:)

  t=t+dt
  !----------------------------------------------------------------------|
  !     boundary condition
  call bnd1d_f_lr(rho,ix)
  call bnd1d_f_lr(vx,ix)
  call bnd1d_f_lr(eps,ix)
  call bnd1d_f_lr(p,ix)

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! artificial viscosity
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  qv=4.0d0

  if(nd == 1) then
    call qv_param(qv)
  endif

  call calc_coef_k(qv,kappa,vx,ix,dx)

  u=rho
  call arvis1d(qv,kappa,u,um,ix,dt,dx)
  rhom=um

  u = rho*vx
  call arvis1d(qv,kappa,u,um,ix,dt,dx)
  vxm = um/rhom

  u = eps
  call arvis1d(qv,kappa,u,um,ix,dt,dx)
  epsm = um

  ! update
  rho(:)=rhom(:)
  vx(:)=vxm(:)
  eps(:) = epsm(:)

  !----------------------------------------------------------------------|
  !     boundary condition
  call bnd1d_f_lr(rho,ix)
  call bnd1d_f_lr(vx,ix)
  call bnd1d_f_lr(eps,ix)

  !----------------------------------------------------------------------|

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
