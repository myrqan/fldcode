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
  double precision, parameter :: pi = 4 * atan(1.d0)

  character(len=100) :: formatoutdata,formatoutmessage,formatstopmessage

  integer,parameter :: margin=1
  integer,parameter :: gr_num_x=1000
  integer,parameter :: ix=2*margin+gr_num_x

  !double precision :: x(0:ix),dx
  double precision :: x(ix),dx
  !double precision :: s(0:ix), sm(0:ix)
  double precision :: s(ix),sm(ix)
  !double precision :: u(0:ix),um(0:ix),un(0:ix),r(0:ix)
  double precision :: u(ix),um(ix),un(ix)
  double precision :: r(ix)
  !double precision :: f(0:ix),f0(0:ix)
  double precision :: f(ix),f0(ix)

  !double precision :: rho(0:ix),rhom(0:ix),rhon(0:ix)
  double precision :: rho(ix),rhom(ix),rhon(ix)
  !double precision :: vx(0:ix),vxm(0:ix),vxn(0:ix)
  double precision :: vx(ix),vxm(ix),vxn(ix)
  !double precision :: p(0:ix),pm(0:ix),pn(0:ix)
  double precision :: p(ix),pm(ix),pn(ix)
  !double precision :: eps(0:ix),epsm(0:ix),epsn(0:ix)
  double precision :: eps(ix),epsm(ix),epsn(ix)
  !double precision :: kappa(0:ix),qv

  double precision :: kappa(ix),qv

  double precision, parameter :: gamma = 5.d0/3.d0

  double precision :: t,dt,tend
  integer :: ns,nd
  double precision :: tout,dtout

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
  call init1d(s); call init1d(sm)
  call init1d(u); call init1d(um); call init1d(un)
  call init1d(r)
  call init1d(f); call init1d(f0)
  call init1d(rho); call init1d(rhom); call init1d(rhon)
  call init1d(vx); call init1d(vxm); call init1d(vxn)
  call init1d(p); call init1d(pm); call init1d(pn)
  call init1d(eps); call init1d(epsm); call init1d(epsn)
  call init1d(kappa)

  !----------------------------------------------------------------------|
  !   time control parameters
  tend=5.d0  ! time for end of calculation
  dtout=0.5d0 ! time spacing for data output
  !----------------------------------------------------------------------|
  !  initialize counters
  t=0.d0
  ns=0      ! number of time steps
  tout=0.d0  ! time for next output
  nd=1        ! number of output data sets
  !----------------------------------------------------------------------|

  ! initialize output files
  call reset(ix)

  !   setup numerical model (grid, initial conditions, etc.)
  call gridx(ix,margin,dx,x)
  call init_shocktube(ix,x,rho,p)
  !call init_sedov1d(ix,x,rho,p)

  do i = 1, ix
  s(i) = x(i)**2
  if (i < ix) then
    sm(i) = (0.5d0*(x(i)+x(i+1)))**2
  endif
  enddo

  !----------------------------------------------------------------------|

  call put1dreal(11,'x.dac',x)
  call put0dreal(10,'t.dac',t)
  call put1dreal(15,'rho.dac',rho)
  call put1dreal(16,'vx.dac',vx)
  call put1dreal(17,'p.dac',p)
  call put1dreal(18,'eps.dac',eps)


  call put1dreal_txt(21,'x.txt',x)
  call put0dreal_txt(20,'t.txt',t)
  call put1dreal_txt(25,'rho.txt',rho)
  call put1dreal_txt(26,'vx.txt',vx)
  call put1dreal_txt(27,'p.txt',p)
  call put1dreal_txt(28,'eps.txt',eps)

    !======================================================================|
    !     time integration 
    !======================================================================|
    do while (t < tend)
    ns=ns+1
    !----------------------------------------------------------------------|
    !     time spacing

    dt = 3.d-4
    dtout = dt
    tend = dt * 10.d0
    !----------------------------------------------------------------------|
    !     solve difference equations
    !     by two-step Lax-Wendroff scheme

    !----------------------------------------------------------------------|
    ! FIRST STEP
    !----------------------------------------------------------------------|

    eps = p / (gamma - 1.d0) + 0.5d0 * rho * vx**2

    f0 = rho * vx * s
    u = rho * s
    call mlw1d1st(u,um,f0,ix,dt,dx)
    rhom = um / sm

    f0 = (rho * vx**2 + p) * s
    u = rho * vx * s
    call mlw1d1st(u,um,f0,ix,dt,dx)
    vxm = um / rhom / sm

    f0 = (eps + p) * vx * s
    u = eps * s
    call mlw1d1st(u,um,f0,ix,dt,dx)
    epsm = um / sm

    !----------------------------------------------------------------------|
    ! SECOND STEP
    !----------------------------------------------------------------------|

    pm = (gamma-1.d0) * (epsm - 0.5d0 * rhom * vxm**2)

    f = rhom*vxm*sm
    u = rho*s
    call mlw1d2nd(u,un,f,ix,dt,dx)
    rhon = un/s

    f = (rhom*vxm**2 + pm) * sm
    u = rho*vx * s
    do i = 1, ix
    r(i) = p(i) * 2*x(i)
    enddo

    call mlw1d2nd(u,un,f,ix,dt,dx)
    call mlw1dsrc(un,r,ix,dt,dx)
    vxn = un/rhon/s

    f = (epsm + pm) * vxm * sm
    u = eps * s
    call mlw1d2nd(u,un,f,ix,dt,dx)
    epsn = un / s

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
    !call bnd1d_f_lr(vx,ix)
    call bnd1d_fix_l(vx, ix)
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

    u=rho*s
    call arvis1d(kappa,u,um,ix,dt,dx)
    rhom=um/s

    u = rho*vx*s
    call arvis1d(kappa,u,um,ix,dt,dx)
    vxm = um/rhom/s

    u = eps*s
    call arvis1d(kappa,u,um,ix,dt,dx)
    epsm = um/s

    pm = (gamma - 1.d0) * (epsm - 0.5d0 * rhom * vxm**2)
    ! update
    rho(:)=rhom(:)
    vx(:)=vxm(:)
    eps(:) = epsm(:)
    p(:) = pm(:)

    !----------------------------------------------------------------------|
    !     boundary condition
    call bnd1d_f_lr(rho,ix)
    !call bnd1d_f_lr(vx,ix)
    call bnd1d_fix_l(vx,ix)
    call bnd1d_f_lr(eps,ix)
    call bnd1d_f_lr(p,ix)

    !----------------------------------------------------------------------|

    !     data output 
    if (t >= tout) then

      call put0dreal(10,'t.dac',t)
    call put1dreal(15,'rho.dac',rho)
    call put1dreal(16,'vx.dac',vx)
    call put1dreal(17,'p.dac',p)
    call put1dreal(18,'eps.dac',eps)

    call put0dreal_txt(20,'t.txt',t)
    call put1dreal_txt(25,'rho.txt',rho)
    call put1dreal_txt(26,'vx.txt',vx)
    call put1dreal_txt(27,'p.txt',p)
    call put1dreal_txt(28,'eps.txt',eps)
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
