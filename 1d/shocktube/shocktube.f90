program main
  use file_output
  use initial_condition
  !======================================================================|
  !     array definitions
  !======================================================================|
  implicit none
  character(len=100) :: formatoutdata,formatoutmessage,formatstopmessage
  character(len=80) :: filename
  character(len=30) :: filenum
  integer :: mfile
  integer,parameter :: ix=1000
  double precision :: x(0:ix),dx
  double precision :: u(0:ix),um(0:ix)
  double precision :: f(0:ix),f0(0:ix)
  double precision :: t,dt,tend
  integer :: ns
  double precision :: tout,dtout
  integer :: nd
  integer :: i
  !double precision :: cs
  double precision, parameter :: pi = 4 * atan(1.d0)
  double precision :: rho(0:ix),rhom(0:ix),rhon(0:ix)
  double precision :: vx(0:ix),vxm(0:ix),vxn(0:ix)
  double precision :: p(0:ix),pm(0:ix),pn(0:ix)
  double precision :: eps(0:ix),epsm(0:ix),epsn(0:ix)
  double precision :: kappa(0:ix),qv
  double precision, parameter :: gamma = 1.4d0
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

  x(:)=0.d0
  u(:)=0.d0; um(:)=0.d0
  f(:)=0.d0; f0(:)=0.d0
  rho(:)=0.d0; rhom(:)=0.d0; rhon(:)=0.d0
  vx(:)=0.d0; vxm(:)=0.d0; vxn(:)=0.d0
  p(:) = 0.d0; pm(:) = 0.d0; pn(:) = 0.d0
  eps(:)=0.d0; epsm(:) = 0.d0; epsn(:) = 0.d0
  kappa(:)=0.d0

  !----------------------------------------------------------------------|
  !   time control parameters
  tend=0.141d0  ! time for end of calculation
  dtout=0.02d0 ! time spacing for data output
  !----------------------------------------------------------------------|
  !  file open
  !mfile=10; filename='out.dat'
  !open(mfile,file=filename,form='formatted')
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

  call gridx(ix,x)
  !dx=1.d0/dble(ix)
  !x(0)=0.d0
  !do i=1,ix
  !x(i)=x(i-1)+dx
  !enddo

  !vx(:)=0.d0

  call init_shocktube(ix,x,rho,p)
  !do i=0,ix
  !if (x(i) <= 0.5d0) then
    !rho(i)=1.d0
  !else
    !rho(i)=0.125d0
  !endif
  !enddo
  !do i = 0, ix
  !if(x(i) <= 0.5d0) then
    !p(i) = 1.d0
  !else
    !p(i) = 0.1d0
  !endif
  !enddo
  !----------------------------------------------------------------------|

  call put1dreal(11,'x.dac',x)

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
  ! FIRST STEP
  eps(:) = p(:) / (gamma - 1.d0) + 0.5d0 * rho(:) * vx(:)**2

  f0(:)=rho(:)*vx(:)
  u(:)=rho(:)
  do i=0,ix-1
  um(i)=0.5d0*(u(i+1)+u(i))-0.5d0*dt/dx*(f0(i+1)-f0(i))
  enddo
  rhom(:)=um(:)

  f0(:)=rho(:)*vx(:)**2 + p(:)
  u(:)=rho(:)*vx(:)
  do i=0,ix-1
  um(i)=0.5d0*(u(i+1)+u(i))-0.5d0*dt/dx*(f0(i+1)-f0(i))
  enddo
  vxm(:)=um(:)/rhom(:)

  f0(:) = (eps(:) + p(:)) * vx(:)
  u(:) = eps(:)
  do i = 0, ix-1
  um(i)=0.5d0*(u(i+1)+u(i))-0.5d0*dt/dx*(f0(i+1)-f0(i))
  enddo
  epsm(:) = um(:)

  ! SECOND STEP

  pm(:) = (gamma - 1.d0) * (epsm(:) - 0.5d0 * rhom(:) * vxm(:)**2)

  f(:)=rhom(:)*vxm(:)
  u(:)=rho(:)
  do i=1,ix-1
  u(i)=u(i)-dt/dx*(f(i)-f(i-1))
  enddo
  rhon(:)=u(:)

  f(:)=rhom(:)*vxm(:)**2+pm(:)
  u(:)=rho(:)*vx(:)
  do i=1,ix-1
  u(i)=u(i)-dt/dx*(f(i)-f(i-1))
  enddo
  vxn(:)=u(:)/rhon(:)

  f(:) = (epsm(:) + pm(:)) * vxm(:)
  u(:) = eps(:)
  do i=1,ix-1
  u(i)=u(i)-dt/dx*(f(i)-f(i-1))
  enddo
  epsn(:) = u(:)

  pn(:) = (gamma - 1.d0) * (epsn(:) - 0.5d0 * rhon(:) * vxn(:)**2)
  !----------------------------------------------------------------------|
  !     update
  rho(:)=rhon(:)
  vx(:)=vxn(:)
  eps(:) = epsn(:)
  p(:) = pn(:)
  t=t+dt
  !----------------------------------------------------------------------|
  !     boundary condition
  rho(0)=rho(1)
  vx(0)=vx(1)
  eps(0) = eps(1)
  p(0) = p(1)
  rho(ix)=rho(ix-1)
  vx(ix)=vx(ix-1)
  eps(ix) = eps(ix-1)
  p(ix) = p(ix-1)
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! artificial viscosity
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  qv=5.0d0
  do i=0,ix-1
  kappa(i)=qv*dx*abs(vx(i+1)-vx(i))
  enddo

  u(:)=rho(:)
  do i=1,ix-1
  um(i)=u(i)+dt/dx*( kappa(i  )/dx*(u(i+1)-u(i)) &
    &                - kappa(i-1)/dx*(u(i)-u(i-1)) )
  enddo
  rhom(:)=um(:)

  u(:)=rho(:)*vx(:)
  do i=1,ix-1
  um(i)=u(i)+dt/dx*( kappa(i  )/dx*(u(i+1)-u(i)) &
    &                - kappa(i-1)/dx*(u(i)-u(i-1)) )
  enddo
  vxm(:)=um(:)/rhom(:)

  u(:) = eps(:)
  do i=1,ix-1
  um(i)=u(i)+dt/dx*( kappa(i  )/dx*(u(i+1)-u(i)) &
    &                - kappa(i-1)/dx*(u(i)-u(i-1)) )
  enddo
  epsm(:) = um(:)


  rho(:)=rhom(:)
  vx(:)=vxm(:)
  eps(:) = epsm(:)
  !----------------------------------------------------------------------|
  !     boundary condition
  rho(0)=rho(1)
  vx(0)=vx(1)
  eps(0) = eps(1)
  rho(ix)=rho(ix-1)
  vx(ix)=vx(ix-1)
  eps(ix) = eps(ix-1)
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
  !     epilogue
  !======================================================================|
  !  ending message
  write(*,formatstopmessage) ns,t
  write(*,*) '  ### normal stop ###'

  stop
end program main

