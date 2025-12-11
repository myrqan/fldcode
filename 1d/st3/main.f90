program main
!===================================================================!
! use modules
!===================================================================!
use init
use model
use datw
use cfl
use mlw
use bnd
use arvis

!===================================================================!
! array definitions
!===================================================================!
implicit none
character(len=100) :: foutmsg,fstpmsg
integer :: mfile
integer,parameter :: margin=1
integer,parameter :: gx = 1200
integer,parameter :: ix = 2*margin+gx
real(8),parameter :: gm = 5.d0/3.d0
real(8) :: xx(ix),xm(ix),dx
real(8) :: ro(ix),vx(ix),pr(ix),ei(ix)
real(8) :: rx(ix)
real(8) :: rom(ix),vxm(ix),prm(ix),eim(ix)
real(8) :: rxm(ix)
real(8) :: dro(ix),drx(ix),dei(ix)

real(8),parameter :: av = 5.0
real(8) :: kx(ix)
real(8) :: time,dt,tend,tout,dtout
integer :: ns,nout ! ns = # of stage, nout = # of output

real(8) :: qu(ix),qfx(ix),qfxm(ix)


integer :: ii, jj


!===================================================================!
! prologue
!===================================================================!
! clean past data
call dataclean()

! set output format
!! terminal message (std output)
foutmsg='(1x, "[w] ", "step= ", i5, " time= ", e10.3, " nout= ", i3)'
fstpmsg='(1x, "[s] ", "step= ", i5, " time= ", e10.3, " nout= ", i3)'


!! initialize
xx(:) = 0.d0; xm(:) = 0.d0; dx = 0.d0
ro(:) = 0.d0; vx(:) = 0.d0; pr(:) = 0.d0; ei(:) = 0.d0
rx(:) = 0.d0
kx(:) = 0.d0

time = 0.d0; dt = 0.d0; tend = 0.d0; tout = 0.d0; dtout = 0.d0


! time control
tend = 0.14d0
dtout = 0.01d0
tout = dtout


! initialize counters
ns = 0; nout = 0


! setup numerlcal model
! (grid & initial condition)
call mkgrd(xx,xm,dx,ix,margin)
call model_shocktube(ro,vx,pr,ei,ix,xx,gm)

! write initial condition
call w0d(time,'dat/t.dat')
call w1d(ix,xx,'dat/x.dat')
call w1d(ix,ro,'dat/ro.dat')
call w1d(ix,vx,'dat/vx.dat')
call w1d(ix,pr,'dat/pr.dat')
call w1d(ix,ei,'dat/ei.dat')

nout= nout + 1
write(*,foutmsg) ns, time, nout

!===================================================================!
! main loop
!===================================================================!
do while (time < tend)
ns = ns + 1

! determine time step, dt
call calc_dt(dt,ix,dx,ro,vx,pr,gm)

! cleaning & setup
dro(:) = 0.d0; drx(:) = 0.d0; dei(:) = 0.d0

rx(:) = ro(:)*vx(:)


!===================================================================!
! 1st step
!===================================================================!

! mass conservation
qu(:) = ro(:)
qfx(:) = rx(:)
call mlw1d1st(rom,qu,qfx,ix,dt,dx)

! x-momentum conservation
qu(:) = rx(:)
qfx(:) = ro(:)*vx(:)**2+pr(:)
call mlw1d1st(rxm,qu,qfx,ix,dt,dx)

! energy conservation
qu(:) = ei(:)
qfx(:) = (ei(:)+pr(:))*vx(:)
call mlw1d1st(eim,qu,qfx,ix,dt,dx)


vxm(:) = rxm(:)/rom(:)
prm(:) = (gm-1.d0) * (eim(:) - 0.5d0*rom(:)*vxm(:)**2)


!===================================================================!
! second step
!===================================================================!
! mass conservation
qfx(:) = rx(:)
qfxm(:) = rxm(:)
call mlw1d2nd(dro,qfx,qfxm,ix,dt,dx)

! x-momentum conservation
qfx(:) = ro(:)*vx(:)**2 + pr(:)
qfxm(:) = rom(:)*vxm(:)**2 + prm(:)
call mlw1d2nd(drx,qfx,qfxm,ix,dt,dx)

! energy conservation
qfx(:) = (ei(:) + pr(:)) * vx(:)
qfxm(:) = (eim(:) + prm(:)) * vxm(:)
call mlw1d2nd(dei,qfx,qfxm,ix,dt,dx)


!===================================================================!
! update
!===================================================================!
ro(:) = ro(:) + dro(:)
rx(:) = rx(:) + drx(:)
ei(:) = ei(:) + dei(:)

vx(:) = rx(:) / ro(:)

! apply boudary condition
call bnd1d(ro,0 ,ix)
call bnd1d(ro,10,ix)
call bnd1d(vx,0 ,ix)
call bnd1d(vx,10,ix)
call bnd1d(ei,0 ,ix)
call bnd1d(ei,10,ix)

pr(:) = (gm-1.d0) * (ei(:) - 0.5d0*ro(:)*vx(:)**2)


!===================================================================!
! apply artificial viscosity
!===================================================================!

call arvis1d(ro,dro,ix,dt,dx,av,vx)
call arvis1d(rx,drx,ix,dt,dx,av,vx)
call arvis1d(ei,dei,ix,dt,dx,av,vx)


!===================================================================!
! update
!===================================================================!
ro(:) = ro(:) + dro(:)
rx(:) = rx(:) + drx(:)
ei(:) = ei(:) + dei(:)

vx(:) = rx(:) / ro(:)

! apply boudary condition
call bnd1d(ro,0 ,ix)
call bnd1d(ro,10,ix)
call bnd1d(vx,0 ,ix)
call bnd1d(vx,10,ix)
call bnd1d(ei,0 ,ix)
call bnd1d(ei,10,ix)

pr(:) = (gm-1.d0) * (ei(:) - 0.5d0*ro(:)*vx(:)**2)



time = time + dt

if (time >= tout) then
!if(.true.) then ! for debug
  call w0d(time,'dat/t.dat')
  call w1d(ix,ro,'dat/ro.dat')
  call w1d(ix,vx,'dat/vx.dat')
  call w1d(ix,pr,'dat/pr.dat')
  call w1d(ix,ei,'dat/ei.dat')
  tout = tout + dtout
  nout = nout + 1
  write(*,foutmsg) ns,time,nout
endif

! for debug
!if(ns >= 10) then
!  exit 
!end if
enddo !main loop


!===================================================================!
! ending message
!===================================================================!
write(*,fstpmsg) ns,time,nout
write(*,*) '### NORMAL STOP ###'


stop
end program main
