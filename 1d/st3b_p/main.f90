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
use mpi

!===================================================================!
! array definitions
!===================================================================!
implicit none
integer :: nprocs,my_rank,ierr
character(len=100) :: foutmsg,fstpmsg
integer :: mfile
integer,parameter :: margin=1
integer,parameter :: gx = 1000
integer,parameter :: mpix = 4
integer,parameter :: lix = 2*margin+gx/mpix
integer,parameter :: gix = 2*margin+gx
!integer,parameter :: ix = 2*margin+gx
real(8),parameter :: gm = 2.d0
!real(8),parameter :: gm = 5.d0/3.d0
real(8) :: xx(lix),xm(lix),dx
real(8) :: ro(lix),vx(lix),vy(lix),bx(lix),by(lix),pr(lix),&
  & etot(lix),ez(lix)
real(8) :: rx(lix),ry(lix)
real(8) :: rom(lix),vxm(lix),vym(lix),bxm(lix),bym(lix),prm(lix),&
  & etotm(lix),ezm(lix)
real(8) :: rxm(lix),rym(lix)
real(8) :: dro(lix),drx(lix),dry(lix),dby(lix),dbz(lix),detot(lix)
real(8) :: bb(lix),vv(lix),bbm(lix),vvm(lix)

real(8),parameter :: av = 2.0 !lapidus artificial viscosity coef.
real(8) :: kx(lix)
real(8) :: time,dt,tend,tout,dtout
integer :: ns,nout ! ns = # of stage, nout = # of output

real(8) :: qu(lix),qfx(lix),qfxm(lix)

real(8),parameter :: pi = 4.d0*atan(1.d0)
real(8) :: pii = 1.d0 / pi

real(8),parameter :: ro_floor=1.d-6
real(8),parameter :: pr_floor=1.d-7


integer :: ii, jj


!===================================================================!
! prologue
!===================================================================!
! ready for MPI
!===================================================================eo 
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)



if(nprocs /= mpix) then
  write(*,*) "number of processes is not consistent"
  write(*,*) "change header in main.f90 or ns in makefile"
  write(*,*) "nprocs= ", nprocs
  write(*,*) "mpix= ", mpix
  stop


! clean past data
call dataclean()

! set output format
!! terminal message (std output)
foutmsg='(1x, "[w] ", "step= ", i5, " time= ", e10.3, " nout= ", i3)'
fstpmsg='(1x, "[s] ", "step= ", i5, " time= ", e10.3, " nout= ", i3)'


!! initialize
xx(:) = 0.d0; xm(:) = 0.d0; dx = 0.d0
ro(:) = 0.d0; vx(:) = 0.d0; vy(:) = 0.d0; 
bx(:) = 0.d0; by(:) = 0.d0; pr(:) = 0.d0; etot(:) = 0.d0
etot(:) = 0.d0; ez(:) = 0.d0
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
call mkgrd(xx,xm,dx,gx,margin)
call model_shocktube(ro,vx,vy,bx,by,pr,ez,etot,ix,xx,gm)

! write initial condition
call w0d(time,'dat/t.dat')
call w1d(ix,xx,'dat/x.dat')
call w1d(ix,ro,'dat/ro.dat')
call w1d(ix,vx,'dat/vx.dat')
call w1d(ix,vy,'dat/vy.dat')
call w1d(ix,bx,'dat/bx.dat')
call w1d(ix,by,'dat/by.dat')
call w1d(ix,pr,'dat/pr.dat')

nout= nout + 1
write(*,foutmsg) ns, time, nout

!===================================================================!
! main loop
!===================================================================!
do while (time < tend)
ns = ns + 1

! determine time step, dt
call calc_dt(dt,ix,dx,ro,vx,vy,bx,by,pr,gm)

! cleaning & setup
dro(:) = 0.d0; drx(:) = 0.d0; dry(:) = 0.d0;
dby(:) = 0.d0; detot(:) = 0.d0;

rx(:) = ro(:)*vx(:)
ry(:) = ro(:)*vy(:)
bb(:) = bx(:)**2 + by(:)**2
vv(:) = vx(:)**2 + vy(:)**2


!===================================================================!
! 1st step
!===================================================================!

! mass conservation
qu(:) = ro(:)
qfx(:) = ro(:) * vx(:)
call mlw1d1st(rom,qu,qfx,ix,dt,dx)

! x-momentum conservation
qu(:) = rx(:)
qfx(:) = ro(:)*vx(:)**2 + pr(:) + bb(:)/(8.d0*pi) &
  & - bx(:)**2/(4.d0*pi)
call mlw1d1st(rxm,qu,qfx,ix,dt,dx)

! y-momentum conservation
qu(:) = ry(:)
qfx(:) = ro(:)*vx(:)*vy(:) - bx(:)*by(:)/(4.d0*pi)
call mlw1d1st(rym,qu,qfx,ix,dt,dx)

! by equation
qu(:) = by(:)
qfx(:) = -ez(:)
call mlw1d1st(bym,qu,qfx,ix,dt,dx)

! energy conservation
qu(:) = etot(:)
qfx(:) = (etot(:) + pr(:) + bb(:)/(8.d0*pi))*vx(:) &
  & - (vx(:)*bx(:)+vy(:)*by(:))*bx(:)/(4.d0*pi)
call mlw1d1st(etotm,qu,qfx,ix,dt,dx)


vxm(:) = rxm(:)/rom(:)
vym(:) = rym(:)/rom(:)
vvm(:) = vxm(:)**2+vym(:)**2
bxm(:) = bx(:)
bbm(:) = bxm(:)**2+bym(:)**2
prm(:) = (etotm(:)-0.5d0*rom(:)*vvm(:) -bbm(:)/(8.d0*pi))*(gm-1.d0)
ezm(:) = -vxm(:)*bym(:) + vym(:)*bxm(:)



!===================================================================!
! second step
!===================================================================!
! mass conservation
qfx(:) = rx(:)
qfxm(:) = rxm(:)
call mlw1d2nd(dro,qfx,qfxm,ix,dt,dx)

! x-momentum conservation
qfx(:) = ro(:)*vx(:)**2 + pr(:) + bb(:)/(8.d0*pi) &
  & - bx(:)**2/(4.d0*pi)
qfxm(:) = rom(:)*vxm(:)**2 + prm(:) + bbm(:)/(8.d0*pi) &
  & - bxm(:)**2/(4.d0*pi)
call mlw1d2nd(drx,qfx,qfxm,ix,dt,dx)

! y-momentum conservation
qfx(:) = ro(:)*vx(:)*vy(:) - bx(:)*by(:)/(4.d0*pi)
qfxm(:)= rom(:)*vxm(:)*vym(:) - bxm(:)*bym(:)/(4.d0*pi)
call mlw1d2nd(dry,qfx,qfxm,ix,dt,dx)

! by equation
qfx(:) = -ez(:)
qfxm(:) = -ezm(:)
call mlw1d2nd(dby,qfx,qfxm,ix,dt,dx)

! energy conservation
qfx(:) = (etot(:) + pr(:) + bb(:)/(8.d0*pi))*vx(:) &
  & - (vx(:)*bx(:)+vy(:)*by(:))*bx(:)/(4.d0*pi)
qfxm(:) = (etotm(:) + prm(:) + bbm(:)/(8.d0*pi))*vxm(:) &
  & -(vxm(:)*bxm(:)+vym(:)*bym(:))*bxm(:)/(4.d0*pi)
call mlw1d2nd(detot,qfx,qfxm,ix,dt,dx)


!===================================================================!
! update
!===================================================================!
ro(:) = ro(:) + dro(:)
rx(:) = rx(:) + drx(:)
ry(:) = ry(:) + dry(:)
by(:) = by(:) + dby(:)
etot(:) = etot(:) + detot(:)


vx(:) = rx(:)/ro(:)
vy(:) = ry(:)/ro(:)
vv(:) = vx(:)**2+vy(:)**2
bb(:) = bx(:)**2+by(:)**2
pr(:) = (etot(:)-0.5d0*ro(:)*vv(:)-bb(:)/(8.d0*pi))*(gm-1.d0)


! apply boudary condition
call bnd1d(ro,0 ,ix)
call bnd1d(ro,10,ix)
call bnd1d(vx,0 ,ix)
call bnd1d(vx,10,ix)
call bnd1d(vy,0 ,ix)
call bnd1d(vy,10,ix)
call bnd1d(by,0 ,ix)
call bnd1d(by,10,ix)
call bnd1d(pr,0 ,ix)
call bnd1d(pr,10,ix)

ez(:) = -vx(:)*by(:) + vy(:)*bx(:)
etot(:) = 0.5d0*ro(:)*(vx(:)**2+vy(:)**2) + pr(:)/(gm-1.d0) &
  & + (bx(:)**2+by(:)**2)/(8.d0*pi)


!===================================================================!
! apply artificial viscosity
!===================================================================!

call arvis1d(ro,dro,ix,dt,dx,av,vx,vy)
call arvis1d(rx,drx,ix,dt,dx,av,vx,vy)
call arvis1d(ry,dry,ix,dt,dx,av,vx,vy)
call arvis1d(by,dby,ix,dt,dx,av,vx,vy)
call arvis1d(etot,detot,ix,dt,dx,av,vx,vy)


!===================================================================!
! update
!===================================================================!

ro(:) = ro(:) + dro(:)
rx(:) = rx(:) + drx(:)
ry(:) = ry(:) + dry(:)
by(:) = by(:) + dby(:)
etot(:) = etot(:) + detot(:)

vx(:) = rx(:)/ro(:)
vy(:) = ry(:)/ro(:)
vv(:) = vx(:)**2+vy(:)**2
bb(:) = bx(:)**2+by(:)**2
pr(:) = (etot(:)-0.5d0*ro(:)*vv(:)-bb(:)/(8.d0*pi))*(gm-1.d0)


! apply boudary condition
call bnd1d(ro,0 ,ix)
call bnd1d(ro,10,ix)
call bnd1d(vx,0 ,ix)
call bnd1d(vx,10,ix)
call bnd1d(vy,0 ,ix)
call bnd1d(vy,10,ix)
call bnd1d(by,0 ,ix)
call bnd1d(by,10,ix)
call bnd1d(pr,0 ,ix)
call bnd1d(pr,10,ix)

ez(:) = -vx(:)*by(:) + vy(:)*bx(:)
etot(:) = 0.5d0*ro(:)*(vx(:)**2+vy(:)**2) + pr(:)/(gm-1.d0) &
  & +(bx(:)**2+by(:)**2)/(8.d0*pi)

!==============================
! if ro, pr < 0 then apply flooring
!==============================
do ii = 1,ix
if(ro(ii) < 0.d0) then
  ro(ii) = ro_floor
endif
if(pr(ii) < 0.d0) then
  pr(ii) = pr_floor
endif
enddo

! re-calculate etot
etot(:) = 0.5d0*ro(:)*(vx(:)**2+vy(:)**2) + pr(:)/(gm-1.d0) &
  & + (bx(:)**2+by(:)**2)/(8.d0*pi)



time = time + dt

if (time >= tout) then
  call w0d(time,'dat/t.dat')
  call w1d(ix,ro,'dat/ro.dat')
  call w1d(ix,vx,'dat/vx.dat')
  call w1d(ix,vy,'dat/vy.dat')
  call w1d(ix,bx,'dat/bx.dat')
  call w1d(ix,by,'dat/by.dat')
  call w1d(ix,pr,'dat/pr.dat')
  tout = tout + dtout
  nout = nout + 1
  write(*,foutmsg) ns,time,nout
endif
enddo !main loop


!===================================================================!
! ending message
!===================================================================!
write(*,fstpmsg) ns,time,nout
write(*,*) '### NORMAL STOP ###'


stop
end program main
