program main
!===================================================================!
! use modules
!===================================================================!
use vars, only : nprocs,my_rank,ierr
use init
use model
use datw
use cfl
use mlw
use bnd
use arvis
use mpi
use exc

!===================================================================!
! array definitions
!===================================================================!
implicit none
character(len=100) :: foutmsg,fstpmsg
integer,parameter :: mg=2
integer,parameter :: gx = 3600
integer,parameter :: mpix = 4
integer,parameter :: lix = 2*mg+gx/mpix
integer,parameter :: gix = 2*mg+gx
!integer,parameter :: ix = 2*mg+gx
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

real(8),parameter :: av = 5.d0 !lapidus artificial viscosity coef.
real(8) :: kx(lix)
real(8) :: time,dt,tend,tout,dtout
integer :: ns,nout ! ns = # of stage, nout = # of output

real(8) :: qu(lix),qfx(lix),qfxm(lix)

real(8),parameter :: pi = 4.d0*atan(1.d0)
!real(8) :: pii = 1.d0 / pi

real(8),parameter :: ro_floor=1.d-6
real(8),parameter :: pr_floor=1.d-7

! for file output
integer :: ifxx,ifro,ifvx,ifvy,ifbx,ifby,ifpr
integer :: write_cnt
integer(kind=MPI_OFFSET_KIND) :: disp,offset,file_end_pos
integer :: status(MPI_STATUS_SIZE)

integer :: ii!, jj


!===================================================================!
! prologue
!===================================================================!
! ready for MPI
!===================================================================!
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)



if(my_rank == 0) then
  if(nprocs /= mpix) then
    write(*,*) "number of processes is not consistent"
    write(*,*) "change mpix in main.f90 or MPINUM in makefile"
    write(*,*) "nprocs= ", nprocs
    write(*,*) "mpix= ", mpix
    call MPI_FINALIZE(ierr)
    stop
  end if

  if (mod(gx,mpix)/=0) then
    write(*,*) "gx/mpix is not integer"
    write(*,*) "gx=", gx
    write(*,*) "mpix= ", mpix
    call MPI_FINALIZE(ierr)
    stop
  end if 
end if

write_cnt = lix-2*mg
disp = int(write_cnt,kind=MPI_OFFSET_KIND)* 8_MPI_OFFSET_KIND
offset=int(my_rank,kind=MPI_OFFSET_KIND)*disp

!===================================================================!


! clean past data
if(my_rank == 0) then
  call dataclean()
end if

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
call mkgrd(xx,xm,dx,gix,lix,mg)
  !subroutine mkgrd(xx,xm,dx,gix,lix,mg)
call model_shocktube(ro,vx,vy,bx,by,pr,ez,etot,lix,xx,gm)


! mpi output


if(my_rank == 0) then
 call w0d(time,'dat/tt.dat')
end if


!write_cnt=lix-2*mg
!disp = int(write_cnt,kind=MPI_OFFSET_KIND)*8_MPI_OFFSET_KIND
!offset=int(my_rank,kind=MPI_OFFSET_KIND)*disp

call MPI_File_Open(MPI_Comm_World,'dat/x.dat',&
  MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifxx,ierr)
call MPI_File_Open(MPI_Comm_World,'dat/ro.dat',&
  MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifro,ierr)
call MPI_File_Open(MPI_Comm_World,'dat/vx.dat',&
  MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifvx,ierr)
call MPI_File_Open(MPI_Comm_World,'dat/vy.dat',&
  MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifvy,ierr)
call MPI_File_Open(MPI_Comm_World,'dat/bx.dat',&
  MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifbx,ierr)
call MPI_File_Open(MPI_Comm_World,'dat/by.dat',&
  MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifby,ierr)
call MPI_File_Open(MPI_Comm_World,'dat/pr.dat',&
  MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifpr,ierr)


call MPI_File_Write_at_All(ifxx,offset,xx(1+mg),write_cnt,&
  MPI_DOUBLE_PRECISION, status, ierr)
call MPI_File_Write_at_All(ifro,offset,ro(1+mg),write_cnt,&
  MPI_DOUBLE_PRECISION, status, ierr)
call MPI_File_Write_at_All(ifvx,offset,vx(1+mg),write_cnt,&
  MPI_DOUBLE_PRECISION, status, ierr)
call MPI_File_Write_at_All(ifvy,offset,vy(1+mg),write_cnt,&
  MPI_DOUBLE_PRECISION, status, ierr)
call MPI_File_Write_at_All(ifbx,offset,bx(1+mg),write_cnt,&
  MPI_DOUBLE_PRECISION, status, ierr)
call MPI_File_Write_at_All(ifby,offset,by(1+mg),write_cnt,&
  MPI_DOUBLE_PRECISION, status, ierr)
call MPI_File_Write_at_All(ifpr,offset,pr(1+mg),write_cnt,&
  MPI_DOUBLE_PRECISION, status, ierr)

call MPI_File_Close(ifxx,ierr)
call MPI_File_Close(ifro,ierr)
call MPI_File_Close(ifvx,ierr)
call MPI_File_Close(ifvy,ierr)
call MPI_File_Close(ifbx,ierr)
call MPI_File_Close(ifby,ierr)
call MPI_File_Close(ifpr,ierr)


! write initial condition
! call w0d(time,'dat/t.dat')
! call w1d(ix,xx,'dat/x.dat')
! call w1d(ix,ro,'dat/ro.dat')
! call w1d(ix,vx,'dat/vx.dat')
! !call w1d(ix,vy,'dat/vy.dat')
! !call w1d(ix,bx,'dat/bx.dat')
! !call w1d(ix,by,'dat/by.dat')
! call w1d(ix,pr,'dat/pr.dat')

nout= nout + 1
if(my_rank == 0) then
  write(*,foutmsg) ns, time, nout
end if

!===================================================================!
! main loop
!===================================================================!
do while (time < tend)
ns = ns + 1

! determine time step, dt
call calc_dt(dt,lix,dx,ro,vx,vy,bx,by,pr,gm)


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
call mlw1d1st(rom,qu,qfx,lix,dt,dx)

! x-momentum conservation
qu(:) = rx(:)
qfx(:) = ro(:)*vx(:)**2 + pr(:) + bb(:)/(8.d0*pi) &
  & - bx(:)**2/(4.d0*pi)
call mlw1d1st(rxm,qu,qfx,lix,dt,dx)

! y-momentum conservation
qu(:) = ry(:)
qfx(:) = ro(:)*vx(:)*vy(:) - bx(:)*by(:)/(4.d0*pi)
call mlw1d1st(rym,qu,qfx,lix,dt,dx)

! by equation
qu(:) = by(:)
qfx(:) = -ez(:)
call mlw1d1st(bym,qu,qfx,lix,dt,dx)

! energy conservation
qu(:) = etot(:)
qfx(:) = (etot(:) + pr(:) + bb(:)/(8.d0*pi))*vx(:) &
  & - (vx(:)*bx(:)+vy(:)*by(:))*bx(:)/(4.d0*pi)
call mlw1d1st(etotm,qu,qfx,lix,dt,dx)


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
call mlw1d2nd(dro,qfx,qfxm,lix,dt,dx)

! x-momentum conservation
qfx(:) = ro(:)*vx(:)**2 + pr(:) + bb(:)/(8.d0*pi) &
  & - bx(:)**2/(4.d0*pi)
qfxm(:) = rom(:)*vxm(:)**2 + prm(:) + bbm(:)/(8.d0*pi) &
  & - bxm(:)**2/(4.d0*pi)
call mlw1d2nd(drx,qfx,qfxm,lix,dt,dx)

! y-momentum conservation
qfx(:) = ro(:)*vx(:)*vy(:) - bx(:)*by(:)/(4.d0*pi)
qfxm(:)= rom(:)*vxm(:)*vym(:) - bxm(:)*bym(:)/(4.d0*pi)
call mlw1d2nd(dry,qfx,qfxm,lix,dt,dx)

! by equation
qfx(:) = -ez(:)
qfxm(:) = -ezm(:)
call mlw1d2nd(dby,qfx,qfxm,lix,dt,dx)

! energy conservation
qfx(:) = (etot(:) + pr(:) + bb(:)/(8.d0*pi))*vx(:) &
  & - (vx(:)*bx(:)+vy(:)*by(:))*bx(:)/(4.d0*pi)
qfxm(:) = (etotm(:) + prm(:) + bbm(:)/(8.d0*pi))*vxm(:) &
  & -(vxm(:)*bxm(:)+vym(:)*bym(:))*bxm(:)/(4.d0*pi)
call mlw1d2nd(detot,qfx,qfxm,lix,dt,dx)


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


! exchange margin data
call exchange_1d(ro,lix,mg)
call exchange_1d(vx,lix,mg)
call exchange_1d(vy,lix,mg)
!call exchange_1d(bx,lix,mg)
call exchange_1d(by,lix,mg)
call exchange_1d(pr,lix,mg)



!===================================================================!
! apply boudary condition
! left side
if(my_rank == 0) then
  call bnd1d(ro,0 ,lix)
  call bnd1d(vx,0 ,lix)
  call bnd1d(vy,0 ,lix)
  call bnd1d(by,0 ,lix)
  call bnd1d(pr,0 ,lix)
end if

!right side
if(my_rank == nprocs-1) then
  call bnd1d(ro,10,lix)
  call bnd1d(vx,10,lix)
  call bnd1d(vy,10,lix)
  call bnd1d(by,10,lix)
  call bnd1d(pr,10,lix)
end if

! wait MPi
call MPI_Barrier(MPI_Comm_World,ierr)
!===================================================================!
ez(:) = -vx(:)*by(:) + vy(:)*bx(:)
etot(:) = 0.5d0*ro(:)*(vx(:)**2+vy(:)**2) + pr(:)/(gm-1.d0) &
  & + (bx(:)**2+by(:)**2)/(8.d0*pi)


!===================================================================!
! apply artificial viscosity
!===================================================================!

call arvis1d(ro,dro,lix,dt,dx,av,vx,vy)
call arvis1d(rx,drx,lix,dt,dx,av,vx,vy)
call arvis1d(ry,dry,lix,dt,dx,av,vx,vy)
call arvis1d(by,dby,lix,dt,dx,av,vx,vy)
call arvis1d(etot,detot,lix,dt,dx,av,vx,vy)


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


! exchange margin data
call exchange_1d(ro,lix,mg)
call exchange_1d(vx,lix,mg)
call exchange_1d(vy,lix,mg)
!call exchange_1d(bx,lix,mg)
call exchange_1d(by,lix,mg)
call exchange_1d(pr,lix,mg)

!=======================================
! apply boudary condition

! Left side
if(my_rank == 0) then
  call bnd1d(ro,0 ,lix)
  call bnd1d(vx,0 ,lix)
  call bnd1d(vy,0 ,lix)
  call bnd1d(by,0 ,lix)
  call bnd1d(pr,0 ,lix)
end if

! Right side
if(my_rank == nprocs-1) then
  call bnd1d(ro,10,lix)
  call bnd1d(vx,10,lix)
  call bnd1d(vy,10,lix)
  call bnd1d(by,10,lix)
  call bnd1d(pr,10,lix)
end if

! wait mpi

call MPI_Barrier(MPI_Comm_World,ierr)

!=======================================




ez(:) = -vx(:)*by(:) + vy(:)*bx(:)
etot(:) = 0.5d0*ro(:)*(vx(:)**2+vy(:)**2) + pr(:)/(gm-1.d0) &
  & +(bx(:)**2+by(:)**2)/(8.d0*pi)

!==============================
! if ro, pr < 0 then apply flooring
!==============================
do ii = 1,lix
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
!if(.true.) then
  if(my_rank == 0) then
    call w0d(time,'dat/tt.dat')
  end if

!===================================================================!
! mpi output
!===================================================================!
  call MPI_File_Open(MPI_Comm_World,'dat/ro.dat',&
    MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifro,ierr)
  call MPI_File_Open(MPI_Comm_World,'dat/vx.dat',&
    MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifvx,ierr)
  call MPI_File_Open(MPI_Comm_World,'dat/vy.dat',&
    MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifvy,ierr)
  call MPI_File_Open(MPI_Comm_World,'dat/bx.dat',&
    MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifbx,ierr)
  call MPI_File_Open(MPI_Comm_World,'dat/by.dat',&
    MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifby,ierr)
  call MPI_File_Open(MPI_Comm_World,'dat/pr.dat',&
    MPI_Mode_Wronly+MPI_Mode_Create,MPI_Info_Null,ifpr,ierr)

  ! calculate file size and determine last position
  call MPI_File_Get_Size(ifro,file_end_pos,ierr)
  offset = file_end_pos + (int(my_rank,kind=MPI_OFFSET_KIND)*disp)
  ! write_cnt and disp is calculated in first output

  call MPI_File_Write_at_All(ifro,offset,ro(1+mg),write_cnt,&
    MPI_DOUBLE_PRECISION, status, ierr)
  call MPI_File_Write_at_All(ifvx,offset,vx(1+mg),write_cnt,&
    MPI_DOUBLE_PRECISION, status, ierr)
  call MPI_File_Write_at_All(ifvy,offset,vy(1+mg),write_cnt,&
    MPI_DOUBLE_PRECISION, status, ierr)
  call MPI_File_Write_at_All(ifbx,offset,bx(1+mg),write_cnt,&
    MPI_DOUBLE_PRECISION, status, ierr)
  call MPI_File_Write_at_All(ifby,offset,by(1+mg),write_cnt,&
    MPI_DOUBLE_PRECISION, status, ierr)
  call MPI_File_Write_at_All(ifpr,offset,pr(1+mg),write_cnt,&
    MPI_DOUBLE_PRECISION, status, ierr)

  call MPI_File_Close(ifro,ierr)
  call MPI_File_Close(ifvx,ierr)
  call MPI_File_Close(ifvy,ierr)
  call MPI_File_Close(ifbx,ierr)
  call MPI_File_Close(ifby,ierr)
  call MPI_File_Close(ifpr,ierr)




  !call w0d(time,'dat/t.dat')
  !call w1d(ix,ro,'dat/ro.dat')
  !call w1d(ix,vx,'dat/vx.dat')
  !call w1d(ix,vy,'dat/vy.dat')
  !call w1d(ix,bx,'dat/bx.dat')
  !call w1d(ix,by,'dat/by.dat')
  !call w1d(ix,pr,'dat/pr.dat')


    tout = tout + dtout
    nout = nout + 1
  if(my_rank == 0) then
    write(*,foutmsg) ns,time,nout
  end if
endif

! for debug
!if(ns >= 10) then
!  exit
!end if

enddo !main loop


!===================================================================!
! ending message
!===================================================================!
if (my_rank == 0) then
  write(*,fstpmsg) ns,time,nout
  write(*,*) '### NORMAL STOP ###'
end if

CALL MPI_FINALIZE(ierr)

stop
end program main
