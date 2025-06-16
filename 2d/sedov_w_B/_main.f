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

  use, intrinsic :: IEEE_ARITHMETIC, only: IEEE_IS_NAN
  !======================================================================|
  !     array definitions
  !======================================================================|
  IMPLICIT NONE
  character(len=100) :: formatoutdata,formatoutmessage,formatstopmessage
  integer :: mfile
  integer,parameter :: margin = 1
  integer,parameter :: grid_x = 200
  integer,parameter :: grid_z = 200
  integer,parameter :: ix=2*margin+grid_x
  integer,parameter :: jx=2*margin+grid_z
  integer :: grid_ary(2)
  double precision :: x(ix,jx),xm(ix,jx),dx
  double precision :: z(ix,jx),zm(ix,jx),dz

  double precision :: u(ix,jx)
  double precision :: fx(ix,jx)
  double precision :: fz(ix,jx)
  double precision :: r(ix,jx)

  double precision :: t,dt,tend,tout,dtout
  integer :: ns, nd
  double precision, parameter :: pi = 4.d0 * atan(1.d0)

  double precision :: rho(ix,jx),rhon(ix,jx)
  double precision :: vx(ix,jx),vxn(ix,jx)
  double precision :: vy(ix,jx),vyn(ix,jx)
  double precision :: vz(ix,jx),vzn(ix,jx)
  double precision :: bx(ix,jx),bxn(ix,jx)
  double precision :: by(ix,jx),byn(ix,jx)
  double precision :: bz(ix,jx),bzn(ix,jx)
  double precision :: ex(ix,jx),exn(ix,jx)
  double precision :: ey(ix,jx),eyn(ix,jx)
  double precision :: ez(ix,jx),ezn(ix,jx)
  double precision :: p(ix,jx),pn(ix,jx)
  double precision :: eps(ix,jx),epsn(ix,jx)
  double precision :: ay(ix,jx),ayn(ix,jx)
  double precision :: b2(ix,jx),v2(ix,jx)
  double precision :: b2n(ix,jx),v2n(ix,jx)

  double precision :: dro(ix,jx),drvx(ix,jx),drvy(ix,jx),drvz(ix,jx),&
    de(ix,jx),dbx(ix,jx),dby(ix,jx),dbz(ix,jx),day(ix,jx)
  double precision :: ron(ix,jx),rvxn(ix,jx),rvyn(ix,jx),rvzn(ix,jx),&
    en(ix,jx)
  double precision :: ro(ix,jx),rvx(ix,jx),rvy(ix,jx),rvz(ix,jx),e(ix,jx)

  double precision :: kx(ix,jx),kz(ix,jx),qv
  double precision, parameter :: gamma = 5.d0/3.d0

  integer :: i,j
  !======================================================================|
  !     prologue
  !======================================================================|
  !----------------------------------------------------------------------|
  !   output format
  formatoutdata= '(1x,e13.5e3,1x,e13.5e3,1x,e13.5e3,1x,e13.5e3,1x,e13.5e3)'
  formatoutmessage= '(1x," write    ","step=",i8," t=",e10.3," nd =",i4)'
  formatstopmessage= '(1x," stop     ","step=",i8," t=",e10.3)'
  !----------------------------------------------------------------------|
  !   initialize

  call init2d(x); call init2d(xm)
  call init2d(u); call init2d(fx); call init2d(fz);
  call init2d(r); 

  call init2d(rho); call init2d(rhon)
  call init2d(vx); call init2d(vxn)
  call init2d(vy); call init2d(vyn)
  call init2d(vz); call init2d(vzn)
  call init2d(bx); call init2d(bxn)
  call init2d(by); call init2d(byn)
  call init2d(bz); call init2d(bzn)
  call init2d(ex); call init2d(exn)
  call init2d(ey); call init2d(eyn)
  call init2d(ez); call init2d(ezn)
  call init2d(p); call init2d(pn)
  call init2d(eps); call init2d(epsn)
  call init2d(ay); call init2d(ayn)
  
  call init2d(dro); call init2d(drvx); call init2d(drvy); call init2d(drvz)
  call init2d(ro);  call init2d(rvx);  call init2d(rvy);  call init2d(rvz)
  call init2d(ron); call init2d(rvxn); call init2d(rvyn); call init2d(rvzn);
  call init2d(e); call init2d(de); call init2d(en)
  call init2d(dbx); call init2d(dby); call init2d(dbz)
  call init2d(day)

  call init2d(kx); CALL init2d(kz)

  !----------------------------------------------------------------------|
  !   time control parameters
  tend=6.0d0  ! time for end of calculation
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

  call sedov2d_w_b(ix,jx,x,z,rho,vx,vy,vz,p,bx,by,bz)

  !! apply boundary conditions
  call bc_free_b(rho,ix,jx); call bc_free_u(rho,ix,jx)
  call bc_free_l(rho,ix,jx); call bc_free_r(rho,ix,jx)

  call bc_free_b(vx,ix,jx); call bc_free_u(vx,ix,jx)
  call bc_fix_l(vx,ix,jx); call bc_free_r(vx,ix,jx)

  call bc_free_b(vy,ix,jx); call bc_free_u(vy,ix,jx)
  call bc_fix_l(vy,ix,jx); call bc_free_r(vy,ix,jx)
  
  call bc_fix_b(vz,ix,jx); call bc_free_u(vz,ix,jx)
  call bc_free_l(vz,ix,jx); call bc_free_r(vz,ix,jx)

  call bc_fix_b(bx,ix,jx); call bc_free_u(bx,ix,jx)
  call bc_fix_l(bx,ix,jx); call bc_free_r(bx,ix,jx)

  call bc_fix_b(by,ix,jx); call bc_free_u(by,ix,jx)
  call bc_fix_l(by,ix,jx); call bc_free_r(by,ix,jx)

  call bc_free_b(bz,ix,jx); call bc_free_u(bz,ix,jx)
  call bc_free_l(bz,ix,jx); call bc_free_r(bz,ix,jx)

  call bc_free_b(ay,ix,jz); call bc_free_u(ay,ix,jx)
  call bc_free_l(ay,ix,jx); call bc_free_r(ay,ix,jx)

  call bc_free_b(p,ix,jx); call bc_free_u(p,ix,jx)
  call bc_free_l(p,ix,jx); call bc_free_r(p,ix,jx)

  !----------------------------------------------------------------------|
  !! BINARY OUTPUT
  ! calc eps
  eps(:,:) = p(:,:)/(gamma-1.d0) &
    + 0.5d0*rho(:,:)*(vx(:,:)**2 + vy(:,:)**2 + vz(:,:)**2)&
    + 0.5d0*(bx(:,:)**2 + by(:,:)**2 + bz(:,:)**2);
  ex(:,:) = -vy(:,:)*bz(:,:) + vz(:,:)*by(:,:)
  ey(:,:) = -vz(:,:)*bx(:,:) + vx(:,:)*bz(:,:)
  ez(:,:) = -vx(:,:)*by(:,:) + vy(:,:)*bx(:,:)


  grid_ary(1)=ix;grid_ary(2)=jx
  CALL put1dint(9,'ixjx.dac',grid_ary)
  call put2dreal(11,'x.dac',x)
  call put2dreal(12,'z.dac',z)
  call put0dreal(10,'t.dac',t)
  call put2dreal(16,'vx.dac',vx)
  call put2dreal(17,'vz.dac',vy)
  call put2dreal(18,'vz.dac',vz)
  call put2dreal(20,'rho.dac',rho)
  call put2dreal(21,'p.dac',p)
  call put2dreal(22,'eps.dac',eps)
  call put2dreal(26,'bx.dac',bx)
  call put2dreal(27,'by.dac',by)
  call put2dreal(28,'bz.dac',bz)
  call put2dreal(30,'ay.dac',ay)

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

  CALL calc_dt(rho,vx,vy,vz,p,bx,by,bz,gamma,ix,jx,dt,dx,dz)
  !---------------------------------------------------------------------|
  ! FIRST STEP
  !----------------------------------------------------------------------|
  ! initialize d** 
  CALL init2d(dro); CALL init2d(drvx); CALL init2d(drvz); CALL init2d(de)
  CALL init2d(ron); CALL init2d(rvxn); CALL init2d(rvzn); CALL init2d(en)
  CALL init2d(rvyn); CALL init2d(bxn); CALL init2d(byn); CALL init2d(bzn)
  CALL init2d(ayn)

  b2(:,:) = bx(:,:)**2 + by(:,:)**2 + bz(:,:)**2
  v2(:,:) = vx(:,:)**2 + vy(:,:)**2 + vz(:,:)**2


  ro(:,:)=rho(:,:)
  fx(:,:)=rho(:,:)*vx(:,:)
  fz(:,:)=rho(:,:)*vz(:,:)
  r(:,:)=-fx(:,:)/x(:,:)
  CALL mlw2d1st(ro,fx,fz,ron,dro,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,ron,dro,ix,jx,dt,dx,dz)

  rvx(:,:)= rho(:,:)*vx(:,:)
  fx(:,:) = rho(:,:)*vx(:,:)**2 + p(:,:) + 0.5d0*b2(:,:) - bx(:,:)**2
  fz(:,:) = rho(:,:)*vx(:,:)*vz(:,:) - bx(:,:)*bx(:,:)
  r(:,:)  =-(rho(:,:)*(vx(:,:)**2-vy(:,:)**2) &
    - (bx(:,:)**2-by(:,:)**2))/x(:,:)
  CALL mlw2d1st(rvx,fx,fz,rvxn,drvx,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,rvxn,drvx,ix,jx,dt,dx,dz)

  rvy(:,:)= rho(:,:)        *vy(:,:)
  fx(:,:) = rho(:,:)*vx(:,:)*vy(:,:) - bx(:,:)*by(:,:)
  fz(:,:) = rho(:,:)*vz(:,:)*vy(:,:) - bz(:,:)*by(:,:)
  r(:,:)  =-2.d0 * fx(:,:) / x(:,:)
  CALL mlw2d1st(rvy,fx,fz,rvyn,drvy,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,rvyn,drvy,ix,jx,dt,dx,dz)

  rvz(:,:) = rho(:,:) * vz(:,:)
  fx(:,:) = rho(:,:)*vx(:,:)*vz(:,:) - bx(:,:)*bz(:,:)
  fz(:,:) = rho(:,:)*vz(:,:)**2 + p(:,:) + 0.5d0 * b2(:,:) - bz(:,:)**2
  r(:,:) = - fx(:,:) / x(:,:)
  CALL mlw2d1st(rvz,fx,fz,rvzn,drvz,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,rvzn,drvz,ix,jx,dt,dx,dz)

  e(:,:)  =  eps(:,:)
  fx(:,:) = (gamma/(gamma-1.d0)*p(:,:) + 0.5d0 * rho(:,:) * v2(:,:)) * vx(:,:)&
    + (ey(:,:)*bz(:,:) - ez(:,:)*by(:,:))
  fy(:,:) = (gamma-(gamma-1.d0)*p(:,:) + 0.5d0 * rho(:,:) * v2(:,:)) * vz(:,:)&
    + (ex(:,:)*by(:,:) - ey(:,:)*bx(:,:))
  r(:,:) = -fx(:,:)/x(:,:)
  CALL mlw2d1st(e,fx,fz,en,de,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,en,de,ix,jx,dt,dx,dz)




  rhon(:,:)=ron(:,:)
  vxn(:,:)=rvxn(:,:)/ron(:,:)
  vzn(:,:)=rvzn(:,:)/ron(:,:)
  epsn(:,:)=en(:,:)


  pn(:,:) = (gamma-1.d0) * (epsn(:,:)-0.5d0*rhon(:,:)*(vxn(:,:)**2+vzn(:,:)**2))


  !----------------------------------------------------------------------|
  ! SECOND STEP
  !----------------------------------------------------------------------|

  fx(:,:)=rhon(:,:)*vxn(:,:)
  fz(:,:)=rhon(:,:)*vzn(:,:)
  r(:,:)=-rhon(:,:)*vxn(:,:)/xm(:,:)
  CALL mlw2d2nd(fx,fz,dro,ix,jx,dt,dx,dz)
  CALL mlw2dsrc2nd(dro,r,ix,jx,dt,dx,dz)

  fx(:,:)=rhon(:,:)*vxn(:,:)**2+pn(:,:)
  fz(:,:)=rhon(:,:)*vxn(:,:)*vzn(:,:)
  r(:,:)=-rhon(:,:)*vxn(:,:)**2/xm(:,:)
  CALL mlw2d2nd(fx,fz,drvx,ix,jx,dt,dx,dz)
  CALL mlw2dsrc2nd(drvx,r,ix,jx,dt,dx,dz)

  fx(:,:)=rhon(:,:)*vxn(:,:)*vzn(:,:)
  fz(:,:)=rhon(:,:)*vzn(:,:)**2+pn(:,:)
  r(:,:)=-rhon(:,:)*vxn(:,:)*vzn(:,:)/xm(:,:)
  CALL mlw2d2nd(fx,fz,drvz,ix,jx,dt,dx,dz)
  CALL mlw2dsrc2nd(drvz,r,ix,jx,dt,dx,dz)

  fx(:,:)=(en(:,:)+pn(:,:))*vxn(:,:)
  fz(:,:)=(en(:,:)+pn(:,:))*vzn(:,:)
  r(:,:)=-(en(:,:)+pn(:,:))*vxn(:,:)/xm(:,:)
  CALL mlw2d2nd(fx,fz,de,ix,jx,dt,dx,dz)
  CALL mlw2dsrc2nd(de,r,ix,jx,dt,dx,dz)



  !=====================
  ! ARTIFICIAL VISCOSITY
  !=====================

  qv = 3.d0

  if(ns == 1) then
    call qv_param(qv)
  endif

  call calc_coef_k(qv,kx,kz,vx,vz,ix,jx,dx,dz)

  ro(:,:)=rho(:,:)
  CALL arvis2d(kx,kz,ro,dro,ix,jx,dt,dx,dz)

  rvx(:,:) = rho(:,:)*vx(:,:)
  CALL arvis2d(kx,kz,rvx,drvx,ix,jx,dt,dx,dz)

  rvz(:,:) = rho(:,:)*vz(:,:)
  CALL arvis2d(kx,kz,rvz,drvz,ix,jx,dt,dx,dz)

  e(:,:)=eps(:,:)
  CALL arvis2d(kx,kz,e,de,ix,jx,dt,dx,dz)

  ro(:,:)=ro(:,:)+dro(:,:)
  rvx(:,:)=rvx(:,:)+drvx(:,:)
  rvz(:,:)=rvz(:,:)+drvz(:,:)
  e(:,:)=e(:,:)+de(:,:)

  rho(:,:)=ro(:,:)
  vx(:,:)=rvx(:,:)/rho(:,:)
  vz(:,:)=rvz(:,:)/rho(:,:)
  eps(:,:)=e(:,:)

  p(:,:) = (gamma - 1.d0) * (eps(:,:) - 0.5d0*rho(:,:)*(vx(:,:)**2+vz(:,:)**2))

  !==============================
  !!     boundary condition
  !==============================

  call bc_free_b(rho,ix,jx); call bc_free_u(rho,ix,jx)
  call bc_free_l(rho,ix,jx); call bc_free_r(rho,ix,jx)

  call bc_fix_b(vx,ix,jx); call bc_free_u(vx,ix,jx)
  call bc_free_b(vx,ix,jx); call bc_free_r(vx,ix,jx)
  
  call bc_fix_b(vz,ix,jx); call bc_free_u(vz,ix,jx)
  call bc_free_l(vz,ix,jx); call bc_free_r(vz,ix,jx)

  call bc_free_b(eps,ix,jx); call bc_free_u(eps,ix,jx)
  call bc_free_l(eps,ix,jx); call bc_free_r(eps,ix,jx)
  
  p(:,:) = (gamma - 1.d0) * (eps(:,:) - 0.5d0*rho(:,:)*(vx(:,:)**2+vz(:,:)**2))

  call bc_free_b(p,ix,jx); call bc_free_u(p,ix,jx)
  call bc_free_l(p,ix,jx); call bc_free_r(p,ix,jx)
 

  t=t+dt

  !----------------------------------------------------------------------|
  !     data output 
  !----------------------------------------------------------------------|
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
