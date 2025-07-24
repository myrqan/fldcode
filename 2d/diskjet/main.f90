PROGRAM main
  !========================================
  !! IMPORT MODULES
  !========================================
  USE arvis
  USE boundary2d
  USE cfl
  USE calculate_variables
  USE file_output
  USE mlw
  USE model
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN

  !========================================
  !! VARIABLES DEFINITION
  !========================================
  IMPLICIT NONE
  DOUBLE PRECISION,PARAMETER :: pi = 4.d0*ATAN(1.d0)
  CHARACTER(100) :: foutmsg,fstopmsg
  INTEGER,PARAMETER :: margin = 2
  INTEGER,PARAMETER :: grid_x = 150
  INTEGER,PARAMETER :: grid_z = 200
  INTEGER,PARAMETER :: ix = 2*margin+grid_x
  INTEGER,PARAMETER :: jx = 2*margin+grid_z
  INTEGER :: ixjx(3) = (/ix,jx,margin/)

  DOUBLE PRECISION :: x(ix,jx),xm(ix,jx),dx
  DOUBLE PRECISION :: z(ix,jx),zm(ix,jx),dz
  DOUBLE PRECISION :: gx(ix,jx),gxm(ix,jx),gz(ix,jx),gzm(ix,jx)

  DOUBLE PRECISION :: u(ix,jx),fx(ix,jx),fz(ix,jx),r(ix,jx)
  DOUBLE PRECISION :: ro(ix,jx),vx(ix,jx),vy(ix,jx),vz(ix,jx),&
    bx(ix,jx),by(ix,jx),bz(ix,jx),ex(ix,jx),ey(ix,jx),ez(ix,jx),&
    p(ix,jx),etot(ix,jx),ay(ix,jx)
  DOUBLE PRECISION :: ron(ix,jx),vxn(ix,jx),vyn(ix,jx),vzn(ix,jx),&
    bxn(ix,jx),byn(ix,jx),bzn(ix,jx),exn(ix,jx),eyn(ix,jx),ezn(ix,jx),&
    pn(ix,jx),etotn(ix,jx)
  DOUBLE PRECISION :: rvx(ix,jx),rvy(ix,jx),rvz(ix,jx)
  DOUBLE PRECISION :: rvxn(ix,jx),rvyn(ix,jx),rvzn(ix,jx)
  DOUBLE PRECISION :: dro(ix,jx),drvx(ix,jx),drvy(ix,jx),drvz(ix,jx),&
    detot(ix,jx),dbx(ix,jx),dby(ix,jx),dbz(ix,jx)
  DOUBLE PRECISION :: v2(ix,jx),v2n(ix,jx),b2(ix,jx),b2n(ix,jx)
  DOUBLE PRECISION :: kx(ix,jx),kz(ix,jx),qv
  DOUBLE PRECISION,PARAMETER :: gm = 5.d0/3.d0

  INTEGER :: ns,nd
  DOUBLE PRECISION :: t,dt,tend,tout,dtout

  INTEGER :: i,j

  !========================================
  !!  setup & initialize
  !========================================
  !!  make stdout format
  !========================================

  foutmsg='(1x," WRITE  ","STEP=  ",i7," time=",e10.3," nd=",i4)'
  fstopmsg='(1x," STOP  ","STEP=  ",i7," time=",e10.3)'

  !========================================
  !! initialize variables
  !========================================

  x(:,:)=0.d0;xm(:,:)=0.d0;dx=0.d0;
  z(:,:)=0.d0;zm(:,:)=0.d0;dz=0.d0;
  gx(:,:)=0.d0;gxm(:,:)=0.d0;gz(:,:)=0.d0;gzm(:,:)=0.d0
  u(:,:)=0.d0;fx(:,:)=0.d0;fz(:,:)=0.d0;r(:,:)=0.d0
  ro(:,:)=0.d0;vx(:,:)=0.d0;vy(:,:)=0.d0;vz(:,:)=0.d0;
  bx(:,:)=0.d0;by(:,:)=0.d0;bz(:,:)=0.d0;
  ex(:,:)=0.d0;ey(:,:)=0.d0;ez(:,:)=0.d0;
  p(:,:)=0.d0;etot(:,:)=0.d0;ay(:,:)=0.d0
  ron(:,:)=0.d0;vxn(:,:)=0.d0;vyn(:,:)=0.d0;vzn(:,:)=0.d0;
  bxn(:,:)=0.d0;byn(:,:)=0.d0;bzn(:,:)=0.d0;
  exn(:,:)=0.d0;eyn(:,:)=0.d0;ezn(:,:)=0.d0;
  pn(:,:)=0.d0;etotn(:,:)=0.d0
  rvx(:,:)=0.d0;rvy(:,:)=0.d0;rvz(:,:)=0.d0;
  rvxn(:,:)=0.d0;rvyn(:,:)=0.d0;rvzn(:,:)=0.d0
  dro(:,:)=0.d0;drvx(:,:)=0.d0;drvy(:,:)=0.d0;drvz(:,:)=0.d0;
  detot(:,:)=0.d0;dbx(:,:)=0.d0;dby(:,:)=0.d0;dbz(:,:)=0.d0
  v2(:,:)=0.d0;b2(:,:)=0.d0;v2n(:,:)=0.d0;b2n(:,:)=0.d0
  kx(:,:)=0.d0;kz(:,:)=0.d0;qv=0.d0


  !========================================
  !!  time parameters
  !========================================
  tend = 6.0d0 !! end of calculation
  dtout = 0.1d-0!! time interval for output
  t = 0.d0
  tout = 0.d0
  ns = 0 !! # of steps
  nd = 0 !! # of printed data

  !========================================
  !! initialize parameter files & data files
  !========================================
  CALL clean()
  CALL put_param_int("margin:",margin)
  CALL put_param_int("ix: ",ix)
  CALL put_param_int("jx: ",jx)
  CALL write_1d_int("ixjx.dat",ixjx)


  !========================================
  !! setup models
  !========================================

  CALL make_grid(x,xm,dx,z,zm,dz,ix,jx,margin)
  CALL diskjet(ro,vx,vy,vz,bx,by,bz,p,gx,gxm,gz,gzm,x,z,ix,jx,dx,dz)
  !CALL sedov_with_b(ro,vx,vy,vz,bx,by,bz,p,x,z,ix,jx)
  CALL calc_ay(x,z,bx,by,bz,ix,jx,ay)
  !gx(:,:)=0.d0;gz(:,:)=0.d0;gxm(:,:)=0.d0;gzm(:,:)=0.d0


  ! Apply boundary conditions
  !
  ! (zmax) |-------------|
  !        |      2      |
  !        |             |
  !        | 3         1 |
  !        |             |
  !        |      0      |
  ! (zmin) |-------------|
  !      (xmin)        (xmax)
  !
  !  how to make grid:
  !  e.g.) margin = 1
  !
  !  (mgn)
  !
  !  (1,3) | (2,3)   (3,3)
  !        |
  !  (1,2) | (2,2)   (3,2)
  !        ---------------   
  !  (1,1)   (2,1)   (3,1)  (margin)
  !
  !  0 to 3 represents boundary surface
  !
  ! name : bc(val,num1,num2,margin,ix,jx)
  ! num1 : surface number (0 to 3)
  ! num2 : symmetric same sign(0) or symmetric different sign(1) or free(2)


  ! (apply BC for ro,v,b,p)
  !0220
  CALL bc(ro,0,0,margin,ix,jx)
  CALL bc(ro,1,2,margin,ix,jx)
  CALL bc(ro,2,2,margin,ix,jx)
  CALL bc(ro,3,0,margin,ix,jx)

  !0221
  CALL bc(vx,0,0,margin,ix,jx)
  CALL bc(vx,1,2,margin,ix,jx)
  CALL bc(vx,2,2,margin,ix,jx)
  CALL bc(vx,3,1,margin,ix,jx)

  !0221
  CALL bc(vy,0,0,margin,ix,jx)
  CALL bc(vy,1,2,margin,ix,jx)
  CALL bc(vy,2,2,margin,ix,jx)
  CALL bc(vy,3,1,margin,ix,jx)

  !1220
  CALL bc(vz,0,1,margin,ix,jx)
  CALL bc(vz,1,2,margin,ix,jx)
  CALL bc(vz,2,2,margin,ix,jx)
  CALL bc(vz,3,0,margin,ix,jx)

  !1221
  CALL bc(bx,0,1,margin,ix,jx)
  CALL bc(bx,1,2,margin,ix,jx)
  CALL bc(bx,2,2,margin,ix,jx)
  CALL bc(bx,3,1,margin,ix,jx)

  !1221
  CALL bc(by,0,1,margin,ix,jx)
  CALL bc(by,1,2,margin,ix,jx)
  CALL bc(by,2,2,margin,ix,jx)
  CALL bc(by,3,1,margin,ix,jx)

  !0220
  CALL bc(bz,0,0,margin,ix,jx)
  CALL bc(bz,1,2,margin,ix,jx)
  CALL bc(bz,2,2,margin,ix,jx)
  CALL bc(bz,3,0,margin,ix,jx)

  !0220
  CALL bc(p,0,0,margin,ix,jx)
  CALL bc(p,1,2,margin,ix,jx)
  CALL bc(p,2,2,margin,ix,jx)
  CALL bc(p,3,0,margin,ix,jx)

 !! apply boudary conditions for g{x,z}(m)
 !0221
  CALL bc(gx,0,0,margin,ix,jx)
  CALL bc(gx,1,2,margin,ix,jx)
  CALL bc(gx,2,2,margin,ix,jx)
  CALL bc(gx,3,1,margin,ix,jx)
  !1220
  CALL bc(gz,0,1,margin,ix,jx)
  CALL bc(gz,1,2,margin,ix,jx)
  CALL bc(gz,2,2,margin,ix,jx)
  CALL bc(gz,3,0,margin,ix,jx)
  !0221
  CALL bc(gxm,0,0,margin,ix,jx)
  CALL bc(gxm,1,2,margin,ix,jx)
  CALL bc(gxm,2,2,margin,ix,jx)
  CALL bc(gxm,3,1,margin,ix,jx)
  !1220 
  CALL bc(gzm,0,1,margin,ix,jx)
  CALL bc(gzm,1,2,margin,ix,jx)
  CALL bc(gzm,2,2,margin,ix,jx)
  CALL bc(gzm,3,0,margin,ix,jx)


  CALL calc_e(bx,by,bz,vx,vy,vz,ex,ey,ez)
  CALL calc_etot(ro,p,vx,vy,vz,bx,by,bz,etot,gm,ix,jx)

  !========================================
  !! data output (1st time)
  !========================================
  CALL write_0d_dble("t.dat",t)
  CALL write_2d_dble("x.dat",x)
  CALL write_2d_dble("z.dat",z)
  CALL write_2d_dble("ro.dat",ro)
  CALL write_2d_dble("vx.dat",vx)
  CALL write_2d_dble("vy.dat",vy)
  CALL write_2d_dble("vz.dat",vz)
  CALL write_2d_dble("bx.dat",bx)
  CALL write_2d_dble("by.dat",by)
  CALL write_2d_dble("bz.dat",bz)
  CALL write_2d_dble("p.dat",p)
  CALL write_2d_dble('etot.dat',etot)
  CALL write_2d_dble("ay.dat",ay)
  write(*,foutmsg) ns,t,nd
  nd = nd + 1
  tout = tout + dtout

  !========================================
  !!  MAIN LOOP (TIME INTEGRATION)
  !========================================
  main_loop:&
    do while(t < tend)
    !do while(ns < 10)
!    if(mod(ns,10)==0) then
!      write(*,*) "ns=", ns, "t=", t
!    endif

  !! Calculate dt using CFL condition
  CALL calc_dt(ro,vx,vy,vz,p,bx,by,bz,gm,ix,jx,dt,dx,dz,margin)
  

  !! initialize variables
  dro(:,:)=0.d0;drvx(:,:)=0.d0;drvy(:,:)=0.d0;drvz(:,:)=0.d0
  detot(:,:)=0.d0;dbx(:,:)=0.d0;dby(:,:)=0.d0;dbz(:,:)=0.d0
  ron(:,:)=0.d0;rvxn(:,:)=0.d0;rvyn(:,:)=0.d0;rvzn(:,:)=0.d0
  etotn(:,:)=0.d0;bxn(:,:)=0.d0;byn(:,:)=0.d0;bzn(:,:)=0.d0
  exn(:,:)=0.d0;eyn(:,:)=0.d0;ezn(:,:)=0.d0

  rvx(:,:) = ro(:,:)*vx(:,:)
  rvy(:,:) = ro(:,:)*vy(:,:)
  rvz(:,:) = ro(:,:)*vz(:,:)
  b2(:,:) = bx(:,:)**2+by(:,:)**2+bz(:,:)**2
  v2(:,:) = vx(:,:)**2+vy(:,:)**2+vz(:,:)**2

  !========================================
  !! 1st STEP
  !========================================
  !----- density-----!
  fx(:,:) = ro(:,:)*vx(:,:)
  fz(:,:) = ro(:,:)*vz(:,:)
  r(:,:) = -fx(:,:) / x(:,:)
  CALL mlw2d1st(ro,fx,fz,ron,dro,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,ron,dro,ix,jx,dt,dx,dz)

  !----- x-momentum -----!
  fx(:,:) = ro(:,:)*vx(:,:)**2+p(:,:)&
    +0.5d0*(-bx(:,:)**2+by(:,:)**2+bz(:,:)**2)
  fz(:,:) = ro(:,:)*vz(:,:)*vx(:,:)-bz(:,:)*bx(:,:)
  r(:,:) = - (ro(:,:)*(vx(:,:)**2-vy(:,:)**2) &
    - (bx(:,:)**2-by(:,:)**2))/x(:,:) + ro(:,:)*gx(:,:)
  CALL mlw2d1st(rvx,fx,fz,rvxn,drvx,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,rvxn,drvx,ix,jx,dt,dx,dz)

  !----- y-momentum -----!
  fx(:,:) = ro(:,:)*vx(:,:)*vy(:,:)-bx(:,:)*by(:,:)
  fz(:,:) = ro(:,:)*vz(:,:)*vy(:,:)-bz(:,:)*by(:,:)
  r(:,:) = -2.d0 * fx(:,:)/x(:,:)
  CALL mlw2d1st(rvy,fx,fz,rvyn,drvy,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,rvyn,drvy,ix,jx,dt,dx,dz)

  !----- z-momentum -----!
  fx(:,:) = ro(:,:)*vx(:,:)*vz(:,:)-bx(:,:)*bz(:,:)
  fz(:,:) = ro(:,:)*vz(:,:)**2+p(:,:)&
    +0.5d0*(bx(:,:)**2+by(:,:)**2-bz(:,:)**2)
  r(:,:) = - fx(:,:)/x(:,:) + ro(:,:)*gz(:,:)
  CALL mlw2d1st(rvz,fx,fz,rvzn,drvz,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,rvzn,drvz,ix,jx,dt,dx,dz)

  !----- energy -----!
  fx(:,:) = (gm/(gm-1.d0)*p(:,:) + 0.5d0*ro(:,:)*v2(:,:))*vx(:,:)& 
    + (ey(:,:)*bz(:,:)-ez(:,:)*by(:,:))
  fz(:,:) = (gm/(gm-1.d0)*p(:,:) + 0.5d0*ro(:,:)*v2(:,:))*vz(:,:)&
    + (ex(:,:)*by(:,:)-ey(:,:)*bx(:,:))
  r(:,:) = -fx(:,:) / x(:,:) &
    + ro(:,:)*gx(:,:)*vx(:,:) + ro(:,:)*gz(:,:)*vz(:,:)
  CALL mlw2d1st(etot,fx,fz,etotn,detot,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,etotn,detot,ix,jx,dt,dx,dz)

  fx(:,:) = 0.d0
  fz(:,:) = -ey(:,:)
  r(:,:) = 0.d0
  CALL mlw2d1st(bx,fx,fz,bxn,dbx,ix,jx,dt,dx,dz)
  ! CALL mlw2dsrc1st()

  fx(:,:) = -ez(:,:)
  fz(:,:) = ex(:,:)
  r(:,:) = 0.d0
  CALL mlw2d1st(by,fx,fz,byn,dby,ix,jx,dt,dx,dz)
  ! CALL mlw2dsrc1st()

  fx(:,:) = ey(:,:)
  fz(:,:) = 0.d0
  r(:,:) = - fx(:,:) / x(:,:)
  CALL mlw2d1st(bz,fx,fz,bzn,dbz,ix,jx,dt,dx,dz)
  CALL mlw2dsrc1st(r,bzn,dbz,ix,jx,dt,dx,dz)

  !  get half step value
  vxn(:,:) = rvxn(:,:) / ron(:,:)
  vyn(:,:) = rvyn(:,:) / ron(:,:)
  vzn(:,:) = rvzn(:,:) / ron(:,:)
  
  ! get pn & e{x,y,z}n
  CALL calc_p(ron,pn,vxn,vyn,vzn,bxn,byn,bzn,etotn,gm,ix,jx)
  CALL calc_e(bxn,byn,bzn,vxn,vyn,vzn,exn,eyn,ezn)


  !!!!!
  do j = 1,jx
  do i = 1,ix
    if(pn(i,j) < 0.d0) then
      !write(*,*) "pressure is negative at (i,j)=", i, j
      !write(*,*) "end of 1st step" 
      !pn(i,j) = 1.d-3
      !stop 
    endif
  enddo
  enddo
  !!!!!

  
  !========================================
  !! 2nd STEP
  !========================================

  !----- density -----!
  fx(:,:) = ron(:,:)*vxn(:,:)
  fz(:,:) = ron(:,:)*vzn(:,:)
  r(:,:) = -fx(:,:) / xm(:,:)
  CALL mlw2d2nd(fx,fz,dro,ix,jx,dt,dx,dz)
  CALL mlw2dsrc2nd(dro,r,ix,jx,dt,dx,dz)

  !----- x-momentum -----!
  fx(:,:) = ron(:,:)*vxn(:,:)**2+pn(:,:)&
    +0.5d0*(-bxn(:,:)**2+byn(:,:)**2+bzn(:,:)**2)
  fz(:,:) = ron(:,:)*vzn(:,:)*vxn(:,:)-bzn(:,:)*bxn(:,:)
  r(:,:) = - (ron(:,:)*(vxn(:,:)**2-vyn(:,:)**2)&
    -(bxn(:,:)**2-byn(:,:)**2)) / xm(:,:) + ron(:,:)*gxm(:,:)
  CALL mlw2d2nd(fx,fz,drvx,ix,jx,dt,dx,dz)
  CALL mlw2dsrc2nd(drvx,r,ix,jx,dt,dx,dz)

  !----- y-momentum -----!
  fx(:,:) = ron(:,:)*vxn(:,:)*vyn(:,:)-bxn(:,:)*byn(:,:)
  fz(:,:) = ron(:,:)*vzn(:,:)*vyn(:,:)-bzn(:,:)*byn(:,:)
  r(:,:) = -2.d0 * fx(:,:)/xm(:,:)
  CALL mlw2d2nd(fx,fz,drvy,ix,jx,dt,dx,dz)
  CALL mlw2dsrc2nd(drvy,r,ix,jx,dt,dx,dz)

  !----- z-momentum -----!
  fx(:,:) = ron(:,:)*vxn(:,:)*vzn(:,:)-bxn(:,:)*bzn(:,:)
  fz(:,:) = ron(:,:)*vzn(:,:)**2+pn(:,:)&
    +0.5d0*(bxn(:,:)**2+byn(:,:)**2-bzn(:,:)**2)
  r(:,:) = - fx(:,:)/xm(:,:) + ron(:,:)*gzm(:,:)
  CALL mlw2d2nd(fx,fz,drvz,ix,jx,dt,dx,dz)
  CALL mlw2dsrc2nd(drvz,r,ix,jx,dt,dx,dz)

  v2n(:,:) = vxn(:,:)**2 + vyn(:,:)**2 + vzn(:,:)**2

  !----- energy -----!
  fx(:,:) = (gm/(gm-1.d0)*pn(:,:) + 0.5d0*ron(:,:)*v2n(:,:))*vxn(:,:)&
    + (eyn(:,:)*bzn(:,:)-ezn(:,:)*byn(:,:))
  fz(:,:) = (gm/(gm-1.d0)*pn(:,:) + 0.5d0*ron(:,:)*v2n(:,:))*vzn(:,:)&
    + (exn(:,:)*byn(:,:)-eyn(:,:)*bxn(:,:))
  r(:,:) = -fx(:,:) / xm(:,:) &
     + ron(:,:)*gxm(:,:)*vxn(:,:) + ron(:,:)*gzm(:,:)*vzn(:,:)
  CALL mlw2d2nd(fx,fz,detot,ix,jx,dt,dx,dz)
  CALL mlw2dsrc2nd(detot,r,ix,jx,dt,dx,dz)

  fx(:,:) = 0.d0
  fz(:,:) = -eyn(:,:)
  r(:,:) = 0.d0
  CALL mlw2d2nd(fx,fz,dbx,ix,jx,dt,dx,dz)
  ! CALL mlw2dsrc2nd()

  fx(:,:) = -ezn(:,:)
  fz(:,:) = exn(:,:)
  r(:,:) = 0.d0
  CALL mlw2d2nd(fx,fz,dby,ix,jx,dt,dx,dz)
  ! CALL mlw2dsrc2nd()

  fx(:,:) = eyn(:,:)
  fz(:,:) = 0.d0
  r(:,:) = - fx(:,:) / xm(:,:)
  CALL mlw2d2nd(fx,fz,dbz,ix,jx,dt,dx,dz)
  CALL mlw2dsrc2nd(dbz,r,ix,jx,dt,dx,dz)


  !!!
  CALL calc_p(ron,pn,vxn,vyn,vzn,bxn,byn,bzn,etotn,gm,ix,jx)
  do j = 1,jx
  do i = 1,ix
    if(pn(i,j) < 0.d0) then
      !.w.rite(*,*) "pressure is negative at (i,j)=", i, j
      !write(*,*) "pressure is negative at (r,z)=", x(i,j), z(i,j)
      !write(*,*) "end of 2nd step" 
      !stop 
      !pn(i,j) = 1.d-3
    endif
  enddo
  enddo
  !!!


  !========================================
  !!  Artificial Viscosity
  !========================================

  qv = 10.d0
  if(ns == 0) then
    CALL put_param_dble("qv:",qv)
  endif

  ! Calculate Artificial Viscosity coefficients (kx,kz)
  CALL calc_coef_k(qv,kx,kz,vx,vz,ix,jx,dx,dz)

  CALL arvis2d(kx,kz,ro,dro,ix,jx,dt,dx,dz)
  CALL arvis2d(kx,kz,rvx,drvx,ix,jx,dt,dx,dz)
  CALL arvis2d(kx,kz,rvy,drvy,ix,jx,dt,dx,dz)
  CALL arvis2d(kx,kz,rvz,drvz,ix,jx,dt,dx,dz)
  CALL arvis2d(kx,kz,etot,detot,ix,jx,dt,dx,dz)
  CALL arvis2d(kx,kz,bx,dbx,ix,jx,dt,dx,dz)
  CALL arvis2d(kx,kz,by,dby,ix,jx,dt,dx,dz)
  CALL arvis2d(kx,kz,bz,dbz,ix,jx,dt,dx,dz)

  !========================================
  !! Update values
  !========================================
  ro(:,:) = ro(:,:) + dro(:,:)
  rvx(:,:) = rvx(:,:) + drvx(:,:)
  rvy(:,:) = rvy(:,:) + drvy(:,:)
  rvz(:,:) = rvz(:,:) + drvz(:,:)
  vx(:,:) = rvx(:,:)/ro(:,:)
  vy(:,:) = rvy(:,:)/ro(:,:)
  vz(:,:) = rvz(:,:)/ro(:,:)
  etot(:,:) = etot(:,:) + detot(:,:)
  bx(:,:) = bx(:,:) + dbx(:,:)
  by(:,:) = by(:,:) + dby(:,:)
  bz(:,:) = bz(:,:) + dbz(:,:)


  !========================================
  !! Apply boundary conditions
  !========================================
  !0220
  CALL bc(ro,0,0,margin,ix,jx)
  CALL bc(ro,1,2,margin,ix,jx)
  CALL bc(ro,2,2,margin,ix,jx)
  CALL bc(ro,3,0,margin,ix,jx)

  !0221
  CALL bc(vx,0,0,margin,ix,jx)
  CALL bc(vx,1,2,margin,ix,jx)
  CALL bc(vx,2,2,margin,ix,jx)
  CALL bc(vx,3,1,margin,ix,jx)

  !0221
  CALL bc(vy,0,0,margin,ix,jx)
  CALL bc(vy,1,2,margin,ix,jx)
  CALL bc(vy,2,2,margin,ix,jx)
  CALL bc(vy,3,1,margin,ix,jx)

  !1220
  CALL bc(vz,0,1,margin,ix,jx)
  CALL bc(vz,1,2,margin,ix,jx)
  CALL bc(vz,2,2,margin,ix,jx)
  CALL bc(vz,3,0,margin,ix,jx)

  !1221
  CALL bc(bx,0,1,margin,ix,jx)
  CALL bc(bx,1,2,margin,ix,jx)
  CALL bc(bx,2,2,margin,ix,jx)
  CALL bc(bx,3,1,margin,ix,jx)

  !1221
  CALL bc(by,0,1,margin,ix,jx)
  CALL bc(by,1,2,margin,ix,jx)
  CALL bc(by,2,2,margin,ix,jx)
  CALL bc(by,3,1,margin,ix,jx)

  !0220
  CALL bc(bz,0,0,margin,ix,jx)
  CALL bc(bz,1,2,margin,ix,jx)
  CALL bc(bz,2,2,margin,ix,jx)
  CALL bc(bz,3,0,margin,ix,jx)

  !0220
  CALL bc(etot,0,0,margin,ix,jx)
  CALL bc(etot,1,2,margin,ix,jx)
  CALL bc(etot,2,2,margin,ix,jx)
  CALL bc(etot,3,0,margin,ix,jx)

  CALL calc_p(ro,p,vx,vy,vz,bx,by,bz,etot,gm,ix,jx)
  CALL calc_e(bx,by,bz,vx,vy,vz,ex,ey,ez)
  
  ay(:,:) = ay(:,:) - ey(:,:) * dt

  !!!!!
 do j = 1,jx
 do i = 1,ix
   if(p(i,j) < 0.d0) then
     !write(*,*) "pressure is negative at (i,j)=", i, j
     !write(*,*) "pressure is negative at (r,z)=", x(i,j), z(i,j)
     !write(*,*) "end of artv step" 
     !write(*,*) "pressure:", p(i,j)j
     !write(*,*) "t=", t, "ns=", ns
     !stop 
     !p(i,j) = 1.d-3
   endif
 enddo
 enddo
  !!!!!


  ns = ns+1
  t = t+dt

  !write(*,*) 'ns:',ns
   !CALL write_0d_dble("dt.dat",dt)

  if (t >= tout) then
  !if(ns > 2) then
  !  stop 1000
  !endif
  !if(ns <= 2) then
    CALL write_0d_dble("t.dat",t)
    CALL write_2d_dble("ro.dat",ro)
    CALL write_2d_dble("vx.dat",vx)
    CALL write_2d_dble("vy.dat",vy)
    CALL write_2d_dble("vz.dat",vz)
    CALL write_2d_dble("bx.dat",bx)
    CALL write_2d_dble("by.dat",by)
    CALL write_2d_dble("bz.dat",bz)
    CALL write_2d_dble("p.dat",p)
    CALL write_2d_dble('etot.dat',etot)
    CALL write_2d_dble("ay.dat",ay)


    write(*,foutmsg) ns,t,nd
    tout = tout + dtout
    nd = nd + 1
  endif
  end do main_loop

  write(*,fstopmsg) ns,t
  write(*,*) 'Normal Termination'

END PROGRAM main
