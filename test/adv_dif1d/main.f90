program main
  use datw
  use basic
  use time_step
  use update
  use sts
  implicit none

  integer,parameter :: margin = 0
  integer,parameter :: gx = 1001
  real(8),parameter :: xmin = 0.d0
  real(8),parameter :: xmax = 1.d0
  real(8),parameter :: pi = 4.d0 * atan(1.d0)
  real(8),parameter :: eta = 1.d-1
  real(8) :: xx(gx),dx
  real(8) :: uu(gx),un(gx)

  integer(8) :: ns
  real(8) :: time, t_end, dt, t_out, dt_out
  real(8) :: teps= 1.d-3
  real(8) :: rr

  integer :: sts_s

  integer :: ii

  !================
  ! make grid
  ! initial condition
  !================
  call dataclean()

  call makegrid(xx,dx,xmax,xmin,gx)

  call initial_condition(uu,xx,gx)



  ns = 0
  time = 0.d0
  t_end = teps
  t_end = t_end + 0.5d0
  dt_out = 0.05d0
  t_out = dt_out



  call w0d(time,'dat/t.dat')
  call w1d(gx,xx,'dat/x.dat')
  call w1d(gx,uu,'dat/uu.dat')
  write(*,*) '[write], time = ', time, 'steps=', ns

  !=======================

  do while (time < t_end)

  call determine_sts_s(sts_s,dx,eta)


  !!!!!!!!!!!!!!!!!!!!!
  !! Strang operator splitting methodを実装
  !!!!!!!!!!!!!!!!!!!!!

  call cfl_ad(dt,dx,gx)

  call advflow(uu,dx,dt,gx)



  !call cfl_tc(dt,dx,gx,eta)
  !rr = eta * dt / dx**2
  !
  !do ii = 2,gx-1
  !un(ii) = rr * uu(ii-1) + (1-2*rr) * uu(ii) + rr * uu(ii+1)
  !enddo
  ! 
  !uu(:) = un(:)

  ! boudnary condition
  !uu(1) = 0.d0
  !uu(gx) = 0.d0


  time = time + dt
  ns = ns+1


  if(time >= t_out) then
    call w0d(time,'dat/t.dat')
    call w1d(gx,uu,'dat/uu.dat')
    t_out = t_out+ dt_out

    write(*,*) '[write], time = ', time, 'steps=', ns
  end if

  end do


  stop

end program main
