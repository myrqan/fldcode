program main
  use datw
  implicit none

  integer,parameter :: margin = 0
  integer,parameter :: gx = 1001
  real(8),parameter :: xmin = 0.d0
  real(8),parameter :: xmax = 1.d0
  real(8),parameter :: pi = 4.d0 * atan(1.d0)
  real(8) :: xx(gx),dx
  real(8) :: uu(gx),un(gx)

  real(8) :: ns
  real(8) :: time, t_end, dt, t_out, dt_out
  real(8) :: teps= 1.d-3
  real(8),parameter :: safety = 0.4d0
  real(8) :: rr

  integer :: ii

  !================
  ! make grid
  ! initial condition
  !================
  call dataclean()

  dx = (xmax-xmin) / (gx-1)
  xx(1) = 0.d0
  do ii = 2, gx
  xx(ii) = xx(ii-1) + dx
  end do

  do ii = 1, gx
  uu(ii) = dsin(pi*xx(ii))
  end do


  ns = 0
  time = 0.d0
  t_end = 0.5d0 + teps
  dt_out = 0.05d0
  t_out = dt_out



  call w0d(time,'dat/t.dat')
  call w1d(gx,xx,'dat/x.dat')
  call w1d(gx,uu,'dat/uu.dat')
  write(*,*) '[write], time = ', time

  !=======================

  do while (time < t_end)

  dt = safety * (0.5d0 * dx**2)
  rr = dt / dx**2
  
  do ii = 2,gx-1
  un(ii) = rr * uu(ii-1) + (1-2*rr) * uu(ii) + rr * uu(ii+1)
  enddo
   
  uu(:) = un(:)

  ! boudnary condition
  uu(1) = 0.d0
  uu(gx) = 0.d0


  if(time >= t_out) then
    call w0d(time,'dat/t.dat')
    call w1d(gx,uu,'dat/uu.dat')
    t_out = t_out+ dt_out

    write(*,*) '[write], time = ', time
  end if

  time = time + dt
  end do


  stop

end program main
