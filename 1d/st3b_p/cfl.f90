module cfl
  implicit none
contains
  subroutine calc_dt(dt,ix,dx,ro,vx,vy,bx,by,pr,gm)
    integer,intent(in) :: ix
    real(8),intent(inout) :: dt
    real(8),intent(in) :: dx,ro(ix),vx(ix),vy(ix),bx(ix),by(ix),&
      & pr(ix),gm

    real(8) :: safety
    real(8) :: dt_min,dt_tmp
    real(8) :: vv,bb
    real(8),parameter :: pi = 4.d0 * atan(1.d0)
    integer :: ii

    dt_min = 1.d-3
    dt = dt_min
    safety = 0.4d0
    !dt = safety * dx * 1.d-3
    do ii = 1,ix
    vv = vx(ii)**2 + vy(ii)**2
    bb = bx(ii)**2 + by(ii)**2
    dt_tmp = dx / sqrt(vv + bb/(4.d0*pi*ro(ii)) + gm*pr(ii)/ro(ii))
    dt = min(dt_tmp,dt_min,dt)
    enddo

    dt = dt * safety

    if (dt < 1.d-10) then
      write(*,*) 'dt is too small < 1e-10'
      stop
    endif

  end subroutine calc_dt
end module cfl

