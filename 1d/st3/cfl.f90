module cfl
  implicit none
contains
  subroutine calc_dt(dt,ix,dx,ro,vx,pr,gm)
    integer,intent(in) :: ix
    real(8),intent(inout) :: dt
    real(8),intent(in) :: dx,ro(ix),vx(ix),pr(ix),gm

    real(8) :: safety
    real(8) :: dt_min,dt_tmp
    integer :: ii

    dt_min = 1.d-3
    safety = 0.4d0
    !dt = safety * dx * 1.d-3
    do ii = 1,ix
    dt = dx / sqrt(vx(ii)**2 + gm*pr(ii)/ro(ii))
    dt = min(dt,dt_min)
    enddo

    dt = dt * safety

    if (dt < 1.d-10) then
      write(*,*) 'dt is too small < 1e-10'
      stop
    endif

  end subroutine calc_dt
end module cfl

