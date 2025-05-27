module cfl
  implicit none
  contains
    subroutine calc_dt(rho,vx,p,gm,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: rho(ix),vx(ix),p(ix)
      double precision,intent(in) :: gm,dx
      double precision,intent(inout) :: dt
      double precision :: dt_cfl(ix), safety = 0.4d0
      INTEGER :: i
      dt_cfl = 1e20
      do i = 2, ix-1
      dt_cfl(i) = dx / (abs(vx(i)) + sqrt(gm*p(i)/rho(i)))
      enddo
      dt = safety * MINVAL(dt_cfl)
    end subroutine calc_dt
end module cfl
