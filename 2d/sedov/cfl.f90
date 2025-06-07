module cfl
  implicit none
  contains
    subroutine calc_dt(rho,vx,vz,p,gm,ix,jx,dt,dx,dz)
      integer,intent(in) :: ix,jx
      double precision,intent(in) :: rho(ix,jx),vx(ix,jx),vz(ix,jx),p(ix,jx)
      double precision,intent(in) :: gm,dx,dz
      double precision,intent(inout) :: dt
      double precision :: dt_cfl,safety = 0.4d0
      DOUBLE PRECISION :: min_val
      INTEGER :: i,j
      INTEGER :: iout, jout
      DOUBLE PRECISION :: dtmin = 1.d-10
      min_val = 1e20
      do j = 2, jx-1
      do i = 2, ix-1
        dt_cfl = min(dx,dz) / sqrt(vx(i,j)**2+vz(i,j)**2 + gm*p(i,j)/rho(i,j))
        min_val = MIN(min_val,dt_cfl)
        if(min_val == dt_cfl) then
          iout = i; jout = j
        endif
      enddo
      enddo
      dt = safety * min_val
      if (dt < dtmin) then
        stop 999
      endif
      write(*, *) "dt:", dt
      write(*, *) "iout, jout", iout, jout
    end subroutine calc_dt
end module cfl

