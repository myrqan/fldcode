module cfl
  implicit none
  contains
    subroutine calc_dt(rho,vx,vy,vz,p,bx,by,bz,gm,ix,jx,dt,dx,dz,margin)
      integer,intent(in) :: ix,jx,margin
      double precision,intent(in) :: rho(ix,jx),vx(ix,jx),vy(ix,jx),&
        vz(ix,jx),p(ix,jx),bx(ix,jx),by(ix,jx),bz(ix,jx)
      double precision,intent(in) :: gm,dx,dz
      double precision,intent(inout) :: dt
      double precision :: dt_cfl,safety = 0.4d0
      DOUBLE PRECISION :: min_val,v2,va2,vs2
      INTEGER :: i,j
      INTEGER :: iout, jout
      DOUBLE PRECISION :: dtmin = 1.d-20

      min_val = 1e20
      do j = margin, jx-margin
      do i = margin, ix-margin
        v2 = vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2
        va2 = (bx(i,j)**2 + by(i,j)**2 + bz(i,j)**2) / rho(i,j)
        vs2 = gm * p(i,j) / rho(i,j)
        dt_cfl = min(dx,dz) / sqrt(v2+va2+vs2)
        min_val = MIN(min_val,dt_cfl)
        if(min_val == dt_cfl) then
          iout = i; jout = j
        endif
      enddo
      enddo
      dt = safety * min_val
      if (dt < dtmin) then
        write(*,*) 'dt=', dt
        stop 349
      endif
      !write(*, *) "dt:", dt
      !write(*, *) "iout, jout", iout, jout
    end subroutine calc_dt
end module cfl

