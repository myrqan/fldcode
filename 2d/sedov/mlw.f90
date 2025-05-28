module mlw
  implicit none
contains
  subroutine mlw1d1st(u,f,um,du,ix,dt,dx)
    INTEGER,INTENT(IN)::ix
    DOUBLE PRECISION,INTENT(IN)::u(ix),f(ix)
    DOUBLE PRECISION,INTENT(IN)::dt,dx
    DOUBLE PRECISION,INTENT(INOUT)::um(ix),du(ix)
    INTEGER::j
    do j = 1,ix-1
    if(j>=2) then
      du(j) = du(j) - 0.25d0*dt/dx * (f(j+1)-f(j-1))
    endif
    um(j) = 0.5d0*(u(j+1)+u(j))-dt/dx*(f(j+1)-f(j))
    enddo
  end subroutine mlw1d1st
  subroutine mlw1dsrc1st(um,du,s,ix,dt)
    INTEGER,INTENT(IN)::ix
    DOUBLE PRECISION,INTENT(IN)::s(ix),dt
    DOUBLE PRECISION,INTENT(INOUT)::um(ix),du(ix)
    INTEGER::j
    do j = 1, ix-1
    du(j) = du(j) + 0.5d0 * dt * s(j)
    um(j) = um(j) + 0.5d0 * dt * (s(j+1)+s(j))
    enddo
  end subroutine mlw1dsrc1st
  subroutine mlw1d2nd(fm,du,ix,dt,dx)
    INTEGER,INTENT(IN)::ix
    DOUBLE PRECISION,INTENT(IN)::fm(ix),dt,dx
    DOUBLE PRECISION,INTENT(INOUT)::du(ix)
    INTEGER::j
    do j = 2, ix-1
    du(j) = du(j) - 0.5d0*dt/dx*(fm(j)-fm(j-1))
    enddo
  end subroutine mlw1d2nd
  subroutine mlw1dsrc2nd(du,sm,ix,dt)
    INTEGER,INTENT(IN)::ix
    DOUBLE PRECISION,INTENT(IN)::sm(ix),dt
    DOUBLE PRECISION,INTENT(INOUT)::du(ix)
    INTEGER::j
    do j = 2, ix-1
    du(j) = du(j) + 0.5d0*dt*sm(j)
    enddo
  end subroutine mlw1dsrc2nd

  subroutine mlw2d1st(u,fx,fz,un,du,ix,jx,dt,dx,dz)
    INTEGER,INTENT(IN)::ix,jx
    DOUBLE PRECISION,INTENT(IN)::dt,dx,dz
    DOUBLE PRECISION,INTENT(IN)::u(ix,jx),fx(ix,jx),fz(ix,jx)
    DOUBLE PRECISION,INTENT(OUT)::un(ix,jx),du(ix,jx)
    INTEGER::i,j
    do j=1,ix-1
    do i=1,ix-1
    if(i>=2 .and. j>=2) then
      du(i,j) = du(i,j) - 0.5d0*dt * 0.5d0/dx * (fx(i+1,j)-fx(i-1,j))&
                        - 0.5d0*dt * 0.5d0/dx * (fz(i,j+1)-fz(i,j-1))
    endif

      un(i,j) = 0.25d0 * (u(i,j)+u(i+1,j)+u(i,j)+u(i+1,j+1))&
                - dt/dx * 0.5d0 * (fx(i+1,j)-fx(i,j)+fx(i+1,j+1)-fx(i,j+1))&
                - dt/dz + 0.5d0 * (fz(i,j+1)-fz(i,j)+fz(i+1,j+1)-fz(i,j+1))
    enddo
    enddo
  end subroutine mlw2d1st

  subroutine mlw2dsrc1st(r,un,du,ix,jx,dt,dx,dz)
    INTEGER,INTENT(in)::ix,jx
    DOUBLE PRECISION,INTENT(IN)::r(ix,jx),dt,dx,dz
    DOUBLE PRECISION,INTENT(INOUT)::un(ix,jx),du(ix,jx)
    INTEGER::i,j
    do j = 1,jx-1
    do i = 1,ix-1
      du(i,j) = du(i,j) + 0.5d0*dt * r(i,j)
      un(i,j) = un(i,j) + 0.25d0*dt*(r(i+1,j)+r(i,j)+r(i+1,j+1)+r(i,j+1))
    enddo
    enddo
  end subroutine mlw2dsrc1st

  subroutine mlw2d2nd(fx,fz,du,ix,jx,dt,dx,dz)
    INTEGER,INTENT(IN)::ix,jx
    DOUBLE PRECISION,INTENT(IN)::fx(ix,jx),fz(ix,jx),dt,dx,dz
    DOUBLE PRECISION,INTENT(INOUT)::du(ix,jx)
    INTEGER::i,j
    do j=2,jx-1
    do i=2,ix-1
      du(i,j) = du(i,j) - 0.5d0*dt/dx*(fx(i,j-1)-fx(i-1,j-1))
    enddo
    enddo
  end subroutine mlw2d2nd

  subroutine mlw2dsrc2nd(du,r,ix,jx,dt,dx,dz)
    INTEGER,INTENT(IN)::ix,jx
    DOUBLE PRECISION,INTENT(IN):: dt,dx,dz,r(ix,jx)
    DOUBLE PRECISION,INTENT(INOUT)::du(ix,jx)
    INTEGER::i,j

    do j=2,ix
    do i=2,ix
      du(i,j) = du(i,j) + 0.5d0 * dt &
                  * 0.25d0 *(r(i-1,j-1)+r(i,j-1)+r(i-1,j)+r(i,j))
    enddo
    enddo
  end subroutine mlw2dsrc2nd
end module mlw
