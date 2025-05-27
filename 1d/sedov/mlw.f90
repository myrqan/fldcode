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
end module mlw
