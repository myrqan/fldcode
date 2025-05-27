module mlw
  implicit none
contains
  subroutine mlw1d1st(u,f,um,du,ix,dt,dx)
    INTEGER,INTENT(IN)::ix
    DOUBLE PRECISION,INTENT(IN)::u(ix),f(ix)
    DOUBLE PRECISION,INTENT(IN)::dt,dx
    DOUBLE PRECISION,INTENT(INOUT)::um(ix),du(ix)
    INTEGER::i
    do i = 1,ix-1
    if(i>=2) then
      du(i) = du(i) - 0.5d0*dt/dx * 0.5d0*(f(i+1)-f(i-1))
    endif
    um(i) = 0.5d0*(u(i+1)+u(i))-dt/dx*(f(i+1)-f(i))
    enddo
  end subroutine mlw1d1st
  subroutine mlw1dsrc1st(um,du,s,ix,dt)
    INTEGER,INTENT(IN)::ix
    DOUBLE PRECISION,INTENT(IN)::s(ix),dt
    DOUBLE PRECISION,INTENT(INOUT)::um(ix),du(ix)
    INTEGER::i
    do i = 1, ix-1
    du(i) = du(i) + 0.5d0 * dt * s(i)
    um(i) = um(i) + 0.5d0 * dt * (s(i+1)+s(i))
    enddo
  end subroutine mlw1dsrc1st
  subroutine mlw1d2nd(fm,du,ix,dt,dx)
    INTEGER,INTENT(IN)::ix
    DOUBLE PRECISION,INTENT(IN)::fm(ix),dt,dx
    DOUBLE PRECISION,INTENT(INOUT)::du(ix)
    INTEGER::i
    do i = 2, ix-1
    du(i) = du(i) - 0.5d0*dt/dx*(fm(i)-fm(i-1))
    enddo
  end subroutine mlw1d2nd
  subroutine mlw1dsrc2nd(du,sm,ix,dt)
    INTEGER,INTENT(IN)::ix
    DOUBLE PRECISION,INTENT(IN)::sm(ix),dt
    DOUBLE PRECISION,INTENT(INOUT)::du(ix)
    INTEGER::i
    do i = 2, ix-1
    du(i) = du(i) + 0.5d0*dt* 0.5d0*(sm(i-1)+sm(i))
    enddo
  end subroutine mlw1dsrc2nd

 ! subroutine mlw1dsrc1st(u,f,um,du,s,ix,dt,dx)
 !   integer,intent(in) :: ix
 !   doubleprecision,intent(in) :: u(ix),f(ix),s(ix)
 !   doubleprecision,intent(in) :: dt,dx
 !   doubleprecision,intent(inout) :: um(ix),du(ix)
 !   !!
 !   !! um : t+dt,n+1/2; umm : t+dt, n
 !   !! du : contribution to next dt
 !   !!
 !   integer :: i
 !   um = 0.d0; umm = 0.d0; du = 0.d0

 !   do i = 1, ix-1
 !   if(i >= 2) then
 !     du(i) = du(i) - 0.5d0*dt/dx * 0.5d0*(f(i+1)-f(i-1))! + 0.5d0*dt*s(i)
 !     umm(i) = 0.5d0*(u(i+1)+u(i-1)) - 0.5d0*dt/dx*(f(i+1)-f(i-1)) + dt*s(i)
 !   endif
 !   um(i) = 0.5d0*(u(i+1)+u(i)) - dt/dx*(f(i+1)-f(i)) + 0.5d0*dt*(s(i+1)+s(i))
 !   enddo
 ! end subroutine mlw1dsrc1st

 ! subroutine mlw1dsrc2nd(fm,du,smm,ix,dt,dx)
 !   integer,intent(in) :: ix
 !   double precision,intent(in) :: fm(ix),smm(ix)
 !   double precision,intent(in) :: dt,dx
 !   double precision,intent(inout) :: du(ix)
 !   integer :: i

 !   do i = 2, ix-1
 !   du(i) = du(i) - 0.5d0*dt/dx*(fm(i)-fm(i-1)) + 0.5d0*dt*smm(i)
 !   enddo
 ! end subroutine mlw1dsrc2nd
end module mlw
