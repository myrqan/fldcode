module mlw
  implicit none
  contains
    subroutine mlw1d1st(uu,uum,ff,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: uu(0:ix),ff(0:ix)
      double precision,intent(inout) :: uum(0:ix)
      double precision,intent(in) :: dt,dx
      integer :: i
      uum(:) = 0.d0
      do i = 0, ix-1
      uum(i) = 0.5d0 * (uu(i+1)+uu(i)) - 0.5d0*dt/dx * (ff(i+1)-ff(i))
      enddo
    end subroutine mlw1d1st

    subroutine mlw1d2nd(uu,uun,ff,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: uu(0:ix),ff(0:ix)
      double precision,intent(inout) :: uun(0:ix)
      double precision,intent(in) :: dt,dx
      integer :: i
      uun = 0.d0
      do i = 1, ix-1
      uun(i) = uu(i) - dt/dx * (ff(i) - ff(i-1))
      enddo
    end subroutine mlw1d2nd 

    subroutine mlw1dsrc(uu,uun,ss,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: uu(0:ix),ss(0:ix)
      double precision,intent(in) :: dt,dx
      double precision,intent(inout) :: uun(0:ix)
      integer :: i
      uun = 0.d0
      do i = 1, ix-1
      uun(i) = uu(i) + ss(i) + dt
      enddo
    end subroutine mlw1dsrc
end module mlw
