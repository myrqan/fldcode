module mlw
  implicit none
  contains
    subroutine mlw1d1st(uu,uum,ff,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: uu(ix),ff(ix)
      double precision,intent(inout) :: uum(ix)
      double precision,intent(in) :: dt,dx
      integer :: i
      uum(:) = 0.d0
      do i = 1, ix-1
      uum(i) = 0.5d0 * (uu(i+1)+uu(i)) - 0.5d0*dt/dx * (ff(i+1)-ff(i))
      enddo
    end subroutine mlw1d1st

    subroutine mlw1d2nd(uu,uun,ff,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: uu(ix),ff(ix)
      double precision,intent(inout) :: uun(ix)
      double precision,intent(in) :: dt,dx
      integer :: i
      uun = 0.d0
      do i = 2, ix
      uun(i) = uu(i) - dt/dx * (ff(i) - ff(i-1))
      enddo
    end subroutine mlw1d2nd 

    subroutine mlw1dsrc(uu,ss,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: ss(ix)
      double precision,intent(in) :: dt,dx
      double precision,intent(inout) :: uu(ix)
      double precision :: uun(ix)
      integer :: i
      uun = 0.d0
      do i = 1, ix
      uun(i) = uu(i) + ss(i) + dt
      enddo
      uu = uun
    end subroutine mlw1dsrc
end module mlw
