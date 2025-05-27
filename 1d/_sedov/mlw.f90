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
      uum(i) = 0.5d0 * (uu(i+1)+uu(i)) - 0.5d0 * dt/dx * (ff(i+1)-ff(i))
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

    subroutine mlw1dsrc1st(uu,ss,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: ss(ix)
      double precision,intent(in) :: dt,dx
      double precision,intent(inout) :: uu(ix)
      double precision :: uun(ix)
      integer :: i
      uun = 0.d0
      do i = 1, ix-1
      uun(i) = uu(i) + 0.25d0 * dt * (ss(i)+ss(i+1))
      enddo
      uu = uun
    end subroutine mlw1dsrc1st

    subroutine mlw1dsrc2nd(uu,ssm,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: ssm(ix)
      double precision,intent(in) :: dt,dx
      double precision,intent(inout) :: uu(ix)
      double precision :: uun(ix)
      integer :: i
      uun = 0.d0
      do i = 2, ix
      uun(i) = uu(i) +  dt * ssm(i)
      enddo
      uu = uun
    end subroutine mlw1dsrc2nd


    subroutine mlw1d1st_wsrc(u,um,f,r,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: u(ix),f(ix),r(ix)
      double precision,intent(in) :: dt,dx
      double precision,intent(inout) :: um(ix)
      integer :: i
      um = 0.d0
      do i = 1, ix-1
      um(i) = 0.5d0 * (u(i+1) + u(i)) - 0.5d0 * dt/dx * (f(i+1)- f(i)) &
        & + 0.5d0 * dt * 0.5d0 * (r(i) + r(i+1))
      enddo
    end subroutine mlw1d1st_wsrc

    subroutine mlw1d2nd_wsrc(u,un,fm,rm,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: u(ix),fm(ix),rm(ix)
      double precision,intent(in) :: dt,dx
      double precision,intent(inout) :: un(ix)
      integer :: i
      un = 0.d0
      do i = 2, ix
      un(i) = u(i) - dt/dx * (fm(i) - fm(i-1)) &
        & + dt * rm(i)
        !& + dt * 0.5d0 * (rm(i) + rm(i-1))
      enddo
    end subroutine mlw1d2nd_wsrc

end module mlw
