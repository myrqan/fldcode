module mlw
  implicit none
contains

  subroutine mlw2d1st(u,fx,fz,un,du,ix,jx,dt,dx,dz)
    integer,intent(in)::ix,jx
    real(8),intent(in)::dt,dx,dz
    real(8),intent(in)::u(ix,jx),fx(ix,jx),fz(ix,jx)
    real(8),intent(out)::un(ix,jx),du(ix,jx)
    integer::i,j
    do j=1,jx-1
    do i=1,ix-1
    if(i>=2 .and. j>=2) then
      du(i,j) = du(i,j) - 0.5d0*dt * 0.5d0/dx * (fx(i+1,j)-fx(i-1,j))&
        - 0.5d0*dt * 0.5d0/dz * (fz(i,j+1)-fz(i,j-1))
    endif

    un(i,j) = 0.25d0 * (u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))&
      - dt/dx * 0.5d0 * (fx(i+1,j)-fx(i,j)+fx(i+1,j+1)-fx(i,j+1))&
      - dt/dz * 0.5d0 * (fz(i,j+1)-fz(i,j)+fz(i+1,j+1)-fz(i+1,j))
    enddo
    enddo
  end subroutine mlw2d1st

  subroutine mlw2dsrc1st(r,un,du,ix,jx,dt,dx,dz)
    integer,intent(in)::ix,jx
    real(8),intent(in)::r(ix,jx),dt,dx,dz
    real(8),intent(inout)::un(ix,jx),du(ix,jx)
    integer::i,j
    do j=1,jx-1
    do i=1,ix-1
    if(i>=2 .and. j>=2) then 
      du(i,j) = du(i,j) + 0.5d0*dt * r(i,j)
    endif
    un(i,j) = un(i,j) + 0.25d0*dt * (r(i+1,j)+r(i,j)+r(i+1,j+1)+r(i,j+1))
    enddo
    enddo
  end subroutine mlw2dsrc1st

  subroutine mlw2d2nd(fx,fz,du,ix,jx,dt,dx,dz)
    integer,intent(in)::ix,jx
    real(8),intent(in)::fx(ix,jx),fz(ix,jx),dt,dx,dz
    real(8),intent(inout)::du(ix,jx)
    integer::i,j
    do j=2,jx-1
    do i=2,ix-1
    du(i,j) = du(i,j) &
      - 0.25d0*dt/dx*(fx(i,j-1)-fx(i-1,j-1)+fx(i,j)-fx(i-1,j))&
      - 0.25d0*dt/dz*(fz(i-1,j)-fz(i-1,j-1)+fz(i,j)-fz(i,j-1))
    enddo
    enddo
  end subroutine mlw2d2nd

  subroutine mlw2dsrc2nd(du,r,ix,jx,dt,dx,dz)
    integer,intent(in)::ix,jx
    real(8),intent(in):: dt,dx,dz,r(ix,jx)
    real(8),intent(inout)::du(ix,jx)
    integer::i,j
    do j=2,jx-1
    do i=2,ix-1
    du(i,j) = du(i,j) + 0.5d0 * dt &
      * 0.25d0 *(r(i-1,j-1)+r(i,j-1)+r(i-1,j)+r(i,j))
    enddo
    enddo
  end subroutine mlw2dsrc2nd
end module mlw
