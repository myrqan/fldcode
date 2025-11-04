module mlw
  implicit none
contains
  subroutine mlw1d1st(qm,uu,fx,ix,dt,dx)
    integer,intent(in) :: ix
    real(8),intent(inout) :: qm(ix)
    real(8),intent(in) :: uu(ix),fx(ix),dt,dx
    real(8) :: dxi
    integer :: ii

    dxi = 1.d0/dx

    do ii = 1, ix-1
    qm(ii) = 0.5d0 * (uu(ii) + uu(ii+1)) &
      & - dt*dxi*(fx(ii+1) - fx(ii))
    enddo
  end subroutine mlw1d1st

  subroutine mlw1d2nd(dq,fx,fxm,ix,dt,dx)
    integer,intent(in) :: ix
    real(8),intent(inout) :: dq(ix)
    real(8),intent(in) :: fx(ix),fxm(ix),dt,dx
    real(8) :: dxi
    integer :: ii

    dxi = 1.d0/dx
    do ii = 2,ix-1
    dq(ii) = - 0.25d0 * dt * dxi * (fx(ii+1) - fx(ii-1) &
      & +2.d0*fxm(ii) - 2.d0*fxm(ii-1))
    enddo
  end subroutine mlw1d2nd
end module mlw
