module arvis
  implicit none
contains
  subroutine arvis1d(qq,dqq,ix,dt,dx,av,vx)
    integer,intent(in) :: ix
    real(8),intent(inout) :: qq(ix),dqq(ix)
    real(8),intent(in) :: vx(ix)
    real(8),intent(in) :: dt,dx,av
    real(8),parameter :: dvmin = 1.d-2
    real(8) :: kx(ix),dxi,kh(ix)
    integer :: ii

    dxi = 1.d0 / dx

    do ii = 2,ix-1
    kx(ii) = av * dx * (max(0.5d0*abs(vx(ii+1) - vx(ii-1)),dvmin)&
      & - dvmin)
    enddo

    do ii = 1, ix-1
    kh(ii) = 0.5d0 * (kx(ii) + kx(ii+1))
    enddo

    do ii = 2,ix-1
    dqq(ii) = dt/dx * (kh(ii) * (qq(ii+1) - qq(ii)) * dxi &
      & - kh(ii-1) * (qq(ii) - qq(ii-1)) * dxi)
    enddo

  end subroutine arvis1d
end module arvis
