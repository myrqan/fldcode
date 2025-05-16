module arvis
  implicit none
  contains
    subroutine calc_coef_k(qv,k,vvx,ix,dx)
      double precision,intent(in) :: qv,vvx(0:ix),dx
      double precision,intent(inout) :: k(0:ix)
      integer,intent(in) :: ix
      integer :: i
      do i = 0, ix-1
      k(i) = qv*dx*abs(vvx(i+1)-vvx(i))
      enddo
    end subroutine calc_coef_k

    subroutine arvis1d(qv,k,uu,uum,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: qv,k(0:ix),uu(0:ix),dt,dx
      double precision,intent(inout) :: uum(0:ix)
      integer :: i
      uum = 0.d0
      do i = 1, ix-1
      uum(i) = uu(i) + dt/dx * &
        & (k(i)/dx * (uu(i+1)-uu(i)) &
        & -k(i-1)/dx * (uu(i) - uu(i-1))) 
      enddo
    end subroutine arvis1d
end module arvis
