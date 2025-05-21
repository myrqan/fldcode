module arvis
  implicit none
  contains
    subroutine calc_coef_k(qv,k,vvx,ix,dx)
      double precision,intent(in) :: qv,vvx(ix),dx
      double precision,intent(inout) :: k(ix)
      integer,intent(in) :: ix
      integer :: i
      do i = 1, ix-1
      k(i) = qv*dx*abs(vvx(i+1)-vvx(i))
      enddo
    end subroutine calc_coef_k

    subroutine arvis1d(k,uu,uum,ix,dt,dx)
      integer,intent(in) :: ix
      double precision,intent(in) :: k(ix),uu(ix),dt,dx
      double precision,intent(inout) :: uum(ix)
      integer :: i
      uum = 0.d0
      do i = 2, ix-1
      uum(i) = uu(i) + dt/dx * &
        & (k(i)/dx * (uu(i+1)-uu(i)) &
        & -k(i-1)/dx * (uu(i)-uu(i-1))) 
      enddo
    end subroutine arvis1d
end module arvis
