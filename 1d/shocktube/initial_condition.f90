module initial_condition
  implicit none
  contains
    subroutine init0d(sc)
      double precision, intent(inout) :: sc
      ! sc : scalar value (double)
      sc = 0.d0
    end subroutine init0d

    subroutine init1d(ary)
      double precision,intent(inout) :: ary(:)
      ! ary : list value (double)
      ary = 0.d0
    end subroutine init1d

    subroutine gridx(ix,dx,x)
      integer, intent(in) :: ix
      double precision, intent(inout):: dx
      double precision,intent(inout) :: x(0:ix)
      ! ix : number of grids
      ! dx : grid size (double)
      ! x(0:ix) : position (double)
      integer :: i
      
      dx = 1.d0 / dble(ix)

      x(0) = 0.d0
      do i = 1, ix
        x(i) = x(i-1) + dx
      enddo
    end subroutine gridx

    subroutine init_shocktube(ix,x,rho,p)
      integer,intent(in) :: ix
      double precision,intent(inout) :: x(0:ix),rho(0:ix),p(0:ix)
      ! ix : number of grids
      ! x, rho, p : array (double)
      integer :: i

      do i = 0, ix
      if (x(i) <= 0.5d0) then
        rho(i)=1.d0
        p(i)=1.d0
      else
        rho(i)=0.125d0
        p(i)=0.1d0
      endif
      enddo
    end subroutine init_shocktube
end module initial_condition
