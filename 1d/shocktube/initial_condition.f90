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

    subroutine gridx(ix,margin,dx,x)
      integer, intent(in) :: ix, margin
      double precision, intent(inout):: dx
      double precision,intent(inout) :: x(ix)
      ! ix : number of grids
      ! dx : grid size (double)
      ! x(0:ix) : position (double)
      integer :: i, iix
      iix = ix - 2*margin
      dx = 1.d0 / dble(iix)
      x(1) = 0.d0
      do i = 2, ix
        x(i) = x(i-1) + dx
      enddo
    end subroutine gridx

    subroutine init_shocktube(ix,x,rho,p)
      integer,intent(in) :: ix
      double precision,intent(inout) :: x(ix),rho(ix),p(ix)
      ! ix : number of grids
      ! x, rho, p : array (double)
      integer :: i

      do i = 1, ix
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
