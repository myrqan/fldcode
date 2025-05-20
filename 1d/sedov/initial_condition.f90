module initial_condition
  implicit none
  contains
!=======================================================
!    INITIALIZE w/ zeros
!=======================================================
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

!=======================================================
!    MAKE GRIDS
!=======================================================
    subroutine gridx(ix,margin,dx,x)
      integer, intent(in) :: ix, margin
      double precision, intent(inout):: dx
      double precision,intent(inout) :: x(ix)
      integer :: ixx
      ! ix : number of grids
      ! dx : grid size (double)
      ! x(0:ix) : position (double)
      integer :: i
      ixx = ix - margin
      
      dx = 1.d0 / dble(ixx)
      x(1) = dx
      do i = 2, ix
        x(i) = x(i-1) + dx
      enddo
    end subroutine gridx

!=======================================================
!    EXAMPLE OF INITIAL CONDITIONS
!=======================================================
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

    subroutine init_sedov1d(ix,x,rho,p)
      integer,intent(in) :: ix
      double precision,intent(inout) :: x(ix),rho(ix),p(ix)
      integer :: i
      double precision :: we

      rho = 1.d0;
      we = 0.02d0;
      do i = 1, ix
        p(i) = (1.d-8 + (1-1.d-8)*exp(-x(i)**2/we**2))
      enddo
    end subroutine init_sedov1d
end module initial_condition
