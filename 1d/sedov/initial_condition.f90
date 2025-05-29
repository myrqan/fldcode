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

    subroutine gridx(ix,margin,dx,x,xm)
      integer, intent(in) :: ix, margin
      double precision,intent(inout):: dx
      double precision,intent(inout) :: x(ix),xm(ix)
      ! ix : number of grids
      ! dx : grid size (double)
      ! x(0:ix) : position (double)
      integer :: i, iix
      double precision :: xmin
      iix = ix - 2*margin
      dx = 1.d0 / dble(iix)
      !xmin = 0.02d0
      xmin = dx
      !x(1) = 2*dx
      x(1) = xmin - 0.5d0 * dx
      xm(1) = xmin
      do i = 2, ix
        x(i) = x(i-1) + dx
        xm(i) = x(i) + 0.5d0 * dx
      enddo
      !xm(ix) = xm(ix-1)
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

    subroutine init_sedov1d(ix,x,rho,vx,p)
      integer,intent(in) :: ix
      double precision,intent(in) :: x(ix)
      double precision,intent(inout) :: rho(ix),vx(ix),p(ix)
      integer :: i
      double precision :: we
      rho = 1.d0;
      we = 0.1d0;
      vx = 0.d0;
      do i = 1, ix
        p(i) = (1.d-8 + (1-1.d-8)*exp(-x(i)**2/(we**2)))
      enddo
    end subroutine init_sedov1d

    subroutine cross_section(ix,s,sm,ds,dsm,x,xm,dx)
      INTEGER,INTENT(IN)::ix
      DOUBLE PRECISION,INTENT(INOUT)::s(ix),sm(ix),ds(ix),dsm(ix)
      DOUBLE PRECISION,INTENT(IN)::x(ix),xm(ix),dx
      INTEGER::j
      do j = 1, ix
      s(j) = x(j)**2
      if(j<ix)  then
        sm(j) = xm(j)**2
      endif

      if(j > 1) then
        ds(j) = (sm(j)-sm(j-1))/dx
      endif
      if(j < ix) then
        dsm(j) = (s(j+1)-s(j)) / dx
      endif
      enddo
      ds(1) = 0.d0
      dsm(1) = 0.d0

    end subroutine cross_section
end module initial_condition
