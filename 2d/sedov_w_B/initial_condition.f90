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
    
    subroutine init2d(ary)
      double precision,intent(inout) :: ary(:,:)
      ! ary : list value (double)
      ary(:,:) = 0.d0
    end subroutine init2d

    subroutine grid2d(ix,jx,margin,dx,dz,x,z,xm,zm)
      integer, intent(in) :: ix,jx,margin
      double precision,intent(inout):: dx,dz
      double precision,intent(inout) :: x(ix,jx),xm(ix,jx),z(ix,jx),zm(ix,jx)
      ! ix : number of grids
      ! dx,dz : grid size (double)
      ! x(ix), z(jx) : position (double)
      integer :: i,j, iix,jjx
      double precision :: xmin, zmin

      iix = ix - 2*margin
      jjx = jx - 2*margin
      dx = 6.d0 / dble(iix)
      dz = 6.d0 / dble(jjx)
      xmin = dx! 0.04d0 ! 0.04d0
      zmin = 0.d0 !dz ! 0.04d0
      !xmin = dx
      !x(1) = 2*dx
      x(1,:) = xmin - 0.5d0 * dx
      z(:,1) = zmin - 0.5d0 * dz
      xm(1,:) = xmin
      zm(:,1) = zmin
      do j = 1, jx
      do i = 1, ix
      if (i>1) then
        x(i,j) = x(i-1, j) + dx
      endif
      if (j>1) then
        z(i,j) = z(i, j-1) + dz
      endif
        xm(i,j) = x(i,j) + 0.5d0 * dx
        zm(i,j) = z(i,j) + 0.5d0 * dz
        enddo
      enddo
      !xm(ix) = xm(ix-1)
    end subroutine grid2d

    subroutine sedov2d_w_b(ix,jx,x,z,rho,vx,vy,vz,p,bx,by,bz)
      integer,intent(in) :: ix,jx
      double precision,intent(in) :: x(ix,jx),z(ix,jx)
      double precision,intent(inout) :: rho(ix,jx),vx(ix,jx),vz(ix,jx),p(ix,jx),&
        vy(ix,jx),bx(ix,jx),by(ix,jx),bz(ix,jx)
      integer :: i,j
      double precision :: we,prism
      we = 0.1d0; prism = 1.d-8
      rho(:,:) = 1.d0;
      vx(:,:) = 0.d0; vz(:,:) = 0.d0; vy(:,:) = 0.d0
      bx(:,:) = 0.d0; by(:,:) = 0.d0; bz(:,:) = 1.d0
      do j = 1,jx
      do i = 1,ix
        p(i,j) = prism + (1.d0-prism)*exp(-(sqrt(x(i,j)**2+z(i,j)**2)/we)**2)
      enddo
      enddo
    end subroutine sedov2d_w_b

    
end module initial_condition
