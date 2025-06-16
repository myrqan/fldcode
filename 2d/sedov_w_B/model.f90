MODULE model
  IMPLICIT NONE
  CONTAINS
    subroutine make_grid(x,xm,dx,z,zm,dz,ix,jx,margin)
      INTEGER,INTENT(IN) :: ix,jx,margin
      DOUBLE PRECISION,INTENT(INOUT) :: x(ix,jx),xm(ix,jx),dx,&
        z(ix,jx),zm(ix,jx),dz
      DOUBLE PRECISION :: xmin,xmax,zmin,zmax
      INTEGER :: ii,jj,iix,jjx

      !------------------------------
      ! these values can be changed
      xmax = 1.0d0; zmax = 1.d0
      xmin = 0.04d0; zmin = 0.d0
      !------------------------------
      
      iix = ix-2*margin
      jjx = jx-2*margin
      dx = xmax / dble(iix)
      dz = zmax / dble(jjx)


      x(1,:) = xmin - 0.5d0*dx
      z(:,1) = zmin - 0.5d0*dz
      xm(1,:) = xmin; zm(:,1) = zmin

      do jj = 1, jx
        do ii = 1, ix
          if(ii>1) then
            x(ii,jj) = x(ii-1,jj)+dx
          endif
          if(jj>1) then
            z(ii,jj) = z(ii,jj-1)+dz
          endif
          xm(ii,jj)=x(ii,jj)+0.5d0*dx
          zm(ii,jj)=z(ii,jj)+0.5d0*dz
        enddo
      enddo
    END subroutine make_grid

    !#################################

    subroutine sedov_with_b(ro,vx,vy,vz,bx,by,bz,p,x,z,ix,jx)
      INTEGER,INTENT(IN) :: ix,jx
      DOUBLE PRECISION,INTENT(IN) :: x(ix,jx),z(ix,jx)
      DOUBLE PRECISION,INTENT(INOUT) :: ro(ix,jx),vx(ix,jx),&
        vy(ix,jx),vz(ix,jx),bx(ix,jx),by(ix,jx),bz(ix,jx),p(ix,jx)
      DOUBLE PRECISION :: we, prism
      INTEGER :: ii,jj

      we = 0.1d0
      prism = 1.d-8

      ro(:,:)=1.d0
      vx(:,:)=0.d0; vy(:,:)=0.d0; vz(:,:)=0.d0
      bx(:,:)=0.d0; by(:,:)=0.d0; bz(:,:)=1.d0

      do jj = 1,jx
        do ii = 1,ix
          p(ii,jj) = prism + (1.d0-prism) &
          * exp(-(sqrt(x(ii,jj)**2+z(ii,jj)**2)/we)**2)
        enddo
      enddo

    END subroutine sedov_with_b

END module model
