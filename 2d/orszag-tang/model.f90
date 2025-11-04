MODULE model
  IMPLICIT NONE
  CONTAINS
    subroutine make_grid(x,xm,dx,z,zm,dz,ix,jx,margin)
      INTEGER,INTENT(IN) :: ix,jx,margin
      DOUBLE PRECISION,INTENT(OUT) :: x(ix,jx),xm(ix,jx),dx,&
        z(ix,jx),zm(ix,jx),dz
      DOUBLE PRECISION :: xmin,xmax,zmin,zmax
      INTEGER :: ii,jj,iix,jjx

      !------------------------------
      ! these values can be changed
      xmax = 1.02d0; zmax = 1.d0
      xmin = 0.02d0; zmin = 0.d0
      !xmax = 1.d0; zmax = 1.d0
      !xmin = 0.d0; zmin = -1.d0
      !------------------------------
      xmax = 1.d0 + xmin
      
      iix = ix-2*margin
      jjx = jx-2*margin
      dx = (xmax-xmin) / dble(iix)
      dz = (zmax-zmin) / dble(jjx)

      ! xmin replace 
      !xmin = 2*margin*dx

      x(margin,:) = xmin - 0.5d0*dx
      z(:,margin) = zmin - 0.5d0*dz
      xm(margin,:) = xmin; zm(:,margin) = zmin

      do ii = margin-1,1,-1
        x(ii,:) = x(ii+1,:)-dx
        xm(ii,jj)=x(ii,jj)+0.5d0*dx
      enddo

      do jj = margin-1,1,-1
        z(:,jj) = z(:,jj+1)-dz
        zm(:,jj)=z(:,jj)+0.5d0*dz
      enddo

      do ii = margin, ix
        if(ii>margin) then
          x(ii,:) = x(ii-1,:)+dx
        endif
        xm(ii,:)=x(ii,:)+0.5d0*dx
      enddo

      do jj = margin, jx
        if(jj>margin) then
          z(:,jj) = z(:,jj-1)+dz
        endif
        zm(:,jj)=z(:,jj)+0.5d0*dz
      enddo
    END subroutine make_grid

    !#################################

    subroutine sedov_with_b(ro,vx,vy,vz,bx,by,bz,p,x,z,ix,jx)
      INTEGER,INTENT(IN) :: ix,jx
      DOUBLE PRECISION,INTENT(IN) :: x(ix,jx),z(ix,jx)
      DOUBLE PRECISION,INTENT(OUT) :: ro(ix,jx),vx(ix,jx),&
        vy(ix,jx),vz(ix,jx),bx(ix,jx),by(ix,jx),bz(ix,jx),p(ix,jx)
      DOUBLE PRECISION :: we, prism
      INTEGER :: ii,jj
      DOUBLE PRECISION, PARAMETER :: pi = 4.d0 * ATAN(1.d0)

      we = 0.1d0
      prism = 1.d-8

      ro(:,:)=1.d0
      vx(:,:)=0.d0; vy(:,:)=0.d0; vz(:,:)=0.d0
      bx(:,:)=0.d0; by(:,:)=0.d0; bz(:,:)=1.d-3

      do jj = 1,jx
        do ii = 1,ix
          p(ii,jj) = prism + (1.d0-prism) &
          * exp(-(x(ii,jj)**2+z(ii,jj)**2)/we**2)
        enddo
      enddo

    END subroutine sedov_with_b

END module model
