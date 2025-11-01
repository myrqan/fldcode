MODULE model
  IMPLICIT NONE
  CONTAINS
    subroutine make_grid(x,xm,dx,z,zm,dz,ix,jx,margin)
      USE file_output
      INTEGER,INTENT(IN) :: ix,jx,margin
      DOUBLE PRECISION,INTENT(OUT) :: x(ix,jx),xm(ix,jx),dx,&
        z(ix,jx),zm(ix,jx),dz
      DOUBLE PRECISION :: xmin,xmax,zmin,zmax
      INTEGER :: ii,jj,iix,jjx

      !------------------------------
      ! these values can be changed
      xmax = 3.d0; zmax = 7.d0
      xmin = 0.07d0; zmin = 0.00d0
      !xmax = 1.d0; zmax = 1.d0
      !xmin = 0.d0; zmin = -1.d0
      !------------------------------
      !xmax = 1.d0 + xmin
      CALL put_param_dble("xmin:",xmin)
      CALL put_param_dble("xmax:",xmax)
      CALL put_param_dble("zmin:",zmin)
      CALL put_param_dble("zmax:",zmax)
      
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

    subroutine diskjet(ro,vx,vy,vz,bx,by,bz,p,gm,gx,gxm,gz,gzm,x,z,ix,jx,dx,dz)
      USE file_output
      INTEGER,INTENT(IN) :: ix,jx
      DOUBLE PRECISION,INTENT(IN) :: x(ix,jx),z(ix,jx),dx,dz,gm
      DOUBLE PRECISION,INTENT(OUT) :: ro(ix,jx),vx(ix,jx),&
        vy(ix,jx),vz(ix,jx),bx(ix,jx),by(ix,jx),bz(ix,jx),p(ix,jx),&
        gx(ix,jx),gxm(ix,jx),gz(ix,jx),gzm(ix,jx)
      DOUBLE PRECISION :: gpot(ix,jx),dis
      DOUBLE PRECISION :: srad,aa,nn,alpha,roc,eth,emg,tec0,tec00
      DOUBLE PRECISION :: psi0,b0,ro_c,ro_d,p_c,p_d,vy_c,vy_d,te
      INTEGER :: i,j
      DOUBLE PRECISION, PARAMETER :: pi = 4.d0 * ATAN(1.d0)

      !! def. parameters
      srad = 0.2d0
      aa = 0.d0
      nn = 3.d0
      alpha = 1.d0
      tec0 = 1.d0
      tec00 = tec0 * gm
      roc = 1.d-3
      eth = 5.d-2
      emg = 5.d-4
      !emg = 4.d-3

      CALL put_param_dble("srad:",srad)
      CALL put_param_dble("aa:",aa)
      CALL put_param_dble("nn:",nn)
      CALL put_param_dble("alpha:",alpha)
      CALL put_param_dble("roc:",roc)
      CALL put_param_dble("eth:",eth)
      CALL put_param_dble("emg:",emg)
      
      !! calculate gravitational potential and gravitational acceration
      !! softening radius

      do j=1,jx
      do i=1,ix
        dis = sqrt(x(i,j)**2 + z(i,j)**2)

        if(dis > srad) then
          gpot(i,j) = -1.d0 / dis
        else if (dis > 0.5d0 * srad) then
          gpot(i,j) = -(2.d0 / srad - dis / srad**2)
        else
          gpot(i,j) = -1.5d0 / srad
        endif
      enddo
      enddo

      do j=1,jx-1
      do i=1,ix-1
        if(i >= 2 .AND. j >= 2) then
          gx(i,j) = -(gpot(i+1,j) - gpot(i-1,j)) / (2.d0 * dx)
          gz(i,j) = -(gpot(i,j+1) - gpot(i,j-1)) / (2.d0 * dz)
        endif
        gxm(i,j) = -(gpot(i+1,j) - gpot(i,j)) / dx
        gzm(i,j) = -(gpot(i,j+1) - gpot(i,j)) / dz
      enddo
      enddo

      !! apply boudary conditions for g{x,z}(m) (in the main.f90)

      ! calculate other variables

      psi0 = -1.d0 + 0.5d0 / (1.d0-aa) + (nn+1.d0) * eth
      b0 = sqrt(emg)

      do j=1,jx
      do i=1,ix
        dis = sqrt(x(i,j)**2+z(i,j)**2)
        ro_c = roc * exp(-1.d0/tec0*(gpot(i,j)+1.d0))
        p_c = ro_c * tec0

        ro_d = 0.d0
        p_d = 0.d0
        vy_d = 0.d0
        if(x(i,j) > srad) then

          te = (psi0 + 1.d0/dis - 0.5d0/(1.d0-aa)*x(i,j)**(2*aa-2)) / (nn+1) * gm

          if(te > 0.d0) then

            ro_d = (te/gm/eth)**nn
            p_d = te*ro_d/gm
            vy_d = x(i,j)**(aa-1)

          endif
        endif
        ro(i,j) = ro_c + ro_d
        p(i,j) = p_c + p_d
        vy(i,j) = vy_c + vy_d

      !--- ---!

!        dis = sqrt(x(i,j)**2 + z(i,j)**2)
!        ro_c = roc * exp(-alpha * (gpot(i,j) + 1.d0))
!        ro_d = ( 1/(eth*(nn+1)) &
!          * ( psi0 - gpot(i,j) - 1.d0/(2*(1-aa)) * x(i,j)**(2*aa-2) ) )**nn
!        
!
!        if(ro_d < 0.d0) then
!          !---- outside disk ----!
!          ro(i,j) = ro_c
!          vy(i,j) = 0.d0
!        else
!          !---- inside disk ----!
!          ro(i,j) = ro_c + ro_d
!          vy(i,j) = x(i,j) ** (aa-1.d0)
!        endif

      enddo
      enddo

!      p(:,:) = eth * ro(:,:) ** (1.d0 + 1.d0/nn)

      vx(:,:) = 0.d0; vz(:,:) = 0.d0
      bx(:,:) = 0.d0; by(:,:) = 0.d0
      bz(:,:) = b0
      !bz(:,:) = 0.d0

    END subroutine diskjet


END module model
