module model
  use vars
  use const
  use datw
  use bnd2d
  implicit none
contains
  subroutine make_grid(xx,xm,dx,zz,zm,dz,gix,lix,giz,liz,mg)
    integer,intent(in):: gix,lix,giz,liz,mg
    real(8),intent(inout):: xx(lix),xm(lix),zz(liz),zm(liz)
    real(8),intent(inout):: dx,dz
    real(8):: xmin,xmax,zmin,zmax
    integer:: gx,gz,lx,lz,g_idx_ofs,g_idz_ofs
    integer:: ii,jj

    !--------------------------------------------------
    ! xmin, xmax, zmin, zmax value
    !--------------------------------------------------
    xmin= 0.07d0
    !xmin=0.d0
    xmax= 3.d0
    zmin= 0.d0
    zmax= 7.d0

    !--------------------------------------------------
    ! determine mpi area of responsibility
    !--------------------------------------------------
    mplx=mod(my_rank,mpx)
    mplz=int(my_rank/mpx)
    ! call MPI_BARRIER(MPI_Comm_World,ierr)
    !--------------------------------------------------
    ! setup grid
    !--------------------------------------------------
    gx=gix-2*mg; lx=lix-2*mg
    gz=giz-2*mg; lz=liz-2*mg

    dx=(xmax-xmin)/real(gx,8)
    dz=(zmax-zmin)/real(gz,8)

    do ii=1,lix
      !g_idx_ofs=mplx*lx+ii-mg-1
      !xx(ii)=xmin+(g_idx_ofs+0.5d0)*dx
      g_idx_ofs = mplx * lx + (ii - mg)
      xx(ii) = xmin + (dble(g_idx_ofs) - 0.5d0) * dx
    end do
    do jj=1,liz
      !g_idz_ofs=mplz*lz+jj-mg-1
      !zz(jj)=zmin+(g_idz_ofs+0.5d0)*dz
      g_idz_ofs = mplz * lz + (jj - mg)
      zz(jj) = zmin + (dble(g_idz_ofs) - 0.5d0) * dz
    end do

    do ii=1,lix-1
      xm(ii)=0.5d0*(xx(ii)+xx(ii+1))
    end do
    do jj=1,liz-1
      zm(jj)=0.5d0*(zz(jj)+zz(jj+1))
    end do
    xm(lix) = xm(lix-1)+dx
    zm(liz) = zm(liz-1)+dz

    !--------------------------------------------------
    ! write parameter on param.txt
    !--------------------------------------------------
    if(my_rank==0) then
      call put_param_real("xmin:",xmin)
      call put_param_real("xmax:",xmax)
      call put_param_real("zmin:",zmin)
      call put_param_real("zmax:",zmax)
      call put_param_real("dx:",dx)
      call put_param_real("dz:",dz)
    end if
    !--------------------------------------------------
    ! write grid data on x.dat and z.dat
    !--------------------------------------------------
    call put_grid_data(lix,liz,mg,xx,zz)
  end subroutine make_grid

  subroutine model_diskjet(&
      ro,vx,vy,vz,bx,by,bz,pr,&
      gm,grx,grxm,grz,grzm,xx,zz,&
      lix,liz,dx,dz,mg)
    integer,intent(in):: lix,liz,mg
    real(8),intent(in):: gm,dx,dz
    real(8),intent(in):: xx(lix),zz(liz)
    real(8),intent(inout):: ro(lix,liz),&
      vx(lix,liz),vy(lix,liz),vz(lix,liz),&
      bx(lix,liz),by(lix,liz),bz(lix,liz),&
      pr(lix,liz),grx(lix,liz),grxm(lix,liz),&
      grz(lix,liz),grzm(lix,liz)

    real(8):: gpot(lix,liz),dis
    real(8):: srad,aa,nn,alpha,roc,eth,emg,tec0,tec00
    real(8):: psi0,b0,ro_c,ro_d,pr_c,pr_d,vy_c,vy_d,te

    integer:: ii,jj
    !--------------------------------------------------
    ! define parameters
    !--------------------------------------------------
    ! softening radius
    srad = 0.2d0
    ! distribution of angular momentum and pressure
    ! L = L_0 r^a, p = K rho^{1+1/n}
    aa = 0.d0
    nn = 3.d0
    ! ratio of Keplerian velocity @ r_0 vs. 
    ! sound speed in the corona
    alpha = 1.d0
    ! temperature
    tec0 = 1.d0
    tec00 = tec0*gm
    ! coronal density at radius r0
    roc =1.d-3
    ! non-dimensional parameters
    ! 1. sound speed vs. Keplarian velocity
    ! 2. Alfven speed vs. Keplarian velocity
    eth = 5.d-2
    emg = 5.d-4

    !--------------------------------------------------
    ! write parameters on param.txt
    !--------------------------------------------------
    if(my_rank==0) then
      call put_param_real("srad:",srad)
      call put_param_real("aa:",aa)
      call put_param_real("nn:",nn)
      call put_param_real("alpha:",alpha)
      call put_param_real("tec0:",tec0)
      call put_param_real("roc:",roc)
      call put_param_real("eth:",eth)
      call put_param_real("emg:",emg)
    end if

    !--------------------------------------------------
    ! calculate gravitational potential and acceleration
    !--------------------------------------------------
    do jj = 1,liz
    do ii = 1,lix
      dis = sqrt(xx(ii)**2+zz(jj)**2)
      if(dis>srad) then
        gpot(ii,jj)=-1.d0/dis
      else if (dis>0.5d0*srad) then
        gpot(ii,jj)=-(2.d0/srad-dis/srad**2)
      else 
        gpot(ii,jj)=-1.5d0/srad
      end if
    end do
    end do

    do jj = 1,liz-1
    do ii = 1,lix-1
      if(ii>=2 .and. jj>=2) then
        grx(ii,jj)=-(gpot(ii+1,jj)-gpot(ii-1,jj))/(2.d0*dx)
        grz(ii,jj)=-(gpot(ii,jj+1)-gpot(ii,jj-1))/(2.d0*dz)
      end if
      grxm(ii,jj)=-(gpot(ii+1,jj)-gpot(ii,jj))/dx
      grzm(ii,jj)=-(gpot(ii,jj+1)-gpot(ii,jj))/dz
    end do
    end do
    
    !--------------------------------------------------
    ! apply boundary condition for grx,grz,grxm,grzm
    !--------------------------------------------------
    call bnd_2d_symm(grx,0,mg,lix,liz)
    call bnd_2d_free(grx,1,mg,lix,liz)
    call bnd_2d_free(grx,2,mg,lix,liz)
    call bnd_2d_asym(grx,3,mg,lix,liz)

    call bnd_2d_asym(grz,0,mg,lix,liz)
    call bnd_2d_free(grz,1,mg,lix,liz)
    call bnd_2d_free(grz,2,mg,lix,liz)
    call bnd_2d_symm(grz,3,mg,lix,liz)

    call bnd_2d_symm(grxm,0,mg,lix,liz)
    call bnd_2d_free(grxm,1,mg,lix,liz)
    call bnd_2d_free(grxm,2,mg,lix,liz)
    call bnd_2d_asym(grxm,3,mg,lix,liz)

    call bnd_2d_asym(grzm,0,mg,lix,liz)
    call bnd_2d_free(grzm,1,mg,lix,liz)
    call bnd_2d_free(grzm,2,mg,lix,liz)
    call bnd_2d_symm(grzm,3,mg,lix,liz)
    !--------------------------------------------------
    ! calculate other variables
    !--------------------------------------------------
    psi0 = -1.d0+0.5d0/(1.d0-aa)+(nn+1.d0)*eth
    b0=sqrt(emg)

    do jj=1,liz
    do ii=1,lix
      dis=sqrt(xx(ii)**2+zz(jj)**2)
      ro_c=roc*exp(-1.d0/tec0*(gpot(ii,jj)+1.d0))
      pr_c=ro_c*tec0
      vy_c=0.d0

      ro_d=0.d0
      pr_d=0.d0
      vy_d=0.d0
      if(xx(ii)>srad) then
        te=(psi0+1.d0/dis-0.5d0/(1.d0-aa)*xx(ii)**(2.d0*aa-2.d0))&
          /(nn+1)*gm
        if(te>0.d0) then
          ro_d=(te/gm/eth)**nn
          pr_d=te*ro_d/gm
          vy_d=xx(ii)**(aa-1)
        end if
      end if
      ro(ii,jj)=ro_c+ro_d
      pr(ii,jj)=pr_c+pr_d
      vy(ii,jj)=vy_c+vy_d
    end do
    end do
    vx(:,:)=0.d0;vz(:,:)=0.d0
    bx(:,:)=0.d0;by(:,:)=0.d0
    bz(:,:)=b0

    !--------------------------------------------------
    ! apply boundary condition(s)
    !--------------------------------------------------
    call bnd_2d_symm(ro,0,mg,lix,liz)
    call bnd_2d_free(ro,1,mg,lix,liz)
    call bnd_2d_free(ro,2,mg,lix,liz)
    call bnd_2d_symm(ro,3,mg,lix,liz)

    call bnd_2d_symm(vx,0,mg,lix,liz)
    call bnd_2d_free(vx,1,mg,lix,liz)
    call bnd_2d_free(vx,2,mg,lix,liz)
    call bnd_2d_asym(vx,3,mg,lix,liz)

    call bnd_2d_symm(vy,0,mg,lix,liz)
    call bnd_2d_free(vy,1,mg,lix,liz)
    call bnd_2d_free(vy,2,mg,lix,liz)
    call bnd_2d_asym(vy,3,mg,lix,liz)

    call bnd_2d_asym(vz,0,mg,lix,liz)
    call bnd_2d_free(vz,1,mg,lix,liz)
    call bnd_2d_free(vz,2,mg,lix,liz)
    call bnd_2d_symm(vz,3,mg,lix,liz)

    call bnd_2d_asym(bx,0,mg,lix,liz)
    call bnd_2d_free(bx,1,mg,lix,liz)
    call bnd_2d_free(bx,2,mg,lix,liz)
    call bnd_2d_asym(bx,3,mg,lix,liz)

    call bnd_2d_asym(by,0,mg,lix,liz)
    call bnd_2d_free(by,1,mg,lix,liz)
    call bnd_2d_free(by,2,mg,lix,liz)
    call bnd_2d_asym(by,3,mg,lix,liz)

    call bnd_2d_symm(bz,0,mg,lix,liz)
    call bnd_2d_free(bz,1,mg,lix,liz)
    call bnd_2d_free(bz,2,mg,lix,liz)
    call bnd_2d_symm(bz,3,mg,lix,liz)

    call bnd_2d_symm(pr,0,mg,lix,liz)
    call bnd_2d_free(pr,1,mg,lix,liz)
    call bnd_2d_free(pr,2,mg,lix,liz)
    call bnd_2d_symm(pr,3,mg,lix,liz)
  end subroutine model_diskjet

end module model
