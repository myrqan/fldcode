module model
  use vars
  use datw
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
    xmax= 3.d0
    zmin= 0.d0
    zmax= 7.d0


    !--------------------------------------------------
    ! determine mpi area of responsibility
    ! todo ここをnprocsに応じて書き換える．modとか使えばいけそう
    !--------------------------------------------------
    if(my_rank==0) then 
      mplx=0;mplz=0
    else if(my_rank==1) then
      mplx=1;mplz=0
    else if(my_rank==2) then
      mplx=0;mplz=1
    else!if(my_rank==3) then
      mplx=1;mplz=1
    end if
    call MPI_BARRIER(MPI_Comm_World,ierr)

    !--------------------------------------------------
    ! write parameter in param.txt
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
    ! setup grid
    !--------------------------------------------------
    gx=gix-2*mg; lx=lix-2*mg
    gz=giz-2*mg; lz=liz-2*mg

    dx=(xmax-xmin)/real(gx-1,8)
    dz=(zmax-zmin)/real(gz-1,8)

    do ii=1,lix
      g_idx_ofs=real(mplx*lx+ii-mg-1,8)
      xx(ii)=xmin+g_idx_ofs*dx
    end do
    do ii=1,liz
      g_idz_ofs=real(mplz*lz+ii-mg-1,8)
      zz(ii)=zmin+g_idz_ofs*dz
    end do
    do ii=1,lix-1
      xm(ii)=0.5d0*(xx(ii)+xx(ii+1))
    end do
    do ii=1,liz-1
      zm(ii)=0.5d0*(zz(ii)+zz(ii+1))
    end do
    xm(lix) = xm(lix-1)+dx
    zm(liz) = zm(liz-1)+dz


    call put_grid_data(lix,liz,mg,xx,zz)


  end subroutine make_grid

end module model
