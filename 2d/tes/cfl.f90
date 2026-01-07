module cfl
  use mpi
  use vars
  use datw
  implicit none
contains
  subroutine calc_dt(&
      ro,vx,vy,vz,pr,bx,by,bz,gm,lix,liz,dt,dx,dz,mg)
    integer,intent(in):: lix,liz,mg
    real(8),intent(in):: dx,dz,gm
    real(8),intent(in):: ro(lix,liz),vx(lix,liz),vy(lix,liz),&
      vz(lix,liz),pr(lix,liz),bx(lix,liz),by(lix,liz),bz(lix,liz)
    real(8),intent(inout):: dt
    real(8):: safety=0.4d0
    real(8):: dtmin=1.d-20
    real(8):: cs2,va2,v2,dtcfl,min_val,gmin_val
    integer:: ii,jj

    if(my_rank==0) then
      call put_param_real("cfl safety:",safety)
    end if

    min_val=1.d20
    do jj = mg,liz-mg
    do ii = mg,lix-mg
      v2 = vx(ii,jj)**2+vy(ii,jj)**2+vz(ii,jj)**2
      va2 = (bx(ii,jj)**2+by(ii,jj)**2+bz(ii,jj)**2)/ro(ii,jj)
      cs2 = gm*pr(ii,jj)/ro(ii,jj)
      dtcfl=min(dx,dz)/sqrt(v2+va2+cs2)
      min_val=min(min_val,dtcfl)
    end do
    end do
    call MPI_Allreduce(min_val,gmin_val,1,MPI_Double_Precision,&
      MPI_MIN,MPI_Comm_World,ierr)
    dt=gmin_val*safety

    ! check if time interval is not too small
    if(dt < dtmin) then
      if(my_rank==0) then
        write(*,*) 'dt is too small: dt<1e-20'
      end if
      call MPI_FINALIZE(ierr)
      stop
    end if
  end subroutine calc_dt
end module cfl
