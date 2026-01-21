module calc
  implicit none
contains
  subroutine calc_ay(ay,bz,xx,zz,lix,liz)
    integer,intent(in):: lix,liz
    real(8),intent(in):: xx(lix),zz(liz),bz(lix,liz)
    real(8),intent(inout):: ay(lix,liz)
    integer:: ii,jj

    do ii = 1,lix
    do jj = 1,liz
      ay(ii,jj) = 0.5d0 * bz(ii,jj) * xx(ii)
    end do
    end do 
  end subroutine calc_ay

  subroutine calc_efield(bx,by,bz,vx,vy,vz,ex,ey,ez)
    real(8),intent(in):: bx(:,:),by(:,:),bz(:,:),&
      vx(:,:),vy(:,:),vz(:,:)
    real(8),intent(inout):: ex(:,:),ey(:,:),ez(:,:)
    ex(:,:)=-vy(:,:)*bz(:,:)+vz(:,:)*by(:,:)
    ey(:,:)=-vz(:,:)*bx(:,:)+vx(:,:)*bz(:,:)
    ez(:,:)=-vx(:,:)*by(:,:)+vy(:,:)*bx(:,:)
  end subroutine calc_efield

  subroutine calc_etot(ro,pr,vx,vy,vz,bx,by,bz,etot,gm)
    real(8),intent(in):: ro(:,:),pr(:,:),vx(:,:),vy(:,:),vz(:,:),&
      bx(:,:),by(:,:),bz(:,:)
    real(8),intent(in):: gm
    real(8),intent(inout):: etot(:,:)

    etot(:,:)=pr(:,:)/(gm-1.d0)+&
      0.5d0*ro(:,:)*(vx(:,:)**2+vy(:,:)**2+vz(:,:)**2)+&
      0.5d0*(bx(:,:)**2+by(:,:)**2+bz(:,:)**2)
  end subroutine

  subroutine calc_pr(ro,pr,vx,vy,vz,bx,by,bz,etot,gm,lix,liz)
    integer,intent(in):: lix,liz
    real(8),intent(in):: ro(lix,liz),vx(lix,liz),vy(lix,liz),vz(lix,liz),&
      bx(lix,liz),by(lix,liz),bz(lix,liz),etot(lix,liz)
    real(8),intent(in):: gm
    real(8),intent(inout):: pr(lix,liz)

    integer :: ii,jj
    do jj=1,liz
    do ii=1,lix
      pr(ii,jj) = (gm-1.d0)&
        *(etot(ii,jj)-0.5d0*ro(ii,jj)*(vx(ii,jj)**2+vy(ii,jj)**2+vz(ii,jj)**2)&
      -0.5d0*(bx(ii,jj)**2+by(ii,jj)**2+bz(ii,jj)**2))
    end do
    end do
  end subroutine calc_pr
end module calc
