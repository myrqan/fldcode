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
end module calc
