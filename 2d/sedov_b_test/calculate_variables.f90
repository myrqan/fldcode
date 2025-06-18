module calculate_variables
implicit none
contains
subroutine calc_e(bx,by,bz,vx,vy,vz,ex,ey,ez)
  DOUBLE PRECISION,INTENT(IN) :: bx(:,:),by(:,:),bz(:,:),&
    vx(:,:),vy(:,:),vz(:,:)
  DOUBLE PRECISION,INTENT(INOUT) :: ex(:,:),ey(:,:),ez(:,:)

  ex(:,:) = -vy(:,:)*bz(:,:) + vz(:,:)*by(:,:)
  ey(:,:) = -vz(:,:)*bx(:,:) + vx(:,:)*bz(:,:)
  ez(:,:) = -vx(:,:)*by(:,:) + vy(:,:)*bx(:,:)
end subroutine calc_e

subroutine calc_etot(ro,p,vx,vy,vz,bx,by,bz,etot,gm,ix,jx)
  DOUBLE PRECISION,INTENT(IN) :: ro(:,:),p(:,:),&
    vx(:,:),vy(:,:),vz(:,:),bx(:,:),by(:,:),bz(:,:),gm
  DOUBLE PRECISION,INTENT(INOUT) :: etot(:,:)
  INTEGER,INTENT(IN) :: ix,jx
  DOUBLE PRECISION :: v2(ix,jx), b2(ix,jx)

  v2(:,:) = vx(:,:)**2+vy(:,:)**2+vz(:,:)**2
  b2(:,:) = bx(:,:)**2+by(:,:)**2+bz(:,:)**2

  etot(:,:) = p(:,:)/(gm-1.d0) + 0.5d0*ro(:,:)*v2(:,:)&
    +0.5d0*b2(:,:)
end subroutine calc_etot

subroutine calc_p(ro,p,vx,vy,vz,bx,by,bz,etot,gm,ix,jx)
  DOUBLE PRECISION,INTENT(IN) :: ro(:,:),vx(:,:),vy(:,:),&
    vz(:,:),bx(:,:),by(:,:),bz(:,:),etot(:,:),gm
  DOUBLE PRECISION,INTENT(INOUT) :: p(:,:)
  INTEGER,INTENT(IN) :: ix,jx
  DOUBLE PRECISION :: v2(ix,jx), b2(ix,jx)

  v2(:,:) = vx(:,:)**2+vy(:,:)**2+vz(:,:)**2
  b2(:,:) = bx(:,:)**2+by(:,:)**2+bz(:,:)**2

  p(:,:) = (gm-1.d0)&
    *(etot(:,:)-0.5d0*ro(:,:)*v2(:,:)-0.5d0*b2(:,:))
end subroutine calc_p
end module calculate_variables
