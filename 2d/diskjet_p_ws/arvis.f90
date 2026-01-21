module arvis
  implicit none
contains
  subroutine calc_coef_k(qv,kx,kz,vx,vz,ix,jx,dx,dz)
    integer,intent(in) :: ix,jx
    double precision,intent(in) :: qv,vx(ix,jx),vz(ix,jx),dx,dz
    double precision,intent(inout) :: kx(ix,jx), kz(ix,jx)
    integer :: i,j
    kx(:,:)=0.d0
    kz(:,:)=0.d0
    do j = 1,jx-1
    do i = 1,ix-1
    kx(i,j) = qv*dx*max(abs(vx(i+1,j)-vx(i,j))-1.0d-4, 0.d0)
    kz(i,j) = qv*dz*max(abs(vz(i,j+1)-vz(i,j))-1.0d-4, 0.d0)
    enddo
    enddo
  end subroutine calc_coef_k

  subroutine arvis2d(kx,kz,u,du,ix,jx,dt,dx,dz)
    integer,intent(in) :: ix,jx
    double precision,intent(in) :: kx(ix,jx),kz(ix,jx),u(ix,jx),dt,dx,dz
    double precision,intent(inout) :: du(ix,jx)
    integer :: i,j
    !du(:,:)=0.d0
    do j=2,jx
    do i=2,ix
    du(i,j)=du(i,j) &
      +dt/dx* ((kx(i  ,j)/dx * (u(i+1,j)-u(i,  j))&
      - kx(i-1,j)/dx * (u(i  ,j)-u(i-1,j))))&
      +dt/dz* ((kz(i  ,j)/dz * (u(i,j+1)-u(i,  j))&
      - kz(i,j-1)/dz * (u(i,  j)-u(i,j-1))))
    enddo
    enddo
  end subroutine arvis2d
end module arvis
