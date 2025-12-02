module model
  implicit none
contains  

  subroutine model_shocktube(ro,vx,vy,bx,by,pr,ez,etot,ix,xx,gm)
    ! setup numerical model (initial condition)
    integer,intent(in) :: ix
    real(8),intent(in) :: xx(ix)
    real(8),intent(in) :: gm
    real(8),intent(inout) :: ro(ix),vx(ix),vy(ix),bx(ix),by(ix),&
      & pr(ix),etot(ix),ez(ix)
    
    integer :: ii
    real(8),parameter :: pi = 4.d0*atan(1.d0)
    real(8),parameter :: pii = 1.d0 / pi
    do ii = 1,ix
    if(xx(ii) <= 0.5d0) then
      ro(ii) = 1.d0
      vx(ii) = 0.d0
      vy(ii) = 0.d0
      !bx(ii) = 0.75d0*sqrt(4.d0*pi)
      bx(ii) = 0.d0
      !by(ii) = 1.d0*sqrt(4.d0*pi)
      by(ii) = 0.d0
      pr(ii) = 1.d0
    else
      ro(ii) = 0.125d0
      vx(ii) = 0.d0
      vy(ii) = 0.d0
      !bx(ii) = 0.75d0*sqrt(4.d0*pi)
      bx(ii) = 0.d0
      !by(ii) = -1.d0*sqrt(4.d0*pi)
      by(ii) = 0.d0
      pr(ii) = 0.1d0
    endif
    enddo
    etot(:) = 0.5d0*ro(:)*(vx(:)**2 + vy(:)**2) + pr(:)/(gm-1.d0) &
      & + (bx(:)**2 + by(:)**2) * 0.125d0 * pii
    ez(:) = -vx(:)*by(:) + vy(:)*bx(:)

  end subroutine model_shocktube
end module
