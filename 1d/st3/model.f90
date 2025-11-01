module model
  implicit none
contains  

  subroutine model_shocktube(ro,vx,pr,ei,ix,xx,gm)
    ! setup numerical model (initial condition)
    integer,intent(in) :: ix
    real(8),intent(in) :: xx(ix)
    real(8),intent(in) :: gm
    real(8),intent(inout) :: ro(ix),vx(ix),pr(ix),ei(ix)
    
    integer :: ii
    do ii = 1,ix
    if(xx(ii) <= 0.5d0) then
      ro(ii) = 1.d0
      vx(ii) = 0.d0
      pr(ii) = 1.d0
    else
      ro(ii) = 0.125d0
      vx(ii) = 0.d0
      pr(ii) = 0.1d0
    endif
    enddo
    ei(:) = 0.5d0*ro(:)*vx(:)**2 + pr(:)/(gm-1.d0)

  end subroutine model_shocktube
end module
