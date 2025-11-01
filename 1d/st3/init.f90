module init
  implicit none
contains
  subroutine mkgrd(xx,xm,dx,ix,mg)
    ! xx = center of cells, xm = edge of cells
    ! xm(i) is (xx(i) + xx(i+1))/2
    ! dx = grid width 
    ! ix = # of cells, mg = margin
    integer,intent(in) :: ix,mg
    real(8),intent(inout) :: xx(ix), xm(ix)
    real(8),intent(inout) :: dx

    integer :: ii, r_ix
    real(8) :: xmin,xmax

    r_ix = ix - 2 * mg ! real ix (w/o margin)
    xmin = 0.d0
    xmax = 1.d0
    dx = (xmax - xmin) / dble(r_ix)
    
    xx(mg+1) = xmin

    do ii = mg+2, ix
    xx(ii) = xx(ii-1) + dx
    enddo

    do ii = mg,1,-1
    xx(ii) = xx(ii+1) - dx
    enddo

    do ii = 1,ix-1
    xm(ii) = 0.5d0 * (xx(ii) + xx(ii+1))
    enddo

    xm(ix) = xm(ix-1) + dx
    
  end subroutine mkgrd
end module init
