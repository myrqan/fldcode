module init
  use vars
  implicit none
contains
  subroutine mkgrd(xx,xm,dx,gix,lix,mg)
    ! xx = center of cells, xm = edge of cells
    ! xm(i) is (xx(i) + xx(i+1))/2  i.e., xx(i+1/2)
    ! dx = grid width 
    ! ix = # of cells, mg = margin
    ! gix = # of GLOBAL cells
    integer,intent(in) :: mg,gix,lix
    real(8),intent(inout) :: xx(lix), xm(lix)
    real(8),intent(inout) :: dx

    integer :: icnt,gx,lx,ii
    real(8) :: xmin,xmax,g_idx_ofs

    gx = gix-2*mg
    lx = lix-2*mg
    xmin = 0.d0
    xmax = 1.d0
    dx = (xmax-xmin)/dble(gx-1)

    do ii=1,lix
      g_idx_ofs = dble(my_rank*lx+ii-mg-1)
      xx(ii)=xmin+g_idx_ofs*dx
    end do

    do ii = 1,lix-1
      xm(ii) = 0.5d0 * (xx(ii) + xx(ii+1))
    enddo

    xm(lix) = xm(lix-1) + dx
    
  end subroutine mkgrd
end module init
