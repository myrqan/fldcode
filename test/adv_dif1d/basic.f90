module basic
  implicit none
  contains
    subroutine makegrid(xx,dx,xmax,xmin,gx)
      integer,intent(in) :: gx
      real(8),intent(in) :: xmax,xmin
      real(8),intent(inout) :: xx(gx),dx
      integer :: ii ! counter

      dx = (xmax - xmin) / (gx - 1)
      xx(1) = xmin
      do ii = 2, gx
      xx(ii) = xx(ii-1) + dx
      end do
    end subroutine makegrid

    subroutine initial_condition(qq,xx,gx)
      integer,intent(in) :: gx
      real(8),intent(in) :: xx(gx)
      real(8),intent(inout) :: qq(gx)
      real(8) :: width = 0.05d0
      real(8) :: mid = 0.5d0
      integer :: ii

      do ii = 1, gx
        qq(ii) = exp(-((xx(ii)-mid)/width)**2)
      end do

    end subroutine initial_condition
end module basic
