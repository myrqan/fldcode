module update
  implicit none
  contains
    subroutine advflow(qq,dx,dt,gx)
      integer,intent(in) :: gx
      real(8),intent(in) :: dx,dt
      real(8),intent(inout) :: qq(gx)

      integer :: ii
      real(8) :: qqh(gx),qqn(gx)

      do ii = 1, gx-1
        qqh(ii) = 0.5d0 * (qq(ii+1) + qq(ii)) &
          &- 0.5d0*dt/dx * (qq(ii+1) - qq(ii))
      end do
      
      !! bnd condition
      qqh(1) = 0.d0
      qqh(gx) = 0.d0

      do ii = 2, gx
        qqn(ii) = qq(ii) - dt/dx * (qqh(ii) - qqh(ii-1))
      end do

      qqn(1) = 0.d0
      qqn(gx) = 0.d0

      qq(:) = qqn(:)


    end subroutine advflow

end module update
