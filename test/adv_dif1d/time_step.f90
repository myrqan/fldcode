module time_step
  implicit none
  real(8) :: safety = 0.4d0
  contains
    subroutine cfl_tc(dt,dx,gx,eta)
      integer,intent(in) :: gx
      real(8),intent(in) :: dx,eta
      real(8),intent(inout) :: dt
      !real(8) :: safety = 0.4d0

      dt = safety * dx**2 / (2.d0 * eta)

    end subroutine cfl_tc
    subroutine cfl_ad(dt,dx,gx)
      integer,intent(in) :: gx
      real(8),intent(in) :: dx
      real(8),intent(inout) :: dt
      !real(8) :: safety = 0.4d0

      dt = safety * dx

    end subroutine cfl_ad
end module time_step
