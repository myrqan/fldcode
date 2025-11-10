module sts
  implicit none
  contains
    subroutine determine_sts_s(ss,dx,eta)
      integer,intent(inout) :: ss
      real(8),intent(in) :: dx,eta

      integer :: trys
      real(8) :: t
      
      t = 2.d0 * eta / dx

      trys = ceiling(0.5d0 * (-1 + sqrt(16*t+9)))
      if(mod(trys,2) == 1) then
        ss = trys
      else
        ss = trys + 1
      end if


    end subroutine determine_sts_s
end module sts
