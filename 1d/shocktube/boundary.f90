module boundary
  implicit none
  contains
    subroutine bnd1d_f_lr(ary,ix)
      integer,intent(in) :: ix
      double precision, intent(inout) :: ary(ix)
      ary(1) = ary(2); ary(ix)= ary(ix-1)
    end subroutine bnd1d_f_lr
end module boundary
