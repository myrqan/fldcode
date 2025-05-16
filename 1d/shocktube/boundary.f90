module boundary
  implicit none
  contains
    subroutine bnd1d_f_lr(ary,ix)
      integer,intent(in) :: ix
      double precision, intent(inout) :: ary(0:ix)
      ary(0) = ary(1); ary(ix)= ary(ix-1)
    end subroutine bnd1d_f_lr
end module boundary
