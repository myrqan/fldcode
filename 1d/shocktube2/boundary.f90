module boundary
  implicit none
  contains
    subroutine bnd1d_f_lr(ary,ix)
      integer,intent(in) :: ix
      double precision, intent(inout) :: ary(ix)
      ary(3) = ary(4); ary(2) = ary(3)
      ary(1) = ary(2); ary(ix)= ary(ix-1)
    end subroutine bnd1d_f_lr

    subroutine bnd1d_f_r(ary,ix)
      integer,intent(in) :: ix
      double precision, intent(inout) :: ary(ix)
      ary(ix)= ary(ix-1)
    end subroutine bnd1d_f_r

    subroutine bnd1d_fix_l(ary,ix)
      integer,intent(in) :: ix
      double precision,intent(inout) :: ary(ix)
      ary(1) = 0.d0; ary(2) = 0.d0
    end subroutine
end module boundary
