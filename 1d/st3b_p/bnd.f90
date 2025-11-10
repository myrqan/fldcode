module bnd
  implicit none
contains
  subroutine bnd1d(qq,pm,ix)
    integer,intent(in) :: pm,ix
    real(8),intent(inout) :: qq(ix)
    if(pm == 0) then 
      ! left fix (00)
      qq(1) = qq(2)
    else if(pm == 10) then
      ! right fix (10)
      qq(ix) = qq(ix-1)
    endif
  end subroutine
end module bnd
