module boundary2d
  implicit none
  contains
!===========================
! boudnary condititon for Left/Right/Upper/Bottom
!===========================

    subroutine bc_free_l(ary,ix,jx)
      INTEGER,INTENT(IN)::ix,jx
      DOUBLE PRECISION,INTENT(OUT)::ary(ix,jx)
      ary(1,:)=ary(2,:)
    end subroutine bc_free_l

    subroutine bc_free_r(ary,ix,jx)
      INTEGER,INTENT(IN)::ix,jx
      DOUBLE PRECISION,INTENT(OUT)::ary(ix,jx)
      ary(ix,:)=ary(ix-1,:)
    end subroutine bc_free_r

    subroutine bc_free_u(ary,ix,jx)
      INTEGER,INTENT(IN)::ix,jx
      DOUBLE PRECISION,INTENT(OUT)::ary(ix,jx)
      ary(:,jx)=ary(:,jx-1)
    end subroutine bc_free_u

    subroutine bc_free_b(ary,ix,jx)
      INTEGER,INTENT(IN)::ix,jx
      DOUBLE PRECISION,INTENT(OUT)::ary(ix,jx)
      ary(:,1)=ary(:,2)
    end subroutine bc_free_b


    subroutine bc_fix_l(ary,ix,jx)
      INTEGER,INTENT(IN)::ix,jx
      DOUBLE PRECISION,INTENT(OUT)::ary(ix,jx)
      ary(1,:)=-ary(2,:)
    end subroutine bc_fix_l

    subroutine bc_fix_r(ary,ix,jx)
      INTEGER,INTENT(IN)::ix,jx
      DOUBLE PRECISION,INTENT(OUT)::ary(ix,jx)
      ary(ix-1,:)=-ary(ix,:)
    end subroutine bc_fix_r

    subroutine bc_fix_u(ary,ix,jx)
      INTEGER,INTENT(IN)::ix,jx
      DOUBLE PRECISION,INTENT(OUT)::ary(ix,jx)
      ary(:,jx)=-ary(:,jx-1)
    end subroutine bc_fix_u

    subroutine bc_fix_b(ary,ix,jx)
      INTEGER,INTENT(IN)::ix,jx
      DOUBLE PRECISION,INTENT(OUT)::ary(ix,jx)
      ary(:,1)=-ary(:,2)
    end subroutine bc_fix_b


end module boundary2d
