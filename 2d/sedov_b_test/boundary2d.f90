module boundary2d
  implicit none
contains
  !
  !
  ! (zmax) |------------------|
  !        |         2        |
  !        |                  |
  !        |                  |
  !        | 3              1 |
  !        |                  |
  !        |                  |
  !        |         0        |
  ! (zmin) |------------------|
  !      (xmin)            (xmax)
  !  0 to 3 represents boundary surface
  !
  ! name: bc(val,num1,num2,margin,ix,jx)
  ! num1 : surface number (0 to 3)
  ! num2 : symmetric same sign(0) or symmetric different sign(1) or free(2)

  subroutine bc(val,num1,num2,margin,ix,jx)
    INTEGER,INTENT(IN) :: ix,jx,margin
    INTEGER,INTENT(IN) :: num1,num2
    DOUBLE PRECISION,INTENT(INOUT) :: val(ix,jx)
    INTEGER :: ii,jj

    if(num1 == 0) then  !! at z = zmin
      condition0:&
        select case(num2)
      case(0) !! same
        do jj = 0,margin-1
        val(:,margin-jj) = val(:,margin+jj+1)
        enddo
      case(1) !! different
        do jj = 0,margin-1
        val(:,margin-jj) = - val(:,margin+jj+1)
        enddo
      case(2) !! free boundary condition
        do jj = margin,1,-1
        val(:,jj) = val(:,jj+1)
        enddo
      case default
        WRITE(*,*) "ILLEGAL BOUNDARY CONDITION : num1=0"
      end select condition0
!----------------------------
    else if(num1 == 1) then !! at x = xmax
      condition1:&
        select case(num2)
      case(0) !! same
        do ii = 1,margin
        val(ix-margin+ii,:) = val(ix-margin-ii+1,:)
        enddo
      case(1) !! different
        do ii = 1,margin
        val(ix-margin+ii,:) = - val(ix-margin-ii+1,:)
        enddo
      case(2) !! free boundary condition
        do ii = 1,margin
        val(ix-margin+ii,:) = val(ix-margin+ii-1,:)
        enddo
      case default
        WRITE(*,*) "ILLEGAL BOUNDARY CONDITION : num1=1"
      end select condition1
!----------------------------
    else if (num1 == 2) then !! at z = zmax
      condition2:&
        select case(num2)
      case(0) !! same
        do jj = 1,margin
        val(:,jx-margin+jj) = val(:,jx-margin-jj+1)
        enddo
      case(1) !! different
        do jj = 1,margin
        val(:,jx-margin+jj) = - val(:,jx-margin-jj+1)
        enddo
      case(2) !! free boundary condition
        do jj = 1,margin
        val(:,jx-margin+jj) = val(:,jx-margin+jj-1)
        enddo
      case default
        WRITE(*,*) "ILLEGAL BOUNDARY CONDITION : num1=2"
      end select condition2
!----------------------------
    else if (num1 == 3) then !! at x = xmin
      condition3:&
        select case(num2)
      case(0) !! same
        do ii = 0,margin-1
        val(margin-ii,:) = val(margin+ii+1,:)
        enddo
      case(1) !! different
        do ii = 0,margin-1
        val(margin-ii,:) = - val(margin+ii+1,:)
        enddo
      case(2) !! free boundary condition
        do ii = margin,1,-1
        val(ii,:) = val(ii+1,:)
        enddo
      case default
        WRITE(*,*) "ILLEGAL BOUNDARY CONDITION : num1=3 "
      end select condition3
!----------------------
    else 
      WRITE(*,*) "ILLEGAL BOUNDARY CONDITION : num1 value illegal"
    end if
  end subroutine bc

end module boundary2d
