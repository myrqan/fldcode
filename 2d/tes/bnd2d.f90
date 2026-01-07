module bnd2d
  use mpi
  use vars
  implicit none
contains
    !--------------------------------------------------
    ! surface number (snum) description
    ! (zmax) |-------------|
    !        |      2      |
    !        |             |
    !        | 3         1 |
    !        |             |
    !        |      0      |
    ! (zmin) |-------------|
    !      (xmin)        (xmax)
    !--------------------------------------------------
  subroutine bnd_2d_free(qq,snum,mg,lix,liz)
    integer,intent(in):: lix,liz,mg
    ! surface number
    integer,intent(in):: snum
    real(8),intent(inout):: qq(lix,liz)
    integer:: ii,jj
    select case(snum)
      case(0)
        if(mplz==0) then
          do jj=mg,1,-1
            qq(:,jj)=qq(:,jj+1)
          end do
        end if
      case(1)
        if(mplx==mpx-1) then
          do ii=lix-mg+1,lix
            qq(ii,:)=qq(ii-1,:)
          end do
        end if
      case(2)
        if(mplz==mpz-1) then
          do jj=liz-mg+1,liz
            qq(:,jj)=qq(:,jj-1)
          end do
        end if
      case(3)
        if(mplx==0) then
          do ii=mg,1,-1
            qq(ii,:)=qq(ii+1,:)
          end do
        end if
    end select
  end subroutine bnd_2d_free

  subroutine bnd_2d_symm(qq,snum,mg,lix,liz)
    integer,intent(in):: lix,liz,mg
    ! surface number
    integer,intent(in):: snum
    real(8),intent(inout):: qq(lix,liz)
    integer:: ii,jj
    select case(snum)
      case(0)
        if(mplz==0) then
          do jj=1,mg
            qq(:,jj)=qq(:,2*mg-jj+1)
          end do
        end if
      case(1)
        if(mplx==mpx-1) then
          do ii=1,mg
            qq(lix-mg+ii,:)=qq(lix-mg-ii+1,:)
          end do
        end if
      case(2)
        if(mplz==mpz-1) then
          do jj=1,mg
            qq(:,liz-mg+jj)=qq(:,lix-mg-jj+1)
          end do
        end if
      case(3)
        if(mplx==0) then
          do ii=1,mg
            qq(ii,:)=qq(2*mg-ii+1,:)
          end do
        end if
    end select
  end subroutine bnd_2d_symm

  subroutine bnd_2d_asym(qq,snum,mg,lix,liz)
    integer,intent(in):: lix,liz,mg
    ! surface number
    integer,intent(in):: snum
    real(8),intent(inout):: qq(lix,liz)
    integer:: ii,jj
    select case(snum)
      case(0)
        if(mplz==0) then
          do jj=1,mg
            qq(:,jj)=-qq(:,2*mg-jj+1)
          end do
        end if
      case(1)
        if(mplx==mpx-1) then
          do ii=1,mg
            qq(lix-mg+ii,:)=-qq(lix-mg-ii+1,:)
          end do
        end if
      case(2)
        if(mplz==mpz-1) then
          do jj=1,mg
            qq(:,liz-mg+jj)=-qq(:,lix-mg-jj+1)
          end do
        end if
      case(3)
        if(mplx==0) then
          do ii=1,mg
            qq(ii,:)=-qq(2*mg-ii+1,:)
          end do
        end if
    end select
  end subroutine bnd_2d_asym
end module bnd2d
