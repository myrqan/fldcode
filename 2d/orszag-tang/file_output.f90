module file_output
  implicit none
contains
  subroutine clean()
    CALL system("rm *.dat param.txt")
  end subroutine clean

  ! #############################

  subroutine put_param_int(char,num)
    CHARACTER(*),INTENT(IN) :: char
    INTEGER,INTENT(IN) :: num
    INTEGER :: unit_num

    OPEN(newunit=unit_num,&
      file='param.txt',&
      form='formatted',&
      position='append')
    WRITE(unit_num,*) char, num
    CLOSE(unit_num)
  end subroutine put_param_int

  ! #############################

  subroutine put_param_dble(char,num)
    CHARACTER(*),INTENT(IN) :: char
    DOUBLE PRECISION,INTENT(IN) :: num
    INTEGER :: unit_num

    OPEN(newunit=unit_num,&
      file='param.txt',&
      form='formatted',&
      position='append')
    WRITE(unit_num,*) char, num
    CLOSE(unit_num)
  end subroutine put_param_dble

  ! #############################

  subroutine write_0d_dble(char,val)
    CHARACTER(*),INTENT(IN) :: char
    DOUBLE PRECISION,INTENT(IN) :: val
    INTEGER :: unit_num

    OPEN(newunit=unit_num,&
      file=char,&
      form='unformatted',&
      position='append')
    WRITE(unit_num) val
    CLOSE(unit_num)

  end subroutine write_0d_dble

  ! #############################
  
  subroutine write_1d_int(char,val)
    CHARACTER(*),INTENT(IN) :: char
    INTEGER,INTENT(IN) :: val(:)
    INTEGER :: unit_num

    OPEN(newunit=unit_num,&
      file=char,&
      form='unformatted',&
      position='append')
    WRITE(unit_num) val
    CLOSE(unit_num)

  end subroutine write_1d_int
  ! #############################

  subroutine write_2d_dble(char,val)
    CHARACTER(*),INTENT(IN) :: char
    DOUBLE PRECISION,INTENT(IN) :: val(:,:)
    INTEGER :: unit_num

    OPEN(newunit=unit_num,&
      file=char,&
      form='unformatted',&
      position='append')
    WRITE(unit_num) val
    CLOSE(unit_num)

  end subroutine write_2d_dble

  ! #############################

end module file_output
