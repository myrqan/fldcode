module file_output
  implicit none
contains
  subroutine reset(ix)
    integer,intent(in) :: ix
    ! ix : grid size
    call system("rm t.dac x.dac rho.dac vx.dac p.dac eps.dac param.txt")!
    call system("rm t.txt x.txt rho.txt vx.txt p.txt eps.txt")!
    open(50, file='param.txt',form='formatted',position='append')
    write(50, *) "ix:", ix
    close(50)
  end subroutine reset

  subroutine time_param(t)
    double precision,intent(in) :: t
    ! t : time
    open(50, file='param.txt',form='formatted',position='append')
    write(50, *) "t:", t
    close(50)
  end subroutine time_param

  subroutine qv_param(qv)
    double precision,intent(in) :: qv
    open(50, file='param.txt',form='formatted',position='append')
    write(50, *) "qv:", qv
    close(50)
  end subroutine qv_param

  subroutine put0dreal(device_num,filename,sc)
    integer,intent(in) :: device_num
    character(len=*),intent(in) :: filename
    double precision,intent(in) :: sc
    ! device num : device number 
    !   use a number (>= 10)
    ! filename : name for data file
    ! sc : scalar (real number)
    open(device_num,file=filename,form='unformatted',position='append')
    write(device_num) sc
    close(device_num)
  end subroutine put0dreal

  subroutine put0dreal_txt(device_num,filename,sc)
    integer,intent(in) :: device_num
    character(len=*),intent(in) :: filename
    double precision,intent(in) :: sc
    ! device num : device number 
    !   use a number (>= 10)
    ! filename : name for data file
    ! sc : scalar (real number)
    open(device_num,file=filename,form='formatted',position='append')
    write(device_num, *) sc
    write(device_num, *)
    close(device_num)
  end subroutine put0dreal_txt


  subroutine put0dint(device_num,filename,sc)
    integer,intent(in) :: device_num
    character(len=*),intent(in) :: filename
    integer,intent(in) :: sc
    ! device num : device number 
    !   use a number (>= 10)
    ! filename : name for data file
    ! sc : scalar (integer)
    open(device_num,file=filename,form='unformatted',position='append')
    write(device_num) sc
    close(device_num)
  end subroutine put0dint

  subroutine put0dint_txt(device_num,filename,sc)
    integer,intent(in) :: device_num
    character(len=*),intent(in) :: filename
    integer,intent(in) :: sc
    ! device num : device number 
    !   use a number (>= 10)
    ! filename : name for data file
    ! sc : scalar (integer)
    open(device_num,file=filename,form='formatted',position='append')
    write(device_num, *) sc
    write(device_num, *)
    close(device_num)
  end subroutine put0dint_txt


  subroutine put1dreal(device_num,filename,ary)
    integer,intent(in) :: device_num
    character(len=*),intent(in) :: filename
    double precision,intent(in) :: ary(:)
    ! device num : device number 
    !   use a number (>= 10)
    ! filename : name for data file
    ! ary : array (double precision)
    open(device_num,file=filename,form='unformatted',position='append')
    write(device_num) ary
    close(device_num)
  end subroutine put1dreal

  subroutine put1dreal_txt(device_num,filename,ary)
    integer,intent(in) :: device_num
    character(len=*),intent(in) :: filename
    double precision,intent(in) :: ary(:)
    ! device num : device number 
    !   use a number (>= 10)
    ! filename : name for data file
    ! ary : array (double precision)
    open(device_num,file=filename,form='formatted',position='append')
    write(device_num, *) ary
    write(device_num, *)
    close(device_num)
  end subroutine put1dreal_txt

  
end module file_output
