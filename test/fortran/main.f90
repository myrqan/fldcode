program main
  implicit none

  integer,parameter :: ix=200
  real(8) :: xx(ix)
  integer,parameter :: pi = 4.d0*ATAN(1.d0)
  integer :: ii

  do ii = 1,ix
    xx(ii) = DSIN(pi*real(ii,8)/real(ix,8))
  end do

  open(unit=10,file='test.dat',access='stream',form='unformatted')
  write(10) xx
  close(10)
end program main
