module sts
  implicit none
  contains
    subroutine determine_sts_s(ss,dx,dt,eta)
      integer,intent(inout) :: ss
      real(8),intent(in) :: dx,eta,dt

      integer :: trys
      real(8) :: t
      

      !t = dt / (dx**2/(2.d0*eta))
      t = 2.d0 * eta * dt / dx**2

      trys = ceiling(0.5d0 * (-1.d0 + sqrt(16.d0*t+9.d0)))

      if(mod(trys,2) == 1) then
        ss = trys
      else
        ss = trys + 1
      end if
    end subroutine determine_sts_s

    subroutine sts_update(qq,dt,dx,gx,eta,ss_max)
      integer,intent(in) :: gx
      real(8),intent(in) :: dt,dx,eta
      real(8),intent(inout) :: qq(gx)
      integer,intent(inout) :: ss_max


      integer :: ss
      real(8),allocatable,dimension(:) :: aa,bb
      real(8),allocatable,dimension(:) :: mu, mut, nu, gmt
      real(8) :: w1
      !real(8) :: w1, mut(0:ss), mu(ss), nu(ss), gmt(ss)

      real(8),allocatable,dimension(:,:) :: qqn
      real(8) :: ly0(gx), lyj1(gx)

      integer :: ii, jj
      real(8) :: jr,sr


      !===============
      ! determine sts steps
      !===============

      call determine_sts_s(ss,dx,dt,eta)
      ss_max = max(ss_max, ss)
      sr = real(ss,8)

      allocate(aa(0:ss),bb(0:ss))
      allocate(mut(ss),mu(ss),nu(ss),gmt(ss))
      allocate(qqn(0:ss,gx))

      w1 = 4.d0 / (sr**2 + sr - 2.d0)
      bb(0) = 1.d0/3.d0
      bb(1) = bb(0)
      aa(0) = 1.d0 - bb(0)
      aa(1) = 1.d0 - bb(1)
      do jj = 2, ss
        jr = real(jj,8)
        bb(jj) = (jr**2+jr-2.d0)/(2.d0*jr*(jr+1.d0))
        aa(jj) = 1.d0 - bb(jj)
      end do

      mut(1) = w1/3.d0
      do jj = 2, ss
        jr = real(jj,8)
        mu(jj) = (2.d0*jr-1.d0)/jr * bb(jj) / bb(jj-1)
        nu(jj) = - (jr-1.d0)/(jr) * bb(jj) / bb(jj-2)
        mut(jj) = w1*mu(jj)
        gmt(jj) = -aa(jj-1) * mut(jj)
      end do

      !==================
      ! sts main update
      !==================
      qqn(0,:) = qq(:)
      call lap(ly0,qqn(0,:),dx,gx,eta)
      qqn(1,:) = qqn(0,:) + mut(1)*dt*ly0

      do jj = 2, ss
        call lap(lyj1,qqn(jj-1,:),dx,gx,eta)
        qqn(jj,:) = mu(jj)*qqn(jj-1,:) + nu(jj)*qqn(jj-2,:) &
          & + (1 - mu(jj) - nu(jj)) * qqn(0,:) + mut(jj) * dt * lyj1(:)&
          & + gmt(jj) * dt * ly0(:)
      end do
      
      qq(:) = qqn(ss,:)
    end subroutine sts_update

    subroutine lap(lqq,qq,dx,gx,eta)
      integer,intent(in) :: gx
      real(8),intent(inout) :: lqq(gx)
      real(8),intent(in) :: qq(gx),dx,eta

      integer :: i
      do i = 2, gx-1
        lqq(i) = eta/(dx**2) * (qq(i+1) - 2 * qq(i) + qq(i-1))
      end do

      lqq(1) = 0.d0
      lqq(gx) = 0.d0
    end subroutine lap
end module sts
