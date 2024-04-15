      
      subroutine pois
      include 'var.cmn'

      integer lbi,ldw,mxmv,info,nn_max,nn,typmv
      integer ind0,ind1,n
      parameter (nn_max=km*jm*im+1,
     $           lbi=1,ldw=nn_max*(5+2*lbi)) 
      real*8  xx(nn_max),rhs(nn_max),workbi(nn_max,5+2*lbi)
      real*8  tolref,tol,delta
      logical okprint,nonzero
 
      okprint=.false.
      nonzero=.true.
      delta=1d-2
      mxmv=5000
      tolref=1d-7

c      compute work      

!$omp parallel do private(i,j,k,ind1)
      do k=1,km
      do j=1,jm
      do i=1,im
       ind1=(k-1)*jm*im+(j-1)*im+i

       if ( typ(i,j,k).eq.1) then
        rhs(ind1)=work(i,j,k)
        xx(ind1)=p(i,j,k)
       else
        rhs(ind1)=0d0
        xx(ind1)=0d0
       end if

      end do
      end do
      end do

c     tol is changed in bicgstab2, so it should be defined as a variable
      tol=tolref
      nn=km*jm*im

      typmv=0

      call bicgstab2(okprint,lbi,nn,xx,rhs,nonzero,tol,delta,
     $              'max',mxmv,workbi,ldw,info,typmv)


       if (time.eq.0d0) then
       open(unit=620,file='p_pois_iter.dat')
       else 
       open(unit=620,file='p_pois_iter.dat',position='append')   
       end if

 620    format(e14.6,1i6,e14.6)

        write(620,620)time,mxmv,tol

       close (620)






      write(6,*)time,mxmv,tol
      if (info.ne.0) then
       if ((info.eq.1).and.(tol.lt.10d0*tolref)) then
         write(6,*)'tol slightly too large'
       else
         write(6,*)'bicg failed',info,tol,time
         stop
       end if
      end if

!$omp parallel do private(i,j,k)
      do k=1,km
      do j=1,jm
      do i=1,im
        p(i,j,k)=xx( (k-1)*jm*im+(j-1)*im+i )
      end do
      end do
      end do

      call bou_p

      return
      end


      subroutine matvec(n,xx,yy,typmv)
      include 'var.cmn'

      integer n,typmv
      integer ind1,ind2,ind3
      real*8 xx(n),yy(n)

!$omp parallel do private(i,j,k,ind1,ind2)
      do k=1,km
        ind1=(k-1)*jm*im
        do j=1,jm
         ind2=ind1+(j-1)*im
         do i=1,im
          p(i,j,k)=xx(ind2+i)
         end do
        end do
      end do
 
      call bou_p

!$omp parallel do private(i,j,k,ind1,ind2,ind3)
      do k=1,km
         ind1=(k-1)*jm*im
         do j=1,jm
          ind2=ind1+(j-1)*im
          do i=1,im
           ind3=ind2+i
           yy(ind3)=(D1(i,j,k)*p(i-1,j,k)+D2(i,j,k)*p(i+1,j,k)
     $             +D3(i,j,k)*p(i,j-1,k)+D4(i,j,k)*p(i,j+1,k)
     $             +D5(i,j,k)*p(i,j,k-1)+D6(i,j,k)*p(i,j,k+1))
     $             +p(i,j,k)
          end do
         end do
      end do    
 
      return
      end 



      subroutine sumprod(n,xx,yy,sum1)
      include 'var.cmn'
      integer n,n1,ind1,ind2,ind3
      real*8 xx(n),yy(n),sum1,parsum(km)

      n1=im*jm
!$omp parallel do private(i,k,ind1,ind2,ind3)
      do k=1,km
       parsum(k)=0d0
       ind1=(k-1)*n1
       if (k.eq.km) then
        ind2=n-ind1
       else
        ind2=n1
       end if

       do i=1,ind2
        ind3=ind1+i
        parsum(k)=parsum(k)+xx(ind3)*yy(ind3)
       end do
      end do

      sum1=0d0
      do k=1,km
       sum1=sum1+parsum(k)
      end do

      return
      end


      subroutine maxnorm(n,xx,maxval1)
      include 'var.cmn'
      integer n,n1,ind1,ind2
      real*8 xx(n),maxval1,parmax(km)

      n1=im*jm
!$omp parallel do private(i,k,ind1,ind2)
      do k=1,km
       parmax(k)=0d0
       ind1=(k-1)*n1
       if (k.eq.km) then
        ind2=n-ind1
       else
        ind2=n1
       end if
       do i=1,ind2
        parmax(k)=max(parmax(k),abs(xx(ind1+i)))
       end do
      end do

      maxval1=0d0
      do k=1,km
       maxval1=max(maxval1,parmax(k))
      end do

      return
      end


 
      subroutine sumprod0(n,xx,yy,sum1)
      include 'var.cmn'
      integer n,n1,ind1,ind2,ind3
      real*8 xx(n),yy(n),sum1,parsum(km)

      sum1=0d0
      do i=1,n
       sum1=sum1+xx(i)*yy(i)
      end do

      return
      end


      subroutine maxnorm0(n,xx,maxval1)
      include 'var.cmn'
      integer n,n1,ind1,ind2
      real*8 xx(n),maxval1,parmax(km)

      maxval1=0d0
      do i=1,n
       maxval1=max(maxval1,abs(xx(i)))
      end do

      return
      end

          
