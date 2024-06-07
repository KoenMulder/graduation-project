      subroutine eval_time
      include 'var.cmn'
      integer n,ilist,mini,maxi,minj,maxj,mink,maxk,
     $        i1,j1,k1
      real*8 maxuabs,u1cen,u2cen,u3cen,uabs
      real*8 maxc(3),minc(3),maxdiv,I_node(0:im,0:jm,0:km),I_ave
      real*8 max_koh,max_h2,cH_node(0:im,0:jm,0:km)
      real*8 sum1,sum2,sum3,FF1(0:jm+1,0:km+1),FF2(0:jm,0:km)


!$omp parallel do private(i,j,k)
       do k=0,km+1
            do j=0,jm+1
                  do i=0,im
                        current(i,j,k)=-conduct(i,j,k)*(phi(i+1,j,k)-phi(i,j,k))
     $                  /dx1_i(i)
                  end do
            end do
      end do

      do k=0,km
            do j=0,jm
                  do i=0,im
                        I_node(i,j,k)=0.25d0*(current(i,j,k)+current(i,j+1,k)
     $                  +current(i,j,k+1)+current(i,j+1,k+1))
                  end do
            end do
      end do

      do k=0,km
            do j=0,jm
                  do i=0,im
                        if ((x1_i(i)-xc)**2+(x2_j(j)-yc)**2+
     $                  (x3_k(k)-zc)**2.le.radius**2) then
                              I_node(i,j,k)=0d0
                        end if
                  end do
            end do
      end do

c      I_ave=0d0

c      do k=0,km
c      do j=0,jm
c      do i=0,im
c          I_ave=I_ave+I_node(i,j,k)
c      end do
c      end do
c      end do

c      I_ave=I_ave/(im+1)/(jm+1)/(km+1)

      do k=0,km
            do j=0,jm
                  I_ave=I_ave+I_node(0,j,k)
            end do
      end do

      I_ave=I_ave/(jm+1)/(km+1)

      if (time.lt.0.5d0*dtime) then
            open(unit=20,file='datt.dat')
            open(unit=211,file='c_H.dat')
c           open(unit=212,file='c_E.dat')
c           open(unit=213,file='c_W.dat')
      else 
            open(unit=20,file='datt.dat',position='append')
            open(unit=211,file='c_H.dat',position='append')
c           open(unit=212,file='c_E.dat',position='append')
c           open(unit=213,file='c_W.dat',position='append') 
      end if

      max_h2=0d0
      do k=1,km
            do j=1,jm
                  do i=1,im
                        max_h2=max(max_h2,c(i,j,k,1))
                        max_koh=max(max_koh,c(i,j,k,2))
                  end do
            end do
      end do
c     These parts are not needed
      maxuabs=-1d30
      maxdiv=-1d30
      do i=1,3
       maxc(i)=-1d30
       minc(i)=1d30
      end do

      do k=1,km
      do j=1,jm
      do i=1,im
           u1cen=0.5d0*(u1(i-1,j,k)+u1(i,j,k))
           u2cen=0.5d0*(u2(i,j-1,k)+u2(i,j,k))
           u3cen=0.5d0*(u3(i,j,k-1)+u3(i,j,k))
           uabs=sqrt(u1cen**2+u2cen**2+u3cen**2)
           if (typ(i,j,k).eq.1) then
            maxc(1)=max(maxc(1),c(i,j,k,1))
            maxc(2)=max(maxc(2),c(i,j,k,2))
            maxc(3)=max(maxc(3),c(i,j,k,3))
            minc(1)=min(minc(1),c(i,j,k,1))
            minc(2)=min(minc(2),c(i,j,k,2))
            minc(3)=min(minc(3),c(i,j,k,3))
            maxdiv=max(maxdiv,abs(div(i,j,k)))
            if (uabs.gt.maxuabs) then
                maxuabs=uabs
                i1=i
                j1=j
                k1=k

             end if


             
            if (uabs.gt.1d29) then
c             write(6,*)'uabs too large',m1,i,j,k
c             write(6,*)u1cen,u2cen,u3cen
c             stop
            end if

           end if


      end do
      end do
      end do

      do k=0,km
      do j=0,jm
      do i=0,im

        cH_node(i,j,k)=(1d0/8d0)*(c(i,j,k,1)+c(i+1,j,k,1)
     $      +c(i,j+1,k,1)+c(i+1,j+1,k,1)+c(i,j,k+1,1)
     $      +c(i+1,j,k+1,1)+c(i,j+1,k+1,1)+c(i+1,j+1,k+1,1))

      end do
      end do
      end do

      sum1=0d0
      sum2=0d0
      sum3=0d0
      

      do k=0,km+1
      do j=0,jm+1


        FF1(j,k)=(c(1,j,k,1)-c(0,j,k,1))/dx1(1)
        sum1=sum1+FF1(j,k)


      end do
      end do

      do k=0,km
      do j=0,jm


        FF2(j,k)=(cH_node(1,j,k)-cH_node(0,j,k))/dx1(1)
        sum2=sum2+FF2(j,k)
        sum3=sum3+cH_node(0,j,k)

      end do
      end do

      


      do i=1,3
       maxc(i)=maxc(i)/c_ref(i) 
       minc(i)=minc(i)/c_ref(i)
      end do
 
      mini=im
      maxi=1
      minj=jm
      maxj=1
      mink=km
      maxk=1
      do ilist=1,change-1
       mini=min(mini,x_IB(ilist))
       maxi=max(maxi,x_IB(ilist))
       minj=min(minj,y_IB(ilist))
       maxj=max(maxj,y_IB(ilist))
       mink=min(mink,z_IB(ilist))
       maxk=max(maxk,z_IB(ilist))
      end do

      write(6,60)time,m1,change-1,cube_number,
     $           mini,maxi,minj,maxj,mink,maxk,
     $           radius,height,total,dxcdt,drdt,maxuabs,
     $           i1,j1,k1,maxdiv,
     $           (minc(i),i=1,3),(maxc(i),i=1,3),
     $           I_ave,max_h2,max_koh,sum1/(jm+2)/(km+2),
     $           sum2/(jm+1)/(km+1),sum3/(jm+1)/(km+1)
     
      write(20,60)time,m1,change-1,cube_number,
     $            mini,maxi,minj,maxj,mink,maxk,
     $            radius,height,total,dxcdt,drdt,maxuabs,
     $            i1,j1,k1,maxdiv,
     $            (minc(i),i=1,3),(maxc(i),i=1,3),
     $            I_ave,max_h2,max_koh,sum1/(jm+2)/(km+2),
     $           sum2/(jm+1)/(km+1),sum3/(jm+1)/(km+1)



 60   format(e14.6,3i9,6i6,6e14.6,3i4,100e14.6)
      close(20)

      if (maxuabs.gt.2d0) then
       write(6,*)'maxuabs too large'
       stop
      end if

      if (1d0.eq.1d0) then


      do k=0,km+1
      do j=0,jm+1
      do i=0,im+1
c       write(211,60)(c(i,j,k,1))
c       write(212,60)(c(i,j,k,2))
c       write(213,60)(c(i,j,k,3))
      end do
      end do
      end do



      do j=0,jm
      do i=0,im
       write(211,60)(cH_node(i,j,km/2))
      end do
      end do


      do j=0,jm
      do i=0,im
       write(211,60)(cH_node(i,j,km/2))
      end do
      end do

      end if
      close (211)
c      close (212)
c      close (213)


      if (total.gt.1d-10) then
       write(6,*)'total is positive'
       write(6,*)'dtime=',dtime
       stop
      end if
      if (radius.gt.R_max) then
       write(6,*)'radius 0.4 lx1 reached'
       write(6,*)'beta range'
       open(unit=61,file='betarange.dat')
       do i=1,36
 61     format(i5,10e14.6)
        write(6,61)i,minbeta(i),maxbeta(i)
        write(61,61)i,minbeta(i),maxbeta(i)
       end do
       close(61)

       open(unit=111,file='datf',form='unformatted')
       write(111)time,p,u1,u2,u3,c,eta,phi,total,radius,
     $          height,xc,ka,kc,phiL,radius0,conduct,
     $          typ,typ1,typ2,typ3,drdt
       close(111) 

      call eval_prof

       stop
      end if

      return
      end 


      subroutine eval_prof
      include 'var.cmn'
      integer n,o


      do k=0,km
      do j=0,jm
      do i=0,im

        u_node(i,j,k)=0.25d0*(u1(i,j,k)+u1(i,j+1,k)
     $      +u1(i,j,k+1)+u1(i,j+1,k+1))
        v_node(i,j,k)=0.25d0*(u2(i,j,k)+u2(i+1,j,k)
     $      +u2(i,j,k+1)+u2(i+1,j,k+1))
        w_node(i,j,k)=0.25d0*(u3(i,j,k)+u3(i+1,j,k)
     $      +u3(i,j+1,k)+u3(i+1,j+1,k))
        p_node(i,j,k)=(1d0/8d0)*(p(i,j,k)+p(i+1,j,k)
     $      +p(i,j+1,k)+p(i+1,j+1,k)+p(i,j,k+1)
     $      +p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1))

      end do
      end do
      end do




 11   format(1000e36.24)
      open(unit=11,file='xfile_NV.dat')
      open(unit=21,file='xfile.dat')
      open(unit=22,file='dc_current.dat')
c      open(unit=23,file='overpotential.dat')
      open(unit=51,file='radial_u.dat')
      open(unit=52,file='radial_v.dat')
      open(unit=53,file='radial_w.dat')
      open(unit=54,file='KOH.dat')


      do o=1,o_u
        write(51,11)thetaS_U(o),phiS_U(o),
     $  un_surf_u(o),ut_surf_u(o),us_surf_u(o),
     $  dc_t_u(o),dc_s_u(o)
      end do

      do o=1,o_v
        write(52,11)thetaS_V(o),phiS_V(o),
     $  un_surf_v(o),ut_surf_v(o),us_surf_v(o),
     $  dc_t_v(o),dc_s_v(o)
      end do

      do o=1,o_w
        write(52,11)thetaS_W(o),phiS_W(o),
     $  un_surf_w(o),ut_surf_w(o),us_surf_w(o),
     $  dc_t_w(o),dc_s_w(o)
      end do

      do o=1,change-1
        write(54,11)thetaS(o),phiS(o),
     $  c_KOH_surf(o)
      end do

      close(51)
      close(52)
      close(53)
      close(54)



      do k=0,km
      do j=0,jm
      do i=0,im

       write(11,11)x1_i(i),x2_j(j),x3_k(k),u_node(i,j,k),
     $      v_node(i,j,k),w_node(i,j,k),p_node(i,j,k)

      end do
      end do
      end do

      do k=0,km+1
      do j=0,jm+1
      do i=0,im+1

       write(21,11)x1(i),x2(j),x3(k),(c(i,j,k,n),n=1,nmax),phi(i,j,k)

      end do
      end do
      end do

      do k=0,km+1
      do j=0,jm+1

c      write(23,11)x2(j),x3(k),eta(j,k)
                     
      end do
      end do

      do k=0,km+1
      do j=0,jm+1                  
      do i=0,im

       write(22,11)current(i,j,k),
     $             ((c(i+1,j,k,n)-c(i,j,k,n))/dx1_i(i),n=1,nmax)

      end do
      end do
      end do

      close(11)
      close(21)
      close(22)
c      close(23)
            
      return 
      end


