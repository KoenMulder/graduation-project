      subroutine u 
      include 'var.cmn'
      real*8  tp,ua,va,wa

      
!$omp parallel do private(ua,va,wa,tp,i,j,k)
       do k=1,km
       do j=1,jm
       do i=1,im-1

        ua=0.5d0*(u1(i,j,k)+u1(i-1,j,k))
        tp=-ua**2/dx1_i(i)
        ua=0.5d0*(u1(i+1,j,k)+u1(i,j,k))
        tp=tp+ua**2/dx1_i(i)

        va=0.5d0*(u2(i+1,j-1,k)+u2(i,j-1,k))
        tp=tp-0.5d0*(u1(i,j-1,k)+u1(i,j,k))*va/dx2(j)
        va=0.5d0*(u2(i+1,j,k)+u2(i,j,k))
        tp=tp+0.5d0*(u1(i,j,k)+u1(i,j+1,k))*va/dx2(j)

        wa=0.5d0*(u3(i+1,j,k-1)+u3(i,j,k-1))
        tp=tp-0.5d0*(u1(i,j,k-1)+u1(i,j,k))*wa/dx3(k)
        wa=0.5d0*(u3(i+1,j,k)+u3(i,j,k))
        tp=tp+0.5d0*(u1(i,j,k)+u1(i,j,k+1))*wa/dx3(k)

        rhu(i,j,k)=u1(i,j,k)-tp*dtime
     $  	-nu*dtime/dx1_i(i)*((u1(i,j,k)-u1(i-1,j,k))/dx1(i)
     $          -(u1(i+1,j,k)-u1(i,j,k))/dx1(i+1))
     $	        -nu*dtime/dx2(j)*((u1(i,j,k)-u1(i,j-1,k))/dx2_j(j-1)
     $          -(u1(i,j+1,k)-u1(i,j,k))/dx2_j(j))
     $	        -nu*dtime/dx3(k)*((u1(i,j,k)-u1(i,j,k-1))/dx3_k(k-1)
     $          -(u1(i,j,k+1)-u1(i,j,k))/dx3_k(k))


       end do
       end do
       end do

 
      return
      end

      subroutine bou_rhu
      include 'var.cmn' 
       
      do k=1,km
      do j=1,jm
      do i=1,im-1
      if (typ1(i,j,k).eq.1) then

       u1(i,j,k)=rhu(i,j,k)

      end if
      end do
      end do
      end do

      do k=1,km
      do j=1,jm

        u1(im,j,k)=0d0
        u1(0,j,k)=0d0

      end do
      end do


      return
      end





      subroutine bou_u
      include 'var.cmn' 
       
c     Simple Dirichlet boundary at left and right side
c     Second order

      u1(0,1:jm,1:km)=0d0
      u1(im,1:jm,1:km)=0d0  


c     Simple Dirichlet boundary at inflow
c     Second order

      u1(0:im,0,1:km)=-u1(0:im,1,1:km)


c     Zero Neumann conditions at outflow boundary

      u1(0:im,jm+1,1:km)=u1(0:im,jm,1:km)


c     Periodic boundary conditions

      u1(0:im,0:jm+1,0)=u1(0:im,0:jm+1,km)
      u1(0:im,0:jm+1,km+1)=u1(0:im,0:jm+1,1)


      return
      end




      subroutine v 
      include 'var.cmn'
      real*8  tp,ua,va,wa
      


!$omp parallel do private(ua,va,wa,tp,i,j,k)
       do k=1,km
       do j=1,jm
       do i=1,im

        ua=0.5d0*(u1(i-1,j,k)+u1(i-1,j+1,k))
        tp=-0.5d0*(u2(i-1,j,k)+u2(i,j,k))*ua/dx1(i)
        ua=0.5d0*(u1(i,j,k)+u1(i,j+1,k))
        tp=tp+0.5d0*(u2(i,j,k)+u2(i+1,j,k))*ua/dx1(i)

        va=0.5d0*(u2(i,j,k)+u2(i,j-1,k))
        tp=tp-va**2/dx2_j(j)
        va=0.5d0*(u2(i,j,k)+u2(i,j+1,k))
        tp=tp+va**2/dx2_j(j)

        wa=0.5d0*(u3(i,j,k-1)+u3(i,j+1,k-1))
        tp=tp-0.5d0*(u2(i,j,k-1)+u2(i,j,k))*wa/dx3(k)
        wa=0.5d0*(u3(i,j,k)+u3(i,j+1,k))
        tp=tp+0.5d0*(u2(i,j,k)+u2(i,j,k+1))*wa/dx3(k)

        rhv(i,j,k)=u2(i,j,k)-tp*dtime
     $  	-nu*dtime/dx1(i)*((u2(i,j,k)-u2(i-1,j,k))/dx1_i(i-1)
     $          -(u2(i+1,j,k)-u2(i,j,k))/dx1_i(i))
     $	        -nu*dtime/dx2_j(j)*((u2(i,j,k)-u2(i,j-1,k))/dx2(j)
     $          -(u2(i,j+1,k)-u2(i,j,k))/dx2(j+1))
     $	        -nu*dtime/dx3(k)*((u2(i,j,k)-u2(i,j,k-1))/dx3_k(k-1)
     $          -(u2(i,j,k+1)-u2(i,j,k))/dx3_k(k))


       end do
       end do
       end do

 
      return
      end

      subroutine bou_rhv
      include 'var.cmn' 
       
      do k=1,km
      do j=1,jm
      do i=1,im

      if (typ2(i,j,k).eq.1) then

       u2(i,j,k)=rhv(i,j,k)

      end if

      end do
      end do
      end do

       do k=1,km
       do i=1,im

          u2(i,0,k)=4d0*umax*(x1(i)/lx1)*(1d0-(x1(i)/lx1))

       end do
       end do


      return
      end


      subroutine bou_v
      include 'var.cmn' 


c      Simple Dirichlet boundary at inflow
c      Second order

      do k=1,km
      do i=1,im

       u2(i,0,k)=4d0*umax*(x1(i)/lx1)*(1d0-(x1(i)/lx1))

      end do
      end do

      u2(1:im,jm+1,1:km)=u2(1:im,jm,1:km)



c     Simple Dirichlet boundary at left and right side
c     Second order

      u2(0,0:jm+1,1:km)=-u2(1,0:jm+1,1:km)
      u2(im+1,0:jm+1,1:km)=-u2(im,0:jm+1,1:km)  



c     Periodic boundary conditions

      u2(0:im+1,0:jm+1,0)=u2(0:im+1,0:jm+1,km)
      u2(0:im+1,0:jm+1,km+1)=u2(0:im+1,0:jm+1,1)

  
      return
      end



      subroutine w    
      include 'var.cmn'
      real*8  tp,ua,va,wa
      
!$omp parallel do private(ua,va,wa,tp,i,j,k)
       do k=1,km
       do j=1,jm
       do i=1,im

        ua=0.5d0*(u1(i-1,j,k)+u1(i-1,j,k+1))
        tp=-0.5d0*(u3(i-1,j,k)+u3(i,j,k))*ua/dx1(i)
        ua=0.5d0*(u1(i,j,k)+u1(i,j,k+1))
        tp=tp+0.5d0*(u3(i,j,k)+u3(i+1,j,k))*ua/dx1(i)

        va=0.5d0*(u2(i,j-1,k)+u2(i,j-1,k+1))
        tp=tp-0.5d0*(u3(i,j-1,k)+u3(i,j,k))*va/dx2(j)
        va=0.5d0*(u2(i,j,k)+u2(i,j,k+1))
        tp=tp+0.5d0*(u3(i,j,k)+u3(i,j+1,k))*va/dx2(j)

        wa=0.5d0*(u3(i,j,k)+u3(i,j,k-1))
        tp=tp-wa**2/dx3_k(k)
        wa=0.5d0*(u3(i,j,k+1)+u3(i,j,k))
        tp=tp+wa**2/dx3_k(k)

        rhw(i,j,k)=u3(i,j,k)-tp*dtime
     $  	-nu*dtime/dx1(i)*((u3(i,j,k)-u3(i-1,j,k))/dx1_i(i-1)
     $          -(u3(i+1,j,k)-u3(i,j,k))/dx1_i(i))
     $	        -nu*dtime/dx2(j)*((u3(i,j,k)-u3(i,j-1,k))/dx2_j(j-1)
     $          -(u3(i,j+1,k)-u3(i,j,k))/dx2_j(j))
     $	        -nu*dtime/dx3(k)*((u3(i,j,k)-u3(i,j,k-1))/dx3(k)
     $          -(u3(i,j,k+1)-u3(i,j,k))/dx3(k+1))

       end do
       end do
       end do



      return
      end



      subroutine bou_rhw
      include 'var.cmn' 
       
      do k=1,km
      do j=1,jm
      do i=1,im

       if(typ3(i,j,k).eq.1) then

       u3(i,j,k)=rhw(i,j,k)

       end if

      end do
      end do
      end do

       do j=1,jm
       do i=1,im

        u3(i,j,0)=rhw(i,j,km)

       end do
       end do


      return
      end


      subroutine bou_w
      include 'var.cmn'

       
c     Simple Dirichlet boundary at left and right side
c     Second order

      u3(0,1:jm,1:km)=-u3(1,1:jm,1:km)
      u3(im+1,1:jm,1:km)=-u3(im,1:jm,1:km) 


c     Simple Dirichlet boundary at inflow
c     Second order

      u3(0:im+1,0,1:km)=-u3(0:im+1,1,1:km)


c     Zero Neumann conditions at outflow boundary

      u3(0:im+1,jm+1,1:km)=u3(0:im+1,jm,1:km)


c     Periodic boundary conditions

      u3(0:im+1,0:jm+1,0)=u3(0:im+1,0:jm+1,km)
      u3(0:im+1,0:jm+1,km+1)=u3(0:im+1,0:jm+1,1)

 
      return
      end


      subroutine pre_p    
      include 'var.cmn'
      integer o

      do k=1,km
      do j=1,jm
      do i=1,im
       D1(i,j,k)=1d0/dx1(i)/dx1_i(i-1)
       D2(i,j,k)=1d0/dx1(i)/dx1_i(i)
       D3(i,j,k)=1d0/dx2(j)/dx2_j(j-1) 
       D4(i,j,k)=1d0/dx2(j)/dx2_j(j) 
       D5(i,j,k)=1d0/dx3(k)/dx3_k(k-1)
       D6(i,j,k)=1d0/dx3(k)/dx3_k(k)

      end do
      end do
      end do


      do o=1,change-1
       i=x_IB(o)
       j=y_IB(o)
       k=z_IB(o)
  
       if (neighbor(o,1).eq.1) then
          D1(i+1,j,k)=0d0
       end if

       if (neighbor(o,2).eq.1) then
          D2(i-1,j,k)=0d0
       end if

       if (neighbor(o,3).eq.1) then
          D3(i,j+1,k)=0d0
       end if

       if (neighbor(o,4).eq.1) then
          D4(i,j-1,k)=0d0
       end if

       if (neighbor(o,5).eq.1) then
          D5(i,j,k+1)=0d0
       end if

       if (neighbor(o,6).eq.1) then
          D6(i,j,k-1)=0d0
       end if

      end do




!$omp parallel do
      do k=1,km
      do j=1,jm
      do i=1,im




       D(i,j,k)=-(D1(i,j,k)+D2(i,j,k)+D3(i,j,k)+
     $          D4(i,j,k)+D5(i,j,k)+D6(i,j,k))



       work(i,j,k)=ro/dtime/D(i,j,k)*
     $            ((u1(i,j,k)-u1(i-1,j,k))/dx1(i)+
     $            (u2(i,j,k)-u2(i,j-1,k))/dx2(j)+
     $            (u3(i,j,k)-u3(i,j,k-1))/dx3(k))



      D1(i,j,k)=D1(i,j,k)/D(i,j,k)
      D2(i,j,k)=D2(i,j,k)/D(i,j,k)
      D3(i,j,k)=D3(i,j,k)/D(i,j,k)
      D4(i,j,k)=D4(i,j,k)/D(i,j,k)
      D5(i,j,k)=D5(i,j,k)/D(i,j,k)
      D6(i,j,k)=D6(i,j,k)/D(i,j,k)








      end do
      end do
      end do




      return
      end




      subroutine pressure    
      include 'var.cmn'
      real*8 relax

      relax=1.5d0
c     Jacobi iteration step
c     Relaxation factor 1 
      do k=1,km
      do j=1,jm
      do i=1,im
      if (typ(i,j,k).eq.1) then
       p(i,j,k)=p(i,j,k)*(1-relax)+relax*(
     $          -D1(i,j,k)*p(i-1,j,k)
     $          -D2(i,j,k)*p(i+1,j,k)
     $          -D3(i,j,k)*p(i,j-1,k)
     $          -D4(i,j,k)*p(i,j+1,k)
     $          -D5(i,j,k)*p(i,j,k-1)
     $          -D6(i,j,k)*p(i,j,k+1)+
     $          work(i,j,k))
      end if
      end do
      end do
      end do

      return
      end


      subroutine velocity
      include 'var.cmn' 

!$omp parallel do
      do k=1,km
      do j=1,jm
      do i=1,im

       if(typ1(i,j,k).eq.1) then
        u1(i,j,k)=u1(i,j,k)-dtime/ro/dx1_i(i)*(
     $          p(i+1,j,k)-p(i,j,k))
       end if

       if(typ2(i,j,k).eq.1) then
        u2(i,j,k)=u2(i,j,k)-dtime/ro/dx2_j(j)*(
     $          p(i,j+1,k)-p(i,j,k))
       end if

       if(typ3(i,j,k).eq.1) then
        u3(i,j,k)=u3(i,j,k)-dtime/ro/dx3_k(k)*(
     $          p(i,j,k+1)-p(i,j,k))
       end if

      end do
      end do
      end do

      return
      end

 

      subroutine bou_p
      include 'var.cmn' 
      
       
c     Simple Dirichlet boundary at left and right side
c     Second order

      p(0,1:jm,1:km)=p(1,1:jm,1:km)
      p(im+1,1:jm,1:km)=p(im,1:jm,1:km) 


c     Simple Dirichlet boundary at inflow
c     Second order

      p(0:im+1,0,1:km)=p(0:im+1,1,1:km)
      p(0:im+1,jm+1,1:km)=-p(0:im+1,jm,1:km)


c     Periodic boundary conditions

      p(0:im+1,0:jm+1,0)=p(0:im+1,0:jm+1,km)
      p(0:im+1,0:jm+1,km+1)=p(0:im+1,0:jm+1,1)

      return
      end


      subroutine beta_function(icall)
      include 'var.cmn' 
      integer icall

      minbeta(icall)=min(beta,minbeta(icall))
      maxbeta(icall)=max(beta,maxbeta(icall))

      if (mod(icall,5).eq.1.and.icall.le.30)  then
       if ((beta.lt.-0.5d0-1d-12).or.(beta.gt.1.0d0+1d-12)) then
        write(6,*)'beta out of range',beta,'icall',icall
        write(6,*)'i,j,k',i,j,k
        stop
       end if
      end if
      if (icall.ge.31) then
       if ((beta.lt.-1d-12).or.(beta.gt.1d0+1d-12)) then
        write(6,*)'beta out of range',beta,'icall',icall
        write(6,*)'i,j,k',i,j,k
        stop
       end if
      end if

      if (beta .gt. 0.5d0) then

       alfa_s=8d0/3d0
       beta_s=2d0-(2d0-beta)*alfa_s 
       gama_s=-1d0+(1d0-beta)*alfa_s
c       beta_s=(-10d0+8d0*beta)/3d0
c       gama_s=(5d0-8d0*beta)/3d0
      else

       alfa_s=2d0/(1d0-beta)/(2d0-beta)
       beta_s=-2d0*beta/(1d0-beta)
       gama_s=beta/(2d0-beta)

      end if

      return
      end

      subroutine beta_function_linear(icall)
      include 'var.cmn' 
      integer icall

      minbeta(icall)=min(beta,minbeta(icall))
      maxbeta(icall)=max(beta,maxbeta(icall))

      if (mod(icall,5).eq.1.and.icall.le.30)  then
       if ((beta.lt.-0.5d0-1d-12).or.(beta.gt.1.0d0+1d-12)) then
        write(6,*)'beta out of range',beta,'icall',icall
        write(6,*)'i,j,k',i,j,k
        stop
       end if
      end if
      if (icall.ge.31) then
       if ((beta.lt.-1d-12).or.(beta.gt.1d0+1d-12)) then
        write(6,*)'beta out of range',beta,'icall',icall
        write(6,*)'i,j,k',i,j,k
        stop
       end if
      end if
      if (beta .gt. 0.5d0) then

       alfa_s=2d0
       beta_s=2d0-(2d0-beta)*alfa_s
       gama_s=-1d0+(1d0-beta)*alfa_s
c       beta_s=-2d0+2d0*beta
c       gama_s=1d0-2d0*beta

      else

       alfa_s=1d0/(1d0-beta)
       beta_s=-beta/(1d0-beta)
       gama_s=0d0

      end if

      return
      end



      subroutine z_plus
      include 'var.cmn'
      real*8  z,u1s,u2s,u3s

      z=zc+sqrt(radius**2-((x1(i)-xc)**2 + (x2(j)-yc)**2))
      beta=(z-x3_k(k))/dx3(k)
      call beta_function(1)
      u3s=drdt*(z-zc)/radius
      u3(i,j,k)=alfa_s*u3s+beta_s*u3(i,j,k+1)+gama_s*u3(i,j,k+2)

      z=zc+sqrt(radius**2-((x1_i(i)-xc)**2 + (x2(j)-yc)**2))
      beta=(z-x3(k))/dx3_k(k)
      call beta_function(2)
      u1s=dxcdt+drdt*(x1_i(i)-xc)/radius
      u1(i,j,k)=alfa_s*u1s+beta_s*u1(i,j,k+1)+gama_s*u1(i,j,k+2)

c     Exclude i=1 for the following point (x1_i(i-1)=0 (on the wall) and  
c     the argument of sqrt is negative)

      if (i.ge.2) then
       z=zc+sqrt(radius**2-((x1_i(i-1)-xc)**2 + (x2(j)-yc)**2))
       beta=(z-x3(k))/dx3_k(k)
       call beta_function(3)
       u1s=dxcdt+drdt*(x1_i(i-1)-xc)/radius
       u1(i-1,j,k)=alfa_s*u1s+beta_s*u1(i-1,j,k+1)+gama_s*u1(i-1,j,k+2)
      end if

      z=zc+sqrt(radius**2-((x1(i)-xc)**2 + (x2_j(j)-yc)**2))
      beta=(z-x3(k))/dx3_k(k)
      call beta_function(4)
      u2s=drdt*(x2_j(j)-yc)/radius
      u2(i,j,k)=alfa_s*u2s+beta_s*u2(i,j,k+1)+gama_s*u2(i,j,k+2)

      z=zc+sqrt(radius**2-((x1(i)-xc)**2 + (x2_j(j-1)-yc)**2))
      beta=(z-x3(k))/dx3_k(k)
      call beta_function(5)
      u2s=drdt*(x2_j(j-1)-yc)/radius
      u2(i,j-1,k)=alfa_s*u2s+beta_s*u2(i,j-1,k+1)+gama_s*u2(i,j-1,k+2)

      return
      end


      subroutine z_minus
      include 'var.cmn'
      real*8  z,u1s,u2s,u3s

      z=zc-sqrt(radius**2-((x1(i)-xc)**2 + (x2(j)-yc)**2))
      beta=-(z-x3_k(k-1))/dx3(k)
      call beta_function(6)
      u3s=drdt*(z-zc)/radius
      u3(i,j,k-1)=alfa_s*u3s+beta_s*u3(i,j,k-2)+gama_s*u3(i,j,k-3)

      z=zc-sqrt(radius**2-((x1_i(i)-xc)**2 + (x2(j)-yc)**2))
      beta=-(z-x3(k))/dx3_k(k)
      call beta_function(7)
      u1s=dxcdt+drdt*(x1_i(i)-xc)/radius
      u1(i,j,k)=alfa_s*u1s+beta_s*u1(i,j,k-1)+gama_s*u1(i,j,k-2)

c     Exclude i=1 for the following point (x1_i(i-1)=0 (on the wall) and  
c     the argument of sqrt is negative)

      if (i.ge.2) then
       z=zc-sqrt(radius**2-((x1_i(i-1)-xc)**2 + (x2(j)-yc)**2))
       beta=-(z-x3(k))/dx3_k(k)
       u1s=dxcdt+drdt*(x1_i(i-1)-xc)/radius
       call beta_function(8)
       u1(i-1,j,k)=alfa_s*u1s+beta_s*u1(i-1,j,k-1)+gama_s*u1(i-1,j,k-2)
      end if

      z=zc-sqrt(radius**2-((x1(i)-xc)**2 + (x2_j(j)-yc)**2))
      beta=-(z-x3(k))/dx3_k(k)
      call beta_function(9)
      u2s=drdt*(x2_j(j)-yc)/radius
      u2(i,j,k)=alfa_s*u2s+beta_s*u2(i,j,k-1)+gama_s*u2(i,j,k-2)

      z=zc-sqrt(radius**2-((x1(i)-xc)**2 + (x2_j(j-1)-yc)**2))
      beta=-(z-x3(k))/dx3_k(k)
      call beta_function(10)
      u2s=drdt*(x2_j(j-1)-yc)/radius
      u2(i,j-1,k)=alfa_s*u2s+beta_s*u2(i,j-1,k-1)+gama_s*u2(i,j-1,k-2)

      return
      end


      subroutine y_plus
      include 'var.cmn'
      real*8  y,u1s,u2s,u3s

      y=yc+sqrt(radius**2-( (x1(i)-xc)**2 + (x3(k)-zc)**2 ))
      beta=(y-x2_j(j))/dx2(j)
      call beta_function(11)
      u2s=drdt*(y-yc)/radius
      u2(i,j,k)=alfa_s*u2s+beta_s*u2(i,j+1,k)+gama_s*u2(i,j+2,k)

      y=yc+sqrt(radius**2-( (x1_i(i)-xc)**2 + (x3(k)-zc)**2 ))
      beta=(y-x2(j))/dx2_j(j)
      call beta_function(12)
      u1s=dxcdt+drdt*(x1_i(i)-xc)/radius
      u1(i,j,k)=alfa_s*u1s+beta_s*u1(i,j+1,k)+gama_s*u1(i,j+2,k)

c     Exclude i=1 for the following point (x1_i(i-1)=0 (on the wall) and  
c     the argument of sqrt is negative)

      if (i.ge.2) then
       y=yc+sqrt(radius**2-( (x1_i(i-1)-xc)**2 + (x3(k)-zc)**2 ))
       beta=(y-x2(j))/dx2_j(j)
       call beta_function(13)
       u1s=dxcdt+drdt*(x1_i(i-1)-xc)/radius
       u1(i-1,j,k)=alfa_s*u1s+beta_s*u1(i-1,j+1,k)+gama_s*u1(i-1,j+2,k)
      end if

      y=yc+sqrt(radius**2-( (x1(i)-xc)**2 + (x3_k(k)-zc)**2 ))
      beta=(y-x2(j))/dx2_j(j)
      call beta_function(14)
      u3s=drdt*(x3_k(k)-zc)/radius
      u3(i,j,k)=alfa_s*u3s+beta_s*u3(i,j+1,k)+gama_s*u3(i,j+2,k)

      y=yc+sqrt(radius**2-( (x1(i)-xc)**2 + (x3_k(k-1)-zc)**2 ))
      beta=(y-x2(j))/dx2_j(j)
      call beta_function(15)
      u3s=drdt*(x3_k(k-1)-zc)/radius
      u3(i,j,k-1)=alfa_s*u3s+beta_s*u3(i,j+1,k-1)+gama_s*u3(i,j+2,k-1)

      return
      end


      subroutine y_minus
      include 'var.cmn'
      real*8  y,u1s,u2s,u3s

      y=yc-sqrt(radius**2-( (x1(i)-xc)**2 + (x3(k)-zc)**2 ))
      beta=-(y-x2_j(j-1))/dx2(j)
      call beta_function(16)
      u2s=drdt*(y-yc)/radius
      u2(i,j-1,k)=alfa_s*u2s+beta_s*u2(i,j-2,k)+gama_s*u2(i,j-3,k)

      y=yc-sqrt(radius**2-( (x1_i(i)-xc)**2 + (x3(k)-zc)**2 ))
      beta=-(y-x2(j))/dx2_j(j)
      call beta_function(17)
      u1s=dxcdt+drdt*(x1_i(i)-xc)/radius
      u1(i,j,k)=alfa_s*u1s+beta_s*u1(i,j-1,k)+gama_s*u1(i,j-2,k)

c     Exclude i=1 for the following point (x1_i(i-1)=0 (on the wall) and  
c     the argument of sqrt is negative)

      if (i.ge.2) then
       y=yc-sqrt(radius**2-( (x1_i(i-1)-xc)**2 + (x3(k)-zc)**2 ))
       beta=-(y-x2(j))/dx2_j(j)
       call beta_function(18)
       u1s=dxcdt+drdt*(x1_i(i-1)-xc)/radius
       u1(i-1,j,k)=alfa_s*u1s+beta_s*u1(i-1,j-1,k)+gama_s*u1(i-1,j-2,k)
      end if

      y=yc-sqrt(radius**2-( (x1(i)-xc)**2 + (x3_k(k)-zc)**2 ))
      beta=-(y-x2(j))/dx2_j(j)
      call beta_function(19)
      u3s=drdt*(x3_k(k)-zc)/radius
      u3(i,j,k)=alfa_s*u3s+beta_s*u3(i,j-1,k)+gama_s*u3(i,j-2,k)

      y=yc-sqrt(radius**2-( (x1(i)-xc)**2 + (x3_k(k-1)-zc)**2 ))
      beta=-(y-x2(j))/dx2_j(j)
      call beta_function(20)
      u3s=drdt*(x3_k(k-1)-zc)/radius
      u3(i,j,k-1)=alfa_s*u3s+beta_s*u3(i,j-1,k-1)+gama_s*u3(i,j-2,k-1)

      return
      end



      subroutine x_plus
      include 'var.cmn'
      real*8  x,u1s,u2s,u3s

      x=xc+sqrt(radius**2-( (x2(j)-yc)**2 + (x3(k)-zc)**2 )) 
      beta=(x-x1_i(i))/dx1(i)
      call beta_function(21)
      u1s=dxcdt+drdt*(x-xc)/radius
      u1(i,j,k)=alfa_s*u1s+beta_s*u1(i+1,j,k)+gama_s*u1(i+2,j,k)

      x=xc+sqrt(radius**2-( (x2_j(j)-yc)**2 + (x3(k)-zc)**2 ))
      beta=(x-x1(i))/dx1_i(i)
      call beta_function(22)
      u2s=drdt*(x2_j(j)-yc)/radius
      u2(i,j,k)=alfa_s*u2s+beta_s*u2(i+1,j,k)+gama_s*u2(i+2,j,k)

      x=xc+sqrt(radius**2-( (x2_j(j-1)-yc)**2 + (x3(k)-zc)**2 ))
      beta=(x-x1(i))/dx1_i(i)
      call beta_function(23)
      u2s=drdt*(x2_j(j-1)-yc)/radius
      u2(i,j-1,k)=alfa_s*u2s+beta_s*u2(i+1,j-1,k)+gama_s*u2(i+2,j-1,k)

      x=xc+sqrt(radius**2-( (x2(j)-yc)**2 + (x3_k(k)-zc)**2 ))
      beta=(x-x1(i))/dx1_i(i)
      call beta_function(24)
      u3s=drdt*(x3_k(k)-zc)/radius
      u3(i,j,k)=alfa_s*u3s+beta_s*u3(i+1,j,k)+gama_s*u3(i+2,j,k)

      x=xc+sqrt(radius**2-( (x2(j)-yc)**2 + (x3_k(k-1)-zc)**2 ))
      beta=(x-x1(i))/dx1_i(i)
      call beta_function(25)
      u3s=drdt*(x3_k(k-1)-zc)/radius
      u3(i,j,k-1)=alfa_s*u3s+beta_s*u3(i+1,j,k-1)+gama_s*u3(i+2,j,k-1)

      return
      end


      subroutine x_minus
      include 'var.cmn'
      real*8  x,u1s,u2s,u3s


      x=xc-sqrt(radius**2-( (x2(j)-yc)**2 + (x3(k)-zc)**2 )) 
      u1s=dxcdt+drdt*(x-xc)/radius


      if (i.le.2) then
c      linear extrapolation using u1s and wall value (zero) 
       u1(i-1,j,k)=u1s*x1_i(i-1)/max(x,0.25d0*dx1(1))  
      else
       beta=-(x-x1_i(i-1))/dx1(i)
       call beta_function(26)
       u1(i-1,j,k)=alfa_s*u1s+beta_s*u1(i-2,j,k)+gama_s*u1(i-3,j,k)
      end if


      x=xc-sqrt(radius**2-( (x2_j(j)-yc)**2 + (x3(k)-zc)**2 )) 
      u2s=drdt*(x2_j(j)-yc)/radius
      if (i.le.2) then
c      linear extrapolation using u2s and wall value (zero) 
       u2(i,j,k)=u2s*x1(i)/max(x,0.25d0*dx1(1))  
      else
       beta=-(x-x1(i))/dx1_i(i)
       call beta_function(27)
       u2(i,j,k)=alfa_s*u2s+beta_s*u2(i-1,j,k)+gama_s*u2(i-2,j,k)
      end if


      x=xc-sqrt(radius**2-( (x2_j(j-1)-yc)**2 + (x3(k)-zc)**2 )) 
      u2s=drdt*(x2_j(j-1)-yc)/radius
      if (i.le.2) then
c      linear extrapolation using u2s and wall value (zero) 
       u2(i,j-1,k)=u2s*x1(i)/max(x,0.25d0*dx1(1))
      else  
       beta=-(x-x1(i))/dx1_i(i)
       call beta_function(28)
       u2(i,j-1,k)=alfa_s*u2s+beta_s*u2(i-1,j-1,k)+gama_s*u2(i-2,j-1,k)
      end if



      x=xc-sqrt(radius**2-( (x2(j)-yc)**2 + (x3_k(k)-zc)**2 )) 
      u3s=drdt*(x3_k(k)-zc)/radius
      if (i.le.2) then
c      linear extrapolation using u3s and wall value (zero) 
       u3(i,j,k)=u3s*x1(i)/max(x,0.25d0*dx1(1))  
      else
       beta=-(x-x1(i))/dx1_i(i)
       call beta_function(29)
       u3(i,j,k)=alfa_s*u3s+beta_s*u3(i-1,j,k)+gama_s*u3(i-2,j,k)
      end if



      x=xc-sqrt(radius**2-( (x2(j)-yc)**2 + (x3_k(k-1)-zc)**2 )) 
      u3s=drdt*(x3_k(k-1)-zc)/radius
      if (i.le.2) then
c      linear extrapolation using u3s and wall value (zero) 
       u3(i,j,k-1)=u3s*x1(i)/max(x,0.25d0*dx1(1))  
      else
      beta=-(x-x1(i))/dx1_i(i)
      call beta_function(30)
      u3(i,j,k-1)=alfa_s*u3s+beta_s*u3(i-1,j,k-1)+gama_s*u3(i-2,j,k-1)
      end if


      

      return
      end


      subroutine bou_phi
      include 'var.cmn' 
      integer oo
      real*8 az,bz,cz,fac,kka,kkc,kkae,kkce

      real*8 Error_phi(0:im+1,0:jm+1,0:km+1),
     $ phi_old(0:im+1,0:jm+1,0:km+1),ER_phi


c     phi1=phi(i=1)
c     phi0=phi(i=0)
c     electrode at x=0: -condL*(phi1-phi0)/dx=ka*exp(aa*eta)+kc*exp(-ac*eta) 
c     condL is conductivity at left boundary, ka>0, kc >0 
c     kka=-0.5*ka*dx/condL < 0 
c     kkc=-0.5*kc*dx/condL < 0
c     0.5*(phi1-phi0)=kka*exp(aa*eta)+kkc*exp(-ac*eta) 
c     (phi0+phi1)/2 = -eta + phiLe   (assume phi=potential - potential
c     electrode). 
c     phi1-phi0 = phi1-(-phi1-2*eta+2*phiLe)=2*phi1+2*eta-2*phiLe  
c     (phi1-phiLe + eta) = kka*exp(aa*eta)+kkc*exp(-ac*eta) 
c     -eta+kkc*exp(-ac*eta)+kka*exp(aa*eta) = phi1-phiLe  
c     solve eta_new using linearization (Newton method) near  previous eta
c     az*(eta_new - eta) +bz = cz with:
c           cz = phi1-phiLe
c           bz = -eta+kkc*exp(-ac*eta)+kka*exp(aa*eta)
c           az = -1 - kkc*ac*exp(-ac*eta) +kka*aa*exp(aa*eta) 
c     eta_new = eta + (cz-bz)/az
 
c     phi0_new = 2*(-eta_new+phiLe)-phi1
      
c     compute aa, ac, ka, kc and ke
!$omp parallel do private(cz,fac,kka,kkc,kkae,kkce,bz,az)
      do k=1,km
      do j=1,jm
       
       do i=0,im
        conduct(i,j,k)=condfac*0.5d0*(c(i,j,k,2)+c(i+1,j,k,2)) 
       end do

       cz=phi(1,j,k)-phiLe

c      If needed perform the next 8 statements (Newton iteration) multiple times
 
       fac = -0.5d0*dx1_i(0)/conduct(0,j,k)
       kka = fac*ka(j,k)
       kkc = fac*kc(j,k)
       kkae = kka*exp(aa*eta(j,k))
       kkce = kkc*exp(-ac*eta(j,k))      
       bz=-eta(j,k)+kkce+kkae 
       az=-1d0-ac*kkce+aa*kkae
       eta(j,k)=eta(j,k) + (cz-bz)/az

       phiL(j,k)=-eta(j,k)+phiLe
c      Impose phiLe at left side, like a Dirichlet boundary condition
c      Second order


       phi(0,j,k)=2d0*phiL(j,k)-phi(1,j,k)

c      Simple Dirichlet boundary at right side
c      Second order
       phi(im+1,j,k)=2d0*phiR-phi(im,j,k)  

      end do
      end do

c     zero Neumann conditions at inflow and outflow boundaries
       phi(0:im+1,0,1:km)=phi(0:im+1,1,1:km)
       phi(0:im+1,jm+1,1:km)=phi(0:im+1,jm,1:km)
       phi(0:im+1,0:jm+1,0)=phi(0:im+1,0:jm+1,km)
       phi(0:im+1,0:jm+1,km+1)=phi(0:im+1,0:jm+1,1)



      if (bubble) then


        iter3=1
        do k=0,km+1
        do j=0,jm+1
        do i=0,im+1
         phi_old(i,j,k)=phi(i,j,k)
        end do
        end do
        end do
        call bulb_phi

        do k=0,km+1
        do j=0,jm+1
        do i=0,im+1
         Error_phi(i,j,k)=abs(phi(i,j,k)-phi_old(i,j,k))              
        end do
        end do
        end do
        Max_Err_phi=maxval(Error_phi)

c        if (iter3.gt.5000) then
c         write(*,*)'Iteration problem in bulb_phi'
c         stop
c        end if


      end if


      conduct(0:im,0,1:km)=conduct(0:im,1,1:km)
      conduct(0:im,jm+1,1:km)=conduct(0:im,jm,1:km)
      conduct(0:im,0:jm+1,0)=conduct(0:im,0:jm+1,1)
      conduct(0:im,0:jm+1,km+1)=conduct(0:im,0:jm+1,km)

      eta(0,1:km)=eta(1,1:km)
      eta(jm+1,1:km)=eta(jm,1:km)
      eta(0:jm+1,0)=eta(0:jm+1,1)
      eta(0:jm+1,km+1)=eta(0:jm+1,km)



      return
      end

      subroutine bulb_phi
      include 'var.cmn' 
      integer o,mm
      real*8  e0,e1,e2,e3,e4,e5,e6,e7
      real*8  a1,a2,a3,phi_x,phi_y,phi_z
      real*8  nx,ny,nz,length_n
      real*8  phi_ref

      mm=1
      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)


          e0=(-f01(o,mm)*
     $         phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)) 
     $        +f02(o,mm)*
     $         phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1)   
     $        +f03(o,mm)*
     $         phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm))   
     $        -f04(o,mm)*
     $         phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1)
     $        +f05(o,mm)*
     $         phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm))
     $        -f06(o,mm)*
     $         phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1)
     $        -f07(o,mm)*
     $         phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm))
     $        +f08(o,mm)*
     $         phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1))
     $             /denom(o,mm)

          e1=(f11(o,mm)*
     $       (phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm))-
     $        phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)))+
     $        f12(o,mm)*
     $       (phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1)-
     $        phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1))+
     $        f13(o,mm)*
     $       (phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm))-
     $        phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)))+
     $        f14(o,mm)*
     $       (phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1)-
     $        phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1)))
     $       /denom(o,mm)

          e2=( f21(o,mm)*
     $       (phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm))-
     $        phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)))+
     $        f22(o,mm)*
     $       (phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1)-
     $        phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1))+
     $        f23(o,mm)*
     $       (phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm))-
     $        phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)))+
     $        f24(o,mm)*
     $       (phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1)-
     $        phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1)))
     $       /denom(o,mm)


          e3=( f31(o,mm)*
     $        (phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm))-
     $         phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1))+
     $        f32(o,mm)*
     $       (phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1)-
     $        phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)))+
     $        f33(o,mm)*
     $       (phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1)-
     $        phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)))+
     $        f34(o,mm)*
     $       (phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm))-
     $        phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1)))
     $        /denom(o,mm)


          e4=(x3(nz_prob(o,mm))
     $      *(phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1)
     $      +phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1)
     $      -phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1)
     $      -phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1))+
     $       x3(nz_prob(o,mm)+1)*
     $      (phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm))
     $      +phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm))
     $      -phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm))
     $      -phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm))))
     $      /denom(o,mm)

          e5=(x2(ny_prob(o,mm))*
     $       (phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm))
     $       +phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1)
     $       -phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1)
     $       -phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)))+
     $        x2(ny_prob(o,mm)+1)*
     $       (phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1)
     $       +phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm))
     $       -phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm))
     $       -phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1)))
     $       /denom(o,mm)

          e6=(x1(nx_prob(o,mm))*
     $       (phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm))
     $       +phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1)
     $       -phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1)
     $       -phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)))+
     $        x1(nx_prob(o,mm)+1)*
     $       (phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1)
     $       +phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm))
     $       -phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm))
     $       -phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1)))
     $       /denom(o,mm)

          e7=(phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm))+
     $        phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1)+
     $        phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1)+
     $        phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm))
     $       -phi(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1)
     $       -phi(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm))
     $       -phi(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm))
     $       -phi(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1))
     $       /denom(o,mm)

         phi_ref=e0+e1*x_prob(o,mm)+e2*y_prob(o,mm)+e3*z_prob(o,mm)+
     $           e4*
     $           x_prob(o,mm)*y_prob(o,mm)+e5*x_prob(o,mm)*z_prob(o,mm)+
     $           e6*y_prob(o,mm)*z_prob(o,mm)+
     $           e7*x_prob(o,mm)*y_prob(o,mm)*z_prob(o,mm) 


         if (direction(o) .eq. 1) then
          phi(i,j,k)=alfa_s_P(o)*phi_ref+beta_s_P(o)*phi(i+1,j,k)+
     $              gama_s_P(o)*phi(i+2,j,k)

         else if (direction(o) .eq. -1) then
          phi(i,j,k)=alfa_s_P(o)*phi_ref+beta_s_P(o)*phi(i-1,j,k)+
     $              gama_s_P(o)*phi(i-2,j,k)

         else if (direction(o) .eq. 2) then
          phi(i,j,k)=alfa_s_P(o)*phi_ref+beta_s_P(o)*phi(i,j+1,k)+
     $              gama_s_P(o)*phi(i,j+2,k)

         else if (direction(o) .eq. -2) then
          phi(i,j,k)=alfa_s_P(o)*phi_ref+beta_s_P(o)*phi(i,j-1,k)+
     $              gama_s_P(o)*phi(i,j-2,k)

         else if (direction(o) .eq. 3) then
          phi(i,j,k)=alfa_s_P(o)*phi_ref+beta_s_P(o)*phi(i,j,k+1)+
     $              gama_s_P(o)*phi(i,j,k+2)

         else if (direction(o) .eq. -3) then
          beta=-(zb(o)-x3(k))/dx3_k(k)
          call beta_function(32)
          phi(i,j,k)=alfa_s_P(o)*phi_ref+beta_s_P(o)*phi(i,j,k-1)+
     $              gama_s_P(o)*phi(i,j,k-2)

         end if          


      end do


      return
      end

      subroutine prep_phi
      include 'var.cmn'
      real*8 coef0_tmp

c     determine matrix coefficients normalized by minus the diagonal element
c     only valid for solution of div(c * grad(phi))=0
!$omp parallel do 
      do k=0,km
      do j=0,jm
      do i=0,im
       work1(i,j,k)=(c(i,j,k,2)+c(i+1,j,k,2))/dx1_i(i)
       work2(i,j,k)=(c(i,j,k,2)+c(i,j+1,k,2))/dx2_j(j)
       work3(i,j,k)=(c(i,j,k,2)+c(i,j,k+1,2))/dx3_k(k)
      end do
      end do
      end do

!$omp parallel do private(coef0_tmp)
      do k=1,km
      do j=1,jm
      do i=1,im
       coef0_tmp=-work1(i-1,j,k)-work1(i,j,k)
     $           -work2(i,j-1,k)-work2(i,j,k)
     $           -work3(i,j,k-1)-work3(i,j,k)

       coef1_phi(i,j,k)=work1(i-1,j,k)/coef0_tmp
       coef2_phi(i,j,k)=work1(i,j,k)  /coef0_tmp
       coef3_phi(i,j,k)=work2(i,j-1,k)/coef0_tmp
       coef4_phi(i,j,k)=work2(i,j,k)  /coef0_tmp
       coef5_phi(i,j,k)=work3(i,j,k-1)/coef0_tmp
       coef6_phi(i,j,k)=work3(i,j,k)  /coef0_tmp
      end do
      end do
      end do

c     compute kinetic prefactors
      do k=1,km
      do j=1,jm
       ka(j,k)=ka0*(0.5d0*(c(0,j,k,2)+c(1,j,k,2))/c_ref(2))*
     $          sqrt(abs(0.5d0*(c(0,j,k,1)+c(1,j,k,1)))/c_ref(1)) 
       kc(j,k)=kc0*0.5d0*(c(0,j,k,3)+c(1,j,k,3))/c_ref(3) 
      end do
      end do

      ka(0,0:km)=ka(1,0:km)
      ka(jm+1,0:km)=ka(jm,0:km)
      ka(0:jm+1,0)=ka(0:jm+1,1)
      ka(0:jm+1,km+1)=ka(0:jm+1,km)

      kc(0,0:km)=kc(1,0:km)
      kc(jm+1,0:km)=kc(jm,0:km)
      kc(0:jm+1,0)=kc(0:jm+1,1)
      kc(0:jm+1,km+1)=kc(0:jm+1,km)


      return
      end 

      subroutine spatial_phi 
      include 'var.cmn'

c     Jacobi iteration step
c     Relaxation factor 1, so phi(i) at the rhs cancels out 
!$omp parallel do
      do k=1,km
      do j=1,jm
      do i=1,im
      if (typ(i,j,k).eq.1) then
       phi(i,j,k)=-(coef1_phi(i,j,k)*phi(i-1,j,k)+
     $              coef2_phi(i,j,k)*phi(i+1,j,k)+
     $              coef3_phi(i,j,k)*phi(i,j-1,k)+
     $              coef4_phi(i,j,k)*phi(i,j+1,k)+
     $              coef5_phi(i,j,k)*phi(i,j,k-1)+
     $              coef6_phi(i,j,k)*phi(i,j,k+1))
      end if
      end do
      end do
      end do

      return
      end

      subroutine over
      include 'var.cmn' 

      do k=0,km+1
      do j=0,jm+1
       
       eta(j,k)=-0.5d0*(phi(0,j,k)+phi(1,j,k))-0.6d0

      end do
      end do


      return
      end



      subroutine bou_c
      include 'var.cmn'
      integer n,oo
      real*8 dphi
      real*8 Error_c(0:im+1,0:jm+1,0:km+1,1:3),Max_Err_c,
     $ c_old(0:im+1,0:jm+1,0:km+1,1:3),ER_c


       Error_c(0:im+1,0:jm+1,0:km+1,1)=0d0

c      current1=10000d0
c      F1=96485d0
      
c      dphi1=current1/F1*dx1_i(0)


c     No values of c and phi are needed in the corner (dummy) points 
c     (0,0),(im+1,0),(0,jm+1),(im+1,jm+1). 

      do k=1,km
      do j=1,jm   
   
        dphi=conduct(0,j,k)*(phi(1,j,k)-phi(0,j,k))/Far           
        c(0,j,k,1)=c(1,j,k,1)+0.5d0*dphi/dif(1)
        c(0,j,k,2)=c(1,j,k,2)+0.5d0*dphi/dif(2)
        c(0,j,k,3)=c(1,j,k,3)-dphi/dif(3)

     
      end do
      end do


      do n=1,nmax
!$omp parallel do 
       do k=1,km
       do j=1,jm 
        c(im+1,j,k,n)=2d0*c_ref(n)-c(im,j,k,n)
       end do
       end do

c      Dirichlet condition at inflow boundary
c      zero Neumann conditions at outflow boundary

       c(0:im+1,0,1:km,n)=2d0*c_ref(n)-c(0:im+1,1,1:km,n)        
       c(0:im+1,jm+1,1:km,n)=c(0:im+1,jm,1:km,n)
       c(0:im+1,0:jm+1,0,n)=c(0:im+1,0:jm+1,km,n)        
       c(0:im+1,0:jm+1,km+1,n)=c(0:im+1,0:jm+1,1,n)

      end do


      iter4=0
      if (bubble) then
       ER_c=1e-3
       Max_Err_c=5d0
       do while (Max_Err_c .gt. ER_c)
        do n=2,3
        do k=0,km+1
        do j=0,jm+1
        do i=0,im+1
         c_old(i,j,k,n)=c(i,j,k,n)
        end do
        end do
        end do
        end do
        call bulb_c
        iter4=iter4+1
        do n=2,3
        do k=0,km+1
        do j=0,jm+1
        do i=0,im+1
         Error_c(i,j,k,n)=abs(c(i,j,k,n)-c_old(i,j,k,n))             
        end do
        end do
        end do
        end do

        Max_Err_c=maxval(Error_c)
c        if (iter4.gt.5000) then
c         write(*,*)'Iteration problem in bulb_c'
c         stop
c        end if

       end do


      end if


      return
      end


      subroutine spatial_c 
      include 'var.cmn'
      integer n



      do n=1,nmax      
!$omp parallel do
       do k=0,km
       do j=0,jm 
       do i=0,im
        work1(i,j,k)=dif(n)*(c(i+1,j,k,n)-c(i,j,k,n))/dx1_i(i)
     $            -u1(i,j,k)*0.5d0*(c(i,j,k,n)+c(i+1,j,k,n)) 
        work2(i,j,k)=dif(n)*(c(i,j+1,k,n)-c(i,j,k,n))/dx2_j(j) 
     $            -u2(i,j,k)*0.5d0*(c(i,j,k,n)+c(i,j+1,k,n)) 
        work3(i,j,k)=dif(n)*(c(i,j,k+1,n)-c(i,j,k,n))/dx3_k(k) 
     $            -u3(i,j,k)*0.5d0*(c(i,j,k,n)+c(i,j,k+1,n)) 
       end do
       end do
       end do

!$omp parallel do
       do k=1,km
       do j=1,jm
       do i=1,im
       if (typ(i,j,k).eq.1) then
        c(i,j,k,n)=c(i,j,k,n)+dtime*(
     $              (work1(i,j,k)-work1(i-1,j,k))/dx1(i)+
     $              (work2(i,j,k)-work2(i,j-1,k))/dx2(j)+ 
     $              (work3(i,j,k)-work3(i,j,k-1))/dx3(k)) 
       end if
       end do
       end do
       end do

      end do
 
      return
      end





      subroutine bulb_c
      include 'var.cmn' 
      integer o,n,mm
      real*8  a1,a2,a3,cx,cy,cz
      real*8  e0,e1,e2,e3,e4,e5,e6,e7
      real*8  nx,ny,nz,length_n
      mm=1
      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)


         do n=2,3

          e0=(-f01(o,mm)*
     $         c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm),n) 
     $        +f02(o,mm)*
     $         c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1,n)   
     $        +f03(o,mm)*
     $         c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm),n)   
     $        -f04(o,mm)*
     $         c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)
     $        +f05(o,mm)*
     $         c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm),n)
     $        -f06(o,mm)*
     $         c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1,n)
     $        -f07(o,mm)*
     $         c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm),n)
     $        +f08(o,mm)*
     $         c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1,n))
     $        /denom(o,mm)

          e1=(f11(o,mm)*
     $       (c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm),n)-
     $        c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm),n))+
     $        f12(o,mm)*
     $       (c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1,n)-
     $        c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1,n))+
     $        f13(o,mm)*
     $       (c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm),n)-
     $        c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm),n))+
     $        f14(o,mm)*
     $       (c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)-
     $        c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)))
     $       /denom(o,mm)

          e2=(f21(o,mm)*
     $       (c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm),n)-
     $        c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm),n))+
     $        f22(o,mm)*
     $       (c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)-
     $        c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1,n))+
     $        f23(o,mm)*
     $       (c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm),n)-
     $        c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm),n))+
     $        f24(o,mm)*
     $       (c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1,n)-
     $        c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)))
     $       /denom(o,mm)

          e3=(f31(o,mm)*
     $       (c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm),n)-
     $        c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1,n))+
     $        f32(o,mm)*
     $       (c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)-
     $        c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm),n))+
     $        f33(o,mm)*
     $       (c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1,n)-
     $        c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm),n))+
     $        f34(o,mm)*
     $       (c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm),n)-
     $        c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)))
     $       /denom(o,mm)


          e4=(x3(nz_prob(o,mm))*
     $       (c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1,n)
     $       +c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)
     $       -c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)
     $       -c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1,n))+
     $        x3(nz_prob(o,mm)+1)*
     $       (c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm),n)
     $       +c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm),n)
     $       -c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm),n)
     $       -c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm),n)))
     $       /denom(o,mm)

          e5=(x2(ny_prob(o,mm))
     $      *(c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm),n)
     $       +c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)
     $       -c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)
     $       -c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm),n))+
     $        x2(ny_prob(o,mm)+1)*
     $       (c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1,n)
     $       +c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm),n)
     $       -c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm),n)
     $       -c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1,n)))
     $       /denom(o,mm)

          e6=(x1(nx_prob(o,mm))*
     $       (c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm),n)
     $       +c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)
     $       -c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1,n)
     $       -c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm),n))+
     $        x1(nx_prob(o,mm)+1)*
     $       (c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1,n)
     $       +c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm),n)
     $       -c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm),n)
     $       -c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)))
     $       /denom(o,mm)

          e7=(c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm),n)+
     $        c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm)+1,n)+
     $        c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm)+1,n)+
     $        c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm),n)
     $       -c(nx_prob(o,mm),ny_prob(o,mm),nz_prob(o,mm)+1,n)
     $       -c(nx_prob(o,mm),ny_prob(o,mm)+1,nz_prob(o,mm),n)
     $       -c(nx_prob(o,mm)+1,ny_prob(o,mm),nz_prob(o,mm),n)
     $      -c(nx_prob(o,mm)+1,ny_prob(o,mm)+1,nz_prob(o,mm)+1,n))
     $       /denom(o,mm)

         c_s(n)=e0+e1*x_prob(o,mm)+e2*y_prob(o,mm)+e3*z_prob(o,mm)+
     $          e4*
     $          x_prob(o,mm)*y_prob(o,mm)+e5*x_prob(o,mm)*z_prob(o,mm)+
     $          e6*y_prob(o,mm)*z_prob(o,mm)+
     $          e7*x_prob(o,mm)*y_prob(o,mm)*z_prob(o,mm)

         end do 

         c_s(1)=c_ref(1)



         do n=1,3
         if (direction(o) .eq. 1) then
          c(i,j,k,n)=alfa_s_P(o)*c_s(n)+beta_s_P(o)*c(i+1,j,k,n)+
     $              gama_s_P(o)*c(i+2,j,k,n)

         else if (direction(o) .eq. -1) then
          c(i,j,k,n)=alfa_s_P(o)*c_s(n)+beta_s_P(o)*c(i-1,j,k,n)
     $              +gama_s_P(o)*c(i-2,j,k,n)

         else if (direction(o) .eq. 2) then
         c(i,j,k,n)=alfa_s_P(o)*c_s(n)+beta_s_P(o)*c(i,j+1,k,n)
     $             +gama_s_P(o)*c(i,j+2,k,n)

         else if (direction(o) .eq. -2) then
          c(i,j,k,n)=alfa_s_P(o)*c_s(n)+beta_s_P(o)*c(i,j-1,k,n)
     $              +gama_s_P(o)*c(i,j-2,k,n)


         else if (direction(o) .eq. 3) then
          c(i,j,k,n)=alfa_s_P(o)*c_s(n)+beta_s_P(o)*c(i,j,k+1,n)
     $              +gama_s_P(o)*c(i,j,k+2,n)


         else if (direction(o) .eq. -3) then
          c(i,j,k,n)=alfa_s_P(o)*c_s(n)+beta_s_P(o)*c(i,j,k-1,n)
     $              +gama_s_P(o)*c(i,j,k-2,n)

         end if
         end do
          


      end do


      return
      end




      subroutine bulb_p
      include 'var.cmn' 
      integer o


      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)

         if (direction(o) .eq. 1) then
            p(i,j,k)=p(i+1,j,k)

         else if (direction(o) .eq. -1) then
            p(i,j,k)=p(i-1,j,k)

         else if (direction(o) .eq. 2) then
            p(i,j,k)=p(i,j+1,k)

         else if (direction(o) .eq. -2) then
            p(i,j,k)=p(i,j-1,k)

         else if (direction(o) .eq. 3) then
            p(i,j,k)=p(i,j,k+1)

         else if (direction(o) .eq. -3) then
            p(i,j,k)=p(i,j,k-1)

         end if

      end do


      return
      end

c     writing to 51 is now in bulb_flux_write, because it is not needed 
c     at each time step

      subroutine bulb_flux
      include 'var.cmn' 
      integer o
      real*8 term1,term2

      total=0d0

      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)

         if (x1_minus(o) .eq. 1) then

            term1=-0.5d0*(c(i,j,k,1)+c(i-1,j,k,1))*
     $               u1(i-1,j,k)*dx2(j)*dx3(k)

            term2=dif(1)*(c(i,j,k,1)-c(i-1,j,k,1))
     $              /dx1_i(i-1)*dx2(j)*dx3(k)

            if (i .eq. nx_start) then
             term1=0d0
             term2=0d0
            end if

            total=total+term2

c 2          format(A4,4i6,1000e14.6)
c            if (term2.gt.0d0) then
c             write(20,2)'x1m',o,i,j,k,term1,term2,
c     $                  c(i,j,k,1),c(i-1,j,k,1)
c            end if

         end if


         if (x1_plus(o) .eq. 1) then

            term1=-(-0.5d0*(c(i+1,j,k,1)+c(i,j,k,1))*
     $               u1(i,j,k)*dx2(j)*dx3(k))

            term2=-(dif(1)*(c(i+1,j,k,1)-c(i,j,k,1))
     $              /dx1_i(i)*dx2(j)*dx3(k))

            total=total+term2

c            if (term2.gt.0d0) then
c             write(20,2)'x1m',o,i,j,k,term1,term2,
c     $                  c(i,j,k,1),c(i+1,j,k,1)
c            end if

         end if

         if (x2_minus(o) .eq. 1) then


            term1=-0.5d0*(c(i,j,k,1)+c(i,j-1,k,1))*
     $               u2(i,j-1,k)*dx1(i)*dx3(k)

            term2=dif(1)*(c(i,j,k,1)-c(i,j-1,k,1))
     $              /dx2_j(j-1)*dx1(i)*dx3(k)

            total=total+term2


c            if (term2.gt.0d0) then
c             write(20,2)'x2m',o,i,j,k,term1,term2,
c     $                  c(i,j,k,1),c(i,j-1,k,1)
c            end if

         end if


         if (x2_plus(o) .eq. 1) then

            term1=-(-0.5d0*(c(i,j+1,k,1)+c(i,j,k,1))*
     $               u2(i,j,k)*dx1(i)*dx3(k))

            term2=-(dif(1)*(c(i,j+1,k,1)-c(i,j,k,1))
     $              /dx2_j(j)*dx1(i)*dx3(k))

            total=total+term2

c            if (term2.gt.0d0) then
c             write(20,2)'x2p',o,i,j,k,term1,term2,
c     $                  c(i,j,k,1),c(i,j+1,k,1)
c            end if

         end if

         if (x3_minus(o) .eq. 1) then

            term1=-0.5d0*(c(i,j,k,1)+c(i,j,k-1,1))*
     $               u3(i,j,k-1)*dx1(i)*dx2(j)

            term2=dif(1)*(c(i,j,k,1)-c(i,j,k-1,1))
     $              /dx3_k(k-1)*dx1(i)*dx2(j)

            total=total+term2

c            if (term2.gt.0d0) then
c             write(20,2)'x3m',o,i,j,k,term1,term2,
c     $                  c(i,j,k,1),c(i,j,k-1,1)
c            end if

         end if

         if (x3_plus(o) .eq. 1) then

            term1=-(-0.5d0*(c(i,j,k+1,1)+c(i,j,k,1))*
     $               u3(i,j,k)*dx1(i)*dx2(j))

            term2=-(dif(1)*(c(i,j,k+1,1)-c(i,j,k,1))
     $              /dx3_k(k)*dx1(i)*dx2(j))

            total=total+term2

c            if (term2.gt.0d0) then
c             write(20,2)'x3p',o,i,j,k,term1,term2,
c     $                  c(i,j,k,1),c(i,j,k+1,1)
c            end if

         end if

      end do

      return
      end


      subroutine bulb_flux_write
      include 'var.cmn' 
      integer o,counter_write
      real*8 term1,term2

 100   format(4I5,1000e14.6)
      total=0d0
      counter_write=0

      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)


         if (x1_minus(o) .eq. 1) then


            term1=-0.5d0*(c(i,j,k,1)+c(i-1,j,k,1))*
     $               u1(i-1,j,k)*dx2(j)*dx3(k)

            term2=dif(1)*(c(i,j,k,1)-c(i-1,j,k,1))
     $              /dx1_i(i-1)*dx2(j)*dx3(k)

            if (i .eq. nx_start) then

             term1=0d0
             term2=0d0

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,-1,term1,term2

             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,-1,term1,term2

             end if

             close(51)

            end if

             
            total=total+term2
 2          format(A4,4i6,1000e14.6)
            if (term2.gt.0d0) then
             write(20,2)'x1m',o,i,j,k,term1,term2,
     $                  c(i,j,k,1),c(i-1,j,k,1)
            end if
         end if


         if (x1_plus(o) .eq. 1) then


            term1=-(-0.5d0*(c(i+1,j,k,1)+c(i,j,k,1))*
     $               u1(i,j,k)*dx2(j)*dx3(k))

            term2=-(dif(1)*(c(i+1,j,k,1)-c(i,j,k,1))
     $              /dx1_i(i)*dx2(j)*dx3(k))

            if (i .eq. nx_start) then

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,1,term1,term2

             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,1,term1,term2

             end if

             close(51)

            end if

              
            total=total+term2
            if (term2.gt.0d0) then
             write(20,2)'x1m',o,i,j,k,term1,term2,
     $                  c(i,j,k,1),c(i+1,j,k,1)
            end if

         end if

         if (x2_minus(o) .eq. 1) then


            term1=-0.5d0*(c(i,j,k,1)+c(i,j-1,k,1))*
     $               u2(i,j-1,k)*dx1(i)*dx3(k)

            term2=dif(1)*(c(i,j,k,1)-c(i,j-1,k,1))
     $              /dx2_j(j-1)*dx1(i)*dx3(k)

            if (i .eq. nx_start) then

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,-2,term1,term2
             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,-2,term1,term2
             end if

             close(51)

            end if

              
            total=total+term2
            if (term2.gt.0d0) then
             write(20,2)'x2m',o,i,j,k,term1,term2,
     $                  c(i,j,k,1),c(i,j-1,k,1)
            end if

         end if


         if (x2_plus(o) .eq. 1) then


            term1=-(-0.5d0*(c(i,j+1,k,1)+c(i,j,k,1))*
     $               u2(i,j,k)*dx1(i)*dx3(k))

            term2=-(dif(1)*(c(i,j+1,k,1)-c(i,j,k,1))
     $              /dx2_j(j)*dx1(i)*dx3(k))

            if (i .eq. nx_start) then

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,2,term1,term2
             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,2,term1,term2
             end if

             close(51)

            end if

              
            total=total+term2
            if (term2.gt.0d0) then
             write(20,2)'x2p',o,i,j,k,term1,term2,
     $                  c(i,j,k,1),c(i,j+1,k,1)
            end if

         end if

         if (x3_minus(o) .eq. 1) then


            term1=-0.5d0*(c(i,j,k,1)+c(i,j,k-1,1))*
     $               u3(i,j,k-1)*dx1(i)*dx2(j)

            term2=dif(1)*(c(i,j,k,1)-c(i,j,k-1,1))
     $              /dx3_k(k-1)*dx1(i)*dx2(j)

            if (i .eq. nx_start) then

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,-3,term1,term2
             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,-3,term1,term2
             end if

             close(51)

            end if

              
            total=total+term2
            if (term2.gt.0d0) then
             write(20,2)'x3m',o,i,j,k,term1,term2,
     $                  c(i,j,k,1),c(i,j,k-1,1)
            end if

         end if



         if (x3_plus(o) .eq. 1) then


            term1=-(-0.5d0*(c(i,j,k+1,1)+c(i,j,k,1))*
     $               u3(i,j,k)*dx1(i)*dx2(j))

            term2=-(dif(1)*(c(i,j,k+1,1)-c(i,j,k,1))
     $              /dx3_k(k)*dx1(i)*dx2(j))

            if (i .eq. nx_start) then

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,3,term1,term2
             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,3,term1,term2
             end if

             close(51)

            end if

              
            total=total+term2
            if (term2.gt.0d0) then
             write(20,2)'x3p',o,i,j,k,term1,term2,
     $                  c(i,j,k,1),c(i,j,k+1,1)
            end if

         end if



      end do


      return
      end



      subroutine domain_flux
      include 'var.cmn' 
      integer o
      real*8 term1,term2,TERM3,TERM4,TERM5,TERM6,
     $       A12(1:jm,1:km),B12(1:jm,1:km),
     $       A21(1:im,1:km),B21(1:im,1:km),
     $       A22(1:im,1:km),B22(1:im,1:km),
     $       A32(1:im,1:jm),B32(1:im,1:jm),
     $       Error


      do k=1,km
      do j=1,jm

         A12(j,k)=-(dif(1)*(c(1,j,k,1)-c(0,j,k,1))/dx1_i(0)
     $          *dx2(j)*dx3(k))

         B12(j,k)=dif(1)*(c(im+1,j,k,1)-c(im,j,k,1))/dx1_i(im)
     $          *dx2(j)*dx3(k)

      end do
      end do

      do o=1,change-1
         i=x_IB(o)
         if (i .eq. nx_start) then
          A12(y_IB(o),z_IB(o))=0d0
         end if
      end do


      term1=sum(A12)
      term2=sum(B12)


      do k=1,km
      do i=1,im

         A21(i,k)=-(-0.5d0*(c(i,1,k,1)+c(i,0,k,1))
     $              *u2(i,0,k)*dx1(i)*dx3(k))

         B21(i,k)=-0.5d0*(c(i,jm+1,k,1)+c(i,jm,k,1))
     $              *u2(i,jm,k)*dx1(i)*dx3(k)


         A22(i,k)=-(dif(1)*(c(i,1,k,1)-c(i,0,k,1))/dx2_j(0)
     $           *dx1(i)*dx3(k))

         B22(i,k)=dif(1)*(c(i,jm+1,k,1)-c(i,jm,k,1))/dx2_j(jm)
     $           *dx1(i)*dx3(k)


      end do
      end do


      term3=sum(A21)+sum(A22)
      term4=sum(b21)+sum(B22)

      do j=1,jm
      do i=1,im


         A32(i,j)=-(dif(1)*(c(i,j,1,1)-c(i,j,0,1))/dx3_k(0)
     $           *dx1(i)*dx2(j))
         B32(i,j)=dif(1)*(c(i,j,km+1,1)-c(i,j,km,1))/dx3_k(km)
     $           *dx1(i)*dx2(j)


      end do
      end do

      term5=sum(A32)
      term6=sum(B32)


      Error=abs(term1+term2+term3+term4+term5+term6+total)
      write(*,*)'Error of hydrogen balance',Error
      write(6,60)term1,term2,term3,term4,term5,term6,total
 60   format(1000e14.6)


      return
      end



      subroutine growing
      include 'var.cmn' 
      integer o,oo
      real*8 Rid,temp,sigma,delta_v,P_IB_fluid,P_outside,
     $       P_bubble,R_new,dr,pi,height0,az,bz




c     average pressure over the IB_fluid cells
      P_IB_fluid=0d0

      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)

         P_IB_fluid=P_IB_fluid+p(i,j,k)

      end do
 
      P_outside=P_IB_fluid/(change-1)

      Rid=8.3142d0
      temp=353d0
      sigma=7.7d-2


c     laplace pressure equation
      P_bubble=P_outside+(2d0/radius)*sigma+1d5


      delta_v=-1d0*total*Rid*temp*dtime/P_bubble

      pi=2d0*asin(1d0)
      radius0=radius

      height0=0d0
      height=0d0

      radius=(3d0/4d0/pi*delta_v+radius**3)**(1d0/3d0)
      drdt=(radius-radius0)/dtime

c      write(*,*)chord,radius,height
      return
      end





      subroutine bulb_flux_force
      include 'var.cmn' 
      integer o,counter_write
      real*8 term1,term2,term3


 100   format(4I5,1000e14.6)
      total=0d0
      counter_write=0
      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)


         if (x1_minus(o) .eq. 1) then

            term1=-(0.5d0*(u1(i-1,j,k)+u1(i,j,k)))**2*dx2(j)*dx3(k)

            term2=-p(i,j,k)/ro*dx2(j)*dx3(k)

            term3=nu*(u1(i,j,k)-u1(i-1,j,k))/dx1(i)*dx2(j)*dx3(k)


            if (i .eq. nx_start) then

             term1=0d0
             term2=0d0
             term3=0d0

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,-1,term1,term2,term3

             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,-1,term1,term2,term3

             end if

             close(51)

            end if


              
            total=total+term1+term2+term3

         end if


         if (x1_plus(o) .eq. 1) then


            term1=-(-(0.5d0*(u1(i+1,j,k)+u1(i,j,k)))**2*dx2(j)*dx3(k))

            term2=-(-p(i+1,j,k)/ro*dx2(j)*dx3(k))
    
            term3=-(nu*(u1(i+1,j,k)-u1(i,j,k))/dx1(i+1)*dx2(j)*dx3(k)) 

            if (i .eq. nx_start) then

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,1,term1,term2,term3

             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,1,term1,term2,term3

             end if

             close(51)

            end if

              
            total=total+term1+term2+term3

         end if

         if (x2_minus(o) .eq. 1) then


            term1=-0.5d0*(u1(i,j,k)+u1(i,j-1,k))
     $     *0.5d0*(u2(i,j-1,k)+u2(i+1,j-1,k))*dx1(i)*dx3(k)

            term2=0d0
    
            term3=nu*(u1(i,j,k)-u1(i,j-1,k))/dx2_j(j-1)*dx1(i)*dx3(k)    

            if (i .eq. nx_start) then

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,-2,term1,term2,term3

             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,-2,term1,term2,term3

             end if
             close(51)

            end if
            
            total=total+term1+term2+term3

         end if


         if (x2_plus(o) .eq. 1) then


            term1=-(-0.5d0*(u1(i,j,k)+u1(i,j+1,k))
     $        *0.5d0*(u2(i,j+1,k)+u2(i+1,j+1,k))*dx1(i)*dx3(k))

            term2=0d0
    
            term3=-(nu*(u1(i,j+1,k)-u1(i,j,k))/dx2_j(j)*dx1(i)*dx3(k)) 

            if (i .eq. nx_start) then

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,2,term1,term2,term3

             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,2,term1,term2,term3

             end if

             close(51)

            end if
       
            total=total+term1+term2+term3

         end if

         if (x3_minus(o) .eq. 1) then


            term1=-0.5d0*(u1(i,j,k)+u1(i,j,k-1))
     $     *0.5d0*(u3(i,j,k-1)+u3(i+1,j,k-1))*dx1(i)*dx2(j)

            term2=0d0
    
            term3=nu*(u1(i,j,k)-u1(i,j,k-1))/dx3_k(k-1)*dx1(i)*dx2(j) 

            if (i .eq. nx_start) then

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,-3,term1,term2,term3

             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,-3,term1,term2,term3

             end if

             close(51)

            end if
              
            total=total+term1+term2+term3

         end if



         if (x3_plus(o) .eq. 1) then


            term1=-(-0.5d0*(u1(i,j,k)+u1(i,j,k+1))
     $     *0.5d0*(u3(i,j,k+1)+u3(i+1,j,k+1))*dx1(i)*dx2(j))

            term2=0d0
    
            term3=-(nu*(u1(i,j,k+1)-u1(i,j,k))/dx3_k(k)*dx1(i)*dx2(j))   

            if (i .eq. nx_start) then

             counter_write=counter_write+1

             if (counter_write .eq. 1) then

              open(unit=51,file='node.dat')
              write(51,100)i,j,k,3,term1,term2,term3

             else 

              open(unit=51,file='node.dat',position='append')
              write(51,100)i,j,k,3,term1,term2,term3

             end if

             close(51)

            end if

            total=total+term1+term2+term3

         end if



      end do


      return
      end



      subroutine domain_flux_force
      include 'var.cmn' 
      integer o       
      real*8 term1,term2,term3,term4,term5,term6,
     $       u_i(1:im,1:km),u_o(1:im,1:km),
     $       w_i(1:im,1:km),w_o(1:im,1:km),
     $       u_f(1:im,1:jm),u_b(1:im,1:jm),
     $       v_f(1:im,1:jm),v_b(1:im,1:jm),
     $       A11(1:jm,1:km),B11(1:jm,1:km),
     $       A12(1:jm,1:km),B12(1:jm,1:km),
     $       A13(1:jm,1:km),B13(1:jm,1:km),
     $       A21(1:im,1:km),B21(1:im,1:km),
     $       A22(1:im,1:km),B22(1:im,1:km),
     $       A23(1:im,1:km),B23(1:im,1:km),
     $       A31(1:im,1:jm),B31(1:im,1:jm),
     $       A32(1:im,1:jm),B32(1:im,1:jm),
     $       A33(1:im,1:jm),B33(1:im,1:jm),Error


      do k=1,km
      do j=1,jm

         A11(j,k)=0d0
         B11(j,k)=0d0

         A12(j,k)=-(-0.5d0*(p(0,j,k)+p(1,j,k))*dx2(j)*dx3(k)/ro)
         B12(j,k)=-0.5d0*(p(im,j,k)+p(im+1,j,k))*dx2(j)*dx3(k)/ro

         A13(j,k)=0d0
         B13(j,k)=0d0


      end do
      end do

      do o=1,change-1
         i=x_IB(o)
         if (i .eq. nx_start) then
          A12(y_IB(o),z_IB(o))=0d0
         end if
      end do


      term1=sum(A11)+sum(A12)+sum(A13)
      term2=sum(B11)+sum(B12)+sum(B13)


      do k=1,km
      do i=1,im

         u_i(i,k)=0.25d0*(u1(i-1,0,k)+u1(i-1,1,k)
     $           +u1(i,0,k)+u1(i,1,k))
         u_o(i,k)=0.25d0*(u1(i-1,jm,k)+u1(i-1,jm+1,k)
     $           +u1(i,jm,k)+u1(i,jm+1,k))
         w_i(i,k)=0.25d0*(u3(i,0,k-1)+u3(i,1,k-1)+u3(i,0,k)+u3(i,1,k))
         w_o(i,k)=0.25d0*(u3(i,jm,k-1)+u3(i,jm+1,k-1)
     $         +u3(i,jm,k)+u3(i,jm+1,k))

      end do
      end do


      do k=1,km
      do i=1,im

         A21(i,k)=-(-u_i(i,k)*u2(i,0,k)*dx1(i)*dx3(k))
         B21(i,k)=-u_o(i,k)*u2(i,jm,k)*dx1(i)*dx3(k)
   
         A22(i,k)=0d0
         B22(i,k)=0d0

         A23(i,k)=-(0.5d0*nu*((u1(i,1,k)+u1(i-1,1,k))
     $            -(u1(i,0,k)+u1(i-1,0,k)))*dx1(i)*dx3(k)/dx2_j(0))
         B23(i,k)=0.5d0*nu*((u1(i,jm+1,k)+u1(i-1,jm+1,k))
     $            -(u1(i,jm,k)+u1(i-1,jm,k)))*dx1(i)*dx3(k)/dx2_j(jm)


      end do
      end do

      term3=sum(A21)+sum(A22)+sum(A23)
      term4=sum(B21)+sum(B22)+sum(B23)


      do j=1,jm
      do i=1,im

         u_f(i,j)=0.25d0*(u1(i-1,j,0)+u1(i-1,j,1)
     $            +u1(i,j,0)+u1(i,j,1))
         u_b(i,j)=0.25d0*(u1(i-1,j,km)+u1(i-1,j,km+1)
     $            +u1(i,j,km)+u1(i,j,km+1))
         v_f(i,j)=0.25d0*(u2(i,j-1,0)+u2(i,j-1,1)
     $            +u2(i,j,0)+u2(i,j,1))
         v_b(i,j)=0.25d0*(u2(i,j-1,km)+u2(i,j-1,km+1)
     $            +u2(i,j,km)+u2(i,j,km+1))

      end do
      end do


      do j=1,jm
      do i=1,im

         A31(i,j)=-(-u_f(i,j)*u3(i,j,0)*dx1(i)*dx2(j))
         B31(i,j)=-u_b(i,j)*u3(i,j,km)*dx1(i)*dx2(j)

         A32(i,j)=0d0
         B32(i,j)=0d0

         A33(i,j)=-(0.5d0*nu*((u1(i,j,1)+u1(i-1,j,1))-(u1(i,j,0)
     $            +u1(i-1,j,0)))*dx1(i)*dx2(j)/dx3_k(0))
         B33(i,j)=0.5d0*nu*((u1(i,j,km+1)+u1(i-1,j,km+1))-(u1(i,j,km)
     $            +u1(i-1,j,km)))*dx1(i)*dx2(j)/dx3_k(km)

      end do
      end do

      term5=sum(A31)+sum(A32)+sum(A33)
      term6=sum(B31)+sum(B32)+sum(B33)


      Error=abs(term1+term2+term3+term4+term5+term6+total)
      write(*,*)'Error of force balance',Error
      write(*,*)term1,term2,term3,term4,term5,term6,total



      return
      end



