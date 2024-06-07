      subroutine u 
      include 'var.cmn'
      real*8  tp,ua,va,wa
  


       if (rung.eq.1) then
          storage_u=u1
       end if

c     Following code is for the descretisation of the u-momentum equation
c     at faces 1 to Im-1 of the interior cells, based upon the indexing 
c     of the cell nodes.    
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

        rhu(i,j,k)=storage_u(i,j,k)-tp*dtime*alpha_rung(rung)
     $  	-nu*dtime*alpha_rung(rung)/dx1_i(i)*
     $       ((u1(i,j,k)-u1(i-1,j,k))/dx1(i)
     $          -(u1(i+1,j,k)-u1(i,j,k))/dx1(i+1))
     $	        -nu*dtime*alpha_rung(rung)/dx2(j)*
     $        ((u1(i,j,k)-u1(i,j-1,k))/dx2_j(j-1)
     $          -(u1(i,j+1,k)-u1(i,j,k))/dx2_j(j))
     $	        -nu*dtime*alpha_rung(rung)/dx3(k)*
     $        ((u1(i,j,k)-u1(i,j,k-1))/dx3_k(k-1)
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
      do i=1,im
      if (typ1(i,j,k).eq.1) then

       u1(i,j,k)=rhu(i,j,k)

      end if
      end do
      end do
      end do

      do k=1,km
      do j=1,jm

        u1(0,j,k)=0d0
        u1(im,j,k)=0d0 

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


      u1(0:im+1,0,1:km)=-u1(0:im+1,1,1:km)

      u1(0:im+1,jm+1,1:km)=u1(0:im+1,jm,1:km)


c     Periodic boundary conditions

      u1(0:im+1,0:jm+1,0)=u1(0:im+1,0:jm+1,km)
      u1(0:im+1,0:jm+1,km+1)=u1(0:im+1,0:jm+1,1)



      return
      end




      subroutine v 
      include 'var.cmn'
      real*8  tp,ua,va,wa

       if (rung.eq.1) then
          storage_v=u2
       end if
      


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

        rhv(i,j,k)=storage_v(i,j,k)-tp*dtime*alpha_rung(rung)
     $  	-nu*dtime*alpha_rung(rung)/dx1(i)*
     $       ((u2(i,j,k)-u2(i-1,j,k))/dx1_i(i-1)
     $          -(u2(i+1,j,k)-u2(i,j,k))/dx1_i(i))
     $	        -nu*dtime*alpha_rung(rung)/dx2_j(j)*
     $         ((u2(i,j,k)-u2(i,j-1,k))/dx2(j)
     $          -(u2(i,j+1,k)-u2(i,j,k))/dx2(j+1))
     $	        -nu*dtime*alpha_rung(rung)/dx3(k)*
     $         ((u2(i,j,k)-u2(i,j,k-1))/dx3_k(k-1)
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


      do k=1,km
      do i=1,im

       u2(i,0,k)=4d0*umax*(x1(i)/lx1)*(1d0-(x1(i)/lx1))

      end do
      end do

      u2(1:im,jm+1,1:km)=u2(1:im,jm,1:km)


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

      if (rung.eq.1) then

        storage_w=u3
      end if
      
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

        rhw(i,j,k)=storage_w(i,j,k)-tp*dtime*alpha_rung(rung)
     $  	-nu*dtime*alpha_rung(rung)/dx1(i)*
     $      ((u3(i,j,k)-u3(i-1,j,k))/dx1_i(i-1)
     $          -(u3(i+1,j,k)-u3(i,j,k))/dx1_i(i))
     $	        -nu*dtime*alpha_rung(rung)/dx2(j)*
     $         ((u3(i,j,k)-u3(i,j-1,k))/dx2_j(j-1)
     $          -(u3(i,j+1,k)-u3(i,j,k))/dx2_j(j))
     $	        -nu*dtime*alpha_rung(rung)/dx3(k)*
     $        ((u3(i,j,k)-u3(i,j,k-1))/dx3(k)
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



      return
      end


      subroutine bou_w
      include 'var.cmn'

       
c     Simple Dirichlet boundary at left and right side
c     Second order

      u3(0,1:jm,1:km)=-u3(1,1:jm,1:km)
      u3(im+1,1:jm,1:km)=-u3(im,1:jm,1:km) 


      u3(0:im+1,0,1:km)=-u3(0:im+1,1,1:km)
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




!$omp parallel do private(i,j,k)
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

!$omp parallel do private(i,j,k)
      do k=1,km
      do j=1,jm
      do i=1,im

       if(typ1(i,j,k).eq.1) then
        u1(i,j,k)=u1(i,j,k)
     $  -dtime*alpha_rung(rung)/ro/dx1_i(i)*(
     $          p(i+1,j,k)-p(i,j,k))
       end if



       if(typ2(i,j,k).eq.1) then
        u2(i,j,k)=u2(i,j,k)
     $  -dtime*alpha_rung(rung)/ro/dx2_j(j)
     $          *(p(i,j+1,k)-p(i,j,k))
       end if



       if(typ3(i,j,k).eq.1) then
        u3(i,j,k)=u3(i,j,k)
     $   -dtime*alpha_rung(rung)/ro/dx3_k(k)*(
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

      return
      end


      subroutine y_plus
      include 'var.cmn'
      real*8  y,u1s,u2s,u3s


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

 
       phi(im+1,j,k)=2d0*phiR-phi(im,j,k) 

      end do
      end do

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
      integer o,mm,n
      real*8  e0,e1,e2,e3,e4,e5,e6,e7
      real*8  a1,a2,a3,phi_x,phi_y,phi_z
      real*8  nx,ny,nz,length_n
      real*8  phi_ref,xd,yd,zd,c0,c1
      real*8  c00,c01,c10,c11

      mm=1

      if (1d0.eq.1d0) then
      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)

         xd=(x_prob_c(o,mm)-x1(nx_prob_C(o,mm)))/(dx1_i(1))
         yd=(y_prob_c(o,mm)-x2(ny_prob_C(o,mm)))/(dx2_j(1))
         zd=(z_prob_c(o,mm)-x3(nz_prob_C(o,mm)))/(dx3_k(1))


         c00=(1d0-xd)*
     $   phi(nx_prob_C(o,mm),ny_prob_C(o,mm),nz_prob_C(o,mm))+
     $   xd*phi(nx_prob_C(o,mm)+1,ny_prob_C(o,mm),nz_prob_C(o,mm))

         c01=(1d0-xd)*
     $   phi(nx_prob_C(o,mm),ny_prob_C(o,mm),nz_prob_C(o,mm)+1)+
     $   xd*phi(nx_prob_C(o,mm)+1,ny_prob_C(o,mm),nz_prob_C(o,mm)+1)

         c10=(1d0-xd)*
     $   phi(nx_prob_C(o,mm),ny_prob_C(o,mm)+1,nz_prob_C(o,mm))+
     $   xd*phi(nx_prob_C(o,mm)+1,ny_prob_C(o,mm)+1,nz_prob_C(o,mm))

         c11=(1d0-xd)*
     $   phi(nx_prob_C(o,mm),ny_prob_C(o,mm)+1,nz_prob_C(o,mm)+1)+
     $   xd*phi(nx_prob_C(o,mm)+1,ny_prob_C(o,mm)+1,nz_prob_C(o,mm)+1)

         c0=c00*(1-yd)+c10*yd
         c1=c01*(1-yd)+c11*yd

         phi(i,j,k)=c0*(1-zd)+c1*zd


         end do

         end if

      if (0d0.eq.1d0) then

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

      end if


      return
      end

      subroutine prep_phi
      include 'var.cmn'
      real*8 coef0_tmp

c     determine matrix coefficients normalized by minus the diagonal element
c     only valid for solution of div(c * grad(phi))=0
!$omp parallel do private(i,j,k)
      do k=0,km
      do j=0,jm
      do i=0,im
       work1(i,j,k)=(c(i,j,k,2)+c(i+1,j,k,2))/dx1_i(1)
       work2(i,j,k)=(c(i,j,k,2)+c(i,j+1,k,2))/dx2_j(1)
       work3(i,j,k)=(c(i,j,k,2)+c(i,j,k+1,2))/dx3_k(1)
      end do
      end do
      end do

!$omp parallel do private(coef0_tmp,i,j,k)
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
      real*8 phi_temp(0:im+1,0:jm+1,0:km+1)

c     Temporary array is used to avoid modifying current values with future values while current values 
c     are still neede amongs the cell centers/nodes.
      phi_temp=phi

c     Jacobi iteration step
c     Relaxation factor 1, so phi(i) at the rhs cancels out 
!$omp parallel do private(i,j,k)
      do k=1,km
      do j=1,jm
      do i=1,im
      if (typ(i,j,k).eq.1) then
       phi_temp(i,j,k)=-(coef1_phi(i,j,k)*phi(i-1,j,k)+
     $              coef2_phi(i,j,k)*phi(i+1,j,k)+
     $              coef3_phi(i,j,k)*phi(i,j-1,k)+
     $              coef4_phi(i,j,k)*phi(i,j+1,k)+
     $              coef5_phi(i,j,k)*phi(i,j,k-1)+
     $              coef6_phi(i,j,k)*phi(i,j,k+1))
      end if
      end do
      end do
      end do

      phi=phi_temp

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
!$omp parallel do private(j,k)
       do k=1,km
       do j=1,jm 
        c(im+1,j,k,n)=2d0*c_ref(n)-c(im,j,k,n)
       end do
       end do


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
      real*8  ar1,ar2,ar3,rhs(1:im,1:jm,1:km)


      if (rung.eq.1) then

         storage_c=c
      end if

      

      do n=1,nmax  
       if (0d0.eq.1d0) then    
!$omp parallel do private(i,j,k,ar1,ar2,ar3)
       do k=0,km
       do j=0,jm 
       do i=0,im

        ar1=0.5d0*abs(u1(i,j,k))*dx1(1)
        ar2=0.5d0*abs(u2(i,j,k))*dx2(1)
        ar3=0.5d0*abs(u3(i,j,k))*dx3(1)

    

        work1(i,j,k)=(dif(n)+ar1)*(c(i+1,j,k,n)-c(i,j,k,n))/dx1_i(1)
     $            -u1(i,j,k)*0.5d0*(c(i,j,k,n)+c(i+1,j,k,n)) 
        work2(i,j,k)=(dif(n)+ar2)*(c(i,j+1,k,n)-c(i,j,k,n))/dx2_j(1) 
     $            -u2(i,j,k)*0.5d0*(c(i,j,k,n)+c(i,j+1,k,n)) 
        work3(i,j,k)=(dif(n)+ar3)*(c(i,j,k+1,n)-c(i,j,k,n))/dx3_k(1) 
     $            -u3(i,j,k)*0.5d0*(c(i,j,k,n)+c(i,j,k+1,n))

       end do
       end do
       end do


!$omp parallel do private(i,j,k)
       do k=1,km
       do j=1,jm
       do i=1,im
       if (typ(i,j,k).eq.1) then



         c(i,j,k,n)=storage_c(i,j,k,n)
     $       +dtime*alpha_rung(rung)*(
     $              (work1(i,j,k)-work1(i-1,j,k))/dx1(i)+
     $              (work2(i,j,k)-work2(i,j-1,k))/dx2(j)+ 
     $              (work3(i,j,k)-work3(i,j,k-1))/dx3(k))


       end if
       end do
       end do
       end do

       end if






        if(1d0.eq.1d0) then 

!$omp parallel do private(i,j,k,ar1,ar2,ar3)
       do k=0,km
       do j=0,jm 
       do i=0,im


        work1(i,j,k)=dif(n)*(c(i+1,j,k,n)-c(i,j,k,n))/dx1_i(1)
             
        work2(i,j,k)=dif(n)*(c(i,j+1,k,n)-c(i,j,k,n))/dx2_j(1) 
            
        work3(i,j,k)=dif(n)*(c(i,j,k+1,n)-c(i,j,k,n))/dx3_k(1) 
 

       end do
       end do
       end do

!$omp parallel do private(ar1,ar2,ar3)
       do k=1,km
       do j=1,jm
       do i=1,im
       if (typ(i,j,k).eq.1) then
        rhs(i,j,k)=
     $              (work1(i,j,k)-work1(i-1,j,k))/dx1(i)+
     $              (work2(i,j,k)-work2(i,j-1,k))/dx2(j)+ 
     $              (work3(i,j,k)-work3(i,j,k-1))/dx3(k)


        

        
        ar1=0.5d0*(u1(i,j,k)+u1(i-1,j,k))
        ar2=0.5d0*(u2(i,j,k)+u2(i,j-1,k))
        ar3=0.5d0*(u3(i,j,k)+u3(i,j,k-1))






        if (ar1.gt.0d0) then
           rhs(i,j,k)=rhs(i,j,k)-ar1*(c(i,j,k,n)
     $                -c(i-1,j,k,n))/dx1_i(0)
        else
           rhs(i,j,k)=rhs(i,j,k)-ar1*(c(i+1,j,k,n)
     $                -c(i,j,k,n))/dx1_i(0)
        end if

        if (ar2.gt.0d0) then
           rhs(i,j,k)=rhs(i,j,k)-ar2*(c(i,j,k,n)
     $               -c(i,j-1,k,n))/dx2_j(0)
        else
           rhs(i,j,k)=rhs(i,j,k)-ar2*(c(i,j+1,k,n)
     $                -c(i,j,k,n))/dx2_j(0)
        end if

        if (ar3.gt.0d0) then
           rhs(i,j,k)=rhs(i,j,k)-ar3*(c(i,j,k,n)
     $               -c(i,j,k-1,n))/dx3_k(0)
        else
           rhs(i,j,k)=rhs(i,j,k)-ar3*(c(i,j,k+1,n)
     $               -c(i,j,k,n))/dx3_k(0)
        end if


        
     
       end if
       end do
       end do
       end do

     


!$omp parallel do private(i,j,k)
       do k=1,km
       do j=1,jm
       do i=1,im
       if (typ(i,j,k).eq.1) then



         c(i,j,k,n)=storage_c(i,j,k,n)
     $       +dtime*alpha_rung(rung)*rhs(i,j,k)

       end if
       end do
       end do
       end do

       end if


      end do


       balance_h2=0d0
       do k=1,km
       do j=1,jm
       do i=1,im
       if (typ(i,j,k).eq.1) then



        balance_h2=balance_h2+(c(i,j,k,1)-storage_c(i,j,k,1))
     $  /dtime*dx1_i(0)*dx2_j(0)*dx3_k(0)

       end if
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
      real*8  nx,ny,nz,length_n,beta_c
      real*8  xd,yd,zd,c00,c01,c10,c11,c0,c1

      mm=1

      if(1d0.eq.1d0) then

       
      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)

         xd=(x_prob_c(o,mm)-x1(nx_prob_C(o,mm)))/(dx1_i(1))
         yd=(y_prob_c(o,mm)-x2(ny_prob_C(o,mm)))/(dx2_j(1))
         zd=(z_prob_c(o,mm)-x3(nz_prob_C(o,mm)))/(dx3_k(1))

         do n=2,3

         c00=(1d0-xd)*
     $   c(nx_prob_C(o,mm),ny_prob_C(o,mm),nz_prob_C(o,mm),n)+
     $   xd*c(nx_prob_C(o,mm)+1,ny_prob_C(o,mm),nz_prob_C(o,mm),n)

         c01=(1d0-xd)*
     $   c(nx_prob_C(o,mm),ny_prob_C(o,mm),nz_prob_C(o,mm)+1,n)+
     $   xd*c(nx_prob_C(o,mm)+1,ny_prob_C(o,mm),nz_prob_C(o,mm)+1,n)

         c10=(1d0-xd)*
     $   c(nx_prob_C(o,mm),ny_prob_C(o,mm)+1,nz_prob_C(o,mm),n)+
     $   xd*c(nx_prob_C(o,mm)+1,ny_prob_C(o,mm)+1,nz_prob_C(o,mm),n)

         c11=(1d0-xd)*
     $   c(nx_prob_C(o,mm),ny_prob_C(o,mm)+1,nz_prob_C(o,mm)+1,n)+
     $   xd*c(nx_prob_C(o,mm)+1,ny_prob_C(o,mm)+1,nz_prob_C(o,mm)+1,n)

         c0=c00*(1-yd)+c10*yd
         c1=c01*(1-yd)+c11*yd

         c_s(n)=c0*(1-zd)+c1*zd

         end do 

         beta_c=rb(o)/dx1(1)
         c(i,j,k,1)=(1+beta_c)*c_ref(1)-(beta_c)*c_s(1)


         c(i,j,k,2)=c_s(2)
         c(i,j,k,3)=c_s(3)

         if (radius.ge.R_max.or.m1.eq.m1_max) then
         c_KOH_surf(o)=c_s(2)
         end if



         end do

         

         end if


      if (0d0.eq.1d0) then
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

         end do



         

         c_s(1)=c_ref(1)


      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)




         do n=1,1
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

         end if
      return
      end



      subroutine bulb_p
      include 'var.cmn' 
      integer o
      real*8 p_temp(0:im+1,0:jm+1,0:km+1)

      p_temp=p
      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)

         if (direction(o) .eq. 1) then
            p_temp(i,j,k)=p(i+1,j,k)

         else if (direction(o) .eq. -1) then
            p_temp(i,j,k)=p(i-1,j,k)

         else if (direction(o) .eq. 2) then
            p_temp(i,j,k)=p(i,j+1,k)

         else if (direction(o) .eq. -2) then
            p_temp(i,j,k)=p(i,j-1,k)

         else if (direction(o) .eq. 3) then
            p_temp(i,j,k)=p(i,j,k+1)

         else if (direction(o) .eq. -3) then
            p_temp(i,j,k)=p(i,j,k-1)

         end if

      end do
      p=p_temp

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
     $       A31(1:im,1:jm),B31(1:im,1:jm),
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


         A31(i,j)=-(-0.5d0*(c(i,j,1,1)+c(i,j,0,1))
     $           *u3(i,j,0)*dx1(i)*dx2(j))
         B31(i,j)=-0.5d0*(c(i,j,km+1,1)+c(i,j,km,1))
     $           *u3(i,j,km)*dx1(i)*dx2(j)

         A32(i,j)=-(dif(1)*(c(i,j,1,1)-c(i,j,0,1))/dx3_k(0)
     $           *dx1(i)*dx2(j))
         B32(i,j)=dif(1)*(c(i,j,km+1,1)-c(i,j,km,1))/dx3_k(km)
     $           *dx1(i)*dx2(j)


      end do
      end do

      term5=sum(A31)+sum(A32)
      term6=sum(B31)+sum(B32)

      write(*,*)sum(A31),sum(B31)
      write(*,*)sum(A32),sum(B32)

      Error=abs(term1+term2+term3+term4+term5+term6+total-balance_h2)


      if (time.lt.0.5d0*dtime) then
       open(unit=32,file='H2_balance.dat')
      else 
       open(unit=32,file='H2_balance.dat',position='append')
      end if

 30   format(e14.6,1i9,1000e36.24)

      write(6,30)time,m1,term1,term2,
     $            term3,term4,term5,term6,total,-balance_h2,Error

      write(32,30)time,m1,term1,term2,
     $            term3,term4,term5,term6,total,-balance_h2,Error
      close(32)


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



      subroutine finding_u_IB
      include 'var.cmn' 
      integer o

      typ_u_IB(0:im+1,0:jm+1,0:km+1)=0

      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)

          
         if (typ(i+1,j,k).eq.1) then

            typ_u_IB(i,j,k)=1

         else

          if (typ_IB(i+1,j,k).eq.1) then

             typ_u_IB(i,j,k)=1
           
          end if

         end if


         if (typ(i-1,j,k).eq.1) then

            typ_u_IB(i-1,j,k)=1

         else

          if (typ_IB(i-1,j,k).eq.1) then

             typ_u_IB(i-1,j,k)=1
           
          end if


         end if

      end do

      
      o_u=0
      do k=0,km+1
      do j=0,jm+1
      do i=0,im+1

         if (typ_u_IB(i,j,k).eq.1) then

          o_u=o_u+1

          x_IB_U(o_u)=i
          y_IB_U(o_u)=j
          z_IB_U(o_u)=k

         end if

      end do
      end do
      end do


      return
      end



      subroutine probe_IB_u
      include 'var.cmn' 
      integer o,mm
      real*8  rnode_U

c        Neumann BC in Spherical  coordinates

c        phiS = angle phi in spherical coordinates.
c        not to be confused with phi the potential
c        if x=0 then atan(y,x)=pi/2. This is not true however,
c        if y is negative, then it should be 3pi/2.

      mm=1
      
      do o=1,o_u
         i=x_IB_U(o)
         j=y_IB_U(o)
         k=z_IB_U(o)

         phiS_U(o)=datan2((x2(j)-yc),(x1_i(i)-xc))

         thetaS_U(o)=acos((x3(k)-zc)/sqrt((x1_i(i)-xc)**2+
     $    (x2(j)-yc)**2+(x3(k)-zc)**2))


         T_ni_U(o)=cos(phiS_U(o))*sin(thetaS_U(o))
         T_nj_U(o)=sin(phiS_U(o))*sin(thetaS_U(o))
         T_nk_U(o)=cos(thetaS_U(o))

         T_ti_U(o)=-sin(phiS_U(o))
         T_tj_U(o)=cos(phiS_U(o))
         T_tk_U(o)=0d0

         T_si_U(o)=cos(phiS_U(o))*cos(thetaS_U(o))
         T_sj_U(o)=sin(phiS_U(o))*cos(thetaS_U(o))
         T_sk_U(o)=-sin(thetaS_U(o))

c        initialise x-direction
c        Normal distance between interface and bubble centre
         rnode_U=sqrt((x1_i(i)-xc)**2+(x2(j)-yc)**2+(x3(k)-zc)**2)


c        normal distance of interface to bubble surface 
c        (can be negative, as it can be outside of the bubble.)
         rb_U(o)=radius-rnode_U

c        Intersection with the bubble interface and the normal u_IB velocity point
         xb_U(o)=x1_i(i)+rb_U(o)*T_ni_U(o)
         yb_U(o)=x2(j)  +rb_U(o)*T_nj_U(o)
         zb_U(o)=x3(k)  +rb_U(o)*T_nk_U(o)

c        Obtain probe point with distance dx1 normal to the bubble interface
         x_prob_U(o,mm)=xb_U(o)+dx1_i(0)*T_ni_U(o)
         y_prob_U(o,mm)=yb_U(o)+dx1_i(0)*T_nj_U(o)
         z_prob_U(o,mm)=zb_U(o)+dx1_i(0)*T_nk_U(o)

c        Marangoni flow part (skip)
         x_prob_mt(o,1)=xb_U(o)+dx1_i(0)*(T_ti_U(o)+dd*T_ni_U(o))
         y_prob_mt(o,1)=yb_U(o)+dx1_i(0)*(T_tj_U(o)+dd*T_nj_U(o))
         z_prob_mt(o,1)=zb_U(o)+dx1_i(0)*(T_tk_U(o)+dd*T_nk_U(o))

         x_prob_mt(o,2)=xb_U(o)-dx1_i(0)*(T_ti_U(o)-dd*T_ni_U(o))
         y_prob_mt(o,2)=yb_U(o)-dx1_i(0)*(T_tj_U(o)-dd*T_nj_U(o))
         z_prob_mt(o,2)=zb_U(o)-dx1_i(0)*(T_tk_U(o)-dd*T_nk_U(o))



         x_prob_ms(o,1)=xb_U(o)+dx1_i(0)*(T_si_U(o)+dd*T_ni_U(o))
         y_prob_ms(o,1)=yb_U(o)+dx1_i(0)*(T_sj_U(o)+dd*T_nj_U(o))
         z_prob_ms(o,1)=zb_U(o)+dx1_i(0)*(T_sk_U(o)+dd*T_nk_U(o))

         x_prob_ms(o,2)=xb_U(o)-dx1_i(0)*(T_si_U(o)-dd*T_ni_U(o))
         y_prob_ms(o,2)=yb_U(o)-dx1_i(0)*(T_sj_U(o)-dd*T_nj_U(o))
         z_prob_ms(o,2)=zb_U(o)-dx1_i(0)*(T_sk_U(o)-dd*T_nk_U(o))


c        Trilinear interpolation of the probe point. Determin i,j,k position. 
c        prob number of cells starts at an interface so only 
c        +0.5d if that direction is on the level of a node
         nx_prob_UU(o,mm)=floor(x_prob_U(o,mm)/dx1_i(0))
         ny_prob_UU(o,mm)=floor(y_prob_U(o,mm)/dx2_j(0)+0.5d0)
         nz_prob_UU(o,mm)=floor(z_prob_U(o,mm)/dx3_k(0)+0.5d0)

         nx_prob_UV(o,mm)=floor(x_prob_U(o,mm)/dx1_i(0)+0.5d0)
         ny_prob_UV(o,mm)=floor(y_prob_U(o,mm)/dx2_j(0))
         nz_prob_UV(o,mm)=floor(z_prob_U(o,mm)/dx3_k(0)+0.5d0)

         nx_prob_UW(o,mm)=floor(x_prob_U(o,mm)/dx1_i(0)+0.5d0)
         ny_prob_UW(o,mm)=floor(y_prob_U(o,mm)/dx2_j(0)+0.5d0)
         nz_prob_UW(o,mm)=floor(z_prob_U(o,mm)/dx3_k(0))

         nx_prob_mt(o,1)=floor(x_prob_mt(o,1)/dx1_i(0)+0.5d0)
         ny_prob_mt(o,1)=floor(y_prob_mt(o,1)/dx2_j(0)+0.5d0)
         nz_prob_mt(o,1)=floor(z_prob_mt(o,1)/dx3_k(0)+0.5d0)

         nx_prob_mt(o,2)=floor(x_prob_mt(o,2)/dx1_i(0)+0.5d0)
         ny_prob_mt(o,2)=floor(y_prob_mt(o,2)/dx2_j(0)+0.5d0)
         nz_prob_mt(o,2)=floor(z_prob_mt(o,2)/dx3_k(0)+0.5d0)

         nx_prob_ms(o,1)=floor(x_prob_ms(o,1)/dx1_i(0)+0.5d0)
         ny_prob_ms(o,1)=floor(y_prob_ms(o,1)/dx2_j(0)+0.5d0)
         nz_prob_ms(o,1)=floor(z_prob_ms(o,1)/dx3_k(0)+0.5d0)

         nx_prob_ms(o,2)=floor(x_prob_ms(o,2)/dx1_i(0)+0.5d0)
         ny_prob_ms(o,2)=floor(y_prob_ms(o,2)/dx2_j(0)+0.5d0)
         nz_prob_ms(o,2)=floor(z_prob_ms(o,2)/dx3_k(0)+0.5d0)


      end do


      return
      end

      subroutine bulb_u
      include 'var.cmn' 
      integer o,mm,test,nn,ii,jj,kk
      integer nx_prob_e,ny_prob_e,nz_prob_e
c     Solving u1  on IB interface

      real*8 u1i_ref,u2i_ref,u3i_ref
      real*8 uni_ref,uti_ref,usi_ref
      real*8 u_star,v_star,w_star
      real*8 unis,utis,usis
      real*8 uni,uti,usi
      real*8 c_tt(1:2),c_ss(1:2)
      real*8 constant_t,constant_s
      real*8 R_IB_U,R_P_U



      real*8  xd,yd,zd,beta_U
      real*8  c00,c10,c01,c11,c0,c1

      real*8  phiS_e,thetaS_e
      real*8  T_ni_e,T_nj_e,T_nk_e
      real*8  rnode_e,rb_e
      real*8  xb_e,yb_e,zb_e,u1_temp(0:im+1,0:jm+1,0:km+1)
      real*8  x_prob_e,y_prob_e,z_prob_e

      mm=1
      u1_temp=u1
      do o=1,o_u

          i=x_IB_U(o)
          j=y_IB_U(o)
          k=z_IB_U(o)

c     Marangoni Flow (skip)
      do nn=1,2

      do kk=nz_prob_mt(o,nn),nz_prob_mt(o,nn)+1
      do jj=ny_prob_mt(o,nn),ny_prob_mt(o,nn)+1
      do ii=nx_prob_mt(o,nn),nx_prob_mt(o,nn)+1

        if(typ(ii,jj,kk)+typ_IB(ii,jj,kk).eq.0.and.ii.ne.0) then


             phiS_e=datan2((x2(jj)-yc),(x1(ii)-xc))
  
             thetaS_e=acos((x3(kk)-zc)/sqrt((x1(ii)-xc)**2+
     $                (x2(jj)-yc)**2+(x3(kk)-zc)**2))

             T_ni_e=cos(phiS_e)*sin(thetaS_e)
             T_nj_e=sin(phiS_e)*sin(thetaS_e)
             T_nk_e=cos(thetaS_e)

c            initialise x-direction
c            Normal distance between interface and bubble centre
         rnode_e=sqrt((x1(ii)-xc)**2+(x2(jj)-yc)**2
     $       +(x3(kk)-zc)**2)

c            normal distance of interface to bubble surface 
c            (can be negative)
             rb_e=radius-rnode_e

             xb_e=x1(ii)+rb_e*T_ni_e
             yb_e=x2(jj)+rb_e*T_nj_e
             zb_e=x3(kk)+rb_e*T_nk_e

             x_prob_e=xb_e+dx1_i(0)*T_ni_e
             y_prob_e=yb_e+dx1_i(0)*T_nj_e
             z_prob_e=zb_e+dx1_i(0)*T_nk_e

             nx_prob_e=floor(x_prob_e/dx1_i(0)+0.5d0)
             ny_prob_e=floor(y_prob_e/dx2_j(0)+0.5d0)
             nz_prob_e=floor(z_prob_e/dx3_k(0)+0.5d0)

             xd=(x_prob_e-x1(nx_prob_e))/(dx1(1))
             yd=(y_prob_e-x2(ny_prob_e))/(dx2(1))
             zd=(z_prob_e-x3(nz_prob_e))/(dx3(1))

             c00=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e,2)

             c01=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e+1,2)

             c10=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1,nz_prob_e,2)
 
             c11=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1
     $       ,nz_prob_e+1,2)

             c0=c00*(1-yd)+c10*yd
             c1=c01*(1-yd)+c11*yd

             

             c(ii,jj,kk,2)=c0*(1-zd)+c1*zd


        end if

           
      end do
      end do
      end do

      

      do kk=nz_prob_ms(o,nn),nz_prob_ms(o,nn)+1
      do jj=ny_prob_ms(o,nn),ny_prob_ms(o,nn)+1
      do ii=nx_prob_ms(o,nn),nx_prob_ms(o,nn)+1

        if(typ(ii,jj,kk)+typ_IB(ii,jj,kk).eq.0.and.ii.ne.0) then



             phiS_e=datan2((x2(jj)-yc),(x1(ii)-xc))
  
             thetaS_e=acos((x3(kk)-zc)/sqrt((x1(ii)-xc)**2+
     $                (x2(jj)-yc)**2+(x3(kk)-zc)**2))

             T_ni_e=cos(phiS_e)*sin(thetaS_e)
             T_nj_e=sin(phiS_e)*sin(thetaS_e)
             T_nk_e=cos(thetaS_e)

c            initialise x-direction
c            Normal distance between interface and bubble centre
         rnode_e=sqrt((x1(ii)-xc)**2+(x2(jj)-yc)**2
     $       +(x3(kk)-zc)**2)

c            normal distance of interface to bubble surface 
c            (can be negative)
             rb_e=radius-rnode_e

             xb_e=x1(ii)+rb_e*T_ni_e
             yb_e=x2(jj)+rb_e*T_nj_e
             zb_e=x3(kk)+rb_e*T_nk_e

             x_prob_e=xb_e+dx1_i(0)*T_ni_e
             y_prob_e=yb_e+dx1_i(0)*T_nj_e
             z_prob_e=zb_e+dx1_i(0)*T_nk_e

             nx_prob_e=floor(x_prob_e/dx1_i(0)+0.5d0)
             ny_prob_e=floor(y_prob_e/dx2_j(0)+0.5d0)
             nz_prob_e=floor(z_prob_e/dx3_k(0)+0.5d0)

             xd=(x_prob_e-x1(nx_prob_e))/(dx1(1))
             yd=(y_prob_e-x2(ny_prob_e))/(dx2(1))
             zd=(z_prob_e-x3(nz_prob_e))/(dx3(1))

             c00=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e,2)

             c01=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e+1,2)

             c10=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1,nz_prob_e,2)
 
             c11=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1
     $       ,nz_prob_e+1,2)

             c0=c00*(1-yd)+c10*yd
             c1=c01*(1-yd)+c11*yd

             c(ii,jj,kk,2)=c0*(1-zd)+c1*zd



        end if

           
      end do
      end do
      end do

      end do

c     Trilinear interpolation for u_velocity
c     x-component u1


      xd=(x_prob_U(o,mm)-x1_i(nx_prob_UU(o,mm)))/(dx1(1))
      yd=(y_prob_U(o,mm)-x2(ny_prob_UU(o,mm)))/(dx2_j(1))
      zd=(z_prob_U(o,mm)-x3(nz_prob_UU(o,mm)))/(dx3_k(1))


      c00=(1d0-xd)*
     $  u1(nx_prob_UU(o,mm),ny_prob_UU(o,mm),nz_prob_UU(o,mm))+
     $  xd*u1(nx_prob_UU(o,mm)+1,ny_prob_UU(o,mm),nz_prob_UU(o,mm))

      c01=(1d0-xd)*
     $  u1(nx_prob_UU(o,mm),ny_prob_UU(o,mm),nz_prob_UU(o,mm)+1)+
     $  xd*u1(nx_prob_UU(o,mm)+1,ny_prob_UU(o,mm),nz_prob_UU(o,mm)+1)

      c10=(1d0-xd)*
     $  u1(nx_prob_UU(o,mm),ny_prob_UU(o,mm)+1,nz_prob_UU(o,mm))+
     $  xd*u1(nx_prob_UU(o,mm)+1,ny_prob_UU(o,mm)+1,nz_prob_UU(o,mm))

      c11=(1d0-xd)*
     $  u1(nx_prob_UU(o,mm),ny_prob_UU(o,mm)+1,nz_prob_UU(o,mm)+1)+
     $  xd*u1(nx_prob_UU(o,mm)+1,ny_prob_UU(o,mm)+1,nz_prob_UU(o,mm)+1)

      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      u1i_ref=c0*(1-zd)+c1*zd


c     y-component u1

      xd=(x_prob_U(o,mm)-x1(nx_prob_UV(o,mm)))/(dx1_i(1))
      yd=(y_prob_U(o,mm)-x2_j(ny_prob_UV(o,mm)))/(dx2(1))
      zd=(z_prob_U(o,mm)-x3(nz_prob_UV(o,mm)))/(dx3_k(1))

      c00=(1d0-xd)*
     $  u2(nx_prob_UV(o,mm),ny_prob_UV(o,mm),nz_prob_UV(o,mm))+
     $  xd*u2(nx_prob_UV(o,mm)+1,ny_prob_UV(o,mm),nz_prob_UV(o,mm))

      c01=(1d0-xd)*
     $  u2(nx_prob_UV(o,mm),ny_prob_UV(o,mm),nz_prob_UV(o,mm)+1)+
     $  xd*u2(nx_prob_UV(o,mm)+1,ny_prob_UV(o,mm),nz_prob_UV(o,mm)+1)

      c10=(1d0-xd)*
     $  u2(nx_prob_UV(o,mm),ny_prob_UV(o,mm)+1,nz_prob_UV(o,mm))+
     $  xd*u2(nx_prob_UV(o,mm)+1,ny_prob_UV(o,mm)+1,nz_prob_UV(o,mm))

      c11=(1d0-xd)*
     $  u2(nx_prob_UV(o,mm),ny_prob_UV(o,mm)+1,nz_prob_UV(o,mm)+1)+
     $  xd*u2(nx_prob_UV(o,mm)+1,ny_prob_UV(o,mm)+1,nz_prob_UV(o,mm)+1)

      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      u2i_ref=c0*(1-zd)+c1*zd

c     z-component u1

      xd=(x_prob_U(o,mm)-x1(nx_prob_UW(o,mm)))/(dx1_i(1))
      yd=(y_prob_U(o,mm)-x2(ny_prob_UW(o,mm)))/(dx2_j(1))
      zd=(z_prob_U(o,mm)-x3_k(nz_prob_UW(o,mm)))/(dx3(1))

      c00=(1d0-xd)*
     $  u3(nx_prob_UW(o,mm),ny_prob_UW(o,mm),nz_prob_UW(o,mm))+
     $  xd*u3(nx_prob_UW(o,mm)+1,ny_prob_UW(o,mm),nz_prob_UW(o,mm))

      c01=(1d0-xd)*
     $  u3(nx_prob_UW(o,mm),ny_prob_UW(o,mm),nz_prob_UW(o,mm)+1)+
     $  xd*u3(nx_prob_UW(o,mm)+1,ny_prob_UW(o,mm),nz_prob_UW(o,mm)+1)

      c10=(1d0-xd)*
     $  u3(nx_prob_UW(o,mm),ny_prob_UW(o,mm)+1,nz_prob_UW(o,mm))+
     $  xd*u3(nx_prob_UW(o,mm)+1,ny_prob_UW(o,mm)+1,nz_prob_UW(o,mm))

      c11=(1d0-xd)*
     $  u3(nx_prob_UW(o,mm),ny_prob_UW(o,mm)+1,nz_prob_UW(o,mm)+1)+
     $  xd*u3(nx_prob_UW(o,mm)+1,ny_prob_UW(o,mm)+1,nz_prob_UW(o,mm)+1)

      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      u3i_ref=c0*(1-zd)+c1*zd

c     Marangoni flow (skip)
      do nn=1,2

         xd=(x_prob_mt(o,nn)-x1(nx_prob_mt(o,nn)))/(dx1_i(1))
         yd=(y_prob_mt(o,nn)-x2(ny_prob_mt(o,nn)))/(dx2_j(1))
         zd=(z_prob_mt(o,nn)-x3(nz_prob_mt(o,nn)))/(dx3_k(1))

         c00=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn),nz_prob_mt(o,nn),2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn),nz_prob_mt(o,nn),2)

         c01=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn),nz_prob_mt(o,nn)+1,2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn),nz_prob_mt(o,nn)+1,2)

         c10=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn)+1,nz_prob_mt(o,nn),2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn)+1,nz_prob_mt(o,nn),2)

         c11=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn)+1,nz_prob_mt(o,nn)+1,2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn)+1
     $   ,nz_prob_mt(o,nn)+1,2)

         c0=c00*(1-yd)+c10*yd
         c1=c01*(1-yd)+c11*yd

         c_tt(nn)=c0*(1-zd)+c1*zd

         xd=(x_prob_ms(o,nn)-x1(nx_prob_ms(o,nn)))/(dx1_i(1))
         yd=(y_prob_ms(o,nn)-x2(ny_prob_ms(o,nn)))/(dx2_j(1))
         zd=(z_prob_ms(o,nn)-x3(nz_prob_ms(o,nn)))/(dx3_k(1))

         c00=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn),nz_prob_ms(o,nn),2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn),nz_prob_ms(o,nn),2)

         c01=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn),nz_prob_ms(o,nn)+1,2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn),nz_prob_ms(o,nn)+1,2)

         c10=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn),2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn),2)

         c11=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn)+1,2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn)+1
     $   ,nz_prob_ms(o,nn)+1,2)

         c0=c00*(1-yd)+c10*yd
         c1=c01*(1-yd)+c11*yd

         c_ss(nn)=c0*(1-zd)+c1*zd

      end do


      constant_t=(c_tt(1)-c_tt(2))/(2d0*dx1_i(0))
     $           *dsigma/(ro*nu)
      constant_s=(c_ss(1)-c_ss(2))/(2d0*dx1_i(0))
     $           *dsigma/(ro*nu)



c     Translate Cartesian velocities to spherical velocities
c     Look at bulb.f for specification on each T (conversion of x,y,z to n,t,s of the 
c     probe point with respect tot the bubble center)
        uni_ref=u1i_ref*T_ni_U(o)+
     $          u2i_ref*T_nj_U(o)+
     $          u3i_ref*T_nk_U(o)
        uti_ref=u1i_ref*T_ti_U(o)+
     $          u2i_ref*T_tj_U(o)+
     $          u3i_ref*T_tk_U(o)
        usi_ref=u1i_ref*T_si_U(o)+
     $          u2i_ref*T_sj_U(o)+
     $          u3i_ref*T_sk_U(o)

 
c     Dirichlet BC for normal direction, beta can be negative as rb_i can be negative
      beta_U=rb_U(o)/dx1(1)
c     Velocity at the bubble surface in cartisian coordinate (the bubble grows and
c     therefore this cannot be set to zero (no-slip B.C))
      u_star=dxcdt+drdt*T_ni_U(o)
      v_star=      drdt*T_nj_U(o)
      w_star=      drdt*T_nk_U(o)

c     Obtain bubble surface velocity in spherical coordinates
        unis=u_star*T_ni_U(o)+
     $       v_star*T_nj_U(o)+
     $       w_star*T_nk_U(o)
        utis=u_star*T_ti_U(o)+
     $       v_star*T_tj_U(o)+
     $       w_star*T_tk_U(o)
        usis=u_star*T_si_U(o)+
     $       v_star*T_sj_U(o)+
     $       w_star*T_sk_U(o)


      R_IB_U=radius-rb_U(o)
      R_P_U=radius+dx1(1)

c     Linear interpolation for obtaining all velocities in spherical coordinates 
c     at the IB_u velocity location 
      uni=(1+beta_U)*unis-(beta_U)*uni_ref
      uti=(1+beta_U)*utis-(beta_U)*uti_ref
      usi=(1+beta_U)*usis-(beta_U)*usi_ref

  
c     Solve u1


c      u1(i,j,k)=uni*T_ni_U(o)+
c     $          uti*T_ti_U(o)+usi*T_si_U(o)
  
c     Convert from spherical coordiante to carthesian coordinate
      u1_temp(i,j,k)=cos(phiS_U(o))*sin(thetaS_U(o))*uni-
     $                           sin(phiS_U(o))*uti+
     $          cos(phiS_U(o))*cos(thetaS_U(o))*usi

c     not part of the computation
      if (radius.ge.R_max.or.m1.eq.m1_max) then
      un_surf_u(o)=unis
      ut_surf_u(o)=utis
      us_surf_u(o)=usis
      dc_t_u(o)=(c_tt(1)-c_tt(2))/(2d0*dx1_i(0))
      dc_s_u(o)=(c_ss(1)-c_ss(2))/(2d0*dx1_i(0))
      end if


      end do
      u1=u1_temp

      return
      end


      subroutine finding_v_IB
      include 'var.cmn' 
      integer o

      typ_v_IB(0:im+1,0:jm+1,0:km+1)=0

c     if cell centre is inside of the bubble typ(i,j,k)=0
c     if cell centre is outside of the bubble typ(i,j,k)=1
c     if typ_IB(i,j,k) is 1 means it is IB cell


      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)


          
         if (typ(i,j+1,k).eq.1) then

            typ_v_IB(i,j,k)=1

          else

           if(typ_IB(i,j+1,k).eq.1) then

             typ_v_IB(i,j,k)=1


           end if


         end if


         if (typ(i,j-1,k).eq.1) then

            typ_v_IB(i,j-1,k)=1

          else

           if(typ_IB(i,j-1,k).eq.1) then

             typ_v_IB(i,j-1,k)=1



           end if

         end if


      end do


      o_v=0
      do k=0,km+1
      do j=0,jm+1
      do i=0,im+1

         if (typ_v_IB(i,j,k).eq.1) then

          o_v=o_v+1

          x_IB_V(o_v)=i
          y_IB_V(o_v)=j
          z_IB_V(o_v)=k

         end if

      end do
      end do
      end do


      return
      end


      subroutine probe_IB_V
      include 'var.cmn' 
      integer o,mm
      real*8  rnode_V

c        Neumann BC in Spherical  coordinates

c        phiS = angle phi in spherical coordinates.
c        not to be confused with phi the potential
c        if x=0 then atan(y,x)=pi/2. This is not true however,
c        if y is negative, then it should be 3pi/2.

      mm=1
      
      do o=1,o_v
         i=x_IB_V(o)
         j=y_IB_V(o)
         k=z_IB_V(o)

         phiS_V(o)=datan2((x2_j(j)-yc),(x1(i)-xc))

         thetaS_V(o)=acos((x3(k)-zc)/sqrt((x1(i)-xc)**2+
     $    (x2_j(j)-yc)**2+(x3(k)-zc)**2))


         T_ni_V(o)=cos(phiS_V(o))*sin(thetaS_V(o))
         T_nj_V(o)=sin(phiS_V(o))*sin(thetaS_V(o))
         T_nk_V(o)=cos(thetaS_V(o))

         T_ti_V(o)=-sin(phiS_V(o))
         T_tj_V(o)=cos(phiS_V(o))
         T_tk_V(o)=0d0

         T_si_V(o)=cos(phiS_V(o))*cos(thetaS_V(o))
         T_sj_V(o)=sin(phiS_V(o))*cos(thetaS_V(o))
         T_sk_V(o)=-sin(thetaS_V(o))

c        initialise x-direction
c        Normal distance between interface and bubble centre
         rnode_V=sqrt((x1(i)-xc)**2+(x2_j(j)-yc)**2+(x3(k)-zc)**2)


c        normal distance of interface to bubble surface 
c        (can be negative)
         rb_V(o)=radius-rnode_V

         xb_V(o)=x1(i)  +rb_V(o)*T_ni_V(o)
         yb_V(o)=x2_j(j)+rb_V(o)*T_nj_V(o)
         zb_V(o)=x3(k)  +rb_V(o)*T_nk_V(o)

         x_prob_V(o,mm)=xb_V(o)+dx1_i(0)*T_ni_V(o)
         y_prob_V(o,mm)=yb_V(o)+dx1_i(0)*T_nj_V(o)
         z_prob_V(o,mm)=zb_V(o)+dx1_i(0)*T_nk_V(o)

         x_prob_mt(o,1)=xb_V(o)+dx1_i(0)*(T_ti_V(o)+dd*T_ni_V(o))
         y_prob_mt(o,1)=yb_V(o)+dx1_i(0)*(T_tj_V(o)+dd*T_nj_V(o))
         z_prob_mt(o,1)=zb_V(o)+dx1_i(0)*(T_tk_V(o)+dd*T_nk_V(o))

         x_prob_mt(o,2)=xb_V(o)-dx1_i(0)*(T_ti_V(o)-dd*T_ni_V(o))
         y_prob_mt(o,2)=yb_V(o)-dx1_i(0)*(T_tj_V(o)-dd*T_nj_V(o))
         z_prob_mt(o,2)=zb_V(o)-dx1_i(0)*(T_tk_V(o)-dd*T_nk_V(o))

         x_prob_ms(o,1)=xb_V(o)+dx1_i(0)*(T_si_V(o)+dd*T_ni_V(o))
         y_prob_ms(o,1)=yb_V(o)+dx1_i(0)*(T_sj_V(o)+dd*T_nj_V(o))
         z_prob_ms(o,1)=zb_V(o)+dx1_i(0)*(T_sk_V(o)+dd*T_nk_V(o))

         x_prob_ms(o,2)=xb_V(o)-dx1_i(0)*(T_si_V(o)-dd*T_ni_V(o))
         y_prob_ms(o,2)=yb_V(o)-dx1_i(0)*(T_sj_V(o)-dd*T_nj_V(o))
         z_prob_ms(o,2)=zb_V(o)-dx1_i(0)*(T_sk_V(o)-dd*T_nk_V(o))



c        prob number of cells starts at an interface so only 
c        +0.5d if that direction is on the level of a node
         nx_prob_VU(o,mm)=floor(x_prob_V(o,mm)/dx1_i(0))
         ny_prob_VU(o,mm)=floor(y_prob_V(o,mm)/dx2_j(0)+0.5d0)
         nz_prob_VU(o,mm)=floor(z_prob_V(o,mm)/dx3_k(0)+0.5d0)

         nx_prob_VV(o,mm)=floor(x_prob_V(o,mm)/dx1_i(0)+0.5d0)
         ny_prob_VV(o,mm)=floor(y_prob_V(o,mm)/dx2_j(0))
         nz_prob_VV(o,mm)=floor(z_prob_V(o,mm)/dx3_k(0)+0.5d0)

         nx_prob_VW(o,mm)=floor(x_prob_V(o,mm)/dx1_i(0)+0.5d0)
         ny_prob_VW(o,mm)=floor(y_prob_V(o,mm)/dx2_j(0)+0.5d0)
         nz_prob_VW(o,mm)=floor(z_prob_V(o,mm)/dx3_k(0))

         nx_prob_mt(o,1)=floor(x_prob_mt(o,1)/dx1_i(0)+0.5d0)
         ny_prob_mt(o,1)=floor(y_prob_mt(o,1)/dx2_j(0)+0.5d0)
         nz_prob_mt(o,1)=floor(z_prob_mt(o,1)/dx3_k(0)+0.5d0)

         nx_prob_mt(o,2)=floor(x_prob_mt(o,2)/dx1_i(0)+0.5d0)
         ny_prob_mt(o,2)=floor(y_prob_mt(o,2)/dx2_j(0)+0.5d0)
         nz_prob_mt(o,2)=floor(z_prob_mt(o,2)/dx3_k(0)+0.5d0)

         nx_prob_ms(o,1)=floor(x_prob_ms(o,1)/dx1_i(0)+0.5d0)
         ny_prob_ms(o,1)=floor(y_prob_ms(o,1)/dx2_j(0)+0.5d0)
         nz_prob_ms(o,1)=floor(z_prob_ms(o,1)/dx3_k(0)+0.5d0)

         nx_prob_ms(o,2)=floor(x_prob_ms(o,2)/dx1_i(0)+0.5d0)
         ny_prob_ms(o,2)=floor(y_prob_ms(o,2)/dx2_j(0)+0.5d0)
         nz_prob_ms(o,2)=floor(z_prob_ms(o,2)/dx3_k(0)+0.5d0)




      end do


      return
      end



      subroutine bulb_v
      include 'var.cmn' 
      integer o,mm,nn,ii,jj,kk
      integer nx_prob_e,ny_prob_e,nz_prob_e
c     Solving u2 on IB interface
      real*8 u1j_ref,u2j_ref,u3j_ref
      real*8 unj_ref,utj_ref,usj_ref
      real*8  xd,yd,zd
      real*8  c00,c10,c01,c11,c0,c1
      real*8 u_star,v_star,w_star
      real*8 unjs,utjs,usjs
      real*8 unj,utj,usj,beta_V

      real*8 R_IB_V,R_P_V

      real*8 c_tt(1:2),c_ss(1:2)
      real*8 constant_t,constant_s

      real*8  phiS_e,thetaS_e
      real*8  T_ni_e,T_nj_e,T_nk_e
      real*8  rnode_e,rb_e,u2_temp(0:im+1,0:jm+1,0:km+1)
      real*8  xb_e,yb_e,zb_e
      real*8  x_prob_e,y_prob_e,z_prob_e




      u2_temp=u2
      mm=1
      do o=1,o_v
          i=x_IB_V(o)
          j=y_IB_V(o)
          k=z_IB_V(o)


          do nn=1,2

          do kk=nz_prob_mt(o,nn),nz_prob_mt(o,nn)+1
          do jj=ny_prob_mt(o,nn),ny_prob_mt(o,nn)+1
          do ii=nx_prob_mt(o,nn),nx_prob_mt(o,nn)+1

          if(typ(ii,jj,kk)+typ_IB(ii,jj,kk).eq.0.and.ii.ne.0) then


             phiS_e=datan2((x2(jj)-yc),(x1(ii)-xc))
  
             thetaS_e=acos((x3(kk)-zc)/sqrt((x1(ii)-xc)**2+
     $                (x2(jj)-yc)**2+(x3(kk)-zc)**2))

             T_ni_e=cos(phiS_e)*sin(thetaS_e)
             T_nj_e=sin(phiS_e)*sin(thetaS_e)
             T_nk_e=cos(thetaS_e)

c            initialise x-direction
c            Normal distance between interface and bubble centre
         rnode_e=sqrt((x1(ii)-xc)**2+(x2(jj)-yc)**2
     $       +(x3(kk)-zc)**2)


c            normal distance of interface to bubble surface 
c            (can be negative)
             rb_e=radius-rnode_e

             xb_e=x1(ii)+rb_e*T_ni_e
             yb_e=x2(jj)+rb_e*T_nj_e
             zb_e=x3(kk)+rb_e*T_nk_e

             x_prob_e=xb_e+dx1_i(0)*T_ni_e
             y_prob_e=yb_e+dx1_i(0)*T_nj_e
             z_prob_e=zb_e+dx1_i(0)*T_nk_e

             nx_prob_e=floor(x_prob_e/dx1_i(0)+0.5d0)
             ny_prob_e=floor(y_prob_e/dx2_j(0)+0.5d0)
             nz_prob_e=floor(z_prob_e/dx3_k(0)+0.5d0)

             xd=(x_prob_e-x1(nx_prob_e))/(dx1(1))
             yd=(y_prob_e-x2(ny_prob_e))/(dx2(1))
             zd=(z_prob_e-x3(nz_prob_e))/(dx3(1))


             c00=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e,2)

             c01=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e+1,2)

             c10=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1,nz_prob_e,2)
 
             c11=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1
     $       ,nz_prob_e+1,2)

             c0=c00*(1-yd)+c10*yd
             c1=c01*(1-yd)+c11*yd
 
             c(ii,jj,kk,2)=c0*(1-zd)+c1*zd


          end if

           
          end do
          end do
          end do

          do kk=nz_prob_ms(o,nn),nz_prob_ms(o,nn)+1
          do jj=ny_prob_ms(o,nn),ny_prob_ms(o,nn)+1
          do ii=nx_prob_ms(o,nn),nx_prob_ms(o,nn)+1

          if(typ(ii,jj,kk)+typ_IB(ii,jj,kk).eq.0.and.ii.ne.0) then


             phiS_e=datan2((x2(jj)-yc),(x1(ii)-xc))
  
             thetaS_e=acos((x3(kk)-zc)/sqrt((x1(ii)-xc)**2+
     $                (x2(jj)-yc)**2+(x3(kk)-zc)**2))

             T_ni_e=cos(phiS_e)*sin(thetaS_e)
             T_nj_e=sin(phiS_e)*sin(thetaS_e)
             T_nk_e=cos(thetaS_e)

c            initialise x-direction
c            Normal distance between interface and bubble centre
         rnode_e=sqrt((x1(ii)-xc)**2+(x2(jj)-yc)**2
     $       +(x3(kk)-zc)**2)


c            normal distance of interface to bubble surface 
c            (can be negative)
             rb_e=radius-rnode_e

             xb_e=x1(ii)+rb_e*T_ni_e
             yb_e=x2(jj)+rb_e*T_nj_e
             zb_e=x3(kk)+rb_e*T_nk_e

             x_prob_e=xb_e+dx1_i(0)*T_ni_e
             y_prob_e=yb_e+dx1_i(0)*T_nj_e
             z_prob_e=zb_e+dx1_i(0)*T_nk_e

             nx_prob_e=floor(x_prob_e/dx1_i(0)+0.5d0)
             ny_prob_e=floor(y_prob_e/dx2_j(0)+0.5d0)
             nz_prob_e=floor(z_prob_e/dx3_k(0)+0.5d0)

             xd=(x_prob_e-x1(nx_prob_e))/(dx1(1))
             yd=(y_prob_e-x2(ny_prob_e))/(dx2(1))
             zd=(z_prob_e-x3(nz_prob_e))/(dx3(1))


             c00=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e,2)

             c01=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e+1,2)

             c10=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1,nz_prob_e,2)
 
             c11=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1
     $       ,nz_prob_e+1,2)

             c0=c00*(1-yd)+c10*yd
             c1=c01*(1-yd)+c11*yd
 
             c(ii,jj,kk,2)=c0*(1-zd)+c1*zd


          end if

           
         end do
         end do
         end do

         end do


c     x-component u2

      xd=(x_prob_V(o,mm)-x1_i(nx_prob_VU(o,mm)))/(dx1(1))
      yd=(y_prob_V(o,mm)-x2(ny_prob_VU(o,mm)))/(dx2_j(1))
      zd=(z_prob_V(o,mm)-x3(nz_prob_VU(o,mm)))/(dx3_k(1))


      c00=(1d0-xd)*
     $  u1(nx_prob_VU(o,mm),ny_prob_VU(o,mm),nz_prob_VU(o,mm))+
     $  xd*u1(nx_prob_VU(o,mm)+1,ny_prob_VU(o,mm),nz_prob_VU(o,mm))

      c01=(1d0-xd)*
     $  u1(nx_prob_VU(o,mm),ny_prob_VU(o,mm),nz_prob_VU(o,mm)+1)+
     $  xd*u1(nx_prob_VU(o,mm)+1,ny_prob_VU(o,mm),nz_prob_VU(o,mm)+1)

      c10=(1d0-xd)*
     $  u1(nx_prob_VU(o,mm),ny_prob_VU(o,mm)+1,nz_prob_VU(o,mm))+
     $  xd*u1(nx_prob_VU(o,mm)+1,ny_prob_VU(o,mm)+1,nz_prob_VU(o,mm))

      c11=(1d0-xd)*
     $  u1(nx_prob_VU(o,mm),ny_prob_VU(o,mm)+1,nz_prob_VU(o,mm)+1)+
     $  xd*u1(nx_prob_VU(o,mm)+1,ny_prob_VU(o,mm)+1,nz_prob_VU(o,mm)+1)

      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      u1j_ref=c0*(1-zd)+c1*zd


c     y-component u2


      xd=(x_prob_V(o,mm)-x1(nx_prob_VV(o,mm)))/(dx1_i(1))
      yd=(y_prob_V(o,mm)-x2_j(ny_prob_VV(o,mm)))/(dx2(1))
      zd=(z_prob_V(o,mm)-x3(nz_prob_VV(o,mm)))/(dx3_k(1))

      c00=(1d0-xd)*
     $  u2(nx_prob_VV(o,mm),ny_prob_VV(o,mm),nz_prob_VV(o,mm))+
     $  xd*u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm),nz_prob_VV(o,mm))

      c01=(1d0-xd)*
     $  u2(nx_prob_VV(o,mm),ny_prob_VV(o,mm),nz_prob_VV(o,mm)+1)+
     $  xd*u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm),nz_prob_VV(o,mm)+1)

      c10=(1d0-xd)*
     $  u2(nx_prob_VV(o,mm),ny_prob_VV(o,mm)+1,nz_prob_VV(o,mm))+
     $  xd*u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm)+1,nz_prob_VV(o,mm))

      c11=(1d0-xd)*
     $  u2(nx_prob_VV(o,mm),ny_prob_VV(o,mm)+1,nz_prob_VV(o,mm)+1)+
     $  xd*u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm)+1,nz_prob_VV(o,mm)+1)

      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      u2j_ref=c0*(1-zd)+c1*zd

      if (i.eq.4.and.k.eq.8) then
      if (j.eq.11.or.j.eq.13) then

c         write(*,*)'here',j,
c     $  u2(nx_prob_VV(o,mm),ny_prob_VV(o,mm),nz_prob_VV(o,mm)+1),
c     $  u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm),nz_prob_VV(o,mm)+1)

      end if
      end if



      if (0d0.eq.1d0) then
      if (i.eq.4.and.k.eq.8) then
      if (j.eq.11) then
         
         write(*,*)'here',j,o,nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm),
     $   nz_prob_VV(o,mm)+1,
     $   u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm),nz_prob_VV(o,mm)+1)



      end if

      if (j.eq.13) then
         
         write(*,*)'here',j,o,
     $   nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm)+1,nz_prob_VV(o,mm)+1,
     $   u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm)+1,nz_prob_VV(o,mm)+1)


      end if

      end if
      end if






c      z-component u2


      xd=(x_prob_V(o,mm)-x1(nx_prob_VW(o,mm)))/(dx1_i(1))
      yd=(y_prob_V(o,mm)-x2(ny_prob_VW(o,mm)))/(dx2_j(1))
      zd=(z_prob_V(o,mm)-x3_k(nz_prob_VW(o,mm)))/(dx3(1))

      c00=(1d0-xd)*
     $  u3(nx_prob_VW(o,mm),ny_prob_VW(o,mm),nz_prob_VW(o,mm))+
     $  xd*u3(nx_prob_VW(o,mm)+1,ny_prob_VW(o,mm),nz_prob_VW(o,mm))

      c01=(1d0-xd)*
     $  u3(nx_prob_VW(o,mm),ny_prob_VW(o,mm),nz_prob_VW(o,mm)+1)+
     $  xd*u3(nx_prob_VW(o,mm)+1,ny_prob_VW(o,mm),nz_prob_VW(o,mm)+1)

      c10=(1d0-xd)*
     $  u3(nx_prob_VW(o,mm),ny_prob_VW(o,mm)+1,nz_prob_VW(o,mm))+
     $  xd*u3(nx_prob_VW(o,mm)+1,ny_prob_VW(o,mm)+1,nz_prob_VW(o,mm))

      c11=(1d0-xd)*
     $  u3(nx_prob_VW(o,mm),ny_prob_VW(o,mm)+1,nz_prob_VW(o,mm)+1)+
     $  xd*u3(nx_prob_VW(o,mm)+1,ny_prob_VW(o,mm)+1,nz_prob_VW(o,mm)+1)

      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      u3j_ref=c0*(1-zd)+c1*zd

      do nn=1,2

         xd=(x_prob_mt(o,nn)-x1(nx_prob_mt(o,nn)))/(dx1_i(1))
         yd=(y_prob_mt(o,nn)-x2(ny_prob_mt(o,nn)))/(dx2_j(1))
         zd=(z_prob_mt(o,nn)-x3(nz_prob_mt(o,nn)))/(dx3_k(1))

         c00=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn),nz_prob_mt(o,nn),2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn),nz_prob_mt(o,nn),2)

         c01=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn),nz_prob_mt(o,nn)+1,2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn),nz_prob_mt(o,nn)+1,2)

         c10=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn)+1,nz_prob_mt(o,nn),2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn)+1,nz_prob_mt(o,nn),2)

         c11=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn)+1,nz_prob_mt(o,nn)+1,2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn)+1
     $   ,nz_prob_mt(o,nn)+1,2)

         c0=c00*(1-yd)+c10*yd
         c1=c01*(1-yd)+c11*yd

         c_tt(nn)=c0*(1-zd)+c1*zd

         xd=(x_prob_ms(o,nn)-x1(nx_prob_ms(o,nn)))/(dx1_i(1))
         yd=(y_prob_ms(o,nn)-x2(ny_prob_ms(o,nn)))/(dx2_j(1))
         zd=(z_prob_ms(o,nn)-x3(nz_prob_ms(o,nn)))/(dx3_k(1))

         c00=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn),nz_prob_ms(o,nn),2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn),nz_prob_ms(o,nn),2)

         c01=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn),nz_prob_ms(o,nn)+1,2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn),nz_prob_ms(o,nn)+1,2)

         c10=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn),2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn),2)

         c11=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn)+1,2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn)+1
     $   ,nz_prob_ms(o,nn)+1,2)

         c0=c00*(1-yd)+c10*yd
         c1=c01*(1-yd)+c11*yd

         c_ss(nn)=c0*(1-zd)+c1*zd

      if (i.eq.9.and.k.eq.8) then
      if (j.eq.9.or.j.eq.15) then
      if (nn.eq.2) then

c         write(*,*)j,o,
c     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn),nz_prob_ms(o,nn)+1,2),
c     $   nx_prob_ms(o,nn),ny_prob_ms(o,nn),nz_prob_ms(o,nn)+1

c         write(*,*)j,o,
c     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn)+1,2),
c     $   nx_prob_ms(o,nn),ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn)+1








c         write(*,*)j,o,c01,yd,c11,nx_prob_ms(o,nn),
c     $ny_prob_ms(o,nn),nz_prob_ms(o,nn)
c         write(*,*)xd,yd,zd
c         write(*,*)x_prob_ms(o,nn),y_prob_ms(o,nn),
c     $z_prob_ms(o,nn)
      end if
      end if
      end if

      end do


       constant_t=(c_tt(1)-c_tt(2))/(2d0*dx1_i(0))
     $           *dsigma/(ro*nu)
       constant_s=(c_ss(1)-c_ss(2))/(2d0*dx1_i(0))
     $           *dsigma/(ro*nu)



      if (i.eq.9.and.k.eq.8) then
      if (j.eq.9.or.j.eq.15) then
c         write(*,*)j,o,c_ss(1),c_ss(2)
      end if
      end if




        unj_ref=u1j_ref*T_ni_V(o)+
     $          u2j_ref*T_nj_V(o)+
     $          u3j_ref*T_nk_V(o)
        utj_ref=u1j_ref*T_ti_V(o)+
     $          u2j_ref*T_tj_V(o)+
     $          u3j_ref*T_tk_V(o)
        usj_ref=u1j_ref*T_si_V(o)+
     $          u2j_ref*T_sj_V(o)+
     $          u3j_ref*T_sk_V(o)

c       Dirichlet BC for normal direction
        beta_V=rb_V(o)/dx2(j)
c       Velocity at the bubble surface in cartisian coordinate
        u_star=dxcdt+drdt*T_ni_V(o)
        v_star=      drdt*T_nj_V(o)
        w_star=      drdt*T_nk_V(o)




        unjs=  u_star*T_ni_V(o)+
     $         v_star*T_nj_V(o)+
     $         w_star*T_nk_V(o)
        utjs=  u_star*T_ti_V(o)+
     $         v_star*T_tj_V(o)+
     $         w_star*T_tk_V(o)
        usjs=  u_star*T_si_V(o)+
     $         v_star*T_sj_V(o)+
     $         w_star*T_sk_V(o)


        R_IB_V=radius-rb_V(o)
        R_P_V=radius+dx1(1)



        unj=(1+beta_V)*unjs-(beta_V)*unj_ref
        utj=(1+beta_V)*utjs-(beta_V)*utj_ref
        usj=(1+beta_V)*usjs-(beta_V)*usj_ref

      if (i.eq.4.and.k.eq.8) then
      if (j.eq.11.or.j.eq.13) then
c         write(*,*)j,o,u1j_ref,u2j_ref,u3j_ref
      end if
      end if

      

c       Solve u2



c          u2(i,j,k)=unj*T_nj_V(o)+
c     $              utj*T_tj_V(o)+
c     $              usj*T_sj_V(o)


      u2_temp(i,j,k)=sin(phiS_V(o))*sin(thetaS_V(o))*unj+
     $                           cos(phiS_V(o))*utj+
     $          sin(phiS_V(o))*cos(thetaS_V(o))*usj


      if (i.eq.4.and.k.eq.8) then
      if (j.eq.11.or.j.eq.13) then
c         write(*,*)j,o,u2(i,j,k),unj,utj,usj
      end if
      end if

      if (radius.ge.R_max.or.m1.eq.m1_max) then
      un_surf_v(o)=unjs
      ut_surf_v(o)=utjs
      us_surf_v(o)=usjs
      dc_t_v(o)=(c_tt(1)-c_tt(2))/(2d0*dx1_i(0))
      dc_s_v(o)=(c_ss(1)-c_ss(2))/(2d0*dx1_i(0))
      end if


      


      end do

      u2=u2_temp

      return
      end

      subroutine finding_w_IB
      include 'var.cmn' 
      integer o
      typ_w_IB(0:im+1,0:jm+1,0:km+1)=0

      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)

          
         if (typ(i,j,k+1).eq.1) then

            typ_w_IB(i,j,k)=1

         else

          if (typ_IB(i,j,k+1).eq.1) then

             typ_w_IB(i,j,k)=1
           
          end if

         end if


         if (typ(i,j,k-1).eq.1) then

            typ_w_IB(i,j,k-1)=1

         else

          if (typ_IB(i,j,k-1).eq.1) then

             typ_w_IB(i,j,k-1)=1
           
          end if


         end if


      end do

      o_w=0
      do k=0,km+1
      do j=0,jm+1
      do i=0,im+1

         if (typ_w_IB(i,j,k).eq.1) then

          o_w=o_w+1

          x_IB_W(o_w)=i
          y_IB_W(o_w)=j
          z_IB_W(o_w)=k

         end if

      end do
      end do
      end do




      
      return
      end



      subroutine probe_IB_W
      include 'var.cmn' 
      integer o,mm
      real*8  rnode_W

c        Neumann BC in Spherical  coordinates

c        phiS = angle phi in spherical coordinates.
c        not to be confused with phi the potential
c        if x=0 then atan(y,x)=pi/2. This is not true however,
c        if y is negative, then it should be 3pi/2.

      mm=1
      
      do o=1,o_w
         i=x_IB_W(o)
         j=y_IB_W(o)
         k=z_IB_W(o)

         phiS_W(o)=datan2((x2(j)-yc),(x1(i)-xc))

         thetaS_W(o)=acos((x3_k(k)-zc)/sqrt((x1(i)-xc)**2+
     $    (x2(j)-yc)**2+(x3_k(k)-zc)**2))


         T_ni_W(o)=cos(phiS_W(o))*sin(thetaS_W(o))
         T_nj_W(o)=sin(phiS_W(o))*sin(thetaS_W(o))
         T_nk_W(o)=cos(thetaS_W(o))

         T_ti_W(o)=-sin(phiS_W(o))
         T_tj_W(o)=cos(phiS_W(o))
         T_tk_W(o)=0d0

         T_si_W(o)=cos(phiS_W(o))*cos(thetaS_W(o))
         T_sj_W(o)=sin(phiS_W(o))*cos(thetaS_W(o))
         T_sk_W(o)=-sin(thetaS_W(o))

c        initialise x-direction
c        Normal distance between interface and bubble centre
         rnode_W=sqrt((x1(i)-xc)**2+(x2(j)-yc)**2+(x3_k(k)-zc)**2)


c        normal distance of interface to bubble surface 
c        (can be negative)
         rb_W(o)=radius-rnode_W

         xb_W(o)=x1(i)  +rb_W(o)*T_ni_W(o)
         yb_W(o)=x2(j)  +rb_W(o)*T_nj_W(o)
         zb_W(o)=x3_k(k)+rb_W(o)*T_nk_W(o)


         

         x_prob_W(o,mm)=xb_W(o)+dx1_i(0)*T_ni_W(o)
         y_prob_W(o,mm)=yb_W(o)+dx1_i(0)*T_nj_W(o)
         z_prob_W(o,mm)=zb_W(o)+dx1_i(0)*T_nk_W(o)

         x_prob_mt(o,1)=xb_W(o)+dx1_i(0)*(T_ti_W(o)+dd*T_ni_w(o))
         y_prob_mt(o,1)=yb_W(o)+dx1_i(0)*(T_tj_W(o)+dd*T_nj_w(o))
         z_prob_mt(o,1)=zb_W(o)+dx1_i(0)*(T_tk_W(o)+dd*T_nk_w(o))

         x_prob_mt(o,2)=xb_W(o)-dx1_i(0)*(T_ti_W(o)-dd*T_ni_w(o))
         y_prob_mt(o,2)=yb_W(o)-dx1_i(0)*(T_tj_W(o)-dd*T_nj_w(o))
         z_prob_mt(o,2)=zb_W(o)-dx1_i(0)*(T_tk_W(o)-dd*T_nk_w(o))

         x_prob_ms(o,1)=xb_W(o)+dx1_i(0)*(T_si_W(o)+dd*T_ni_w(o))
         y_prob_ms(o,1)=yb_W(o)+dx1_i(0)*(T_sj_W(o)+dd*T_nj_w(o))
         z_prob_ms(o,1)=zb_W(o)+dx1_i(0)*(T_sk_W(o)+dd*T_nk_w(o))

         x_prob_ms(o,2)=xb_W(o)-dx1_i(0)*(T_si_W(o)-dd*T_ni_w(o))
         y_prob_ms(o,2)=yb_W(o)-dx1_i(0)*(T_sj_W(o)-dd*T_nj_w(o))
         z_prob_ms(o,2)=zb_W(o)-dx1_i(0)*(T_sk_W(o)-dd*T_nk_w(o))



c        prob number of cells starts at an interface so only 
c        +0.5d if that direction is on the level of a node
         nx_prob_WU(o,mm)=floor(x_prob_W(o,mm)/dx1_i(0))
         ny_prob_WU(o,mm)=floor(y_prob_W(o,mm)/dx2_j(0)+0.5d0)
         nz_prob_WU(o,mm)=floor(z_prob_W(o,mm)/dx3_k(0)+0.5d0)

         nx_prob_WV(o,mm)=floor(x_prob_W(o,mm)/dx1_i(0)+0.5d0)
         ny_prob_WV(o,mm)=floor(y_prob_W(o,mm)/dx2_j(0))
         nz_prob_WV(o,mm)=floor(z_prob_W(o,mm)/dx3_k(0)+0.5d0)

         nx_prob_WW(o,mm)=floor(x_prob_W(o,mm)/dx1_i(0)+0.5d0)
         ny_prob_WW(o,mm)=floor(y_prob_W(o,mm)/dx2_j(0)+0.5d0)
         nz_prob_WW(o,mm)=floor(z_prob_W(o,mm)/dx3_k(0))

         nx_prob_mt(o,1)=floor(x_prob_mt(o,1)/dx1_i(0)+0.5d0)
         ny_prob_mt(o,1)=floor(y_prob_mt(o,1)/dx2_j(0)+0.5d0)
         nz_prob_mt(o,1)=floor(z_prob_mt(o,1)/dx3_k(0)+0.5d0)

         nx_prob_mt(o,2)=floor(x_prob_mt(o,2)/dx1_i(0)+0.5d0)
         ny_prob_mt(o,2)=floor(y_prob_mt(o,2)/dx2_j(0)+0.5d0)
         nz_prob_mt(o,2)=floor(z_prob_mt(o,2)/dx3_k(0)+0.5d0)

         nx_prob_ms(o,1)=floor(x_prob_ms(o,1)/dx1_i(0)+0.5d0)
         ny_prob_ms(o,1)=floor(y_prob_ms(o,1)/dx2_j(0)+0.5d0)
         nz_prob_ms(o,1)=floor(z_prob_ms(o,1)/dx3_k(0)+0.5d0)

         nx_prob_ms(o,2)=floor(x_prob_ms(o,2)/dx1_i(0)+0.5d0)
         ny_prob_ms(o,2)=floor(y_prob_ms(o,2)/dx2_j(0)+0.5d0)
         nz_prob_ms(o,2)=floor(z_prob_ms(o,2)/dx3_k(0)+0.5d0)




      end do


      return
      end


      subroutine bulb_w
      include 'var.cmn' 
      integer o,mm,nn,ii,jj,kk
      integer nx_prob_e,ny_prob_e,nz_prob_e
c     solving u3 IB interface
      real*8  u1k_ref,u2k_ref,u3k_ref
      real*8  unk_ref,utk_ref,usk_ref
      real*8  xd,yd,zd
      real*8  c00,c10,c01,c11,c0,c1
      real*8  u_star,v_star,w_star
      real*8  unks,utks,usks,u3_temp(0:im+1,0:jm+1,0:km+1)
      real*8  unk,utk,usk,beta_W

      real*8  R_IB_W,R_P_W

      real*8 c_tt(1:2),c_ss(1:2)
      real*8 constant_t,constant_s

      real*8  phiS_e,thetaS_e
      real*8  T_ni_e,T_nj_e,T_nk_e
      real*8  rnode_e,rb_e
      real*8  xb_e,yb_e,zb_e
      real*8  x_prob_e,y_prob_e,z_prob_e


      u3_temp=u3
      mm=1
      do o=1,o_w
          i=x_IB_W(o)
          j=y_IB_W(o)
          k=z_IB_W(o)

          do nn=1,2

          do kk=nz_prob_mt(o,nn),nz_prob_mt(o,nn)+1
          do jj=ny_prob_mt(o,nn),ny_prob_mt(o,nn)+1
          do ii=nx_prob_mt(o,nn),nx_prob_mt(o,nn)+1

          if(typ(ii,jj,kk)+typ_IB(ii,jj,kk).eq.0.and.ii.ne.0) then

             phiS_e=datan2((x2(jj)-yc),(x1(ii)-xc))
  
             thetaS_e=acos((x3(kk)-zc)/sqrt((x1(ii)-xc)**2+
     $                (x2(jj)-yc)**2+(x3(kk)-zc)**2))

             T_ni_e=cos(phiS_e)*sin(thetaS_e)
             T_nj_e=sin(phiS_e)*sin(thetaS_e)
             T_nk_e=cos(thetaS_e)

c            initialise x-direction
c            Normal distance between interface and bubble centre
         rnode_e=sqrt((x1(ii)-xc)**2+(x2(jj)-yc)**2
     $       +(x3(kk)-zc)**2)


c            normal distance of interface to bubble surface 
c            (can be negative)
             rb_e=radius-rnode_e

             xb_e=x1(ii)+rb_e*T_ni_e
             yb_e=x2(jj)+rb_e*T_nj_e
             zb_e=x3(kk)+rb_e*T_nk_e

             x_prob_e=xb_e+dx1_i(0)*T_ni_e
             y_prob_e=yb_e+dx1_i(0)*T_nj_e
             z_prob_e=zb_e+dx1_i(0)*T_nk_e

             nx_prob_e=floor(x_prob_e/dx1_i(0)+0.5d0)
             ny_prob_e=floor(y_prob_e/dx2_j(0)+0.5d0)
             nz_prob_e=floor(z_prob_e/dx3_k(0)+0.5d0)

             xd=(x_prob_e-x1(nx_prob_e))/(dx1(1))
             yd=(y_prob_e-x2(ny_prob_e))/(dx2(1))
             zd=(z_prob_e-x3(nz_prob_e))/(dx3(1))


             c00=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e,2)

             c01=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e+1,2)

             c10=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1,nz_prob_e,2)
 
             c11=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1
     $       ,nz_prob_e+1,2)

             c0=c00*(1-yd)+c10*yd
             c1=c01*(1-yd)+c11*yd
 
             c(ii,jj,kk,2)=c0*(1-zd)+c1*zd

          end if

           
          end do
          end do
          end do

          do kk=nz_prob_ms(o,nn),nz_prob_ms(o,nn)+1
          do jj=ny_prob_ms(o,nn),ny_prob_ms(o,nn)+1
          do ii=nx_prob_ms(o,nn),nx_prob_ms(o,nn)+1

          if(typ(ii,jj,kk)+typ_IB(ii,jj,kk).eq.0.and.ii.ne.0) then


             phiS_e=datan2((x2(jj)-yc),(x1(ii)-xc))
  
             thetaS_e=acos((x3(kk)-zc)/sqrt((x1(ii)-xc)**2+
     $                (x2(jj)-yc)**2+(x3(kk)-zc)**2))

             T_ni_e=cos(phiS_e)*sin(thetaS_e)
             T_nj_e=sin(phiS_e)*sin(thetaS_e)
             T_nk_e=cos(thetaS_e)

c            initialise x-direction
c            Normal distance between interface and bubble centre
         rnode_e=sqrt((x1(ii)-xc)**2+(x2(jj)-yc)**2
     $       +(x3(kk)-zc)**2)


c            normal distance of interface to bubble surface 
c            (can be negative)
             rb_e=radius-rnode_e

             xb_e=x1(ii)+rb_e*T_ni_e
             yb_e=x2(jj)+rb_e*T_nj_e
             zb_e=x3(kk)+rb_e*T_nk_e

             x_prob_e=xb_e+dx1_i(0)*T_ni_e
             y_prob_e=yb_e+dx1_i(0)*T_nj_e
             z_prob_e=zb_e+dx1_i(0)*T_nk_e

             nx_prob_e=floor(x_prob_e/dx1_i(0)+0.5d0)
             ny_prob_e=floor(y_prob_e/dx2_j(0)+0.5d0)
             nz_prob_e=floor(z_prob_e/dx3_k(0)+0.5d0)

             xd=(x_prob_e-x1(nx_prob_e))/(dx1(1))
             yd=(y_prob_e-x2(ny_prob_e))/(dx2(1))
             zd=(z_prob_e-x3(nz_prob_e))/(dx3(1))


             c00=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e,2)

             c01=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e,nz_prob_e+1,2)

             c10=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1,nz_prob_e,2)
 
             c11=(1d0-xd)*
     $       c(nx_prob_e,ny_prob_e+1,nz_prob_e+1,2)+
     $       xd*c(nx_prob_e+1,ny_prob_e+1
     $       ,nz_prob_e+1,2)

             c0=c00*(1-yd)+c10*yd
             c1=c01*(1-yd)+c11*yd
 
             c(ii,jj,kk,2)=c0*(1-zd)+c1*zd



          end if

           
         end do
         end do
         end do

         end do



c     x-component u3

      xd=(x_prob_W(o,mm)-x1_i(nx_prob_WU(o,mm)))/(dx1(1))
      yd=(y_prob_W(o,mm)-x2(ny_prob_WU(o,mm)))/(dx2_j(1))
      zd=(z_prob_W(o,mm)-x3(nz_prob_WU(o,mm)))/(dx3_k(1))


      c00=(1d0-xd)*
     $  u1(nx_prob_WU(o,mm),ny_prob_WU(o,mm),nz_prob_WU(o,mm))+
     $  xd*u1(nx_prob_WU(o,mm)+1,ny_prob_WU(o,mm),nz_prob_WU(o,mm))

      c01=(1d0-xd)*
     $  u1(nx_prob_WU(o,mm),ny_prob_WU(o,mm),nz_prob_WU(o,mm)+1)+
     $  xd*u1(nx_prob_WU(o,mm)+1,ny_prob_WU(o,mm),nz_prob_WU(o,mm)+1)

      c10=(1d0-xd)*
     $  u1(nx_prob_WU(o,mm),ny_prob_WU(o,mm)+1,nz_prob_WU(o,mm))+
     $  xd*u1(nx_prob_WU(o,mm)+1,ny_prob_WU(o,mm)+1,nz_prob_WU(o,mm))

      c11=(1d0-xd)*
     $  u1(nx_prob_WU(o,mm),ny_prob_WU(o,mm)+1,nz_prob_WU(o,mm)+1)+
     $  xd*u1(nx_prob_WU(o,mm)+1,ny_prob_WU(o,mm)+1,nz_prob_WU(o,mm)+1)

      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      u1k_ref=c0*(1-zd)+c1*zd


c       y-component u3


      xd=(x_prob_W(o,mm)-x1(nx_prob_WV(o,mm)))/(dx1_i(1))
      yd=(y_prob_W(o,mm)-x2_j(ny_prob_WV(o,mm)))/(dx2(1))
      zd=(z_prob_W(o,mm)-x3(nz_prob_WV(o,mm)))/(dx3_k(1))

      c00=(1d0-xd)*
     $  u2(nx_prob_WV(o,mm),ny_prob_WV(o,mm),nz_prob_WV(o,mm))+
     $  xd*u2(nx_prob_WV(o,mm)+1,ny_prob_WV(o,mm),nz_prob_WV(o,mm))

      c01=(1d0-xd)*
     $  u2(nx_prob_WV(o,mm),ny_prob_WV(o,mm),nz_prob_WV(o,mm)+1)+
     $  xd*u2(nx_prob_WV(o,mm)+1,ny_prob_WV(o,mm),nz_prob_WV(o,mm)+1)

      c10=(1d0-xd)*
     $  u2(nx_prob_WV(o,mm),ny_prob_WV(o,mm)+1,nz_prob_WV(o,mm))+
     $  xd*u2(nx_prob_WV(o,mm)+1,ny_prob_WV(o,mm)+1,nz_prob_WV(o,mm))

      c11=(1d0-xd)*
     $  u2(nx_prob_WV(o,mm),ny_prob_WV(o,mm)+1,nz_prob_WV(o,mm)+1)+
     $  xd*u2(nx_prob_WV(o,mm)+1,ny_prob_WV(o,mm)+1,nz_prob_WV(o,mm)+1)

      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      u2k_ref=c0*(1-zd)+c1*zd

c     z-component u3

      xd=(x_prob_W(o,mm)-x1(nx_prob_WW(o,mm)))/(dx1_i(1))
      yd=(y_prob_W(o,mm)-x2(ny_prob_WW(o,mm)))/(dx2_j(1))
      zd=(z_prob_W(o,mm)-x3_k(nz_prob_WW(o,mm)))/(dx3(1))

      c00=(1d0-xd)*
     $  u3(nx_prob_WW(o,mm),ny_prob_WW(o,mm),nz_prob_WW(o,mm))+
     $  xd*u3(nx_prob_WW(o,mm)+1,ny_prob_WW(o,mm),nz_prob_WW(o,mm))

      c01=(1d0-xd)*
     $  u3(nx_prob_WW(o,mm),ny_prob_WW(o,mm),nz_prob_WW(o,mm)+1)+
     $  xd*u3(nx_prob_WW(o,mm)+1,ny_prob_WW(o,mm),nz_prob_WW(o,mm)+1)

      c10=(1d0-xd)*
     $  u3(nx_prob_WW(o,mm),ny_prob_WW(o,mm)+1,nz_prob_WW(o,mm))+
     $  xd*u3(nx_prob_WW(o,mm)+1,ny_prob_WW(o,mm)+1,nz_prob_WW(o,mm))

      c11=(1d0-xd)*
     $  u3(nx_prob_WW(o,mm),ny_prob_WW(o,mm)+1,nz_prob_WW(o,mm)+1)+
     $  xd*u3(nx_prob_WW(o,mm)+1,ny_prob_WW(o,mm)+1,nz_prob_WW(o,mm)+1)

      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      u3k_ref=c0*(1-zd)+c1*zd

      do nn=1,2

         xd=(x_prob_mt(o,nn)-x1(nx_prob_mt(o,nn)))/(dx1_i(1))
         yd=(y_prob_mt(o,nn)-x2(ny_prob_mt(o,nn)))/(dx2_j(1))
         zd=(z_prob_mt(o,nn)-x3(nz_prob_mt(o,nn)))/(dx3_k(1))

         c00=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn),nz_prob_mt(o,nn),2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn),nz_prob_mt(o,nn),2)

         c01=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn),nz_prob_mt(o,nn)+1,2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn),nz_prob_mt(o,nn)+1,2)

         c10=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn)+1,nz_prob_mt(o,nn),2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn)+1,nz_prob_mt(o,nn),2)

         c11=(1d0-xd)*
     $   c(nx_prob_mt(o,nn),ny_prob_mt(o,nn)+1,nz_prob_mt(o,nn)+1,2)+
     $   xd*c(nx_prob_mt(o,nn)+1,ny_prob_mt(o,nn)+1
     $   ,nz_prob_mt(o,nn)+1,2)

         c0=c00*(1-yd)+c10*yd
         c1=c01*(1-yd)+c11*yd

         c_tt(nn)=c0*(1-zd)+c1*zd

         xd=(x_prob_ms(o,nn)-x1(nx_prob_ms(o,nn)))/(dx1_i(1))
         yd=(y_prob_ms(o,nn)-x2(ny_prob_ms(o,nn)))/(dx2_j(1))
         zd=(z_prob_ms(o,nn)-x3(nz_prob_ms(o,nn)))/(dx3_k(1))

         c00=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn),nz_prob_ms(o,nn),2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn),nz_prob_ms(o,nn),2)

         c01=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn),nz_prob_ms(o,nn)+1,2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn),nz_prob_ms(o,nn)+1,2)

         c10=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn),2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn),2)

         c11=(1d0-xd)*
     $   c(nx_prob_ms(o,nn),ny_prob_ms(o,nn)+1,nz_prob_ms(o,nn)+1,2)+
     $   xd*c(nx_prob_ms(o,nn)+1,ny_prob_ms(o,nn)+1
     $   ,nz_prob_ms(o,nn)+1,2)

         c0=c00*(1-yd)+c10*yd
         c1=c01*(1-yd)+c11*yd

         c_ss(nn)=c0*(1-zd)+c1*zd

      end do

      constant_t=(c_tt(1)-c_tt(2))/(2d0*dx1_i(0))
     $           *dsigma/(ro*nu)
      constant_s=(c_ss(1)-c_ss(2))/(2d0*dx1_i(0))
     $           *dsigma/(ro*nu)


        unk_ref=u1k_ref*T_ni_W(o)+
     $          u2k_ref*T_nj_W(o)+
     $          u3k_ref*T_nk_W(o)
        utk_ref=u1k_ref*T_ti_W(o)+
     $          u2k_ref*T_tj_W(o)+
     $          u3k_ref*T_tk_W(o)
        usk_ref=u1k_ref*T_si_W(o)+
     $          u2k_ref*T_sj_W(o)+
     $          u3k_ref*T_sk_W(o)


c     Dirichlet BC for normal direction
      beta_W=rb_W(o)/dx3(k)
      u_star=dxcdt+drdt*T_ni_W(o)
      v_star=      drdt*T_nj_W(o)
      w_star=      drdt*T_nk_W(o)

        unks=  u_star*T_ni_W(o)+
     $         v_star*T_nj_W(o)+
     $         w_star*T_nk_W(o)
        utks=  u_star*T_ti_W(o)+
     $         v_star*T_tj_W(o)+
     $         w_star*T_tk_W(o)
        usks=  u_star*T_si_W(o)+
     $         v_star*T_sj_W(o)+
     $         w_star*T_sk_W(o)


      R_IB_W=radius-rb_W(o)
      R_P_W=radius+dx1(1)



      unk=(1+beta_W)*unks-(beta_W)*unk_ref
      utk=(1+beta_W)*utks-(beta_W)*utk_ref
      usk=(1+beta_W)*usks-(beta_W)*usk_ref



c     Solve u3


c       u3(i,j,k)=unk*T_nk_W(o)+
c     $   utk*T_tk_W(o)+usk*T_sk_W(o)

      u3_temp(i,j,k)=cos(thetaS_W(o))*unk-
     $          sin(thetaS_W(o))*usk


      if (radius.ge.R_max.or.m1.eq.m1_max) then
      un_surf_w(o)=unks
      ut_surf_w(o)=utks
      us_surf_w(o)=usks
      dc_t_w(o)=(c_tt(1)-c_tt(2))/(2d0*dx1_i(0))
      dc_s_w(o)=(c_ss(1)-c_ss(2))/(2d0*dx1_i(0))
      end if




      end do

      u3=u3_temp

      return
      end


      subroutine finding_vv_IB
      include 'var.cmn' 
      integer o

      typ_v_IB(0:im+1,0:jm+1,0:km+1)=0

c     if cell centre is inside of the bubble typ(i,j,k)=0
c     if cell centre is outside of the bubble typ(i,j,k)=1
c     if typ_IB(i,j,k) is 1 means it is IB cell

      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)

          
         if (typ(i,j+1,k).eq.1) then

            typ_v_IB(i,j,k)=1

          else

            if(typ_IB(i,j+1,k).eq.1) then

             typ_v_IB(i,j,k)=1

            end if

         end if


         if (typ(i,j-1,k).eq.1) then

            typ_v_IB(i,j-1,k)=1

          else

           if(typ_IB(i,j-1,k).eq.1) then

             typ_v_IB(i,j-1,k)=1

           end if

         end if


      end do



      o_v=0
      do k=0,km+1
      do j=0,jm+1
      do i=0,im+1

         if (typ_v_IB(i,j,k).eq.1) then

          o_v=o_v+1

          x_IB_V(o_v)=i
          y_IB_V(o_v)=j
          z_IB_V(o_v)=k

         end if

      end do
      end do
      end do

     

 55   format(7i9)
      open(unit=55,file='vIB.dat')

      do o=1,o_v
       write(55,55)o,x_IB_V(o),
     $   y_IB_V(o),z_IB_V(o)
      end do
      close(55)


      return
      end



      subroutine bulb_force_x
      include 'var.cmn' 
      integer o,counter_write,typ_force_x(0:im+1,0:jm+1,0:km+1)
      real*8 term1,term2,term3,total_up,total_bellow

      total=0d0
      total_bellow=0d0
      total_up=0d0
      typ_force_x(0:im+1,0:jm+1,0:km+1)=0
c     Label the velocity points that are outside of the bubble, and set their type to 1. Label the velocities inside of
c     the bubble are set to type zero
      do k=0,km+1
      do j=0,jm+1
      do i=0,im+1

         if ((x1_i(i)-xc)**2+(x2(j)-yc)**2
     $   +(x3(k)-zc)**2.gt.radius**2) then

            typ_force_x(i,j,k)=1

         end if


      end do
      end do
      end do

      do o=1,o_u
         i=x_IB_u(o)
         j=y_IB_u(o)
         k=z_IB_u(o)

         typ_force_x(i,j,k)=0

      end do
c     Label all IB_u-velocities to type zero.
      do o=1,o_u
         i=x_IB_u(o)
         j=y_IB_u(o)
         k=z_IB_u(o)
c        Select only the surfaces of the cube based on the velocity next to this surface. That velocity should be outside 
c        of the bubble. Do not include the ones inside the bubble!
c        The equation used is the steady state NS-equation, transformed to a surface interagral of  to calculate the force.
         if(typ_force_x(i+1,j,k).eq.1) then

           term1=-(-(0.5d0*(u1(i+1,j,k)+u1(i,j,k)))**2*dx2(j)*dx3(k))
           term2=-(-p(i+1,j,k)/ro*dx2(j)*dx3(k))
           term3=-(nu*(u1(i+1,j,k)-u1(i,j,k))/dx1(i+1)*dx2(j)*dx3(k))    
           total=total+term1+term2+term3
c        The below if-statement is not necessary for computation.
           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if


         end if

         if(typ_force_x(i-1,j,k).eq.1) then

           term1=-(0.5d0*(u1(i-1,j,k)+u1(i,j,k)))**2*dx2(j)*dx3(k)
           term2=-p(i,j,k)/ro*dx2(j)*dx3(k)
           term3=nu*(u1(i,j,k)-u1(i-1,j,k))/dx1(i)*dx2(j)*dx3(k)
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if


         end if

         if(typ_force_x(i,j+1,k).eq.1) then

           term1=-(-0.5d0*(u1(i,j,k)+u1(i,j+1,k))
     $     *0.5d0*(u2(i,j,k)+u2(i+1,j,k))*dx1(i)*dx3(k))
           term2=0d0 
           term3=-(nu*(u1(i,j+1,k)-u1(i,j,k))/dx2_j(j)*dx1(i)*dx3(k)) 
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if


         end if

         if(typ_force_x(i,j-1,k).eq.1) then

           term1=-0.5d0*(u1(i,j,k)+u1(i,j-1,k))
     $     *0.5d0*(u2(i,j-1,k)+u2(i+1,j-1,k))*dx1(i)*dx3(k)
           term2=0d0
           term3=nu*(u1(i,j,k)-u1(i,j-1,k))/dx2_j(j-1)*dx1(i)*dx3(k)    
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if


         end if

         if(typ_force_x(i,j,k+1).eq.1) then

           term1=-(-0.5d0*(u1(i,j,k)+u1(i,j,k+1))
     $    *0.5d0*(u3(i,j,k)+u3(i+1,j,k))*dx1(i)*dx2(j))
           term2=0d0    
           term3=-(nu*(u1(i,j,k+1)-u1(i,j,k))/dx3_k(k)*dx1(i)*dx2(j))   
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if


         end if

         if(typ_force_x(i,j,k-1).eq.1) then

           term1=-0.5d0*(u1(i,j,k)+u1(i,j,k-1))
     $     *0.5d0*(u3(i,j,k-1)+u3(i+1,j,k-1))*dx1(i)*dx2(j)
           term2=0d0    
           term3=nu*(u1(i,j,k)-u1(i,j,k-1))/dx3_k(k-1)*dx1(i)*dx2(j)              
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if


         end if

      end do

      if (time.lt.0.5d0*dtime) then
       open(unit=79,file='Fx_top_bottom.dat')
      else 
       open(unit=79,file='Fx_top_bottom.dat',position='append')
      end if

 30   format(e14.6,1i9,1000e36.24)

      write(6,30)time,m1,total_up,total_bellow,total_up+total_bellow

      write(79,30)time,m1,total_up,total_bellow,total_up+total_bellow



      close(79)


      return
      end

      subroutine domain_flux_force_x
      include 'var.cmn' 
      integer o       
      real*8 term1,term2,term3,term4,term5,term6,
     $       A11(1:jm,1:km),B11(1:jm,1:km),
     $       A12(1:jm,1:km),B12(1:jm,1:km),
     $       A13(1:jm,1:km),B13(1:jm,1:km),
     $       A21(1:im,1:km),B21(1:im,1:km),
     $       A22(1:im,1:km),B22(1:im,1:km),
     $       A23(1:im,1:km),B23(1:im,1:km),
     $       A31(1:im,1:jm),B31(1:im,1:jm),
     $       A32(1:im,1:jm),B32(1:im,1:jm),
     $       A33(1:im,1:jm),B33(1:im,1:jm),Error

      a11=0d0
      b11=0d0
      a12=0d0
      b12=0d0
      a13=0d0
      b13=0d0

      a21=0d0
      b21=0d0
      a22=0d0
      b22=0d0
      a23=0d0
      b23=0d0

      a31=0d0
      b31=0d0
      a32=0d0
      b32=0d0
      a33=0d0
      b33=0d0


      do k=1,km
      do j=1,jm
c        A11 = first term of the electrode A1, A12 is the second term of the electrode, ...
c        B11 = first term of the membrame B1, B12 is the second term of the membrame ...
c        These are part of the force balance of the domain

         A11(j,k)=-(-(0.5d0*(u1(0,j,k)+u1(1,j,k)))**2
     $             *dx2(j)*dx3(k))
         B11(j,k)=-(0.5d0*(u1(im,j,k)+u1(im-1,j,k)))**2
     $             *dx2(j)*dx3(k) 

         A12(j,k)=-(-p(1,j,k)*dx2(j)*dx3(k)/ro)
         B12(j,k)=-p(im,j,k)*dx2(j)*dx3(k)/ro

         A13(j,k)=-(nu*(u1(1,j,k)-u1(0,j,k))
     $            /dx1(1)*dx2(j)*dx3(k))
         B13(j,k)=nu*(u1(im,j,k)-u1(im-1,j,k))
     $            /dx1(1)*dx2(i)*dx3(k)


      end do
      end do

      term1=sum(A11)+sum(A12)+sum(A13)
      term2=sum(B11)+sum(B12)+sum(B13)

      do k=1,km
      do i=1,im-1

         A21(i,k)=-(-0.5d0*(u1(i,0,k)+u1(i,1,k))
     $   *0.5d0*(u2(i,0,k)+u2(i+1,0,k))*dx1(i)*dx3(k))
         B21(i,k)=-0.5d0*(u1(i,jm,k)+u1(i,jm+1,k))
     $   *0.5d0*(u2(i,jm,k)+u2(i+1,jm,k))*dx1(i)*dx3(k)
  
         A22(i,k)=0d0
         B22(i,k)=0d0

         A23(i,k)=-nu*(u1(i,1,k)-u1(i,0,k))/dx2_j(0) 
     $   *dx1(i)*dx3(k) 
         B23(i,k)=0d0

      end do
      end do


      term3=sum(A21)+sum(A22)+sum(A23)
      term4=sum(B21)+sum(B22)+sum(B23)

      do j=1,jm
      do i=1,im-1

         A31(i,j)=-(-0.5d0*(u1(i,j,0)+u1(i,j,1))
     $   *0.5d0*(u3(i,j,0)+u3(i+1,j,0))*dx1(i)*dx2(j))
         B31(i,j)=-0.5d0*(u1(i,j,km)+u1(i,j,km+1))
     $   *0.5d0*(u3(i,j,km)+u3(i+1,j,km))*dx1(i)*dx2(j)

         A32(i,j)=0d0
         B32(i,j)=0d0

         A33(i,j)=-nu*(u1(i,j,1)-u1(i,j,0))
     $            *dx1(i)*dx2(j)/dx3_k(0)            
         B33(i,j)=nu*(u1(i,j,km+1)-u1(i,j,km))
     $            *dx1(i)*dx2(j)/dx3_k(km)

      end do
      end do


 

      term5=sum(A31)+sum(A32)+sum(A33)
      term6=sum(B31)+sum(B32)+sum(B33)

      Error=abs(term1+term2+term3+term4+term5+term6+total)

      if (time.lt.0.5d0*dtime) then
       open(unit=31,file='Fxbalance.dat')
      else 
       open(unit=31,file='Fxbalance.dat',position='append')
      end if

 30   format(e14.6,1i9,1000e36.24)

      write(6,30)time,m1,term1,term2,
     $            term3,term4,term5,term6,total,Error

      write(31,30)time,m1,term1,term2,
     $            term3,term4,term5,term6,total,Error
      close(31)

      return
      end

      subroutine bulb_force_y
      include 'var.cmn' 
      integer o,counter_write,typ_force_y(0:im+1,0:jm+1,0:km+1)
      real*8 term1,term2,term3,total_up,total_bellow

      total=0d0
      total_up=0d0
      total_bellow=0d0
      typ_force_y(0:im+1,0:jm+1,0:km+1)=0

      do k=0,km+1
      do j=0,jm+1
      do i=0,im+1


        if ((x1(i)-xc)**2+(x2_j(j)-yc)**2
     $  +(x3(k)-zc)**2.gt.radius**2) then

            typ_force_y(i,j,k)=1

        end if


      end do
      end do
      end do


      do o=1,o_v
         i=x_IB_V(o)
         j=y_IB_V(o)
         k=z_IB_V(o)

         typ_force_y(i,j,k)=0

      end do

      do o=1,o_v
         i=x_IB_V(o)
         j=y_IB_V(o)
         k=z_IB_V(o)


         if(typ_force_y(i+1,j,k).eq.1) then

           term1=-(-0.5d0*(u2(i,j,k)+u2(i+1,j,k))*
     $     0.5d0*(u1(i,j,k)+u1(i,j+1,k))*dx2(j)*dx3(k))
           term2=0d0    
           term3=-(nu*(u2(i+1,j,k)-u2(i,j,k))/dx1_i(i)*dx2(j)*dx3(k))          
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if
            
         end if

         if(typ_force_y(i-1,j,k).eq.1) then

           term1=-0.5d0*(u2(i,j,k)+u2(i-1,j,k))*
     $     0.5d0*(u1(i-1,j,k)+u1(i-1,j+1,k))*dx2(j)*dx3(k)
           term2=0d0
           term3=nu*(u2(i,j,k)-u2(i-1,j,k))/dx1_i(i-1)*dx2(j)*dx3(k)  
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if

                
         end if

         if(typ_force_y(i,j+1,k).eq.1) then

           term1=-(-(0.5d0*(u2(i,j,k)+u2(i,j+1,k)))**2*dx1(i)*dx3(k))
           term2=-(-p(i,j+1,k)/ro*dx1(i)*dx3(k))    
           term3=-(nu*(u2(i,j+1,k)-u2(i,j,k))/dx2_j(j)*dx1(i)*dx3(k))       
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if

               
         end if

         if(typ_force_y(i,j-1,k).eq.1) then

           term1=-(0.5d0*(u2(i,j,k)+u2(i,j-1,k)))**2*dx1(i)*dx3(k)
           term2=-p(i,j,k)/ro*dx1(i)*dx3(k)
           term3=nu*(u2(i,j,k)-u2(i,j-1,k))/dx2_j(j-1)*dx1(i)*dx3(k)     
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if

                
         end if

         if(typ_force_y(i,j,k+1).eq.1) then

           term1=-(-0.5d0*(u2(i,j,k)+u2(i,j,k+1))
     $     *0.5d0*(u3(i,j,k)+u3(i,j+1,k))*dx1(i)*dx2(j))
           term2=0d0    
           term3=-(nu*(u2(i,j,k+1)-u2(i,j,k))/dx3_k(k)*dx1(i)*dx2(j))   
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if


         end if

         if(typ_force_y(i,j,k-1).eq.1) then

           term1=-0.5d0*(u2(i,j,k)+u2(i,j,k-1))*
     $     0.5d0*(u3(i,j,k-1)+u3(i,j+1,k-1))*dx1(i)*dx2(j)
           term2=0d0    
           term3=nu*(u2(i,j,k)-u2(i,j,k-1))/dx3_k(k-1)*dx1(i)*dx2(j)               
           total=total+term1+term2+term3

           if (x2_j(j).gt.yc) then
              total_up=total_up+term1+term2+term3
           else
              total_bellow=total_bellow+term1+term2+term3
           end if

                
         end if

      end do  

      if (time.lt.0.5d0*dtime) then
       open(unit=69,file='Fy_top_bottom.dat')
      else 
       open(unit=69,file='Fy_top_bottom.dat',position='append')
      end if

 30   format(e14.6,1i9,1000e36.24)

      write(6,30)time,m1,total_up,total_bellow,total_up+total_bellow
      write(69,30)time,m1,total_up,total_bellow,total_up+total_bellow

      close(69)
    
      
      return
      end



      subroutine domain_flux_force_y
      include 'var.cmn' 
      integer o       
      real*8 term1,term2,term3,term4,term5,term6,
     $       A11(1:jm,1:km),B11(1:jm,1:km),
     $       A12(1:jm,1:km),B12(1:jm,1:km),
     $       A13(1:jm,1:km),B13(1:jm,1:km),
     $       A21(1:im,1:km),B21(1:im,1:km),
     $       A22(1:im,1:km),B22(1:im,1:km),
     $       A23(1:im,1:km),B23(1:im,1:km),
     $       A31(1:im,1:jm),B31(1:im,1:jm),
     $       A32(1:im,1:jm),B32(1:im,1:jm),
     $       A33(1:im,1:jm),B33(1:im,1:jm),Error


      a11=0d0
      b11=0d0
      a12=0d0
      b12=0d0
      a13=0d0
      b13=0d0

      a21=0d0
      b21=0d0
      a22=0d0
      b22=0d0
      a23=0d0
      b23=0d0

      a31=0d0
      b31=0d0
      a32=0d0
      b32=0d0
      a33=0d0
      b33=0d0

      do k=1,km
      do j=1,jm-1

         A11(j,k)=0d0
         B11(j,k)=-0.5d0*(u2(im,j,k)+u2(im+1,j,k))
     $   *0.5d0*(u1(im,j,k)+u1(im,j+1,k))*dx2(j)*dx3(k)

         A12(j,k)=0d0
         B12(j,k)=0d0

         A13(j,k)=(-nu*(u2(1,j,k)-u2(0,j,k))*
     $            dx2(j)*dx3(k)/dx1_i(0))
         B13(j,k)=nu*(u2(im+1,j,k)-u2(im,j,k))*
     $           dx2(j)*dx3(k)/dx1_i(im)

      end do
      end do

      term1=sum(A11)+sum(A12)+sum(A13) 
      term2=sum(B11)+sum(B12)+sum(B13)

      do k=1,km
      do i=1,im

 

	A21(i,k)=-(-(0.5d0*(u2(i,0,k)+u2(i,1,k)))**2
     $           *dx1(i)*dx3(k))
        B21(i,k)=-(0.5d0*(u2(i,jm,k)+u2(i,jm-1,k)))**2
     $           *dx1(i)*dx3(k) 

        A22(i,k)=-(-(p(i,1,k))*dx1(i)*dx3(k)/ro)
        B22(i,k)=-p(i,jm,k)*dx1(i)*dx3(k)/ro

        A23(i,k)=-(nu*(u2(i,1,k)-u2(i,0,k))/dx2(1)*dx1(i)*dx3(k))
        B23(i,k)=0d0

      end do
      end do

      term3=sum(A21)+sum(A22)+sum(A23)
      term4=sum(B21)+sum(B22)+sum(B23)


      do j=1,jm
      do i=1,im

         A31(i,j)=-(-0.5d0*(u2(i,j,0)+u2(i,j,1))
     $            *0.5d0*(u3(i,j,0)+u3(i,j+1,0))*dx1(i)*dx2(j))
         B31(i,j)=-0.5d0*(u2(i,j,km)+u2(i,j,km+1))
     $            *0.5d0*(u3(i,j,km)+u3(i,j+1,km))*dx1(i)*dx2(j)

         A32(i,j)=0d0
         B32(i,j)=0d0

         A33(i,j)=-(nu*(u2(i,j,1)-u2(i,j,0))
     $            *dx1(i)*dx2(j)/dx3_k(0))
         B33(i,j)=nu*(u2(i,j,km+1)-u2(i,j,km))
     $            *dx1(i)*dx2(j)/dx3_k(km)

      end do
      end do

      term5=sum(A31)+sum(A32)+sum(A33)
      term6=sum(B31)+sum(B32)+sum(B33)
      
      Error=abs(term1+term2+term3+term4+term5+term6+total)

      if (time.lt.0.5d0*dtime) then
       open(unit=30,file='Fybalance.dat')
      else 
       open(unit=30,file='Fybalance.dat',position='append')
      end if

 30   format(e14.6,1i9,1000e36.24)

      write(6,30)time,m1,term1,term2,
     $            term3,term4,term5,term6,total,Error

      write(30,30)time,m1,term1,term2,
     $            term3,term4,term5,term6,total,Error
      close(30)

      return
      end
