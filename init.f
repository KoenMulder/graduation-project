      subroutine initialize
      include 'var.cmn'
      integer n
      real*8 Rid,temp,frt,j0,alphac

c     x1=x1L=0 at the electrode (left side of domain)
c     x1=x1R=lx1 at a location in the bulk flow (right side of domain)
c     domain length is lx1

      lx1=1d-4
      lx2=1d-4
      lx3=1d-4
 
      do i=-1,im+1
       x1_i(i)=dfloat(i)*lx1/dfloat(im)
      end do 
 
      do i=0,im+1
       x1(i)=0.5d0*(x1_i(i-1)+x1_i(i))
       dx1(i)=x1_i(i)-x1_i(i-1)
      end do

      do i=0,im
       dx1_i(i)=x1(i+1)-x1(i)
      end do
 
      do j=-1,jm+1
       x2_j(j)=dfloat(j)*lx2/dfloat(jm)
      end do 

      do j=0,jm+1
       x2(j)=0.5d0*(x2_j(j-1)+x2_j(j))
       dx2(j)=x2_j(j)-x2_j(j-1)
      end do

      do j=0,jm
       dx2_j(j)=x2(j+1)-x2(j)
      end do

      do k=-1,km+1
       x3_k(k)=dfloat(k)*lx3/dfloat(km)
      end do 
 
      do k=0,km+1
       x3(k)=0.5d0*(x3_k(k-1)+x3_k(k))
       dx3(k)=x3_k(k)-x3_k(k-1)
      end do

      do k=0,km
       dx3_k(k)=x3(k+1)-x3(k)
      end do



c     bubble
      bubble=.true.
c     free ('f') or no ('n')
      slip='n'


      height=0d0
      chord=0d0
      radius=3d0*x1_i(1)
      radius0=radius

      drdt=0d0
  
      xc=radius+x1_i(1)
      yc=lx1/2d0
      zc=lx1/2d0


c     phiLe = potential of electrode plus standard potential 
c     phiL  = phiLe - eta 
c     phiR  = potential at right hand side

      phiLe=-0.6d0
      phiR=0d0

c     exchange current density
      j0=1d0
     
c     other kinetic constants
      kc0=-j0
      ka0=j0      
 
      Rid=8.3142d0
      Far=96485d0
      temp=353d0
      frt=Far/Rid/temp
      alphac=0.5d0

      ac=alphac*frt
      aa=(1d0-alphac)*frt



c     1: dissolved H2 concentration 
c     2: electrolyte concentration (1m: OHmin, 1p, Kplus)
c     3: H2O concentration


c     reference conditions to be applied in the kinetics
c     and also used as bulk flow boundary condition
 
      c_ref(1)=1d0*0.16d0
      c_ref(2)=1d0*6700d0
      c_ref(3)=1d0*49000d0
 
      c_s(1)=c_ref(1)
c     Diffusivities of K+ and OH- are assumed to be the same. 
c     Electrolyte concentration is high and nonequal Fickian diffusivities 
c     are not realistic and consistent with the mass averaged velocity.   

      dif(1)=5.8d-9
      dif(2)=3.2d-9
      dif(3)=dif(2)

      condfac=Far*frt*2d0*dif(2)

c     set initial conditions

      do k=0,km+1
      do j=0,jm+1 
       eta(j,k)=0d0
       phiL(j,k)=phiLe
      end do
      end do


      umax=1d-3

!$omp parallel do private(i,j,k)
      do k=0,km+1
      do j=0,jm+1
      do i=0,im+1
       phi(i,j,k)=phiR
       p(i,j,k)=0d0
       u1(i,j,k)=0d0
       u2(i,j,k)=4d0*umax*(x1(i)/lx1)*(1d0-(x1(i)/lx1))
       u3(i,j,k)=0d0
      end do
      end do
      end do

      do n=1,nmax
!$omp parallel do 
       do k=0,km+1
       do j=0,jm+1
       do i=0,im+1
        c(i,j,k,n)=c_ref(n)
       end do
       end do
       end do
      end do


 
     
      time=0d0
      ro=1258d0
      nu=6.7d-7
      c_w=1d-30+10d0*umax
      co=1d0/3d0

c     set the explicit time step, we assume 
c     that the maximum diffusivity is dif(1) 
      dtime=0.5d0*min(co*dx1_i(0)/c_w,0.25d0*(dx1_i(0)**2)/nu,   
     $      0.25d0*(dx1_i(0)**2)/dif(1))

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

      

      do i=1,36
       minbeta(i)=1d30
       maxbeta(i)=-1d30
      end do
 
      return
      end
     

