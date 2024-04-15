      subroutine IBM
      include 'var.cmn' 
      integer counter,xx_plus,xx_minus,
     $        yy_plus,yy_minus,zz_plus,zz_minus,
     $        nr,nyc,nzc,absdir,o,mm
      real*8  nx,ny,nz,xc0,udummy
      real*8  length_n

      change=1
      cube_number=0
      neighbor(1:im*jm*km/2,1:6)=0

      xc0=xc
      xc=radius+x1_i(1)

      dxcdt=(xc-xc0)/dtime

c     dangerous, do you realize that nr, yc and zc are equated to 
c     a real, and that Fortrans rounds of or takes just int? 
c     perhaps better to take a wider loop and to create a dist2 array
c     to make it faster

      nr=radius/dx1(0)
      nyc=yc/dx1(0)
      nzc=zc/dx1(0)

      nx_start=1
      nx_end=2*nr
      ny_start=nyc-nr+1
      ny_end=nyc+nr
      nz_start=nzc-nr+1
      nz_end=nzc+nr

      udummy=0d0

!$omp parallel do private(i,j,k) 
      do k=0,km+1
      do j=0,jm+1
      do i=0,im+1
       typ(i,j,k)=0
       typ1(i,j,k)=0
       typ2(i,j,k)=0
       typ3(i,j,k)=0
       typ1a(i,j,k)=0
       typ2a(i,j,k)=0
       typ3a(i,j,k)=0

      end do
      end do
      end do

!$omp parallel do private(i,j,k) 
      do k=1,km
      do j=1,jm
      do i=1,im
       if ((x1(i)-xc)**2+(x2(j)-yc)**2+(x3(k)-zc)**2.gt.radius**2) then
        typ(i,j,k)=1
       else
        u1(i-1,j,k)=udummy
        u1(i,j,k)=udummy
        u2(i,j-1,k)=udummy
        u2(i,j,k)=udummy
        u3(i,j,k-1)=udummy
        u3(i,j,k)=udummy
       end if
 
      end do
      end do
      end do
        
      do k=1,km
      do j=1,jm
      do i=1,im

       counter=0

       x1_plus(change)=0
       x1_minus(change)=0
       x2_plus(change)=0
       x2_minus(change)=0
       x3_plus(change)=0
       x3_minus(change)=0

       if (typ(i,j,k).eq.0) then

        cube_number=cube_number+1

        if (typ(i+1,j,k).eq.1) then
         counter=counter+1
         x1_plus(change)=1
        end if
        if (typ(i-1,j,k).eq.1) then
         counter=counter+1
         x1_minus(change)=1
        end if
        if (typ(i,j+1,k).eq.1) then
         counter=counter+1
         x2_plus(change)=1
        end if
        if (typ(i,j-1,k).eq.1) then
         counter=counter+1
         x2_minus(change)=1
        end if
        if (typ(i,j,k+1).eq.1) then
         counter=counter+1
         x3_plus(change)=1
        end if
        if (typ(i,j,k-1).eq.1) then
         counter=counter+1
         x3_minus(change)=1
        end if

        if (counter .ge. 1) then

         if (x1_plus(change).eq.1) then
          neighbor(change,1)=1
         end if

         if (x1_minus(change).eq.1) then
          neighbor(change,2)=1
         end if


         if (x2_plus(change).eq.1) then
          neighbor(change,3)=1
         end if

         if (x2_minus(change).eq.1) then
          neighbor(change,4)=1
         end if


         if (x3_plus(change).eq.1) then
          neighbor(change,5)=1
         end if

         if (x3_minus(change).eq.1) then
          neighbor(change,6)=1
         end if




         x_IB(change)=i
         y_IB(change)=j
         z_IB(change)=k

         nx=x1(i)-xc
         ny=x2(j)-yc
         nz=x3(k)-zc

         if (abs(nx).ge.max(abs(ny),abs(nz))) then 
          absdir=1
         else
          if (abs(ny).ge.abs(nz)) then
           absdir=2
          else
           absdir=3
          end if
         end if

         if (absdir.eq.3) then 
          if (nz .gt. 0d0) then

           if (k.le.km-1) then
            direction(change)=3
            direction_v(change)=3
           else 
            write(6,*)'IB problem z_plus'
            stop
           end if

          else 

           if (k .ge. 2) then
            direction(change)=-3
            direction_v(change)=-3
           else
            write(6,*)'IB problem z_minus'
            stop
           end if
 
          end if
         end if

         if (absdir.eq.2) then
          if (ny .gt. 0d0) then

           if (j .le. jm-1) then
            direction(change)=2
            direction_v(change)=2
           else
            write(6,*)'IB problem y_plus'
            stop
           end if

          else 

           if (j .ge. 2) then
            direction(change)=-2
            direction_v(change)=-2

           else
            write(6,*)'IB problem y_minus'
            stop
           end if

          end if
         end if

         if (absdir.eq.1) then
          if (nx .gt. 0d0) then

           if (i .le. im-1) then
            direction(change)=1
            direction_v(change)=1

           else
            write(6,*)'IB problem x_plus'
            stop
           end if

          else

           if (i .ge. 2) then
            direction(change)=-1
            direction_v(change)=-1

           else

            if (i.eq.1) then

             direction_v(change)=-1


             if (abs(nz).gt.abs(ny)) then

              if (nz.gt.0d0) then
               if (k.le.km-1) then
                direction(change)=3
               else
                write(6,*)'IB problem z_plus at i=1'
                stop
               end if
              else
               if (k.ge.2) then
                direction(change)=-3
               else
                write(6,*)'IB problem z_minus at i=1'
                stop
               end if
              end if

             else

              if (ny.gt.0d0) then
               if (j.le.jm-1) then
                direction(change)=2
               else
                write(6,*)'IB problem y_plus at i=1'
                stop
               end if
              else
               if (j.ge.2) then
                direction(change)=-2
               else
                write(6,*)'IB problem y_minus at i=1'
                stop
               end if
              end if
     
             end if

            else

             write(6,*)'IB problem x_minus for i>=2'
             stop

            end if

           end if

          end if

         end if

         change=change+1
        end if

       end if

      end do
      end do
      end do

 
      do k=1,km
      do j=1,jm
      do i=1,im

       if(typ(i,j,k)*typ(i+1,j,k).eq.1) then
         typ1(i,j,k)=1
       end if

       if(typ(i,j,k)*typ(i,j+1,k).eq.1) then
         typ2(i,j,k)=1
       end if

       if(typ(i,j,k)*typ(i,j,k+1).eq.1) then
         typ3(i,j,k)=1
       end if

       if(typ(i,j,k)+typ(i+1,j,k).ge.1) then
         typ1a(i,j,k)=1
       end if

       if(typ(i,j,k)+typ(i,j+1,k).ge.1) then
         typ2a(i,j,k)=1
       end if

       if(typ(i,j,k)+typ(i,j,k+1).ge.1) then
         typ3a(i,j,k)=1
       end if


      end do
      end do
      end do

      typ1(0,1:jm,1:km)=0
      typ1(im,1:jm,1:km)=0

      typ2(1:im,0,1:km)=0
      typ2(1:im,jm,1:km)=1

      typ3(1:im,1:jm,0)=0
      typ3(1:im,1:jm,km)=1



      mm=1
      do o=1,change-1
         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)

         if (direction(o) .eq. 1) then
            xb(o)=xc+sqrt(radius**2-( (x2(j)-yc)**2 + (x3(k)-zc)**2 )) 
            yb(o)=x2(j)
            zb(o)=x3(k)
            beta=(xb(o)-x1(i))/dx1_i(i)
            call beta_function(35)
            alfa_s_P(o)=alfa_s
            beta_s_P(o)=beta_s
            gama_s_P(o)=gama_s

         else if (direction(o) .eq. -1) then
            xb(o)=xc-sqrt(radius**2-( (x2(j)-yc)**2 + (x3(k)-zc)**2))
            yb(o)=x2(j)
            zb(o)=x3(k)
            beta=-(xb(o)-x1(i))/dx1_i(i)
            call beta_function(36)
            alfa_s_P(o)=alfa_s
            beta_s_P(o)=beta_s
            gama_s_P(o)=gama_s

         else if (direction(o) .eq. 2) then
            xb(o)=x1(i)
            yb(o)=yc+sqrt(radius**2-( (x1(i)-xc)**2 + (x3(k)-zc)**2 ))
            zb(o)=x3(k)
            beta=(yb(o)-x2(j))/dx2_j(j)
            call beta_function(33)
            alfa_s_P(o)=alfa_s
            beta_s_P(o)=beta_s
            gama_s_P(o)=gama_s


         else if (direction(o) .eq. -2) then
            xb(o)=x1(i)
            yb(o)=yc-sqrt(radius**2-( (x1(i)-xc)**2 + (x3(k)-zc)**2 ))
            zb(o)=x3(k)
            beta=-(yb(o)-x2(j))/dx2_j(j)
            call beta_function(34)
            alfa_s_P(o)=alfa_s
            beta_s_P(o)=beta_s
            gama_s_P(o)=gama_s


         else if (direction(o) .eq. 3) then
            xb(o)=x1(i)
            yb(o)=x2(j)
            zb(o)=zc+sqrt(radius**2-((x1(i)-xc)**2 + (x2(j)-yc)**2))
            beta=(zb(o)-x3(k))/dx3_k(k)
            call beta_function(31)
            alfa_s_P(o)=alfa_s
            beta_s_P(o)=beta_s
            gama_s_P(o)=gama_s


         else if (direction(o) .eq. -3) then
            xb(o)=x1(i)
            yb(o)=x2(j)
            zb(o)=zc-sqrt(radius**2-((x1(i)-xc)**2 + (x2(j)-yc)**2))
            beta=-(zb(o)-x3(k))/dx3_k(k)
            call beta_function(32)
            alfa_s_P(o)=alfa_s
            beta_s_P(o)=beta_s
            gama_s_P(o)=gama_s

         end if
     

          nx=2d0*(xb(o)-xc)
          ny=2d0*(yb(o)-yc)
          nz=2d0*(zb(o)-zc)

          length_n=sqrt(nx**2+ny**2+nz**2)

          nx_n=nx/length_n
          ny_n=ny/length_n
          nz_n=nz/length_n

          x_prob(o,mm)=xb(o)+dx1_i(0)*nx_n
          y_prob(o,mm)=yb(o)+dx1_i(0)*ny_n
          z_prob(o,mm)=zb(o)+dx1_i(0)*nz_n

          nx_prob(o,mm)=floor(x_prob(o,mm)/dx1_i(0)+0.5d0)
          ny_prob(o,mm)=floor(y_prob(o,mm)/dx2_j(0)+0.5d0)
          nz_prob(o,mm)=floor(z_prob(o,mm)/dx3_k(0)+0.5d0)

         if (nx_prob(o,mm).lt.0) then
         if (mm.eq.1) then
            nx_prob(o,mm)=0
            x_prob(o,mm)=x1(nx_prob(o,mm))
         end if
         end if

         denom(o,mm)=(x1(nx_prob(o,mm))-x1(nx_prob(o,mm)+1))
     $      *(x2(ny_prob(o,mm))-x2(ny_prob(o,mm)+1))
     $      *(x3(nz_prob(o,mm))-x3(nz_prob(o,mm)+1))


         f01(o,mm)=x1(nx_prob(o,mm)+1)*x2(ny_prob(o,mm)+1)
     $            *x3(nz_prob(o,mm)+1)
         f02(o,mm)=x1(nx_prob(o,mm)+1)*x2(ny_prob(o,mm)+1)
     $            *x3(nz_prob(o,mm))
         f03(o,mm)=x1(nx_prob(o,mm)+1)*x2(ny_prob(o,mm))
     $            *x3(nz_prob(o,mm)+1)
         f04(o,mm)=x1(nx_prob(o,mm)+1)*x2(ny_prob(o,mm))
     $            *x3(nz_prob(o,mm))
         f05(o,mm)=x1(nx_prob(o,mm))*x2(ny_prob(o,mm)+1)
     $            *x3(nz_prob(o,mm)+1)      
         f06(o,mm)=x1(nx_prob(o,mm))*x2(ny_prob(o,mm)+1)
     $            *x3(nz_prob(o,mm))  
         f07(o,mm)=x1(nx_prob(o,mm))*x2(ny_prob(o,mm))
     $            *x3(nz_prob(o,mm)+1)  
         f08(o,mm)=x1(nx_prob(o,mm))*x2(ny_prob(o,mm))
     $            *x3(nz_prob(o,mm))   

         f11(o,mm)=x2(ny_prob(o,mm)+1)*x3(nz_prob(o,mm)+1)
         f12(o,mm)=x2(ny_prob(o,mm)+1)*x3(nz_prob(o,mm))
         f13(o,mm)=x2(ny_prob(o,mm))*x3(nz_prob(o,mm)+1)
         f14(o,mm)=x2(ny_prob(o,mm))*x3(nz_prob(o,mm))

         f21(o,mm)=x1(nx_prob(o,mm)+1)*x3(nz_prob(o,mm)+1)
         f22(o,mm)=x1(nx_prob(o,mm)+1)*x3(nz_prob(o,mm))
         f23(o,mm)=x1(nx_prob(o,mm))*x3(nz_prob(o,mm)+1)
         f24(o,mm)=x1(nx_prob(o,mm))*x3(nz_prob(o,mm))

         f31(o,mm)=x1(nx_prob(o,mm)+1)*x2(ny_prob(o,mm)+1)
         f32(o,mm)=x1(nx_prob(o,mm)+1)*x2(ny_prob(o,mm))
         f33(o,mm)=x1(nx_prob(o,mm))*x2(ny_prob(o,mm)+1)
         f34(o,mm)=x1(nx_prob(o,mm))*x2(ny_prob(o,mm))

      end do


      return
      end




      subroutine bulb_velocity
      include 'var.cmn' 
      integer o

      do o=1,change-1

         i=x_IB(o)
         j=y_IB(o)
         k=z_IB(o)


         if (direction_v(o) .eq. 1) then
           call x_plus

         else if (direction_v(o) .eq. -1) then
          call x_minus

         else if (direction_v(o) .eq. 2) then
          call y_plus

         else if (direction_v(o) .eq. -2) then
          call y_minus

         else if (direction_v(o) .eq. 3) then
          call z_plus

         else if (direction_v(o) .eq. -3) then
          call z_minus

         end if         


      end do


      return
      end
