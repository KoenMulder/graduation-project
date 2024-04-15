      program elec
      use omp_lib
      include 'var.cmn'

      integer restart,i1,j1,k1,do_artcom

      real*8 wtime,maxuabs,u1cen,u2cen,u3cen,uabs,usou2

      real*8 Error_v(0:im+1,0:jm+1,0:km+1),Max_Err_v,
     $ v_old(0:im+1,0:jm+1,0:km+1),ER_v

      real*8 tmax,
     $ res(1:im,1:jm,1:km)



      call omp_set_num_threads(48)     
      wtime=omp_get_wtime()

      call initialize


      restart=0
      if (restart.eq.1) then

        open(unit=111,file='datf1',form='unformatted')
        read(111)time,p,u1,u2,u3,c,eta,phi,total,radius,
     $          height,xc,ka,kc,phiL,radius0,conduct
     $          ,typ,typ1,typ2,typ3,drdt

        close(111)
      end if

      do_artcom=1
      if (do_artcom.eq.1) then
       usou2=0.05d0*(dx1(0)/dtime)**2
       write(6,*)'usou',sqrt(usou2)
      end if

      call eval_time



      call bou_u
      call bou_w
      call bou_v


      do m1=1,1000
       do m2=1,700

       

          call IBM

          call bulb_velocity

          call u
          call v
          call w


          call bou_rhu
          call bou_rhw
          call bou_rhv
 
          if (do_artcom.eq.0) then
           call pre_P
           call pois
          end if

          if (0d0.eq.1d0) then

          ER_v=1e-8
          Max_Err_v=5d0
          iter1=0
          do while (Max_Err_v .gt. ER_v)
           do k=0,km+1
           do j=0,jm+1
           do i=0,im+1
             v_old(i,j,k)=p(i,j,k)
           end do
           end do
           end do

           call pressure
           call bou_P
           call bulb_p
           iter1=iter1+1
           do k=0,km+1
           do j=0,jm+1
           do i=0,im+1

             Error_v(i,j,k)=abs(p(i,j,k)-v_old(i,j,k))       
               
           end do
           end do
           end do


           Max_Err_v=maxval(Error_v)
c           if (iter1.gt.50000) then
c            write(*,*)'Iteration problem in P',Max_Err_v
c            stop
c           end if

          end do
          end if



      if (0d0.eq.1d0) then
      do k=1,km
      do j=1,jm
      do i=1,im
       
       res(i,j,k)=-p(i,j,k)
     $          -D1(i,j,k)*p(i-1,j,k)
     $          -D2(i,j,k)*p(i+1,j,k)
     $          -D3(i,j,k)*p(i,j-1,k)
     $          -D4(i,j,k)*p(i,j+1,k)
     $          -D5(i,j,k)*p(i,j,k-1)
     $          -D6(i,j,k)*p(i,j,k+1)+
     $          work(i,j,k)

      end do
      end do
      end do


          do k=1,km
          do j=1,jm 
          do i=1,im

           if (typ(i,j,k).eq.0) then
            res(i,j,k)=0d0
           end if

          end do
          end do
          end do

          tmax=0d0

          do k=1,km
          do j=1,jm 
          do i=1,im
  
           tmax=max(tmax,abs(res(i,j,k)))
           if (tmax.eq.abs(div(i,j,k))) then
            i1=i
            j1=j
            k1=k
           end if

          end do
          end do
          end do


          write(*,*)'Res',tmax,i1,j1,k1
       end if



          call velocity

          call bou_u
          call bou_w
          call bou_v

!omp parallel do
          do k=1,km
          do j=1,jm 
          do i=1,im
           div(i,j,k)=(u1(i,j,k)-u1(i-1,j,k))/dx1(i)+
     $             (u2(i,j,k)-u2(i,j-1,k))/dx2(j)+ 
     $             (u3(i,j,k)-u3(i,j,k-1))/dx3(k) 
          end do
          end do
          end do

          if (do_artcom.eq.1) then
           call bou_p
!omp parallel do
           do k=1,km+1
           do j=1,jm+1
           do i=1,im+1
            if (typ(i,j,k).eq.1) then
             p(i,j,k)=p(i,j,k)-ro*usou2*dtime*div(i,j,k)
            end if
           end do
           end do
           end do
          end if


           if (1d0.eq.1d0) then
          call prep_phi

          ER_v=1e-7
          Max_Err_v=5d0
          Max_Err_phi=5d0
          iter2=0
          do while (Max_Err_v.gt.ER_v .or. Max_Err_phi.gt.1.01d0*ER_v)
           do k=0,km+1
           do j=0,jm+1
           do i=0,im+1
             v_old(i,j,k)=phi(i,j,k)
           end do
           end do
           end do
           call spatial_phi
           call bou_phi
           iter2=iter2+1
           do k=0,km+1
           do j=0,jm+1
           do i=0,im+1

             Error_v(i,j,k)=abs(phi(i,j,k)-v_old(i,j,k))       
               
           end do
           end do
           end do


           Max_Err_v=maxval(Error_v)
c           if (mod(iter2,100).eq.0) then
c            write(6,30)'iter2',time,Max_Err_phi,Max_Err_v,iter2
c 30         format(a6,e16.8,2e12.4,i8)
c           end if
c           write(*,*)Max_Err_phi,Max_Err_v
c           if (iter2.gt.5000) then
c            write(*,*)'Iteration problem in Phi'
c            stop
c           end if


          end do
          end if



c          call over
          call spatial_c


          call bou_c
   
          time=time+dtime

          call bulb_flux

         call growing


       end do

       call eval_time

       if (1d0.eq.1d0) then

       do k=1,km
       do j=1,jm 
       do i=1,im
        div(i,j,k)=(u1(i,j,k)-u1(i-1,j,k))/dx1(i)+
     $             (u2(i,j,k)-u2(i,j-1,k))/dx2(j)+ 
     $             (u3(i,j,k)-u3(i,j,k-1))/dx3(k) 
       end do
       end do
       end do


       do k=1,km
       do j=1,jm 
       do i=1,im

        if (typ(i,j,k).eq.0) then
         div(i,j,k)=0d0
        end if

       end do
       end do
       end do



       tmax=0d0

       do k=1,km
       do j=1,jm 
       do i=1,im
  
        tmax=max(tmax,abs(div(i,j,k)))

        if (tmax.eq.abs(div(i,j,k))) then
         i1=i
         j1=j
         k1=k
        end if

       end do
       end do
       end do

      if(1d0.eq.1d0) then
       if (m1.eq.1) then
        open(unit=31,file='Div.dat')
       else 
        open(unit=31,file='Div.dat',position='append')
       end if
 31   format(e14.6,e14.6,3i6)
       write(31,31)time,tmax,i1,j1,k1
      close (31)


       if (m1.eq.1) then
       open(unit=320,file='p_iter.dat')
       open(unit=321,file='phi_iter.dat')
       open(unit=322,file='phi_bulb_iter.dat')
       open(unit=323,file='c_bulb_iter.dat')

       else 
       open(unit=320,file='p_iter.dat',position='append')
       open(unit=321,file='phi_iter.dat',position='append')
       open(unit=322,file='phi_bulb_iter.dat',position='append')
       open(unit=323,file='c_bulb_iter.dat',position='append')
   
       end if

 320    format(e14.6,3i6,3i6)

        write(320,320)time,m1,iter1
        write(321,320)time,m1,iter2
        write(322,320)time,m1,iter3
        write(323,320)time,m1,iter4

       close (320)
       close (321)
       close (322)
       close (323)
       end if
       end if




      end do  

       open(unit=111,file='datf',form='unformatted')
       write(111)time,p,u1,u2,u3,c,eta,phi,total,radius,
     $          height,xc,ka,kc,phiL,radius0,conduct,
     $          typ,typ1,typ2,typ3,drdt
       close(111) 



    
 

      call eval_prof

c      call bulb_flux_write
c      call domain_flux

c      call bulb_flux_force
c      call domain_flux_force




      wtime=omp_get_wtime()-wtime
      write(*,*)'Wtime=',wtime

      stop
      end 


