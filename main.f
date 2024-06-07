      program elec
      use omp_lib
      include 'var.cmn'

      integer restart,i1,j1,k1,do_artcom,i_num
      integer iii,jjj,kkk,o,k_num,max_i,max_j

      real*8 wtime,maxuabs,u1cen,u2cen,u3cen,uabs,usou2

      real*8 Error_v(0:im+1,0:jm+1,0:km+1),Max_Err_v,
     $ v_old(0:im+1,0:jm+1,0:km+1),ER_v

      real*8 tmax,test4,
     $ res(1:im,1:jm,1:km),test1,test2,test3



      call omp_set_num_threads(24)     
      wtime=omp_get_wtime()

c     Initialize domain arrays with initial values
c     - x,y,z grid created
c     - potential electrode (phiL) and (eta)
c     - pressure, velocity and species concentration
c     - determine kinetic prefactors

      call initialize
      dd=0d0

c     Restart-check function. If the code breaks/does not converge,
c     use the last known converged values and start from there.
c     --> it either starts the simulation from time zero, or uses the data from the most recent saved values
      restart=0
      if (restart.eq.1) then
        open(unit=111,file='datf1',form='unformatted')
        read(111)time,p,u1,u2,u3,c,eta,phi,total,radius,
     $          height,xc,ka,kc,phiL,radius0,conduct
     $          ,typ,typ1,typ2,typ3,drdt
        close(111)
      end if

c     Use the artificial compressibility method to decouple the pressure term of the Navier Stokes equation
c     usou2 is the speed of sound?
      do_artcom=1
      if (do_artcom.eq.1) then
        usou2=0.05d0*(dx1(0)/dtime)**2
        write(6,*)'usou',sqrt(usou2)
      end if

c     Evaluate the elapsed time for saving and post-processing. Not needed for computation.
      call eval_time

c     Apply boundary conditions for flow velocity in x,y,z direction
      call bou_u
      call bou_w
      call bou_v

c     Loop over number of time steps m1 and m2, and save results at each m1 (eval_time subroutine) 
      m1_max=2000
      do m1=1,m1_max
        do m2=1,700         
c         Use Inverse Boundary Method to probe point normal to bubble surface
          call IBM
          do rung=1,1
            call bulb_c
            call finding_u_IB
            call finding_v_IB
            call finding_w_IB

            do k_num=1,1
              call probe_IB_u
              call  bulb_u
              call probe_IB_v 
              call bulb_v
              call probe_IB_w
              call bulb_w
            end do

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
c           if we want to change the approach the velocity with the projection method
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
              end do
            end if
c           call velocity = add pressure to momentum equation
            call velocity
            call bou_u
            call bou_w
            call bou_v
c           Compute the divergence of the velocity == continuity equation
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
c           Call artificial compressibility subroutines for solving the navier stokes equation. Two ways of solving:
c           1. Artificial compressibility
c           2. Projection method
            if (do_artcom.eq.1) then
              if (rung.eq.1) then
                storage_p=p
              end if
c             usou2 == velocity of sound squared
c             ro == rho (density)
!omp parallel do
              do k=1,km
                do j=1,jm
                  do i=1,im
                    if (typ(i,j,k).eq.1) then
                      p(i,j,k)=storage_p(i,j,k)
     $                -ro*usou2*dtime*alpha_rung(rung)*div(i,j,k)
                    end if
                  end do
                end do
              end do
              call bou_p
c             call bulb_p
            end if
c           Solving potential equation iterative method
c           Look at the subroutine for the potential
            if (1d0.eq.1d0) then
              call prep_phi
c           Critertia error of the potential
              ER_v=1e-7
              Max_Err_v=5d0
              Max_Err_phi=5d0
              iter2=0
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
            end if
c           call concentration. Different time integration schemes are implemented: Euler method and the upwind method.
c           For current situation (no-slip) only consider the Euler method and disable the upwind method.
            if(1d0.eq.1d0) then
              call spatial_c
              call bou_c
            end if
          end do
          time=time+dtime
c         Compute bubble hydrogen flux to determine the new radius of the bubble.
          call bulb_flux
          call growing
        end do
c     C~~~~~~~~~~~~~~~~~~~~~~~~~~ End of computing time of m2-amount of time steps
c     C~~~~~~~~~~~~~~~~~~~~~~~~~~ Storing data m1-loop
c       Each m1 timestep (roughly) is saved 
c       eval_time is the save-to-file subroutine (storing subroutine)
        call bulb_flux
        call domain_flux
        call bulb_force_x
        call domain_flux_force_x
        call bulb_force_y
        call domain_flux_force_y
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
c     C~~~~~~~~~~~~~~~~~~~~~~~~~~ End of m1-loop
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
      write(*,*)dtime
      wtime=omp_get_wtime()-wtime
      write(*,*)'Wtime=',wtime

      stop
      end 


