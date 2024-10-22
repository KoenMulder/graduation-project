program elec
    use var2
    use omp_lib
    implicit none
    real(8)  :: wtime, comptime

    call omp_set_num_threads(nThreads)     
    wtime = omp_get_wtime()

    call initialize

    call write_inputParam

    call bou_u
    call bou_w
    call bou_v

    do m1=1,m1_max
      comptime = omp_get_wtime()
      do m2=1,m2_max

        call write_output

        call IBM

        call bulb_c
        call finding_u_IB
        call finding_v_IB
        call finding_w_IB

        call probe_IB_u
        call bulb_u
        
        call probe_IB_v 
        call bulb_v

        call probe_IB_w
        call bulb_w

        ! ---> contained in NSmomentum in F90 version 
        call u
        call v
        call w

        call bou_rhu
        call bou_rhw
        call bou_rhv

        call velocity

        call bou_u
        call bou_w
        call bou_v

!omp parallel do
        do k=1,km
          do j=1,jm
            do i=1,im
              div(i,j,k)=(u1(i,j,k)-u1(i-1,j,k))/dx1(i)+ &
                      (u2(i,j,k)-u2(i,j-1,k))/dx2(j)+ &
                      (u3(i,j,k)-u3(i,j,k-1))/dx3(k) 
            end do
          end do
        end do
        storage_p=p
!omp parallel do
        do k=1,km
          do j=1,jm
            do i=1,im
              if (typ(i,j,k).eq.1) then
                p(i,j,k) = storage_p(i,j,k) - rhoKOH*usou2*dtime*div(i,j,k)
              end if
            end do
          end do
        end do
        call bou_p
        ! ---> end NSmomentum
        

        ! ---> start species mass transport
        call prep_phi
        call spatial_phi
        call bou_phi
        call spatial_c
        call bou_c
        ! ---> end species mass transport

        time = time + dtime

        call bubble_Main

      end do

      comptime = omp_get_wtime() - comptime
1       format(1A30,1I5,1A11,1F11.5,1A4)
      write(*,1) 'Completed m1: ', m1, ' at wtime: ', comptime, ' [s]'

      if ( time.gt.0.125d0 ) then
        stop 'Simulation completed'
      end if

    end do
    write(*,*) 'End of computation'
    wtime=omp_get_wtime()-wtime
    write(*,*)'Wtime=',wtime

    stop
end program elec

