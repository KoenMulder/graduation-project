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
        call finding_u_IB
        call finding_v_IB
        call finding_w_IB
        call probe_IB_u
        call probe_IB_v
        call probe_IB_w

        do rungIter = 1,rungIterEnd
          ! Apply bubble boundary conditions
          call bulb_c
          call bulb_u
          call bulb_v
          call bulb_w

          ! Solve momentum eq. without pressure
          call u
          call v
          call w
 ! compute divergence here 
          ! Apply momentum eq. boundary conditions
          call bou_rhu
          call bou_rhw
          call bou_rhv

          ! Compute fluid velocity with pressure p^n
          call velocity

          ! Apply fluid boundary conditions
          call bou_u
          call bou_w
          call bou_v

          ! Solve pressure eq. artificial comp.
          call p_solve_ArtComp
          call bou_p

          ! ---> start species mass transport
          SELECT CASE (FlowCondPreset)
            CASE ('Stokes')
              ! no computation for species required
              ! only flow
            CASE DEFAULT
              call prep_phi
              call spatial_phi
              call bou_phi
              call spatial_c
              call bou_c
          END SELECT
          ! ---> end species mass transport
        end do

        call bubble_Main

        time = time + dtime

      end do

      comptime = omp_get_wtime() - comptime
1     format(1A25,1I5,1A11,1F11.5,1A4)
      write(*,1) 'Completed m1: ', m1, ' at wtime: ', comptime, ' [s]'
    end do

    write(*,*) 'End of computation: m1 = ', m1
    wtime = omp_get_wtime() - wtime
    write(*,*) 'Wtime=', wtime
    
    stop
end program elec

