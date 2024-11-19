
subroutine write_inputParam()
    !!------------------------------------------------------------------------------
    !! Write input parameters to screen and file inputParam.dat
    !!>  The following parameters are lised in the inputParam.dat file at the start
    !!>  of the simulation:
    !! - Universal speed of sound
    !! - Density of hydrogen gas and KOH-solution
    !! - Kinematic viscosity
    use var2
    implicit none
    real(8) :: u2_inf, Reb

    u2_inf = 2d0*(rhoH2/rhoKOH - 1d0)*R_max**2d0*g/(9d0*nuKOH)
    Reb    = 2d0*radius*abs(u2_inf)/nuKOH

1   format(1A30, 1E20.10, 1A5)
    write(*,*) ! new line
    write(*,'(A55)') '------------------ Grid  settings ------------------'
    write(*,1) 'Grid resolution: DeltaX = ', lx1/dble(im), ' [m]'
    write(*,1) 'Grid resolution: DeltaY = ', lx2/dble(jm), ' [m]'
    write(*,1) 'Grid resolution: DeltaZ = ', lx3/dble(km), ' [m]'
    write(*,*) ! new line
    write(*,'(A55)') '--------------- Simulation  settings ---------------'
    write(*,1) 'Time step: ', dtime, ' [s]'
    write(*,1) 'Saving each: ', dtime*m2_max, ' [s]'
    write(*,1) 'Max. Simulated time: ', dtime*m2_max*m1_max, ' [s]'
    write(*,1) 'Universal Speed of sound: ', sqrt(usou2), ' [m/s]'
    write(*,1) 'Density ratio: ', rhoH2/rhoKOH, ' [-]'
    write(*,1) 'Galilei number: ', sqrt(abs(rhoH2/rhoKOH - 1d0)*abs(g)*(2d0*R_max)**3d0)/nuKOH, '[-]'
    write(*,1) 'Virtual Force Coeff. Cv: ', Cv, '[-]'

    SELECT CASE (FlowCondPreset)
        CASE ('Stokes')
            write(*,1) 'Stokes terminal velocity: ', u2_inf, ' [m/s]'
            write(*,1) 'Stokes Bubble-Reynolds: ', Reb, ' [-]'
    END SELECT 

2   format(1A30, 1L5)
    if ( useRungeKutta.eqv..true.) then
        write(*,2) 'Low storage Runge-Kutta: ', useRungeKutta
    end if

    write(*,2) 'Simulate bubble: ', bubble
    if ( bubble ) then
        write(*,*) ! new line
        write(*,'(A55)') '----------------- Bubble  settings -----------------'
        write(*,1) 'Critical radius: ', R_max, '[m]'
        write(*,'(1A30, 1I20, 1A5)') 'Number of bubbles: ', nBubbles, '[-]'
4       format(1A30,5A20)
        write(*,4) 'Bubble info: ', 'Bubble ID', 'R-initial [m]', 'X-initial [m]', 'Y-initial [m]', 'Z-initial [m]'
6       format(1A30, 1I20, 4E20.10)
        write(*,6) '', 0, radius, xc, yc, zc
    end if

    write(*,*) ! new line
    write(*,'(A55)') '---------------- Running Simulation ----------------'

end subroutine write_inputParam


subroutine write_output()
    !!------------------------------------------------------------------------------
    !! Write simulation results for bubble only to bubbleInfo.dat
    !!>  The following parameters are lised in the bubbleInfo.dat file
    !!>  of the simulation:
    !! - time and time step, 
    !! - Bubble location and velocity
    !! - Net forces acting on the bubble
    !! - Mass flux and bubble growth rate
    use var2
    implicit none

666 format(19E36.24,1L36)!,5I10)
    if ((m1.eq.1).and.(m2.eq.1)) then
        open(unit=999,file='./output/BubbleInfo.dat')
333     format(20A36,10A6)
        write(999,333) 'time', 'dtime','xc', 'yc', 'zc', 'massflux',    &
        'drdt', 'radius', 'dxcdt', 'dycdt', 'dzcdt', 'H2_mass',         &
        'Fsx', 'Fsy', 'Fsz', 'Fvx', 'Fvy', 'Fvz', 'Fby',                &
        'isDetached' !, 'nIBc', 'nIBFx', 'nIBFy', 'nIBFz', 'nBtot'

        write(999,666) time, dtime, xc, yc, zc, bubbleMassflux, drdt,   &
        radius, dxcdt, dycdt, dzcdt, balance_h2, fsx, fsy, fsz, fvx,    &
        fvy, fvz, fby, isDetached !, change-1, o_u, o_v, o_w, cube_number
        close (999,status='keep')
    end if
    if ((m1.gt.1).and.(m2.eq.m2_max)) then
        open(unit=999,file='./output/BubbleInfo.dat',position='append')
        write(999,666) time, dtime, xc, yc, zc, bubbleMassflux, drdt,   &
        radius, dxcdt, dycdt, dzcdt, balance_h2, fsx, fsy, fsz, fvx,    &
        fvy, fvz, fby, isDetached !, change-1, o_u, o_v, o_w, cube_number
        close (999,status='keep')
    end if

    if ( FullDomainOutput.eqv..true.) then
        ! Only save at end of m2-loop
        if ( m2.eq.m2_max ) then
            ! Write output of extrapolated node points parameters
            call write_output_CellNodePoint()

            ! Write output of all scalar parameters
            call write_output_CellCenterPoints()

            ! Write output of all velocity parameters
            call write_output_CellFacePoints()
        end if
    end if

    if ( WriteToBinaryOutput.eqv..true. ) then
        if ( m2.eq.m2_max ) then
            call write_output_UnformattedBinary()
        end if
    end if
end subroutine write_output




subroutine write_output_CellNodePoint()
    use var2
    implicit none
    character(len=60) :: filePath
   
    ! Storing all cell and face values at the node points
    do k=0,km
        do j=0,jm
            do i=0,im
              u_node(i,j,k) = 0.25d0*(u1(i,j,k)   + u1(i,j+1,k)   &
                                    + u1(i,j,k+1) + u1(i,j+1,k+1) )
              v_node(i,j,k) = 0.25d0*(u2(i,j,k)   + u2(i+1,j,k)   &
                                    + u2(i,j,k+1) + u2(i+1,j,k+1) )
              w_node(i,j,k) = 0.25d0*(u3(i,j,k)   + u3(i+1,j,k)   &
                                     +u3(i,j+1,k) + u3(i+1,j+1,k) )
              p_node(i,j,k) = (1d0/8d0)*(p(i,j,k) + p(i+1,j,k)    &
                                     + p(i,j+1,k) + p(i+1,j+1,k)  &
                                     + p(i,j,k+1) + p(i+1,j,k+1)  &
                                     + p(i,j+1,k+1) + p(i+1,j+1,k+1))
            end do
        end do
    end do

    ! Define filename
1   format(1A,1I4.4,1A4)
    write(filePath,1) './output/NodeData_', m1, '.dat'

    ! Write data to file
    open(111,file=trim(filePath))
2   format(7A36)
    write(111,2) 'x1_i', 'x2_j', 'x3_k', 'u_node', 'v_node', 'w_node', 'p_node'
3   format(7E36.24)
    do k = 0,km
        do j = 0,jm
            do i = 0,im
                write(111,3) x1_i(i), x2_j(j), x3_k(k), u_node(i,j,k),  &
                            v_node(i,j,k), w_node(i,j,k), p_node(i,j,k)
            end do
        end do
    end do
    close(111,status='keep')

end subroutine write_output_CellNodePoint




subroutine write_output_CellCenterPoints()
    use var2
    implicit none
    character(len=60) :: filePath
    integer :: n

    ! Define filename
1   format(1A,1I4.4,1A4)
    write(filePath,1) './output/CenterData_', m1, '.dat'

    ! Write data to file
    open(111,file=trim(filePath))
2   format(8A36)
    write(111,2) 'x1', 'x2', 'x3', 'cH2', 'cKOH', 'cH2O', 'phi', 'p'

3   format(8E36.24)
    do k=0,km+1
        do j=0,jm+1
            do i=0,im+1
                write(111,3) x1(i), x2(j), x3(k), (c(i,j,k,n), n=1,nmax), phi(i,j,k), p(i,j,k)
            end do
        end do
    end do
    close(111,status='keep')

end subroutine write_output_CellCenterPoints




subroutine write_output_CellFacePoints()
    use var2
    implicit none
    character(len=60) :: filePath
    ! Velocities u, v, and w are saved subsequently in the same array
    ! with 'direction' 1,2 and 3 respectively. 

    ! Define filename
1   format(1A,1I4.4,1A4)
    write(filePath,1) './output/FaceData_', m1, '.dat'

    ! Write data to file
    open(111,file=trim(filePath))
2   format(5A36)
    write(111,2) 'direction', 'xLocation', 'yLocation', 'zLocation', 'magnitude'

3   format(I36,4E36.24)
    do k = 0,km+1
        do j = 0,jm+1
            do i = 0,im+1
                write(111,3) 1, x1_i(i), x2(j), x3(k), u1(i,j,k)
                write(111,3) 2, x1(i), x2_j(j), x3(k), u2(i,j,k)
                write(111,3) 3, x1(i), x2(j), x3_k(k), u3(i,j,k)
            end do
        end do
    end do

    close(111,status='keep')
    
end subroutine write_output_CellFacePoints




subroutine write_output_UnformattedBinary()
    use var2
    implicit none
    character(len=60) :: filePath

    ! Define filename
1   format(1A,1I4.4,1A4)
    write(filePath,1) './output/simData_', m1, '.bin'

    open(111,file=trim(filePath),form='unformatted')
    write(111) time, x1, x1_i, x2, x2_j, x3, x3_k, p, u1, u2,    &
               u3, c, eta, phi, ka, kc, phiL, conduct, xc, yc,   &
               zc, radius, drdt, dxcdt, dycdt, dzcdt, ddxcdttFv, &
               ddycdttFv, ddzcdttFv, isDetached, fsx, fsy, fsz,  &
               bubbleMassflux
    close(111, status='keep')

    
end subroutine write_output_UnformattedBinary