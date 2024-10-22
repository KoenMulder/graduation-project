subroutine write_output()
    !!------------------------------------------------------------------------------
    !! Write simulation results for bubble only to bubbleInfo.dat
    !!------------------------------------------------------------------------------
    !! DESCRIPTION:
    !!>  The following parameters are lised in the bubbleInfo.dat file
    !!>  of the simulation:
    !! - time and time step, 
    !! - Bubble location and velocity
    !! - Net forces acting on the bubble
    !! - Mass flux and bubble growth rate
    use var2
    implicit none

666 format(19E36.24,1L10,5I10)
    if ((m1.eq.1).and.(m2.eq.1)) then
        open(unit=999,file='./output/BubbleInfo.dat')
333     format(20A36,10A6)
        write(999,333) 'time', 'dtime','xc', 'yc', 'zc', 'massflux',    &
        'drdt', 'radius', 'dxcdt', 'dycdt', 'dzcdt', 'H2_mass',         &
        'Fsx', 'Fsy', 'Fsz', 'Fvx', 'Fvy', 'Fvz', 'Fby',                &
        'isDetached', 'nIBc', 'nIBFx', 'nIBFy', 'nIBFz', 'nBtot'

        write(999,666) time, dtime, xc, yc, zc, bubbleMassflux, drdt,   &
        radius, dxcdt, dycdt, dzcdt, balance_h2, fsx, fsy, fsz, fvx,    &
        fvy, fvz, fby, isDetached, change-1, o_u, o_v, o_w, cube_number
        close (999,status='keep')
    end if
    if ((m1.gt.1).and.(m2.eq.m2_max)) then
        open(unit=999,file='./output/BubbleInfo.dat',access='append')
        write(999,666) time, dtime, xc, yc, zc, bubbleMassflux, drdt,   &
        radius, dxcdt, dycdt, dzcdt, balance_h2, fsx, fsy, fsz, fvx,    &
        fvy, fvz, fby, isDetached, change-1, o_u, o_v, o_w, cube_number
        close (999,status='keep')
    end if
end subroutine write_output


subroutine write_inputParam()
    !!------------------------------------------------------------------------------
    !! Write input parameters to screen and file inputParam.dat
    !!------------------------------------------------------------------------------
    !! DESCRIPTION:
    !!>  The following parameters are lised in the inputParam.dat file at the start
    !!>  of the simulation:
    !! - Universal speed of sound
    !! - Density of hydrogen gas and KOH-solution
    !! - Kinematic viscosity
    use var2
    implicit none
    real(8) :: u2_inf, Reb
    
    usou2 = 0.05d0*(dx1(0)/dtime)**2

    u2_inf = 2*(1 - rhoH2/rhoKOH)*R_max**2*g/(9*nuKOH)
    Reb    = 2*R_max*abs(u2_inf)/nuKOH

    ! Write to screen
    write(*,'(A55)') '------------------ Current  settings ------------------'

1   format(1A30, 1E20.10, 1A5)
    write(*,1) 'Time step: ', dtime, ' [s]'
    write(*,1) 'Saving each: ', dtime*m2_max, ' [s]'
    write(*,1) 'Max. Simulated time: ', dtime*m2_max*m1_max, ' [s]'
    write(*,1) 'Universal Speed of sound: ', sqrt(usou2), ' [m/s]'
    write(*,1) 'Bubble Reynolds number: ', Reb, ' [-]'
    write(*,1) 'Virtual Force Coeff. Cv: ', Cv

2   format(1A30, 1L5)
    write(*,2) 'Simulate bubble: ', bubble
    if ( bubble ) then
        write(*,1)  'R-critical: ', R_max, '[m]'
4       format(1A30,5A20)
        write(*,4) '', 'Bubble ID', 'R-initial [m]', 'X-initial [m]', 'Y-initial [m]', 'Z-initial [m]'
6       format(1A30, 1I20, 4E20.10)
        write(*,6) 'Bubble info ', 0, radius, xc, yc, zc
    end if

!     write(*,2) 'Use previous simdata: ', loadSimdata
!     if ( loadSimdata ) then
!         write(pathDataDump,'(3A,I4.4,A4)') trim(dataPath), trim(fileName),  &
!             '_', m1CompStart,'.dat'
! 7   format(1A30, A)
!         write(*,7) 'Loading initial data from: ', trim(pathDataDump)
!         call initLoadPreviousSimulation
!     end if
    write(*,*) '------------------ Running Simulation -----------------'

    ! Write to file
    ! open(unit=1, file='./output/inputParam.dat')
    ! write(1,)

    
end subroutine write_inputParam