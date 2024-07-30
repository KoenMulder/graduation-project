program ElectrolyserSimulation
    !!------------------------------------------------------------
    !!  **Purpose:**
    !!
    !!      CFD-simulation of electrolyser with (multiple) bubble 
    !!      growth and bubble detachment using Immersed Boundary
    !!      Method and Artificial Compressibility for pressure.
    !!
    !!  **Info:**
    !!      - Input parameters are listed in the paramdata-module.
    !!      - Bubbles can be turned off, resulting in only the 
    !!      concentration transport in the domain.
    !!      - Uses openMP for parallelization.
    !!
    !!  **Credits:**
    !!
    !!      Code based on the implementation of MSc. F. 
    !!      Khalighi.
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    use omp_lib
    implicit none
    real(8) :: wtime

    call omp_set_num_threads(48)
    wtime = omp_get_wtime()

    ! Initialize data and apply initical conditions
    call initialize

    ! Apply boundary conditions
    call boundaryConditions

    ! Define bubble(s)
    call bubbleInitialize

    ! Write info to screen
    call writeStart

    ! Main loop
1   format(1A30,1I5)
    do m1 = m1Start,mComp
        do m2 = 1,mSave
            call writeOutput(m1,m2,0)

            ! Identify cells and faces with corresponding bubble
            call IBMidentify

            ! Apply BC to IB-cells (Immersed Boundary Method)
            call IBMboucond

            ! Central differencing on velocity momentum equation
            call NSmomentum

            ! potential and species concentration (mass) transport 
            call Speciestransport

            ! Compute mass flux H2 into bubble and determine new size
            call bubbleMain

            ! Bubble detachment model

            totaltime = totaltime + dtime

! <<<<--------- Verify the above before continuing --------->>>> !
            ! Todo: calculate forces on bubble
            ! Todo: check detachment criterion
            ! Todo: if detached: 
                ! Todo: Check 1 cell distance between walls condition
                ! Todo: if 1 cell distance present -> do nothing
        end do
        
        write(*,1) 'm1 iteration status: ', m1

        if ( m1.gt.m1Start ) then
            call writeOutput_SimData(m1)
        end if
        

        if ( (m1.eq.mComp).or.(BubbleBlock(1)%radius.ge.radiusCrit) ) then
            call writeOutput_IBM
            call writeOuput_Scalar
            call writeOutput_Velocity
            call writeOutput_NodePoint

2           format(1A30,1E20.10,1A10)
            write(*,2) 'Stopping at time: ', totaltime, ' [s]'
            wtime = omp_get_wtime() - wtime
            write(*,2) 'Total computing time: ', wtime, ' [s]'

            stop
        end if
    end do
end program ElectrolyserSimulation