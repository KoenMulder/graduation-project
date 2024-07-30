subroutine writeStart()
    !!--------------------------------------------------------------
    !!  **Purpose:**
    !!      Write information to screen prior to start of main loop.
    !!
    !!--------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: n
    character(len=40) :: pathDataDump

    write(*,'(A55)') '------------------ Current  settings ------------------'

1   format(1A30, 1E20.10, 1A5)
    write(*,1) 'Time step: ', dtime, ' [s]'
    write(*,1) 'Saving each: ', dtime*mSave, ' [s]'
    write(*,1) 'Max. Simulated time: ', dtime*mSave*mComp, ' [s]'
    write(*,1) 'Universal Speed of sound: ', sqrt(usou2), ' [m/s]'
    write(*,1) 'Reynolds number: ', 2d0/3d0*umax*lx1/nuKOH, ' [-]'
    write(*,1) 'Peclet number: ', umax*lx1/DifKOH, '[-]'


2   format(1A30, 1L5)
    write(*,2) 'Simulate bubble: ', SimBubble
    if ( SimBubble ) then
3   format(1A30, 1i5)
        write(*,3) 'Number of bubbles: ', nBubbles

4       format(1A30,5A20)
        write(*,4) '', 'R-initial [m]', 'X-initial [m]', 'Y-initial [m]', 'Z-initial [m]', 'R-critical [m]'
5       format(1A30, 5E20.10)
        write(*,5) 'Bubble parameters: ', radius0, xCenter0, yCenter0, zCenter0, radiusCrit
6       format(1A30, 4E20.10)
        do n = 1, nBubbles
            write(*,6) 'BubbleBlock info:  ', BubbleBlock(n)%radius, BubbleBlock(n)%xCenter,  &
                                                BubbleBlock(n)%yCenter, BubbleBlock(n)%zCenter
        end do
        if ( SimDetach ) then
            write(*,2) 'Simulate detachment: ', SimDetach
        end if
    end if

    write(*,2) 'Use previous simdata: ', loadSimdata
    if ( loadSimdata ) then
        write(pathDataDump,'(3A,I4.4,A4)') trim(dataPath), trim(fileName),  &
            '_', m1CompStart,'.dat'
7   format(1A30, A)
        write(*,7) 'Loading initial data from: ', trim(pathDataDump)
        call initLoadPreviousSimulation
    end if
    write(*,*) '------------------ Running Simulation -----------------'

end subroutine writeStart




subroutine writeOutput(mIter1, mIter2, IterConverge)
    use commondata
    use paramdata
    implicit none
    integer, intent(in) :: mIter1, mIter2, IterConverge
    integer :: n, i, j, k
    integer, dimension(4) :: ErrorLocationIndex
    real(8) :: MaxError_c

    if ( SimBubble ) then
        ! Store Bubble parameter information
10      format(i10,13E36.24)
        if ( (mIter1.eq.1).and.(mIter2.eq.1) ) then
9           format(1A10,13A36)
            open(999,file='./output/BubbleInfo.dat')
            write(999,9) 'BubbleID', 'xCenter [m]', 'yCenter [m]', 'zCenter [m]', 'massFlux',   &
                        'drdt [m t-1]', 'radius [m]', 'u-COM [m s-1]', 'v-COM [m s-1]',         &
                        'w-COM [m s-1]', 'Force x []', 'Force y []', 'Force z []', 'time [s]'

            do n = 1, nBubbles
                write(999,10) BubbleBlock(n)%iBubble, BubbleBlock(n)%xCenter, BubbleBlock(n)%yCenter,   &
                                BubbleBlock(n)%zCenter, BubbleBlock(n)%massFlux, BubbleBlock(n)%drdt,   &
                                BubbleBlock(n)%radius, BubbleBlock(n)%u, BubbleBlock(n)%v,              &
                                BubbleBlock(n)%w, BubbleBlock(n)%ForceX, BubbleBlock(n)%ForceY,         &
                                BubbleBlock(n)%ForceZ, totaltime
            end do
            close(999,status='keep')
        end if
        if ( (mIter1.gt.1).and.(mIter2.eq.mSave) ) then
            open(999,file='./output/BubbleInfo.dat',access='append')
            do n = 1, nBubbles
                write(999,10) BubbleBlock(n)%iBubble, BubbleBlock(n)%xCenter, BubbleBlock(n)%yCenter,   &
                                BubbleBlock(n)%zCenter, BubbleBlock(n)%massFlux, BubbleBlock(n)%drdt,   &
                                BubbleBlock(n)%radius, BubbleBlock(n)%u, BubbleBlock(n)%v,              &
                                BubbleBlock(n)%w, BubbleBlock(n)%ForceX, BubbleBlock(n)%ForceY,         &
                                BubbleBlock(n)%ForceZ, totaltime
            end do
            close(999,status='keep')
        end if
    end if

    ! Store iterative convergence information
    MaxError_c = maxval(ErrIter_c)
    ErrorLocationIndex = maxloc(ErrIter_c)
    i = ErrorLocationIndex(1)
    j = ErrorLocationIndex(2)
    k = ErrorLocationIndex(3)
    n = ErrorLocationIndex(4)

8   format(1E36.24,6I12)
    if ( (mIter1.eq.1).and.(mIter2.eq.1) ) then
7       format(1A36,6A12)
        open(777,file='./output/IterConvergenceInfo.dat')
        write(777,7) 'Max. Error c', 'numb. Iter', 'Index i', 'Index j', 'Index k', 'IBtype', 'Species n'
        write(777,8)  MaxError_c, IterConverge, i, j, k, IBtype(i,j,k), n
        close(777,status='keep')
    end if

    if ( (mIter1.gt.1).and.(mIter2.eq.mSave) ) then
        open(777,file='./output/IterConvergenceInfo.dat', access='append')
        write(777,8)  MaxError_c, IterConverge, i, j, k, IBtype(i,j,k), n
        close(777,status='keep')
    end if

end subroutine writeOutput




subroutine writeOuput_Scalar()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k

    open(666,file='./output/initialScalar.dat')
6   format(8E36.24)
    do k=0,km+1
        do j=0,jm+1
            do i=0,im+1
                write(666,6) x1(i), x2(j), x3(k), c(i,j,k,1), c(i,j,k,2), c(i,j,k,3), phi(i,j,k), p(i,j,k)
            end do
        end do
    end do
    close(666,status='keep')
    
end subroutine writeOuput_Scalar




subroutine writeOutput_Velocity()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k

    ! Velocity components
    open(666,file='./output/u_velocity.dat')
6   format(4E36.24)
    do k = 0,km+1
        do j = 0,jm+1
            do i = 0,im+1
                write(666,6) x1_i(i), x2(j), x3(k), u1(i,j,k)
            end do
        end do
    end do
    close(666,status='keep')

    open(777,file='./output/v_velocity.dat')
    do k = 0,km+1
        do j = 0,jm+1
            do i = 0,im+1
                write(777,6) x1(i), x2_j(j), x3(k), u2(i,j,k)
            end do
        end do
    end do
    close(777,status='keep')

    open(888,file='./output/w_velocity.dat')
    do k = 0,km+1
        do j = 0,jm+1
            do i = 0,im+1
                write(888,6) x1(i), x2(j), x3_k(k), u3(i,j,k)
            end do
        end do
    end do
    close(888,status='keep')

    
end subroutine writeOutput_Velocity




subroutine writeOutput_IBM()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k, nIB, n

    ! Store bubble and IBM properties 
    if ( SimBubble ) then
        open(888,file='./output/IBMtype.dat')
8      format(3E36.24, I10, 3E36.24, I10)
        do k= 0,km+1
            do j=0,jm+1
                do i=0,im+1
                    write(888,8) x1(i), x2(j), x3(k), IBtype(i,j,k), dx1(i), dx2(j), dx3(k), IBbubble(i,j,k)
                end do
            end do
        end do
        close(888,status='keep')
        
        ! Storing cell faces IBM-types; center of the face along with the face size (width and lengths)
        ! Location of faces is different for u1, u2, and u3
9      format(3E36.24, I10, 2E36.24, I10)
        open(999,file='./output/IBMtypeU1.dat')
        do k = 0,km
            do j = 0,jm
                do i = 0,jm
                    write(999,9) x1_i(i), x2(j), x3(k), IBtypeU1(i,j,k), dx2_j(j), dx3_k(k), IBbubbleU1(i,j,k)
                end do
            end do
        end do
        close(999,status='keep')

        open(999,file='./output/IBMtypeU2.dat')
        do k = 0,km
            do j = 0,jm
                do i = 0,jm
                    write(999,9) x1(i), x2_j(j), x3(k), IBtypeU2(i,j,k), dx1_i(i), dx3_k(k), IBbubbleU2(i,j,k)
                end do
            end do
        end do
        close(999,status='keep')

        open(999,file='./output/IBMtypeU3.dat')
        do k = 0,km
            do j = 0,jm
                do i = 0,jm
                    write(999,9) x1(i), x2(j), x3_k(k), IBtypeU3(i,j,k), dx1_i(i), dx2_j(j), IBbubbleU3(i,j,k)
                end do
            end do
        end do
        close(999,status='keep')

        ! Store probe points
10      format(3E36.24,I4)
        open(1000,file='./output/probeXYZScalar.dat')
        do n = 1, nBubbles
            do nIB = 1, BubbleBlock(n)%nIB
                write(1000,10) probeCell(1,nIB,n), probeCell(2,nIB,n), probeCell(3,nIB,n), BubbleBlock(n)%iBubble
            end do
        end do
        close(1000,status='keep')
    end if
    
end subroutine writeOutput_IBM




subroutine writeOutput_NodePoint()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k
    real(8), allocatable, dimension(:,:,:):: u_node, v_node, w_node, p_node
   
    allocate(u_node(0:im,0:jm,0:km), v_node(0:im,0:jm,0:km),  &
             w_node(0:im,0:jm,0:km), p_node(0:im,0:jm,0:km)   )

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

    open(111,file='./output/NodeData.dat')
1   format(7E36.24)
    do k = 0,km
        do j = 0,jm
            do i = 0,im
                write(111,1) x1_i(i), x2_j(j), x3_k(k), u_node(i,j,k),  &
                             v_node(i,j,k), w_node(i,j,k), p_node(i,j,k)
            end do
        end do
    end do

end subroutine writeOutput_NodePoint




subroutine writeOutput_SimData(m1Iter)
    use commondata
    use paramdata
    implicit none
    integer, intent(in) :: m1Iter
    character(len=60) :: pathDataDump
    integer :: n

2   format(3A,1I4.4,1A4)
    write(pathDataDump,2) trim(dataPath), trim(fileName), '_', m1Iter, '.dat'
    
1   format(1A30,1A40)
    write(*,1) 'Saving unformatted data to: ', trim(pathDataDump)
    open(111,file=trim(pathDataDump),form='unformatted')
    write(111) totaltime, p, u1, u2, u3, c, eta, phi,   &
            KA, KC, phiL, conduct,                      &
            (BubbleBlock(n)%xCenter, n=1,nBubbles),     &
            (BubbleBlock(n)%yCenter, n=1,nBubbles),     &
            (BubbleBlock(n)%zCenter, n=1,nBubbles),     &
            (BubbleBlock(n)%radius, n=1,nBubbles),      &
            (BubbleBlock(n)%drdt, n=1,nBubbles),        &
            (BubbleBlock(n)%u, n=1,nBubbles),           &
            (BubbleBlock(n)%v, n=1,nBubbles),           &
            (BubbleBlock(n)%w, n=1,nBubbles),           &
            (BubbleBlock(n)%massFlux, n=1,nBubbles),    &
            (BubbleBlock(n)%hasDeparted, n=1,nBubbles), &
            (BubbleBlock(n)%ForceX, n=1,nBubbles),      &
            (BubbleBlock(n)%ForceY, n=1,nBubbles),      &
            (BubbleBlock(n)%ForceZ, n=1,nBubbles)
    close(111)    
end subroutine writeOutput_SimData