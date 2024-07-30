subroutine initialize()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      - Allocate memory for (global) parameters.
    !!      - Initialize domain geometry and create (uniform) mesh
    !!      - Apply Initial Conditions (IC)
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    real(8) :: c_w

    call AllocArrays

    call initDomain

    call initInitialConditions

    ! Determine explicit time-step
    c_w     =  1d-30 + 10d0*umax
    dtime   = 0.5d0*min( CFL*dx1_i(0)/c_w,              &
                         0.25d0*(dx1_i(0)**2/nuKOH),    &
                         0.25d0*(dx1_i(0)**2/dif(1))    )
    usou2   = 0.05d0*(dx1(0)/dtime)**2
    totaltime   = 0d0
    m1Start     = 1

1   format(1A30)
    write(*,1) 'Initialization complete.'
end subroutine initialize



subroutine AllocArrays()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      Allocate memory for (global) parameters arrays
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: istat
    character(len=60) :: errMsg = 'Insufficient memory during initialization. Stopping program:'

    ! Arrays of rank 1
    allocate(   x1(0:im+1),     x2(0:jm+1),     x3(0:km+1),     &
                x1_i(-1:im+1),  x2_j(-1:jm+1),  x3_k(-1:km+1),  &
                dx1(0:im+1),    dx2(0:jm+1),    dx3(0:km+1),    &
                dx1_i(0:im+1),  dx2_j(0:jm+1),  dx3_k(0:km+1),  &
                c_ref(nspec),   dif(nspec),     c_s(nspec),     &
                massBalance(3), stat = istat                    )

    ! Arrays of rank 2
    allocate(   eta(0:jm+1,0:km+1),     phiL(0:jm+1,0:km+1),    &
                KA(0:jm+1,0:km+1),      KC(0:jm+1,0:km+1),      &
                IBcellLocation(3,1:im*jm*km,nBubbles),          &
                IBfaceXLocation(3,1:im*jm*km/2,nBubbles),       &
                IBfaceYLocation(3,1:im*jm*km/2,nBubbles),       &
                IBfaceZLocation(3,1:im*jm*km/2,nBubbles),       &
                ForceVector(4,nBubbles), stat = istat           )


    ! Array of rank 3
    allocate(   u1(0:im+1,0:jm+1,0:km+1),   u2(0:im+1,0:jm+1,0:km+1),       &
                u3(0:im+1,0:jm+1,0:km+1),   p(0:im+1,0:jm+1,0:km+1),        &
                phi(0:im+1,0:jm+1,0:km+1),  div(1:im,1:jm,1:km),            &
                rhu(0:im,0:jm,0:km),        rhv(0:im,0:jm,0:km),            &
                rhw(0:im,0:jm,0:km),        rhs(0:im,0:jm,0:km),            &
                store_u(0:im+1,0:jm+1,0:km+1),                              &
                store_v(0:im+1,0:jm+1,0:km+1),                              &
                store_w(0:im+1,0:jm+1,0:km+1),                              &
                store_p(0:im+1,0:jm+1,0:km+1),                              &
                store_phi(0:im+1,0:jm+1,0:km+1),                            &
                coef1_phi(1:im,1:jm,1:km),  coef2_phi(1:im,1:jm,1:km),      &
                coef3_phi(1:im,1:jm,1:km),  coef4_phi(1:im,1:jm,1:km),      &
                coef5_phi(1:im,1:jm,1:km),  coef6_phi(1:im,1:jm,1:km),      &
                work1(0:im,0:jm,0:km),      work2(0:im,0:jm,0:km),          &
                work3(0:im,0:jm,0:km),      conduct(0:im,0:jm+1,0:km+1),    &
                CellIsFlagged(0:im+1,0:jm+1,0:km+1),                        &
                IBtype(0:im+1,0:jm+1,0:km+1),                               &
                IBtypeU1(0:im+1,0:jm+1,0:km+1),                             &
                IBtypeU2(0:im+1,0:jm+1,0:km+1),                             &
                IBtypeU3(0:im+1,0:jm+1,0:km+1),                             &
                IBbubble(0:im+1,0:jm+1,0:km+1),                             &
                IBbubbleU1(0:im+1,0:jm+1,0:km+1),                           &
                IBbubbleU2(0:im+1,0:jm+1,0:km+1),                           &
                IBbubbleU3(0:im+1,0:jm+1,0:km+1),                           &
                probeCell(3,1:im*jm*km/2,nBubbles),                         &
                probeU(3,1:im*jm*km/2,nBubbles),                            &
                probeV(3,1:im*jm*km/2,nBubbles),                            &
                probeW(3,1:im*jm*km/2,nBubbles),    stat = istat            )

    ! Array of rank 4
    allocate(   c(0:im+1,0:jm+1,0:km+1,nspec),                              &
                store_c(0:im+1,0:jm+1,0:km+1,nspec),                        &
                ErrIter_c(0:im+1,0:jm+1,0:km+1,nspec),  stat = istat        )

    ! Store (reference) values in an array for easy acces.
    c_ref(1:3)  = (/crefH2, crefKOH, crefH2O/)
    c_s(1)      = c_ref(1)
    dif(1:3)    = (/DifH2,  DifKOH, DifH2O/)

    
    call checkAllocError(istat, errMsg)
    
end subroutine AllocArrays



subroutine initDomain()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      Initialize domain spatial coordinates with
    !!      corresponding cell index
    !!
    !!      *Spatial coordinates:*
    !!          - cell centers : x1(i), x2(j), x3(k)
    !!          - cell faces : x1_i(i), x2_j(j), x3_k(k)
    !!          - cell center distance : dx1(i), dx2(j), dx3(k)
    !!          - cell face distance : dx1_i(i), dx2_j(j), 
    !!                                  dx3_k(k)
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none

    integer :: i,j,k

    ! x-cells
    do i = -1,im+1
        x1_i(i) = dble(i)*lx1/dble(im)       
    end do

    do i = 0,im+1
        x1(i) = 0.5d0*(x1_i(i-1) + x1_i(i))
        dx1(i) = x1_i(i) - x1_i(i-1)        
    end do

    do i = 0,im
        dx1_i(i) = x1(i+1) - x1(i)
    end do

    ! y-cells
    do j = -1,jm+1
        x2_j(j) = dble(j)*lx2/dble(jm)
    end do

    do j = 0,jm+1
        x2(j) = 0.5d0*(x2_j(j-1) + x2_j(j))
        dx2(j) = x2_j(j) - x2_j(j-1)        
    end do

    do j=0,jm
        dx2_j(j) = x2(j+1) - x2(j)
    end do

    ! z-cells
    do k = -1,km+1
        x3_k(k) = dble(k)*lx3/dble(km)
    end do

    do k = 0,km+1
        x3(k) = 0.5d0*(x3_k(k-1) + x3_k(k))
        dx3(k) = x3_k(k) - x3_k(k-1)
    end do
    
    do k = 0,km
        dx3_k(k) = x3(k+1) - x3(k)
    end do

end subroutine initDomain



subroutine initInitialConditions()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      Apply Initial Conditions to the domain.
    !!
    !!      *Scalar values:*
    !!          - potential : phi(i,j,k) = phiR
    !!          - pressure: p(i,j,k) = 0
    !!          - species concentration : c(i,j,k,n) = c_ref(n)
    !!
    !!      *Velocity values:
    !!          - u1(i,j,k) = u3(i,j,k) = 0
    !!          - u2(i,j,k) = parabolic profile
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k,nConc


    ! Potential IC at electrode
    do j = 0,jm+1
        do k = 0,km+1
            eta(j,k)    = 0.0d0
            phiL(j,k)   = phiLe 
        end do
    end do

    ! Domain IC for potential, pressure and flow
!$OMP PARALLEL DO PRIVATE(i,j,k)
    do i = 0,im+1
        do j = 0,jm+1
            do k = 0,km+1
                phi(i,j,k)  = phiR
                p(i,j,k)    = 0.0d0
                u1(i,j,k)   = 0.0d0 
                u2(i,j,k)   = 4.0d0*umax*(x1(i)/lx1)*(1.0d0 - (x1(i))/lx1)
                u3(i,j,k)   = 0.0d0
            end do
        end do
    end do

    ! Species concentration IC over entire domain
    do nConc = 1,nspec
!$OMP PARALLEL DO PRIVATE(i,j,k)
        do i = 0,im+1
            do j = 0,jm+1
                do k = 0,km+1
                    c(i,j,k,nConc) = c_ref(nConc)
                end do
            end do
        end do
    end do
    
end subroutine initInitialConditions




subroutine initLoadPreviousSimulation()
    use commondata
    use paramdata
    implicit none
    character(len=60) :: pathDataDump
    integer :: n

    m1Start = m1CompStart

    write(pathDataDump,'(3A,I4.4,A4)') trim(dataPath), trim(fileName),  &
            '_', m1CompStart,'.dat'

    open(111,file=trim(pathDataDump),form='unformatted')
    read(111) totaltime, p, u1, u2, u3, c, eta, phi,    &
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
            (BubbleBlock(n)%hasDeparted, n=1,nBubbles)
    close(111)
1   format(A30,A31)
    write(*,1) '','Completed loading previous data'

end subroutine initLoadPreviousSimulation