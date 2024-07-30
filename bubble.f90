subroutine bubbleInitialize()
    !!--------------------------------------------------------------
    !!  **Purpose:**
    !!      Initializes bubbles parameters defined in commondata
    !!      and stores those in the BubbleBlock of type(bubble).
    !!      Data present in the BubbleBlock are:
    !!      - iBubble
    !!      - xCenter
    !!      - yCenter
    !!      - zCenter
    !!      - radius
    !!      - drdt
    !!      - massFlux
    !!      - Velocity: u,v,w
    !!
    !!--------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: n, istat
    character(len=60) :: errMsg = 'Insufficient memory for Bubbles. Stopping program:'
    save

    if ( SimBubble ) then
        allocate(BubbleBlock(nBubbles), stat=istat)

        call checkAllocError(istat, errMsg)

        ! Assign initial values to bubble i
        do n = 1,nBubbles
            BubbleBlock(n)%iBubble  = n
            BubbleBlock(n)%xCenter  = xCenter0 + 4*radius0*(n-1)
            BubbleBlock(n)%yCenter  = yCenter0
            BubbleBlock(n)%zCenter  = zCenter0
            BubbleBlock(n)%radius   = radius0
            BubbleBlock(n)%drdt     = 0.0d0
            BubbleBlock(n)%massFlux = 0.0d0
            BubbleBlock(n)%u        = 0.0d0
            BubbleBlock(n)%v        = 0.0d0
            BubbleBlock(n)%w        = 0.0d0
            BubbleBlock(n)%hasDeparted = .false.
            BubbleBlock(n)%nFacesX  = 0
            BubbleBlock(n)%nFacesY  = 0
            BubbleBlock(n)%nFacesZ  = 0
            BubbleBlock(n)%nIB      = 0
            BubbleBlock(n)%ForceX   = 0.0d0
            BubbleBlock(n)%ForceY   = 0.0d0
            BubbleBlock(N)%ForceZ   = 0.0d0
        end do

        call bubbleCheckDistance
    end if

end subroutine bubbleInitialize




subroutine bubbleMain()
    !!------------------------------------------------------------
    !!  **Main Subroutine for bubble:** 
    !!      - Calls *bubbleFlux* to determine the hydrogen
    !!      diffusion into the bubble.
    !!      - Calls *bubbleGrowth* to determine the new size of
    !!      the bubble with the added mass.
    !!      - Calls *bubbleForces* to compute the hydrodynamic 
    !!      forces acting on the bubble.
    !!
    !!------------------------------------------------------------
    use paramdata
    implicit none

    if ( SimBubble ) then
        ! Compute bubble flux and store it in
        ! bubbleblock(n)%massflux
        call bubbleFlux

        ! Compute bubble growth
        call bubbleGrowth

        ! Verify distsance electrode and bubble surface
        call bubbleCheckDistance

        ! Compute forces acting on bubble
        call bubbleForces
        
    end if


end subroutine BubbleMain




subroutine bubbleFlux()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      - Compute the total flux of hydrogen through the 
    !!      surface into the bubble.
    !!      - Assumes forward staggered grid for velocity
    !!      component
    !!      - Goes through IB-ghost cells (id=-1) and checks
    !!      whether which direction the neighbouring IB-fluid 
    !!      (id=1) is located for the bubble inward normal.
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: i, j, k, n
    real(8) :: term1, term2, massFluxSum
    character(len=20) :: errNaNval, errNaNsub

    write(errNaNval,'(A)') 'massFluxSum'
    write(errNaNsub,'(A)') 'bubbleFlux()' 

    do n = 1,nBubbles
        ! Reset bubble mass flux
        BubbleBlock(n)%massFlux = 0d0
        massFluxSum             = 0d0

        do i = 0,im+1
            do j = 0,jm+1
                do k = 0,km+1
                    if ( (IBtype(i,j,k).eq.IB_label).and.  &
                    (IBbubble(i,j,k).eq.BubbleBlock(n)%iBubble) ) then

                        ! Mass flux in x-direction
                        if ( IBtype(i-1,j,k).eq.fluid_label ) then
                            ! x1_minus: outward normal negative x-direction
                            term1 = -0.5d0*(c(i,j,k,1) + c(i-1,j,k,1))* &
                                    u1(i-1,j,k)*dx2(j)*dx3(k)
                            term2 = dif(1)*(c(i,j,k,1) - c(i-1,j,k,1))/ &
                                    dx1_i(i-1)*dx2(j)*dx3(k)

                            ! if ( term2.gt.0d0 ) then
                            !     write(*,*) 'term2 > 0 in x-minus-direction '
                            ! end if

                            massFluxSum = massFluxSum + term2
                        end if
                        
                        if ( IBtype(i+1,j,k).eq.fluid_label ) then
                            ! x1_plus: outward normal positive x-direction
                            term1 = -(-0.5d0*(c(i+1,j,k,1) + c(i,j,k,1))*   &
                                    u1(i,j,k)*dx2(j)*dx3(k))
                            term2 = -(dif(1)*(c(i+1,j,k,1) - c(i,j,k,1))/   &
                                    dx1_i(i)*dx2(j)*dx3(k))
                            
                            ! if ( term2.gt.0d0 ) then
                            !     write(*,*) 'term2 > 0 in x-plus-direction '
                            !     ! error stop 'term2 > 0 in x-plus-direction '
                            ! end if

                            massFluxSum = massFluxSum + term2
                        end if

                        ! Mass flux in y-direction
                        if ( (IBtype(i,j-1,k).eq.fluid_label) ) then
                            ! x2_minus: outward normal negative y-direction
                            term1 = -0.5d0*(c(i,j,k,1) + c(i,j-1,k,1))* &
                                    u2(i,j-1,k)*dx1(i)*dx3(k)
                            term2 = dif(1)*(c(i,j,k,1) - c(i,j-1,k,1))/ &
                                    dx2_j(j-1)*dx1(i)*dx3(k)

                            ! if ( term2.gt.0d0 ) then
                            !     write(*,*) 'term2 > 0 in y-minus-direction '
                            !     ! error stop 'term2 > 0 in y-minus-direction '
                            ! end if
                            
                            massFluxSum = massFluxSum + term2
                        end if

                        if ( (IBtype(i,j+1,k).eq.fluid_label) ) then
                            ! x2_plus: outward normal positive y-direction
                            term1 = -(-0.5d0*(c(i,j+1,k,1) + c(i,j,k,1))*   &
                                    u2(i,j,k)*dx1(i)*dx3(k))
                            term2 = -(dif(1)*(c(i,j+1,k,1) - c(i,j,k,1))/   &
                                    dx2_j(j)*dx1(i)*dx3(k))

                            ! if ( term2.gt.0d0 ) then
                            !     write(*,*) 'term2 > 0 in y-plus-direction '
                            !     ! error stop 'term2 > 0 in y-plus-direction '
                            ! end if
                            
                            massFluxSum = massFluxSum + term2
                        end if

                        ! Mass flux in z-direction
                        if ( (IBtype(i,j,k-1).eq.fluid_label) ) then
                            ! x3_minus: outward normal negative z-direction
                            term1 = -0.5d0*(c(i,j,k,1) + c(i,j,k-1,1))* &
                                    u3(i,j,k-1)*dx1(i)*dx2(j)
                            term2 = dif(1)*(c(i,j,k,1) - c(i,j,k-1,1))/ &
                                    dx3_k(k-1)*dx1(i)*dx2(j)

                            ! if ( term2.gt.0d0 ) then
                            !     write(*,*) 'term2 > 0 in z-minus-direction '
                            !     ! error stop 'term2 > 0 in z-minus-direction '
                            ! end if

                            massFluxSum = massFluxSum + term2
                            
                        end if

                        if ( (IBtype(i,j,k+1).eq.fluid_label) ) then
                            ! x3_plus: outward normal positive z-direction
                            term1 = -(-0.5d0*(c(i,j,k+1,1) + c(i,j,k,1))*   &
                                    u3(i,j,k)*dx1(i)*dx2(j))
                            term2 = -(dif(1)*(c(i,j,k+1,1) - c(i,j,k,1))/   &
                                    dx3_k(k)*dx1(i)*dx2(j))

                            ! if ( term2.gt.0d0 ) then
                            !     write(*,*) 'term2 > 0 in z-plus-direction '
                            !     ! error stop 'term2 > 0 in z-plus-direction '
                            ! end if

                            massFluxSum = massFluxSum + term2
                        end if
                    end if

                    ! Check for NaN
                    call checkNaNvalue(massFluxSum, errNaNval, errNaNsub)
                end do
            end do
        end do

        BubbleBlock(n)%massFlux = massFluxSum
    end do
end subroutine bubbleFlux




subroutine bubbleGrowth()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      - Compute bubble size via the ideal gas law, Laplace
    !!      pressure, and the total flux of hydrogen through the 
    !!      surface into the bubble.
    !!      - Average pressure over all the IB-ghost cells of a 
    !!      bubble is taken
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: i, j, k, n, countIB
    real(8) :: radiusOld, radiusNew, p_IB_ghost, p_outside,        &
                p_bubble, p_atmospheric, deltaV
    character(len=20) :: errNaNval, errNaNsub
    write(errNaNsub,'(A)') 'bubbleGrowth()' 

    p_atmospheric = 1d5

    do n = 1,nBubbles
        p_IB_ghost = 0d0
        countIB = 0

        do i = 0,im+1
            do j = 0,jm+1
                do k = 0,km+1
                    if ( (IBtype(i,j,k).eq.IB_label).and.  &
                    (IBbubble(i,j,k).eq.BubbleBlock(n)%iBubble) ) then
                        p_IB_ghost  = p_IB_ghost + p(i,j,k)
                        countIB     = countIB + 1
                        
                    end if
                end do
            end do
        end do

        ! Pressure outside the bubble
        p_outside = p_IB_ghost/dble(countIB)

        ! Laplace pressure over the bubble
        p_bubble = p_outside + (2d0/BubbleBlock(n)%radius)*gammaKOH &
                 + p_atmospheric

        ! Ideal gas law
        deltaV = -1d0*BubbleBlock(n)%massFlux*Ru*temp*dtime/p_bubble

        radiusOld = BubbleBlock(n)%radius

        radiusNew = (3d0/4d0/pi*deltaV + radiusOld**3)**(1d0/3d0)

        ! Store values
        BubbleBlock(n)%radius = radiusNew
        BubbleBlock(n)%drdt   = (radiusNew - radiusOld)/dtime

        ! check for NaN values
        write(errNaNval,'(A)') 'p_bubble'
        call checkNaNvalue(radiusNew,errNaNval,errNaNsub)

    end do
end subroutine bubbleGrowth




subroutine bubbleCheckDistance()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      - If the bubble has not departed, check whether there
    !!      is 1 cell distance between bubble surface and the
    !!      electrode. 
    !!      - If the bubble has departed, check whether there is
    !!      1 cell distance between surfaces of other bubbles,
    !!      the electrode and the membrame.
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: n
    real(8) :: xC_old, xC_new, dxcdt, radius, delta

    do n = 1,nBubbles
        xC_old = BubbleBlock(n)%xCenter
        radius = BubbleBlock(n)%radius
        delta = 0

        ! Check electrode distance
        if( (xC_old - radius).lt.x1_i(1)) then
            delta  = x1_i(1) - (xC_old - radius)
            xC_new = xC_old + delta
            dxcdt  = (xC_new - xC_old)/dtime

            ! store values
            BubbleBlock(n)%xCenter = xC_new
            BubbleBlock(n)%u       = dxcdt
        else
            dxcdt            = 0.0d0
            BubbleBlock(n)%u = dxcdt
        end if

        ! Check membrame distance

        ! Check distance between bubbles
    end do
    
end subroutine bubbleCheckDistance





subroutine bubbleForces()
    use commondata
    use paramdata
    implicit none
    integer :: id,jd,kd, n, nFace
    real(8) :: term1, term2, term3, SumForceX, SumForceY,   &
                SumForceYbellow, SumForceYabove

    SumForceX       = 0.0d0
    SumForceY       = 0.0d0
    SumForceYbellow = 0.0d0
    SumForceYabove  = 0.0d0
    do n = 1,nBubbles
        ! Check forces in x-direction
        do nFace = 1,BubbleBlock(n)%nFacesX
            id = IBfaceXLocation(1,nFace,n)
            jd = IBfaceXLocation(2,nFace,n)
            kd = IBfaceXLocation(3,nFace,n)
            
            if ( IBtypeU1(id+1,jd,kd).eq.0 ) then
                term1 = -(-(0.5d0*(u1(id+1,jd,kd) + u1(id,jd,kd)))**2*&
                        dx2(jd)*dx3(kd))
                term2 = -(-p(id+1,jd,kd)/rhoKOH*dx2(jd)*dx3(kd))
                term3 = -(nuKOH*(u1(id+1,jd,kd) - u1(id,jd,kd))/    &
                        dx1(id+1)*dx2(jd)*dx3(kd))  
                SumForceX = SumForceX + term1 + term2 + term3

                ! ???
                ! if (x2_j(j).gt.yc) then
                !     total_up=total_up+term1+term2+term3
                !  else
                !     total_bellow=total_bellow+term1+term2+term3
                !  end if
            end if

            if ( IBtypeU1(id-1,jd,kd).eq.0 ) then
                term1 = -(0.5d0*(u1(id-1,jd,kd) + u1(id,jd,kd)))**2*&
                        dx2(jd)*dx3(kd)
                term2 = -p(id,jd,kd)/rhoKOH*dx2(jd)*dx3(kd)
                term3 = nuKOH*(u1(id,jd,kd) - u1(id-1,jd,kd))/      &
                        dx1(id)*dx2(jd)*dx3(kd)
                SumForceX = SumForceX + term1 + term2 + term3

                ! ???
                ! if (x2_j(j).gt.yc) then
                !     total_up=total_up+term1+term2+term3
                !  else
                !     total_bellow=total_bellow+term1+term2+term3
                !  end if
            end if

            if ( IBtypeU1(id,jd+1,kd).eq.0 ) then
                term1 = -(-0.5d0*(u1(id,jd,kd) + u1(id,jd+1,kd))*   &
                           0.5d0*(u2(id,jd,kd) + u2(id+1,jd,kd))*   &
                           dx1(id)*dx3(kd))
                term2 = 0d0 
                term3 = -(nuKOH*(u1(id,jd+1,kd) - u1(id,jd,kd))/    &
                        dx2_j(jd)*dx1(id)*dx3(kd))
                SumForceX = SumForceX + term1 + term2 + term3

                ! ???
                ! if (x2_j(j).gt.yc) then
                !     total_up=total_up+term1+term2+term3
                !  else
                !     total_bellow=total_bellow+term1+term2+term3
                !  end if
            end if

            if ( IBtypeU1(id,jd-1,kd).eq.0 ) then
                term1 = -0.5d0*(u1(id,jd,kd) + u1(id,jd-1,kd))*     &
                         0.5d0*(u2(id,jd-1,kd) + u2(id+1,jd-1,kd))* &
                         dx1(id)*dx3(kd)
                term2 = 0d0
                term3 = nuKOH*(u1(id,jd,kd) - u1(id,jd-1,kd))/      &
                        dx2_j(jd-1)*dx1(id)*dx3(kd)
                SumForceX = SumForceX + term1 + term2 + term3
                
                !???
                ! if (x2_j(j).gt.yc) then
                !     total_up=total_up+term1+term2+term3
                ! else
                !     total_bellow=total_bellow+term1+term2+term3
                ! end if
            end if

            if ( IBtypeU1(id,jd,kd+1).eq.0 ) then
                term1 = -(-0.5d0*(u1(id,jd,kd) + u1(id,jd,kd+1))*   &
                           0.5d0*(u3(id,jd,kd) + u3(id+1,jd,kd))*   &
                           dx1(id)*dx2(jd))
                term2 = 0d0    
                term3 = -(nuKOH*(u1(id,jd,kd+1) - u1(id,jd,kd))/    &
                        dx3_k(kd)*dx1(id)*dx2(jd))  
                SumForceX = SumForceX + term1 + term2 + term3

                ! ???
                ! if (x2_j(j).gt.yc) then
                !     total_up=total_up+term1+term2+term3
                ! else
                !     total_bellow=total_bellow+term1+term2+term3
                ! end if
            end if

            if ( IBtypeU1(id,jd,kd-1).eq.0 ) then
                term1 = -0.5d0*(u1(id,jd,kd) + u1(id,jd,kd-1))*     &
                         0.5d0*(u3(id,jd,kd-1) + u3(id+1,jd,kd-1))* &
                         dx1(id)*dx2(jd)
                term2 = 0d0
                term3 = nuKOH*(u1(id,jd,kd) - u1(id,jd,kd-1))/      &
                        dx3_k(kd-1)*dx1(id)*dx2(jd)
                SumForceX = SumForceX + term1 + term2 + term3

                !???
                ! if (x2_j(j).gt.yc) then
                !     total_up=total_up+term1+term2+term3
                !  else
                !     total_bellow=total_bellow+term1+term2+term3
                !  end if
            end if
        end do

        ! Check forces in y-direction
        do nFace = 1,BubbleBlock(n)%nFacesY
            id = IBfaceXLocation(1,nFace,n)
            jd = IBfaceXLocation(2,nFace,n)
            kd = IBfaceXLocation(3,nFace,n)
            
            if ( IBtypeU2(id+1,jd,kd).eq.0 ) then
                term1 = -(-0.5d0*(u2(id,jd,kd) + u2(id+1,jd,kd))*   &
                           0.5d0*(u1(id,jd,kd) + u1(id,jd+1,kd))*   &
                           dx2(jd)*dx3(kd))
                term2 = 0d0
                term3 = -(nuKOH*(u2(id+1,jd,kd) - u2(id,jd,kd))/    &
                        dx1_i(id)*dx2(jd)*dx3(kd)) 
                SumForceY = SumForceY + term1 + term2 +term3

                if (x2_j(jd).gt.BubbleBlock(n)%yCenter) then
                    SumForceYabove = SumForceYabove + term1 + term2 + term3
                else
                    SumForceYbellow = SumForceYbellow + term1 + term2 + term3
                end if
            end if

            if ( IBtypeU2(id-1,jd,kd).eq.0 ) then
                term1 = -0.5d0*(u2(id,jd,kd) + u2(id-1,jd,kd))*     &
                        0.5d0*(u1(id-1,jd,kd) + u1(id-1,jd+1,kd))*  &
                        dx2(jd)*dx3(kd)
                term2 = 0d0
                term3 = nuKOH*(u2(id,jd,kd) - u2(id-1,jd,kd))/      &
                        dx1_i(id-1)*dx2(jd)*dx3(kd)
                SumForceY = SumForceY + term1 + term2 + term3

                if (x2_j(jd).gt.BubbleBlock(n)%nFacesY) then
                    SumForceYabove = SumForceYabove + term1 + term2 + term3
                else
                    SumForceYbellow = SumForceYbellow + term1 + term2 + term3
                end if
            end if

            if ( IBtypeU2(id,jd+1,kd).eq.0 ) then
                term1 = -(-(0.5d0*(u2(id,jd,kd) + u2(id,jd+1,kd)))**2*  &
                        dx1(id)*dx3(kd))
                term2 = -(-p(id,jd+1,kd)/rhoKOH*dx1(id)*dx3(kd))
                term3 = -(nuKOH*(u2(id,jd+1,kd) - u2(id,jd,kd))/        &
                        dx2_j(jd)*dx1(id)*dx3(kd))
                SumForceY = SumForceY + term1 + term2 + term3

                if (x2_j(jd).gt.BubbleBlock(n)%yCenter) then
                    SumForceYabove = SumForceYabove + term1 + term2 + term3
                else
                    SumForceYbellow = SumForceYbellow + term1 + term2 + term3
                end if
            end if

            if ( IBtypeU2(id,jd-1,kd).eq.0 ) then
                term1 = -(0.5d0*(u2(id,jd,kd) + u2(id,jd-1,kd)))**2*    &
                        dx1(id)*dx3(kd)
                term2 = -p(id,jd,kd)/rhoKOH*dx1(id)*dx3(kd)
                term3 = nuKOH*(u2(id,jd,kd) - u2(id,jd-1,kd))/          &
                        dx2_j(jd-1)*dx1(id)*dx3(kd) 
                SumForceY = SumForceY + term1 + term2 + term3

                if (x2_j(jd).gt.BubbleBlock(n)%yCenter) then
                    SumForceYabove = SumForceYabove + term1 + term2 + term3
                else
                    SumForceYbellow = SumForceYbellow + term1 + term2 + term3
                end if
            end if

            if ( IBtypeU2(id,jd,kd+1).eq.0 ) then
                term1 = -(-0.5d0*(u2(id,jd,kd) + u2(id,jd,kd+1))*   &
                           0.5d0*(u3(id,jd,kd) + u3(id,jd+1,kd))*   &
                           dx1(id)*dx2(jd))
                term2 = 0d0
                term3 = -(nuKOH*(u2(id,jd,kd+1) - u2(id,jd,kd))/    &
                        dx3_k(kd)*dx1(id)*dx2(jd))
                SumForceY = SumForceY + term1 + term2 + term3

                if (x2_j(jd).gt.BubbleBlock(n)%yCenter) then
                    SumForceYabove = SumForceYabove + term1 + term2 + term3
                else
                    SumForceYbellow = SumForceYbellow + term1 + term2 + term3
                end if
            end if

            if ( IBtypeU2(id,jd,kd-1).eq.0 ) then
                term1 = -0.5d0*(u2(id,jd,kd) + u2(id,jd,kd-1))*     &
                         0.5d0*(u3(id,jd,kd-1) + u3(id,jd+1,kd-1))* &
                         dx1(id)*dx2(jd)
                term2 = 0d0
                term3 = nuKOH*(u2(id,jd,kd) - u2(id,jd,kd-1))/      &
                        dx3_k(kd-1)*dx1(id)*dx2(jd)
                SumForceY = SumForceY + term1 + term2 + term3

                if (x2_j(jd).gt.BubbleBlock(n)%yCenter) then
                    SumForceYabove = SumForceYabove + term1 + term2 + term3
                else
                    SumForceYbellow = SumForceYbellow + term1 + term2 + term3
                end if
            end if
        end do

        ! Check forces in z-direction
        ! do nFace = 1, BubbleBlock(n)%nFacesZ
        ! end do

        ForceVector(1,n) = SumForceX
        ForceVector(2,n) = SumForceY
        ForceVector(3,n) = SumForceYabove
        ForceVector(4,n) = SumForceYbellow

        BubbleBlock(n)%ForceX = SumForceX
        BubbleBlock(n)%ForceY = SumForceY
    end do

    
end subroutine bubbleForces
