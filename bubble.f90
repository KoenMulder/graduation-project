subroutine bubble_Main()
    !! ---------------------------------
    !! Governs the translational kinematics 
    !! of the bubble after departure and
    !! the growing of the bubble before
    !! departure.
    !! ---------------------------------
    use var2
    implicit none
    real(8) ::  xc_old

    if ( bubble ) then
        SELECT CASE (FlowCondPreset)
            CASE ('Stokes')
                ! No growth in Stokes flow

            CASE DEFAULT
                ! Determine bubble growth
                if ( radius.lt.R_max ) then
                    call bulb_flux
                    call growing
        
                    ! Check if bubble is at least 1 cell away from electrode and membrame
                    ! only do this during growth
                    if ( (xc - radius).lt.x1_i(1) ) then
                        xc_old = xc
                        xc = radius + 1d0*x1_i(1)
                        dxcdt = (xc - xc_old)/dtime
                    end if
                end if
        END SELECT

        ! These parts are always calculated
        call domain_flux

        call bulb_force_x
        call domain_flux_force_x

        call bulb_force_y
        call domain_flux_force_y

        call bulb_force_z

        !<<<<<<<<<<< Bubble departure
        if ( radius.ge.R_max) then
            ! Set bubble properties to zero for better post-processing
            bubbleMassflux = 0d0
            drdt = 0d0
            ! Assign detachment label
            if ( isDetached.eqv..false. ) then
                isDetached = .true.
            end if
            
            ! Compute bubble position
            call bubble_COM
        end if
    end if
    
end subroutine bubble_Main




subroutine bubble_COM()
    !!------------------------------------------
    !! Computes the new velocity and position of
    !! the bubble's center of mass.
    !!
    !! Check 
    !! Select a scheme for the virtual force
    !! acceleration term:
    !!  - Explicit single step
    !!  - Adams-Bashforth two-step
    !!------------------------------------------
    use var2
    implicit none
    real(8) :: xc_old, deltaX

    ! If virtual force method is toggled on or off
    if ( VirtForcMethod.eqv..true. ) then

        call bubble_COM_VirtualForceMethod
        
    elseif ( VirtForcMethod.eqv..false.) then

        call bubble_COM_StandardMethod

    else
        error stop 'Invalid VirtualForceMethod setting'
    end if

    ! Check if bubble is at least 1 cell away from electrode and membrame
    if ( (xc - radius).lt.x1_i(1) ) then
        xc_old = xc
        deltaX = x1_i(1) - (xc_old - radius)
        xc = xc_old + deltaX
        dxcdt_old = dxcdt
        dxcdt = dxcdt_old + (xc - xc_old)/dtime
    end if
end subroutine bubble_COM




subroutine bubble_COM_VirtualForceMethod()
    !!------------------------------------------
    !! Computes the new velocity and position of
    !! the bubble's center of mass explicitly 
    !! using the virtual force method.
    !!
    !! Select a scheme for the virtual force
    !! acceleration term specified in the 
    !! module:
    !!  - Explicit single step
    !!  - Adams-Bashforth two-step
    !!------------------------------------------
    use var2
    implicit none
    real(8) :: xc_old, yc_old, zc_old, cu, ca, cg, &
    Vbubble
    
    Vbubble = 4d0/3d0*pi*radius**3d0

    ! Virtual force coefficients with dimensions [m-3], [-]
    ! and [-] respectively:
    cu = rhoKOH/((rhoH2 + Cv*rhoKOH)*Vbubble)
    ca = Cv*rhoKOH/(rhoH2 + Cv*rhoKOH)
    cg = (rhoH2 - rhoKOH)/(rhoH2 + Cv*rhoKOH)

    ! Compute acceleration by bubble surface force: [m/s2]
    ! A negative sign due to change of reference frame
    fsx = cu*(-1d0*FxCOM)
    fsy = cu*(-1d0*FyCOM)
    fsz = cu*(-1d0*FzCOM)
    
    ! Select integration scheme
    select case (timeDiscrScheme)
    CASE('FOBD')
        ! First order explicit time scheme for virtual force acceleration
        ! (du/dt/)_n = (u_n - u_n{n-1})/dt
        ddxcdttFv = (dxcdt - dxcdt_old)/dtime
        ddycdttFv = (dycdt - dycdt_old)/dtime
        ddzcdttFv = (dzcdt - dzcdt_old)/dtime

    CASE('SOBD')
        ! Second order implicit time scheme for virtual force acceleration
        ! (du/dt)_n = (3*u_n - 4*u_{n-1} + n_{n-2})/dt
        ddxcdttFv = (3d0*dxcdt - 4d0*dxcdt_old + dxcdt_oldold)/dtime
        ddycdttFv = (3d0*dycdt - 4d0*dycdt_old + dycdt_oldold)/dtime
        ddzcdttFv = (3d0*dzcdt - 4d0*dzcdt_old + dzcdt_oldold)/dtime

        ! store n-2 velocity
        dxcdt_oldold = dxcdt_old
        dycdt_oldold = dycdt_old
        dzcdt_oldold = dzcdt_old

    CASE DEFAULT
        error stop 'Invalid temporal integration scheme for bubble motion'
    end select

    ! Compute accelertation by virtual force using backwards
    ! finite differencing [m/s2]
    fvx = ca*ddxcdttFv
    fvy = ca*ddycdttFv
    fvz = ca*ddzcdttFv

    ! Compute buoyancy acceleration (only present for bubble)
    ! Direction: downwards= -g,  upwards= +g. Specified in module
    fby = cg*g

    ! Explicit temporal scheme for COM velocity
    ! and store n-1 velocity
    dxcdt_old = dxcdt
    dxcdt = dxcdt_old + dtime*(fsx + fvx)
    dycdt_old = dycdt
    dycdt = dycdt_old + dtime*(fsy + fvy + fby)
    dzcdt_old = dzcdt
    dzcdt = dzcdt_old + dtime*(fsz + fvz)


    ! Explicit temporal scheme for COM position
    xc_old = xc
    xc = xc_old + dtime*dxcdt
    yc_old = yc
    yc = yc_old + dtime*dycdt
    zc_old = zc
    zc = zc_old + dtime*dzcdt
    
end subroutine bubble_COM_VirtualForceMethod




subroutine bubble_COM_StandardMethod()
    !!------------------------------------------
    !! Computes the new velocity and position of
    !! the bubble's center of mass explicitly 
    !! using the standard bubble equation of 
    !! motion.
    !!
    !!------------------------------------------
    use var2
    real(8) :: xc_old, yc_old, zc_old, cu, cg, &
    Vbubble

    Vbubble = 4d0/3d0*pi*radius**3d0 

    ! Virtual force coefficients with dimensions [m-3], [-]
    ! and [-] respectively:
    cu = rhoKOH/(rhoH2*Vbubble)
    cg = (rhoH2 - rhoKOH)/rhoH2

    ! Compute acceleration by bubble surface force: [m/s2]
    ! A negative sign due to change of reference frame
    fsx = cu*(-1d0*FxCOM)
    fsy = cu*(-1d0*FyCOM)
    fsz = cu*(-1d0*FzCOM)

    ! Compute buoyancy acceleration (only present for bubble)
    fby = cg*g

    ! Explicit temporal scheme for COM velocity
    dxcdt_old = dxcdt
    dxcdt = dxcdt_old + dtime*(fsx)
    dycdt_old = dycdt
    dycdt = dycdt_old + dtime*(fsy + fby)
    dzcdt_old = dzcdt
    dzcdt = dzcdt_old + dtime*(fsz)

    ! Explicit temporal scheme for COM position
    xc_old = xc
    xc = xc_old + dtime*dxcdt
    yc_old = yc
    yc = yc_old + dtime*dycdt
    zc_old = zc
    zc = zc_old + dtime*dzcdt

end subroutine bubble_COM_StandardMethod