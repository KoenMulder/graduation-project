subroutine bubbleDetachment()
    !!------------------------------------------------------------
    !!  **Main Subroutine for bubble detachment:** 
    !!
    !!      Implementation of the bubble detachment criterion.
    !!      Options are:
    !!      - *Critical Diameter* defined by means of experimental
    !!      data.
    !!      - *Force Balance* with a resulting force in the
    !!      x-direction larger than a threshold value.
    !!
    !!------------------------------------------------------------
    use paramdata
    use commondata
    implicit none
    integer :: n
    real(8) :: uVelocity, vVelocity, wVelocity

    if ( SimBubble ) then

        ! Simulate detachment?
        if ( SimDetach ) then        
            detachCase: SELECT CASE (detachModel)

            case DEFAULT
                write(*,*) 'Error: undefined model, detachModel=', detachModel
                write(*,*) '--> Check detachModel in module,f90. Stopping simulation'
                stop

            case (1)
                ! Critical Bubble
                do n = 1, nBubbles
                    if ( BubbleBlock(n)%radius.ge.radiusCrit ) then
                        ! Explicit time step:
                        uVelocity = BubbleBlock(n)%ForceX*dtime
                        vVelocity = BubbleBlock(n)%ForceY*dtime
                        
                    end if
                end do

            case (2)
                ! Force balance

            end SELECT detachCase
        end if
    end if
   
end subroutine bubbleDetachment



