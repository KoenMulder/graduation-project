subroutine boundaryConditions()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      Applies Boundary Conditions for
    !!      - Kinetic prefactors KA and KC
    !!      - Velocity u1, u2, and u3
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata

    call BCkinPref

    call BCu1
    call BCu2
    call BCu3

    write(*,*) 'Boundary Conditions complete.'
    
end subroutine boundaryConditions



subroutine BCu1()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      Applies Boundary Conditions for u-velocity (u1)
    !!      - Dirichlet at the walls (u_wall = 0)
    !!      - Neumann at inflow and outflow (du/dy = 0)
    !!      - Periodic in z-direction (u_front = u_back)
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none

    u1(0,1:jm,1:km)         = 0d0                   ! Electrode
    u1(im,1:jm,1:km)        = 0d0                   ! Membrame

    u1(0:im+1,0,1:km)       = -u1(0:im+1,1,1:km)    ! Bottom (inflow)
    u1(0:im+1,jm+1,1:km)    = u1(0:im+1,jm,1:km)    ! Top (outflow)

    u1(0:im+1,0:jm+1,0)     = u1(0:im+1,0:jm+1,km)  ! Front = back
    u1(0:im+1,0:jm+1,km+1)  = u1(0:im+1,0:jm+1,1)   ! Back = front

end subroutine BCu1



subroutine BCu2()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      Applies Boundary Conditions for v-velocity (u2)
    !!      - Prescribed parabolic flow at inflow
    !!      - Neumann at outflow (dv/dy = 0)
    !!      - Periodic in z-direction (v_front = v_back)
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: i,k

    
    do i = 1,im
        do k = 1,km
            u2(i,0,k) = 4.0d0*umax * (x1(i)/lx1) * (1.0d0 - (x1(i)/lx1)) 
        end do
    end do
    u2(1:im,jm+1,1:km)      = u2(1:im,jm,1:km)      ! Top (Neumann)

    u2(0,0:jm+1,1:km)       = -u2(1,0:jm+1,1:km)    ! Electrode (Zero Neumann) ?
    u2(im+1,0:jm+1,1:km)    = -u2(im,0:jm+1,1:km)   ! Membrame (Zero Neumann) ?


    u2(0:im+1,0:jm+1,0)     = u2(0:im+1,0:jm+1,km)  ! Front = back
    u2(0:im+1,0:jm+1,km+1)  = u2(0:im+1,0:jm+1,1)   ! Back = front

end subroutine BCu2



subroutine BCu3()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      Applies Boundary Conditions for w-velocity (u3)
    !!      - Neumann at outflow (dw/dy = 0)
    !!      - Periodic in z-direction (w_front = w_back)
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    
    u3(0,1:jm,1:km)         = -u3(1,1:jm,1:km)      ! Electrode (Dirichlet)
    u3(im+1,1:jm,1:km)      = -u3(im,1:jm,1:km)     ! Membrame 

    u3(0:im+1,0,1:km)       = -u3(0:im+1,1,1:km)    ! Bottom
    u3(0:im+1,jm+1,1:km)    = u3(0:im+1,jm,1:km)    ! Top

    u3(0:im+1,0:jm+1,0)     = u3(0:im+1,0:jm+1,km)  ! Front = back
    u3(0:im+1,0:jm+1,km+1)  = u3(0:im+1,0:jm+1,1)   ! Back = front
    
end subroutine BCu3




subroutine BCkinPref()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      Applies Boundary Conditions Kinetic prefactors KA
    !!      and KC
    !!
    !!  **Requires:** mathFunctions.f90
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: j,k
    ! real(8) :: kinPrefKA, kinPrefKC

    ! Compute kinetic prefactors
    ! do j = 1,jm
    !     do k = 1,km
    !         KA(j,k) = KinPrefKA(j0, c(0,j,k,2), c(1,j,k,2), c_ref(2),   &
    !                              c(0,j,k,1), c(1,j,k,1), c_ref(1))
    !         KC(j,k) = KinPrefKC(-j0, c(0,j,k,3), c(1,j,k,3), c_ref(3))
    !     end do
    ! end do

    ! Compute kinetic prefactors
    do k=1,km
        do j=1,jm
            ka(j,k) = ka0*(0.5d0*(c(0,j,k,2)+c(1,j,k,2))/c_ref(2))*     &
                    sqrt(abs(0.5d0*(c(0,j,k,1)+c(1,j,k,1)))/c_ref(1)) 
            kc(j,k) = kc0*0.5d0*(c(0,j,k,3)+c(1,j,k,3))/c_ref(3) 
        end do
    end do

    Ka(0,0:km)      = Ka(1,0:km)    ! Electrode Bottom
    Ka(jm+1,0:km)   = Ka(jm,0:km)   ! Electrode Top
    Ka(0:jm+1,0)    = Ka(0:jm+1,1)  ! Electrode Front
    Ka(0:jm+1,km+1) = Ka(0:jm+1,km) ! Electrode Back

    Kc(0,0:km)      = Kc(1,0:km)    ! Electrode Bottom
    Kc(jm+1,0:km)   = Kc(jm,0:km)   ! Electrode Top
    Kc(0:jm+1,0)    = Kc(0:jm+1,1)  ! Electrode Front
    Kc(0:jm+1,km+1) = Kc(0:jm+1,km) ! Electrode Back
    
end subroutine BCkinPref




subroutine BCpressure()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      - Applies (periodic) Boundary Conditions for the 
    !!      pressure in the fluid domain (not bubble)
    !!      - Dirichlet boundary condition
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none

    ! Simple Dirichlet boundary at left and right side
    ! Second order

    p(0,1:jm,1:km)=p(1,1:jm,1:km)
    p(im+1,1:jm,1:km)=p(im,1:jm,1:km) 

    p(0:im+1,0,1:km)=p(0:im+1,1,1:km)
    p(0:im+1,jm+1,1:km)=-p(0:im+1,jm,1:km)

    ! Periodic boundary conditions

    p(0:im+1,0:jm+1,0)=p(0:im+1,0:jm+1,km)
    p(0:im+1,0:jm+1,km+1)=p(0:im+1,0:jm+1,1)
    
end subroutine BCpressure

