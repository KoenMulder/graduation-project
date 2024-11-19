subroutine u()
    use var2
    implicit none
    real(8) ::  tp, ua, va, wa

    if ( rungIter.eq.1 ) then
        storage_u = u1
    end if
    

!$omp parallel do private(ua,va,wa,tp,i,j,k)
    do k=1,km
        do j=1,jm
            do i=1,im!-1 <<--- aparantly this is needed? If code breaks check this
                ua = 0.5d0*(u1(i,j,k) + u1(i-1,j,k))
                tp = -ua**2d0/dx1_i(i)
                ua = 0.5d0*(u1(i+1,j,k) + u1(i,j,k))
                tp= tp + ua**2d0/dx1_i(i)
    
                va = 0.5d0*(u2(i+1,j-1,k) + u2(i,j-1,k))
                tp = tp - 0.5d0*(u1(i,j-1,k) + u1(i,j,k))*va/dx2(j)
                va = 0.5d0*(u2(i+1,j,k) + u2(i,j,k))
                tp = tp + 0.5d0*(u1(i,j,k) + u1(i,j+1,k))*va/dx2(j)
    
                wa = 0.5d0*(u3(i+1,j,k-1) + u3(i,j,k-1))
                tp = tp - 0.5d0*(u1(i,j,k-1) + u1(i,j,k))*wa/dx3(k)
                wa = 0.5d0*(u3(i+1,j,k) + u3(i,j,k))
                tp = tp + 0.5d0*(u1(i,j,k) + u1(i,j,k+1))*wa/dx3(k)
    
                rhu(i,j,k) = storage_u(i,j,k) - tp*dtime*RK_alpha(rungIter)    &
                            - nuKOH*dtime/dx1_i(i)*RK_alpha(rungIter)*         &
                            ((u1(i,j,k) - u1(i-1,j,k))/dx1(i) - (u1(i+1,j,k) - u1(i,j,k))/dx1(i+1))     &
                            - nuKOH*dtime/dx2(j)*RK_alpha(rungIter)*           &
                            ((u1(i,j,k) - u1(i,j-1,k))/dx2_j(j-1) - (u1(i,j+1,k) - u1(i,j,k))/dx2_j(j)) &
                            - nuKOH*dtime/dx3(k)*RK_alpha(rungIter)*           &
                            ((u1(i,j,k) - u1(i,j,k-1))/dx3_k(k-1) - (u1(i,j,k+1) - u1(i,j,k))/dx3_k(k)) 
            end do
        end do
    end do  
end subroutine u



subroutine bou_rhu()
    use var2
    implicit none
    
    do k=1,km
        do j=1,jm
            do i=1,im
                if (typ1(i,j,k).eq.1) then
                    u1(i,j,k) = rhu(i,j,k)
                end if
            end do
        end do
    end do

    SELECT CASE (FlowCondPreset)
        CASE ('Stokes')
            ! periodidc, do nothing
        CASE ('free-slip')
            ! No mass transport at walls
            do k=1,km
                do j=1,jm
                    u1(0,j,k)  = 0d0 ! Electrode
                    u1(im,j,k) = 0d0 ! Membrame
                end do
            end do
        CASE DEFAULT
            ! No mass transport at walls
            do k=1,km
                do j=1,jm
                    u1(0,j,k)  = 0d0 ! Electrode
                    u1(im,j,k) = 0d0 ! Membrame
                end do
            end do
    END SELECT

end subroutine bou_rhu


subroutine bou_u()
    !!----------------------------------------------
    !! Boundary conditions for u-velocity on domain.
    !! Boundaries applied to:
    !! - Outlet = dirichlet
    !! - Inlet  = dirichlet
    !! - Electrode  = dirichlet
    !! - Membrame = dirichlet
    !!---------------------------------------------
    use var2
    implicit none

    SELECT CASE (FlowCondPreset)
        CASE('Stokes')
            ! Periodic in width (x)
            u1(0,   0:jm+1,0:km+1) = u1(im,0:jm+1,0:km+1)
            u1(im+1,0:jm+1,0:km+1) = u1(1, 0:jm+1,0:km+1)

            ! Periodic in height (y)
            u1(0:im+1,0,   0:km+1) = u1(0:im+1,jm,0:km+1)
            u1(0:im+1,jm+1,0:km+1) = u1(0:im+1,1, 0:km+1)

            ! Periodic in depth (z)
            u1(0:im+1,0:jm+1,0)    = u1(0:im+1,0:jm+1,km)
            u1(0:im+1,0:jm+1,km+1) = u1(0:im+1,0:jm+1,1)

        CASE ('free-slip')
            ! No flow normal to wall surface
            u1(0,1:jm,1:km)      = 0d0  ! Electrode
            u1(im,1:jm,1:km)     = 0d0  ! Membrame

            ! Nonzero gradient? at inlet
            u1(0:im+1,0,1:km)    = -u1(0:im+1,1,1:km)   ! Inlet
            ! Zero gradient at oulet
            u1(0:im+1,jm+1,1:km) = u1(0:im+1,jm,1:km)   ! Outlet

            ! Periodic in depth (z)
            u1(0:im+1,0:jm+1,0)    = u1(0:im+1,0:jm+1,km)
            u1(0:im+1,0:jm+1,km+1) = u1(0:im+1,0:jm+1,1)

        CASE DEFAULT
            ! Simple Dirichlet boundary at left and right side
            ! Second order
            ! No flow normal to wall surface
            u1(0,1:jm,1:km)      = 0d0                  ! Electrode
            u1(im,1:jm,1:km)     = 0d0                  ! Membrame

            ! Nonzero gradient? at inlet
            u1(0:im+1,0,1:km)    = -u1(0:im+1,1,1:km)   ! Inlet
            ! Zero gradient at oulet
            u1(0:im+1,jm+1,1:km) = u1(0:im+1,jm,1:km)   ! Outlet

            ! Periodic in depth
            u1(0:im+1,0:jm+1,0)    = u1(0:im+1,0:jm+1,km)
            u1(0:im+1,0:jm+1,km+1) = u1(0:im+1,0:jm+1,1)
    END SELECT
end subroutine bou_u



subroutine v()
    use var2
    implicit none
    real(8) ::  tp,ua,va,wa

    if ( rungIter.eq.1 ) then
        storage_v = u2
    end if

!$omp parallel do private(ua,va,wa,tp,i,j,k)
    do k=1,km
        do j=1,jm
            do i=1,im
                ua = 0.5d0*(u1(i-1,j,k) + u1(i-1,j+1,k))
                tp = -0.5d0*(u2(i-1,j,k) + u2(i,j,k))*ua/dx1(i)
                ua = 0.5d0*(u1(i,j,k) + u1(i,j+1,k))
                tp = tp + 0.5d0*(u2(i,j,k) + u2(i+1,j,k))*ua/dx1(i)

                va = 0.5d0*(u2(i,j,k) + u2(i,j-1,k))
                tp = tp - va**2d0/dx2_j(j)
                va = 0.5d0*(u2(i,j,k) + u2(i,j+1,k))
                tp = tp + va**2d0/dx2_j(j)

                wa = 0.5d0*(u3(i,j,k-1) + u3(i,j+1,k-1))
                tp = tp - 0.5d0*(u2(i,j,k-1) + u2(i,j,k))*wa/dx3(k)
                wa = 0.5d0*(u3(i,j,k) + u3(i,j+1,k))
                tp = tp + 0.5d0*(u2(i,j,k) + u2(i,j,k+1))*wa/dx3(k)

                rhv(i,j,k) = storage_v(i,j,k) - tp*dtime*RK_alpha(rungIter)    &
                        - nuKOH*dtime/dx1(i)*RK_alpha(rungIter)*               &
                        ((u2(i,j,k) - u2(i-1,j,k))/dx1_i(i-1) - (u2(i+1,j,k) - u2(i,j,k))/dx1_i(i)) &
                        - nuKOH*dtime/dx2_j(j)*RK_alpha(rungIter)*             &
                        ((u2(i,j,k) - u2(i,j-1,k))/dx2(j) - (u2(i,j+1,k) - u2(i,j,k))/dx2(j+1))     &
                        - nuKOH*dtime/dx3(k)*RK_alpha(rungIter)*               &
                        ((u2(i,j,k) - u2(i,j,k-1))/dx3_k(k-1) - (u2(i,j,k+1) - u2(i,j,k))/dx3_k(k))
            end do
        end do
    end do
end subroutine v



subroutine bou_rhv()
    use var2
    implicit none
    do k=1,km
        do j=1,jm
            do i=1,im
                if (typ2(i,j,k).eq.1) then
                    u2(i,j,k) = rhv(i,j,k)
                end if
            end do
        end do
    end do

    ! Apply quiescent flow when StokesFlow is selected
    SELECT CASE (FlowCondPreset)
        CASE ('Stokes','free-slip')
            ! Do nothing
            ! u2(1:im,0,1:km) = 0d0 <-- this preset is wrong
        CASE DEFAULT
            do k=1,km
                do i=1,im
                    u2(i,0,k) = 4d0*umax*(x1(i)/lx1)*(1d0 - (x1(i)/lx1))
                end do
            end do
    END SELECT
end subroutine bou_rhv



subroutine bou_v()
    use var2
    implicit none

    SELECT CASE (FlowCondPreset)
        CASE ('Stokes')
            ! Periodic in width (x)
            u2(0,   0:jm+1,0:km+1) = u2(im,0:jm+1,0:km+1)
            u2(im+1,0:jm+1,0:km+1) = u2(1, 0:jm+1,0:km+1)
    
            ! Periodic in height (y)
            u2(0:im+1,0,   0:km+1) = u2(0:im+1,jm,0:km+1)
            u2(0:im+1,jm+1,0:km+1) = u2(0:im+1,1, 0:km+1)
    
            ! Periodic in depth (z)
            u2(0:im+1,0:jm+1,0   ) = u2(0:im+1,0:jm+1,km)
            u2(0:im+1,0:jm+1,km+1) = u2(0:im+1,0:jm+1,1 )

        CASE ('free-slip')
            ! Quiescient flow
            u2(1:im,0,1:km) = 0d0
            u2(1:im,jm+1,1:km) = u2(1:im,jm,1:km)
    
            ! Free slip at wall (zero gradient) 
            u2(0,0:jm+1,1:km)    = u2(1,0:jm+1,1:km)    ! Electrode
            u2(im+1,0:jm+1,1:km) = u2(im,0:jm+1,1:km)   ! Memrame
    
            ! Periodic boundary conditions
            u2(0:im+1,0:jm+1,0)    = u2(0:im+1,0:jm+1,km)
            u2(0:im+1,0:jm+1,km+1) = u2(0:im+1,0:jm+1,1)
        CASE DEFAULT
            ! Prescribed inlet flow velocity
            do k=1,km
                do i=1,im
                    u2(i,0,k) = 4d0*umax*(x1(i)/lx1)*(1d0 - (x1(i)/lx1))
                end do
            end do
            u2(1:im,jm+1,1:km) = u2(1:im,jm,1:km)   ! Outlet

            ! Dirichlet boundary conditions
            u2(0,0:jm+1,1:km)    = -u2(1,0:jm+1,1:km)   ! Electrode
            u2(im+1,0:jm+1,1:km) = -u2(im,0:jm+1,1:km)  ! Memrame

            ! Periodic boundary conditions
            u2(0:im+1,0:jm+1,0)    = u2(0:im+1,0:jm+1,km)
            u2(0:im+1,0:jm+1,km+1) = u2(0:im+1,0:jm+1,1)
    END SELECT
end subroutine bou_v



subroutine w()
    use var2
    implicit none
    real(8) ::  tp,ua,va,wa

    if ( rungIter.eq.1 ) then
        storage_w = u3
    end if
    
!$omp parallel do private(ua,va,wa,tp,i,j,k)
     do k=1,km
        do j=1,jm
            do i=1,im
                ua =  0.5d0*(u1(i-1,j,k) + u1(i-1,j,k+1))
                tp = -0.5d0*(u3(i-1,j,k) + u3(i,j,k))*ua/dx1(i)
                ua =  0.5d0*(u1(i,j,k)   + u1(i,j,k+1))
                tp = tp + 0.5d0*(u3(i,j,k) + u3(i+1,j,k))*ua/dx1(i)

                va = 0.5d0*(u2(i,j-1,k) + u2(i,j-1,k+1))
                tp = tp - 0.5d0*(u3(i,j-1,k) + u3(i,j,k))*va/dx2(j)
                va = 0.5d0*(u2(i,j,k) + u2(i,j,k+1))
                tp = tp + 0.5d0*(u3(i,j,k) + u3(i,j+1,k))*va/dx2(j)

                wa = 0.5d0*(u3(i,j,k) + u3(i,j,k-1))
                tp = tp - wa**2d0/dx3_k(k)
                wa = 0.5d0*(u3(i,j,k+1)+u3(i,j,k))
                tp = tp + wa**2d0/dx3_k(k)

                rhw(i,j,k) = storage_w(i,j,k) - tp*dtime*RK_alpha(rungIter)    &
                            - nuKOH*dtime/dx1(i)*RK_alpha(rungIter)*           &
                            ((u3(i,j,k) - u3(i-1,j,k))/dx1_i(i-1) - (u3(i+1,j,k) - u3(i,j,k))/dx1_i(i)) &
                            - nuKOH*dtime/dx2(j)*RK_alpha(rungIter)*           &
                            ((u3(i,j,k) - u3(i,j-1,k))/dx2_j(j-1) - (u3(i,j+1,k) - u3(i,j,k))/dx2_j(j)) &
                            - nuKOH*dtime/dx3(k)*RK_alpha(rungIter)*           &
                            ((u3(i,j,k) - u3(i,j,k-1))/dx3(k) - (u3(i,j,k+1) - u3(i,j,k))/dx3(k+1))
            end do
        end do
     end do
  
end subroutine w



subroutine bou_rhw()
    use var2
    implicit none
     
    do k=1,km
        do j=1,jm
            do i=1,im
                if(typ3(i,j,k).eq.1) then
                    u3(i,j,k) = rhw(i,j,k)
                end if
            end do
        end do
    end do
  
end subroutine bou_rhw



subroutine bou_w()
    use var2
    implicit none

    SELECT CASE (FlowCondPreset)
        CASE ('Stokes')
            ! Periodic in width
            u3(0,   0:jm+1,0:km+1) = u3(im,0:jm+1,0:km+1)
            u3(im+1,0:jm+1,0:km+1) = u3(1, 0:jm+1,0:km+1) 

            ! Periodic in height
            u3(0:im+1,0,   0:km+1) = u3(0:im+1,jm,0:km+1) 
            u3(0:im+1,jm+1,0:km+1) = u3(0:im+1,1, 0:km+1) 

            ! Periodic in depth
            u3(0:im+1,0:jm+1,0)    = u3(0:im+1,0:jm+1,km)
            u3(0:im+1,0:jm+1,km+1) = u3(0:im+1,0:jm+1,1)

        CASE ('free-slip')
            ! Free slip at wall (zero gradient)
            u3(0,1:jm,1:km)      = u3(1,1:jm,1:km)     ! Electrode
            u3(im+1,1:jm,1:km)   = u3(im,1:jm,1:km)    ! Memrame

            ! No swirl at inlet
            u3(0:im+1,0,1:km)    = -u3(0:im+1,1,1:km)   ! Inlet

            ! Zero Neumann at outlet (zero gradient)
            u3(0:im+1,jm+1,1:km) = u3(0:im+1,jm,1:km)   ! Outlet

            ! Periodic boundary conditions
            u3(0:im+1,0:jm+1,0)    = u3(0:im+1,0:jm+1,km)
            u3(0:im+1,0:jm+1,km+1) = u3(0:im+1,0:jm+1,1)

        CASE DEFAULT
            ! Simple Dirichlet boundary at left and right side
            ! Second order
            ! Dirichlet conditions
            u3(0,1:jm,1:km)      = -u3(1,1:jm,1:km)     ! Electrode
            u3(im+1,1:jm,1:km)   = -u3(im,1:jm,1:km)    ! Memrame

            u3(0:im+1,0,1:km)    = -u3(0:im+1,1,1:km)   ! Inlet

            ! Neumann conditions
            u3(0:im+1,jm+1,1:km) = u3(0:im+1,jm,1:km)   ! Outlet 

            ! Periodic boundary conditions
            u3(0:im+1,0:jm+1,0)    = u3(0:im+1,0:jm+1,km)
            u3(0:im+1,0:jm+1,km+1) = u3(0:im+1,0:jm+1,1)
    END SELECT
end subroutine bou_w



subroutine velocity()
    use var2
    implicit none   

!$omp parallel do private(i,j,k)
    do k=1,km
        do j=1,jm
            do i=1,im
                if(typ1(i,j,k).eq.1) then
                    u1(i,j,k) = u1(i,j,k) - dtime/rhoKOH/dx1_i(i)*(p(i+1,j,k) - p(i,j,k))*RK_alpha(rungIter)
                end if

                if(typ2(i,j,k).eq.1) then
                    u2(i,j,k) = u2(i,j,k) - dtime/rhoKOH/dx2_j(j)*(p(i,j+1,k) - p(i,j,k))*RK_alpha(rungIter)
                end if

                if(typ3(i,j,k).eq.1) then
                    u3(i,j,k) = u3(i,j,k) - dtime/rhoKOH/dx3_k(k)*(p(i,j,k+1) - p(i,j,k))*RK_alpha(rungIter)
                end if
            end do
        end do
    end do

end subroutine velocity



subroutine p_solve_ArtComp()
    use var2
    implicit none

!omp parallel do
    do k=1,km
        do j=1,jm
            do i=1,im
                div(i,j,k) = (u1(i,j,k) - u1(i-1,j,k))/dx1(i) + &
                             (u2(i,j,k) - u2(i,j-1,k))/dx2(j) + &
                             (u3(i,j,k) - u3(i,j,k-1))/dx3(k) 
            end do
        end do
    end do

    if ( rungIter.eq.1 ) then
        storage_p = p
    end if
    
!omp parallel do
    do k=1,km
        do j=1,jm
            do i=1,im
                if (typ(i,j,k).eq.1) then
                    p(i,j,k) = storage_p(i,j,k) - rhoKOH*usou2*dtime*div(i,j,k)*RK_alpha(rungIter)
                end if
            end do
        end do
    end do
    
end subroutine p_solve_ArtComp



subroutine bou_p()
    use var2
    implicit none

    SELECT CASE (FlowCondPreset)
    CASE ('Stokes')
        ! Periodic in height (y)
        p(0,   0:jm+1,0:km+1) = p(im,0:jm+1,0:km+1)
        p(im+1,0:jm+1,0:km+1) = p(1, 0:jm+1,0:km+1)

        ! Periodic in height (y)
        p(0:im+1,0,   0:km+1) = p(0:im+1,jm,0:km+1)
        p(0:im+1,jm+1,0:km+1) = p(0:im+1,1, 0:km+1)

        ! Periodic in depth (z)
        p(0:im+1,0:jm+1,0   ) = p(0:im+1,0:jm+1,km)
        p(0:im+1,0:jm+1,km+1) = p(0:im+1,0:jm+1,1 )

    CASE ('free-slip')
        ! Simple Dirichlet boundary at left and right side
        ! Second order
        p(0,1:jm,1:km)      =  p(1,1:jm,1:km)
        p(im+1,1:jm,1:km)   =  p(im,1:jm,1:km) 
        p(0:im+1,0,1:km)    =  p(0:im+1,1,1:km)
        p(0:im+1,jm+1,1:km) = -p(0:im+1,jm,1:km)

        ! Periodic boundary conditions
        p(0:im+1,0:jm+1,0)    = p(0:im+1,0:jm+1,km)
        p(0:im+1,0:jm+1,km+1) = p(0:im+1,0:jm+1,1)

    CASE DEFAULT
        ! Simple Dirichlet boundary at left and right side
        ! Second order
        p(0,1:jm,1:km)      =  p(1,1:jm,1:km)
        p(im+1,1:jm,1:km)   =  p(im,1:jm,1:km)

        ! 
        p(0:im+1,0,1:km)    =  p(0:im+1,1,1:km)
        p(0:im+1,jm+1,1:km) = -p(0:im+1,jm,1:km)

        ! Periodic boundary conditions
        p(0:im+1,0:jm+1,0)    = p(0:im+1,0:jm+1,km)
        p(0:im+1,0:jm+1,km+1) = p(0:im+1,0:jm+1,1)
    END SELECT
end subroutine bou_p





subroutine bou_phi()
    use var2
    implicit none   
    real(8) :: az,bz,cz,fac,kka,kkc,kkae,kkce
    ! Error_phi and phi_old are now also allocated in memory of init_new
    ! real(8) :: Error_phi(0:im+1,0:jm+1,0:km+1), &
    !        phi_old(0:im+1,0:jm+1,0:km+1)

! c     phi1 = phi(i=1)
! c     phi0 = phi(i=0)
! c     electrode at x=0: -condL*(phi1-phi0)/dx=ka*exp(aa*eta)+kc*exp(-ac*eta) 
! c     condL is conductivity at left boundary, ka>0, kc >0 
! c     kka=-0.5*ka*dx/condL < 0 
! c     kkc=-0.5*kc*dx/condL < 0
        
        ! eq.1: Butler-Volmer and Ohm's law
        ! 0.5*(phi1 - phi0) = kka*exp(aa*eta) + kkc*exp(-ac*eta)
        ! eq.2: electrode potenital plus a standard potential
        ! phiLe = 0.5*(phi0 + phi1) + E^0 = constant
        ! eq.3: overpotential at electrode
        ! (phi0 + phi1)/2 = -eta + phiLe   (assume phi = potential - potential electrode). 

! c     phi1-phi0 = phi1-(-phi1 - 2*eta + 2*phiLe) = 2*phi1 + 2*eta - 2*phiLe  
! c     (phi1 - phiLe + eta) = kka*exp(aa*eta) + kkc*exp(-ac*eta) 
! c     -eta + kkc*exp(-ac*eta) + kka*exp(aa*eta) = phi1 - phiLe  
! c     Solve eta_new using linearization (Newton method) near  previous eta
! c     az*(eta_new - eta) +bz = cz with:
! c           cz = phi1-phiLe
! c           bz = -eta+kkc*exp(-ac*eta)+kka*exp(aa*eta)
! c           az = -1 - kkc*ac*exp(-ac*eta) +kka*aa*exp(aa*eta) 
! c     eta_new = eta + (cz-bz)/az

! c     phi0_new = 2*(-eta_new+phiLe)-phi1
    
! c     compute aa, ac, ka, kc and ke
!$omp parallel do private(cz,fac,kka,kkc,kkae,kkce,bz,az)
    do k=1,km
        do j=1,jm
            do i=0,im
                conduct(i,j,k) = condfac*0.5d0*(c(i,j,k,2) + c(i+1,j,k,2)) 
            end do

            cz=phi(1,j,k)-phiLe
            ! If needed perform the next 8 statements (Newton iteration) multiple times
            fac = -0.5d0*dx1_i(0)/conduct(0,j,k)
            kka = fac*ka(j,k)
            kkc = fac*kc(j,k)
            kkae = kka*exp(aa*eta(j,k))
            kkce = kkc*exp(-ac*eta(j,k))      
            bz   = -eta(j,k) + kkce + kkae 
            az   = -1d0 - ac*kkce + aa*kkae
            eta(j,k)  = eta(j,k) + (cz-bz)/az
            phiL(j,k) = -eta(j,k) + phiLe

            ! Impose phiLe at left side, like a Dirichlet boundary condition
            ! Second order
            phi(0,j,k)    = 2d0*phiL(j,k) - phi(1,j,k)
            phi(im+1,j,k) = 2d0*phiR - phi(im,j,k) 
        end do
    end do

     phi(0:im+1,0,1:km)=phi(0:im+1,1,1:km)
     phi(0:im+1,jm+1,1:km)=phi(0:im+1,jm,1:km)
     phi(0:im+1,0:jm+1,0)=phi(0:im+1,0:jm+1,km)
     phi(0:im+1,0:jm+1,km+1)=phi(0:im+1,0:jm+1,1)


    if (bubble) then
        iter3=1
        do k=0,km+1
            do j=0,jm+1
                do i=0,im+1
                    storage_phi(i,j,k)=phi(i,j,k)
                end do
            end do
        end do

        call bulb_phi

        do k=0,km+1
            do j=0,jm+1
                do i=0,im+1
                    Error_phi(i,j,k) = abs(phi(i,j,k) - storage_phi(i,j,k))              
                end do
            end do
        end do

        Max_Err_phi=maxval(Error_phi)

! c        if (iter3.gt.5000) then
! c         write(*,*)'Iteration problem in bulb_phi'
! c         stop
! c        end if
    end if
    conduct(0:im,0,1:km)=conduct(0:im,1,1:km)
    conduct(0:im,jm+1,1:km)=conduct(0:im,jm,1:km)
    conduct(0:im,0:jm+1,0)=conduct(0:im,0:jm+1,1)
    conduct(0:im,0:jm+1,km+1)=conduct(0:im,0:jm+1,km)

    eta(0,1:km)=eta(1,1:km)
    eta(jm+1,1:km)=eta(jm,1:km)
    eta(0:jm+1,0)=eta(0:jm+1,1)
    eta(0:jm+1,km+1)=eta(0:jm+1,km)

end subroutine bou_phi



subroutine bulb_phi()
    use var2
    implicit none
    integer ::   o,mm
    real(8) ::  c00,c01,c10,c11,xd,yd,zd,c0,c1

    mm=1

    do o=1,change-1
       i=x_IB(o)
       j=y_IB(o)
       k=z_IB(o)

       xd=(x_prob_c(o,mm)-x1(nx_prob_C(o,mm)))/(dx1_i(1))
       yd=(y_prob_c(o,mm)-x2(ny_prob_C(o,mm)))/(dx2_j(1))
       zd=(z_prob_c(o,mm)-x3(nz_prob_C(o,mm)))/(dx3_k(1))


       c00=(1d0-xd)* &
      phi(nx_prob_C(o,mm),ny_prob_C(o,mm),nz_prob_C(o,mm))+ &
      xd*phi(nx_prob_C(o,mm)+1,ny_prob_C(o,mm),nz_prob_C(o,mm))

       c01=(1d0-xd)* &
      phi(nx_prob_C(o,mm),ny_prob_C(o,mm),nz_prob_C(o,mm)+1)+ &
      xd*phi(nx_prob_C(o,mm)+1,ny_prob_C(o,mm),nz_prob_C(o,mm)+1)

       c10=(1d0-xd)*&
      phi(nx_prob_C(o,mm),ny_prob_C(o,mm)+1,nz_prob_C(o,mm))+&
      xd*phi(nx_prob_C(o,mm)+1,ny_prob_C(o,mm)+1,nz_prob_C(o,mm))

       c11=(1d0-xd)* &
      phi(nx_prob_C(o,mm),ny_prob_C(o,mm)+1,nz_prob_C(o,mm)+1)+ &
      xd*phi(nx_prob_C(o,mm)+1,ny_prob_C(o,mm)+1,nz_prob_C(o,mm)+1)

       c0=c00*(1d0-yd)+c10*yd
       c1=c01*(1d0-yd)+c11*yd

       phi(i,j,k)=c0*(1d0-zd)+c1*zd


    end do

end subroutine bulb_phi




subroutine prep_phi()
    use var2
    implicit none   
    real(8) :: coef0_tmp

! determine matrix coefficients normalized by minus the diagonal element
! only valid for solution of div(c * grad(phi))=0
!$omp parallel do private(i,j,k)
    do k=0,km
        do j=0,jm
            do i=0,im
                work1(i,j,k) = (c(i,j,k,2) + c(i+1,j,k,2))/dx1_i(1)
                work2(i,j,k) = (c(i,j,k,2) + c(i,j+1,k,2))/dx2_j(1)
                work3(i,j,k) = (c(i,j,k,2) + c(i,j,k+1,2))/dx3_k(1)
            end do
        end do
    end do
    
!$omp parallel do private(coef0_tmp,i,j,k)
    do k=1,km
        do j=1,jm
            do i=1,im
                coef0_tmp = -work1(i-1,j,k) - work1(i,j,k) &
                            -work2(i,j-1,k) - work2(i,j,k) &
                            -work3(i,j,k-1) - work3(i,j,k)

                coef1_phi(i,j,k) = work1(i-1,j,k)/coef0_tmp
                coef2_phi(i,j,k) = work1(i,j,k)  /coef0_tmp
                coef3_phi(i,j,k) = work2(i,j-1,k)/coef0_tmp
                coef4_phi(i,j,k) = work2(i,j,k)  /coef0_tmp
                coef5_phi(i,j,k) = work3(i,j,k-1)/coef0_tmp
                coef6_phi(i,j,k) = work3(i,j,k)  /coef0_tmp
            end do
        end do
    end do
    
! compute kinetic prefactors
    do k=1,km
        do j=1,jm
            ka(j,k) = ka0*(0.5d0*(c(0,j,k,2)+c(1,j,k,2))/c_ref(2))*&
                        sqrt(abs(0.5d0*(c(0,j,k,1)+c(1,j,k,1)))/c_ref(1)) 
            kc(j,k) = kc0*0.5d0*(c(0,j,k,3)+c(1,j,k,3))/c_ref(3) 
        end do
    end do

    ka(0,0:km)=ka(1,0:km)
    ka(jm+1,0:km)=ka(jm,0:km)
    ka(0:jm+1,0)=ka(0:jm+1,1)
    ka(0:jm+1,km+1)=ka(0:jm+1,km)

    kc(0,0:km)=kc(1,0:km)
    kc(jm+1,0:km)=kc(jm,0:km)
    kc(0:jm+1,0)=kc(0:jm+1,1)
    kc(0:jm+1,km+1)=kc(0:jm+1,km)
end subroutine prep_phi




subroutine spatial_phi()
    use var2
    implicit none   
    ! real(8) :: phi_temp(0:im+1,0:jm+1,0:km+1)

    ! phi_temp = phi
    storage_phi = phi

    ! Jacobi iteration step
    ! Relaxation factor 1, so phi(i) at the rhs cancels out 
!$omp parallel do private(i,j,k)
    do k=1,km
        do j=1,jm
            do i=1,im
                if (typ(i,j,k).eq.1) then
                    storage_phi(i,j,k) = -(coef1_phi(i,j,k)*phi(i-1,j,k) + &
                                            coef2_phi(i,j,k)*phi(i+1,j,k) + &
                                            coef3_phi(i,j,k)*phi(i,j-1,k) + &
                                            coef4_phi(i,j,k)*phi(i,j+1,k) + &
                                            coef5_phi(i,j,k)*phi(i,j,k-1) + &
                                            coef6_phi(i,j,k)*phi(i,j,k+1))
                end if
            end do
        end do
    end do

    phi = storage_phi

    ! phi=phi_temp
end subroutine spatial_phi



subroutine bou_c()
    use var2
    implicit none   
    integer ::   n
    real(8) :: dphi, ER_c, Max_Err_c

    Error_c(0:im+1,0:jm+1,0:km+1,1:3)=0d0

! c     No values of c and phi are needed in the corner (dummy) points 
! c     (0,0),(im+1,0),(0,jm+1),(im+1,jm+1). 
    do k=1,km
        do j=1,jm   
            dphi = conduct(0,j,k)*(phi(1,j,k) - phi(0,j,k))/Far           
            c(0,j,k,1) = c(1,j,k,1) + 0.5d0*dphi/dif(1)
            c(0,j,k,2) = c(1,j,k,2) + 0.5d0*dphi/dif(2)
            c(0,j,k,3) = c(1,j,k,3) - dphi/dif(3)
        end do
    end do

    do n=1,nmax
!$omp parallel do private(j,k)
        do k=1,km
            do j=1,jm 
                c(im+1,j,k,n) = 2d0*c_ref(n) - c(im,j,k,n)
            end do
        end do
        c(0:im+1,0,1:km,n)=2d0*c_ref(n)-c(0:im+1,1,1:km,n)        
        c(0:im+1,jm+1,1:km,n)=c(0:im+1,jm,1:km,n)
        c(0:im+1,0:jm+1,0,n)=c(0:im+1,0:jm+1,km,n)        
        c(0:im+1,0:jm+1,km+1,n)=c(0:im+1,0:jm+1,1,n)
    end do

    iter4=0
    if (bubble) then
        ER_c=1.0d-3
        Max_Err_c=5d0
        do while (Max_Err_c .gt. ER_c)
            do n=2,3
                do k=0,km+1
                    do j=0,jm+1
                        do i=0,im+1
                            c_old(i,j,k,n)=c(i,j,k,n)
                        end do
                    end do
                end do
            end do
                            
            call bulb_c
            iter4=iter4+1
            do n=2,3
                do k=0,km+1
                    do j=0,jm+1
                        do i=0,im+1
                            Error_c(i,j,k,n)=abs(c(i,j,k,n)-c_old(i,j,k,n))             
                        end do
                    end do
                end do
            end do
            Max_Err_c=maxval(Error_c)
        end do
    end if
end subroutine bou_c



subroutine spatial_c()
    use var2
    implicit none   
    integer ::   n
    real(8) ::  ar1,ar2,ar3
    
    if ( rungIter.eq.1 ) then
        storage_c = c
    end if

    do n=1,nmax

!$omp parallel do private(i,j,k,ar1,ar2,ar3)
        do k=0,km
            do j=0,jm 
                do i=0,im
                    work1(i,j,k)=dif(n)*(c(i+1,j,k,n)-c(i,j,k,n))/dx1_i(1)
                    work2(i,j,k)=dif(n)*(c(i,j+1,k,n)-c(i,j,k,n))/dx2_j(1) 
                    work3(i,j,k)=dif(n)*(c(i,j,k+1,n)-c(i,j,k,n))/dx3_k(1) 
                end do
            end do
        end do

!$omp parallel do private(ar1,ar2,ar3)
        do k=1,km
            do j=1,jm
                do i=1,im
                    if (typ(i,j,k).eq.1) then
                        rhs(i,j,k)= &
                            (work1(i,j,k)-work1(i-1,j,k))/dx1(i)+ &
                            (work2(i,j,k)-work2(i,j-1,k))/dx2(j)+ &
                            (work3(i,j,k)-work3(i,j,k-1))/dx3(k)
                        ar1=0.5d0*(u1(i,j,k)+u1(i-1,j,k))
                        ar2=0.5d0*(u2(i,j,k)+u2(i,j-1,k))
                        ar3=0.5d0*(u3(i,j,k)+u3(i,j,k-1))
                        if (ar1.gt.0d0) then
                            rhs(i,j,k)=rhs(i,j,k)-ar1*(c(i,j,k,n) &
                                    -c(i-1,j,k,n))/dx1_i(0)
                        else
                            rhs(i,j,k)=rhs(i,j,k)-ar1*(c(i+1,j,k,n) &
                                    -c(i,j,k,n))/dx1_i(0)
                        end if

                        if (ar2.gt.0d0) then
                            rhs(i,j,k)=rhs(i,j,k)-ar2*(c(i,j,k,n)&
                                    -c(i,j-1,k,n))/dx2_j(0)
                        else
                            rhs(i,j,k)=rhs(i,j,k)-ar2*(c(i,j+1,k,n) &
                                    -c(i,j,k,n))/dx2_j(0)
                        end if

                        if (ar3.gt.0d0) then
                            rhs(i,j,k)=rhs(i,j,k)-ar3*(c(i,j,k,n) &
                                    -c(i,j,k-1,n))/dx3_k(0)
                        else
                            rhs(i,j,k)=rhs(i,j,k)-ar3*(c(i,j,k+1,n) &
                                    -c(i,j,k,n))/dx3_k(0)
                        end if
                    end if
                end do
            end do
        end do

!$omp parallel do private(i,j,k)
        do k=1,km
            do j=1,jm
                do i=1,im
                    if (typ(i,j,k).eq.1) then
                        c(i,j,k,n)=storage_c(i,j,k,n) + dtime*rhs(i,j,k)*RK_alpha(rungIter)
                    end if
                end do
            end do
        end do
    end do


    balance_h2=0d0
    do k=1,km
        do j=1,jm
            do i=1,im
                if (typ(i,j,k).eq.1) then
                    balance_h2=balance_h2+(c(i,j,k,1)-storage_c(i,j,k,1)) &
                    /dtime*dx1_i(0)*dx2_j(0)*dx3_k(0)
                end if
            end do
        end do
    end do
end subroutine spatial_c




subroutine bulb_c()
    use var2
    implicit none   
    integer ::   o,n,mm
    real(8) ::  xd,yd,zd,c00,c01,c10,c11,c0,c1,beta_c

    mm=1

    do o=1,change-1
       i=x_IB(o)
       j=y_IB(o)
       k=z_IB(o)

       xd=(x_prob_c(o,mm)-x1(nx_prob_C(o,mm)))/(dx1_i(1))
       yd=(y_prob_c(o,mm)-x2(ny_prob_C(o,mm)))/(dx2_j(1))
       zd=(z_prob_c(o,mm)-x3(nz_prob_C(o,mm)))/(dx3_k(1))

       do n=1,3

       c00=(1d0-xd)* &
     c(nx_prob_C(o,mm),ny_prob_C(o,mm),nz_prob_C(o,mm),n)+ &
      xd*c(nx_prob_C(o,mm)+1,ny_prob_C(o,mm),nz_prob_C(o,mm),n)

       c01=(1d0-xd)* &
      c(nx_prob_C(o,mm),ny_prob_C(o,mm),nz_prob_C(o,mm)+1,n)+ &
      xd*c(nx_prob_C(o,mm)+1,ny_prob_C(o,mm),nz_prob_C(o,mm)+1,n)

       c10=(1d0-xd)* &
      c(nx_prob_C(o,mm),ny_prob_C(o,mm)+1,nz_prob_C(o,mm),n)+ &
      xd*c(nx_prob_C(o,mm)+1,ny_prob_C(o,mm)+1,nz_prob_C(o,mm),n)

       c11=(1d0-xd)* &
      c(nx_prob_C(o,mm),ny_prob_C(o,mm)+1,nz_prob_C(o,mm)+1,n)+ &
      xd*c(nx_prob_C(o,mm)+1,ny_prob_C(o,mm)+1,nz_prob_C(o,mm)+1,n)

       c0=c00*(1d0-yd)+c10*yd
       c1=c01*(1d0-yd)+c11*yd

       c_s(n)=c0*(1d0-zd)+c1*zd

       end do 

       beta_c=rb(o)/dx1(1)
       c(i,j,k,1)=(1d0 + beta_c)*c_ref(1)-(beta_c)*c_s(1)


       c(i,j,k,2)=c_s(2)
       c(i,j,k,3)=c_s(3)

       end do
end subroutine bulb_c




subroutine bulb_flux()
    use var2
    implicit none   
    integer ::   o
    real(8) :: term1,term2

    bubbleMassflux = 0d0

    do o=1,change-1
       i=x_IB(o)
       j=y_IB(o)
       k=z_IB(o)

       if (x1_minus(o) .eq. 1) then
          term1=-0.5d0*(c(i,j,k,1)+c(i-1,j,k,1))* &
                  u1(i-1,j,k)*dx2(j)*dx3(k)
          term2=dif(1)*(c(i,j,k,1)-c(i-1,j,k,1)) &
                 /dx1_i(i-1)*dx2(j)*dx3(k)
          if (i .eq. nx_start) then
           term1=0d0
           term2=0d0
          end if
          bubbleMassflux = bubbleMassflux + term2
       end if
       if (x1_plus(o) .eq. 1) then
          term1=-(-0.5d0*(c(i+1,j,k,1)+c(i,j,k,1))* &
                  u1(i,j,k)*dx2(j)*dx3(k))
          term2=-(dif(1)*(c(i+1,j,k,1)-c(i,j,k,1)) &
                 /dx1_i(i)*dx2(j)*dx3(k))
          bubbleMassflux = bubbleMassflux + term2
       end if

       if (x2_minus(o) .eq. 1) then
          term1=-0.5d0*(c(i,j,k,1)+c(i,j-1,k,1))* &
                  u2(i,j-1,k)*dx1(i)*dx3(k)
          term2=dif(1)*(c(i,j,k,1)-c(i,j-1,k,1)) &
                 /dx2_j(j-1)*dx1(i)*dx3(k)
          bubbleMassflux=bubbleMassflux+term2
       end if
       if (x2_plus(o) .eq. 1) then
          term1=-(-0.5d0*(c(i,j+1,k,1)+c(i,j,k,1))*&
                  u2(i,j,k)*dx1(i)*dx3(k))
          term2=-(dif(1)*(c(i,j+1,k,1)-c(i,j,k,1)) &
                 /dx2_j(j)*dx1(i)*dx3(k))
          bubbleMassflux=bubbleMassflux+term2
       end if

       if (x3_minus(o) .eq. 1) then
          term1=-0.5d0*(c(i,j,k,1)+c(i,j,k-1,1))* &
                  u3(i,j,k-1)*dx1(i)*dx2(j)
          term2=dif(1)*(c(i,j,k,1)-c(i,j,k-1,1)) &
                /dx3_k(k-1)*dx1(i)*dx2(j)
          bubbleMassflux=bubbleMassflux+term2
       end if
       if (x3_plus(o) .eq. 1) then
          term1=-(-0.5d0*(c(i,j,k+1,1)+c(i,j,k,1))* &
                  u3(i,j,k)*dx1(i)*dx2(j))
          term2=-(dif(1)*(c(i,j,k+1,1)-c(i,j,k,1)) &
                 /dx3_k(k)*dx1(i)*dx2(j))
          bubbleMassflux=bubbleMassflux+term2
       end if
    end do

end subroutine bulb_flux




subroutine domain_flux()
    use var2
    implicit none   
    integer ::   o
    real(8) :: term1,term2,TERM3,TERM4,TERM5,TERM6, &
        !   A12(1:jm,1:km),B12(1:jm,1:km), &
        !   A21(1:im,1:km),B21(1:im,1:km), &
        !   A22(1:im,1:km),B22(1:im,1:km), &
        !   A32(1:im,1:jm),B32(1:im,1:jm), &
        !   A31(1:im,1:jm),B31(1:im,1:jm), &
          Error

    do k=1,km
    do j=1,jm
       A12(j,k)=-(dif(1)*(c(1,j,k,1)-c(0,j,k,1))/dx1_i(0)&
            *dx2(j)*dx3(k))
       B12(j,k)=dif(1)*(c(im+1,j,k,1)-c(im,j,k,1))/dx1_i(im) &
             *dx2(j)*dx3(k)
    end do
    end do

    do o=1,change-1
       i=x_IB(o)
       if (i .eq. nx_start) then
        A12(y_IB(o),z_IB(o))=0d0
       end if
    end do


    term1=sum(A12)
    term2=sum(B12)


    do k=1,km
    do i=1,im

       A21(i,k)=-(-0.5d0*(c(i,1,k,1)+c(i,0,k,1)) &
                 *u2(i,0,k)*dx1(i)*dx3(k))
       B21(i,k)=-0.5d0*(c(i,jm+1,k,1)+c(i,jm,k,1)) &
                 *u2(i,jm,k)*dx1(i)*dx3(k)
       A22(i,k)=-(dif(1)*(c(i,1,k,1)-c(i,0,k,1))/dx2_j(0) &
              *dx1(i)*dx3(k))
       B22(i,k)=dif(1)*(c(i,jm+1,k,1)-c(i,jm,k,1))/dx2_j(jm) &
              *dx1(i)*dx3(k)
    end do
    end do


    term3=sum(A21)+sum(A22)
    term4=sum(b21)+sum(B22)

    do j=1,jm
    do i=1,im
       A31(i,j)=-(-0.5d0*(c(i,j,1,1)+c(i,j,0,1)) &
              *u3(i,j,0)*dx1(i)*dx2(j))
       B31(i,j)=-0.5d0*(c(i,j,km+1,1)+c(i,j,km,1)) &
              *u3(i,j,km)*dx1(i)*dx2(j)
       A32(i,j)=-(dif(1)*(c(i,j,1,1)-c(i,j,0,1))/dx3_k(0) &
              *dx1(i)*dx2(j))
       B32(i,j)=dif(1)*(c(i,j,km+1,1)-c(i,j,km,1))/dx3_k(km) &
              *dx1(i)*dx2(j)
    end do
    end do

    term5=sum(A31)+sum(A32)
    term6=sum(B31)+sum(B32)

    Error=abs(term1+term2+term3+term4+term5+term6+bubbleMassflux-balance_h2)


1   format(14E36.24)
    if ( (m1.eq.1).and.(m2.eq.1) ) then
        open(unit=1,file='./output/H2_balance.dat')
2       format(14A36)
        write(1,2) 'time', 'term1', 'term2', 'term3', 'term4', &
                    'term5', 'term6', 'totBubbleFlux', 'H2_balance', &
                    'H2_balanceError', 'A31', 'A32', 'B31', 'B32'
        write(1,1) time, term1, term2, term3, term4, term5, term6,  &
                    bubbleMassflux, -balance_h2, Error, sum(A31), sum(A32),  &
                    sum(B31), sum(B32)
        close(1,status='keep')
    end if
    if ( (m1.gt.1).and.(m2.eq.m2_max) ) then
        open(unit=1,file='./output/H2_balance.dat',position='append')
        write(1,1) time, term1, term2, term3, term4, term5, term6,  &
                    bubbleMassflux, -balance_h2, Error, sum(A31), sum(A32),  &
                    sum(B31), sum(B32)
        close(1,status='keep')
    end if

end subroutine domain_flux




subroutine growing()
    use var2
    implicit none   
    integer ::   o
    real(8) :: delta_v,P_IB_fluid,P_outside, &
          P_bubble
          
! c     average pressure over the IB_fluid cells
    P_IB_fluid=0d0
    do o=1,change-1
       i=x_IB(o)
       j=y_IB(o)
       k=z_IB(o)
       P_IB_fluid=P_IB_fluid+p(i,j,k)
    end do

    P_outside=P_IB_fluid/(change-1)

! c     laplace pressure equation
    P_bubble = P_outside + (2d0/radius)*sigma + 1d5


    delta_v = -1d0*bubbleMassflux*Rid*temp*dtime/P_bubble

    radius0=radius

    radius=(3d0/4d0/pi*delta_v+radius**3d0)**(1d0/3d0)
    drdt=(radius-radius0)/dtime

end subroutine growing





subroutine bulb_flux_force()
    use var2
    implicit none   
    integer ::   o,counter_write
    real(8) :: term1,term2,term3

100   format(4I5,1000e14.6)
    total=0d0
    counter_write=0
    do o=1,change-1
       i=x_IB(o)
       j=y_IB(o)
       k=z_IB(o)

        if (x1_minus(o) .eq. 1) then
            term1=-(0.5d0*(u1(i-1,j,k)+u1(i,j,k)))**2d0*dx2(j)*dx3(k)
            term2=-p(i,j,k)/rhoKOH*dx2(j)*dx3(k)
            term3=nuKOH*(u1(i,j,k)-u1(i-1,j,k))/dx1(i)*dx2(j)*dx3(k)

            if (i .eq. nx_start) then
                term1=0d0
                term2=0d0
                term3=0d0
                counter_write=counter_write+1
                if (counter_write .eq. 1) then
                    open(unit=51,file='node.dat')
                    write(51,100)i,j,k,-1,term1,term2,term3
                else 
                    open(unit=51,file='node.dat',position='append')
                    write(51,100)i,j,k,-1,term1,term2,term3
                end if
                close(51)
            end if
            total=total+term1+term2+term3
        end if

        if (x1_plus(o) .eq. 1) then
            term1=-(-(0.5d0*(u1(i+1,j,k)+u1(i,j,k)))**2d0*dx2(j)*dx3(k))
            term2=-(-p(i+1,j,k)/rhoKOH*dx2(j)*dx3(k))
            term3=-(nuKOH*(u1(i+1,j,k)-u1(i,j,k))/dx1(i+1)*dx2(j)*dx3(k)) 
            
            if (i .eq. nx_start) then
                counter_write=counter_write+1
                if (counter_write .eq. 1) then
                    open(unit=51,file='node.dat')
                    write(51,100)i,j,k,1,term1,term2,term3
                else 
                    open(unit=51,file='node.dat',position='append')
                    write(51,100)i,j,k,1,term1,term2,term3
                end if
                close(51)
            end if
            total=total+term1+term2+term3
        end if

        if (x2_minus(o) .eq. 1) then
            term1=-0.5d0*(u1(i,j,k)+u1(i,j-1,k)) &
                    *0.5d0*(u2(i,j-1,k)+u2(i+1,j-1,k))*dx1(i)*dx3(k)
            term2=0d0
            term3=nuKOH*(u1(i,j,k)-u1(i,j-1,k))/dx2_j(j-1)*dx1(i)*dx3(k)    
            if (i .eq. nx_start) then
                counter_write=counter_write+1
                if (counter_write .eq. 1) then
                    open(unit=51,file='node.dat')
                    write(51,100)i,j,k,-2,term1,term2,term3
                else 
                    open(unit=51,file='node.dat',position='append')
                    write(51,100)i,j,k,-2,term1,term2,term3
                end if
                close(51)
            end if
            total=total+term1+term2+term3
        end if

        if (x2_plus(o) .eq. 1) then
            term1=-(-0.5d0*(u1(i,j,k)+u1(i,j+1,k)) &
                    *0.5d0*(u2(i,j+1,k)+u2(i+1,j+1,k))*dx1(i)*dx3(k))
            term2=0d0  
            term3=-(nuKOH*(u1(i,j+1,k)-u1(i,j,k))/dx2_j(j)*dx1(i)*dx3(k)) 
            if (i .eq. nx_start) then
                counter_write=counter_write+1
                if (counter_write .eq. 1) then
                    open(unit=51,file='node.dat')
                    write(51,100)i,j,k,2,term1,term2,term3
                else 
                    open(unit=51,file='node.dat',position='append')
                    write(51,100)i,j,k,2,term1,term2,term3
                end if
                close(51)
            end if
            total=total+term1+term2+term3
        end if

       if (x3_minus(o) .eq. 1) then
          term1=-0.5d0*(u1(i,j,k)+u1(i,j,k-1)) &
        *0.5d0*(u3(i,j,k-1)+u3(i+1,j,k-1))*dx1(i)*dx2(j)

          term2=0d0
  
          term3=nuKOH*(u1(i,j,k)-u1(i,j,k-1))/dx3_k(k-1)*dx1(i)*dx2(j) 

          if (i .eq. nx_start) then

           counter_write=counter_write+1

           if (counter_write .eq. 1) then
            open(unit=51,file='node.dat')
            write(51,100)i,j,k,-3,term1,term2,term3
           else 
            open(unit=51,file='node.dat',position='append')
            write(51,100)i,j,k,-3,term1,term2,term3
           end if
           close(51)
          end if
          total=total+term1+term2+term3
       end if

       if (x3_plus(o) .eq. 1) then
          term1=-(-0.5d0*(u1(i,j,k)+u1(i,j,k+1)) &
        *0.5d0*(u3(i,j,k+1)+u3(i+1,j,k+1))*dx1(i)*dx2(j))
          term2=0d0
          term3=-(nuKOH*(u1(i,j,k+1)-u1(i,j,k))/dx3_k(k)*dx1(i)*dx2(j))   
          if (i .eq. nx_start) then
           counter_write=counter_write+1
           if (counter_write .eq. 1) then
            open(unit=51,file='node.dat')
            write(51,100)i,j,k,3,term1,term2,term3
           else 
            open(unit=51,file='node.dat',position='append')
            write(51,100)i,j,k,3,term1,term2,term3
           end if
           close(51)
          end if
          total=total+term1+term2+term3
       end if
    end do

end subroutine bulb_flux_force





subroutine domain_flux_force()
    use var2
    implicit none   
    integer ::   o       
    real(8) :: term1,term2,term3,term4,term5,term6, &
        !   u_i(1:im,1:km),u_o(1:im,1:km),   &
        !   w_i(1:im,1:km),w_o(1:im,1:km), &
        !   u_f(1:im,1:jm),u_b(1:im,1:jm), &
        !   v_f(1:im,1:jm),v_b(1:im,1:jm), &
        !   A11(1:jm,1:km),B11(1:jm,1:km), &
        !   A12(1:jm,1:km),B12(1:jm,1:km), &
        !   A13(1:jm,1:km),B13(1:jm,1:km), &
        !   A21(1:im,1:km),B21(1:im,1:km), &
        !   A22(1:im,1:km),B22(1:im,1:km), &
        !   A23(1:im,1:km),B23(1:im,1:km), &
        !   A31(1:im,1:jm),B31(1:im,1:jm), &
        !   A32(1:im,1:jm),B32(1:im,1:jm), &
        !   A33(1:im,1:jm),B33(1:im,1:jm),
          Error


    do k=1,km
    do j=1,jm

       A11(j,k)=0d0
       B11(j,k)=0d0

       A12(j,k)=-(-0.5d0*(p(0,j,k)+p(1,j,k))*dx2(j)*dx3(k)/rhoKOH)
       B12(j,k)=-0.5d0*(p(im,j,k)+p(im+1,j,k))*dx2(j)*dx3(k)/rhoKOH

       A13(j,k)=0d0
       B13(j,k)=0d0


    end do
    end do

    do o=1,change-1
       i=x_IB(o)
       if (i .eq. nx_start) then
        A12(y_IB(o),z_IB(o))=0d0
       end if
    end do


    term1=sum(A11)+sum(A12)+sum(A13)
    term2=sum(B11)+sum(B12)+sum(B13)


    do k=1,km
    do i=1,im

       u_i(i,k)=0.25d0*(u1(i-1,0,k)+u1(i-1,1,k) &
              +u1(i,0,k)+u1(i,1,k))
       u_o(i,k)=0.25d0*(u1(i-1,jm,k)+u1(i-1,jm+1,k) &
              +u1(i,jm,k)+u1(i,jm+1,k))
       w_i(i,k)=0.25d0*(u3(i,0,k-1)+u3(i,1,k-1)+u3(i,0,k)+u3(i,1,k))
       w_o(i,k)=0.25d0*(u3(i,jm,k-1)+u3(i,jm+1,k-1) &
            +u3(i,jm,k)+u3(i,jm+1,k))

    end do
    end do


    do k=1,km
    do i=1,im

       A21(i,k)=-(-u_i(i,k)*u2(i,0,k)*dx1(i)*dx3(k))
       B21(i,k)=-u_o(i,k)*u2(i,jm,k)*dx1(i)*dx3(k)
 
       A22(i,k)=0d0
       B22(i,k)=0d0

       A23(i,k)=-(0.5d0*nuKOH*((u1(i,1,k)+u1(i-1,1,k)) &
               -(u1(i,0,k)+u1(i-1,0,k)))*dx1(i)*dx3(k)/dx2_j(0))
       B23(i,k)=0.5d0*nuKOH*((u1(i,jm+1,k)+u1(i-1,jm+1,k)) &
               -(u1(i,jm,k)+u1(i-1,jm,k)))*dx1(i)*dx3(k)/dx2_j(jm)


    end do
    end do

    term3=sum(A21)+sum(A22)+sum(A23)
    term4=sum(B21)+sum(B22)+sum(B23)


    do j=1,jm
    do i=1,im

       u_f(i,j)=0.25d0*(u1(i-1,j,0)+u1(i-1,j,1) &
               +u1(i,j,0)+u1(i,j,1))
       u_b(i,j)=0.25d0*(u1(i-1,j,km)+u1(i-1,j,km+1) &
               +u1(i,j,km)+u1(i,j,km+1))
       v_f(i,j)=0.25d0*(u2(i,j-1,0)+u2(i,j-1,1) &
               +u2(i,j,0)+u2(i,j,1))
       v_b(i,j)=0.25d0*(u2(i,j-1,km)+u2(i,j-1,km+1) &
               +u2(i,j,km)+u2(i,j,km+1))

    end do
    end do


    do j=1,jm
    do i=1,im

       A31(i,j)=-(-u_f(i,j)*u3(i,j,0)*dx1(i)*dx2(j))
       B31(i,j)=-u_b(i,j)*u3(i,j,km)*dx1(i)*dx2(j)

       A32(i,j)=0d0
       B32(i,j)=0d0

       A33(i,j)=-(0.5d0*nuKOH*((u1(i,j,1)+u1(i-1,j,1))-(u1(i,j,0) &
               +u1(i-1,j,0)))*dx1(i)*dx2(j)/dx3_k(0))
       B33(i,j)=0.5d0*nuKOH*((u1(i,j,km+1)+u1(i-1,j,km+1))-(u1(i,j,km) &
               +u1(i-1,j,km)))*dx1(i)*dx2(j)/dx3_k(km)

    end do
    end do

    term5=sum(A31)+sum(A32)+sum(A33)
    term6=sum(B31)+sum(B32)+sum(B33)


    Error=abs(term1+term2+term3+term4+term5+term6+total)
    write(*,*)'Error of force balance',Error
    write(*,*)term1,term2,term3,term4,term5,term6,total

end subroutine domain_flux_force




subroutine finding_u_IB()
    use var2
    implicit none
    integer ::   o

!   Everything is marked as bubble face
    typ_u_IB(0:im+1,0:jm+1,0:km+1)=0
    do o=1,change-1
       i=x_IB(o)
       j=y_IB(o)
       k=z_IB(o)
!   If neighbouring cell is fluid-cell (id=1), mark the face as IB-face (id=1)
        if (typ(i+1,j,k).eq.1) then
            typ_u_IB(i,j,k)=1
        else
!   If neighbouring cell is an IB-cell (id=1), mark the face as IB-face (id=1)
            if (typ_IB(i+1,j,k).eq.1) then
                typ_u_IB(i,j,k)=1
            end if
        end if
       if (typ(i-1,j,k).eq.1) then
            typ_u_IB(i-1,j,k)=1
        else
            if (typ_IB(i-1,j,k).eq.1) then
                typ_u_IB(i-1,j,k)=1
            end if
        end if
    end do
  
!   Count number of u-faces and store the coordinates
    o_u=0
    do k=0,km+1
        do j=0,jm+1
            do i=0,im+1
                if (typ_u_IB(i,j,k).eq.1) then
                    o_u = o_u+1
                    x_IB_U(o_u)=i
                    y_IB_U(o_u)=j
                    z_IB_U(o_u)=k
                end if
            end do
        end do
    end do
end subroutine finding_u_IB


subroutine probe_IB_u()
    use var2
    implicit none
    integer ::   o,mm
    real(8) ::  rnode_U

! c        Neumann BC in Spherical  coordinates

! c        phiS = angle phi in spherical coordinates.
! c        not to be confused with phi the potential
! c        if x=0 then atan(y,x)=pi/2. This is not true however,
! c        if y is negative, then it should be 3pi/2.

    mm=1
    
    do o=1,o_u
       i=x_IB_U(o)
       j=y_IB_U(o)
       k=z_IB_U(o)

       phiS_U(o)=datan2((x2(j)-yc),(x1_i(i)-xc))

       thetaS_U(o)=acos((x3(k)-zc)/sqrt((x1_i(i)-xc)**2d0+ &
       (x2(j)-yc)**2d0+(x3(k)-zc)**2d0))


       T_ni_U(o)=cos(phiS_U(o))*sin(thetaS_U(o))
       T_nj_U(o)=sin(phiS_U(o))*sin(thetaS_U(o))
       T_nk_U(o)=cos(thetaS_U(o))

       T_ti_U(o)=-sin(phiS_U(o))
       T_tj_U(o)=cos(phiS_U(o))
       T_tk_U(o)=0d0

       T_si_U(o)=cos(phiS_U(o))*cos(thetaS_U(o))
       T_sj_U(o)=sin(phiS_U(o))*cos(thetaS_U(o))
       T_sk_U(o)=-sin(thetaS_U(o))

! c        initialise x-direction
! c        Normal distance between interface and bubble centre
       rnode_U=sqrt((x1_i(i)-xc)**2d0+(x2(j)-yc)**2d0+(x3(k)-zc)**2d0)


! c        normal distance of interface to bubble surface 
! c        (can be negative)
       rb_U(o)=radius-rnode_U

       xb_U(o)=x1_i(i)+rb_U(o)*T_ni_U(o)
       yb_U(o)=x2(j)  +rb_U(o)*T_nj_U(o)
       zb_U(o)=x3(k)  +rb_U(o)*T_nk_U(o)


       x_prob_U(o,mm)=xb_U(o)+dx1_i(0)*T_ni_U(o)
       y_prob_U(o,mm)=yb_U(o)+dx1_i(0)*T_nj_U(o)
       z_prob_U(o,mm)=zb_U(o)+dx1_i(0)*T_nk_U(o)


! c        prob number of cells starts at an interface so only 
! c        +0.5d if that direction is on the level of a node
       nx_prob_UU(o,mm)=floor(x_prob_U(o,mm)/dx1_i(0))
       ny_prob_UU(o,mm)=floor(y_prob_U(o,mm)/dx2_j(0)+0.5d0)
       nz_prob_UU(o,mm)=floor(z_prob_U(o,mm)/dx3_k(0)+0.5d0)

        ! Check whether the probe is within the domain
       call chckProbeBounds(nx_prob_UU(o,mm), 0, im+1,'x-UU')
       call chckProbeBounds(ny_prob_UU(o,mm), 0, jm+1,'y-UU')
       call chckProbeBounds(nz_prob_UU(o,mm), 0, km+1,'z-UU')

       nx_prob_UV(o,mm)=floor(x_prob_U(o,mm)/dx1_i(0)+0.5d0)
       ny_prob_UV(o,mm)=floor(y_prob_U(o,mm)/dx2_j(0))
       nz_prob_UV(o,mm)=floor(z_prob_U(o,mm)/dx3_k(0)+0.5d0)

       ! Check whether the probe is within the domain
       call chckProbeBounds(nx_prob_UV(o,mm),0,im+1,'x-UV')
       call chckProbeBounds(ny_prob_UV(o,mm),0,jm+1,'y-UV')
       call chckProbeBounds(nz_prob_UV(o,mm),0,km+1,'z-UV')

       nx_prob_UW(o,mm)=floor(x_prob_U(o,mm)/dx1_i(0)+0.5d0)
       ny_prob_UW(o,mm)=floor(y_prob_U(o,mm)/dx2_j(0)+0.5d0)
       nz_prob_UW(o,mm)=floor(z_prob_U(o,mm)/dx3_k(0))

        ! Check whether the probe is within the domain
       call chckProbeBounds(nx_prob_UW(o,mm),0,im+1,'x-UW')
       call chckProbeBounds(ny_prob_UW(o,mm),0,jm+1,'y-UW')
       call chckProbeBounds(nz_prob_UW(o,mm),0,km+1,'z-UW')

    end do
end subroutine probe_IB_u



subroutine bulb_u()
    use var2
    implicit none
    integer ::   o,mm
! c     Solving u1  on IB interface

    real(8) :: u1i_ref,u2i_ref,u3i_ref
    real(8) :: uni_ref,uti_ref,usi_ref
    real(8) :: u_star,v_star,w_star
    real(8) :: unis,utis,usis
    real(8) :: uni,uti,usi
    real(8) :: R_IB_U,R_P_U
    real(8) ::  xd,yd,zd,beta_U
    real(8) ::  c00,c10,c01,c11,c0,c1
    ! real(8) :: u1_temp(0:im+1,0:jm+1,0:km+1)

    mm=1
    u1_temp = u1
    do o=1,o_u

        i=x_IB_U(o)
        j=y_IB_U(o)
        k=z_IB_U(o)
! c     x-component u1
    xd=(x_prob_U(o,mm)-x1_i(nx_prob_UU(o,mm)))/(dx1(1))
    yd=(y_prob_U(o,mm)-x2(ny_prob_UU(o,mm)))/(dx2_j(1))
    zd=(z_prob_U(o,mm)-x3(nz_prob_UU(o,mm)))/(dx3_k(1))


    c00=(1d0-xd)* &
     u1(nx_prob_UU(o,mm),ny_prob_UU(o,mm),nz_prob_UU(o,mm))+ &
     xd*u1(nx_prob_UU(o,mm)+1,ny_prob_UU(o,mm),nz_prob_UU(o,mm))

    c01=(1d0-xd)* &
     u1(nx_prob_UU(o,mm),ny_prob_UU(o,mm),nz_prob_UU(o,mm)+1)+ &
     xd*u1(nx_prob_UU(o,mm)+1,ny_prob_UU(o,mm),nz_prob_UU(o,mm)+1)

    c10=(1d0-xd)* &
     u1(nx_prob_UU(o,mm),ny_prob_UU(o,mm)+1,nz_prob_UU(o,mm))+ &
     xd*u1(nx_prob_UU(o,mm)+1,ny_prob_UU(o,mm)+1,nz_prob_UU(o,mm))

    c11=(1d0-xd)* &
     u1(nx_prob_UU(o,mm),ny_prob_UU(o,mm)+1,nz_prob_UU(o,mm)+1)+ &
     xd*u1(nx_prob_UU(o,mm)+1,ny_prob_UU(o,mm)+1,nz_prob_UU(o,mm)+1)

    c0=c00*(1d0-yd)+c10*yd
    c1=c01*(1d0-yd)+c11*yd

    u1i_ref=c0*(1d0-zd)+c1*zd


    ! c     y-component u1

    xd=(x_prob_U(o,mm)-x1(nx_prob_UV(o,mm)))/(dx1_i(1))
    yd=(y_prob_U(o,mm)-x2_j(ny_prob_UV(o,mm)))/(dx2(1))
    zd=(z_prob_U(o,mm)-x3(nz_prob_UV(o,mm)))/(dx3_k(1))

    c00=(1d0-xd)* &
     u2(nx_prob_UV(o,mm),ny_prob_UV(o,mm),nz_prob_UV(o,mm))+ &
     xd*u2(nx_prob_UV(o,mm)+1,ny_prob_UV(o,mm),nz_prob_UV(o,mm))

    c01=(1d0-xd)* &
     u2(nx_prob_UV(o,mm),ny_prob_UV(o,mm),nz_prob_UV(o,mm)+1)+ &
     xd*u2(nx_prob_UV(o,mm)+1,ny_prob_UV(o,mm),nz_prob_UV(o,mm)+1)

    c10=(1d0-xd)* &
     u2(nx_prob_UV(o,mm),ny_prob_UV(o,mm)+1,nz_prob_UV(o,mm))+ &
     xd*u2(nx_prob_UV(o,mm)+1,ny_prob_UV(o,mm)+1,nz_prob_UV(o,mm))

    c11=(1d0-xd)* &
     u2(nx_prob_UV(o,mm),ny_prob_UV(o,mm)+1,nz_prob_UV(o,mm)+1)+ &
     xd*u2(nx_prob_UV(o,mm)+1,ny_prob_UV(o,mm)+1,nz_prob_UV(o,mm)+1)

    c0=c00*(1d0-yd)+c10*yd
    c1=c01*(1d0-yd)+c11*yd

    u2i_ref=c0*(1d0-zd)+c1*zd

! c     z-component u1

    xd=(x_prob_U(o,mm)-x1(nx_prob_UW(o,mm)))/(dx1_i(1))
    yd=(y_prob_U(o,mm)-x2(ny_prob_UW(o,mm)))/(dx2_j(1))
    zd=(z_prob_U(o,mm)-x3_k(nz_prob_UW(o,mm)))/(dx3(1))

    c00=(1d0-xd)* &
     u3(nx_prob_UW(o,mm),ny_prob_UW(o,mm),nz_prob_UW(o,mm))+ &
     xd*u3(nx_prob_UW(o,mm)+1,ny_prob_UW(o,mm),nz_prob_UW(o,mm))

    c01=(1d0-xd)* &
     u3(nx_prob_UW(o,mm),ny_prob_UW(o,mm),nz_prob_UW(o,mm)+1)+ &
     xd*u3(nx_prob_UW(o,mm)+1,ny_prob_UW(o,mm),nz_prob_UW(o,mm)+1)

    c10=(1d0-xd)* &
     u3(nx_prob_UW(o,mm),ny_prob_UW(o,mm)+1,nz_prob_UW(o,mm))+ &
     xd*u3(nx_prob_UW(o,mm)+1,ny_prob_UW(o,mm)+1,nz_prob_UW(o,mm))

    c11=(1d0-xd)* &
     u3(nx_prob_UW(o,mm),ny_prob_UW(o,mm)+1,nz_prob_UW(o,mm)+1)+ &
     xd*u3(nx_prob_UW(o,mm)+1,ny_prob_UW(o,mm)+1,nz_prob_UW(o,mm)+1)

    c0=c00*(1d0-yd)+c10*yd
    c1=c01*(1d0-yd)+c11*yd

    u3i_ref=c0*(1d0-zd)+c1*zd

! c     Translate Cartesian velocities to spherical velocities
! c     Look at bulb.f for specification on each T 
    uni_ref = (u1i_ref - dxcdt)*T_ni_U(o) + &
              (u2i_ref - dycdt)*T_nj_U(o) + &
              (u3i_ref - dzcdt)*T_nk_U(o)

    uti_ref = (u1i_ref - dxcdt)*T_ti_U(o) + &
              (u2i_ref - dycdt)*T_tj_U(o) + &
              (u3i_ref - dzcdt)*T_tk_U(o)

    usi_ref = (u1i_ref - dxcdt)*T_si_U(o)+ &
              (u2i_ref - dycdt)*T_sj_U(o)+ &
              (u3i_ref - dzcdt)*T_sk_U(o)


! c     Dirichlet BC for normal direction, beta can be negative as rb_i can be negative
    beta_U=rb_U(o)/dx1(1)
! c     Velocity at the bubble surface in cartisian coordinate
    u_star = dxcdt + drdt*T_ni_U(o)
    v_star = dycdt + drdt*T_nj_U(o)
    w_star = dzcdt + drdt*T_nk_U(o)


      unis = (u_star - dxcdt)*T_ni_U(o) + &
             (v_star - dycdt)*T_nj_U(o) + &
             (w_star - dzcdt)*T_nk_U(o)

      utis = (u_star - dxcdt)*T_ti_U(o) + &
             (v_star - dycdt)*T_tj_U(o) + &
             (w_star - dzcdt)*T_tk_U(o)

      usis = (u_star - dxcdt)*T_si_U(o)+ &
             (v_star - dycdt)*T_sj_U(o)+ &
             (w_star - dzcdt)*T_sk_U(o)


    R_IB_U=radius-rb_U(o)
    R_P_U=radius+dx1(1)


    uni=(1+beta_U)*unis-(beta_U)*uni_ref
    uti=(1+beta_U)*utis-(beta_U)*uti_ref
    usi=(1+beta_U)*usis-(beta_U)*usi_ref


! c     Solve u1
    u1_temp(i,j,k)=T_ni_U(o)*uni+ &
                  T_ti_U(o)*uti+ &
                  T_si_U(o)*usi+dxcdt

    end do
    u1=u1_temp

end subroutine bulb_u



subroutine finding_v_IB()
    use var2
    implicit none
    integer ::   o

    typ_v_IB(0:im+1,0:jm+1,0:km+1)=0
! c     if cell centre is inside of the bubble typ(i,j,k)=0
! c     if cell centre is outside of the bubble typ(i,j,k)=1
! c     if typ_IB(i,j,k) is 1 means it is IB cell
    do o=1,change-1
        i=x_IB(o)
        j=y_IB(o)
        k=z_IB(o)

        if (typ(i,j+1,k).eq.1) then
            typ_v_IB(i,j,k)=1
        else
            if(typ_IB(i,j+1,k).eq.1) then
                typ_v_IB(i,j,k)=1
            end if
        end if
        if (typ(i,j-1,k).eq.1) then
            typ_v_IB(i,j-1,k)=1
        else
            if(typ_IB(i,j-1,k).eq.1) then
                typ_v_IB(i,j-1,k)=1
            end if
        end if
    end do

    o_v=0
    do k=0,km+1
        do j=0,jm+1
            do i=0,im+1
                if (typ_v_IB(i,j,k).eq.1) then
                    o_v=o_v+1
                    x_IB_V(o_v)=i
                    y_IB_V(o_v)=j
                    z_IB_V(o_v)=k
                end if
            end do
        end do
    end do
end subroutine finding_v_IB




subroutine probe_IB_V()
    use var2
    implicit none
    integer ::   o,mm
    real(8) ::  rnode_V

! c        Neumann BC in Spherical  coordinates

! c        phiS = angle phi in spherical coordinates.
! c        not to be confused with phi the potential
! c        if x=0 then atan(y,x)=pi/2. This is not true however,
! c        if y is negative, then it should be 3pi/2.

    mm=1
    
    do o=1,o_v
       i=x_IB_V(o)
       j=y_IB_V(o)
       k=z_IB_V(o)

       phiS_V(o)=datan2((x2_j(j)-yc),(x1(i)-xc))

       thetaS_V(o)=acos((x3(k)-zc)/sqrt((x1(i)-xc)**2d0+ (x2_j(j)-yc)**2d0+(x3(k)-zc)**2d0))


       T_ni_V(o)=cos(phiS_V(o))*sin(thetaS_V(o))
       T_nj_V(o)=sin(phiS_V(o))*sin(thetaS_V(o))
       T_nk_V(o)=cos(thetaS_V(o))

       T_ti_V(o)=-sin(phiS_V(o))
       T_tj_V(o)=cos(phiS_V(o))
       T_tk_V(o)=0d0

       T_si_V(o)=cos(phiS_V(o))*cos(thetaS_V(o))
       T_sj_V(o)=sin(phiS_V(o))*cos(thetaS_V(o))
       T_sk_V(o)=-sin(thetaS_V(o))

! c        initialise x-direction
! c        Normal distance between interface and bubble centre
       rnode_V=sqrt((x1(i)-xc)**2d0+(x2_j(j)-yc)**2d0+(x3(k)-zc)**2d0)


! c        normal distance of interface to bubble surface 
! c        (can be negative)
       rb_V(o)=radius-rnode_V

       xb_V(o)=x1(i)  +rb_V(o)*T_ni_V(o)
       yb_V(o)=x2_j(j)+rb_V(o)*T_nj_V(o)
       zb_V(o)=x3(k)  +rb_V(o)*T_nk_V(o)

       x_prob_V(o,mm)=xb_V(o)+dx1_i(0)*T_ni_V(o)
       y_prob_V(o,mm)=yb_V(o)+dx1_i(0)*T_nj_V(o)
       z_prob_V(o,mm)=zb_V(o)+dx1_i(0)*T_nk_V(o)

! c        prob number of cells starts at an interface so only 
! c        +0.5d if that direction is on the level of a node
       nx_prob_VU(o,mm)=floor(x_prob_V(o,mm)/dx1_i(0))
       ny_prob_VU(o,mm)=floor(y_prob_V(o,mm)/dx2_j(0)+0.5d0)
       nz_prob_VU(o,mm)=floor(z_prob_V(o,mm)/dx3_k(0)+0.5d0)

        ! Check whether the probe is within the domain
       call chckProbeBounds(nx_prob_VU(o,mm),0,im+1,'x-VU')
       call chckProbeBounds(ny_prob_VU(o,mm),0,jm+1,'y-VU')
       call chckProbeBounds(nz_prob_VU(o,mm),0,km+1,'z-VU')

       nx_prob_VV(o,mm)=floor(x_prob_V(o,mm)/dx1_i(0)+0.5d0)
       ny_prob_VV(o,mm)=floor(y_prob_V(o,mm)/dx2_j(0))
       nz_prob_VV(o,mm)=floor(z_prob_V(o,mm)/dx3_k(0)+0.5d0)

        ! Check whether the probe is within the domain
       call chckProbeBounds(nx_prob_VV(o,mm),0,im+1,'x-VV')
       call chckProbeBounds(ny_prob_VV(o,mm),0,jm+1,'y-VV')
       call chckProbeBounds(nz_prob_VV(o,mm),0,km+1,'z-VV')

       nx_prob_VW(o,mm)=floor(x_prob_V(o,mm)/dx1_i(0)+0.5d0)
       ny_prob_VW(o,mm)=floor(y_prob_V(o,mm)/dx2_j(0)+0.5d0)
       nz_prob_VW(o,mm)=floor(z_prob_V(o,mm)/dx3_k(0))

        ! Check whether the probe is within the domain
       call chckProbeBounds(nx_prob_VW(o,mm),0,im+1,'x-VW')
       call chckProbeBounds(ny_prob_VW(o,mm),0,jm+1,'y-VW')
       call chckProbeBounds(nz_prob_VW(o,mm),0,km+1,'z-VW')

    end do
end subroutine probe_IB_V




subroutine bulb_v()
    use var2
    implicit none
    integer ::   o,mm
! c     Solving u2 on IB interface
    real(8) :: u1j_ref,u2j_ref,u3j_ref
    real(8) :: unj_ref,utj_ref,usj_ref
    real(8) ::  xd,yd,zd
    real(8) ::  c00,c10,c01,c11,c0,c1
    real(8) :: u_star,v_star,w_star
    real(8) :: unjs,utjs,usjs
    real(8) :: unj,utj,usj,beta_V
    real(8) :: R_IB_V,R_P_V
    ! real(8) :: u2_temp(0:im+1,0:jm+1,0:km+1)

    u2_temp=u2
    mm=1
    do o=1,o_v
        i=x_IB_V(o)
        j=y_IB_V(o)
        k=z_IB_V(o)

! c     x-component u2

    xd=(x_prob_V(o,mm)-x1_i(nx_prob_VU(o,mm)))/(dx1(1))
    yd=(y_prob_V(o,mm)-x2(ny_prob_VU(o,mm)))/(dx2_j(1))
    zd=(z_prob_V(o,mm)-x3(nz_prob_VU(o,mm)))/(dx3_k(1))


    c00=(1d0-xd)* &
     u1(nx_prob_VU(o,mm),ny_prob_VU(o,mm),nz_prob_VU(o,mm))+ &
     xd*u1(nx_prob_VU(o,mm)+1,ny_prob_VU(o,mm),nz_prob_VU(o,mm))

    c01=(1d0-xd)* &
     u1(nx_prob_VU(o,mm),ny_prob_VU(o,mm),nz_prob_VU(o,mm)+1)+ &
     xd*u1(nx_prob_VU(o,mm)+1,ny_prob_VU(o,mm),nz_prob_VU(o,mm)+1)

    c10=(1d0-xd)*&
     u1(nx_prob_VU(o,mm),ny_prob_VU(o,mm)+1,nz_prob_VU(o,mm))+ &
     xd*u1(nx_prob_VU(o,mm)+1,ny_prob_VU(o,mm)+1,nz_prob_VU(o,mm))

    c11=(1d0-xd)* &
     u1(nx_prob_VU(o,mm),ny_prob_VU(o,mm)+1,nz_prob_VU(o,mm)+1)+ &
     xd*u1(nx_prob_VU(o,mm)+1,ny_prob_VU(o,mm)+1,nz_prob_VU(o,mm)+1)

    c0=c00*(1d0-yd)+c10*yd
    c1=c01*(1d0-yd)+c11*yd

    u1j_ref=c0*(1d0-zd)+c1*zd


! c     y-component u2


    xd=(x_prob_V(o,mm)-x1(nx_prob_VV(o,mm)))/(dx1_i(1))
    yd=(y_prob_V(o,mm)-x2_j(ny_prob_VV(o,mm)))/(dx2(1))
    zd=(z_prob_V(o,mm)-x3(nz_prob_VV(o,mm)))/(dx3_k(1))

    c00=(1d0-xd)* &
     u2(nx_prob_VV(o,mm),ny_prob_VV(o,mm),nz_prob_VV(o,mm))+ &
     xd*u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm),nz_prob_VV(o,mm))

    c01=(1d0-xd)* &
     u2(nx_prob_VV(o,mm),ny_prob_VV(o,mm),nz_prob_VV(o,mm)+1)+ &
     xd*u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm),nz_prob_VV(o,mm)+1)

    c10=(1d0-xd)* &
     u2(nx_prob_VV(o,mm),ny_prob_VV(o,mm)+1,nz_prob_VV(o,mm))+ &
     xd*u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm)+1,nz_prob_VV(o,mm))

    c11=(1d0-xd)* &
     u2(nx_prob_VV(o,mm),ny_prob_VV(o,mm)+1,nz_prob_VV(o,mm)+1)+ &
     xd*u2(nx_prob_VV(o,mm)+1,ny_prob_VV(o,mm)+1,nz_prob_VV(o,mm)+1)

    c0=c00*(1d0-yd)+c10*yd
    c1=c01*(1d0-yd)+c11*yd

    u2j_ref=c0*(1d0-zd)+c1*zd

! c      z-component u2


    xd=(x_prob_V(o,mm)-x1(nx_prob_VW(o,mm)))/(dx1_i(1))
    yd=(y_prob_V(o,mm)-x2(ny_prob_VW(o,mm)))/(dx2_j(1))
    zd=(z_prob_V(o,mm)-x3_k(nz_prob_VW(o,mm)))/(dx3(1))

    c00=(1d0-xd)*&
     u3(nx_prob_VW(o,mm),ny_prob_VW(o,mm),nz_prob_VW(o,mm))+&
     xd*u3(nx_prob_VW(o,mm)+1,ny_prob_VW(o,mm),nz_prob_VW(o,mm))

    c01=(1d0-xd)*&
     u3(nx_prob_VW(o,mm),ny_prob_VW(o,mm),nz_prob_VW(o,mm)+1)+&
     xd*u3(nx_prob_VW(o,mm)+1,ny_prob_VW(o,mm),nz_prob_VW(o,mm)+1)

    c10=(1d0-xd)*&
     u3(nx_prob_VW(o,mm),ny_prob_VW(o,mm)+1,nz_prob_VW(o,mm))+&
     xd*u3(nx_prob_VW(o,mm)+1,ny_prob_VW(o,mm)+1,nz_prob_VW(o,mm))

    c11=(1d0-xd)*&
     u3(nx_prob_VW(o,mm),ny_prob_VW(o,mm)+1,nz_prob_VW(o,mm)+1)+ &
     xd*u3(nx_prob_VW(o,mm)+1,ny_prob_VW(o,mm)+1,nz_prob_VW(o,mm)+1)

    c0=c00*(1d0-yd)+c10*yd
    c1=c01*(1d0-yd)+c11*yd

    u3j_ref=c0*(1d0-zd)+c1*zd

      unj_ref = (u1j_ref - dxcdt)*T_ni_V(o) + &
                (u2j_ref - dycdt)*T_nj_V(o) + &
                (u3j_ref - dzcdt)*T_nk_V(o)

      utj_ref = (u1j_ref - dxcdt)*T_ti_V(o) + &
                (u2j_ref - dycdt)*T_tj_V(o) + &
                (u3j_ref - dzcdt)*T_tk_V(o)

      usj_ref = (u1j_ref - dxcdt)*T_si_V(o) + &
                (u2j_ref - dycdt)*T_sj_V(o) + &
                (u3j_ref - dzcdt)*T_sk_V(o)

! c       Dirichlet BC for normal direction
      beta_V=rb_V(o)/dx2(j)
! c       Velocity at the bubble surface in cartisian coordinate
      u_star = dxcdt + drdt*T_ni_V(o)
      v_star = dycdt + drdt*T_nj_V(o)
      w_star = dzcdt + drdt*T_nk_V(o)

      unjs = (u_star - dxcdt)*T_ni_V(o) + &
             (v_star - dycdt)*T_nj_V(o) + &
             (w_star - dzcdt)*T_nk_V(o)

      utjs = (u_star - dxcdt)*T_ti_V(o) + &
             (v_star - dycdt)*T_tj_V(o) + &
             (w_star - dzcdt)*T_tk_V(o)

     usjs = (u_star - dxcdt)*T_si_V(o) + &
            (v_star - dycdt)*T_sj_V(o) + &
            (w_star - dzcdt)*T_sk_V(o)

      R_IB_V=radius-rb_V(o)
      R_P_V=radius+dx1(1)

      unj=(1+beta_V)*unjs-(beta_V)*unj_ref
      utj=(1+beta_V)*utjs-(beta_V)*utj_ref
      usj=(1+beta_V)*usjs-(beta_V)*usj_ref

! c       Solve u2


    u2_temp(i,j,k) = T_nj_V(o)*unj + &
                     T_tj_V(o)*utj + &
                     T_sj_V(o)*usj + &
                     dycdt
    end do

    u2=u2_temp
end subroutine bulb_v



subroutine finding_w_IB()
    use var2
    implicit none
    integer ::   o
    typ_w_IB(0:im+1,0:jm+1,0:km+1)=0

    do o=1,change-1
        i=x_IB(o)
        j=y_IB(o)
        k=z_IB(o)
        if (typ(i,j,k+1).eq.1) then
            typ_w_IB(i,j,k)=1
        else
            if (typ_IB(i,j,k+1).eq.1) then
                typ_w_IB(i,j,k)=1
            end if
        end if
        if (typ(i,j,k-1).eq.1) then
            typ_w_IB(i,j,k-1)=1
        else
            if (typ_IB(i,j,k-1).eq.1) then
                typ_w_IB(i,j,k-1)=1
            end if
        end if
    end do

    o_w=0
    do k=0,km+1
        do j=0,jm+1
            do i=0,im+1
                if (typ_w_IB(i,j,k).eq.1) then
                    o_w=o_w+1
                    x_IB_W(o_w)=i
                    y_IB_W(o_w)=j
                    z_IB_W(o_w)=k
                end if
            end do
        end do
    end do
end subroutine finding_w_IB




subroutine probe_IB_W()
    use var2
    implicit none
    integer ::   o,mm
    real(8) ::  rnode_W

! c        Neumann BC in Spherical  coordinates

! c        phiS = angle phi in spherical coordinates.
! c        not to be confused with phi the potential
! c        if x=0 then atan(y,x)=pi/2. This is not true however,
! c        if y is negative, then it should be 3pi/2.

    mm=1
    
    do o=1,o_w
       i=x_IB_W(o)
       j=y_IB_W(o)
       k=z_IB_W(o)

       phiS_W(o)=datan2((x2(j)-yc),(x1(i)-xc))

       thetaS_W(o)=acos((x3_k(k)-zc)/sqrt((x1(i)-xc)**2d0+(x2(j)-yc)**2d0+(x3_k(k)-zc)**2d0))


       T_ni_W(o)=cos(phiS_W(o))*sin(thetaS_W(o))
       T_nj_W(o)=sin(phiS_W(o))*sin(thetaS_W(o))
       T_nk_W(o)=cos(thetaS_W(o))

       T_ti_W(o)=-sin(phiS_W(o))
       T_tj_W(o)=cos(phiS_W(o))
       T_tk_W(o)=0d0

       T_si_W(o)=cos(phiS_W(o))*cos(thetaS_W(o))
       T_sj_W(o)=sin(phiS_W(o))*cos(thetaS_W(o))
       T_sk_W(o)=-sin(thetaS_W(o))

! c        initialise x-direction
! c        Normal distance between interface and bubble centre
       rnode_W=sqrt((x1(i)-xc)**2d0+(x2(j)-yc)**2d0+(x3_k(k)-zc)**2d0)


! c        normal distance of interface to bubble surface 
! c        (can be negative)
       rb_W(o)=radius-rnode_W

       xb_W(o)=x1(i)  +rb_W(o)*T_ni_W(o)
       yb_W(o)=x2(j)  +rb_W(o)*T_nj_W(o)
       zb_W(o)=x3_k(k)+rb_W(o)*T_nk_W(o)


       

       x_prob_W(o,mm)=xb_W(o)+dx1_i(0)*T_ni_W(o)
       y_prob_W(o,mm)=yb_W(o)+dx1_i(0)*T_nj_W(o)
       z_prob_W(o,mm)=zb_W(o)+dx1_i(0)*T_nk_W(o)

! c        prob number of cells starts at an interface so only 
! c        +0.5d if that direction is on the level of a node
       nx_prob_WU(o,mm)=floor(x_prob_W(o,mm)/dx1_i(0))
       ny_prob_WU(o,mm)=floor(y_prob_W(o,mm)/dx2_j(0)+0.5d0)
       nz_prob_WU(o,mm)=floor(z_prob_W(o,mm)/dx3_k(0)+0.5d0)

        ! Check whether the probe is within the domain
       call chckProbeBounds(nx_prob_WU(o,mm),0,im+1,'x-WU')
       call chckProbeBounds(ny_prob_WU(o,mm),0,jm+1,'y-WU')
       call chckProbeBounds(nz_prob_WU(o,mm),0,km+1,'z-WU')

       nx_prob_WV(o,mm)=floor(x_prob_W(o,mm)/dx1_i(0)+0.5d0)
       ny_prob_WV(o,mm)=floor(y_prob_W(o,mm)/dx2_j(0))
       nz_prob_WV(o,mm)=floor(z_prob_W(o,mm)/dx3_k(0)+0.5d0)

        ! Check whether the probe is within the domain
       call chckProbeBounds(nx_prob_WV(o,mm),0,im+1,'x-WV')
       call chckProbeBounds(ny_prob_WV(o,mm),0,jm+1,'y-WV')
       call chckProbeBounds(nz_prob_WV(o,mm),0,km+1,'z-WV')

       nx_prob_WW(o,mm)=floor(x_prob_W(o,mm)/dx1_i(0)+0.5d0)
       ny_prob_WW(o,mm)=floor(y_prob_W(o,mm)/dx2_j(0)+0.5d0)
       nz_prob_WW(o,mm)=floor(z_prob_W(o,mm)/dx3_k(0))

        ! Check whether the probe is within the domain
       call chckProbeBounds(nx_prob_WW(o,mm),0,im+1,'x-WW')
       call chckProbeBounds(ny_prob_WW(o,mm),0,jm+1,'y-WW')
       call chckProbeBounds(nz_prob_WW(o,mm),0,km+1,'z-WW')

    end do
end subroutine probe_IB_W




subroutine bulb_w()
    use var2
    implicit none
    integer ::   o,mm
! c     solving u3 IB interface
    real(8) ::  u1k_ref,u2k_ref,u3k_ref
    real(8) ::  unk_ref,utk_ref,usk_ref
    real(8) ::  xd,yd,zd
    real(8) ::  c00,c10,c01,c11,c0,c1
    real(8) ::  u_star,v_star,w_star
    real(8) ::  unks,utks,usks
    real(8) ::  unk,utk,usk,beta_W
    real(8) ::  R_IB_W,R_P_W
    ! real(8) :: u3_temp(0:im+1,0:jm+1,0:km+1)

    u3_temp=u3
    mm=1
    do o=1,o_w
        i=x_IB_W(o)
        j=y_IB_W(o)
        k=z_IB_W(o)

! c     x-component u3

    xd=(x_prob_W(o,mm)-x1_i(nx_prob_WU(o,mm)))/(dx1(1))
    yd=(y_prob_W(o,mm)-x2(ny_prob_WU(o,mm)))/(dx2_j(1))
    zd=(z_prob_W(o,mm)-x3(nz_prob_WU(o,mm)))/(dx3_k(1))


    c00=(1d0-xd)*&
      u1(nx_prob_WU(o,mm),ny_prob_WU(o,mm),nz_prob_WU(o,mm))+&
      xd*u1(nx_prob_WU(o,mm)+1,ny_prob_WU(o,mm),nz_prob_WU(o,mm))

    c01=(1d0-xd)*&
     u1(nx_prob_WU(o,mm),ny_prob_WU(o,mm),nz_prob_WU(o,mm)+1)+&
     xd*u1(nx_prob_WU(o,mm)+1,ny_prob_WU(o,mm),nz_prob_WU(o,mm)+1)

    c10=(1d0-xd)*&
     u1(nx_prob_WU(o,mm),ny_prob_WU(o,mm)+1,nz_prob_WU(o,mm))+&
     xd*u1(nx_prob_WU(o,mm)+1,ny_prob_WU(o,mm)+1,nz_prob_WU(o,mm))

    c11=(1d0-xd)*&
     u1(nx_prob_WU(o,mm),ny_prob_WU(o,mm)+1,nz_prob_WU(o,mm)+1)+&
     xd*u1(nx_prob_WU(o,mm)+1,ny_prob_WU(o,mm)+1,nz_prob_WU(o,mm)+1)

    c0=c00*(1d0-yd)+c10*yd
    c1=c01*(1d0-yd)+c11*yd

    u1k_ref=c0*(1d0-zd)+c1*zd


! c       y-component u3


    xd=(x_prob_W(o,mm)-x1(nx_prob_WV(o,mm)))/(dx1_i(1))
    yd=(y_prob_W(o,mm)-x2_j(ny_prob_WV(o,mm)))/(dx2(1))
    zd=(z_prob_W(o,mm)-x3(nz_prob_WV(o,mm)))/(dx3_k(1))

    c00=(1d0-xd)*&
     u2(nx_prob_WV(o,mm),ny_prob_WV(o,mm),nz_prob_WV(o,mm))+&
     xd*u2(nx_prob_WV(o,mm)+1,ny_prob_WV(o,mm),nz_prob_WV(o,mm))

    c01=(1d0-xd)*&
     u2(nx_prob_WV(o,mm),ny_prob_WV(o,mm),nz_prob_WV(o,mm)+1)+&
     xd*u2(nx_prob_WV(o,mm)+1,ny_prob_WV(o,mm),nz_prob_WV(o,mm)+1)

    c10=(1d0-xd)*&
     u2(nx_prob_WV(o,mm),ny_prob_WV(o,mm)+1,nz_prob_WV(o,mm))+&
     xd*u2(nx_prob_WV(o,mm)+1,ny_prob_WV(o,mm)+1,nz_prob_WV(o,mm))

    c11=(1d0-xd)*&
     u2(nx_prob_WV(o,mm),ny_prob_WV(o,mm)+1,nz_prob_WV(o,mm)+1)+&
     xd*u2(nx_prob_WV(o,mm)+1,ny_prob_WV(o,mm)+1,nz_prob_WV(o,mm)+1)

    c0=c00*(1d0-yd)+c10*yd
    c1=c01*(1d0-yd)+c11*yd

    u2k_ref=c0*(1d0-zd)+c1*zd

! c     z-component u3

    xd=(x_prob_W(o,mm)-x1(nx_prob_WW(o,mm)))/(dx1_i(1))
    yd=(y_prob_W(o,mm)-x2(ny_prob_WW(o,mm)))/(dx2_j(1))
    zd=(z_prob_W(o,mm)-x3_k(nz_prob_WW(o,mm)))/(dx3(1))

    c00=(1d0-xd)*&
     u3(nx_prob_WW(o,mm),ny_prob_WW(o,mm),nz_prob_WW(o,mm))+&
     xd*u3(nx_prob_WW(o,mm)+1,ny_prob_WW(o,mm),nz_prob_WW(o,mm))

    c01=(1d0-xd)*&
     u3(nx_prob_WW(o,mm),ny_prob_WW(o,mm),nz_prob_WW(o,mm)+1)+&
     xd*u3(nx_prob_WW(o,mm)+1,ny_prob_WW(o,mm),nz_prob_WW(o,mm)+1)

    c10=(1d0-xd)*&
     u3(nx_prob_WW(o,mm),ny_prob_WW(o,mm)+1,nz_prob_WW(o,mm))+&
     xd*u3(nx_prob_WW(o,mm)+1,ny_prob_WW(o,mm)+1,nz_prob_WW(o,mm))

    c11=(1d0-xd)*&
     u3(nx_prob_WW(o,mm),ny_prob_WW(o,mm)+1,nz_prob_WW(o,mm)+1)+&
     xd*u3(nx_prob_WW(o,mm)+1,ny_prob_WW(o,mm)+1,nz_prob_WW(o,mm)+1)

    c0=c00*(1d0-yd)+c10*yd
    c1=c01*(1d0-yd)+c11*yd

    u3k_ref=c0*(1d0-zd)+c1*zd

    unk_ref = (u1k_ref - dxcdt)*T_ni_W(o) + &
              (u2k_ref - dycdt)*T_nj_W(o) + &
              (u3k_ref - dzcdt)*T_nk_W(o)

    utk_ref = (u1k_ref - dxcdt)*T_ti_W(o) + &
              (u2k_ref - dycdt)*T_tj_W(o) + &
              (u3k_ref - dzcdt)*T_tk_W(o)

    usk_ref = (u1k_ref - dxcdt)*T_si_W(o) + &
              (u2k_ref - dycdt)*T_sj_W(o) + &
              (u3k_ref - dzcdt)*T_sk_W(o)


! c     Dirichlet BC for normal direction
    beta_W=rb_W(o)/dx3(k)
    u_star = dxcdt + drdt*T_ni_W(o)
    v_star = dycdt + drdt*T_nj_W(o)
    w_star = dzcdt + drdt*T_nk_W(o)

    unks = (u_star - dxcdt)*T_ni_W(o) + &
           (v_star - dycdt)*T_nj_W(o) + &
           (w_star - dzcdt)*T_nk_W(o)

    utks = (u_star - dxcdt)*T_ti_W(o) + &
           (v_star - dycdt)*T_tj_W(o) + &
           (w_star - dzcdt)*T_tk_W(o)

    usks = (u_star - dxcdt)*T_si_W(o) + &
           (v_star - dycdt)*T_sj_W(o) + &
           (w_star - dzcdt)*T_sk_W(o)


    R_IB_W=radius-rb_W(o)
    R_P_W=radius+dx1(1)

    unk=(1+beta_W)*unks-(beta_W)*unk_ref
    utk=(1+beta_W)*utks-(beta_W)*utk_ref
    usk=(1+beta_W)*usks-(beta_W)*usk_ref

! c     Solve u3
    u3_temp(i,j,k) = T_nk_W(o)*unk + &
                     T_tk_W(o)*utk + &
                     T_sk_W(o)*usk + &
                     dzcdt

    end do

    u3=u3_temp
end subroutine bulb_w



subroutine bulb_force_x()
    use var2
    implicit none
    integer :: o
    real(8) :: term1,term2,term3,total_up,total_bellow
    ! integer :: typ_force_x(0:im+1,0:jm+1,0:km+1)

    total=0d0
    total_bellow=0d0
    total_up=0d0
    typ_force_x(0:im+1,0:jm+1,0:km+1)=0
    FxCOM = 0d0
    term1 = 0d0
    term2 = 0d0
    term3 = 0d0

    do k=0,km+1
    do j=0,jm+1
    do i=0,im+1

       if ((x1_i(i)-xc)**2d0+(x2(j)-yc)**2d0+(x3(k)-zc)**2d0.gt.radius**2d0) then
          typ_force_x(i,j,k)=1
       end if
    end do
    end do
    end do

    do o=1,o_u
       i=x_IB_u(o)
       j=y_IB_u(o)
       k=z_IB_u(o)
       typ_force_x(i,j,k)=0
    end do


    do o=1,o_u
       i=x_IB_u(o)
       j=y_IB_u(o)
       k=z_IB_u(o)

       if(typ_force_x(i+1,j,k).eq.1) then
        !  term1=-(-(0.5d0*(u1(i+1,j,k)+u1(i,j,k)))**2d0*dx2(j)*dx3(k))
         term2=-(-p(i+1,j,k)/rhoKOH*dx2(j)*dx3(k))
         term3=-(nuKOH*(u1(i+1,j,k)-u1(i,j,k))/dx1(i+1)*dx2(j)*dx3(k))    
         total=total+term1+term2+term3
         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if

       if(typ_force_x(i-1,j,k).eq.1) then
        !  term1=-(0.5d0*(u1(i-1,j,k)+u1(i,j,k)))**2d0*dx2(j)*dx3(k)
         term2=-p(i,j,k)/rhoKOH*dx2(j)*dx3(k)
         term3=nuKOH*(u1(i,j,k)-u1(i-1,j,k))/dx1(i)*dx2(j)*dx3(k)
         total=total+term1+term2+term3
         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if

       if(typ_force_x(i,j+1,k).eq.1) then
        !  term1=-(-0.5d0*(u1(i,j,k)+u1(i,j+1,k))&
        ! *0.5d0*(u2(i,j,k)+u2(i+1,j,k))*dx1(i)*dx3(k))
         term2=0d0 
         term3=-(nuKOH*(u1(i,j+1,k)-u1(i,j,k))/dx2_j(j)*dx1(i)*dx3(k)) 
         total=total+term1+term2+term3
         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if

       if(typ_force_x(i,j-1,k).eq.1) then
        !  term1=-0.5d0*(u1(i,j,k)+u1(i,j-1,k)) &
        ! *0.5d0*(u2(i,j-1,k)+u2(i+1,j-1,k))*dx1(i)*dx3(k)
         term2=0d0
         term3=nuKOH*(u1(i,j,k)-u1(i,j-1,k))/dx2_j(j-1)*dx1(i)*dx3(k)    
         total=total+term1+term2+term3
         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if

       if(typ_force_x(i,j,k+1).eq.1) then
    !      term1=-(-0.5d0*(u1(i,j,k)+u1(i,j,k+1)) &
    !    *0.5d0*(u3(i,j,k)+u3(i+1,j,k))*dx1(i)*dx2(j))
         term2=0d0    
         term3=-(nuKOH*(u1(i,j,k+1)-u1(i,j,k))/dx3_k(k)*dx1(i)*dx2(j))   
         total=total+term1+term2+term3
         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if

       if(typ_force_x(i,j,k-1).eq.1) then
        !  term1=-0.5d0*(u1(i,j,k)+u1(i,j,k-1)) &
        ! *0.5d0*(u3(i,j,k-1)+u3(i+1,j,k-1))*dx1(i)*dx2(j)
         term2=0d0    
         term3=nuKOH*(u1(i,j,k)-u1(i,j,k-1))/dx3_k(k-1)*dx1(i)*dx2(j)            
         total=total+term1+term2+term3
         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if
    end do

    FxCOM = total

    if ( DragForceOutput ) then
1       format(4E36.24)
        if ((m1.eq.1).and.(m2.eq.1)) then
            open(unit=1,file='./output/FxBubble_top_bottom.dat')
2          format(4A36)
            write(1,2) 'time', 'total_up', 'total_bellow', 'FxCOM'
            write(1,1) time, total_up, total_bellow, FxCOM
            close(1,status='keep')
        end if
        if ( (m1.gt.1).and.(m2.eq.m2_max) ) then
            open(unit=1,file='./output/FxBubble_top_bottom.dat',position='append')
            write(1,1) time, total_up, total_bellow, FxCOM
            close(1,status='keep')
        end if
    end if


end subroutine bulb_force_x




subroutine domain_flux_force_x()
    use var2   
    implicit none   
    real(8) :: term1,term2,term3,term4,term5,term6, &
        !   A11(1:jm,1:km),B11(1:jm,1:km), &
        !   A12(1:jm,1:km),B12(1:jm,1:km), &
        !   A13(1:jm,1:km),B13(1:jm,1:km), &
        !   A21(1:im,1:km),B21(1:im,1:km), &
        !   A22(1:im,1:km),B22(1:im,1:km), &
        !   A23(1:im,1:km),B23(1:im,1:km), &
        !   A31(1:im,1:jm),B31(1:im,1:jm), &
        !   A32(1:im,1:jm),B32(1:im,1:jm), &
        !   A33(1:im,1:jm),B33(1:im,1:jm),
          Error

    a11=0d0
    b11=0d0
    a12=0d0
    b12=0d0
    a13=0d0
    b13=0d0

    a21=0d0
    b21=0d0
    a22=0d0
    b22=0d0
    a23=0d0
    b23=0d0

    a31=0d0
    b31=0d0
    a32=0d0
    b32=0d0
    a33=0d0
    b33=0d0


    do k=1,km
    do j=1,jm

       A11(j,k)=-(-(0.5d0*(u1(0,j,k)+u1(1,j,k)))**2d0&
                *dx2(j)*dx3(k))
       B11(j,k)=-(0.5d0*(u1(im,j,k)+u1(im-1,j,k)))**2d0&
                *dx2(j)*dx3(k) 

       A12(j,k)=-(-p(1,j,k)*dx2(j)*dx3(k)/rhoKOH)
       B12(j,k)=-p(im,j,k)*dx2(j)*dx3(k)/rhoKOH

       A13(j,k)=-(nuKOH*(u1(1,j,k)-u1(0,j,k)) &
               /dx1(1)*dx2(j)*dx3(k))
       B13(j,k)=nuKOH*(u1(im,j,k)-u1(im-1,j,k)) &
               /dx1(1)*dx2(j)*dx3(k)
               !!          ^ was first 'i', but that resulted in an error
    end do
    end do

    term1=sum(A11)+sum(A12)+sum(A13)
    term2=sum(B11)+sum(B12)+sum(B13)

    do k=1,km
    do i=1,im-1
       A21(i,k)=-(-0.5d0*(u1(i,0,k)+u1(i,1,k)) &
      *0.5d0*(u2(i,0,k)+u2(i+1,0,k))*dx1(i)*dx3(k))
       B21(i,k)=-0.5d0*(u1(i,jm,k)+u1(i,jm+1,k)) &
      *0.5d0*(u2(i,jm,k)+u2(i+1,jm,k))*dx1(i)*dx3(k)

       A22(i,k)=0d0
       B22(i,k)=0d0

       A23(i,k)=-nuKOH*(u1(i,1,k)-u1(i,0,k))/dx2_j(0)*dx1(i)*dx3(k) 
       B23(i,k)=0d0
    end do
    end do

    term3=sum(A21)+sum(A22)+sum(A23)
    term4=sum(B21)+sum(B22)+sum(B23)

    do j=1,jm
    do i=1,im-1
       A31(i,j)=-(-0.5d0*(u1(i,j,0)+u1(i,j,1))&
      *0.5d0*(u3(i,j,0)+u3(i+1,j,0))*dx1(i)*dx2(j))
       B31(i,j)=-0.5d0*(u1(i,j,km)+u1(i,j,km+1))&
        *0.5d0*(u3(i,j,km)+u3(i+1,j,km))*dx1(i)*dx2(j)

       A32(i,j)=0d0
       B32(i,j)=0d0

       A33(i,j)=-nuKOH*(u1(i,j,1)-u1(i,j,0))*dx1(i)*dx2(j)/dx3_k(0)            
       B33(i,j)=nuKOH*(u1(i,j,km+1)-u1(i,j,km))*dx1(i)*dx2(j)/dx3_k(km)

    end do
    end do

    term5=sum(A31)+sum(A32)+sum(A33)
    term6=sum(B31)+sum(B32)+sum(B33)

    Error=abs(term1+term2+term3+term4+term5+term6+total)

if ( .false. ) then
1   format(9E36.24)
    if ((m1.eq.1).and.(m2.eq.1)) then
        open(unit=1,file='./output/FxDomain.dat')
2       format(9A36)
        write(1,2) 'time', 'term1', 'term2', 'term3', 'term4',  &
                   'term5', 'term6', 'total', 'Error'
        write(1,1) time, term1, term2, term3, term4, term5,     &
                    term6, total, Error
        close(1,status='keep')
    end if
    if ( (m1.gt.1).and.(m2.eq.m2_max) ) then
        open(unit=1,file='./output/FxDomain.dat',position='append')
        write(1,1) time, term1, term2, term3, term4, term5,     &
                    term6, total, Error
        close(1,status='keep')
    end if
end if


end subroutine domain_flux_force_x




subroutine bulb_force_y()
    !!------------------------------------------------------------------------------
    !! Compute the net acceleration in y-direction
    !!------------------------------------------------------------------------------
    !! DESCRIPTION:
    !!> Computes the net acceleration by evaluating the discrete surface integral 
    !!> over the bubble surface and stores it in FyCOM.
    !!
    !!> Unit is [m4/s2], as it is divided by the fluid density.
    !!
    !!
    use var2  
    implicit none    
    integer :: o
    real(8) :: term1,term2,term3,total_up,total_bellow
    ! integer :: typ_force_y(0:im+1,0:jm+1,0:km+1)

    total=0d0
    total_up=0d0
    total_bellow=0d0
    typ_force_y(0:im+1,0:jm+1,0:km+1)=0
    FyCOM = 0d0
    term1 = 0d0
    term2 = 0d0
    term3 = 0d0

    do k=0,km+1
    do j=0,jm+1
    do i=0,im+1


      if ((x1(i)-xc)**2d0+(x2_j(j)-yc)**2d0 &
      +(x3(k)-zc)**2d0.gt.radius**2d0) then
          typ_force_y(i,j,k)=1
      end if
    end do
    end do
    end do


    do o=1,o_v
       i=x_IB_V(o)
       j=y_IB_V(o)
       k=z_IB_V(o)

       typ_force_y(i,j,k)=0

    end do

    do o=1,o_v
       i=x_IB_V(o)
       j=y_IB_V(o)
       k=z_IB_V(o)


       if(typ_force_y(i+1,j,k).eq.1) then

        !  term1=-(-0.5d0*(u2(i,j,k)+u2(i+1,j,k))*&
        ! 0.5d0*(u1(i,j,k)+u1(i,j+1,k))*dx2(j)*dx3(k))
         term2=0d0    
         term3=-(nuKOH*(u2(i+1,j,k)-u2(i,j,k))/dx1_i(i)*dx2(j)*dx3(k))          
         total=total+term1+term2+term3

         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if

       if(typ_force_y(i-1,j,k).eq.1) then

        !  term1=-0.5d0*(u2(i,j,k)+u2(i-1,j,k))*&
        !  0.5d0*(u1(i-1,j,k)+u1(i-1,j+1,k))*dx2(j)*dx3(k)
         term2=0d0
         term3=nuKOH*(u2(i,j,k)-u2(i-1,j,k))/dx1_i(i-1)*dx2(j)*dx3(k)  
         total=total+term1+term2+term3

         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if

              
       end if

       if(typ_force_y(i,j+1,k).eq.1) then

        !  term1=-(-(0.5d0*(u2(i,j,k)+u2(i,j+1,k)))**2d0*dx1(i)*dx3(k))
         term2=-(-p(i,j+1,k)/rhoKOH*dx1(i)*dx3(k))    
         term3=-(nuKOH*(u2(i,j+1,k)-u2(i,j,k))/dx2_j(j)*dx1(i)*dx3(k))       
         total=total+term1+term2+term3

         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if            
       end if

       if(typ_force_y(i,j-1,k).eq.1) then

        !  term1=-(0.5d0*(u2(i,j,k)+u2(i,j-1,k)))**2d0*dx1(i)*dx3(k)
         term2=-p(i,j,k)/rhoKOH*dx1(i)*dx3(k)
         term3=nuKOH*(u2(i,j,k)-u2(i,j-1,k))/dx2_j(j-1)*dx1(i)*dx3(k)     
         total=total+term1+term2+term3

         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if

       if(typ_force_y(i,j,k+1).eq.1) then

        !  term1=-(-0.5d0*(u2(i,j,k)+u2(i,j,k+1))&
        !  *0.5d0*(u3(i,j,k)+u3(i,j+1,k))*dx1(i)*dx2(j))
         term2=0d0    
         term3=-(nuKOH*(u2(i,j,k+1)-u2(i,j,k))/dx3_k(k)*dx1(i)*dx2(j))   
         total=total+term1+term2+term3

         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if

       if(typ_force_y(i,j,k-1).eq.1) then
        !  term1=-0.5d0*(u2(i,j,k)+u2(i,j,k-1))*&
        ! 0.5d0*(u3(i,j,k-1)+u3(i,j+1,k-1))*dx1(i)*dx2(j)
         term2=0d0    
         term3=nuKOH*(u2(i,j,k)-u2(i,j,k-1))/dx3_k(k-1)*dx1(i)*dx2(j)               
         total=total+term1+term2+term3
         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if
    end do  

    FyCOM = total

    if ( DragForceOutput ) then
1       format(4E36.24)
        if ((m1.eq.1).and.(m2.eq.1)) then
            open(unit=1,file='./output/FyBubble_top_bottom.dat')
2       format(4A36)
            write(1,2) 'time', 'total_up', 'total_bellow', 'FyCOM'
            write(1,1) time, total_up, total_bellow, FyCOM
            close(1,status='keep')
        end if
        if ( (m1.gt.1).and.(m2.eq.m2_max) ) then
            open(unit=1,file='./output/FyBubble_top_bottom.dat',position='append')
            write(1,1) time, total_up, total_bellow, FyCOM
            close(1,status='keep')
        end if
    end if
    
end subroutine bulb_force_y




subroutine domain_flux_force_y()
    use var2
    implicit none
    real(8) :: term1,term2,term3,term4,term5,term6, &
        !   A11(1:jm,1:km),B11(1:jm,1:km), &
        !   A12(1:jm,1:km),B12(1:jm,1:km), &
        !   A13(1:jm,1:km),B13(1:jm,1:km), &
        !   A21(1:im,1:km),B21(1:im,1:km), &
        !   A22(1:im,1:km),B22(1:im,1:km), &
        !   A23(1:im,1:km),B23(1:im,1:km), &
        !   A31(1:im,1:jm),B31(1:im,1:jm), &
        !   A32(1:im,1:jm),B32(1:im,1:jm), &
        !   A33(1:im,1:jm),B33(1:im,1:jm),
          Error


    a11=0d0
    b11=0d0
    a12=0d0
    b12=0d0
    a13=0d0
    b13=0d0

    a21=0d0
    b21=0d0
    a22=0d0
    b22=0d0
    a23=0d0
    b23=0d0

    a31=0d0
    b31=0d0
    a32=0d0
    b32=0d0
    a33=0d0
    b33=0d0

    do k=1,km
    do j=1,jm-1
       A11(j,k)=0d0
       B11(j,k)=-0.5d0*(u2(im,j,k)+u2(im+1,j,k)) &
       *0.5d0*(u1(im,j,k)+u1(im,j+1,k))*dx2(j)*dx3(k)
       A12(j,k)=0d0
       B12(j,k)=0d0
       A13(j,k)=(-nuKOH*(u2(1,j,k)-u2(0,j,k))*&
               dx2(j)*dx3(k)/dx1_i(0))
       B13(j,k)=nuKOH*(u2(im+1,j,k)-u2(im,j,k))*&
              dx2(j)*dx3(k)/dx1_i(im)
    end do
    end do

    term1=sum(A11)+sum(A12)+sum(A13) 
    term2=sum(B11)+sum(B12)+sum(B13)

    do k=1,km
    do i=1,im
      A21(i,k)=-(-(0.5d0*(u2(i,0,k)+u2(i,1,k)))**2d0 &
              *dx1(i)*dx3(k))
      B21(i,k)=-(0.5d0*(u2(i,jm,k)+u2(i,jm-1,k)))**2d0 &
              *dx1(i)*dx3(k) 

      A22(i,k)=-(-(p(i,1,k))*dx1(i)*dx3(k)/rhoKOH)
      B22(i,k)=-p(i,jm,k)*dx1(i)*dx3(k)/rhoKOH

      A23(i,k)=-(nuKOH*(u2(i,1,k)-u2(i,0,k))/dx2(1)*dx1(i)*dx3(k))
      B23(i,k)=0d0

    end do
    end do

    term3=sum(A21)+sum(A22)+sum(A23)
    term4=sum(B21)+sum(B22)+sum(B23)


    do j=1,jm
    do i=1,im

       A31(i,j)=-(-0.5d0*(u2(i,j,0)+u2(i,j,1)) &
               *0.5d0*(u3(i,j,0)+u3(i,j+1,0))*dx1(i)*dx2(j))
       B31(i,j)=-0.5d0*(u2(i,j,km)+u2(i,j,km+1)) &
               *0.5d0*(u3(i,j,km)+u3(i,j+1,km))*dx1(i)*dx2(j)

       A32(i,j)=0d0
       B32(i,j)=0d0

       A33(i,j)=-(nuKOH*(u2(i,j,1)-u2(i,j,0)) &
               *dx1(i)*dx2(j)/dx3_k(0))
       B33(i,j)=nuKOH*(u2(i,j,km+1)-u2(i,j,km)) &
               *dx1(i)*dx2(j)/dx3_k(km)

    end do
    end do

    term5=sum(A31)+sum(A32)+sum(A33)
    term6=sum(B31)+sum(B32)+sum(B33)
    
    Error=abs(term1+term2+term3+term4+term5+term6+total)


if ( .false. ) then
1   format(9E36.24)
    if ((m1.eq.1).and.(m2.eq.1)) then
        open(unit=1,file='./output/FyDomain.dat')
2       format(9A36)
        write(1,2) 'time', 'term1', 'term2', 'term3', 'term4',  &
                   'term5', 'term6', 'total', 'Error'
        write(1,1) time, term1, term2, term3, term4, term5,     &
                    term6, total, Error
        close(1,status='keep')
    end if
    if ( (m1.gt.1).and.(m2.eq.m2_max) ) then
        open(unit=1,file='./output/FyDomain.dat',position='append')
        write(1,1) time, term1, term2, term3, term4, term5,     &
                    term6, total, Error
        close(1,status='keep')
    end if
end if


end subroutine domain_flux_force_y




subroutine bulb_force_z()
    use var2
    implicit none
    integer :: o
    real(8) :: term1,term2,term3,total_up,total_bellow
    ! integer :: typ_force_z(0:im+1,0:jm+1,0:km+1)

    total = 0d0
    total_up = 0d0
    total_bellow = 0d0
    typ_force_z(0:im+1,0:jm+1,0:km+1) = 0
    FzCOM = 0d0
    term1 = 0d0
    term2 = 0d0
    term3 = 0d0

    do k=0,km+1
        do j=0,jm+1
            do i=0,im+1
                if ((x1(i) - xc)**2d0 + (x2(j) - yc)**2d0 + (x3_k(k) - zc)**2d0.gt.radius**2d0) then
                    typ_force_z(i,j,k)=1
                end if
            end do
        end do
    end do

    do o=1,o_w
       i=x_IB_W(o)
       j=y_IB_W(o)
       k=z_IB_W(o)
       typ_force_z(i,j,k)=0
    end do

    do o=1,o_w
       i=x_IB_W(o)
       j=y_IB_W(o)
       k=z_IB_W(o)

        if(typ_force_z(i+1,j,k).eq.1) then
            !  term1=-(-0.5d0*(u2(i,j,k)+u2(i+1,j,k))*&
            ! 0.5d0*(u1(i,j,k)+u1(i,j+1,k))*dx2(j)*dx3(k))
            term2 = 0d0    
            term3 = -(nuKOH*(u3(i+1,j,k)-u3(i,j,k))/dx1_i(i)*dx2(j)*dx3(k))          
            total = total + term1 + term2 + term3
            if (x2_j(j).gt.yc) then
                total_up = total_up + term1 + term2 + term3
            else
                total_bellow = total_bellow + term1 + term2 + term3
            end if
        end if

        if(typ_force_z(i-1,j,k).eq.1) then
            !  term1=-0.5d0*(u2(i,j,k)+u2(i-1,j,k))*&
            !  0.5d0*(u1(i-1,j,k)+u1(i-1,j+1,k))*dx2(j)*dx3(k)
            term2 = 0d0
            term3 = nuKOH*(u3(i,j,k) - u3(i-1,j,k))/dx1_i(i-1)*dx2(j)*dx3(k)  
            total = total + term1 + term2 + term3
            if (x2_j(j).gt.yc) then
                total_up=total_up+term1+term2+term3
            else
                total_bellow=total_bellow+term1+term2+term3
            end if              
        end if

        if(typ_force_z(i,j+1,k).eq.1) then
            !  term1=-(-(0.5d0*(u2(i,j,k)+u2(i,j+1,k)))**2d0*dx1(i)*dx3(k))
            term2 = 0d0
            term3=-(nuKOH*(u3(i,j+1,k) - u3(i,j,k))/dx2_j(j)*dx1(i)*dx3(k))       
            total=total+term1+term2+term3
            if (x2_j(j).gt.yc) then
                total_up=total_up+term1+term2+term3
            else
                total_bellow=total_bellow+term1+term2+term3
            end if            
        end if

       if(typ_force_z(i,j-1,k).eq.1) then

        !  term1=-(0.5d0*(u2(i,j,k)+u2(i,j-1,k)))**2d0*dx1(i)*dx3(k)
         term2= 0d0
         term3=nuKOH*(u3(i,j,k) - u3(i,j-1,k))/dx2_j(j-1)*dx1(i)*dx3(k)     
         total=total+term1+term2+term3

         if (x2_j(j).gt.yc) then
            total_up=total_up+term1+term2+term3
         else
            total_bellow=total_bellow+term1+term2+term3
         end if
       end if

       if(typ_force_z(i,j,k+1).eq.1) then
            !  term1=-(-0.5d0*(u2(i,j,k)+u2(i,j,k+1))&
            !  *0.5d0*(u3(i,j,k)+u3(i,j+1,k))*dx1(i)*dx2(j))
            term2= -(-p(i,j,k+1)/rhoKOH*dx1(i)*dx2(j)) 
            term3= -(nuKOH*(u3(i,j,k+1) - u3(i,j,k))/dx3(k+1)*dx1(i)*dx2(j))   
            total=total+term1+term2+term3

            if (x2_j(j).gt.yc) then
                total_up=total_up+term1+term2+term3
            else
                total_bellow=total_bellow+term1+term2+term3
            end if
        end if

        if(typ_force_z(i,j,k-1).eq.1) then
            !  term1=-0.5d0*(u2(i,j,k)+u2(i,j,k-1))*&
            ! 0.5d0*(u3(i,j,k-1)+u3(i,j+1,k-1))*dx1(i)*dx2(j)
            term2 = -p(i,j,k)/rhoKOH*dx1(i)*dx2(j)   
            term3 = nuKOH*(u3(i,j,k) - u3(i,j,k-1))/dx3(k)*dx1(i)*dx2(j)               
            total = total+term1+term2+term3
            if (x2_j(j).gt.yc) then
                total_up=total_up+term1+term2+term3
            else
                total_bellow=total_bellow+term1+term2+term3
            end if
        end if
    end do  

    FzCOM = total

if ( .false. ) then
1   format(4E36.24)
    if ((m1.eq.1).and.(m2.eq.1)) then
        open(unit=1,file='./output/FzBubble_top_bottom.dat')
2       format(4A36)
        write(1,2) 'time', 'total_up', 'total_bellow', 'FyCOM'
        write(1,1) time, total_up, total_bellow, FzCOM
        close(1,status='keep')
    end if
    if ( (m1.gt.1).and.(m2.eq.m2_max) ) then
        open(unit=1,file='./output/FzBubble_top_bottom.dat',position='append')
        write(1,1) time, total_up, total_bellow, FzCOM
        close(1,status='keep')
    end if
end if
   
end subroutine bulb_force_z