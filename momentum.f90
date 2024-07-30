subroutine NSmomentum()
    !!------------------------------------------------------------
    !!  **Main subroutine for Momentum Equation**
    !!  
    !!  Applies central differencing scheme to u, v, and w
    !!  momentum equation (complete form), excluding the pressure
    !!  term (for now).
    !!
    !!
    !!------------------------------------------------------------
    use commondata
    implicit none


    ! Determine Right Hand Side of the momentum equation solved
    ! with Euler explicit timestep
    call u_momentum
    call v_momentum
    call w_momentum

    ! Apply Boundary Conditions to RHS
    call u_momentum_BC
    call v_momentum_BC
    call w_momentum_BC

    ! Add pressure term u,v,w velocity
    call velocity

    ! Apply boundary conditions to u,v,w velocity
    call BCu1
    call BCu2
    call BCu3

    ! Use Artificial Compressibility Method for pressure correction
    call ArtCompres

    
end subroutine NSmomentum




subroutine u_momentum()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k
    real(8) :: tp, ua, va, wa

    store_u = u1

!$OMP PARALLEL DO PRIVATE(ua,va,wa,tp,i,j,k)
    do k = 1,km
        do j = 1,jm
            do i = 1,im-1
                ua =      0.5d0*(u1(i,j,k) + u1(i-1,j,k))
                tp =    -ua**2/dx1_i(i)
                ua =      0.5d0*(u1(i+1,j,k) + u1(i,j,k))
                tp = tp + ua**2/dx1_i(i)
        
                va =      0.5d0*(u2(i+1,j-1,k) + u2(i,j-1,k))
                tp = tp - 0.5d0*(u1(i,j-1,k) + u1(i,j,k))*va/dx2(j)
                va =      0.5d0*(u2(i+1,j,k) + u2(i,j,k))
                tp = tp + 0.5d0*(u1(i,j,k) + u1(i,j+1,k))*va/dx2(j)
        
                wa =      0.5d0*(u3(i+1,j,k-1) + u3(i,j,k-1))
                tp = tp - 0.5d0*(u1(i,j,k-1) + u1(i,j,k))*wa/dx3(k)
                wa =      0.5d0*(u3(i+1,j,k) + u3(i,j,k))
                tp = tp + 0.5d0*(u1(i,j,k) + u1(i,j,k+1))*wa/dx3(k)

                rhu(i,j,k) = store_u(i,j,k) - tp*dtime                  &                    
                            - nuKOH*dtime/dx1_i(i)*                     &
                            ((u1(i,j,k)    - u1(i-1,j,k))/dx1(i)        &
                            - (u1(i+1,j,k) - u1(i,j,k))/dx1(i+1))       &
                            - nuKOH*dtime/dx2(j)*                       &
                            ((u1(i,j,k)    -  u1(i,j-1,k))/dx2_j(j-1)   &
                            - (u1(i,j+1,k) - u1(i,j,k))/dx2_j(j))       &
                            - nuKOH*dtime/dx3(k)*                       &
                            ((u1(i,j,k)    - u1(i,j,k-1))/dx3_k(k-1)    &
                            - (u1(i,j,k+1) - u1(i,j,k))/dx3_k(k))
            end do
        end do
    end do

end subroutine u_momentum



subroutine v_momentum()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k
    real(8) :: tp, ua, va, wa

    store_v = u2

!$OMP PARALLEL DO PRIVATE(ua,va,wa,tp,i,j,k)
    do k=1,km
        do j=1,jm
            do i=1,im
                ua =      0.5d0*(u1(i-1,j,k) + u1(i-1,j+1,k))
                tp =     -0.5d0*(u2(i-1,j,k) + u2(i,j,k))*ua/dx1(i)
                ua =      0.5d0*(u1(i,j,k) + u1(i,j+1,k))
                tp = tp + 0.5d0*(u2(i,j,k) + u2(i+1,j,k))*ua/dx1(i)

                va =      0.5d0*(u2(i,j,k) + u2(i,j-1,k))
                tp = tp - va**2/dx2_j(j)
                va =      0.5d0*(u2(i,j,k) + u2(i,j+1,k))
                tp = tp + va**2/dx2_j(j)

                wa =      0.5d0*(u3(i,j,k-1) + u3(i,j+1,k-1))
                tp = tp - 0.5d0*(u2(i,j,k-1) + u2(i,j,k))*wa/dx3(k)
                wa =      0.5d0*(u3(i,j,k) + u3(i,j+1,k))
                tp = tp + 0.5d0*(u2(i,j,k) + u2(i,j,k+1))*wa/dx3(k)

                rhv(i,j,k) = store_v(i,j,k) - tp*dtime              &
                            -nuKOH*dtime/dx1(i)*                    &
                            ((u2(i,j,k) - u2(i-1,j,k))/dx1_i(i-1)   &
                            - (u2(i+1,j,k) - u2(i,j,k))/dx1_i(i))   &
                            - nuKOH*dtime/dx2_j(j)*                 &
                            ((u2(i,j,k) - u2(i,j-1,k))/dx2(j)       &
                            - (u2(i,j+1,k) - u2(i,j,k))/dx2(j+1))   &
                            - nuKOH*dtime/dx3(k)*                   &
                            ((u2(i,j,k) - u2(i,j,k-1))/dx3_k(k-1)   &
                            - (u2(i,j,k+1) - u2(i,j,k))/dx3_k(k))
            end do
        end do
    end do

end subroutine v_momentum



subroutine w_momentum()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k
    real(8) :: tp, ua, va, wa

    store_w = u3

!$OMP PARALLEL DO PRIVATE(ua,va,wa,tp,i,j,k)
    do k=1,km
        do j=1,jm
            do i=1,im
                ua =      0.5d0*(u1(i-1,j,k) + u1(i-1,j,k+1))
                tp =     -0.5d0*(u3(i-1,j,k) + u3(i,j,k))*ua/dx1(i)
                ua =      0.5d0*(u1(i,j,k) + u1(i,j,k+1))
                tp = tp + 0.5d0*(u3(i,j,k) + u3(i+1,j,k))*ua/dx1(i)
        
                va =      0.5d0*(u2(i,j-1,k) + u2(i,j-1,k+1))
                tp = tp - 0.5d0*(u3(i,j-1,k) + u3(i,j,k))*va/dx2(j)
                va =      0.5d0*(u2(i,j,k) + u2(i,j,k+1))
                tp = tp + 0.5d0*(u3(i,j,k) + u3(i,j+1,k))*va/dx2(j)
        
                wa =      0.5d0*(u3(i,j,k)+u3(i,j,k-1))
                tp = tp - wa**2/dx3_k(k)
                wa =      0.5d0*(u3(i,j,k+1)+u3(i,j,k))
                tp = tp + wa**2/dx3_k(k)
        
                rhw(i,j,k)= store_w(i,j,k) - tp*dtime               &
                            -nuKOH*dtime/dx1(i)*                    &
                            ((u3(i,j,k) - u3(i-1,j,k))/dx1_i(i-1)   &
                            - (u3(i+1,j,k) - u3(i,j,k))/dx1_i(i))   &
                            - nuKOH*dtime/dx2(j)*                   &
                            ((u3(i,j,k) - u3(i,j-1,k))/dx2_j(j-1)   &
                            -(u3(i,j+1,k) - u3(i,j,k))/dx2_j(j))    &
                            -nuKOH*dtime/dx3(k)*                    &
                            ((u3(i,j,k) - u3(i,j,k-1))/dx3(k)       &
                            -(u3(i,j,k+1) - u3(i,j,k))/dx3(k+1))
            end do
        end do
    end do
    
end subroutine w_momentum



subroutine u_momentum_BC()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k

    do k=1,km
        do j=1,jm
            do i=1,im
                if (IBtypeU1(i,j,k).eq.0) then
                    u1(i,j,k) = rhu(i,j,k)
               end if
            end do
        end do
    end do
    
    do k=1,km
        do j=1,jm
            u1(0,j,k)=0d0
            u1(im,j,k)=0d0 
        end do
    end do

end subroutine u_momentum_BC



subroutine v_momentum_BC()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k

    do k=1,km
        do j=1,jm
            do i=1,im
                if (IBtypeU2(i,j,k).eq.0) then
                    u2(i,j,k) = rhv(i,j,k)
                end if
            end do
        end do
    end do
  
    do k=1,km
         do i=1,im
            u2(i,0,k) = 4d0*umax*(x1(i)/lx1)*(1d0 - (x1(i)/lx1))
         end do
    end do
    
end subroutine v_momentum_BC



subroutine w_momentum_BC()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k

    do k=1,km
        do j=1,jm
            do i=1,im
                if(IBtypeU3(i,j,k).eq.0) then
                    u3(i,j,k) = rhw(i,j,k)
                end if
            end do
        end do
    end do
    
end subroutine w_momentum_BC



subroutine velocity()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k

!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=1,km
        do j=1,jm
            do i=1,im
                if(IBtypeU1(i,j,k).eq.0) then
                    u1(i,j,k) = u1(i,j,k) - dtime/rhoKOH/dx1_i(i)*  &
                                (p(i+1,j,k) - p(i,j,k))
                end if
        
                if(IBtypeU2(i,j,k).eq.0) then
                    u2(i,j,k) = u2(i,j,k) - dtime/rhoKOH/dx2_j(j)*  &
                                (p(i,j+1,k) - p(i,j,k))
                end if

                if(IBtypeU3(i,j,k).eq.0) then
                    u3(i,j,k) = u3(i,j,k) - dtime/rhoKOH/dx3_k(k)*  &
                                (p(i,j,k+1) - p(i,j,k))
                end if
            end do
        end do
    end do
    
end subroutine velocity



subroutine ArtCompres()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k

!$OMP PARALLEL DO
    do k=1,km
        do j=1,jm 
            do i=1,im
                div(i,j,k) = (u1(i,j,k) - u1(i-1,j,k))/dx1(i)   &
                           + (u2(i,j,k) - u2(i,j-1,k))/dx2(j)   &
                           + (u3(i,j,k) - u3(i,j,k-1))/dx3(k) 
            end do
        end do
    end do

    store_p = p

!$OMP PARALLEL DO
    do k=1,km
        do j=1,jm
            do i=1,im
                if (IBtype(i,j,k).eq.0) then
                    p(i,j,k) = store_p(i,j,k) - rhoKOH*usou2*dtime*div(i,j,k)
                end if
            end do
        end do
    end do

    call BCpressure

end subroutine ArtCompres