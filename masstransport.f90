subroutine Speciestransport()
    !!------------------------------------------------------------
    !!  **Main subroutine for Species transport**
    !!  
    !!  Determines the distribution of potential and concentration
    !!  in the domain including the bubble.
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    
    call prep_phi
    
    call spatial_phi

    call BCphi

    call spatial_c

    call BCc

end subroutine Speciestransport




subroutine prep_phi()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k
    real(8) :: coef0_tmp

!$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=0,km
        do j=0,jm
            do i=0,im
                work1(i,j,k) = (c(i,j,k,2) + c(i+1,j,k,2))/dx1_i(1)
                work2(i,j,k) = (c(i,j,k,2) + c(i,j+1,k,2))/dx2_j(1)
                work3(i,j,k) = (c(i,j,k,2) + c(i,j,k+1,2))/dx3_k(1)
            end do
        end do
    end do

!$OMP PARALLEL DO PRIVATE(coef0_tmp,i,j,k)
    do k=1,km
        do j=1,jm
            do i=1,im
                coef0_tmp = -work1(i-1,j,k) - work1(i,j,k)  &
                           - work2(i,j-1,k) - work2(i,j,k)  &
                           - work3(i,j,k-1) - work3(i,j,k)

                coef1_phi(i,j,k) = work1(i-1,j,k)/coef0_tmp
                coef2_phi(i,j,k) = work1(i,j,k)  /coef0_tmp
                coef3_phi(i,j,k) = work2(i,j-1,k)/coef0_tmp
                coef4_phi(i,j,k) = work2(i,j,k)  /coef0_tmp
                coef5_phi(i,j,k) = work3(i,j,k-1)/coef0_tmp
                coef6_phi(i,j,k) = work3(i,j,k)  /coef0_tmp
            end do
        end do
    end do

    do k=1,km
        do j=1,jm
            ka(j,k) = ka0*(0.5d0*(c(0,j,k,2) + c(1,j,k,2))/c_ref(2))*   &
                     sqrt(abs(0.5d0*(c(0,j,k,1) + c(1,j,k,1)))/c_ref(1)) 
            kc(j,k) = kc0*0.5d0*(c(0,j,k,3) + c(1,j,k,3))/c_ref(3) 
        end do
    end do

    ka(0,0:km)      = ka(1,0:km)
    ka(jm+1,0:km)   = ka(jm,0:km)
    ka(0:jm+1,0)    = ka(0:jm+1,1)
    ka(0:jm+1,km+1) = ka(0:jm+1,km)

    kc(0,0:km)      = kc(1,0:km)
    kc(jm+1,0:km)   = kc(jm,0:km)
    kc(0:jm+1,0)    = kc(0:jm+1,1)
    kc(0:jm+1,km+1) = kc(0:jm+1,km)
    
end subroutine prep_phi




subroutine spatial_phi()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k

    store_phi = phi

    do k=1,km
        do j=1,jm
            do i=1,im
                if (IBtype(i,j,k).eq.0) then
                    store_phi(i,j,k) = -(coef1_phi(i,j,k)*phi(i-1,j,k)  &
                                       + coef2_phi(i,j,k)*phi(i+1,j,k)  &
                                       + coef3_phi(i,j,k)*phi(i,j-1,k)  &
                                       + coef4_phi(i,j,k)*phi(i,j+1,k)  &
                                       + coef5_phi(i,j,k)*phi(i,j,k-1)  &
                                       + coef6_phi(i,j,k)*phi(i,j,k+1))
                end if
            end do
        end do
    end do

    phi = store_phi

end subroutine spatial_phi




subroutine BCphi()
    use commondata
    use paramdata
    implicit none
    integer :: i, j, k, n, nIB
    real(8) :: az,bz,cz,fac,kka,kkc,kkae,kkce

!$OMP PARALLEL DO PRIVATE(cz,fac,kka,kkc,kkae,kkce,bz,az)
    do k=1,km
        do j=1,jm
            do i=0,im
                conduct(i,j,k) = condfac*0.5d0*(c(i,j,k,2) + c(i+1,j,k,2)) 
            end do
  
         cz   = phi(1,j,k) - phiLe
         fac  = -0.5d0*dx1_i(0)/conduct(0,j,k)
         kka  = fac*ka(j,k)
         kkc  = fac*kc(j,k)
         kkae = kka*exp(aa*eta(j,k))
         kkce = kkc*exp(-ac*eta(j,k))      
         bz   = -eta(j,k) + kkce + kkae 
         az   = -1d0 - ac*kkce + aa*kkae

         eta(j,k)  = eta(j,k) + (cz-bz)/az
         phiL(j,k) = -eta(j,k) + phiLe

         phi(0,j,k)    = 2d0*phiL(j,k) - phi(1,j,k)   
         phi(im+1,j,k) = 2d0*phiR - phi(im,j,k) 
  
        end do
    end do

    phi(0:im+1,0,1:km)      = phi(0:im+1,1,1:km)
    phi(0:im+1,jm+1,1:km)   = phi(0:im+1,jm,1:km)
    phi(0:im+1,0:jm+1,0)    = phi(0:im+1,0:jm+1,km)
    phi(0:im+1,0:jm+1,km+1) = phi(0:im+1,0:jm+1,1)

    if ( SimBubble ) then
        do n = 1,nBubbles
            do nIB = 1,BubbleBlock(n)%nIB
                call IBM_potential(nIB,n)
            end do
        end do
    end if

    conduct(0:im,0,1:km)      = conduct(0:im,1,1:km)
    conduct(0:im,jm+1,1:km)   = conduct(0:im,jm,1:km)
    conduct(0:im,0:jm+1,0)    = conduct(0:im,0:jm+1,1)
    conduct(0:im,0:jm+1,km+1) = conduct(0:im,0:jm+1,km)

    eta(0,1:km)      = eta(1,1:km)
    eta(jm+1,1:km)   = eta(jm,1:km)
    eta(0:jm+1,0)    = eta(0:jm+1,1)
    eta(0:jm+1,km+1) = eta(0:jm+1,km)

end subroutine BCphi




subroutine IBM_potential(nIB, bubbleID)
    use commondata
    use paramdata
    use IBM_functions
    implicit none
    integer, intent(in) :: nIB, bubbleID
    real(8), dimension(4) :: probeXYZsdist
    integer :: id, jd, kd, c000I, c000J, c000K
    real(8) ::  c000, c001, c011, c010, &
                c110, c100, c101, c111, &
                probeVal

    id = IBcellLocation(1,nIB,bubbleID)
    jd = IBcellLocation(2,nIB,bubbleID)
    kd = IBcellLocation(3,nIB,bubbleID)

    ! get probe location
    probeXYZsdist = IBM_get_probe_location(                 &
                            x1(id),    x2(jd),    x3(kd),   &
                            dx1_i(id), dx2_j(jd), dx3_k(kd),&
                            BubbleBlock(bubbleID)%xCenter,  &
                            BubbleBlock(bubbleID)%yCenter,  &
                            BubbleBlock(bubbleID)%zCenter,  &
                            BubbleBlock(bubbleID)%radius    )

    ! Store probe location for future use or debugging
    probeCell(1,nIB,bubbleID) = probeXYZsdist(1)
    probeCell(2,nIB,bubbleID) = probeXYZsdist(2)
    probeCell(3,nIB,bubbleID) = probeXYZsdist(3)

    ! Get c000 index location for trilinear interpolation
    c000I = floor(probeXYZsdist(1)/dx1_i(id) + 0.5d0)
    c000J = floor(probeXYZsdist(2)/dx2_j(jd) + 0.5d0)
    c000K = floor(probeXYZsdist(3)/dx3_k(kd) + 0.5d0)

    ! Perform trilinear interpolation phi
    c000 = phi(c000I,     c000J,     c000K    )
    c001 = phi(c000I,     c000J,     c000K + 1)
    c011 = phi(c000I,     c000J + 1, c000K + 1)
    c010 = phi(c000I,     c000J + 1, c000K    )
    c110 = phi(c000I + 1, c000J + 1, c000K    )
    c100 = phi(c000I + 1, c000J,     c000K    )
    c101 = phi(c000I + 1, c000J,     c000K + 1)
    c111 = phi(c000I + 1, c000J + 1, c000K + 1)

    probeVal = IBM_probe_interp_trilinear(            & 
                    c000, c001, c011, c010,             &
                    c110, c100, c101, c111,             &
                    dx1_i(id), dx2_j(jd), dx3_k(kd),    &
                    probeXYZsdist(1), probeXYZsdist(2), &
                    probeXYZsdist(3),                   &
                    x1(c000I), x2(c000J), x3(c000K)     )

    ! Apply zero neumann condition
    phi(id,jd,kd) = probeVal
    
end subroutine IBM_potential




subroutine spatial_c()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k,nConc
    real(8) ::ar1,ar2,ar3, balance_KOH, balance_H2, balance_H2O

    store_c = c

    do nConc = 1,nspec
!$OMP PARALLEL DO PRIVATE(i,j,k,ar1,ar2,ar3)
        do k=0,km
            do j=0,jm 
                do i=0,im
                    work1(i,j,k) = dif(nConc)*(c(i+1,j,k,nConc) - c(i,j,k,nConc))/dx1_i(1) 
                    work2(i,j,k) = dif(nConc)*(c(i,j+1,k,nConc) - c(i,j,k,nConc))/dx2_j(1)  
                    work3(i,j,k) = dif(nConc)*(c(i,j,k+1,nConc) - c(i,j,k,nConc))/dx3_k(1) 
                end do
            end do
        end do

!$OMP PARALLEL DO PRIVATE(ar1,ar2,ar3)
        do k=1,km
            do j=1,jm
                do i=1,im
                    if (IBtype(i,j,k).eq.fluid_label) then
                        rhs(i,j,k) = (work1(i,j,k) - work1(i-1,j,k))/dx1(i) &
                                   + (work2(i,j,k) - work2(i,j-1,k))/dx2(j) &
                                   + (work3(i,j,k) - work3(i,j,k-1))/dx3(k)
                    
                        ar1 = 0.5d0*(u1(i,j,k) + u1(i-1,j,k))
                        ar2 = 0.5d0*(u2(i,j,k) + u2(i,j-1,k))
                        ar3 = 0.5d0*(u3(i,j,k) + u3(i,j,k-1))
            
                        if (ar1.gt.0d0) then
                            rhs(i,j,k) = rhs(i,j,k) - ar1*(c(i,j,k,nConc)   &
                                        - c(i-1,j,k,nConc))/dx1_i(0)
                        else
                            rhs(i,j,k) = rhs(i,j,k) - ar1*(c(i+1,j,k,nConc) &
                                        - c(i,j,k,nConc))/dx1_i(0)
                        end if
                
                        if (ar2.gt.0d0) then
                            rhs(i,j,k) = rhs(i,j,k) - ar2*(c(i,j,k,nConc)   &
                                        - c(i,j-1,k,nConc))/dx2_j(0)
                        else
                            rhs(i,j,k) = rhs(i,j,k) - ar2*(c(i,j+1,k,nConc) &
                                        - c(i,j,k,nConc))/dx2_j(0)
                        end if
                
                        if (ar3.gt.0d0) then
                            rhs(i,j,k) = rhs(i,j,k) - ar3*(c(i,j,k,nConc)   &
                                        - c(i,j,k-1,nConc))/dx3_k(0)
                        else
                            rhs(i,j,k) = rhs(i,j,k) - ar3*(c(i,j,k+1,nConc) &
                                        - c(i,j,k,nConc))/dx3_k(0)
                        end if
                    end if
                end do
            end do
        end do

!$OMP PARALLEL DO PRIVATE(i,j,k)
        do k=1,km
            do j=1,jm
                do i=1,im
                    if (IBtype(i,j,k).eq.fluid_label) then
                        c(i,j,k,nConc) = store_c(i,j,k,nConc) + dtime*rhs(i,j,k)
                    end if
                end do
            end do
        end do
    end do
    

    balance_H2  = 0.0d0
    balance_KOH = 0.0d0
    balance_H2O = 0.0d0

    do k=1,km
        do j=1,jm
            do i=1,im
                if (IBtype(i,j,k).eq.fluid_label) then
                    balance_H2 = balance_H2 + (c(i,j,k,1) - store_c(i,j,k,1))/dtime*dx1_i(0)*dx2_j(0)*dx3_k(0)
                end if
            end do
        end do
    end do

    massBalance(1) = balance_H2
    massBalance(2) = balance_KOH
    massBalance(3) = balance_H2O

end subroutine spatial_c




subroutine BCc()
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k,nConc, n, nIB
    real(8) :: dphi, ER_c, Max_Err_c

    do k=1,km
        do j=1,jm   
            dphi = conduct(0,j,k)*(phi(1,j,k) - phi(0,j,k))/Fa       

            c(0,j,k,1) = c(1,j,k,1) + 0.5d0*dphi/dif(1)
            c(0,j,k,2) = c(1,j,k,2) + 0.5d0*dphi/dif(2)
            c(0,j,k,3) = c(1,j,k,3) - dphi/dif(3)
        end do
    end do
! 1   format(1A30,1E20.10)
    ! write(*,1) 'dphi: ', dphi

    do nConc = 1,nspec
!$OMP PARALLEL DO PRIVATE(j,k)
        do k=1,km
            do j=1,jm 
                c(im+1,j,k,nConc) = 2d0*c_ref(nConc) - c(im,j,k,nConc)
            end do
        end do

        c(0:im+1,0,1:km,nConc)      = 2d0*c_ref(nConc) - c(0:im+1,1,1:km,nConc)        
        c(0:im+1,jm+1,1:km,nConc)   = c(0:im+1,jm,1:km,nConc)
        c(0:im+1,0:jm+1,0,nConc)    = c(0:im+1,0:jm+1,km,nConc)        
        c(0:im+1,0:jm+1,km+1,nConc) = c(0:im+1,0:jm+1,1,nConc)
    end do

    iterC = 0
    if ( SimBubble ) then
        ER_c      = 1e-3    ! Target value
        Max_Err_c = 5d0     ! initialization of error
        do while (Max_Err_c.gt.ER_c)
            do nConc=2,3
                do k=0,km+1
                    do j=0,jm+1
                        do i=0,im+1
                            store_c(i,j,k,nConc) = c(i,j,k,nConc)
                        end do
                    end do
                end do
            end do

            do n = 1,nBubbles
                do nIB = 1,BubbleBlock(n)%nIB
                    call IBM_concentration(nIB,n)                    
                end do
            end do

            iterC = iterC + 1

            do nConc=2,3
                do k=0,km+1
                    do j=0,jm+1
                        do i=0,im+1
                            ErrIter_c(i,j,k,nConc) = abs(c(i,j,k,nConc) - store_c(i,j,k,nConc))             
                        end do
                    end do
                end do
            end do

            Max_Err_c=maxval(ErrIter_c)
            if ( iterC.eq.500 ) then
                write(*,*) 'Warning! iterC > 500'
                write(*,*) 'Concentration value: ', Max_Err_c
            end if
            if ( iterC.eq.750 ) then
                write(*,*) 'Warning! iterC > 750'
                write(*,*) 'Concentration value: ', Max_Err_c
            end if
            if (iterC.gt.1000) then
                write(*,*) 'Stopping Simulation ---> Iteration problem in IBM_concentration.'
                write(*,*) 'Concentration value: ', Max_Err_c
                write(*,*) 'Number of iterations: ', iterC
                call writeOutput_IBM
                stop
            end if
       end do
    end if
    
end subroutine BCc


