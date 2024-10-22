subroutine initialize()
      !!------------------------------------------
      !! Allocate memory for the arrays and load
      !! initial settings for the domain and the
      !! bubble. 
      !!
      !!------------------------------------------
      use var2
      implicit none

      ! Allocate arrays into memory.
      call initArrays()

      ! Load the domain into the arrays.
      call initLoadDomain()

      ! Load bubble initial conditions
      call initLoadBubble()

      ! Check the parameters defined in the module for possible conflicts.
      call initCheckConfig()

      !<------------- Time step settings
      time = 0d0
      c_w  = 1d-30 + 10d0*umax

      ! Set the explicit time step, we assume that the maximum diffusivity is dif(1) 

      !------------- Artificial Compressibility parameter
      SELECT CASE(artCompresMethod)
            CASE('modified')
                  dtime = 0.5d0*min( co*dx1_i(0)/c_w, 0.25d0*(dx1_i(0)**2)/nuKOH,     &
                                     0.25d0*(dx1_i(0)**2)/dif(1) )
                  dtime = dtime*ACmod
                  usou2 = 0.05d0*(dx1(0)/dtime)**2
            CASE('auto')
                  dtime = 0.5d0*min( co*dx1_i(0)/c_w, 0.25d0*(dx1_i(0)**2)/nuKOH,     &
                  0.25d0*(dx1_i(0)**2)/dif(1) )
                  usou2 = 0.05d0*(dx1(0)/dtime)**2
            CASE DEFAULT
                  ! Same case as auto
                  dtime = 0.5d0*min( co*dx1_i(0)/c_w, 0.25d0*(dx1_i(0)**2)/nuKOH,     &
                                     0.25d0*(dx1_i(0)**2)/dif(1) )
                  usou2 = 0.05d0*(dx1(0)/dtime)**2
      END SELECT
      
end subroutine initialize


subroutine initArrays()
      !!------------------------------------------
      !! Initialize all the arrays specified in
      !! the module and check if any errors 
      !! occur.
      !!
      !!------------------------------------------
      use var2
      integer :: istat

      write(*,*) 'Allocating memory' 
      ! Rank 1 arrays:
      allocate( &
            x1(0:im+1),x1_i(-1:im+1), &
            dx1(0:im+1),dx1_i(0:im),  &
            x2(0:jm+1),x2_j(-1:jm+1), &
            dx2(0:jm+1),dx2_j(0:jm),  &
            x3(0:km+1),x3_k(-1:km+1), &
            dx3(0:km+1),dx3_k(0:km),  &
            c_s(1:nmax), c_ref(nmax), dif(nmax), &
            xb(1:im*jm*km/2),yb(1:im*jm*km/2),zb(1:im*jm*km/2), &
            phiS(1:im*jm*km/2),thetaS(1:im*jm*km/2),    &
            T_ni(1:im*jm*km/2),T_nj(1:im*jm*km/2),  &
            T_nk(1:im*jm*km/2), &
            T_ti(1:im*jm*km/2),T_tj(1:im*jm*km/2),  &
            T_tk(1:im*jm*km/2), &
            T_si(1:im*jm*km/2),T_sj(1:im*jm*km/2),  &
            T_sk(1:im*jm*km/2), &
            rb(1:im*jm*km/2),   &
            xb_c(1:im*jm*km/2),yb_c(1:im*jm*km/2),  &
            zb_c(1:im*jm*km/2), &
            phiS_U(1:im*jm*km/2),thetaS_U(1:im*jm*km/2),    &
            phiS_V(1:im*jm*km/2),thetaS_V(1:im*jm*km/2),    &
            phiS_W(1:im*jm*km/2),thetaS_W(1:im*jm*km/2),    &

            T_ni_U(1:im*jm*km/2),T_nj_U(1:im*jm*km/2),  &
            T_nk_U(1:im*jm*km/2),   &
            T_ti_U(1:im*jm*km/2),T_tj_U(1:im*jm*km/2),  &
            T_tk_U(1:im*jm*km/2),   &
            T_si_U(1:im*jm*km/2),T_sj_U(1:im*jm*km/2),  &
            T_sk_U(1:im*jm*km/2),   &

            T_ni_V(1:im*jm*km/2),T_nj_V(1:im*jm*km/2),  &
            T_nk_V(1:im*jm*km/2),   &
            T_ti_V(1:im*jm*km/2),T_tj_V(1:im*jm*km/2),  &
            T_tk_V(1:im*jm*km/2),   &
            T_si_V(1:im*jm*km/2),T_sj_V(1:im*jm*km/2),  &
            T_sk_V(1:im*jm*km/2),   &

            T_ni_W(1:im*jm*km/2),T_nj_W(1:im*jm*km/2),  &
            T_nk_W(1:im*jm*km/2),   &
            T_ti_W(1:im*jm*km/2),T_tj_W(1:im*jm*km/2),  &
            T_tk_W(1:im*jm*km/2),   &
            T_si_W(1:im*jm*km/2),T_sj_W(1:im*jm*km/2),  &
            T_sk_W(1:im*jm*km/2),   &

            rb_U(1:im*jm*km/2), xb_U(1:im*jm*km/2),  &
            yb_U(1:im*jm*km/2), zb_U(1:im*jm*km/2),  &
            rb_V(1:im*jm*km/2), xb_v(1:im*jm*km/2),  &
            yb_V(1:im*jm*km/2), zb_V(1:im*jm*km/2),  &
            rb_W(1:im*jm*km/2), xb_W(1:im*jm*km/2),  &
            yb_W(1:im*jm*km/2), zb_W(1:im*jm*km/2),  &
            x_IB(1:im*jm*km/2),y_IB(1:im*jm*km/2),  &
            z_IB(1:im*jm*km/2),direction(1:im*jm*km/2), &
            x1_plus(1:im*jm*km/2),x1_minus(1:im*jm*km/2),   &
            x2_plus(1:im*jm*km/2),x2_minus(1:im*jm*km/2),   &
            x3_plus(1:im*jm*km/2),x3_minus(1:im*jm*km/2),   &
            x_IB_U(1:im*jm*km/2),y_IB_U(1:im*jm*km/2),  &
            z_IB_U(1:im*jm*km/2),   &
            x_IB_V(1:im*jm*km/2),y_IB_V(1:im*jm*km/2),  &
            z_IB_V(1:im*jm*km/2),   &
            x_IB_W(1:im*jm*km/2),y_IB_W(1:im*jm*km/2),  &
            z_IB_W(1:im*jm*km/2), stat=istat )

      ! Rank 2 arrays:
      allocate( &
            eta(0:jm+1,0:km+1),phiL(0:jm+1,0:km+1), &
            ka(0:jm+1,0:km+1),kc(0:jm+1,0:km+1),    &
            x_prob(1:im*jm*km/2,1),y_prob(1:im*jm*km/2,1),  &
            z_prob(1:im*jm*km/2,1),   &

            x_prob_c(1:im*jm*km/2,1:2), &
            y_prob_c(1:im*jm*km/2,1:2), &
            z_prob_c(1:im*jm*km/2,1:2), &

            x_prob_U(1:im*jm*km/2,1:2),   &
            y_prob_U(1:im*jm*km/2,1:2),   &
            z_prob_U(1:im*jm*km/2,1:2),   &
            x_prob_V(1:im*jm*km/2,1:2),   &
            y_prob_V(1:im*jm*km/2,1:2),   &
            z_prob_V(1:im*jm*km/2,1:2),   &
            x_prob_W(1:im*jm*km/2,1:2),     &
            y_prob_W(1:im*jm*km/2,1:2),     &
            z_prob_W(1:im*jm*km/2,1:2),     &
            nx_prob(1:im*jm*km/2,1),        &
            ny_prob(1:im*jm*km/2,1),        &
            nz_prob(1:im*jm*km/2,1),        &
            neighbor(1:im*jm*km/2,1:6),     &

            nx_prob_c(1:im*jm*km/2,1:2),    &
            ny_prob_c(1:im*jm*km/2,1:2),    &
            nz_prob_c(1:im*jm*km/2,1:2),    &
            nx_prob_UU(1:im*jm*km/2,1:2),   &
            ny_prob_UU(1:im*jm*km/2,1:2),   &
            nz_prob_UU(1:im*jm*km/2,1:2),   &
            nx_prob_UV(1:im*jm*km/2,1:2),   &
            ny_prob_UV(1:im*jm*km/2,1:2),   &
            nz_prob_UV(1:im*jm*km/2,1:2),   &   
            nx_prob_UW(1:im*jm*km/2,1:2),   &
            ny_prob_UW(1:im*jm*km/2,1:2),   &
            nz_prob_UW(1:im*jm*km/2,1:2),   &

            nx_prob_VU(1:im*jm*km/2,1:2),   &
            ny_prob_VU(1:im*jm*km/2,1:2),   &
            nz_prob_VU(1:im*jm*km/2,1:2),   &
            nx_prob_VV(1:im*jm*km/2,1:2),   &
            ny_prob_VV(1:im*jm*km/2,1:2),   &
            nz_prob_VV(1:im*jm*km/2,1:2),   &
            nx_prob_VW(1:im*jm*km/2,1:2),   &
            ny_prob_VW(1:im*jm*km/2,1:2),   &
            nz_prob_VW(1:im*jm*km/2,1:2),   &

            nx_prob_WU(1:im*jm*km/2,1:2),   &
            ny_prob_WU(1:im*jm*km/2,1:2),   &
            nz_prob_WU(1:im*jm*km/2,1:2),   &
            nx_prob_WV(1:im*jm*km/2,1:2),   &
            ny_prob_WV(1:im*jm*km/2,1:2),   &
            nz_prob_WV(1:im*jm*km/2,1:2),   &
            nx_prob_WW(1:im*jm*km/2,1:2),   &
            ny_prob_WW(1:im*jm*km/2,1:2),   &
            nz_prob_WW(1:im*jm*km/2,1:2),   &
            stat=istat)

      ! Rank 3 arrays:
      allocate( &
            u1(0:im+1,0:jm+1,0:km+1),u2(0:im+1,0:jm+1,0:km+1), &
            u3(0:im+1,0:jm+1,0:km+1),p(0:im+1,0:jm+1,0:km+1),  &
            rhu(0:im,0:jm,0:km),rhv(0:im,0:jm,0:km),           &
            rhw(0:im,0:jm,0:km), &
            div(1:im,1:jm,1:km), &
            u_node(0:im,0:jm,0:km),v_node(0:im,0:jm,0:km), &
            w_node(0:im,0:jm,0:km),p_node(0:im,0:jm,0:km), &
            work1(0:im,0:jm,0:km), work2(0:im,0:jm,0:km),&
            work3(0:im,0:jm,0:km),  &
            phi(0:im+1,0:jm+1,0:km+1),  &
            coef1_phi(1:im,1:jm,1:km),coef2_phi(1:im,1:jm,1:km),    &
            coef3_phi(1:im,1:jm,1:km),coef4_phi(1:im,1:jm,1:km),    &
            coef5_phi(1:im,1:jm,1:km),coef6_phi(1:im,1:jm,1:km),    &
            current(0:im,0:jm+1,0:km+1),conduct(0:im,0:jm+1,0:km+1),&
            storage_u(0:im+1,0:jm+1,0:km+1),    &
            storage_v(0:im+1,0:jm+1,0:km+1),    &
            storage_w(0:im+1,0:jm+1,0:km+1),    &
            storage_p(0:im+1,0:jm+1,0:km+1),    &
            typ(0:im+1,0:jm+1,0:km+1),  &
            typ1(0:im+1,0:jm+1,0:km+1), &
            typ2(0:im+1,0:jm+1,0:km+1), &
            typ3(0:im+1,0:jm+1,0:km+1), &
            typ1a(0:im+1,0:jm+1,0:km+1),    &
            typ2a(0:im+1,0:jm+1,0:km+1),    &
            typ3a(0:im+1,0:jm+1,0:km+1),    &
            typ_IB(0:im+1,0:jm+1,0:km+1),   &
            typ_u_IB(0:im+1,0:jm+1,0:km+1), &
            typ_v_IB(0:im+1,0:jm+1,0:km+1), &
            typ_w_IB(0:im+1,0:jm+1,0:km+1), stat=istat  )
      
      ! Rank 4 arrays:
      allocate( &
            c(0:im+1,0:jm+1,0:km+1,nmax), &
            storage_c(0:im+1,0:jm+1,0:km+1,1:3), stat=istat)

      if ( istat.eq.0 ) then
            write(*,*) 'Allocating memory complete.'
      else
            write(*,*) 'Problem in allocating memory.' 
      end if
end subroutine initArrays



subroutine initLoadDomain()
      !!------------------------------------------
      !! Load grid coordinates and correspondig 
      !! values into the domain.
      !!
      !! The following parameters are loaded into
      !! the arrays:
      !!  - grid Carthesian coordinates
      !!  - potential distribution at walls
      !!  - initial velocity in the domain
      !!  - initial pressure in the domain
      !!  - initial concentration in the domain
      !!------------------------------------------
      use var2
      implicit none
      integer ::  n
      real(8) ::  frt,alphac

      ! Define
      ! x1=x1L=0 at the electrode (left side of domain)
      ! x1=x1R=lx1 at a location in the bulk flow (right side of domain)
      ! domain length is lx1

      ! Cell Face x
      do i=-1,im+1
            x1_i(i)=dfloat(i)*lx1/dfloat(im)
      end do 
      
      ! Cell Center x
      do i=0,im+1
            x1(i)=0.5d0*(x1_i(i-1)+x1_i(i))
            dx1(i)=x1_i(i)-x1_i(i-1)
      end do
     
      ! Distance between cell faces \deltax
      do i=0,im
            dx1_i(i)=x1(i+1)-x1(i)
      end do
      
      do j=-1,jm+1
            x2_j(j)=dfloat(j)*lx2/dfloat(jm)
      end do 
     
      do j=0,jm+1
            x2(j)=0.5d0*(x2_j(j-1)+x2_j(j))
            dx2(j)=x2_j(j)-x2_j(j-1)
      end do
     
      do j=0,jm
            dx2_j(j)=x2(j+1)-x2(j)
      end do
     
      do k=-1,km+1
            x3_k(k)=dfloat(k)*lx3/dfloat(km)
      end do 
      
      do k=0,km+1
            x3(k)=0.5d0*(x3_k(k-1)+x3_k(k))
            dx3(k)=x3_k(k)-x3_k(k-1)
      end do
     
      do k=0,km
            dx3_k(k)=x3(k+1)-x3(k)
      end do


            
!     phiLe = potential of electrode plus standard potential 
!     phiL  = phiLe - eta 
!     phiR  = potential at right hand side

      frt = Far/Rid/temp
      alphac = 0.5d0

      ac = alphac*frt
      aa = (1d0 - alphac)*frt

!     1: dissolved H2 concentration 
!     2: electrolyte concentration (1m: OHmin, 1p, Kplus)
!     3: H2O concentration

!     reference conditions to be applied in the kinetics
!     and also used as bulk flow boundary condition

      c_ref(1:3) = (/C0H2, C0KOH, C0H20/)
      ! c_s(1)=c_ref(1)

!     Diffusivities of K+ and OH- are assumed to be the same. 
!     Electrolyte concentration is high and nonequal Fickian diffusivities 
!     are not realistic and consistent with the mass averaged velocity. 
      dif(1:3) = (/DH2, DKOH, DH2O/)

      condfac=Far*frt*2d0*dif(2)

      ! Load initial conditions into the domain arrays.
      do k=0,km+1
            do j=0,jm+1 
                  eta(j,k)  = 0d0
                  phiL(j,k) = phiLe
            end do
      end do


      ! Load quiescient flow profile for Stokes flow
      SELECT CASE (FlowCondPreset)
            CASE('Stokes')
!$omp parallel do private(i,j,k)
                  do k=0,km+1
                        do j=0,jm+1
                              do i=0,im+1
                                    phi(i,j,k)= phiR
                                    p(i,j,k)  = 0d0
                                    u1(i,j,k) = 0d0
                                    u2(i,j,k) = 0d0
                                    u3(i,j,k) = 0d0
                              end do
                        end do
                  end do

            CASE DEFAULT
!$omp parallel do private(i,j,k)
                  do k=0,km+1
                        do j=0,jm+1
                              do i=0,im+1
                                    phi(i,j,k)= phiR
                                    p(i,j,k)  = 0d0
                                    u1(i,j,k) = 0d0
                                    u2(i,j,k) = 4d0*umax*(x1(i)/lx1)*(1d0 - (x1(i)/lx1))
                                    u3(i,j,k) = 0d0
                              end do
                        end do
                  end do
      END SELECT

      do n=1,nmax
!$omp parallel do private(i,j,k)
            do k=0,km+1
                  do j=0,jm+1
                        do i=0,im+1
                              c(i,j,k,n) = c_ref(n)
                        end do
                  end do
            end do
      end do

      ! Compute kinetic prefactors
      do k=1,km
            do j=1,jm
                  ka(j,k) = ka0*(0.5d0*(c(0,j,k,2) + c(1,j,k,2))/c_ref(2))*   &
                        sqrt(abs(0.5d0*(c(0,j,k,1) + c(1,j,k,1)))/c_ref(1)) 
                  kc(j,k) = kc0*0.5d0*(c(0,j,k,3) + c(1,j,k,3))/c_ref(3) 
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
end subroutine initLoadDomain




subroutine initLoadBubble()
      !!------------------------------------------
      !! Loads bubble initial conditions.
      !! The following parameters are loaded into
      !! the arrays:
      !!  - radius 
      !!  - position
      !!  - forces over bubble surface 
      !!  - velocity
      !!------------------------------------------
      use var2
      implicit none

      ! Bubble surface forces
      FxCOM = 0d0
      FyCOM = 0d0
      FzCOM = 0d0

      ! Bubble growth rate
      drdt = 0d0

      ! Bubble COM-velocity
      dxcdt = 0d0       ! Current time step
      dycdt = 0d0       !
      dzcdt = 0d0       !
      dxcdt_old = 0d0   ! Previous time step
      dycdt_old = 0d0
      dzcdt_old = 0d0   ! Previous Previous time step
      dxcdt_oldold = 0d0
      dycdt_oldold = 0d0
      dzcdt_oldold = 0d0

      ! Check if Stokes flow conditions are required
      SELECT CASE (FlowCondPreset)
            CASE ('Stokes')
                  ! Disable bubble growth by r>r_max
                  radius = 1.0005d0*R_max
                  
                  ! Initial location of bubble in fluid domain
                  xc = lx1/2d0
                  yc = lx2*1d0/3d0
                  zc = lx3/2d0

                  ! In Stokes flow bubble is detached
                  isDetached = .true.
            CASE DEFAULT
                  ! Initial bubble radius
                  radius  = 3d0*x1_i(1)
                  radius0 = radius
      
                  ! Initial location of bubble in fluid domain
                  xc = radius + 1d0*x1_i(1)
                  yc = lx2/2d0
                  zc = lx1/2d0
      
                  ! During growth bubble is attached
                  isDetached = .false.
      END SELECT
     
      ! For time integration scheme using virtual force
      ddxcdttFv = 0d0     ! acceleration Fv x-direction
      ddycdttFv = 0d0     ! acceleration Fv y-direction
      ddzcdttFv = 0d0     ! acceleration Fv z-direction
      fvx = 0d0   ! acceleration virtual force x-direction
      fvy = 0d0   ! acceleration virtual force y-direction
      fvz = 0d0   ! acceleration virtual force z-direction

      ! For bubble equation of motion (always needed)
      fsx = 0d0   ! acceleration exerted force x-direction
      fsy = 0d0   ! acceleration exerted force y-direction
      fsz = 0d0   ! acceleration exerted force z-direction
      fby = 0d0   ! acceleartion due to bouyancy
      
end subroutine initLoadBubble



subroutine initCheckConfig()
      !!------------------------------------------
      !! Check if the parameters defined in the
      !! module does not creat possible conflicts
      !!
      !! The following checks are performed:
      !!  - Virtual force and density ratio
      !!  - Resolution of domain and bubble
      !!------------------------------------------
      use var2
      implicit none

1     format(1A30,1A60)
      ! Density ratio check from Schwarz et al. (2015)
      if ((rhoH2/rhoKOH.le.1.2).and.(VirtForcMethod.eqv..false.)) then
            write(*,1) '!! WARNING: ', 'Density ratio < 1.2 and Virtual force method is dissabled.'
      end if

      ! Bubble resolution required dX <= R/12
      if ( (lx1/dble(im).gt.radius/12d0) ) then
            write(*,1) '!! Warning: ', 'Domain resolution dX>R/12; issue to resolve bubble properly.'
      end if
      if ( (lx2/dble(jm).gt.radius/12d0) ) then
            write(*,1) '!! Warning: ', 'Domain resolution dY>R/12; issue to resolve bubble properly.'
      end if
      if ( (lx3/dble(km).gt.radius/12d0) ) then
            write(*,1) '!! Warning: ', 'Domain resolution dZ>R/12, issue to resolve bubble properly.'
      end if

      ! Check if resolution is equidistant: dx = dy = dz
      ! if ( (lx1/dble(im).ge.lx2/dble(jm)).and.(lx1/dble(im).ge.lx3/dble(km)) ) then
      !       ! do nothing
      ! else
      !        write(*,1) '!! Warning: ', 'Domain resolution is not equisistand: dx =/= dy =/= dz.'
      ! end if

end subroutine initCheckConfig
     

