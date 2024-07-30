module paramdata
    !!--------------------------------------------------------------
    !!  **Purpose:**
    !!      Global parameters are defined here.
    !!
    !!      These values CANNOT be modified during the program (save).
    !!      
    !!      *Parameter defined are:*
    !!      - Universal constants
    !!      - Domain size (length and number of cells)
    !!      - Physical properties
    !!      - Electrochemical properties
    !!      - Simulation conditions
    !!      - Bubble parameters
    !!
    !!--------------------------------------------------------------
    implicit none
    save ! Guarantees that all data values declared in the module 
         ! will be preseved between references in different 
         ! procedures

    integer, parameter  ::  nspec   = 3,        &! num. species
                            im      = 96,       &! num. x-cells
                            jm      = 96,       &! num. y-cells
                            km      = 96,       &! num. z-cells
                            mSave   = 700,      &! Saving frequency
                            mComp   = 2000       ! Computing runs


    real(8), parameter  ::  pi = 3.141592d0,    &! [-]
                            Ru = 8.314462d0,    &! [J/mol K]
                            Fa = 96485.33d0,    &! [C/mol]
    ! Domain properties:
                            lx1 = 1d-4,         &! [m]
                            lx2 = 1d-4,         &! [m]
                            lx3 = 1d-4,         &! [m]
    ! Physical properties:
                            temp    = 353d0,    &! temp. [K]
                            umax    = 1d-3,     &! max. flow [m/s]
                            DifH2   = 5.8d-9,   &! diffusion [m2/s]
                            DifKOH  = 3.2d-9,   &! diffusion [m2/s]
                            DifH2O  = 3.2d-9,   &! diffusion [m2/s]
                            crefH2  = 0.16d0,   &! concntr. [mol/L]
                            crefKOH = 6700d0,   &! concntr. [mol/L]
                            crefH2O = 49000d0,  &! concntr. [mol/L]
                            rhoKOH  = 1258d0,   &! density [kg/m3]
                            rhoH2   = 8.375d-2, &! density [kg/m3]
                            nuKOH   = 6.7d-7,   &! kin. visc. [Pa s]
                            gammaKOH= 7.7d-2,   &! surf. tens. [N/m]

    ! Electrochemical properties:
                            phiLe   = -0.45d0,  &! elec.+std.pot.[V]
                            phiR    = 0.0d0,    &! memb. pot. [V]
                            j0      = 1.0d0,    &! exc.cur.den.[A/m2]
                            ka0     = j0,       &! Buttler-Volmer eq.
                            kc0     = -j0,      &! Buttler-Volmer eq.
                            alphac  = 0.5d0,    &! Cath. charg. [-]
                            alphaa  = 0.5d0,    &! Anod. charg. [-]
    ! Simplified constants:
                            FRT     = Fa/Ru/temp,        &! [C/J] = [1/V]
                            condfac = Fa*FRT*2d0*DifKOH, &! [C2*m2/mol/J/s]
                            aa      = alphac*FRT,        &! exponent coefficient
                            ac      = alphaa*FRT,        &! exponent coefficient
    ! Simulation stability criterion:
                            CFL     = 1d0/3d0,  &! Courant num. [-]
                            ErrMax_c= 5d0        ! Maximum error [-]
    
    ! Bubble parameters and IBM settings:
    logical, parameter  ::  SimBubble   = .true.,   &! Simulate bubble
                            SimDetach   = .false.     ! Simulate detachment

    integer, parameter  ::  nBubbles    = 1,        &! Number of bubbles
                            ghost_label = -1,       &! identifyer for IBM-method 
                            IB_label    = 1,        &! identifyer for IBM-method
                            fluid_label = 0,        &! identifyer for IBM-method
                            detachModel = 1          ! Select detachment model,
                                                     ! 1 = critical diameter,
                                                     ! 2 = force balance

    real(8), parameter  ::  radius0     = 3.0d0*lx1/dble(im),   &! initial bubble radius
                            xCenter0    = radius0+lx1/dble(im), &! initial bubble x-COM
                            yCenter0    = lx2/2.0d0,            &! initial bubble y-COM
                            zCenter0    = lx3/2.0d0,            &! initial bubble z-COM
                            radiusCrit  = 0.25d0*lx1             ! bubble critical radius

    ! Load previous simulation data
    logical, parameter  :: loadSimdata  = .true.       ! Load previous simulation data
    character(len=40),  &                              ! Location of previous simulation 
    parameter           :: dataPath     = './output/simData/'
    character(len=40),  &
    parameter           :: fileName     = 'simData'    ! filename (stored as a .dat)
                                                       ! NOTE: PARAMDATA SETTINGS MUST
                                                       ! MATCH BETWEEN SIMULATIONS!
    integer, parameter  :: m1CompStart  = 0105         ! Starting m1 loop
contains
    ! subroutine can be defined here but must use the parameters
    ! declared in this module only.
end module paramdata



module commondata
    !!--------------------------------------------------------------
    !!  **Purpose:**
    !!      Global arrays for parameters are defined here.
    !!      
    !!      These CAN be modified during the program.
    !!      
    !!      *Spatial coordinates:*
    !!          - cell centers : x1(i), x2(j), x3(k)
    !!          - cell faces : x1_i(i), x2_j(j), x3_k(k)
    !!          - cell center distance : dx1(i), dx2(j), dx3(k)
    !!          - cell face distance : dx1_i(i), dx2_j(j), 
    !!                                  dx3_k(k)
    !!
    !!      *Non-scalar paremeters:*
    !!      - velocity: u(i,j,k), v(i,j,k), w(i,j,k)
    !!
    !!      *Scalar parameters:*
    !!      - pressure: p(i,j,k)
    !!      - species concentration: c(i,j,k,n)
    !!      - potential: eta(j,k), phiL(j,k)
    !!      - Kinetic Prefactors: KA(j,k), KC(j,k)
    !!  
    !!      *Simulation parameter:*
    !!      - time step: dtime
    !!      - elapsed time: totaltime
    !!      - universal speed of sound: usou2
    !!
    !!--------------------------------------------------------------
    implicit none

    ! Important is to define the save-attribute as otherwise these
    ! arrays wont be stored correctly when being allocated in a 
    ! subroutine

    ! <<<< ------------- removed:  save-attribute ------------- >>>>

    integer :: m1, m2, m1Start, iterC
    real(8) :: dtime, totaltime, usou2

    real(8), allocatable, dimension(:)        ::              &
                                            x1,     x2,     x3,     &
                                            x1_i,   x2_j,   x3_k,   &
                                            dx1,    dx2,    dx3,    &
                                            dx1_i,  dx2_j,  dx3_k,  &
                                            c_ref,  dif,    c_s,    &
                                            massBalance

    real(8), allocatable, dimension(:,:)      ::              &
                                            eta,    phiL,   KA,     &
                                            KC,     ForceVector

    real(8), allocatable, dimension(:,:,:)    ::              &
                                            u1,     u2,     u3,     &
                                            p,      phi,    div,    &
                                            rhu,    rhv,    rhw,    &
                                            store_u,store_v,store_w,&
                                            store_p,    store_phi,  &
                                            coef1_phi,  coef2_phi,  &
                                            coef3_phi,  coef4_phi,  &
                                            coef5_phi,  coef6_phi,  &
                                            work1,  work2,  work3,  &
                                            nodeu1, nodeu2, nodeu3, &
                                            conduct, probeCell,     &
                                            probeU, probeV, probeW, &
                                            rhs

    real(8), allocatable, dimension(:,:,:,:)  ::              &
                                            c,  store_c, ErrIter_c

    integer, allocatable, dimension(:,:,:)    ::                    &
                                            IBtype, CellIsFlagged,  &
                                            IBtypeU1,  IBtypeU2,    &
                                            IBtypeU3,  IBbubble,    &
                                            IBbubbleU1,             &
                                            IBbubbleU2,             &
                                            IBbubbleU3,             &
                                            IBcellLocation,         &
                                            IBfaceXLocation,        &
                                            IBfaceYLocation,        &
                                            IBfaceZLocation

    ! Bubble parameters (Object Oriented Coding)
    type, public :: bubble
        integer :: iBubble
        integer :: nFacesX
        integer :: nFacesY
        integer :: nFacesZ
        integer :: nIB
        real(8) :: xCenter
        real(8) :: yCenter
        real(8) :: zCenter
        real(8) :: radius
        real(8) :: drdt
        real(8) :: massFlux
        real(8) :: u
        real(8) :: v
        real(8) :: w
        real(8) :: ForceX
        real(8) :: ForceY
        real(8) :: ForceZ
        logical :: hasDeparted
    end type bubble
    
    type(bubble), allocatable, dimension(:)  :: BubbleBlock

    ! Immersed Boundary Method (could be used to replace u1,u2 and u3 arrays)
    ! type, public :: cellFaceVal
    !     integer :: IBtype
    !     real(8) :: magnitude
    !     real(8) :: xLocation
    !     real(8) :: yLocation
    !     real(8) :: zLocation
    !     real(8) :: faceNormal
    ! end type cellFaceVal
    ! type(velocity), allocatable, dimension(:,:,:) :: u1, u2, u3

    ! Immersed Boundary Method (could be used to replace scalar variables arrays p, phi, c )
    ! type, public :: cellCenterVal
    !     integer :: IBtype
    !     real(8) :: value1
    !     real(8) :: value2
    !     real(8) :: value3
    ! end type cellCenterVal
    !type(cellCenterVal), allocatable, dimension(:,:,:) :: p, phi, c
    
contains
    subroutine checkAllocError(iVal,  errMsg)
        !!--------------------------------------------------------------
        !!  **Purpose:**
        !!      Check memory allocation for array's
        !!      - if succesfull: do nothing
        !!      - if failed: stop program and display error
        !!--------------------------------------------------------------
        implicit none
        integer, intent(in)             :: iVal
        character(len=60), intent(in)   :: errMsg

        if ( iVal.ne.0 ) then
1           format(1A60,I2)
            write(*,1) trim(errMsg), iVal

            error stop 'RAM-error, aborting program.'
        end if

    end subroutine checkAllocError

    subroutine checkNaNvalue(value, variable, subName)
        !!--------------------------------------------------------------
        !!  **Purpose:**
        !!      Check if value is NaN
        !!      - if true: abort program, dispaly location of NaN value,
        !!      and save data (without overwriting simData.dat!)
        !!      - if false: do nothing
        !!--------------------------------------------------------------
        real(8), intent(in) :: value
        character(len=20), intent(in) :: variable, subName
        character(len=60) :: errMsg

        if ( isnan(value) ) then
1           format(4A)
            write(errMsg,1) 'NaN detected for variable "', trim(variable), '" in ', trim(subName)
            write(*,1) errMsg
            ! Save IBM and bubble info data
            call writeOutput_IBM
            call writeOutput(2,1,0)
            error stop 'NaN-error occured, aborting program.'
        end if
    
    end subroutine checkNaNvalue
end module commondata



module IBM_functions
    implicit none
    
contains
    function IBM_get_probe_location(xIB, yIB, zIB, deltaX, deltaY,  &
                                    deltaZ, xBC, yBC, zBC, radius  )&
                                    result(probeLocationAbsolute)
        !!------------------------------------------------------------
        !!  **Purpose:** 
        !!      - Obtain probe location by computing the intersection
        !!      between IB-ghost suface normal vector and bubble
        !!      surface.
        !!  **Output:**
        !!      - Array containing probe point x-,y- and z-coordinate
        !!      along with distance between IB and surface
        !!
        !!------------------------------------------------------------
        implicit none
        real(8), intent(in) :: xIB      ! IB-point x-coordinate
        real(8), intent(in) :: yIB      ! IB-point y-coordinate
        real(8), intent(in) :: zIB      ! IB-point z-coordinate
        real(8), intent(in) :: deltaX   ! distance surface and probe
        real(8), intent(in) :: deltaY   ! distance surface and probe
        real(8), intent(in) :: deltaZ   ! distance surface and probe
        real(8), intent(in) :: xBC      ! Bubble center x-coordinate
        real(8), intent(in) :: yBC      ! Bubble center y-coordinate
        real(8), intent(in) :: zBC      ! Bubble center z-coordinate
        real(8), intent(in) :: radius   ! Bubble radius

        real(8), dimension(4) :: probeLocationAbsolute

        real(8) :: rdist, sdist, phi, theta, rVec_i, rVec_j, rVec_k
        real(8) :: xInter, yInter, zInter
        ! For reference, see the rotation matrix section on Wikipedia
        ! discussing Cartesian to Spherical coordinate conversion.

        ! Distance between IB-ghost and bubble center
        rdist   = sqrt((xIB - xBC)**2 + (yIB - yBC)**2 + (zIB - zBC)**2)

        ! Determine azimuth angle phi and inclination angle theta
        phi     = datan2((yIB - yBC), (xIB - xBC))
        theta   = acos((zIB - zBC)/rdist)

        ! Obtain radial (a.k.a. surface outward normal) direction
        rVec_i  = cos(phi)*sin(theta)
        rVec_j  = sin(phi)*sin(theta)
        rVec_k  = cos(theta)

        ! ! Surface tangential direction: inclination and azimuth
        ! phiVec_i  = -sin(phi)
        ! phiVec_j  = cos(phi)
        ! phiVec_k  = 0.0d0

        ! thetaVec_i  = cos(phi)*cos(theta)
        ! thetaVec_j  = sin(phi)*cos(theta)
        ! thetaVec_k  =-sin(theta)

        ! Distance between IB-ghost and bubble surface
        sdist   = radius - rdist

        ! Intersection between IB-ghost radial vector and bubble surface
        xInter  = xIB + sdist*rVec_i
        yInter  = yIB + sdist*rVec_j
        zInter  = zIB + sdist*rVec_k

        ! Calculate probe point at a cell deistance away from the surface
        probeLocationAbsolute(1)    = xInter + deltaX*rVec_i
        probeLocationAbsolute(2)    = yInter + deltaY*rVec_j
        probeLocationAbsolute(3)    = zInter + deltaZ*rVec_k
        probeLocationAbsolute(4)    = sdist

    end function IBM_get_probe_location
    

    

    function IBM_probe_interp_trilinear(c000, c001, c011, c010, &
                                        c110, c100, c101, c111, &
                                        deltaX, deltaY, deltaZ, &
                                        xp, yp, zp, x, y, z    )&
                                        result(probeValue)
        !!------------------------------------------------------------
        !!  **Purpose:** 
        !!      - Apply trilinear interpolation to probe (xp, yp, zp)
        !!      with c000 (x,y,z)
        !!
        !!  **Output:**
        !!      - Real(8) value of the probe point
        !!
        !!------------------------------------------------------------
        implicit none
        real(8), intent(in) ::  c000, c001, c011, c010,  &  ! Values at the corners
                                c110, c100, c101, c111       
        real(8), intent(in) ::  deltaX, deltaY, deltaZ      ! Distance between corners
        real(8), intent(in) ::  xp, yp, zp  ! Absolute coordinates probe point
        real(8), intent(in) ::  x, y, z     ! Absolute coordinate C000
        real(8)             ::  probeValue  ! Value of the probepoint 
    
        real(8) :: xd, yd, zd
        real(8) :: c00, c01, c10, c11
        real(8) :: C0, C1

        ! Compute normalized distance from c000 to probe point
        xd = (xp - x)/deltaX
        yd = (yp - y)/deltaY
        zd = (zp - z)/deltaZ

        ! Compute intermediate values in x-direction:
        C00 = C000*(1-xd) + C100*xd
        C01 = C001*(1-xd) + C101*xd
        C10 = C010*(1-xd) + C110*xd
        C11 = C011*(1-xd) + C111*xd

        ! Compute intermediate values in y-direction:
        C0  = C00*(1-yd) + C10*yd
        C1  = C01*(1-yd) + C11*yd

        ! Compute intermediate value in z-direction:
        probeValue = C0*(1-zd) + C1*zd

    end function IBM_probe_interp_trilinear




    function IBM_Cartesian_to_Spherical(xIB, yIB, zIB,  &
                                        xBC, yBC, zBC  )&
                                        result(SphericalVector)
        !!------------------------------------------------------------
        !!  **Purpose:** 
        !!      - Uses the IB coordinates (xIB, yIB, zIB) and the
        !!      bubble coordinates (xBC, yBC, zBC) to create spherical
        !!      vector: r-vec, theta-vec, phi-vec.
        !!
        !!  **Output:**
        !!      - Real(8) value of the spherical vectors stored in a
        !!      3-by-3 array:
        !!          | r_i       r_j       r_k     |
        !!          | theta_i   theta_j   theta_k |
        !!          | phi_i     phi_j     phi_k   |
        !!
        !!------------------------------------------------------------
        implicit none
        real(8), intent(in) :: xIB      ! IB-point x-coordinate
        real(8), intent(in) :: yIB      ! IB-point y-coordinate
        real(8), intent(in) :: zIB      ! IB-point z-coordinate
        real(8), intent(in) :: xBC      ! Bubble center x-coordinate
        real(8), intent(in) :: yBC      ! Bubble center y-coordinate
        real(8), intent(in) :: zBC      ! Bubble center z-coordinate
        
        ! array containing i,j,k for r, theta and phi direction
        real(8), dimension(3,3) :: SphericalVector 

        real(8) :: rdist,      phi,        theta,       &
                   rVec_i,     rVec_j,     rVec_k,      &
                   phiVec_i,   phiVec_j,   phiVec_k,    &
                   thetaVec_i, thetaVec_j, thetaVec_k
    
        ! Distance between IB-ghost and bubble center
        rdist   = sqrt((xIB - xBC)**2 + (yIB - yBC)**2 + (zIB - zBC)**2)

        ! Determine azimuth angle phi and inclination angle theta
        phi     = datan2((yIB - yBC), (xIB - xBC))
        theta   = acos((zIB - zBC)/rdist)

        ! Obtain radial (a.k.a. surface outward normal) direction
        rVec_i  = cos(phi)*sin(theta)
        rVec_j  = sin(phi)*sin(theta)
        rVec_k  = cos(theta)

        ! Surface tangential direction: inclination and azimuth
        phiVec_i  = -sin(phi)
        phiVec_j  = cos(phi)
        phiVec_k  = 0.0d0

        thetaVec_i  = cos(phi)*cos(theta)
        thetaVec_j  = sin(phi)*cos(theta)
        thetaVec_k  =-sin(theta)

        ! Store in output array:
        SphericalVector(1,1) = rVec_i
        SphericalVector(1,2) = rVec_j 
        SphericalVector(1,3) = rVec_k

        SphericalVector(2,1) = thetaVec_i
        SphericalVector(2,2) = thetaVec_j
        SphericalVector(2,3) = thetaVec_k

        SphericalVector(3,1) = phiVec_i
        SphericalVector(3,2) = phiVec_j
        SphericalVector(3,3) = phiVec_k

    end function IBM_Cartesian_to_Spherical

end module IBM_functions