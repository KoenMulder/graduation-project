module var2
    !!------------------------------------------------------------------------------
    !! TU/e ET Graduation Project 
    !!------------------------------------------------------------------------------
    !!
    !! MODULE:  var2
    !!
    !!> @author
    !!> Koen Mulder}
    !!
    !! DESCRIPTION: 
    !!>  Contains all the parameters and initial conditions for the
    !!>  the simulation.
    !!
    !! REVISION HISTORY:
    !! dd Mmm yyyy - Initial Version
    !! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
    !!------------------------------------------------------------------------------
    implicit none
    save
    integer, parameter ::       &
                im = 48,        & ! Number of cells in x-direction
                jm = 96,        & ! Number of cells in y-direction
                km = 48,        & ! Number of cells in z-direction
                nThreads = 24,  & ! Number of threads proportional to half the grid
                m1_max = 2000,  & ! Number of m2-iterations
                m2_max = 250,   & ! Saving frequency 
                nmax = 3,       & ! Number of concentration species
                nBubbles = 1      ! Number of bubbles in the simulation

    logical, parameter ::                       &
                 bubble             = .true.,   & ! Toggle simulation of bubble
                 VirtForcMethod     = .true.,   & ! Toggle virtual force method 
                 DragForceOutput    = .false.,  & ! Toggle force debugging ouput files
                 WriteToBinaryOutput= .true.,  & ! Toggle unformatted binary output file
                 FullDomainOutput   = .false.,  & ! Toggle scalar and velocity output files
                 LoadBinaryData     = .false.,  & ! Toggle loading previous simulation data
                 useRungeKutta      = .false.      ! Toggle Fourth order low storage Runge-Kutta

    integer, parameter ::                       &
                LoadBinaryDataM2iter= 8           ! When LoadBinaryData is toggled on specify
                                                  ! the previous simulation m2-iteration

    real(8), parameter ::               &
                ! Domain properties
                lx1     = 1d-4,         & ! width   [m]
                lx2     = 2d-4,         & ! height  [m]
                lx3     = 1d-4,         & ! depth   [m]
                ! Fluid and bubble properties   
                rhoKOH  = 1258d0,       & ! fluid density       [kg/m3]
                nuKOH   = 6.7d-7,       & ! kinmatic viscosity  [m2/s]
                rhoH2   = 0.0686d0,     & ! bubble density      [kg/m3]
                DH2     = 5.8d-9,       & ! diffusivity         [m2/s]
                DKOH    = 3.2d-9,       & ! diffusivity         [m2/s]
                DH2O    = 3.2d-9,       & ! diffusivity         [m2/s]
                C0H2    = 0.16d0,       & ! ref. concentration  [mol/m3]
                C0KOH   = 6700d0,       & ! ref. concentration  [mol/m3]
                C0H20   = 49000d0,      & ! ref. concentration  [mol/m3]
                umax    = 0d0,          & ! maximum velocity    [m/s]
                R_max   = 0.25d-4,      & ! detachment radius   [m]
                sigma   = 7.7d-2,       & ! surface tension     [N/m]
                ! Electrode and Membrame properties
                phiLe   = -0.45d0,      & ! phi_electr. + E^0   [V]
                phiR    = 0d0,          & ! potential Membrame  [V]
                j0      = 1d0,          & ! exchnge curnt. dens.[A/m2]
                kc0     = -j0,          & ! BV-kinetic constant [A/m2]
                ka0     = j0,           & ! BV-kinetic constant [A/m2]
                ! Universal constants
                Rid     = 8.3142d0,     & ! ideal gas constant  [J/K mol]
                Far     = 96485d0,      & ! Faraday's constant  [A/mol]
                temp    = 353d0,        & ! temperature         [K]
                pi      = 3.141592653d0,& ! pi                  [-]
                ! Direction: downwards= -g,  upwards= +g, 
                g       = -9.81d0,      & ! Gravitational accel.[m/s2]
                ! Simulation properties
                Co      = 1d0/3d0,      & ! Courant number      [-]
                Cv      = 1.0d0,        & ! virtual force coef. [-]
                ACmod   = 1.0d-2           ! AC modifyer coef.   [-]


    character(len=10), parameter ::         &
                timeDiscrScheme  = 'FOBD',  &   ! Select modified time integration scheme 
                                                ! for bubble equation of motion with 
                                                ! Virtual force method:
                                                ! 'FOBD' = First Order Backwards (O1)
                                                ! 'SOBD' = Second Order Backwards (O2)

                artCompresMethod =  'auto', &   ! Select between modified or automatic
                                                ! calculation of the Artificial 
                                                ! Compressibility (AC) parameter:
                                                ! 'modified' = add a modifier to the 
                                                !              automatically selected time step  
                                                !              and universal speed of sound
                                                ! 'auto'     = automatically select temporal
                                                !              appropriate conditions

                FlowCondPreset   = 'default'     ! Select boundary and flow conditions based
                                                ! based on the preset:
                                                ! 'Stokes'   = periodic BC at all domain 
                                                !              boundaries. Uniform flow of
                                                !              zero velocity (quiescent).
                                                !              initial bubble position at domain
                                                !              center
                                                ! 'free-slip'= free-slip BC at walls (electrode 
                                                !               + membrame). No flow (quiescent)
                                                !              initial bubble position at domain
                                                !              center
                                                ! 'default'  = electrode configuration, no-slip
                                                !              at the walls, parabolic flow 
                                                !              profile at inlet, bubble initial
                                                !              position is at the electrode.

    real(8) ::  c_w, time, dtime, usou2 ,beta, total,       &
                ac, aa, condfac, radius0, Max_Err_phi,      &
                balance_h2

    ! Bubble variables
    real(8) ::  xc, yc, zc, dxcdt, dycdt, dzcdt, drdt,      &
                radius, FxCOM, FyCOM, FzCOM, bubbleMassflux,&
                fvx, fvy, fvz, fsx, fsy, fsz, fby,          &
                ddxcdttFv, ddycdttFv, ddzcdttFv,            &
                dycdt_old,  dzcdt_old, dxcdt_old,           &
                dxcdt_oldold, dycdt_oldold, dzcdt_oldold

    integer :: cube_number, i, j, k, change, iter1, iter2,  &
               iter3, iter4, o_u, o_v, o_w, m1, m2,         &
               nx_start, rungIter, rungIterEnd

    logical :: isDetached
    
    ! Rank 1 arrays
    real(8), allocatable, dimension(:) :: x1, x1_i, dx1, dx1_i,  &
                                          x2, x2_j, dx2, dx2_j,  &
                                          x3, x3_k, dx3, dx3_k,  &
                                          c_s, dif,         &
                                          xb, yb, zb, c_ref,&
                                          phiS, thetaS,     &
                                          T_ni, T_nj, T_nk, &
                                          T_ti, T_tj, T_tk, &
                                          T_si, T_sj, T_sk, &
                                          rb, xb_c, yb_c,   &
                                          zb_c,             &
                                          phiS_U, thetaS_U, &
                                          phiS_V, thetaS_V, &
                                          phiS_W, thetaS_W, &

                                          T_ni_U, T_nj_U, T_nk_U,   &
                                          T_ti_U, T_tj_U, T_tk_U,   &
                                          T_si_U, T_sj_U,  T_sk_U,  &
                                          
                                          T_ni_V, T_nj_V, T_nk_V,   &
                                          T_ti_V, T_tj_V, T_tk_V,   &
                                          T_si_V, T_sj_V, T_sk_V,   &

                                          T_ni_W, T_nj_W, T_nk_W,   &
                                          T_ti_W, T_tj_W, T_tk_W,   &
                                          T_si_W, T_sj_W, T_sk_W,   &

                                          rb_U, xb_U, yb_U, zb_U,   &
                                          rb_V, xb_v, yb_V, zb_V,   &
                                          rb_W, xb_W, yb_W, zb_W,   &

                                          RK_alpha

    integer, allocatable, dimension(:) ::   &
                                        x_IB, y_IB, z_IB,       &
                                        direction,              &
                                        x1_plus, x1_minus,      &
                                        x2_plus, x2_minus,      &
                                        x3_plus, x3_minus,      &
                                        x_IB_U, y_IB_U, z_IB_U, &
                                        x_IB_V, y_IB_V, z_IB_V, &
                                        x_IB_W, y_IB_W, z_IB_W


    ! Rank 2 arrays
    real(8), allocatable, dimension(:,:) :: &
                                x_prob, y_prob, z_prob,         &
                                x_prob_c, y_prob_c, z_prob_c,   &
                                x_prob_U, y_prob_U, z_prob_U,   &
                                x_prob_V, y_prob_V, z_prob_V,   &
                                x_prob_W, y_prob_W, z_prob_W,   &
                                eta, phiL, ka, kc,              &
                                ! Domain flux (force) arrays
                                A11, A12, A13,   B11, B12, B13, &
                                A21, A22, A23,   B21, B22, B23, &
                                A31, A32, A33,   B31, B32, B33, &
                                u_i, u_o, u_f, u_b, w_i, w_o,   &
                                v_f, v_b

    integer, allocatable, dimension(:,:) ::   &
                                neighbor, &
                                nx_prob, ny_prob, nz_prob,          &
                                nx_prob_c, ny_prob_c, nz_prob_c,    &

                                nx_prob_UU, ny_prob_UU, nz_prob_UU, &
                                nx_prob_UV, ny_prob_UV, nz_prob_UV, &
                                nx_prob_UW, ny_prob_UW, nz_prob_UW, &

                                nx_prob_VU, ny_prob_VU, nz_prob_VU, &
                                nx_prob_VV, ny_prob_VV, nz_prob_VV, &
                                nx_prob_VW, ny_prob_VW, nz_prob_VW, &

                                nx_prob_WU, ny_prob_WU, nz_prob_WU, &
                                nx_prob_WV, ny_prob_WV, nz_prob_WV, &
                                nx_prob_WW, ny_prob_WW, nz_prob_WW


    ! Rank 3 arrays
    real(8), allocatable, dimension(:,:,:) ::   &
                                        u1, u2, u3, p, rhu, rhv, &
                                        rhw, div, u_node, v_node,&
                                        w_node, p_node, work1,   &
                                        work2, work3, rhs,       &
                                        coef1_phi, phi,          &
                                        coef2_phi, coef3_phi,    &
                                        coef4_phi, coef5_phi,    &
                                        coef6_phi, current,      &
                                        conduct,                 &
                                        storage_u, storage_v,    &
                                        storage_w, storage_p,    &
                                        Error_phi, storage_phi,  &
                                        u1_temp, u2_temp, u3_temp
                                        
    integer, allocatable, dimension(:,:,:) ::   &
                            typ, typ1, typ2, typ3, typ1a, typ2a, typ3a, &
                            typ_IB, typ_u_IB, typ_v_IB, typ_w_IB,       &
                            typ_force_x, typ_force_y, typ_force_z

    ! Rank 4 arrays
    real(8), allocatable, dimension(:,:,:,:) ::   &
                                        c, storage_c, c_old, Error_c


    type :: bubble_properties
        real(8) :: xc, yc, zc, radius
        real(8) :: dxcdt, dycdt, dzcdt, drdt
    end type

contains
    subroutine chckProbeBounds(probeIdx,minIdx,maxIdx,msg)
        integer, intent(in)           :: probeIdx, minIdx, maxIdx
        character(len=*), intent(in) :: msg

        if ( (probeIdx.lt.minIdx).or.(probeIdx.gt.maxIdx) ) then
            print*, 'Probe index out of bounds', msg, ' with probeIdx ', probeIdx
            stop
        end if
    end subroutine chckProbeBounds
    
end module var2