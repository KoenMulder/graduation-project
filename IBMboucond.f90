subroutine IBMboucond()
    !!------------------------------------------------------------
    !!  **Main subroutine for IBM boundary conditions**
    !!  
    !!  Requires *IBMidentify* to be run first.
    !!
    !!  Only runs when *simBubble = .true.*
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: n, nIB, nFace

    if ( SimBubble ) then

        ! Set all ghost cells velocity to zero
        call IBM_velocity_ghost

        ! Store current values for interpolation
        store_u = u1
        store_v = u2
        store_w = u3

        do n = 1, nBubbles

            ! Concentration boundary conditions
            do nIB = 1,BubbleBlock(n)%nIB
                call IBM_concentration(nIB,n)
            end do

            ! Surface velocity boundary conditions
            do nFace = 1,BubbleBlock(n)%nFacesX
                call IBM_velocity_u(nFace, n)
            end do

            do nFace = 1,BubbleBlock(n)%nFacesY
                call IBM_velocity_v(nFace, n)
            end do

            do nFace = 1,BubbleBlock(n)%nFacesZ
                call IBM_velocity_w(nFace, n)
            end do
        end do

        ! Save interpolated values
        u1 = store_u
        u2 = store_v
        u3 = store_w
    end if
end subroutine IBMboucond




subroutine IBM_velocity_ghost()
    !!------------------------------------------------------------
    !!  **Purpose:**
    !!
    !!  Applies zero velocity to faces of ghost cells in u1, u2
    !!  and u3 direction.
    !!
    !!  Assumes forward staggered grid.
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k
    
    do i = 1,im
        do j = 1,jm
            do k = 1,km
                if ( IBtype(i,j,k).eq.ghost_label ) then
                    u1(i-1,j,k) = 0d0
                    u1(i,j,k)   = 0d0
                    u2(i,j-1,k) = 0d0
                    u2(i,j,k)   = 0d0
                    u3(i,j,k-1) = 0d0
                    u3(i,j,k)   = 0d0
                end if
            end do
        end do
    end do

end subroutine IBM_velocity_ghost




subroutine IBM_concentration(nIB, bubbleID)
    !!------------------------------------------------------------
    !!  **Purpose:**
    !!
    !!  Applies zero Neumann boundary condition on the surface
    !!  of the bubble for KOH and H2O. 
    !!
    !!  Applies a Dirichlet boundary condition on the surface 
    !!  of the bubble for H2.   
    !!
    !!  Uses the index location of the IB cells from 
    !!  IBcellLocation and performs trilinear interpolation to
    !!  obtain the value for the IB cell for a given boundary
    !!  condition.
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    use IBM_functions
    implicit none
    integer, intent(in) :: nIB, bubbleID
    real(8), dimension(4) :: probeXYZsdist
    integer :: id, jd, kd, c000I, c000J, c000K, nConc
    real(8) ::  c000, c001, c011, c010, &
                c110, c100, c101, c111, &
                beta
    
    ! IBcellLocation(1,nIB,nBubbles) = i-index of IB cell for bubble n
    ! IBcellLocation(2,nIB,nBubbles) = j-index of IB cell for bubble n
    ! IBcellLocation(3,nIB,nBubbles) = k-index of IB cell for bubble n
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

    ! Perform trilinear interpolation cH2, cKOH, cH2O
    do nConc = 1, nspec
        c000 = c(c000I,     c000J,     c000K,     nConc)
        c001 = c(c000I,     c000J,     c000K + 1, nConc)
        c011 = c(c000I,     c000J + 1, c000K + 1, nConc)
        c010 = c(c000I,     c000J + 1, c000K,     nConc)
        c110 = c(c000I + 1, c000J + 1, c000K,     nConc)
        c100 = c(c000I + 1, c000J,     c000K,     nConc)
        c101 = c(c000I + 1, c000J,     c000K + 1, nConc)
        c111 = c(c000I + 1, c000J + 1, c000K + 1, nConc)

        c_s(nConc) = IBM_probe_interp_trilinear(            &
                        c000, c001, c011, c010,             &
                        c110, c100, c101, c111,             &
                        dx1_i(id), dx2_j(jd), dx3_k(kd),    &
                        probeXYZsdist(1), probeXYZsdist(2), &
                        probeXYZsdist(3),                   &
                        x1(c000I), x2(c000J), x3(c000K)     )
    end do

    ! Dirichlet for cH2
    beta = probeXYZsdist(4)/dx1(id)
    c(id, jd, kd, 1) = (1 + beta)*c_ref(1) - beta*c_s(1)
    
    ! Neumann for cKOH and cH2O (flux = 0)
    c(id, jd, kd, 2) = c_s(2)
    c(id, jd, kd, 3) = c_s(3)

end subroutine IBM_concentration




subroutine IBM_velocity_u(nFace, bubbleID )
    !!------------------------------------------------------------
    !!  **Purpose:**
    !!
    !!  Applies a Dirichlet boundary condition on the surface 
    !!  of the bubble for the u-velocity component.   
    !!
    !!  Uses the index location of the IB cells from 
    !!  IBcellLocation and performs trilinear interpolation to
    !!  obtain the value for the IB cell for a given boundary
    !!  condition.
    !!
    !!------------------------------------------------------------
    use commondata
    use IBM_functions
    implicit none
    integer, intent(in) :: nFace, bubbleID
    real(8), dimension(4) :: probeXYZsdist
    real(8), dimension(3,3) :: spherical_coord
    integer ::  id, jd, kd,                 &
                c000UUi, c000UUj, c000UUk,  &
                c000UVi, c000UVj, c000UVk,  &
                c000UWi, c000UWj, c000UWk
    real(8) ::  c000, c001, c011, c010,     &
                c110, c100, c101, c111,     &
                beta, uuVal, uvVal, uwVal,  &
                uni_ref, uti_ref, usi_ref,  &
                u_star, v_star, w_star,     &
                uni_star, uti_star, usi_star, &
                unIB, utIB, usIB 

    ! Get IB-u face index location
    id = IBfaceXLocation(1,nFace,bubbleID)
    jd = IBfaceXLocation(2,nFace,bubbleID)
    kd = IBfaceXLocation(3,nFace,bubbleID)

    ! Get probe point location
    probeXYZsdist = IBM_get_probe_location(                 &
                            x1_i(id),  x2(jd),    x3(kd),   &
                            dx1_i(id), dx2_j(jd), dx3_k(kd),&
                            BubbleBlock(bubbleID)%xCenter,  &
                            BubbleBlock(bubbleID)%yCenter,  &
                            BubbleBlock(bubbleID)%zCenter,  &
                            BubbleBlock(bubbleID)%radius    )

    ! Store probe location for future use or debugging
    probeU(1,nFace,bubbleID) = probeXYZsdist(1)
    probeU(2,nFace,bubbleID) = probeXYZsdist(2)
    probeU(3,nFace,bubbleID) = probeXYZsdist(3)

    ! Get c000 location of probes point for UU, UV, and UW
    c000UUi = floor(probeXYZsdist(1)/dx1_i(id))
    c000UUj = floor(probeXYZsdist(2)/dx2_j(jd) + 0.5d0)
    c000UUk = floor(probeXYZsdist(3)/dx3_k(kd) + 0.5d0)

    c000UVi = floor(probeXYZsdist(1)/dx1_i(id) + 0.5d0)
    c000UVj = floor(probeXYZsdist(2)/dx2_j(jd))
    c000UVk = floor(probeXYZsdist(3)/dx3_k(kd) + 0.5d0)

    c000UWi = floor(probeXYZsdist(1)/dx1_i(id) + 0.5d0)
    c000UWj = floor(probeXYZsdist(2)/dx2_j(jd) + 0.5d0)
    c000UWk = floor(probeXYZsdist(3)/dx3_k(kd))

    ! Perform trilinear interpolation on UU, UV and UW
    ! x-component of u1
    c000 = u1(c000UUi,     c000UUj,     c000UUk    )
    c001 = u1(c000UUi,     c000UUj,     c000UUk + 1)
    c011 = u1(c000UUi,     c000UUj + 1, c000UUk + 1)
    c010 = u1(c000UUi,     c000UUj + 1, c000UUk    )
    c110 = u1(c000UUi + 1, c000UUj + 1, c000UUk    )
    c100 = u1(c000UUi + 1, c000UUj,     c000UUk    )
    c101 = u1(c000UUi + 1, c000UUj,     c000UUk + 1)
    c111 = u1(c000UUi + 1, c000UUj + 1, c000UUk + 1)

    uuVal = IBM_probe_interp_trilinear(                     &
                        c000, c001, c011, c010,             &
                        c110, c100, c101, c111,             &
                        dx1(id), dx2_j(jd), dx3_k(kd),      &
                        probeXYZsdist(1), probeXYZsdist(2), &
                        probeXYZsdist(3),                   &
                        x1_i(c000UUi), x2(c000UUj),         &
                        x3(c000UUk)                         )

    ! y-component of u1
    c000 = u1(c000UVi,     c000UVj,     c000UVk    )
    c001 = u1(c000UVi,     c000UVj,     c000UVk + 1)
    c011 = u1(c000UVi,     c000UVj + 1, c000UVk + 1)
    c010 = u1(c000UVi,     c000UVj + 1, c000UVk    )
    c110 = u1(c000UVi + 1, c000UVj + 1, c000UVk    )
    c100 = u1(c000UVi + 1, c000UVj,     c000UVk    )
    c101 = u1(c000UVi + 1, c000UVj,     c000UVk + 1)
    c111 = u1(c000UVi + 1, c000UVj + 1, c000UVk + 1)

    uvVal = IBM_probe_interp_trilinear(                     &
                        c000, c001, c011, c010,             &
                        c110, c100, c101, c111,             &
                        dx1_i(id), dx2(jd), dx3_k(kd),      &
                        probeXYZsdist(1), probeXYZsdist(2), &
                        probeXYZsdist(3),                   &
                        x1(c000UVi), x2_j(c000UVj),         &
                        x3(c000UVk)                         )

    ! z-component of u1
    c000 = u1(c000UWi,     c000UWj,     c000UWk    )
    c001 = u1(c000UWi,     c000UWj,     c000UWk + 1)
    c011 = u1(c000UWi,     c000UWj + 1, c000UWk + 1)
    c010 = u1(c000UWi,     c000UWj + 1, c000UWk    )
    c110 = u1(c000UWi + 1, c000UWj + 1, c000UWk    )
    c100 = u1(c000UWi + 1, c000UWj,     c000UWk    )
    c101 = u1(c000UWi + 1, c000UWj,     c000UWk + 1)
    c111 = u1(c000UWi + 1, c000UWj + 1, c000UWk + 1)

    uwVal = IBM_probe_interp_trilinear(                     &
                        c000, c001, c011, c010,             &
                        c110, c100, c101, c111,             &
                        dx1_i(id), dx2_j(jd), dx3(kd),      &
                        probeXYZsdist(1), probeXYZsdist(2), &
                        probeXYZsdist(3),                   &
                        x1(c000UWi), x2(c000UWj),           &
                        x3_k(c000UWk)                       )

    ! Translate cartesian to sperhical coordinates for
    ! probe point velocity u1.
    spherical_coord = IBM_Cartesian_to_Spherical(           &
                        x1_i(id), x2(jd), x3(kd),           &
                        BubbleBlock(bubbleID)%xCenter,      &
                        BubbleBlock(bubbleID)%yCenter,      &
                        BubbleBlock(bubbleID)%zCenter       )
    ! radial/surface normal direction
    uni_ref = (uuVal - BubbleBlock(bubbleID)%u)*spherical_coord(1,1) +  &
              (uvVal - BubbleBlock(bubbleID)%v)*spherical_coord(1,2) +  &
              (uwVal - BubbleBlock(bubbleID)%w)*spherical_coord(1,3)
    ! inclination/surface tangential direction
    uti_ref = (uuVal - BubbleBlock(bubbleID)%u)*spherical_coord(2,1) +  &
              (uvVal - BubbleBlock(bubbleID)%v)*spherical_coord(2,2) +  &
              (uwVal - BubbleBlock(bubbleID)%w)*spherical_coord(2,3)
    ! azimuth/surface tangential direction
    usi_ref = (uuVal - BubbleBlock(bubbleID)%u)*spherical_coord(3,1) +  &
              (uvVal - BubbleBlock(bubbleID)%v)*spherical_coord(3,2) +  &
              (uwVal - BubbleBlock(bubbleID)%w)*spherical_coord(3,3)

    ! Translate cartesian to spherical coordinates for IB-u1 velocity
    u_star = BubbleBlock(bubbleID)%u + BubbleBlock(bubbleID)%drdt*spherical_coord(1,1)
    v_star = BubbleBlock(bubbleID)%v + BubbleBlock(bubbleID)%drdt*spherical_coord(1,2)
    w_star = BubbleBlock(bubbleID)%w + BubbleBlock(bubbleID)%drdt*spherical_coord(1,3)

    ! radial/surface normal direction
    uni_star = (u_star - BubbleBlock(bubbleID)%u)*spherical_coord(1,1) +  &
               (v_star - BubbleBlock(bubbleID)%v)*spherical_coord(1,2) +  &
               (w_star - BubbleBlock(bubbleID)%w)*spherical_coord(1,3)
    ! inclination/surface tangential direction
    uti_star = (u_star - BubbleBlock(bubbleID)%u)*spherical_coord(2,1) +  &
               (v_star - BubbleBlock(bubbleID)%v)*spherical_coord(2,2) +  &
               (w_star - BubbleBlock(bubbleID)%w)*spherical_coord(2,3)
    ! azimuth/surface tangential direction
    usi_star = (u_star - BubbleBlock(bubbleID)%u)*spherical_coord(3,1) +  &
               (v_star - BubbleBlock(bubbleID)%v)*spherical_coord(3,2) +  &
               (w_star - BubbleBlock(bubbleID)%w)*spherical_coord(3,3)


    ! Dirichlet Boundary condition for normal/radial direction by means of 
    ! linear interpolation.
    beta = probeXYZsdist(4)/dx1(id)
    unIB = (1 + beta)*uni_star - beta*uni_ref
    utIB = (1 + beta)*uti_star - beta*uti_ref
    usIB = (1 + beta)*usi_star - beta*usi_ref

    ! Translate spherical boundary condition back to cartesian expression for
    ! the IB-point (only take bubble u-velocity)
    store_u(id, jd, kd) = unIB*spherical_coord(1,1) + &
                          utIB*spherical_coord(1,2) + &
                          usIB*spherical_coord(1,3) + &
                          BubbleBlock(bubbleID)%u 
    
end subroutine IBM_velocity_u



subroutine IBM_velocity_v(nFace, bubbleID)
    !!------------------------------------------------------------
    !!  **Purpose:**
    !!
    !!  Applies a Dirichlet boundary condition on the surface 
    !!  of the bubble for the v-velocity component.   
    !!
    !!  Uses the index location of the IB cells from 
    !!  IBcellLocation and performs trilinear interpolation to
    !!  obtain the value for the IB cell for a given boundary
    !!  condition.
    !!
    !!------------------------------------------------------------
    use commondata
    use IBM_functions
    implicit none
    integer, intent(in) :: nFace, bubbleID
    real(8), dimension(4) :: probeXYZsdist
    real(8), dimension(3,3) :: spherical_coord
    integer ::  id, jd, kd,                 &
                c000VUi, c000VUj, c000VUk,  &
                c000VVi, c000VVj, c000VVk,  &
                c000VWi, c000VWj, c000VWk
    real(8) ::  c000, c001, c011, c010,     &
                c110, c100, c101, c111,     &
                beta, vuVal, vvVal, vwVal,  &
                vni_ref, vti_ref, vsi_ref,  &
                u_star, v_star, w_star,     &
                vni_star, vti_star, vsi_star, &
                vnIB, vtIB, vsIB 

    ! Get IB-u face index location
    id = IBfaceYLocation(1,nFace,bubbleID)
    jd = IBfaceYLocation(2,nFace,bubbleID)
    kd = IBfaceYLocation(3,nFace,bubbleID)

    ! Get probe point location
    probeXYZsdist = IBM_get_probe_location(                 &
                            x1(id),    x2_j(jd),  x3(kd),   &
                            dx1_i(id), dx2_j(jd), dx3_k(kd),&
                            BubbleBlock(bubbleID)%xCenter,  &
                            BubbleBlock(bubbleID)%yCenter,  &
                            BubbleBlock(bubbleID)%zCenter,  &
                            BubbleBlock(bubbleID)%radius    )

    ! Store probe location for future use or debugging
    probeV(1,nFace,bubbleID) = probeXYZsdist(1)
    probeV(2,nFace,bubbleID) = probeXYZsdist(2)
    probeV(3,nFace,bubbleID) = probeXYZsdist(3)

    ! Get c000 location of probes point for VU, VV, and VW
    c000VUi = floor(probeXYZsdist(1)/dx1_i(id))
    c000VUj = floor(probeXYZsdist(2)/dx2_j(jd) + 0.5d0)
    c000VUk = floor(probeXYZsdist(3)/dx3_k(kd) + 0.5d0)

    c000VVi = floor(probeXYZsdist(1)/dx1_i(id) + 0.5d0)
    c000VVj = floor(probeXYZsdist(2)/dx2_j(jd))
    c000VVk = floor(probeXYZsdist(3)/dx3_k(kd) + 0.5d0)

    c000VWi = floor(probeXYZsdist(1)/dx1_i(id) + 0.5d0)
    c000VWj = floor(probeXYZsdist(2)/dx2_j(jd) + 0.5d0)
    c000VWk = floor(probeXYZsdist(3)/dx3_k(kd))

    ! Perform trilinear interpolation on VU, VV and VW
    ! x-component of u2
    c000 = u2(c000VUi,     c000VUj,     c000VUk    )
    c001 = u2(c000VUi,     c000VUj,     c000VUk + 1)
    c011 = u2(c000VUi,     c000VUj + 1, c000VUk + 1)
    c010 = u2(c000VUi,     c000VUj + 1, c000VUk    )
    c110 = u2(c000VUi + 1, c000VUj + 1, c000VUk    )
    c100 = u2(c000VUi + 1, c000VUj,     c000VUk    )
    c101 = u2(c000VUi + 1, c000VUj,     c000VUk + 1)
    c111 = u2(c000VUi + 1, c000VUj + 1, c000VUk + 1)

    vuVal = IBM_probe_interp_trilinear(                     &
                        c000, c001, c011, c010,             &
                        c110, c100, c101, c111,             &
                        dx1(id), dx2_j(jd), dx3_k(kd),      &
                        probeXYZsdist(1), probeXYZsdist(2), &
                        probeXYZsdist(3),                   &
                        x1_i(c000VUi), x2(c000VUj),         &
                        x3(c000VUk)                         )

    ! y-component of u2
    c000 = u2(c000VVi,     c000VVj,     c000VVk    )
    c001 = u2(c000VVi,     c000VVj,     c000VVk + 1)
    c011 = u2(c000VVi,     c000VVj + 1, c000VVk + 1)
    c010 = u2(c000VVi,     c000VVj + 1, c000VVk    )
    c110 = u2(c000VVi + 1, c000VVj + 1, c000VVk    )
    c100 = u2(c000VVi + 1, c000VVj,     c000VVk    )
    c101 = u2(c000VVi + 1, c000VVj,     c000VVk + 1)
    c111 = u2(c000VVi + 1, c000VVj + 1, c000VVk + 1)

    vvVal = IBM_probe_interp_trilinear(                     &
                        c000, c001, c011, c010,             &
                        c110, c100, c101, c111,             &
                        dx1_i(id), dx2(jd), dx3_k(kd),      &
                        probeXYZsdist(1), probeXYZsdist(2), &
                        probeXYZsdist(3),                   &
                        x1(c000VVi), x2_j(c000VVj),         &
                        x3(c000VVk)                         )

    ! z-component of u2
    c000 = u2(c000VWi,     c000VWj,     c000VWk    )
    c001 = u2(c000VWi,     c000VWj,     c000VWk + 1)
    c011 = u2(c000VWi,     c000VWj + 1, c000VWk + 1)
    c010 = u2(c000VWi,     c000VWj + 1, c000VWk    )
    c110 = u2(c000VWi + 1, c000VWj + 1, c000VWk    )
    c100 = u2(c000VWi + 1, c000VWj,     c000VWk    )
    c101 = u2(c000VWi + 1, c000VWj,     c000VWk + 1)
    c111 = u2(c000VWi + 1, c000VWj + 1, c000VWk + 1)

    vwVal = IBM_probe_interp_trilinear(                     &
                        c000, c001, c011, c010,             &
                        c110, c100, c101, c111,             &
                        dx1_i(id), dx2_j(jd), dx3(kd),      &
                        probeXYZsdist(1), probeXYZsdist(2), &
                        probeXYZsdist(3),                   &
                        x1(c000VWi), x2(c000VWj),           &
                        x3_k(c000VWk)                       )

    ! Translate cartesian to sperhical coordinates for
    ! probe point velocity u2.
    spherical_coord = IBM_Cartesian_to_Spherical(           &
                        x1(id), x2_j(jd), x3(kd),           &
                        BubbleBlock(bubbleID)%xCenter,      &
                        BubbleBlock(bubbleID)%yCenter,      &
                        BubbleBlock(bubbleID)%zCenter       )

    ! radial/surface normal direction
    vni_ref = (vuVal - BubbleBlock(bubbleID)%u)*spherical_coord(1,1) +  &
              (vvVal - BubbleBlock(bubbleID)%v)*spherical_coord(1,2) +  &
              (vwVal - BubbleBlock(bubbleID)%w)*spherical_coord(1,3)
    ! inclination/surface tangential direction
    vti_ref = (vuVal - BubbleBlock(bubbleID)%u)*spherical_coord(2,1) +  &
              (vvVal - BubbleBlock(bubbleID)%v)*spherical_coord(2,2) +  &
              (vwVal - BubbleBlock(bubbleID)%w)*spherical_coord(2,3)
    ! azimuth/surface tangential direction
    vsi_ref = (vuVal - BubbleBlock(bubbleID)%u)*spherical_coord(3,1) +  &
              (vvVal - BubbleBlock(bubbleID)%v)*spherical_coord(3,2) +  &
              (vwVal - BubbleBlock(bubbleID)%w)*spherical_coord(3,3)

    ! Translate cartesian to spherical coordinates for IB-u1 velocity
    u_star = BubbleBlock(bubbleID)%u + BubbleBlock(bubbleID)%drdt*spherical_coord(1,1)
    v_star = BubbleBlock(bubbleID)%v + BubbleBlock(bubbleID)%drdt*spherical_coord(1,2)
    w_star = BubbleBlock(bubbleID)%w + BubbleBlock(bubbleID)%drdt*spherical_coord(1,3)

    ! radial/surface normal direction
    vni_star = (u_star - BubbleBlock(bubbleID)%u)*spherical_coord(1,1) +  &
               (v_star - BubbleBlock(bubbleID)%v)*spherical_coord(1,2) +  &
               (w_star - BubbleBlock(bubbleID)%w)*spherical_coord(1,3)
    ! inclination/surface tangential direction
    vti_star = (u_star - BubbleBlock(bubbleID)%u)*spherical_coord(2,1) +  &
               (v_star - BubbleBlock(bubbleID)%v)*spherical_coord(2,2) +  &
               (w_star - BubbleBlock(bubbleID)%w)*spherical_coord(2,3)
    ! azimuth/surface tangential direction
    vsi_star = (u_star - BubbleBlock(bubbleID)%u)*spherical_coord(3,1) +  &
               (v_star - BubbleBlock(bubbleID)%v)*spherical_coord(3,2) +  &
               (w_star - BubbleBlock(bubbleID)%w)*spherical_coord(3,3)

    ! Dirichlet Boundary condition for normal/radial direction by means of 
    ! linear interpolation.
    beta = probeXYZsdist(4)/dx2(jd)
    vnIB = (1 + beta)*vni_star - beta*vni_ref
    vtIB = (1 + beta)*vti_star - beta*vti_ref
    vsIB = (1 + beta)*vsi_star - beta*vsi_ref

    ! Translate spherical boundary condition back to cartesian expression for
    ! the IB-point (only take bubble u-velocity)
    store_v(id, jd, kd) = vnIB*spherical_coord(1,1) + &
                          vtIB*spherical_coord(1,2) + &
                          vsIB*spherical_coord(1,3) + &
                          BubbleBlock(bubbleID)%v
    
end subroutine IBM_velocity_v




subroutine IBM_velocity_w(nFace, bubbleID)
    !!------------------------------------------------------------
    !!  **Purpose:**
    !!
    !!  Applies a Dirichlet boundary condition on the surface 
    !!  of the bubble for the w-velocity component.   
    !!
    !!  Uses the index location of the IB cells from 
    !!  IBcellLocation and performs trilinear interpolation to
    !!  obtain the value for the IB cell for a given boundary
    !!  condition.
    !!
    !!------------------------------------------------------------
    use commondata
    use IBM_functions
    implicit none
    integer, intent(in) :: nFace, bubbleID
    real(8), dimension(4) :: probeXYZsdist
    real(8), dimension(3,3) :: spherical_coord
    integer ::  id, jd, kd,                 &
                c000WUi, c000WUj, c000WUk,  &
                c000WVi, c000WVj, c000WVk,  &
                c000WWi, c000WWj, c000WWk
    real(8) ::  c000, c001, c011, c010,     &
                c110, c100, c101, c111,     &
                beta, wuVal, wvVal, wwVal,  &
                wni_ref, wti_ref, wsi_ref,  &
                u_star, v_star, w_star,     &
                wni_star, wti_star, wsi_star, &
                wnIB, wtIB, wsIB 

    ! Get IB-u face index location
    id = IBfaceZLocation(1,nFace,bubbleID)
    jd = IBfaceZLocation(2,nFace,bubbleID)
    kd = IBfaceZLocation(3,nFace,bubbleID)

    ! Get probe point location
    probeXYZsdist = IBM_get_probe_location(                 &
                            x1(id),    x2(jd),  x3_k(kd),   &
                            dx1_i(id), dx2_j(jd), dx3_k(kd),&
                            BubbleBlock(bubbleID)%xCenter,  &
                            BubbleBlock(bubbleID)%yCenter,  &
                            BubbleBlock(bubbleID)%zCenter,  &
                            BubbleBlock(bubbleID)%radius    )

    ! Store probe location for future use or debugging
    probeW(1,nFace,bubbleID) = probeXYZsdist(1)
    probeW(2,nFace,bubbleID) = probeXYZsdist(2)
    probeW(3,nFace,bubbleID) = probeXYZsdist(3)

    ! Get c000 location of probes point for VU, VV, and VW
    c000WUi = floor(probeXYZsdist(1)/dx1_i(id))
    c000WUj = floor(probeXYZsdist(2)/dx2_j(jd) + 0.5d0)
    c000WUk = floor(probeXYZsdist(3)/dx3_k(kd) + 0.5d0)

    c000WVi = floor(probeXYZsdist(1)/dx1_i(id) + 0.5d0)
    c000WVj = floor(probeXYZsdist(2)/dx2_j(jd))
    c000WVk = floor(probeXYZsdist(3)/dx3_k(kd) + 0.5d0)

    c000WWi = floor(probeXYZsdist(1)/dx1_i(id) + 0.5d0)
    c000WWj = floor(probeXYZsdist(2)/dx2_j(jd) + 0.5d0)
    c000WWk = floor(probeXYZsdist(3)/dx3_k(kd))

    ! Perform trilinear interpolation on WU, WV and WW
    ! x-component of u3
    c000 = u3(c000WUi,     c000WUj,     c000WUk    )
    c001 = u3(c000WUi,     c000WUj,     c000WUk + 1)
    c011 = u3(c000WUi,     c000WUj + 1, c000WUk + 1)
    c010 = u3(c000WUi,     c000WUj + 1, c000WUk    )
    c110 = u3(c000WUi + 1, c000WUj + 1, c000WUk    )
    c100 = u3(c000WUi + 1, c000WUj,     c000WUk    )
    c101 = u3(c000WUi + 1, c000WUj,     c000WUk + 1)
    c111 = u3(c000WUi + 1, c000WUj + 1, c000WUk + 1)

    wuVal = IBM_probe_interp_trilinear(                     &
                        c000, c001, c011, c010,             &
                        c110, c100, c101, c111,             &
                        dx1(id), dx2_j(jd), dx3_k(kd),      &
                        probeXYZsdist(1), probeXYZsdist(2), &
                        probeXYZsdist(3),                   &
                        x1_i(c000WUi), x2(c000WUj),         &
                        x3(c000WUk)                         )

    ! y-component of u3
    c000 = u3(c000WVi,     c000WVj,     c000WVk    )
    c001 = u3(c000WVi,     c000WVj,     c000WVk + 1)
    c011 = u3(c000WVi,     c000WVj + 1, c000WVk + 1)
    c010 = u3(c000WVi,     c000WVj + 1, c000WVk    )
    c110 = u3(c000WVi + 1, c000WVj + 1, c000WVk    )
    c100 = u3(c000WVi + 1, c000WVj,     c000WVk    )
    c101 = u3(c000WVi + 1, c000WVj,     c000WVk + 1)
    c111 = u3(c000WVi + 1, c000WVj + 1, c000WVk + 1)

    wvVal = IBM_probe_interp_trilinear(                     &
                        c000, c001, c011, c010,             &
                        c110, c100, c101, c111,             &
                        dx1_i(id), dx2(jd), dx3_k(kd),      &
                        probeXYZsdist(1), probeXYZsdist(2), &
                        probeXYZsdist(3),                   &
                        x1(c000WVi), x2_j(c000WVj),         &
                        x3(c000WVk)                         )

    ! z-component of u3
    c000 = u3(c000WWi,     c000WWj,     c000WWk    )
    c001 = u3(c000WWi,     c000WWj,     c000WWk + 1)
    c011 = u3(c000WWi,     c000WWj + 1, c000WWk + 1)
    c010 = u3(c000WWi,     c000WWj + 1, c000WWk    )
    c110 = u3(c000WWi + 1, c000WWj + 1, c000WWk    )
    c100 = u3(c000WWi + 1, c000WWj,     c000WWk    )
    c101 = u3(c000WWi + 1, c000WWj,     c000WWk + 1)
    c111 = u3(c000WWi + 1, c000WWj + 1, c000WWk + 1)

    wwVal = IBM_probe_interp_trilinear(                     &
                        c000, c001, c011, c010,             &
                        c110, c100, c101, c111,             &
                        dx1_i(id), dx2_j(jd), dx3(kd),      &
                        probeXYZsdist(1), probeXYZsdist(2), &
                        probeXYZsdist(3),                   &
                        x1(c000WWi), x2(c000WWj),           &
                        x3_k(c000WWk)                       )
    
    ! Translate cartesian to sperhical coordinates for
    ! probe point velocity u2.
    spherical_coord = IBM_Cartesian_to_Spherical(           &
                        x1(id), x2(jd), x3_k(kd),           &
                        BubbleBlock(bubbleID)%xCenter,      &
                        BubbleBlock(bubbleID)%yCenter,      &
                        BubbleBlock(bubbleID)%zCenter       )

    ! radial/surface normal direction
    wni_ref = (wuVal - BubbleBlock(bubbleID)%u)*spherical_coord(1,1) +  &
              (wvVal - BubbleBlock(bubbleID)%v)*spherical_coord(1,2) +  &
              (wwVal - BubbleBlock(bubbleID)%w)*spherical_coord(1,3)
    ! inclination/surface tangential direction
    wti_ref = (wuVal - BubbleBlock(bubbleID)%u)*spherical_coord(2,1) +  &
              (wvVal - BubbleBlock(bubbleID)%v)*spherical_coord(2,2) +  &
              (wwVal - BubbleBlock(bubbleID)%w)*spherical_coord(2,3)
    ! azimuth/surface tangential direction
    wsi_ref = (wuVal - BubbleBlock(bubbleID)%u)*spherical_coord(3,1) +  &
              (wvVal - BubbleBlock(bubbleID)%v)*spherical_coord(3,2) +  &
              (wwVal - BubbleBlock(bubbleID)%w)*spherical_coord(3,3)

    ! Translate cartesian to spherical coordinates for IB-u1 velocity
    u_star = BubbleBlock(bubbleID)%u + BubbleBlock(bubbleID)%drdt*spherical_coord(1,1)
    v_star = BubbleBlock(bubbleID)%v + BubbleBlock(bubbleID)%drdt*spherical_coord(1,2)
    w_star = BubbleBlock(bubbleID)%w + BubbleBlock(bubbleID)%drdt*spherical_coord(1,3)

    ! radial/surface normal direction
    wni_star = (u_star - BubbleBlock(bubbleID)%u)*spherical_coord(1,1) +  &
               (v_star - BubbleBlock(bubbleID)%v)*spherical_coord(1,2) +  &
               (w_star - BubbleBlock(bubbleID)%w)*spherical_coord(1,3)
    ! inclination/surface tangential direction
    wti_star = (u_star - BubbleBlock(bubbleID)%u)*spherical_coord(2,1) +  &
               (v_star - BubbleBlock(bubbleID)%v)*spherical_coord(2,2) +  &
               (w_star - BubbleBlock(bubbleID)%w)*spherical_coord(2,3)
    ! azimuth/surface tangential direction
    wsi_star = (u_star - BubbleBlock(bubbleID)%u)*spherical_coord(3,1) +  &
               (v_star - BubbleBlock(bubbleID)%v)*spherical_coord(3,2) +  &
               (w_star - BubbleBlock(bubbleID)%w)*spherical_coord(3,3)

    ! Dirichlet Boundary condition for normal/radial direction by means of 
    ! linear interpolation.
    beta = probeXYZsdist(4)/dx3(kd)
    wnIB = (1 + beta)*wni_star - beta*wni_ref
    wtIB = (1 + beta)*wti_star - beta*wti_ref
    wsIB = (1 + beta)*wsi_star - beta*wsi_ref

    ! Translate spherical boundary condition back to cartesian expression for
    ! the IB-point (only take bubble u-velocity)
    store_w(id, jd, kd) = wnIB*spherical_coord(1,1) + &
                          wtIB*spherical_coord(1,2) + &
                          wsIB*spherical_coord(1,3) + &
                          BubbleBlock(bubbleID)%w

end subroutine IBM_velocity_w