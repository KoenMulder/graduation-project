subroutine IBMidentify()
    !!------------------------------------------------------------
    !!  **Main subroutine of IBM cell & face identification**
    !!  
    !!  Requires IBtype, IBtypeU1, IBtypeU2, IBtypeU3, IBbubble,
    !!  IBbubbleU1, IBbubbleU2, IBbubbleU3, and CellIsFlagged to
    !!  be allocated in memory.
    !!
    !!  Only runs when *simBubble = .true.*
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer i,j,k
    ! integer :: n, nIB, nFaceX, nFaceY, nFaceZ
    ! integer :: iCellCenter, jCellCenter, kCellCenter

    ! Identify Immersed Boundary cells
    if ( SimBubble ) then
        ! Reset IBM-arrays
        do i = 0,im+1
            do j= 0,jm+1
                do k= 0,km+1
                    IBtype(i,j,k)   = fluid_label ! IB-type for scalar parameters
                    IBtypeU1(i,j,k) = 0   ! IB-face type for u1-velocity
                    IBtypeU2(i,j,k) = 0   ! IB-face type for u2-velocity
                    IBtypeU3(i,j,k) = 0   ! IB-face type for u3-velocity
                    IBbubble(i,j,k)   = 0 ! Labels scalar IB-type with bubbleID
                    IBbubbleU1(i,j,k) = 0 ! Labels u1 IB-face with bubbleID
                    IBbubbleU2(i,j,k) = 0 ! Labels u2 IB-face with bubbleID
                    IBbubbleU3(i,j,k) = 0 ! Labels u3 IB-face with bubbleID
                end do
            end do
        end do

        IBcellLocation = 0
        

        ! Recursive form cell and face identification
        ! if ( .false. ) then
        !     CellIsFlagged(i,j,k)= 0 ! Required for recursive subroutine "IBM_IBghostCell_FloodFill"
        !     do n = 1, nBubbles
        !         ! Assign bubble ID
        !         BubbleBlock(n)%iBubble = n

        !         ! Idenfity closes cell center to bubble center
        !         ! coordinates:
        !         iCellCenter = nint(BubbleBlock(n)%xCenter/lx1*dble(im))
        !         jCellCenter = nint(BubbleBlock(n)%yCenter/lx2*dble(jm))
        !         kCellCenter = nint(BubbleBlock(n)%zCenter/lx3*dble(km))

        !         ! Recursively identify all ghost and IB-cells
        !         nIB = 0
        !         call IBM_celltype_floodfill( iCellCenter, jCellCenter,  &
        !                                     kCellCenter, n,            &
        !                                     BubbleBlock(n)%xCenter,    &
        !                                     BubbleBlock(n)%yCenter,    &
        !                                     BubbleBlock(n)%zCenter,    &
        !                                     BubbleBlock(n)%radius,     &
        !                                     nIB                        )
        !         BubbleBlock(n)%nIB = nIB
                
        !         ! Recursively identify all faces between IB
        !         ! and fluid cells
        !         nFaceX = 0
        !         nFaceY = 0
        !         nFaceZ = 0
        !         call IBM_facetype_floodfill( iCellCenter, jCellCenter,  &
        !                                     kCellCenter, n,            &
        !                                     BubbleBlock(n)%xCenter,    &
        !                                     BubbleBlock(n)%yCenter,    &
        !                                     BubbleBlock(n)%zCenter,    &
        !                                     BubbleBlock(n)%radius,     &
        !                                     nFaceX, nFaceY, nFaceZ     )
        !         BubbleBlock(n)%nFacesX = nFaceX
        !         BubbleBlock(n)%nFacesY = nFaceY
        !         BubbleBlock(n)%nFacesZ = nFaceZ
        !     end do
        ! end if

        ! Array form cell and face identification
        if ( .true. ) then
            call IBM_cell_face_type()
        end if
    end if    
end subroutine IBMidentify




recursive subroutine IBM_celltype_floodfill(id, jd, kd, bubbleID, xCenter,  &
                                            yCenter, zCenter, radius, nIB   )
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      - Identifies ghost (id=-1), IB (id=1), and fluid (id=0)
    !!      cells and store their label in IBtype-array
    !!      -  Cells are also marked with corresponding bubble 
    !!      label/id
    !!      - IB cell index location is stored in
    !!      IBcellLocation(xyz,nIB,nBubbles)
    !!      - The algorithm works in 6 directions; positive x,y, and
    !!      z and negative x,y, and z. A cell center is marked 
    !!      when it is within the radius of the bubble using 
    !!      the circle equation.
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer, intent(in) :: id, jd, kd, bubbleID
    integer, intent(inout) :: nIB
    real(8), intent(in) :: xCenter, yCenter, zCenter, radius

    ! Apply label
    IBtype(id,jd,kd)    = ghost_label
    IBbubble(id,jd,kd)  = bubbleID

    ! While inside the bubble, apply ghost labels
    ! Search in x, y, and z direction (6 directions)
    ! x-direction: positive and negative
    if ( (((x1(id+1) - xCenter)**2 + (x2(jd) - yCenter)**2 + (x3(kd) - zCenter)**2).lt.radius**2)&
        .and.(IBtype(id+1,jd,kd).eq.fluid_label) ) then
            call IBM_celltype_floodfill(id+1, jd, kd, bubbleID, xCenter, yCenter, zCenter,       &
                                        radius, nIB)
    end if

    if ( (((x1(id-1) - xCenter)**2 + (x2(jd) - yCenter)**2 + (x3(kd) - zCenter)**2).lt.radius**2)&
        .and.(IBtype(id-1,jd,kd).eq.fluid_label) ) then
            call IBM_celltype_floodfill(id-1, jd, kd, bubbleID, xCenter, yCenter, zCenter,      &
                                        radius, nIB)
    end if

    ! y-drection: positive and negative
    if ( (((x1(id) - xCenter)**2 + (x2(jd+1) - yCenter)**2 + (x3(kd) - zCenter)**2).lt.radius**2)&
        .and.(IBtype(id,jd+1,kd).eq.fluid_label) ) then
            call IBM_celltype_floodfill(id, jd+1, kd, bubbleID, xCenter, yCenter, zCenter,      &
                                        radius, nIB)
    end if

    if ( (((x1(id) - xCenter)**2 + (x2(jd-1) - yCenter)**2 + (x3(kd) - zCenter)**2).lt.radius**2)&
        .and.(IBtype(id,jd-1,kd).eq.fluid_label) ) then
            call IBM_celltype_floodfill(id, jd-1, kd, bubbleID, xCenter, yCenter, zCenter,      &
                                        radius, nIB)
    end if

    ! z-direction: positive and negative
    if ( (((x1(id) - xCenter)**2 + (x2(jd) - yCenter)**2 + (x3(kd+1) - zCenter)**2).lt.radius**2)&
        .and.(IBtype(id,jd,kd+1).eq.fluid_label) ) then
            call IBM_celltype_floodfill(id, jd, kd+1, bubbleID, xCenter, yCenter, zCenter,      &
                                        radius, nIB)
    end if

    if ( (((x1(id) - xCenter)**2 + (x2(jd) - yCenter)**2 + (x3(kd-1) - zCenter)**2).lt.radius**2)&
        .and.(IBtype(id,jd,kd-1).eq.fluid_label) ) then
            call IBM_celltype_floodfill(id, jd, kd-1, bubbleID, xCenter, yCenter, zCenter,      &
                                        radius, nIB)
    end if   

    ! At the edge identify IB-cells, which can only be a neighbour when sharing a face.
    ! in x-direction: positive and negative
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id+1,jd,kd).eq.fluid_label).and.  &
        ((x1(id+1) - xCenter)**2 + (x2(jd) - yCenter)**2 + (x3(kd) - zCenter)**2).gt.radius**2) then
            IBtype(id,jd,kd)      = IB_label
            ! Count number of IB-cells
            nIB = nIB + 1
            IBcellLocation(1,nIB,bubbleID) = id
            IBcellLocation(2,nIB,bubbleID) = jd
            IBcellLocation(3,nIB,bubbleID) = kd
    end if

    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id-1,jd,kd).eq.fluid_label).and.   &
        ((x1(id-1) - xCenter)**2 + (x2(jd) - yCenter)**2 + (x3(kd) - zCenter)**2).gt.radius**2) then
            IBtype(id,jd,kd)      = IB_label
            ! Count number of IB-cells
            nIB = nIB + 1
            IBcellLocation(1,nIB,bubbleID) = id
            IBcellLocation(2,nIB,bubbleID) = jd
            IBcellLocation(3,nIB,bubbleID) = kd
    end if

    ! y-direction: positive and negative
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd+1,kd).eq.fluid_label).and.   &
        ((x1(id) - xCenter)**2 + (x2(jd+1) - yCenter)**2 + (x3(kd) - zCenter)**2).gt.radius**2) then
            IBtype(id,jd,kd)      = IB_label
            ! Count number of IB-cells
            nIB = nIB + 1
            IBcellLocation(1,nIB,bubbleID) = id
            IBcellLocation(2,nIB,bubbleID) = jd
            IBcellLocation(3,nIB,bubbleID) = kd
    end if

    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd-1,kd).eq.fluid_label).and.   &
        ((x1(id) - xCenter)**2 + (x2(jd-1) - yCenter)**2 + (x3(kd) - zCenter)**2).gt.radius**2) then
            IBtype(id,jd,kd)      = IB_label
            ! Count number of IB-cells
            nIB = nIB + 1
            IBcellLocation(1,nIB,bubbleID) = id
            IBcellLocation(2,nIB,bubbleID) = jd
            IBcellLocation(3,nIB,bubbleID) = kd
    end if

    ! z-direction: positive and negative
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd,kd+1).eq.fluid_label).and.   &
        ((x1(id) - xCenter)**2 + (x2(jd) - yCenter)**2 + (x3(kd+1) - zCenter)**2).gt.radius**2) then
            IBtype(id,jd,kd+1)      = IB_label
            ! Count number of IB-cells
            nIB = nIB + 1
            IBcellLocation(1,nIB,bubbleID) = id
            IBcellLocation(2,nIB,bubbleID) = jd
            IBcellLocation(3,nIB,bubbleID) = kd
    end if

    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd,kd-1).eq.fluid_label).and.   &
        ((x1(id) - xCenter)**2 + (x2(jd) - yCenter)**2 + (x3(kd-1) - zCenter)**2).gt.radius**2) then
            IBtype(id,jd,kd-1)      = IB_label
            ! Count number of IB-cells
            nIB = nIB + 1
            IBcellLocation(1,nIB,bubbleID) = id
            IBcellLocation(2,nIB,bubbleID) = jd
            IBcellLocation(3,nIB,bubbleID) = kd
    end if
end subroutine IBM_celltype_floodfill




recursive subroutine IBM_facetype_floodfill(id, jd, kd, bubbleID, xCenter,  &
                                            yCenter, zCenter, radius,       &
                                            nFaceX, nFaceY, nFaceZ          )
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      - Identifies IB-faces (id = 1) of the bubble interface
    !!      for u, v, and w-velocity.
    !!
    !!  **Method:**
    !!      - Identification criterion: if a IB-cell has a 
    !!      neighbouring fluid cell, mark all IB-faces for u, v,
    !!      and w-velocity using *IBM_facetype_checkneighbour()*.
    !!      - Stores marked u-IB, v-IB and w-IB faces in IBtypeU1,
    !!      IBtypeU2 and IBtypU3 respectively.
    !!      - The algorithm works in 14 directions; positive x,y,
    !!      and z and negative x,y, and z (6 directions), and
    !!      in diagonal (8 directions); +x,-y,+z, +x,+y,-z,
    !!       +x,-y,-z, +x,-y,+z, and -x,-y,+z, -x,+y,-z, -x,-y,-z,
    !!      -x,-y,+z.
    !!      - Once a cell is analyzed, the cell is flagged to not 
    !!      be repeated.
    !!
    !!  **Requires:**
    !!      - *IBM_facetype_checkneighbour()* to label all faces
    !!      shared between the IB-cell and the fluid cell.
    !!      - Array CellIsFlagged(i,j,k) to be allocated
    !!
    !!------------------------------------------------------------ 
    use commondata
    use paramdata
    implicit none
    integer, intent(in) :: id, jd, kd, bubbleID
    real(8), intent(in) :: xCenter, yCenter, zCenter, radius
    integer, intent(inout) :: nFaceX, nFaceY, nFaceZ

    CellIsFlagged(id,jd,kd) = 1

    ! Search in x,y and z (6 directions)
    ! x-direction: positive and negative
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id+1,jd,kd).eq.ghost_label)  &
    .and.(CellIsFlagged(id+1,jd,kd).eq.0) ) then
        call IBM_facetype_floodfill( id+1, jd, kd, bubbleID, xCenter, yCenter,      &
                                     zCenter, radius, nFaceX, nFaceY, nFaceZ        )
    end if
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id+1,jd,kd).eq.IB_label)     &
    .and.(CellIsFlagged(id+1,jd,kd).eq.0) ) then
        call IBM_facetype_checkneighbour(id+1, jd, kd, bubbleID, nFaceX, nFaceY,    &
                                         nFaceZ)
        CellIsFlagged(id+1,jd,kd) = 1
    end if

    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id-1,jd,kd).eq.ghost_label)  &
    .and.(CellIsFlagged(id-1,jd,kd).eq.0) ) then
        call IBM_facetype_floodfill( id-1, jd, kd, bubbleID, xCenter, yCenter,      &
                                     zCenter, radius, nFaceX, nFaceY, nFaceZ        )
    end if
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id-1,jd,kd).eq.IB_label)     &
    .and.(CellIsFlagged(id-1,jd,kd).eq.0) ) then
        call IBM_facetype_checkneighbour(id-1, jd, kd, bubbleID, nFaceX, nFaceY,    &
                                         nFaceZ)
        CellIsFlagged(id-1,jd,kd) = 1
    end if


    ! y-direction: positive and negative
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd+1,kd).eq.ghost_label)  &
    .and.(CellIsFlagged(id,jd+1,kd).eq.0) ) then
        call IBM_facetype_floodfill( id, jd+1, kd, bubbleID, xCenter, yCenter,      &
                                     zCenter, radius, nFaceX, nFaceY, nFaceZ        )
    end if
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd+1,kd).eq.IB_label)     &
    .and.(CellIsFlagged(id,jd+1,kd).eq.0) ) then
        call IBM_facetype_checkneighbour(id, jd+1, kd, bubbleID, nFaceX, nFaceY,    &
                                         nFaceZ)
        CellIsFlagged(id,jd+1,kd) = 1
    end if

    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd-1,kd).eq.ghost_label)  &
    .and.(CellIsFlagged(id,jd-1,kd).eq.0) ) then
        call IBM_facetype_floodfill( id, jd-1, kd, bubbleID, xCenter, yCenter,      &
                                     zCenter, radius, nFaceX, nFaceY, nFaceZ        )
    end if
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd-1,kd).eq.IB_label)     &
    .and.(CellIsFlagged(id,jd-1,kd).eq.0) ) then
        call IBM_facetype_checkneighbour(id, jd-1, kd, bubbleID, nFaceX, nFaceY,    &
                                         nFaceZ)
        CellIsFlagged(id,jd-1,kd) = 1
    end if


    ! z-direction: positive and negative
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd,kd+1).eq.ghost_label)  &
    .and.(CellIsFlagged(id,jd,kd+1).eq.0) ) then
        call IBM_facetype_floodfill( id, jd, kd+1, bubbleID, xCenter, yCenter,      &
                                     zCenter, radius, nFaceX, nFaceY, nFaceZ        )
    end if
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd,kd+1).eq.IB_label)     &
    .and.(CellIsFlagged(id,jd,kd+1).eq.0) ) then
        call IBM_facetype_checkneighbour(id, jd, kd+1, bubbleID, nFaceX, nFaceY,    &
                                         nFaceZ)
        CellIsFlagged(id,jd,kd+1) = 1
    end if

    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd,kd-1).eq.ghost_label)  &
    .and.(CellIsFlagged(id,jd,kd-1).eq.0) ) then
        call IBM_facetype_floodfill( id, jd, kd-1, bubbleID, xCenter, yCenter,      &
                                     zCenter, radius, nFaceX, nFaceY, nFaceZ        )
    end if
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id,jd,kd-1).eq.IB_label)     &
    .and.(CellIsFlagged(id,jd,kd-1).eq.0) ) then
        call IBM_facetype_checkneighbour(id, jd, kd-1, bubbleID, nFaceX, nFaceY,    &
                                         nFaceZ)
        CellIsFlagged(id,jd,kd-1) = 1
    end if


    ! This part is neccessary to ensure that all cell faces are correctly identified 
    ! and no gaps occur.
    ! Diagonal: +x,-y,+z
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id+1,jd-1,kd+1).eq.ghost_label)  &
    .and.(CellIsFlagged(id+1,jd-1,kd+1).eq.0) ) then
        call IBM_facetype_floodfill(id+1, jd-1, kd+1, bubbleID, xCenter, yCenter,       &
                                    zCenter, radius, nFaceX, nFaceY, nFaceZ             )
    end if
    ! Diagonal: +x,-y,-z
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id+1,jd-1,kd-1).eq.ghost_label)  &
    .and.(CellIsFlagged(id+1,jd-1,kd-1).eq.0) ) then
        call IBM_facetype_floodfill(id+1, jd-1, kd-1, bubbleID, xCenter, yCenter,       &
                                    zCenter, radius, nFaceX, nFaceY, nFaceZ             )
    end if
    ! Diagonal: +x,+y,-z
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id+1,jd+1,kd-1).eq.ghost_label)  &
    .and.(CellIsFlagged(id+1,jd+1,kd-1).eq.0) ) then
        call IBM_facetype_floodfill(id+1, jd+1, kd-1, bubbleID, xCenter, yCenter,       &
                                     zCenter, radius, nFaceX, nFaceY, nFaceZ            )
    end if
    ! Diagonal:  +x,+y,+z
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id+1,jd+1,kd+1).eq.ghost_label)  &
    .and.(CellIsFlagged(id+1,jd+1,kd+1).eq.0) ) then
        call IBM_facetype_floodfill(id+1, jd+1, kd+1, bubbleID, xCenter, yCenter,       &
                                    zCenter, radius, nFaceX, nFaceY, nFaceZ             )
    end if
    ! Diagonal: -x,-y,+z
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id-1,jd-1,kd+1).eq.ghost_label)  &
    .and.(CellIsFlagged(id-1,jd-1,kd+1).eq.0) ) then
        call IBM_facetype_floodfill(id-1, jd-1, kd+1, bubbleID, xCenter, yCenter,       &
                                    zCenter, radius, nFaceX, nFaceY, nFaceZ             )
    end if
    ! Diagonal: -x,-y,-z
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id-1,jd-1,kd-1).eq.ghost_label)  &
    .and.(CellIsFlagged(id-1,jd-1,kd-1).eq.0) ) then
        call IBM_facetype_floodfill(id-1, jd-1, kd-1, bubbleID, xCenter, yCenter,       &
                                    zCenter, radius, nFaceX, nFaceY, nFaceZ             )
    end if
    ! Diagonal: -x,+y,-z
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id-1,jd+1,kd-1).eq.ghost_label)  &
    .and.(CellIsFlagged(id-1,jd+1,kd-1).eq.0) ) then
        call IBM_facetype_floodfill(id-1, jd+1, kd-1, bubbleID, xCenter, yCenter,       &
                                    zCenter, radius, nFaceX, nFaceY, nFaceZ             )
    end if
    ! Diagonal: -x,+y,+z
    if ( (IBtype(id,jd,kd).eq.ghost_label).and.(IBtype(id-1,jd+1,kd+1).eq.ghost_label)  &
    .and.(CellIsFlagged(id-1,jd+1,kd+1).eq.0) ) then
        call IBM_facetype_floodfill(id-1, jd+1, kd+1, bubbleID, xCenter, yCenter,       &
                                    zCenter, radius, nFaceX, nFaceY, nFaceZ             )
    end if
end subroutine IBM_facetype_floodfill




subroutine IBM_facetype_checkneighbour(id, jd, kd, bubbleID, nFaceX, nFaceY, nFaceZ)
    !!------------------------------------------------------------
    !!  **Purpose**:
    !!  - check if neighbouring cell is fluid cell and assign 
    !!  label to shared face if applicable.
    !!  - Stores face location in IBfaceXlocation(xyz,nFaceX,BubbleId)
    !!
    !!  **Method:**
    !!  - Assumes forward staggered grid for cells and faces
    !!
    !!  **Modification:**
    !!  - Added check of IBtypeU1,2,3 if it is already marked. 
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer, intent(in) :: id, jd, kd, bubbleID
    integer, intent(inout) :: nFaceX, nFaceY, nFaceZ

    ! x-direction: positive and negative
    if ( (IBtype(id,jd,kd).eq.IB_label).and.(IBtype(id+1,jd,kd).eq.fluid_label) &
    .and.(IBtypeU1(id,jd,kd).eq.0) ) then
        IBtypeU1(id,jd,kd)   = 1
        IBbubbleU1(id,jd,kd) = bubbleID
        nFaceX = nFaceX + 1
        IBfaceXLocation(1,nFaceX,bubbleID) = id
        IBfaceXLocation(2,nFaceX,bubbleID) = jd
        IBfaceXLocation(3,nFaceX,bubbleID) = kd
    end if
    if ( (IBtype(id,jd,kd).eq.IB_label).and.(IBtype(id-1,jd,kd).eq.fluid_label) &
    .and.(IBtypeU1(id-1,jd,kd).eq.0) ) then
        IBtypeU1(id-1,jd,kd)   = 1
        IBbubbleU1(id-1,jd,kd) = bubbleID
        nFaceX = nFaceX + 1
        IBfaceXLocation(1,nFaceX,bubbleID) = id
        IBfaceXLocation(2,nFaceX,bubbleID) = jd
        IBfaceXLocation(3,nFaceX,bubbleID) = kd
    end if

    ! y-direction: positive and negative
    if ( (IBtype(id,jd,kd).eq.IB_label).and.(IBtype(id,jd+1,kd).eq.fluid_label) &
    .and.(IBtypeU2(id,jd,kd).eq.0) ) then
        IBtypeU2(id,jd,kd)   = 1
        IBbubbleU2(id,jd,kd) = bubbleID
        nFaceY = nFaceY + 1
        IBfaceYLocation(1,nFaceY,bubbleID) = id
        IBfaceYLocation(2,nFaceY,bubbleID) = jd
        IBfaceYLocation(3,nFaceY,bubbleID) = kd
    end if
    if ( (IBtype(id,jd,kd).eq.IB_label).and.(IBtype(id,jd-1,kd).eq.fluid_label) &
    .and.(IBtypeU2(id,jd-1,kd).eq.0) ) then
        IBtypeU2(id,jd-1,kd)   = 1
        IBbubbleU2(id,jd-1,kd) = bubbleID
        nFaceY = nFaceY + 1
        IBfaceYLocation(1,nFaceY,bubbleID) = id
        IBfaceYLocation(2,nFaceY,bubbleID) = jd
        IBfaceYLocation(3,nFaceY,bubbleID) = kd
    end if

    ! z-direction: positive and negative
    if ( (IBtype(id,jd,kd).eq.IB_label).and.(IBtype(id,jd,kd+1).eq.fluid_label) &
    .and.(IBtypeU3(id,jd,kd).eq.0) ) then
        IBtypeU3(id,jd,kd)   = 1
        IBbubbleU3(id,jd,kd) = bubbleID
        nFaceZ = nFaceZ + 1
        IBfaceZLocation(1,nFaceZ,bubbleID) = id
        IBfaceZLocation(2,nFaceZ,bubbleID) = jd
        IBfaceZLocation(3,nFaceZ,bubbleID) = kd
    end if
    if ( (IBtype(id,jd,kd).eq.IB_label).and.(IBtype(id,jd,kd-1).eq.fluid_label) &
    .and.(IBtypeU3(id,jd,kd-1).eq.0) ) then
        IBtypeU3(id,jd,kd-1)   = 1
        IBbubbleU3(id,jd,kd-1) = bubbleID
        nFaceZ = nFaceZ + 1
        IBfaceZLocation(1,nFaceZ,bubbleID) = id
        IBfaceZLocation(2,nFaceZ,bubbleID) = jd
        IBfaceZLocation(3,nFaceZ,bubbleID) = kd
    end if
end subroutine IBM_facetype_checkneighbour




subroutine IBM_cell_face_type()
    !!------------------------------------------------------------
    !!  **Purpose:** 
    !!      - Identifies ghost (id=-1), IB (id=1), and fluid (id=0)
    !!      cells and store their label in IBtype-array
    !!      -  Cells are also marked with corresponding bubble 
    !!      label/id
    !!      - IB cell index location is stored in
    !!      IBcellLocation(xyz,nIB,nBubbles)
    !!
    !!  **Method:**
    !!      - Uses array indexing instead of recursion.
    !!
    !!------------------------------------------------------------
    use commondata
    use paramdata
    implicit none
    integer :: i,j,k,n, nIB, nFaceX, nFaceY, nFaceZ
    real(8) :: xCenter, yCenter, zCenter, radius

    ! Mark all cells (including IB-cells) inside bubble as ghost cell
    do n = 1, nBubbles
        xCenter = BubbleBlock(n)%xCenter
        yCenter = BubbleBlock(n)%yCenter
        zCenter = BubbleBlock(n)%zCenter
        radius  = BubbleBlock(n)%radius
!$OMP PARALLEL DO PRIVATE(i,j,k) 
        do i = 0,im+1
            do j = 0,jm+1
                do k = 0,km+1
                    if ( ((x1(i) - xCenter)**2 + (x2(j) - yCenter)**2 + &
                          (x3(k) - zCenter)**2).lt.(radius**2) ) then
                        IBtype(i,j,k)   = ghost_label
                        IBbubble(i,j,k) = n
                    end if
                end do
            end do
        end do
    end do

    ! Mark IB-cells inside bubble as ghost cell
    do n = 1, nBubbles
        nIB = 0
!$OMP PARALLEL DO PRIVATE(i,j,k) 
        do i = 1,im
            do j = 1,jm
                do k = 1,km
                    ! Check postive and negative x-direction
                    if ( (IBtype(i,j,k).eq.ghost_label).and.(IBtype(i+1,j,k).eq.fluid_label)   &
                    .and.(IBbubble(i,j,k).eq.n) ) then

                        IBtype(i,j,k) = IB_label
                        nIB = nIB + 1
                        IBcellLocation(1,nIB,n) = i
                        IBcellLocation(2,nIB,n) = j
                        IBcellLocation(3,nIB,n) = k
                    end if

                    if ( (IBtype(i,j,k).eq.ghost_label).and.(IBtype(i-1,j,k).eq.fluid_label)   &
                    .and.(IBbubble(i,j,k).eq.n) ) then

                        IBtype(i,j,k) = IB_label
                        nIB = nIB + 1
                        IBcellLocation(1,nIB,n) = i
                        IBcellLocation(2,nIB,n) = j
                        IBcellLocation(3,nIB,n) = k
                    end if

                    ! Check positive and negative y-direction
                    if ( (IBtype(i,j,k).eq.ghost_label).and.(IBtype(i,j+1,k).eq.fluid_label)   &
                    .and.(IBbubble(i,j,k).eq.n) ) then

                        IBtype(i,j,k) = IB_label
                        nIB = nIB + 1
                        IBcellLocation(1,nIB,n) = i
                        IBcellLocation(2,nIB,n) = j
                        IBcellLocation(3,nIB,n) = k
                    end if

                    if ( (IBtype(i,j,k).eq.ghost_label).and.(IBtype(i,j-1,k).eq.fluid_label)   &
                    .and.(IBbubble(i,j,k).eq.n) ) then

                        IBtype(i,j,k) = IB_label
                        nIB = nIB + 1
                        IBcellLocation(1,nIB,n) = i
                        IBcellLocation(2,nIB,n) = j
                        IBcellLocation(3,nIB,n) = k
                    end if

                    ! Check positive and negative z-direciton
                    if ( (IBtype(i,j,k).eq.ghost_label).and.(IBtype(i,j,k+1).eq.fluid_label)   &
                    .and.(IBbubble(i,j,k).eq.n) ) then

                        IBtype(i,j,k) = IB_label
                        nIB = nIB + 1
                        IBcellLocation(1,nIB,n) = i
                        IBcellLocation(2,nIB,n) = j
                        IBcellLocation(3,nIB,n) = k
                    end if

                    if ( (IBtype(i,j,k).eq.ghost_label).and.(IBtype(i,j,k-1).eq.fluid_label)   &
                    .and.(IBbubble(i,j,k).eq.n) ) then

                        IBtype(i,j,k) = IB_label
                        nIB = nIB + 1
                        IBcellLocation(1,nIB,n) = i
                        IBcellLocation(2,nIB,n) = j
                        IBcellLocation(3,nIB,n) = k
                    end if
                end do
            end do
        end do
    end do

    ! Mark all IB-faces
    do n = 1,nBubbles
        nFaceX = 0
        nFaceY = 0
        nFaceZ = 0

!$OMP PARALLEL DO PRIVATE(i,j,k) 
        do i = 1,im
            do j = 1,jm
                do k = 1,km
                    if ( IBbubble(i,j,k).eq.n ) then
                        call IBM_facetype_checkneighbour(i, j,  &
                                     k, n, nFaceX, nFaceY, nFaceZ)  
                    end if
                end do
            end do
        end do

        BubbleBlock(n)%nFacesX = nFaceX
        BubbleBlock(n)%nFacesY = nFaceY
        BubbleBlock(n)%nFacesZ = nFaceZ
    end do

end subroutine IBM_cell_face_type