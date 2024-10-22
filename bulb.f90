subroutine IBM()
   use var2
   implicit none
   integer :: counter,o,mm
   ! real(8) :: xc_old
   real(8) :: rnode, udummy

   change=1
   cube_number=0
   neighbor(1:im*jm*km/2,1:6)=0

   ! Check distance with electrode and correct if necessary
   ! xc0=xc
   ! xc=radius+1d0*x1_i(1)
   ! dxcdt=(xc-xc0)/dtime
  
   nx_start=1
   udummy=0d0

   ! Check if the entire bubble is inside the domain
   if ( ((xc - radius).lt.x1(1)).or.((xc + radius).gt.x1(im)) ) then
      print*, 'xc: ', xc, ', yc: ', yc, ', zc: ', zc, ', r: ', radius
      error stop "IBM: bubble is outside the domain in x-direction"
   end if
   if ( ((yc - radius).lt.x2(1)).or.((yc + radius).gt.x2(jm)) ) then
      print*, 'xc: ', xc, ', yc: ', yc, ', zc: ', zc, ', r: ', radius
      error stop "IBM: bubble is outside the domain in y-direction"
   end if
   if ( ((zc - radius).lt.x3(1)).or.((zc + radius).gt.x3(km)) ) then
      print*, 'xc: ', xc, ', yc: ', yc, ', zc: ', zc, ', r: ', radius
      error stop "IBM: bubble is outside the domain in z-direction"
   end if

!$omp parallel do private(i,j,k) 
   do k=0,km+1
      do j=0,jm+1
         do i=0,im+1
            typ(i,j,k)=0
            typ1(i,j,k)=0
            typ2(i,j,k)=0
            typ3(i,j,k)=0
            typ_IB(i,j,k)=0
         end do
      end do
   end do

!$omp parallel do private(i,j,k) 
   do k=1,km
      do j=1,jm
         do i=1,im
            if ((x1(i)-xc)**2+(x2(j)-yc)**2+(x3(k)-zc)**2.gt.radius**2) then
               typ(i,j,k)=1
            else
               u1(i-1,j,k)=udummy
               u1(i,j,k)=udummy
               u2(i,j-1,k)=udummy
               u2(i,j,k)=udummy
               u3(i,j,k-1)=udummy
               u3(i,j,k)=udummy
            end if
         end do
      end do
   end do
       
   do k=1,km
      do j=1,jm
         do i=1,im
            counter=0
            x1_plus(change)=0
            x1_minus(change)=0
            x2_plus(change)=0
            x2_minus(change)=0
            x3_plus(change)=0
            x3_minus(change)=0

! Checks if the cell is a bubble cell (id=0)
            if (typ(i,j,k).eq.0) then
               cube_number = cube_number + 1
               ! Checks if neighbouring cells are fluid cells
               if (typ(i+1,j,k).eq.1) then
                  counter=counter+1
                  x1_plus(change)=1
               end if
               if (typ(i-1,j,k).eq.1) then
                  counter=counter+1
                  x1_minus(change)=1
               end if
               if (typ(i,j+1,k).eq.1) then
                  counter=counter+1
                  x2_plus(change)=1
               end if
               if (typ(i,j-1,k).eq.1) then
                  counter=counter+1
                  x2_minus(change)=1
               end if
               if (typ(i,j,k+1).eq.1) then
                  counter=counter+1
                  x3_plus(change)=1
               end if
               if (typ(i,j,k-1).eq.1) then
                  counter=counter+1
                  x3_minus(change)=1
               end if

!              IB-cells are stored in x,y,z _IB and type_IB with id=1 
               if (counter .ge. 1) then
                  x_IB(change)=i
                  y_IB(change)=j
                  z_IB(change)=k

                  typ_IB(i,j,k) = 1
                  
                  change= change + 1
               end if
            end if
         end do
      end do
   end do

!  If fluid cell has a neighboring Fluid-cell (id=1), then mark cell face 
   ! as non-IB.
   do k=1,km
      do j=1,jm
         do i=1,im
            if(typ(i,j,k)*typ(i+1,j,k).eq.1) then
               typ1(i,j,k)=1
            end if
            if(typ(i,j,k)*typ(i,j+1,k).eq.1) then
               typ2(i,j,k)=1
            end if
            if(typ(i,j,k)*typ(i,j,k+1).eq.1) then
               typ3(i,j,k)=1
            end if
         end do
      end do
   end do

   ! Set the following faces to bubble and fluid id
   typ1(0,1:jm,1:km)=0
   typ1(im,1:jm,1:km)=0

   typ2(1:im,0,1:km)=0
   typ2(1:im,jm,1:km)=1

   typ3(1:im,1:jm,0)=0
   typ3(1:im,1:jm,km)=1

   ! Identify probe points for scalar Boundary condition for 
   mm=1
   do o=1,change-1
      i=x_IB(o)
      j=y_IB(o)
      k=z_IB(o)

      phiS(o)=datan2((x2(j)-yc),(x1(i)-xc))

      thetaS(o)=acos((x3(k)-zc)/sqrt((x1(i)-xc)**2+ &
      (x2(j)-yc)**2+(x3(k)-zc)**2))


      T_ni(o)=cos(phiS(o))*sin(thetaS(o))
      T_nj(o)=sin(phiS(o))*sin(thetaS(o))
      T_nk(o)=cos(thetaS(o))

      T_ti(o)=-sin(phiS(o))
      T_tj(o)=cos(phiS(o))
      T_tk(o)=0d0

      T_si(o)=cos(phiS(o))*cos(thetaS(o))
      T_sj(o)=sin(phiS(o))*cos(thetaS(o))
      T_sk(o)=-sin(thetaS(o))

! c         initialise x-direction
! c         Normal distance between interface and bubble centre
      rnode=sqrt((x1(i)-xc)**2+(x2(j)-yc)**2+(x3(k)-zc)**2)


! c        normal distance of interface to bubble surface 
! c        (can be negative)
      rb(o)=radius-rnode

      xb_c(o)=x1(i)+rb(o)*T_ni(o)
      yb_c(o)=x2(j)+rb(o)*T_nj(o)
      zb_c(o)=x3(k)+rb(o)*T_nk(o)

      x_prob_c(o,mm)=xb_c(o)+dx1_i(0)*T_ni(o)
      y_prob_c(o,mm)=yb_c(o)+dx1_i(0)*T_nj(o)
      z_prob_c(o,mm)=zb_c(o)+dx1_i(0)*T_nk(o)


! c         prob number of cells starts at an interface so only 
! c         +0.5d if that direction is on the level of a node
      nx_prob_c(o,mm)=floor(x_prob_c(o,mm)/dx1_i(0)+0.5d0)
      ny_prob_c(o,mm)=floor(y_prob_c(o,mm)/dx2_j(0)+0.5d0)
      nz_prob_c(o,mm)=floor(z_prob_c(o,mm)/dx3_k(0)+0.5d0)

   end do
end subroutine IBM