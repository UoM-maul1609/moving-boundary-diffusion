	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>diffusion solver routines 
    module diffusion
    
    use nrtype
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>variables and types for the moving boundary diffusion model

    implicit none
    real(sp), parameter :: rhow=1000._sp, mw=18.e-3_sp
    
    !>@brief
    !>main model prognostic variables
    type grid
    
        ! variables for grid
        integer(i4b) :: kp, ntim, kp_cur, kp_cur_old
        real(sp) :: dt,rad,t,p,rh,mwsol,rhosol,d_coeff,da_dt, rad_old, rad_min,rad_max
        real(sp), dimension(:), allocatable :: r,u,r05,d,d05, dr,dr05, vol, vol_old, &
                    r_old, r05_old, dr_old, dr05_old
        real(sp), dimension(:,:), allocatable :: c, cold
                                                 
    end type grid



                                            
            

    !>@brief
    !>variables for namelist input
    type namelist_input
        character (len=200) :: inputfile='input'
        character (len=200) :: outputfile='output'
        integer(i4b) :: kp
        real(sp) :: runtime,dt,rad,rad_min,rad_max,t,p,rh,mwsol,rhosol,nu,d_coeff
    end type namelist_input





    
    !>@brief
    !>variables for NetCDF file output
    type io
        ! variables for io
        integer(i4b) :: ncid, varid, x_dimid, y_dimid, z_dimid, &
                        dimids(2), a_dimid, xx_dimid, yy_dimid, &
                        zz_dimid, i_dimid, j_dimid, k_dimid, nq_dimid, nprec_dimid
        integer(i4b) :: icur=1
        logical :: new_file=.true.
    end type io




    ! declare a namelist type
    type(namelist_input) :: nmd
    ! declare a grid type
    type(grid) :: gridd
    ! declare an io type
    type(io) :: iod


    private
    public :: backward_euler, move_boundary, allocate_and_set_diff, diffusion_driver, &
        gridd, nmd,iod, grid
    contains
       
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>solves one step using backward euler method
	!>@param[inout] kp: number of grid points
	!>@param[inout] kpp: number of grid points to solve over
	!>@param[inout] dt: time-step
	!>@param[inout] r: radius array
	!>@param[inout] r05: radius array - half levels
	!>@param[inout] u: velocity of boundary
	!>@param[inout] d: diffusion coefficient
	!>@param[inout] d05: ditto - half levels
	!>@param[inout] dr: grid spacing array
	!>@param[inout] dr05: grid spacing array - half levels
	!>@param[inout] c: concentration array
	!>@param[inout] cold: old concentration array
	!>@param[inout] flux: flux on outer boundary
    subroutine backward_euler(kp,kpp,dt, r,r05,u,d,d05,dr,dr05,c,cold,flux)
		use nrtype
        use nr, only : tridag
        implicit none
		integer(i4b), intent(in) :: kp, kpp
		real(sp), intent(inout) :: flux
		real(sp), intent(in) :: dt

		real(sp), intent(inout), dimension(0:kp+1) :: r, r05, u, d, d05, &
		                                                    dr, dr05
		real(sp), intent(inout), dimension(1:kp+1,1:2) :: c, cold
		

		
		real(sp), dimension(1:kpp) :: cd, b, x, beta, alpha, gamma
		real(sp), dimension(kpp-1) :: ud, ld
		integer(i4b) :: i
		
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! solve tridiagonal matrix:									    	     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ck+1
        !u=-u
        flux=0._sp
        u=0. !1.e-10_sp
        alpha=-(dt*d05(1:kpp)*r05(1:kpp)**2 / &
                (r(1:kpp)**2*dr05(0:kpp-1)*dr(1:kpp)) )
        ! ck
        beta=1._sp+dt*d05(1:kpp)*r05(1:kpp)**2/(r(1:kpp)**2*dr05(0:kpp-1)*dr(1:kpp))+&
            dt*d05(0:kpp-1)*r05(0:kpp-1)**2/(r(1:kpp)**2*dr05(0:kpp-1)*dr(0:kpp-1))
        ! ck-1
        gamma=- dt*d05(0:kpp-1)*r05(0:kpp-1)**2 / &
            (r(1:kpp)**2*dr05(0:kpp-1)*dr(0:kpp-1))
        
        
        do i=1,2
            b=c(1:kpp,i)  ! solution vector
!             if (i .eq. 1) then
!                 b(kpp)=b(kpp) +flux*(dr(kpp-1)+dr(kpp))/d(kpp)*alpha(kpp)
!             endif
            ud=alpha(1:kpp-1)
            ud(1)=ud(1)+gamma(1) ! added to calculate flux across center==0
            cd=beta
            ld=gamma(2:kpp)
            ld(kpp-1)=ld(kpp-1)+alpha(kpp) ! added to calculate outer flux
            call tridag(ld,cd,ud,b,x)
            c(1:kpp,i)=x
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
    end subroutine backward_euler
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>moves the boundary and calculates their velocity
	!>@param[inout] kp: number of grid points
    subroutine move_boundary(kp,kp_cur,dt,radiusold,radius,r,r05,dr,dr05,vol,u,c,flux, &
        rad_min, rad_max, mwsol, rhosol, deltaV)
		use nrtype
		use nr, only : locate
        implicit none
		integer(i4b), intent(inout) :: kp, kp_cur
		real(sp), intent(inout) :: radiusold, radius, flux, mwsol, rhosol
		real(sp), intent(in) :: dt,rad_min,rad_max, deltaV

		real(sp), intent(inout), dimension(0:kp+1) :: r, r05, u,dr, dr05
		real(sp), intent(inout), dimension(1:kp) :: vol
		real(sp), intent(inout), dimension(1:kp+1,1:2) :: c
		
		real(sp), dimension(0:kp+1) :: rold, r05u
		real(sp), dimension(1:kp+1,1:2) :: ctmp
		real(sp) :: da_dt, max_flux, rnew, mf, deltaV2, volwtot, v, volw_first_layer, &
		            volo_outer, deltaVo, volstot
		real(sp), dimension(1:kp+1,1:2) :: moles
		real(sp), dimension(2) :: moles_outer
		integer(i4b) :: i, j, k, kp_new, iter
		

        do i=0,kp+1
            r05u(i)=10._sp**( (log10(rad_max)-log10(rad_min)) &
                /real(kp,sp)*real(i,sp)+ &
                log10(rad_min) )
        enddo

		
		if(deltaV .gt. 0._sp) then
		
		
		
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    ! add water to aerosol particle:                                             !
		    ! assumes it add so the solute of the outer shell is distributed over the    !
		    ! new shells, and mole fraction is the same across all new shells            !
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    rnew=(3._sp*deltaV/(4._sp*pi)+r05(kp_cur)**3)**(1._sp/3._sp)
		    
		    do i=1,2
    		    ! total number of moles
    		    moles(1:kp,i)=c(1:kp,i)*vol
    		    ! moles in outer layer
    		    moles_outer(i)=moles(kp_cur,i)
    		enddo
            ! add water to outer layer
            moles_outer(1)=moles_outer(1)+deltaV*rhow/mw    		    
    		
            ! mole fraction in the outer layer
            mf=moles_outer(1) / sum(moles_outer)
            
            ! concept here is to set the mole fraction in the new layers to mf
            call set_nodes(kp,kp_new,rnew,rad_min,rad_max, r,r05,dr,dr05,vol,c)
            do j=kp_cur,kp_new
                ! water:
                c(j,1)=1._sp/(mw/rhow+(1._sp/mf-1._sp)*mwsol/rhosol)
                ! other component:
                c(j,2)=1._sp/ ( mw/rhow/(1._sp/mf-1._sp)+mwsol/rhosol )
            enddo
!             radius=rnew
            kp_cur=kp_new


            ! set the values of radius / volume and the boundary
            call set_nodes(kp,kp_cur,rnew,rad_min,rad_max, r,r05,dr,dr05,vol,c)
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            radius=rnew
            
            c(kp_cur+1:kp,:)=0._sp
    	else if(deltaV .lt. 0._sp) then
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    ! shrink the particle:                                                       !
		    ! takes off the water from outer layers, then fill outer with solute         !
		    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	    deltaV2=-deltaV !-1e-35_sp ! volume to take off
		    do i=1,2
    		    ! total number of moles
    		    moles(1:kp,i)=c(1:kp,i)*vol
    		enddo
    	    volstot=sum(moles(1:kp_cur,2)) *mwsol/rhosol ! total volume of water
    	    volwtot=sum(moles(1:kp_cur,1)) *mw/rhow ! total volume of water
    	    if (volwtot .le. 0._sp) return
    	    deltaV2=min(volwtot*0.999_sp,deltaV2) ! total volume to remove
    	    
    	    
    	    
    	    
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! take off the water                                                     !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            k=kp_cur
    	    v=0._sp
            do while (v .lt. deltaV2)
                v=v+moles(k,1)*mw/rhow ! volume of water in k integrate inwards
                if(v .ge. deltaV2) exit
                k=k-1 ! this is the index of the new outer shell
            enddo
            ! the new "first" / outer layer has v-deltaV2 water in it
            volw_first_layer = v-deltaV2
!             volw_first_layer= deltaV2*min(volw_first_layer/deltaV2,1._sp)
            ! new radius - with just the water
            rnew=(3._sp*(volw_first_layer) / (4._sp*pi)+r05(k-1)**3)**(1._sp/3._sp) 
            rnew=max(r05(k-1)+1.e-15_sp,rnew)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! assertion:                                                             !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(.not. ((rnew .gt. r05u(k-1)) .and. (rnew .le. r05u(k)))) then
                print *,'assert1: ',rnew, r05u(k-1),r05u(k),r05(k)
                call exit(1)
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! set moles of water:
            moles(k,1)=(volw_first_layer)*rhow/mw
            moles(k+1:kp_cur,1)=0._sp ! no water above last layer
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! move the outer solute mass inwards                                     !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            volo_outer=sum(moles(k:kp_cur,2)) *mwsol/rhosol ! total volume of solute
            iter=0
            do while ((volo_outer .gt. 0._sp) )
                iter=iter+1
                ! this is how much shell volume the solute,in this search, occupies:
                deltaVo=min(4._sp*pi/3._sp*(r05u(k)**3-r05u(k-1)**3)- &
                    moles(k,1)*mw/rhow, volo_outer)
                if(deltaVo .lt. 0._sp) then 
                    print *,'stopping'
                    stop
                endif
                ! corresponding number of moles:
                moles(k,2)=rhosol/mwsol*deltaVo 
                volo_outer=volo_outer-deltaVo
                if(volo_outer .le. 0._sp) exit
                k=k+1
            enddo
!             print *,iter 
            ! find new volume of both components:
            v=moles(k,1)*mw/rhow
            v=v+moles(k,2)*mwsol/rhosol
            ! new radius
            rnew=(3._sp*v/(4._sp*pi)+r05u(k-1)**3)**(1._sp/3._sp)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! assertion:                                                             !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(.not. ((rnew .ge. r05u(k-1)) .and. (rnew .lt. r05u(k)))) then
                print *,'assert2: ',rnew, r05u(k-1),r05u(k)
                call exit(1)
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



            ! set the values of radius / volume and the boundary
            call set_nodes(kp,kp_cur,rnew,rad_min,rad_max, r,r05,dr,dr05,vol,c)

            c(1:kp,1)=moles(1:kp,1)/vol
            c(1:kp,2)=moles(1:kp,2)/vol
            c(kp_cur+1:kp,:)=0._sp
!             print *,volwtot, sum(moles(1:kp_cur,1)) *mw/rhow+deltaV2
!             print *, volstot, sum(moles(1:kp_cur,2)) *mwsol/rhosol
!     	    radius=rnew
		endif
		

    end subroutine move_boundary
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>moves the boundary and calculates their velocity
	!>@param[inout] kp: number of grid points
	!>@param[inout] kp_cur: number of grid points actual
	!>@param[inout] radius: current radius
	!>@param[inout] rad_min: min radius
	!>@param[inout] rad_max: max radius
	!>@param[inout] r: radius array
	!>@param[inout] r05: radius array for boundary
	!>@param[inout] dr: grid spacing
	!>@param[inout] dr05: grid spacing from boundary
	!>@param[inout] vol: volume of grid points
	!>@param[inout] c: concentration of components
    subroutine set_nodes(kp,kp_cur,radius,rad_min, rad_max,r,r05,dr,dr05,vol,c)
		use nrtype
		use nr, only : locate
        implicit none
		integer(i4b), intent(inout) :: kp, kp_cur
		real(sp), intent(inout) :: radius
		real(sp), intent(in) :: rad_min,rad_max

		real(sp), intent(inout), dimension(0:kp+1) :: r, r05, dr, dr05
		real(sp), intent(inout), dimension(1:kp) :: vol
		real(sp), intent(inout), dimension(1:kp+1,1:2) :: c
		
		integer(i4b) :: i,k
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! set arrays                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=0,kp+1
            r05(i)=10._sp**( (log10(rad_max)-log10(rad_min)) &
                /real(kp,sp)*real(i,sp)+ &
                log10(rad_min) )
        enddo
        r(1:kp+1)=(r05(0:kp)+r05(1:kp+1))/2._sp
        r(0)=0._sp
        
        kp_cur=locate(r05(0:kp+1),radius)
        kp_cur=kp_cur
        
        ! set the grid boundary to be the radius
        r05(kp_cur)=radius
        r(kp_cur)=(r05(kp_cur)+r05(kp_cur-1))/2._sp
        ! set the grid point to be equidistant, outside of drop radius
        r(kp_cur+1)=2._sp*radius-r(kp_cur)

        r05(kp_cur+1)=(r(kp_cur+1)+r(kp_cur+2))/2._sp
        
        r(1:kp+1)=(r05(0:kp)+r05(1:kp+1))/2._sp

        
        dr05(0:kp)=r05(1:kp+1)-r05(0:kp)
        dr05(kp+1)=dr05(kp)
        dr(0:kp)=r(1:kp+1)-r(0:kp)
        dr(kp+1)=dr(kp)
        vol=4._sp*pi/3._sp *(r05(1:kp)**3-r05(0:kp-1)**3)

        
    end subroutine set_nodes
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate and set arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays, and initialise them
	!>@param[in] nm_kp: number of grid points
	!>@param[in] nm_dt,nm_runtime: dt, runtime
	!>@param[in] nm_rad: radius of particle
	!>@param[in] nm_rad_min: min radius of particle
	!>@param[in] nm_rad_max: max radius of particle
	!>@param[in] nm_t: temperaturepressurenumber of grid points
	!>@param[in] nm_rh: rh
	!>@param[in] nm_mwsol: molecular weight of solute
	!>@param[in] nm_rhosol: density of solute
	!>@param[in] nm_d_coeff: diffusion coefficient
	!>@param[inout] kp, kp_cur: number of grid points, where to solve
	!>@param[inout] ntim: number of time levels
	!>@param[inout] rad: radius of particle
	!>@param[inout] rad_min: min radius of particle
	!>@param[inout] rad_max: max radius of particle
	!>@param[inout] t: temperaturepressurenumber of grid points
	!>@param[inout] rh: rh
	!>@param[inout] mwsol: molecular weight of solute
	!>@param[inout] nm_nu: van hoff factor
	!>@param[inout] rhosol: density of solute
	!>@param[inout] d_coeff: diffusion coefficient
	!>@param[inout] r, r_old: radius array
	!>@param[inout] r05, r05_old: radius array - half levels
	!>@param[inout] u: velocity of boundary
	!>@param[inout] d: diffusion coefficient
	!>@param[inout] d05: ditto - half levels
	!>@param[inout] dr, dr_old: grid spacing array
	!>@param[inout] dr05, dr05_old: grid spacing array - half levels
	!>@param[inout] vol, vol_old: volume of shell
	!>@param[inout] c: concentration array
	!>@param[inout] cold: old concentration array
	subroutine allocate_and_set_diff(nm_kp,nm_dt,nm_runtime,nm_rad, &
	        nm_rad_min, nm_rad_max, nm_t,nm_p,&
	         nm_rh, nm_mwsol, &
            nm_rhosol, nm_nu, nm_d_coeff, &
            kp,kp_cur,ntim,dt, rad,rad_min,rad_max,t,p,rh, &
            mwsol,rhosol,d_coeff,r,r_old,r05,r05_old, &
            u,d,d05,dr,dr_old,dr05,dr05_old,vol,vol_old, c,cold)
		use nrtype
		use nr, only : locate
		implicit none
		integer(i4b), intent(in) :: nm_kp
		real(sp), intent(in) :: nm_dt,nm_runtime,nm_rad, &
		                        nm_rad_min, nm_rad_max, nm_t,nm_p,nm_rh,&
		                        nm_mwsol,nm_rhosol,nm_nu, nm_d_coeff
		integer(i4b), intent(inout) :: kp,kp_cur, ntim
		real(sp), intent(inout) :: dt,rad,t,p,rh,mwsol,rhosol,d_coeff, &
		            rad_min,rad_max
		real(sp), intent(inout), dimension(:), allocatable :: r, r05, u, d, d05, &
		                                                    dr, dr05, vol, vol_old, &
		                                                    r_old, r05_old, dr_old, &
		                                                    dr05_old
		real(sp), intent(inout), dimension(:,:), allocatable :: c, cold
		
		integer(i4b) :: AllocateStatus, i
		real(sp), dimension(:), allocatable :: nw,na
		


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! set scalars                                                                    !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        kp=nm_kp
        ntim=ceiling(nm_runtime/nm_dt)
        dt=nm_dt
        rad=nm_rad
        t=nm_t
        p=nm_p
        rh=nm_rh
        mwsol=nm_mwsol
        rhosol=nm_rhosol
        d_coeff=nm_d_coeff
        rad_min=nm_rad_min
        rad_max=nm_rad_max
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! allocate arrays                                                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate( r(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( r05(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( u(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( d(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( d05(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( dr(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( dr05(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( vol(1:kp), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        

        allocate( r_old(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( r05_old(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( dr_old(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( dr05_old(0:kp+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( vol_old(1:kp), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        

        allocate( c(1:kp+1,1:2), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( cold(1:kp+1,1:2), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        

        allocate( nw(1:kp), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        allocate( na(1:kp), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call set_nodes(kp,kp_cur,rad,nm_rad_min, nm_rad_max,r,r05,dr,dr05,vol,c)

        d(1:kp+1) = d_coeff
        d05(1:kp+1) = d_coeff
        
        
        
        ! number of moles of water:
        nw=vol / (mw/rhow + (1._sp-rh)/(rh*nm_nu)*mwsol/rhosol)
        na=(1._sp-rh)/(rh*nm_nu)*nw
        
        ! put some more in the outer shell - for testing
!         nw(200:210)=vol(200:210) / (mw/rhow + (1._sp-0.9_sp)/0.9_sp*mwsol/rhosol)
!         na(200:210)=(1._sp-0.9_sp)/0.9_sp*nw(200:210)
        
        
        c(1:kp,1)=nw/vol
        c(1:kp,2)=na/vol
        
        c(kp_cur+1:kp,:)=0._sp
        
        
        ! set the growth rate to zero
        u=0._sp
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! deallocate arrays                                                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(nw)
        deallocate(na)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        
	end subroutine allocate_and_set_diff
	
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calls IO and runs the model
	!>@param[in] outputfile: name of output file
	!>@param[inout] new_file: whether it is a new file or not
	!>@param[inout] kp, kp_cur: number of grid points, up to boundary
	!>@param[inout] ntim: number of time levels
	!>@param[in] runtime: model runtime
	!>@param[inout] rad, rad_min, rad_max: radius of particle, min, max
	!>@param[inout] t: temperaturepressurenumber of grid points
	!>@param[inout] rh: rh
	!>@param[inout] mwsol: molecular weight of solute
	!>@param[inout] rhosol: density of solute
	!>@param[inout] d_coeff: diffusion coefficient
	!>@param[inout] r: radius array
	!>@param[inout] r05: radius array - half levels
	!>@param[inout] u: velocity of boundary
	!>@param[inout] d: diffusion coefficient
	!>@param[inout] d05: ditto - half levels
	!>@param[inout] dr: grid spacing array
	!>@param[inout] dr05: grid spacing array - half levels
	!>@param[inout] vol: volume of shell
	!>@param[inout] c: concentration array
	!>@param[inout] cold: old concentration array
    subroutine diffusion_driver(outputfile, new_file, kp,kp_cur, &
            ntim,dt, runtime, &
            rad,rad_min,rad_max, t,p,rh, &
            mwsol,rhosol,d_coeff,r,r05, &
            u,d,d05,dr,dr05,vol,c,cold)
		use nrtype
        implicit none
        character(len=*), intent(in) :: outputfile
        logical, intent(inout) :: new_file
		integer(i4b), intent(inout) :: kp,kp_cur, ntim
		real(sp), intent(inout) :: dt,rad, &
		                        t,p,rh,mwsol,rhosol,d_coeff, runtime
		real(sp), intent(in) :: rad_min, rad_max
		real(sp), intent(inout), dimension(0:kp+1) :: r, r05, u, d, d05, &
		                                                    dr, dr05
		real(sp), intent(inout), dimension(1:kp) :: vol
		real(sp), intent(inout), dimension(1:kp+1,1:2) :: c, cold
		
		integer(i4b) :: n
		real(sp) :: time, radius, radiusold,flux=0._sp, deltaV
		
		radius=rad
		radiusold=radius
		time=0._sp
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! output                                                                         !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call output(outputfile, new_file, kp, time, radius, r,r05,dr,dr05,c,vol)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! time-loop                                                                      !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do n=1,ntim	
			time=real(n,sp)*dt
			


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! calculate drop growth rate at ambient humidity                             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 			radius=radius+(1.e-9_sp*2._sp*pi/200._sp*cos(2._sp*pi/200._sp*time)*dt)
!             radius=radius-2.e-9_sp
 			deltaV=4._sp*pi/3._sp*(radius**3-radiusold**3)
 			
 			deltaV=(1.e-22*sin(2._sp*pi/200._sp*time))
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! shift radii and calculate the velocity of boundaries                       !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 			call move_boundary(kp,kp_cur,dt,radiusold,radius,r,r05,dr,dr05,vol,u,c,flux, &
 			    rad_min,rad_max, mwsol, rhosol, deltaV)
            radiusold=radius
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Set diffusion coefficient to zero at boundary                              !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			d05(:)=d_coeff
			d05(kp_cur:kp) = 0._sp
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! solve diffusion equation                                                   !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call backward_euler(kp,kp_cur,dt,r,r05,u,d,d05,dr,dr05,c,cold,flux)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		
 
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! output                                                                     !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			call output(outputfile, new_file, kp, time, radius,r,r05,dr,dr05,c,vol)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 		enddo
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








		
	end subroutine diffusion_driver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calls output data from model
	!>@param[in] outputfile: name of output file
	!>@param[inout] new_file: whether it is a new file or not
	!>@param[in] kp: number of grid points
	!>@param[in] time: model time
	!>@param[in] radius: radius 
	!>@param[in] r: radius array
	!>@param[in] r05: radius array - half levels
	!>@param[in] dr: grid spacing array
	!>@param[in] dr05: grid spacing array - half levels
	!>@param[in] c: concentration array
	!>@param[in] vol: volume of each shell
	subroutine output(outputfile, new_file, kp, time, radius, r,r05,dr,dr05,c,vol)
		use nrtype
		use netcdf
        implicit none
        character(len=*), intent(in) :: outputfile
        logical, intent(inout) :: new_file
		integer(i4b), intent(in) :: kp
		real(sp), intent(in) :: time, radius
		real(sp), intent(in), dimension(:) :: r, r05, dr, dr05, vol
		real(sp), intent(in), dimension(:,:) :: c
	
        ! output to netcdf file
        if(new_file) then
            iod%icur=1
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! open / create the netcdf file                                        !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call check( nf90_create(outputfile, NF90_CLOBBER, iod%ncid) )

            ! define dimensions (netcdf hands back a handle)
            call check( nf90_def_dim(iod%ncid, "times", NF90_UNLIMITED, iod%x_dimid) )
            call check( nf90_def_dim(iod%ncid, "kp", kp, iod%y_dimid) )
            call check( nf90_def_dim(iod%ncid, "ncomp", 2, iod%z_dimid) )


            ! close the file, freeing up any internal netCDF resources
            ! associated with the file, and flush any buffers
            call check( nf90_close(iod%ncid) )


            ! now define some variables, units, etc
            call check( nf90_open(outputfile, NF90_WRITE, iod%ncid) )
            ! define mode
            call check( nf90_redef(iod%ncid) )




            ! define variable: time
            call check( nf90_def_var(iod%ncid, "time", NF90_DOUBLE, &
                        (/iod%x_dimid/), iod%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(iod%ncid, "time", iod%a_dimid) )
            ! units
            call check( nf90_put_att(iod%ncid, iod%a_dimid, &
                       "units", "seconds") )

            ! define variable: radius
            call check( nf90_def_var(iod%ncid, "radius", NF90_DOUBLE, &
                        (/iod%x_dimid/), iod%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(iod%ncid, "radius", iod%a_dimid) )
            ! units
            call check( nf90_put_att(iod%ncid, iod%a_dimid, &
                       "units", "m") )

            ! define variable: c
            call check( nf90_def_var(iod%ncid, "c", NF90_DOUBLE, &
                        (/iod%y_dimid, iod%z_dimid, iod%x_dimid/), iod%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(iod%ncid, "c", iod%a_dimid) )
            ! units
            call check( nf90_put_att(iod%ncid, iod%a_dimid, &
                       "units", "moles per m3") )

            ! define variable: r
            call check( nf90_def_var(iod%ncid, "r", NF90_DOUBLE, &
                        (/iod%y_dimid, iod%x_dimid/), iod%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(iod%ncid, "r", iod%a_dimid) )
            ! units
            call check( nf90_put_att(iod%ncid, iod%a_dimid, &
                       "units", "m") )


            ! define variable: vol
            call check( nf90_def_var(iod%ncid, "vol", NF90_DOUBLE, &
                        (/iod%y_dimid, iod%x_dimid/), iod%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(iod%ncid, "vol", iod%a_dimid) )
            ! units
            call check( nf90_put_att(iod%ncid, iod%a_dimid, &
                       "units", "m3") )



            call check( nf90_enddef(iod%ncid) )
            call check( nf90_close(iod%ncid) )

            new_file=.false.

        endif
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! write data to file                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call check( nf90_open(outputfile, NF90_WRITE, iod%ncid) )
        ! write variable: time
        call check( nf90_inq_varid(iod%ncid, "time", iod%varid ) )
        call check( nf90_put_var(iod%ncid, iod%varid, time, &
                    start = (/iod%icur/)))

        ! write variable: radius
        call check( nf90_inq_varid(iod%ncid, "radius", iod%varid ) )
        call check( nf90_put_var(iod%ncid, iod%varid, radius, &
                    start = (/iod%icur/)))

        ! write variable: c
        call check( nf90_inq_varid(iod%ncid, "c", iod%varid ) )
        call check( nf90_put_var(iod%ncid, iod%varid, c(1:kp,1:2), &
                    start = (/1,1,iod%icur/)))

        ! write variable: r
        call check( nf90_inq_varid(iod%ncid, "r", iod%varid ) )
        call check( nf90_put_var(iod%ncid, iod%varid, r(1:kp), &
                    start = (/1,iod%icur/)))

        ! write variable: vol
        call check( nf90_inq_varid(iod%ncid, "vol", iod%varid ) )
        call check( nf90_put_var(iod%ncid, iod%varid, vol(1:kp), &
                    start = (/1,iod%icur/)))

        call check( nf90_close(iod%ncid) )


        iod%icur=iod%icur+1
       
        
        
        
        
	
	end subroutine output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	


	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! HELPER ROUTINE                                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine check(status)
		use netcdf
		use nrtype
		integer(i4b), intent ( in) :: status

		if(status /= nf90_noerr) then
			print *, trim(nf90_strerror(status))
			stop "Stopped"
		end if
	end subroutine check
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end module diffusion
	