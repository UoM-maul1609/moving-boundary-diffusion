	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Moving Boundary Diffusion (MBD): 
	!>Solve Fickian diffusion with a moving boundary
    !>
    !>
	!> <br><br>
	!> compile using the Makefile (note requires netcdf) and then run using: <br>
	!> ./main.exe namelist.in
	!> <br><br>
	!> (namelist used for initialisation).
	!> <br><br>



	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main programme reads in information, allocates arrays, then calls the model driver

    program main
        use nrtype
        use diffusion, only : allocate_and_set_diff, diffusion_driver, gridd, &
                            iod, nmd
        
        implicit none
        character (len=200) :: nmlfile = ' '
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelist for run variables                                           !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /run_vars/ nmd
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists													   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=run_vars)
        close(8)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Allocate and initialise arrays for diffusion model         		   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call allocate_and_set_diff(nmd%kp,nmd%dt,nmd%runtime, &
            nmd%rad, nmd%rad_min, nmd%rad_max, nmd%t,nmd%p, nmd%rh, nmd%mwsol, &
            nmd%rhosol, nmd%nu, nmd%d_coeff, &
            gridd%kp,gridd%kp_cur, &
            gridd%ntim,gridd%dt, gridd%rad,gridd%rad_min, gridd%rad_max, &
            gridd%t,gridd%p,gridd%rh, &
            gridd%mwsol,gridd%rhosol,gridd%d_coeff,gridd%r,gridd%r_old,gridd%r05, &
            gridd%r05_old, &
            gridd%u,gridd%d,gridd%d05,gridd%dr,gridd%dr_old,&
            gridd%dr05,gridd%dr05_old,gridd%vol,gridd%vol_old, &
            gridd%c,gridd%cold)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Driver code: time-loop, advance solution, output	   				   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call diffusion_driver(nmd%outputfile, iod%new_file, &
            gridd%kp,gridd%kp_cur,gridd%ntim,gridd%dt, nmd%runtime, &
            gridd%rad, nmd%rad_min, nmd%rad_max, &
            gridd%t,gridd%p,gridd%rh, &
            gridd%mwsol,gridd%rhosol,gridd%d_coeff,gridd%r,gridd%r05, &
            gridd%u,gridd%d,gridd%d05,gridd%dr,gridd%dr05,gridd%vol, &
            gridd%c,gridd%cold)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






    end program main



