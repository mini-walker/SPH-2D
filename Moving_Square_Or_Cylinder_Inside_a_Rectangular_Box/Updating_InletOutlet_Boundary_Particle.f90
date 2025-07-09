!**************************************************************************************************************
!
!  SUBROUTINE : Updating_InletOutlet_Boundary_Particle
!
!  PURPOSE    : Distribute in/outlet boundary particles on the desired boundary direction
!
!  Programer  : Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location   : MUN
!
!  Time       : 2017.3.18
!
!  Copyright  : Memorial University
!
!  Version    : 1.0
!
!  Note       : MPI version: mpich-3.2
!               Fortran version: Inter Fortran (ifort) 18.0.0
!**************************************************************************************************************

subroutine  Updating_InletOutlet_Boundary_Particle(i_time_step)

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the mother subroutien
    integer,intent(in)::i_time_step                                     ! Current time step 

    ! Variables in program
    integer::i,j,k,L,m                                                  ! Variables for looping
    
    ! Variables for file opeartion
    integer::File_index                                                 ! File index
    character(len=100)::File_name                                       ! File name

    !==========================================================================================================


    ! Body of subroutine Distribute_InletOutlet_Boundary_Particle
   

    !==========================================================================================================

    !----------------------------------------------------------------------------------------------------------
    ! Attention: The variable, old 'Ture_total_particle_number', is the sum of (1) Particle_ture_number;
    !                                                                          (2) Wave_maker_particle_number;
    !                                                                          (3) Fix_ghost_particle_number;
    !                                                                          (4) InletOutlet_Boundary_Particle_Number;
    !                                                                          (5) New_Fix_ghost_particle_number.
    !----------------------------------------------------------------------------------------------------------

    Ture_total_particle_number=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number+New_Fix_ghost_particle_number
    

    !----------------------------------------------------------------------------------------------------------
    ! Moving the inlet and outlet boundary particle
    New_Inlet_Particle_Number=0

    do i=1,Ture_total_particle_number

        ! For the particle moving into the outlet buffer domain
        if (particle_type(i)==Water_Particle_Label) then

            if ( particle_position(i,1)>Outlet_Boundary_X ) then

                particle_type(i)=Outlet_Particle_Label
                
            endif

        elseif ( particle_type(i)==Inlet_Particle_Label ) then

            ! For the particle moving into the interior domain, Assign it as the fluid particle
            ! And added a new inlet boundary particle
            if (particle_position(i,1)>Inlet_Boundary_X) then
               
                particle_type(i)=Water_Particle_Label                                     ! i is the new fluid particle

                ! New inlet boundary particle
                ! Attentation: The start index of the new inlet particle should keep same with the loop Variables
                ! The new inlet particle is saved after old in/outlet boundary particles
                New_Inlet_Particle_Number=New_Inlet_Particle_Number+1
                L=Ture_total_particle_number+New_Inlet_Particle_Number                    ! L is the new inlet boundary particle

                Particle_Position(L,1)=particle_position(i,1)-InletOutlet_Boundary_Zone_Width
                Particle_Position(L,2)=particle_position(i,2)

                particle_type(L)=Inlet_Particle_Label
                particle_initial_type(L)=Inlet_Particle_Label
                
                ! Define the velocity
                ! ***Very important***
                ! For inlet, the velocity are defined by the inlet definition; 
                ! And the pressure and velocity are exterpolation from the fluid domain;
                particle_velocity(L,:)=Ramping_Inlet_Velocity
                

                particle_mass(L)=water_rho_0*intial_volume_0                                                                               
                particle_rho(L)=water_rho_0                                              
                particle_smooth_lengh(L)=smooth_length                   
                particle_c(L)=sqrt(square_c_0)                               
                particle_press(L)=0.0d0
                particle_volume(L)=intial_volume_0                          
                particle_energy(L)=e_0                       

                free_surface_type(L)=free_surface_type(i)                        
                Boundary_particle_type(L)=Boundary_particle_type(i)
                
                particle_smoke_value(L)=particle_smoke_value(i)                           ! Particle smoke line value

            endif

        elseif ( particle_type(i)==Outlet_Particle_Label ) then

            ! Slow the velocity near the boundary
            if (particle_position(i,dim)<=0.25*interior_dx .or. (boundary_size_y-particle_position(i,dim))<=0.25*interior_dx) then
                particle_velocity(i,dim)=0.0d0
            endif

            ! For the particle moving out the interior domain, Delete the outlet bufer particle
            if ( particle_position(i,1)>=(Outlet_Boundary_X+InletOutlet_Boundary_Zone_Width) ) then
                particle_type(i)=Ill_Particle_Label
            endif

            ! ! For the particle return back to the interior domain, Delete it 
            ! if( particle_velocity(i,1)<=Outlet_Boundary_X ) then
            !     particle_type(i)=Ill_Particle_Label
            ! end if 

            ! For the particle return back to the interior domain, transfer it to the fluid particle
            if( particle_position(i,1)<=Outlet_Boundary_X ) then
                particle_type(i)=Water_Particle_Label
            end if 


            ! Out the y boundary 
            if(particle_position(i,dim)<=0.0d0 .or. particle_position(i,dim)>=boundary_size_y) then
                particle_type(i)=Ill_Particle_Label
            end if

        endif

    enddo
    !----------------------------------------------------------------------------------------------------------

    !==========================================================================================================




end subroutine  Updating_InletOutlet_Boundary_Particle