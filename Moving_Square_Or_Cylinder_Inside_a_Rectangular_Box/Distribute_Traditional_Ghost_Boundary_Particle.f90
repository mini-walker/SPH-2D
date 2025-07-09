!**************************************************************************************************************
!
!  SUBROUTINE : Distribute_Traditional_Ghost_Boundary_Particle
!
!  PURPOSE    : Distribute traditional ghost boundary particles on the desired boundary direction
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

subroutine  Distribute_Traditional_Ghost_Boundary_Particle(i_time_step)

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the superior subroutien
    integer,intent(in)::i_time_step                                                  ! Current time step 

    ! Variables in local subroutine
    integer::i,j,k,L,m                                                               ! Variables for loop
    
    ! Variables for file opeartion
    integer::File_index                                                              ! File index
    character(len=100)::File_name                                                    ! File name

    integer::Temp_traditional_ghost_particle_number                                  ! Temp traditional ghost particle number 


    ! Variables for non-reflection boundary (Martin Lastiwka, Mihai Basa paper)
    real(kind=8)::rho_reference,pressure_reference,u_reference                       ! Reference value of the inlet
    real(kind=8)::Temp_c
    real(kind=8)::J_1,J_2,J_3
    real(kind=8),dimension(dim)::Temp_Mirror_Particle_Velocity
    real(kind=8)::Temp_Mirror_Particle_Pressure,Temp_Mirror_Particle_Rho

    !==========================================================================================================



    ! Body of subroutine Distribute_Period_Boundary_Particle


    !----------------------------------------------------------------------------------------------------------
    ! Attention: The variable, 'Ture_total_particle_number', is the sum of (1) Particle_ture_number;
    !                                                                      (2) Wave_maker_particle_number;
    !                                                                      (3) Fix_ghost_particle_number;
    !                                                                      (4) InletOutlet_Boundary_Particle_Number;
    !                                                                      (5) Period_Boundary_Particle_Number.
    !----------------------------------------------------------------------------------------------------------

    Ture_total_particle_number=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number+Period_Boundary_Particle_Number


    !==========================================================================================================
    ! Period boundary in x direction

    Traditional_ghost_particle_number=0

    if ( Traditional_Ghost_X_Direction==1 ) then

        ! Left and right boundary mirror opeartion
        do i=1,Ture_total_particle_number

            if ( particle_type(i)==Ill_Particle_Label ) cycle                                        ! No error particles

            !--------------------------------------------------------------------------------------------------
            ! Mirror the particle near the left boundary 
            if ( Traditional_Ghost_X_Left==1 .and. particle_position(i,1) <= ( Traditional_Ghost_Boundary_X_Min+Traditional_Ghost_Boundary_Zone_Width) .and. particle_position(i,1) >= Traditional_Ghost_Boundary_X_Min ) then

                Traditional_ghost_particle_number=Traditional_ghost_particle_number+1                ! Increase Traditional_ghost_particle_number
                j=Ture_total_particle_number+Traditional_ghost_particle_number                       ! Particle index

                particle_position(j,1)=2*Traditional_Ghost_Boundary_X_Min-particle_position(i,1)
                particle_position(j,2)=particle_position(i,2)

            endif
            !--------------------------------------------------------------------------------------------------


            !--------------------------------------------------------------------------------------------------
            ! Mirror the particle near the right boundary 
            if ( Traditional_Ghost_X_Right==1 .and. particle_position(i,1) >= ( Traditional_Ghost_Boundary_X_Max-Traditional_Ghost_Boundary_Zone_Width) .and. particle_position(i,1) <= Traditional_Ghost_Boundary_X_Max ) then

                Traditional_ghost_particle_number=Traditional_ghost_particle_number+1                ! Increase Traditional_ghost_particle_number
                j=Ture_total_particle_number+Traditional_ghost_particle_number                       ! Particle index

                particle_position(j,1)=2*Traditional_Ghost_Boundary_X_Max-particle_position(i,1)
                particle_position(j,2)=particle_position(i,2)

            endif
            !--------------------------------------------------------------------------------------------------




            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Assign the other physical attributes of the mirror particles
            ! Mirror fluid particle - fluid particle
            position_difference(:) = particle_position(j,:)-particle_position(i,:)
            
            ! Interpolation
            Temp_Mirror_Particle_Pressure = particle_press(i) + particle_rho(i)*DOT_PRODUCT(Gravity_grad,position_difference)
            Temp_Mirror_Particle_Rho      = ( Temp_Mirror_Particle_Pressure-Background_Pressure )/square_c_0+water_rho_0     ! Particle density

            Temp_Mirror_Particle_Velocity = particle_velocity(i,:)

            !--------------------------------------------------------------------------------------------------
            ! Permeable and non-reflecting boundary
            ! Reference: Permeable and non-reflecting boundary conditions in SPH (Martin Lastiwka, Mihai Basa)
            ! Define the reference value
            ! For the inlet particles, the velocity are defined by the boundary conditions
            ! For the Pressure and density, they were exterpolate from the interior domain
            u_reference        = 0.0d0                                  ! Horzontial velocity
            rho_reference      = water_rho_0
            pressure_reference = 0.0d0
            
            Temp_c=sqrt(square_c_0)*(Temp_Mirror_Particle_Rho/water_rho_0)**3

            J_1= -Temp_c**2*(Temp_Mirror_Particle_Rho-rho_reference)                  +(Temp_Mirror_Particle_Pressure-pressure_reference)
            J_2=  Temp_Mirror_Particle_Rho*Temp_c*(particle_velocity(i,1)-u_reference)+(Temp_Mirror_Particle_Pressure-pressure_reference) 
            J_3= -Temp_Mirror_Particle_Rho*Temp_c*(particle_velocity(i,1)-u_reference)+(Temp_Mirror_Particle_Pressure-pressure_reference) 

            Temp_Mirror_Particle_Velocity(1)=u_reference+(J_2-J_3)/(2*Temp_Mirror_Particle_Rho*Temp_c)
            Temp_Mirror_Particle_Pressure=pressure_reference+0.5d0*(J_2+J_3)

            !--------------------------------------------------------------------------------------------------


            ! Slip boundary condition on external boundary 
            if ( ExternalBoundaryCondition==1 ) then

                particle_velocity(j,1)=-Temp_Mirror_Particle_Velocity(1)
                particle_velocity(j,2)= Temp_Mirror_Particle_Velocity(2)
            
            ! No-slip boundary condition on external boundary 
            elseif ( ExternalBoundaryCondition==0 ) then

                particle_velocity(j,1)=-Temp_Mirror_Particle_Velocity(1)
                particle_velocity(j,2)=-Temp_Mirror_Particle_Velocity(2)

            else

                write(*,'(A)') " The 'ExternalBoundaryCondition' variable is not defined successful. (Traditional Ghost)" 
                
            endif

            particle_rho(j)=( Temp_Mirror_Particle_Pressure-Background_Pressure )/square_c_0+water_rho_0     ! Particle density
            particle_mass(j)=particle_mass(i)                                                    ! Particle mass
            particle_volume(j)=particle_mass(j)/particle_rho(j)                                  ! Particle volume
            particle_type(j)=Ghost_Water_Particle_Label                                          ! Particle type
            particle_initial_type(j)=Ghost_Water_Particle_Label                                  ! Particle initial type
            particle_smooth_lengh(j)=particle_smooth_lengh(i)                                    ! Particle smooth lengh
            particle_c(j)=particle_c(i)                                                          ! Particle c
            particle_energy(j)=particle_energy(i)                                                ! Particle energy
            particle_smoke_value(j)=particle_smoke_value(i)                                      ! Particle smoke line value
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



        enddo

    endif
    !==========================================================================================================



    !==========================================================================================================
    if ( Traditional_Ghost_Y_Direction==1 ) then

        ! Have circulation in x direction or not
        if ( Traditional_Ghost_X_Direction==1 ) then

            Temp_traditional_ghost_particle_number=Traditional_ghost_particle_number

        else

            Temp_traditional_ghost_particle_number=0
            
        endif

        ! Top and bottom mirror opeartion
        do i=1,Ture_total_particle_number+Temp_traditional_ghost_particle_number

            if(particle_type(i)==Ill_Particle_Label) cycle        !No error particles
            
            !-------------------------------------------------------------------------------------------------
            ! Mirror the particle near the bottom boundary 
            if( Traditional_Ghost_Y_Front==1 .and. particle_position(i,2) <= ( Traditional_Ghost_Boundary_Y_Min+Traditional_Ghost_Boundary_Zone_Width) .and. particle_position(i,2) >= Traditional_Ghost_Boundary_Y_Min ) then

                Traditional_ghost_particle_number=Traditional_ghost_particle_number+1                ! Increase Traditional_ghost_particle_number
                j=Ture_total_particle_number+Traditional_ghost_particle_number                       ! Particle index

                particle_position(j,1)=particle_position(i,1)
                particle_position(j,2)=2*Traditional_Ghost_Boundary_Y_Min-particle_position(i,2)

            endif
            !-------------------------------------------------------------------------------------------------


            !-------------------------------------------------------------------------------------------------
            ! Mirror the particle near the top boundary 
            if( Traditional_Ghost_Y_Back==1 .and. particle_position(i,2) >= ( Traditional_Ghost_Boundary_Y_Max-Traditional_Ghost_Boundary_Zone_Width) .and. particle_position(i,2) <= Traditional_Ghost_Boundary_Y_Max ) then

                Traditional_ghost_particle_number=Traditional_ghost_particle_number+1                ! Increase Traditional_ghost_particle_number
                j=Ture_total_particle_number+Traditional_ghost_particle_number                       ! Particle index

                particle_position(j,1)=particle_position(i,1)
                particle_position(j,2)=2*Traditional_Ghost_Boundary_Y_Max-particle_position(i,2)

            endif
            !---------------------------------------------------------------------------------------------------



            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Assign the other physical attributes of the mirror particles
            ! Mirror fluid particle - fluid particle
            position_difference(:) = particle_position(j,:)-particle_position(i,:)
            
            ! Interpolation
            Temp_Mirror_Particle_Pressure = particle_press(i) + particle_rho(i)*DOT_PRODUCT(Gravity_grad,position_difference)
            Temp_Mirror_Particle_Rho      = ( Temp_Mirror_Particle_Pressure-Background_Pressure )/square_c_0+water_rho_0     ! Particle density

            Temp_Mirror_Particle_Velocity = particle_velocity(i,:)

            !--------------------------------------------------------------------------------------------------
            ! Permeable and non-reflecting boundary
            ! Reference: Permeable and non-reflecting boundary conditions in SPH (Martin Lastiwka, Mihai Basa)
            ! Define the reference value
            ! For the inlet particles, the velocity are defined by the boundary conditions
            ! For the Pressure and density, they were exterpolate from the interior domain
            u_reference        = 0.0d0                                  ! Horzontial velocity
            rho_reference      = water_rho_0
            pressure_reference = 0.0d0
            
            Temp_c=sqrt(square_c_0)*(Temp_Mirror_Particle_Rho/water_rho_0)**3

            J_1= -Temp_c**2*(Temp_Mirror_Particle_Rho-rho_reference)                  +(Temp_Mirror_Particle_Pressure-pressure_reference)
            J_2=  Temp_Mirror_Particle_Rho*Temp_c*(particle_velocity(i,2)-u_reference)+(Temp_Mirror_Particle_Pressure-pressure_reference) 
            J_3= -Temp_Mirror_Particle_Rho*Temp_c*(particle_velocity(i,2)-u_reference)+(Temp_Mirror_Particle_Pressure-pressure_reference) 

            Temp_Mirror_Particle_Velocity(2)=u_reference+(J_2-J_3)/(2*Temp_Mirror_Particle_Rho*Temp_c)
            Temp_Mirror_Particle_Pressure=pressure_reference+0.5d0*(J_2+J_3)

            !--------------------------------------------------------------------------------------------------


            ! Slip boundary condition on external boundary 
            if ( ExternalBoundaryCondition==1 ) then

                particle_velocity(j,1)= Temp_Mirror_Particle_Velocity(1)
                particle_velocity(j,2)=-Temp_Mirror_Particle_Velocity(2)
            
            ! No-slip boundary condition on external boundary 
            elseif ( ExternalBoundaryCondition==0 ) then

                particle_velocity(j,1)=-Temp_Mirror_Particle_Velocity(1)
                particle_velocity(j,2)=-Temp_Mirror_Particle_Velocity(2)

            else

                write(*,'(A)') " The 'ExternalBoundaryCondition' variable is not defined successful. (Traditional Ghost)" 
                
            endif


            particle_rho(j)=( Temp_Mirror_Particle_Pressure-Background_Pressure )/square_c_0+water_rho_0     ! Particle density
            particle_mass(j)=particle_mass(i)                                                    ! Particle mass
            particle_volume(j)=particle_mass(j)/particle_rho(j)                                  ! Particle volume
            particle_type(j)=Ghost_Water_Particle_Label                                          ! Particle type
            particle_initial_type(j)=Ghost_Water_Particle_Label                                  ! Particle initial type
            particle_smooth_lengh(j)=particle_smooth_lengh(i)                                    ! Particle smooth lengh
            particle_c(j)=particle_c(i)                                                          ! Particle c
            particle_energy(j)=particle_energy(i)                                                ! Particle energy
            particle_smoke_value(j)=particle_smoke_value(i)                                      ! Particle smoke line value
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


        enddo

    endif
    !===========================================================================================================
    



    !----------------------------------------------Very important-----------------------------------------------
    ! Renew the 'Ture_total_particle_number' which contains : (1) Particle_ture_number;
    !                                                         (2) Wave_maker_particle_number;
    !                                                         (3) Fix_ghost_particle_number;
    !                                                         (4) InletOutlet_Boundary_Particle_Number;
    !                                                         (5) Period_Boundary_Particle_Number;
    !                                                         (6) Traditional_ghost_particle_number
    Ture_total_particle_number=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number+Period_Boundary_Particle_Number+Traditional_ghost_particle_number

    !-----------------------------------------------------------------------------------------------------------






    !***********************************************************************************************************
    if ( i_time_step==Actual_start_time_step .and. Current_Processor_ID==Main_Processor ) then

        ! Open the Output file
        File_index = General_File_Port
        File_name  = './Initial_Data/Check_Traditional_Ghost_Particle_Distribution_Result.dat'

        call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
    
        if( IOERROR==0 ) then

            ! Input tecplot header
            write(File_index,*) 'TITLE="DISTRIBUTION"'   
            write(File_index,*) 'VARIABLES= "X" "Y" "VX" "VY" "RHO" "P" '
            
            ! Output fluid particle
            write(File_index,*) "ZONE I=",Particle_ture_number," F=POINT"
            do i=1,Particle_ture_number
                 write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
100              format(8F10.4) 
            enddo
            
            ! Output wave maker particle
            if ( Make_wave_or_not==1 ) then
                write(File_index,*) "ZONE I=",Wave_maker_particle_number," F=POINT"
                do k=1,Wave_maker_particle_number
                    i=k+Particle_ture_number
                    write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                enddo
            endif

            ! Output fix ghost boundary particle
            if ( DistributeFixedBoundaryParticleOrNot==1 ) then
                write(File_index,*) "ZONE I=",Fix_ghost_particle_number," F=POINT"
                do k=1,Fix_ghost_particle_number
                    i=k+Particle_ture_number+Wave_maker_particle_number
                    write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                enddo
            endif

            ! Inlet and outlet boundary particles
            if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then

                write(File_index,*) "ZONE I=",InletOutlet_Boundary_Particle_Number," F=POINT"
                do k=1,InletOutlet_Boundary_Particle_Number
                    i=k+Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number
                    write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                enddo
                
            endif

            ! Output period boundary particle
            if ( DistributePeriodBoundaryParticleOrNot==1 ) then
                write(File_index,*) "ZONE I=",Period_Boundary_Particle_Number," F=POINT"
                do k=1,Period_Boundary_Particle_Number
                    i=k+Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number
                    write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                enddo
            endif

            ! Output traditional ghost boundary particle
            if ( DistributeTraditionalGhostBoundaryParticleOrNot==1 ) then
                write(File_index,*) "ZONE I=",Traditional_ghost_particle_number," F=POINT"
                do k=1,Traditional_ghost_particle_number
                    i=k+Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number+Period_Boundary_Particle_Number
                    write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                enddo
            endif 
            
        endif

        close(File_index)

    endif
    !***********************************************************************************************************


end subroutine  Distribute_Traditional_Ghost_Boundary_Particle
