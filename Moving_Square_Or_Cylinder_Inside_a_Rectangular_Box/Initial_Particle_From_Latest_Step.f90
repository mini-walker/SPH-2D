!**************************************************************************************************************
!  Subroutine: Initial_Particle_From_Latest_Step
!
!  PURPOSE: Set the initial information from Latest domain data for SPH calculation
!
!  Programer: Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location: MUN
!
!  Time：2017.3.18
!
!  Copyright: Memorial University
!
!  Version: 1.0
!
!  Input: domain size and packing information
!
!  Output: all particles information
!
!  Note: MPI version: mpich-3.2
!        Fortran version: Inter Fortran (ifort) 18.0.0
!**************************************************************************************************************

subroutine Initial_Particle_From_Latest_Step()

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI

    implicit none

    !==========================================================================================================
    ! Variables in local subroutine
    integer::i,j,k,L,m                                                   ! Loop Variables

    ! Read file varibales
    character(len=20)::char_time                                         ! Time character type
    character(len=100)::Temp_character                                   ! Temp characters for reading
    real(kind=8)::Temp_real                                              ! Temp real for reading
    integer::File_index                                                  ! File index 
    character(len=100)::File_name                                        ! File name ( The file name character lenghth should be 100 )
    
    ! Variables for latest data flie reading
    integer::nan_value_number                                            ! NAN value number
    integer::skip_line                                                   ! Skip line for reading

    integer::char_start_index,char_end_index                             ! The skip lines number

    !==========================================================================================================

    ! Body of Distribute_Latest_Particle
    
    !Call MPI functions
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )


    !**********************************************************************************************************
    
    ! Transfer the restart time format from real to character 
    write(char_time,'(F10.6)') Restart_Time

    Actual_start_time_step = Int( Restart_Time/Fixed_dt )+1                ! Restart time step

    ! Open the Restart Result File
    File_name  = './Restart_Data/'//trim(adjustl(char_time))//'/'//trim(adjustl(char_time))//'.dat'
    File_index = Tecplot_General_File_Port

    call Initialziting_Reading_File( File_index,File_name,IOERROR )        ! Input variables : File_index,File_name,IOERROR 
     
    ! Synchronize all processors calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

    if( IOERROR==0 ) then

        ! Skip the initial two header titles
        skip_line=2
        do k=1,skip_line
            read(File_index,*) temp_character
        end do

        !------------------------------------------------------------------------------------------------------
        ! Check the particle number of each type and check the results
        read(File_index,'(A100)') temp_character

        ! Get the particle number from 'temp_character'
        do i=1,len( trim(adjustl(temp_character)) )

            if ( temp_character(i:i+1)=='I=' ) then
                char_start_index=i+2
            endif

            if ( temp_character(i:i+1)=='F=' ) then
                char_end_index=i-1
            endif

        enddo

        temp_character = trim(adjustl( temp_character(char_start_index:char_end_index) )) 

        read(temp_character,'(I8)') particle_ture_number
        !------------------------------------------------------------------------------------------------------
        nan_value_number=0
        do i=1,particle_ture_number

            read(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_mass(i),particle_energy(i),particle_smoke_value(i),particle_type(i),particle_initial_type(i),Boundary_particle_type(i),wave_maker_particle_layer(i),particle_division_degree(i)          
100         format(8F16.8,5I)
        
            ! Check the reading data have NAN or not
            do j=1,dim
                if(isnan(particle_position(i,j)) .or. isnan(particle_velocity(i,j)) .or. isnan(particle_rho(i)) .or. isnan(particle_mass(i)) .or. isnan(particle_energy(i)) .or. isnan(particle_smoke_value(i)) ) then
                    nan_value_number=nan_value_number+1
                endif
            enddo
        
        end do

        ! Read wave maker particle
        wave_maker_particle_number=0
        if (make_wave_or_not==1) then

            !--------------------------------------------------------------------------------------------------
            ! Check the particle number of each type and check the results
            read(File_index,'(A100)') temp_character

            ! Get the particle number from 'temp_character'
            do i=1,len( trim(adjustl(temp_character)) )

                if ( temp_character(i:i+1)=='I=' ) then
                    char_start_index=i+2
                endif

                if ( temp_character(i:i+1)=='F=' ) then
                    char_end_index=i-1
                endif

            enddo

            temp_character = trim(adjustl( temp_character(char_start_index:char_end_index) )) 

            read(temp_character,'(I8)') wave_maker_particle_number
            !--------------------------------------------------------------------------------------------------

            do k=1,wave_maker_particle_number
                
                i=k+particle_ture_number
                read(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_mass(i),particle_energy(i),particle_smoke_value(i),particle_type(i),particle_initial_type(i),Boundary_particle_type(i),wave_maker_particle_layer(i),particle_division_degree(i)          

                ! Check the reading data have NAN or not
                do j=1,dim
                    if(isnan(particle_position(i,j)) .or. isnan(particle_velocity(i,j)) .or. isnan(particle_rho(i)) .or. isnan(particle_mass(i)) .or. isnan(particle_energy(i)) .or. isnan(particle_smoke_value(i)) ) then
                        nan_value_number=nan_value_number+1
                    endif
                enddo

            enddo

        endif

        ! Output fix ghost boundary particle
        fix_ghost_particle_number=0
        if (DistributeFixedBoundaryParticleOrNot==1) then
            
            !--------------------------------------------------------------------------------------------------
            ! Check the particle number of each type and check the results
            read(File_index,'(A100)') temp_character

            ! Get the particle number from 'temp_character'
            do i=1,len( trim(adjustl(temp_character)) )

                if ( temp_character(i:i+1)=='I=' ) then
                    char_start_index=i+2
                endif

                if ( temp_character(i:i+1)=='F=' ) then
                    char_end_index=i-1
                endif

            enddo

            temp_character = trim(adjustl( temp_character(char_start_index:char_end_index) )) 

            read(temp_character,'(I8)') fix_ghost_particle_number
            !--------------------------------------------------------------------------------------------------

            do k=1,fix_ghost_particle_number

                i=k+particle_ture_number+wave_maker_particle_number
                read(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_mass(i),particle_energy(i),particle_smoke_value(i),particle_type(i),particle_initial_type(i),Boundary_particle_type(i),wave_maker_particle_layer(i),particle_division_degree(i)          

                ! Check the reading data have NAN or not
                do j=1,dim
                    if(isnan(particle_position(i,j)) .or. isnan(particle_velocity(i,j)) .or. isnan(particle_rho(i)) .or. isnan(particle_mass(i)) .or. isnan(particle_energy(i)) .or. isnan(particle_smoke_value(i)) ) then
                        nan_value_number=nan_value_number+1
                    endif
                enddo

            enddo

        endif

        ! Inlet and outlet boundary particles
        InletOutlet_Boundary_Particle_Number=0
        if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then
            
            !--------------------------------------------------------------------------------------------------
            ! Check the particle number of each type and check the results
            read(File_index,'(A100)') temp_character

            ! Get the particle number from 'temp_character'
            do i=1,len( trim(adjustl(temp_character)) )

                if ( temp_character(i:i+1)=='I=' ) then
                    char_start_index=i+2
                endif

                if ( temp_character(i:i+1)=='F=' ) then
                    char_end_index=i-1
                endif

            enddo

            temp_character = trim(adjustl( temp_character(char_start_index:char_end_index) )) 

            read(temp_character,'(I8)') InletOutlet_Boundary_Particle_Number
            !--------------------------------------------------------------------------------------------------

            do k=1,InletOutlet_Boundary_Particle_Number
                
                i=k+particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
                read(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_mass(i),particle_energy(i),particle_smoke_value(i),particle_type(i),particle_initial_type(i),Boundary_particle_type(i),wave_maker_particle_layer(i),particle_division_degree(i)          

                ! Check the reading data have NAN or not
                do j=1,dim
                    if(isnan(particle_position(i,j)) .or. isnan(particle_velocity(i,j)) .or. isnan(particle_rho(i)) .or. isnan(particle_mass(i)) .or. isnan(particle_energy(i)) .or. isnan(particle_smoke_value(i)) ) then
                        nan_value_number=nan_value_number+1
                    endif
                enddo

            enddo
            
        endif
            
        !------------------------------------------------------------------------------------------------------
        !IF the data is available, exit the loop
        if( nan_value_number/=0 .and. Current_Processor_ID==Main_Processor ) then

            write(*,'(1X,A)') "The restart file is not acceptable, please change the restart time ! "

        endif
        !------------------------------------------------------------------------------------------------------

    endif

    ! Synchronize all processors calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

    close(File_index)
    !**********************************************************************************************************




    !**********************************************************************************************************
    ! Particle volume and mass ( The fluid particle number(particle_ture_number) is defined in the Initial fluid particle subroutine in Subroutine_Module )
    !write(*,*) particle_ture_number    
    !intial_volume_0=tank_square/particle_ture_number           
    intial_volume_0=interior_dx*interior_dy                                        ! Intial volume
    
    !Attentation: the initial water mass is used during the calculation
    !I cost one night to find this bug, donot forget!
    water_particle_mass_0=water_rho_0*intial_volume_0

    ! Ture total particle number 
    Ture_total_particle_number=particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number
     

    ! Define other physical attributes of these particles
    do i=1,Ture_total_particle_number
        
        call status_equation(particle_rho(i),&                                     ! Particle Density(in)
                             particle_energy(i),&                                  ! Particle Energy(in)
                             particle_type(i),&                                    ! Particle Type(in)
                             particle_press(i),&                                   ! Particle Pressure(out)
                             particle_c(i)&                                        ! Particle Sound Velocity(out)
                             ) 
        
        !initialization other information for fixed ghost boundary particles                           
        particle_volume(i)=particle_mass(i)/particle_rho(i)                        ! Particle volume
        particle_smooth_lengh(i)=smooth_length                                     ! Particle smooth lengh
                    
    end do
    !**********************************************************************************************************







    !**********************************************************************************************************
    if (Current_Processor_ID==Main_Processor) then
    
        ! Open the Output file
        File_index = General_File_Port
        File_name  = "./Initial_Data/Initialization_Particle_Restart_"//trim(adjustl(char_time))//".dat"
        call Initialziting_Writing_File( File_index,File_name,IOERROR )            ! Input variables : File_index,File_name,IOERROR
  
        if(IOERROR==0) then

            ! Input tecplot header
            write(File_index,*) "TITLE='DISTRIBUTION'"
            write(File_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
            ! Output fluid particle
            write(File_index,*) "ZONE I=",particle_ture_number," F=POINT"
            do i=1,particle_ture_number
                 write(File_index,500) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
500              format(8F10.4) 
            end do
            
            ! Output wave maker particle
            if (make_wave_or_not==1) then
                write(File_index,*) "ZONE I=",wave_maker_particle_number," F=POINT"
                do k=1,wave_maker_particle_number
                    i=k+particle_ture_number
                    write(File_index,500) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                end do
            endif

            ! Output fix ghost boundary particle
            if (DistributeFixedBoundaryParticleOrNot==1) then
                write(File_index,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
                do k=1,fix_ghost_particle_number
                    i=k+particle_ture_number+wave_maker_particle_number
                    write(File_index,500) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                end do
            endif

            ! Inlet and outlet boundary particles
            if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then

                write(File_index,*) "ZONE I=",InletOutlet_Boundary_Particle_Number," F=POINT"
                do k=1,InletOutlet_Boundary_Particle_Number
                    i=k+particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
                    write(File_index,500) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                end do
                
            endif

        endif
        
        close(File_index)

    end if
    !**********************************************************************************************************

    ! Synchronize all processors calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)


end subroutine Initial_Particle_From_Latest_Step