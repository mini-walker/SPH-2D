!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  Subroutine: Distribute_Latest_Particle
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
!****************************************************************************

subroutine Distribute_Latest_Particle()

    use information_module
    use Public_variable_module
    use MPI

    implicit none

    ! Variables in subroutine
    
    !==========================================================================================================
    integer::i,j,k,L,m                                                  !loop Variables

    !Read file varibales
    integer::ioerror=0                                                  !open file return value
    integer::stat                                                       !read file data return value
    integer::status                                                     !allocate memory status
    integer::skip_lines                                                 !The skip lines number
    character(len=20)::temp_character                                   !temp characters for reading
    integer::file_index                                                 !file index 

    !Variables for latest data flie reading
    integer::Finished_run_step                                           !The steps have finished
    integer::start_time_step                                             !Desired time step
    integer::current_run_step                                            !current time step
    integer::total_input_particle_number                                 !total input particle number (fliud+wave maker+fixed boundary, not contain the ghost particle of fixed boundary particles)
    integer::nan_value_number                                            !NAN value number
    integer::initial_save_time_step                                      !latsted results file first_save_time_step
    integer::end_save_time_step                                          !latsted results file last_save_time_step
    integer::saved_run_step                                              !saved run step
    integer::skip_line                                                   !skip line for reading

    !==========================================================================================================

    ! Body of Distribute_Latest_Particle
    
    !Call MPI functions
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )


    !**********************************************************************************************************
    if (Current_Processor_ID==Main_Processor) then

        !Open the Latest_Time_Results
        file_index=110+Current_Processor_ID
        open(unit=file_index,file='Latest_Time_Results.dat',status="old",position="rewind",action="read",iostat=ioerror) 
         
        if(ioerror==0) then
            
            write(*,*) "The open action of the Latest_Time_Results.dat file is well done! "
             
            ! Skip the initial three header titles
            skip_line=3
            do k=1,skip_line
                read(file_index,*) temp_character
            end do

            !----------------------------------------------------------------------------------------------
            !Get the first and last saved step
            !Get the first time step have saved
            read(file_index,*) temp_character,temp_character,temp_character,temp_character,temp_character,temp_character,temp_character,temp_character,temp_character,temp_character,temp_character,initial_save_time_step

            write(*,*) "******************************************************************"
            write(*,'(1X,A)',advance='NO') "The initial time step saved is :"
            write(*,*) initial_save_time_step
            

            !Rewind to the initial, calculate the number of saved time steps
            rewind(file_index)                                               ! rewind to the initial position of the files
            saved_run_step=0
            do
                read(file_index,*,iostat=status) temp_character
                if(status/=0) exit
                if(trim(adjustl(temp_character))=='ZONE') then
                   saved_run_step=saved_run_step+1
                end if
                !write(*,*) temp_character
            end do
            end_save_time_step=(saved_run_step-1)*save_time_step_for_debug+initial_save_time_step


            write(*,*) "******************************************************************"
            write(*,'(1X,A)') "The time step have runned and saved is :"
            write(*,'(1X,A,I,A,I)') 'From',initial_save_time_step,' to',end_save_time_step
            write(*,*) "******************************************************************"

            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            ! Inquire the user, Do they want to input the start time step or not
            ! If the user don't want, the program will start from the last available time information 
            ! ( Control variable is: "Change_last_start_time_or_not" in information model )

            if(Change_last_start_time_or_not=='y') then
                !(For windows)
                !write(*,'(1X,A,I4,A)',advance='NO') "The step (interger times", save_time_step_for_debug,") you want to start is:" 
                !(For Linux)
                write(*,'(1X,A,I4,A)') "The step (interger times", save_time_step_for_debug,") you want to start is:"               
                read(*,*) start_time_step
                
                ! check the user input is right or not
                do 
                    if(start_time_step>end_save_time_step .or. start_time_step<initial_save_time_step) then
                        write(*,*) "The input start time step is larger than the finished step !" 
                        write(*,'(1X,A)',advance='NO') "The step you can input should :"
                        write(*,'(A,I,A,I)') 'from',initial_save_time_step,'to',end_save_time_step
                        write(*,*) "******************************************************************"
                        write(*,'(1X,A)',advance='NO') "Please input start time step again:"
                        read(*,*) start_time_step
                        write(*,*) "******************************************************************"
                    else
                        exit                                                                                ! exit the loop
                    end if
                end do
                
                ! transfer the input real time step to the save time step order!
                start_time_step=(start_time_step-initial_save_time_step)/save_time_step_for_debug+1
                
            else
                write(*,*) "The program will start from the latest available time step !"
                start_time_step=saved_run_step
            end if
            !-------------------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------------------
            ! Check the data in the saved time step is available or not.
            ! Does the data have 'NAN' value or not.
            ! If data is not acceptable, then start form the latest available data.
            actual_start_time_step=start_time_step
            do
                
                !----------------------------------------------------------------------------------------------
                !Run the the start position of actual_start_time_step data
                rewind(file_index)                ! rewind to the initial position of the files
                current_run_step=0
                do  
                    read(file_index,*,iostat=status) temp_character
                    if(status/=0) exit
                    if(trim(adjustl(temp_character))=='ZONE') then
                       current_run_step=current_run_step+1
                    end if
                
                    if (current_run_step==actual_start_time_step) exit

                    !write(*,*) temp_character
                end do
                !----------------------------------------------------------------------------------------------
            
                !----------------------------------------------------------------------------------------------
                ! backspace one line and start the reading
                backspace(file_index)
                read(file_index,*) temp_character,temp_character,total_input_particle_number
                write(*,'(1X,A,I6,A)',advance='NO') "The actual total input particle number at",(current_run_step-1)*save_time_step_for_debug+initial_save_time_step," is :"
                write(*,*) total_input_particle_number

                nan_value_number=0
                do i=1,total_input_particle_number
                    read(file_index,500) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_press(i),particle_rho(i),particle_mass(i),particle_energy(i),particle_c(i),particle_initial_type(i),wave_maker_particle_layer(i)
    500          format(9F20.8,2I) 
                
                    ! check the reading data have NAN or not
                    do k=1,dim
                        if(isnan(particle_position(i,k)) .or. isnan(particle_velocity(i,k)) .or. isnan(particle_press(i)) .or. isnan(particle_rho(i)) .or. isnan(particle_mass(i)) .or. isnan(particle_energy(i)) .or. isnan(particle_c(i))) then
                            nan_value_number=nan_value_number+1
                        end if
                    end do
                
                    !write(3,500) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_press(i),particle_rho(i),particle_mass(i),particle_energy(i),particle_c(i),particle_initial_type(i)
                end do
            
                !----------------------------------------------------------------------------------------------
                !IF the data is available, exit the loop
                if(nan_value_number/=0) then
                    actual_start_time_step=actual_start_time_step-1
                else
                    exit                         ! exit the loop
                end if
                !----------------------------------------------------------------------------------------------

            
            end do
            
            write(*,*) "******************************************************************"
            write(*,'(1X,A)',advance='NO') "The actual start time step is :"
            write(*,*) (current_run_step-1)*save_time_step_for_debug+initial_save_time_step
            write(*,*) "******************************************************************"

         else
             write(*,*) " The Latest_Time_Results file is not exist! "
         end if
         close(file_index)

         ! Actual start time step (transfer the save step to actual calculation data)
         actual_start_time_step=(current_run_step-1)*save_time_step_for_debug+initial_save_time_step

    end if

    !-----------------------------------------------------------------------------------------------------------
    !Transfer the data to other Processors
    call MPI_BCAST( particle_position,n_total*dim,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_velocity,n_total*dim,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_press,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_rho,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_mass,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_energy,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_c,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    call MPI_BCAST( wave_maker_particle_layer,n_total,MPI_INTEGER,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_initial_type,n_total,MPI_INTEGER,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( actual_start_time_step,1,MPI_INTEGER,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( total_input_particle_number,1,MPI_INTEGER,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    
    !-----------------------------------------------------------------------------------------------------------
     

    !-----------------------------------------------------------------------------------------------------------
     
    !Classify the particles according the particle_initial_type
    particle_ture_number=0
    wave_maker_particle_number=0
    fix_ghost_particle_number=0

    do i=1,total_input_particle_number
        
        if(particle_initial_type(i)==2) then
            particle_ture_number=particle_ture_number+1                 !Fluid particle number
        else if(particle_initial_type(i)==130) then
            wave_maker_particle_number=wave_maker_particle_number+1     !Wave maker particle number
        else if(particle_initial_type(i)==-3) then
            fix_ghost_particle_number=fix_ghost_particle_number+1       !Fix ghost particle number
        end if
        
        !initialization other information for fixed ghost boundary particles
        particle_type(i)=particle_initial_type(i)                           
        particle_smooth_lengh(i)=smooth_length                             
    
    end do

    
    intial_volume_0=tank_square/particle_ture_number                        !intial particle volume

    !***********************************************************************************************************
    if (Current_Processor_ID==Main_Processor) then
    
        !Open the Output file
        open(unit=5,file="Latest_particle_information.dat",status="replace",position="rewind",action="write",iostat=ioerror)    
        if(ioerror==0) then
            write(*,*) "The open action of the Latest_particle_information.dat is fine!"
            
            !input tecplot header
            write(5,*) "TITLE='DISTRIBUTION'"
            write(5,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
            !output fluid particle
            write(5,*) "ZONE I=",particle_ture_number," F=POINT"
            do i=1,particle_ture_number
                 write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
100              format(8F20.10) 
            end do
            
            !output wave maker particle
            write(5,*) "ZONE I=",wave_maker_particle_number," F=POINT"
            do k=1,wave_maker_particle_number
                i=k+particle_ture_number
                write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
            end do

            !output fix ghost boundary particle
            write(5,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
            do k=1,fix_ghost_particle_number
                i=k+particle_ture_number+wave_maker_particle_number
                write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
            end do
            
        else
            write(*,*) "The open action of the particle_information.dat is not right!"
        end if
        close(5)

    end if
    !***********************************************************************************************************
    
    !Ture total particle number 
    ture_total_particle_number=particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number

 
end subroutine Distribute_Latest_Particle