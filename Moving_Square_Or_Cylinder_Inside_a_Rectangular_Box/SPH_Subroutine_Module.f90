!**************************************************************************************************************
!  SUBROUTINE : SPH_Subroutine_Module
!
!  PURPOSE    : Subroutine module for SPH calculation 
!               (contains: the subroutine such as: initialization subroutines for SPH iteration and so on)
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


module SPH_Subroutine_Module
       
    USE Information_module
    USE Public_variable_module
    USE Function_module
    USE MPI

    contains

    !==========================================================================================================
    ! Subroutine Public_Variable_Initialization
    ! Initialze the public variables 
    subroutine Public_Variable_Initialization()
    
        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine

        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Public_Variable_Initialization
        
        particle_position=0.0d0
        particle_mass=0.0d0
        particle_rho=0.0d0
        particle_smooth_lengh=0.0d0
        particle_energy=0.0d0
        particle_c=0.0d0
        particle_press=0.0d0

        particle_type=0
        free_surface_type=0
        particle_initial_type=0
        Boundary_particle_type=0
        particle_division_degree=0


        !Initialize the arraies for iteration
        !------------------------------------------------------------------------------------------------------
        !Temp Variables for SPH Iteration
        before_rho=0.0                               !density in last step
        before_energy=0.0                            !energy in last step
        before_velocity=0.0                          !velocity in last step
        before_position=0.0                          !position in last step

        temp_particle_velocity=0.0                   !velocity of Prediction step in Eular iteration
        temp_drhodt=0.0                              !density rate of Prediction step in Eular iteration
        temp_dedt=0.0                                !energy rate of Prediction step in Eular iteration
        temp_particle_acceleration=0.0               !acceleration of Prediction step in Eular iteration

        k_1_temp_particle_velocity=0.0               !velocity of Prediction step in Runge–Kutta iteration
        k_1_temp_drhodt=0.0                          !density rate of Prediction step in Runge–Kutta iteration
        k_1_temp_dedt=0.0                            !energy rate of Prediction step in Runge–Kutta iteration
        k_1_temp_particle_acceleration=0.0           !acceleration of Prediction step in Runge–Kutta iteration

        k_2_temp_particle_velocity=0.0               !velocity of Prediction step in Runge–Kutta iteration
        k_2_temp_drhodt=0.0                          !density rate of Prediction step in Runge–Kutta iteration
        k_2_temp_dedt=0.0                            !energy rate of Prediction step in Runge–Kutta iteration
        k_2_temp_particle_acceleration=0.0           !acceleration of Prediction step in Runge–Kutta iteration

        k_3_temp_particle_velocity=0.0               !velocity of Prediction step in Runge–Kutta iteration
        k_3_temp_drhodt=0.0                          !density rate of Prediction step in Runge–Kutta iteration
        k_3_temp_dedt=0.0                            !energy rate of Prediction step in Runge–Kutta iteration
        k_3_temp_particle_acceleration=0.0           !acceleration of Prediction step in Runge–Kutta iteration
        !------------------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------------------
        !Temp Variables for SPH-ALE Iteration
        before_mass=0.0                              !mass in last step
        before_volume=0.0                            !volume in last step
        before_momentum=0.0                          !momentum in last step

        temp_dmassdt=0.0                             !mass rate of change of Prediction step in Eular iteration
        temp_dvolumedt=0.0                           !volume rate of change of Prediction step in Eular iteration
        temp_dmomentumdt=0.0                         !momentum rate of change of Prediction step in Eular iteration

        k_1_temp_dmassdt=0.0                         !mass rate of Prediction step in Runge–Kutta iteration
        k_1_temp_dvolumedt=0.0                       !volume rate of Prediction step in Runge–Kutta iteration
        k_1_temp_dmomentumdt=0.0                     !momentum rate of Prediction step in Runge–Kutta iteration

        k_2_temp_dmassdt=0.0                         !mass rate of Prediction step in Runge–Kutta iteration
        k_2_temp_dvolumedt=0.0                       !volume rate of Prediction step in Runge–Kutta iteration
        k_2_temp_dmomentumdt=0.0                     !momentum rate of Prediction step in Runge–Kutta iteration

        k_3_temp_dmassdt=0.0                         !mass rate of Prediction step in Runge–Kutta iteration
        k_3_temp_dvolumedt=0.0                       !volume rate of Prediction step in Runge–Kutta iteration
        k_3_temp_dmomentumdt=0.0                     !momentum rate of Prediction step in Runge–Kutta iteration
        !------------------------------------------------------------------------------------------------------

    end subroutine Public_Variable_Initialization
    !==========================================================================================================




    !==========================================================================================================
    ! Subroutine Directory_and_File_Initialization
    ! Initialize the directories and files for data saving
    subroutine Directory_and_File_Initialization()
        
        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                                  ! Loop Variables

        character(len=20)::Char_File_Index                                  ! Character file index
        integer::File_index                                                 ! File index
        character(len=100)::File_name                                       ! File name ( The file name character lenghth should be 100 )
        
        character(len=100)::current_folder_name                             ! Current folder name
        character(len=100)::RemoveCommand,BuildCommand                      ! Remove and Build Command
        Logical::dirExists                                                  ! Variables for folder opeartion

        !Variables for file opeartion
        character(len=4)::char_Current_probe_order                          ! Character of the probe index
        character(len=4)::char_Current_Processor_ID                         ! Character of the processor ID
        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Directory_and_File_Initialization
        
        !******************************************************************************************************
        ! Directory opeartion
        !------------------------------------------------------------------------------------------------------
        ! Initialize the Directory Name
        Directory_Name(1) = 'Subdomain'
        Directory_Name(2) = 'Sampling'
        Directory_Name(3) = 'Sampling/Press'
        Directory_Name(4) = 'Sampling/WaveHeight'
        Directory_Name(5) = 'Paraview'
        Directory_Name(6) = 'Binary_Tecplot'
        Directory_Name(7) = 'Binary_Paraview'
        Directory_Name(8) = 'Binary_Tecplot_Grid'
        Directory_Name(9) = 'Binary_Paraview_Grid'
        Directory_Name(10)= 'Initial_Data'
        Directory_Name(11)= 'Restart_Data'

        !------------------------------------------------------------------------------------------------------
        ! Creat new directory with 'main processor' 
        Make_Directory: if ( Current_Processor_ID==Main_Processor ) then

            do i=1,Directory_Number

                current_folder_name=Directory_Name(i)

                ! Windows system
                if ( Linux_or_Windows==0 ) then                                        ! Windows command

                    RemoveCommand='rd /s /q '//trim(adjustl(current_folder_name))
                    BuildCommand='md '//trim(adjustl(current_folder_name))

                ! Linux system
                elseif ( Linux_or_Windows==1 ) then                                    ! Linux command

                    RemoveCommand='rm -r '//trim(adjustl(current_folder_name))
                    BuildCommand ='mkdir -p '//trim(adjustl(current_folder_name))

                else
                    write(*,'(A)') " The value of 'Linux_or_Windows' in 'Information_Module' is not right! ( 0 or 1)" 
                endif

                ! Inquire the directory exist or not
                inquire(directory=trim(adjustl(current_folder_name)), exist=dirExists )
                
                ! If the SPH calculation start from the first step, we need remove all the existed directories and the files in them,
                ! we just keep the old folders and results when start from the last time step.
                if     ( StartFrom == 0 ) then                                         ! Satrt from initial step

                    if(dirExists) then

                        call system (trim(adjustl(RemoveCommand)))                     ! Remove the Subdomain folder first
                        call system (trim(adjustl(BuildCommand)))                      ! Make the new folder

                    else

                        call system (trim(adjustl(BuildCommand)))                      ! Make the new folder

                    endif                    

                elseif ( StartFrom == 1 ) then                                         ! Satrt from Latest step

                    if(dirExists) then

                    else

                        call system (trim(adjustl(BuildCommand)))                      ! Make the new folder

                    endif                          

                else
                    write(*,'(A)') " The 'startFrom' variable is not right. ('0---InitialTime' or '1---LatestTime')"
                end if

            enddo

        endif Make_Directory
        !------------------------------------------------------------------------------------------------------
        
        call sleep( 5 )                                                                ! Sleep 5 seconds
        call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)                                    ! Synchronize all processors calculation 


        !------------------------------------------------------------------------------------------------------
        ! Check the directory creation opeartions successful or not
        j=0
        do i=1,Directory_Number

            current_folder_name=Directory_Name(i)

            inquire(directory=trim(adjustl(current_folder_name)), exist=dirExists )
            if(dirExists) then

                j=j+1

            else

                write(*,14) trim(adjustl(current_folder_name)),Current_Processor_ID
    14          Format(" The ' ",A," ' folder for Processor",I5," is Wrong, the program stop running!")
                
                ! Stop the MPI running
                call MPI_FINALIZE(ierror_MPI)

            endif

        enddo

        ! Check the existed directory number equal to the desired directory number or not
        if (j==Directory_Number) then

            write(*,15) Current_Processor_ID
    15      Format(" All folders are ready for the SPH simulation in Processor",I5," !")

        else

            write(*,'(A)') ' Folders opeartions are not successful, the program stop running!**'
            
            ! Finish the MPI running
            call MPI_FINALIZE(ierror_MPI)
          
        endif
        !------------------------------------------------------------------------------------------------------

        call sleep( 2 )                                                                ! Sleep 2 seconds
        call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)                                    ! Synchronize all processors calculation
        
        !******************************************************************************************************





        !******************************************************************************************************
        ! File opeartion
        !------------------------------------------------------------------------------------------------------
        

        !**************************************************************************************************
        !This is for linux, on Windows the foler opeartion is not same
        if ( startFrom == 0 ) then

            !----------------------------------------------------------------------------------------------
            !Wave sampling files
            do i=1,Wave_Probe_number

                Wave_Probe_File_index(i)=Wave_Probe_File_Initial_index+i     !Wave_Probe_File_index=Wave_Probe_File_Initial_index+wave_probe_order
              
                !transfer the Current_Processor_ID from integer to character
                write(char_Current_probe_order,'(I4)') i

                File_name="./Sampling/WaveHeight/wave_gauge_"//trim(adjustl(char_Current_probe_order))//"_result.dat"
                call Initialziting_Writing_File( Wave_Probe_File_index(i),File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
          
            end do
            !----------------------------------------------------------------------------------------------


            !----------------------------------------------------------------------------------------------
            !Pressure sampling files
            do i=1,Sampling_Point_number

                Normalized_Pressure_Sampling_File_index(i) = Normalized_Pressure_Sampling_File_Initial_index+i        !Normalized_Pressure_Sampling_File_index=Normalized_Pressure_Sampling_File_Initial_index+sampling_probe_order
                Average_Pressure_Sampling_File_index(i)    = Average_Pressure_Sampling_File_Initial_index+i           !Average_Pressure_Sampling_File_index=Average_Pressure_Sampling_File_Initial_index+sampling_probe_order
                Kalman_Pressure_Sampling_File_index(i)     = Kalman_Pressure_Sampling_File_Initial_index+i            !Kalman_Pressure_Sampling_File_index=Kalman_Pressure_Sampling_File_Initial_index+sampling_probe_order

                !transfer the Current_Processor_ID from integer to character
                write(char_Current_probe_order,'(I4)') i

                File_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Normalized_Pressure.dat"
                call Initialziting_Writing_File( Normalized_Pressure_Sampling_File_index(i),File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
          
                File_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Average_Pressure.dat"
                call Initialziting_Writing_File( Average_Pressure_Sampling_File_index(i),File_name,IOERROR )       ! Input variables : File_index,File_name,IOERROR
          
                File_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Kalman_Pressure.dat"
                call Initialziting_Writing_File( Kalman_Pressure_Sampling_File_index(i),File_name,IOERROR )        ! Input variables : File_index,File_name,IOERROR
          
            end do
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            ! ASCII animation results
            File_name  = "Animation_Results.dat"
            File_index = Animation_Result_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then
                !Tecplot header
                write(File_index,*) 'TITLE="DISTRIBUTION"'
                write(File_index,*) 'VARIABLES= "X" "Y" "VX" "VY" "P" "U" "Vor" '
            endif
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            ! ASCII animation results for grid
            File_name  = "Animation_Grid_Results.dat"
            File_index = Animation_Grid_Result_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then
                !Tecplot header
                write(File_index,*) 'TITLE="DISTRIBUTION"'
                write(File_index,*) 'VARIABLES= "X" "Y" "VX" "VY" "P" "U" "Vor" "Fai" '
            endif
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            ! Free surface results
            File_name  = "Free_Surface.dat"
            File_index = Free_Surface_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then
                !Tecplot header
                write(File_index,*) ' TITLE="DISTRIBUTION" '
                write(File_index,*) ' VARIABLES= "X" "Y" "Nor_X" "Nor_Y" '
            endif
            !----------------------------------------------------------------------------------------------
            
            !----------------------------------------------------------------------------------------------
            ! Body force results
            File_name  = "Body_Force.dat"
            File_index = Body_Force_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then
                !Tecplot header
                write(File_index,*) "# Time Pressure(1--dim) Viscous(1--dim) Total(1--dim)"
            endif
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            ! Free surface results
            File_name  = "Smoke_Line.dat"
            File_index = Smoke_Line_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then
                !Tecplot header
                write(File_index,*) ' TITLE="DISTRIBUTION" '
                write(File_index,*) ' VARIABLES= "X" "Y" "Smoke" '
            endif
            !----------------------------------------------------------------------------------------------


            !----------------------------------------------------------------------------------------------
            !Wave front results
            File_name  = "Wave_front_results.dat"
            File_index = Wave_Front_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            !Kinetic energy results
            File_name  = "Kinetic_energy.dat"
            File_index = Kinetic_Energy_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            !Maximum velocity magnitude
            File_name  = "Maximum_velocity_magnitude.dat"
            File_index = Maximum_velocity_magnitude_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            !----------------------------------------------------------------------------------------------
            
              
            !Normalized results for taylor-green vortex
            !----------------------------------------------------------------------------------------------
            File_name  = "Norm_Horizontal_velocity.dat"
            File_index = Norm_Horizontal_velocity_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            File_name  = "Norm_Vertical_velocity.dat"
            File_index = Norm_Vertical_velocity_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR

            File_name  = "Norm_Velocity_magnitude.dat"
            File_index = Norm_Velocity_magnitude_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR

            File_name  = "Norm_Pressure.dat"
            File_index = Norm_Pressure_File_Port
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
            !----------------------------------------------------------------------------------------------

      
        elseif( startFrom == 1 ) then
          

            !----------------------------------------------------------------------------------------------
            !Wave sampling files
            do i=1,Wave_Probe_number

                Wave_Probe_File_index(i)=Wave_Probe_File_Initial_index+i     !Wave_Probe_File_index=Wave_Probe_File_Initial_index+wave_probe_order
              
                !transfer the Current_Processor_ID from integer to character
                write(char_Current_probe_order,'(I4)') i

                File_name="./Sampling/WaveHeight/wave_gauge_"//trim(adjustl(char_Current_probe_order))//"_result.dat"
                call Initialziting_Existed_Writing_File( Wave_Probe_File_index(i),File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
          
            end do
            !----------------------------------------------------------------------------------------------


            !----------------------------------------------------------------------------------------------
            !Pressure sampling files
            do i=1,Sampling_Point_number

                Normalized_Pressure_Sampling_File_index(i) = Normalized_Pressure_Sampling_File_Initial_index+i        !Normalized_Pressure_Sampling_File_index=Normalized_Pressure_Sampling_File_Initial_index+sampling_probe_order
                Average_Pressure_Sampling_File_index(i)    = Average_Pressure_Sampling_File_Initial_index+i           !Average_Pressure_Sampling_File_index=Average_Pressure_Sampling_File_Initial_index+sampling_probe_order
                Kalman_Pressure_Sampling_File_index(i)     = Kalman_Pressure_Sampling_File_Initial_index+i            !Kalman_Pressure_Sampling_File_index=Kalman_Pressure_Sampling_File_Initial_index+sampling_probe_order

                !transfer the Current_Processor_ID from integer to character
                write(char_Current_probe_order,'(I4)') i

                File_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Normalized_Pressure.dat"
                call Initialziting_Existed_Writing_File( Normalized_Pressure_Sampling_File_index(i),File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
          
                File_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Average_Pressure.dat"
                call Initialziting_Existed_Writing_File( Average_Pressure_Sampling_File_index(i),File_name,IOERROR )       ! Input variables : File_index,File_name,IOERROR
          
                File_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Kalman_Pressure.dat"
                call Initialziting_Existed_Writing_File( Kalman_Pressure_Sampling_File_index(i),File_name,IOERROR )        ! Input variables : File_index,File_name,IOERROR
          
            end do
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            ! ASCII animation results
            File_name  = "Animation_Results.dat"
            File_index = Animation_Result_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then
                !Tecplot header
                write(File_index,*) 'TITLE="DISTRIBUTION"'
                write(File_index,*) 'VARIABLES= "X" "Y" "VX" "VY" "P" "U" "Vor" '
            endif
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            ! ASCII animation results for grid
            File_name  = "Animation_Grid_Results.dat"
            File_index = Animation_Grid_Result_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then
                !Tecplot header
                write(File_index,*) 'TITLE="DISTRIBUTION"'
                write(File_index,*) 'VARIABLES= "X" "Y" "VX" "VY" "P" "U" "Vor" "Fai" '
            endif
            !----------------------------------------------------------------------------------------------


            !----------------------------------------------------------------------------------------------
            ! Free surface results
            File_name  = "Free_Surface.dat"
            File_index = Free_Surface_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then
                !Tecplot header
                write(File_index,*) ' TITLE="DISTRIBUTION" '
                write(File_index,*) ' VARIABLES= "X" "Y" "Nor_X" "Nor_Y" '
            endif
            !----------------------------------------------------------------------------------------------
            
            !----------------------------------------------------------------------------------------------
            ! Body force results
            File_name  = "Body_Force.dat"
            File_index = Body_Force_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then
                !Tecplot header
                write(File_index,*) "# Time Pressure(1--dim) Viscous(1--dim) Total(1--dim)"
            endif
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            ! Free surface results
            File_name  = "Smoke_Line.dat"
            File_index = Smoke_Line_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then
                !Tecplot header
                write(File_index,*) ' TITLE="DISTRIBUTION" '
                write(File_index,*) ' VARIABLES= "X" "Y" "Smoke" '
            endif
            !----------------------------------------------------------------------------------------------


            !----------------------------------------------------------------------------------------------
            !Wave front results
            File_name  = "Wave_front_results.dat"
            File_index = Wave_Front_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            !Kinetic energy results
            File_name  = "Kinetic_energy.dat"
            File_index = Kinetic_Energy_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            !----------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------
            !Maximum velocity magnitude
            File_name  = "Maximum_velocity_magnitude.dat"
            File_index = Maximum_velocity_magnitude_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            !----------------------------------------------------------------------------------------------
            
              
            !Normalized results for taylor-green vortex
            !----------------------------------------------------------------------------------------------

            File_name  = "Norm_Horizontal_velocity.dat"
            File_index = Norm_Horizontal_velocity_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            File_name  = "Norm_Vertical_velocity.dat"
            File_index = Norm_Vertical_velocity_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR

            File_name  = "Norm_Velocity_magnitude.dat"
            File_index = Norm_Velocity_magnitude_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR

            File_name  = "Norm_Pressure.dat"
            File_index = Norm_Pressure_File_Port
            call Initialziting_Existed_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
            !----------------------------------------------------------------------------------------------

        else

            write(*,'(A)') " The startFrom variable is not right. ('0---InitialTime' or '1---latestTime')"
              
        endif
        !------------------------------------------------------------------------------------------------------



        !------------------------------------------------------------------------------------------------------
        !open the result checking file for each processor

        !transfer the Current_Processor_ID from integer to character
        write(char_Current_Processor_ID,'(I4)') Current_Processor_ID
        !write(*,*) char_Current_Processor_ID

        File_name="./Subdomain/Results_in_subdomain_"//trim(adjustl(char_Current_Processor_ID))//".dat"
        File_index=Initial_File_Port_For_Each_Processor+Current_Processor_ID
        Current_Processor_File_Index=File_index
        call Initialziting_Writing_File( File_index,File_name,IOERROR )              ! Input variables : File_index,File_name,IOERROR
  
        if(IOERROR==0) then
            !tecplot header
            write(File_index,*) 'TITLE="DISTRIBUTION"' 
            write(File_index,*) 'VARIABLES= "X" "Y" "VX" "VY" "RHO" "P"' 
        endif
        !------------------------------------------------------------------------------------------------------



        !******************************************************************************************************

    end subroutine Directory_and_File_Initialization
    !==========================================================================================================









































    !==========================================================================================================
    !                                                                _
    !                                        _._ _..._ .-',     _.._(`))
    !                                       '-. `     '  /-._.-'    ',/
    !                                          )         \            '.
    !                                         / _    _    |             \
    !                                        |  a    a    /              |
    !                                        \   .-.                     ;
    !                                         '-('' ).-'       ,'       ;
    !                                            '-;           |      .'
    !                                               \           \    /
    !                                               | 7  .__  _.-\   \
    !                                               | |  |  ``/  /`  /
    !                                              /,_|  |   /,_/   /
    !                                                 /,_/      '`-'
    !==========================================================================================================




    !==========================================================================================================
    ! Subroutine Initial_Fluid_Particle_From_Tecplot
    ! Initialze the fluid particle from Tecplot file (Particle packing result) 
    subroutine Initial_Fluid_Particle_From_Tecplot()
    
        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                                             ! loop Variables
        
        !Read file varibales
        integer::skip_lines                                                            ! The skip lines number
        character(len=20)::temp_character                                              ! Temp characters for reading
        integer::File_index                                                            ! File index
        character(len=100)::File_name                                                  ! File name ( The file name character lenghth should be 100 )
            
        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Initial_Fluid_Particle_From_Tecplot
        

        !******************************************************************************************************
        ! Open the input Tecplot data files
        File_index=General_File_Port
        File_name =Trim(adjustl(Tecplot_Initial_Particle_File_name))  
        call Initialziting_Reading_File( File_index,File_name,ioerror )

        ! Read the particle data from the file
        if( ioerror==0 ) then

            ! Skip the initial three head titles
            skip_lines=3
            do i=1,skip_lines
               read(File_index,*) temp_character
            end do
            
            k=0
            do
                k=k+1                                                     ! Calculator plus 1
                read(File_index,*,iostat=status) (particle_position(k,j),j=1,dim)
                if(status/=0) exit
                
                !particle_press(k)=rho_0*g*(y-particle_position(k,2))                                 
                !particle_rho(k)=rho_0*((1+7*g*(y-particle_position(k,2))/square_c_0)**(0.142857142)) 

                particle_rho(k)=water_rho_0*(-6.0*g*(particle_position(k,2)-domain_size_y)/square_c_0+1.0)**(1.0/6.0)
                particle_press(k)=square_c_0*water_rho_0*((particle_rho(k)/water_rho_0)**7.0-1.0)/7.0
  
            end do

            particle_ture_number=k-1                                                      ! Actual partcile number

        endif

        ! Synchronize all processors calculation           
        call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

        !close the file index
        close(File_index)
        !******************************************************************************************************

    end subroutine Initial_Fluid_Particle_From_Tecplot
    !==========================================================================================================






    !==========================================================================================================
    ! Subroutine Initial_Fluid_Particle_From_Gambit
    ! Initialze the fluid particle from Gambit file ( Gambit grid file ) 
    subroutine Initial_Fluid_Particle_From_Gambit()
    
        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                                             ! Loop Variables
        
        !Read file varibales
        integer::skip_lines                                                            ! The skip lines number
        character(len=20)::temp_character                                              ! Temp characters for reading
        integer::File_index                                                            ! File index
        character(len=100)::File_name                                                  ! File name ( The file name character lenghth should be 100 )
    
        !Variables for gambit initialization 
        real(kind=8),allocatable,dimension(:,:)::node                                  ! Node coordinate(node index, x/y/z)
        integer,allocatable,dimension(:,:)::element                                    ! Nodes in gambit element(node index, 1/2/3)
        integer::node_number,element_number                                            ! Node and element number
        integer::temp_int                                                              ! Temp integer value
        real(kind=8)::position                                                         ! Position value
        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Initial_Fluid_Particle_From_Gambit
        

        !******************************************************************************************************
        !open the input Gambit data files
        File_index=General_File_Port
        File_name =Trim( adjustl(Gambit_Initial_Particle_File_name) )  
        call Initialziting_Reading_File( File_index,File_name,ioerror )
        
        if ( ioerror==0 ) then

            !Skip the initial five head titles
            skip_lines=5
            do i=1,skip_lines
               read(File_index,*) temp_character
            end do

            !Read node and element number
            read(File_index,*) node_number,element_number

            !Allocate memory for node, element, and center points
            allocate(node(node_number,dim))                                                              !node coordinate
            allocate(element(element_number,gambit_element_node_number))                                 !element point index

            !Skip two head titles
            skip_lines=2
            do i=1,skip_lines
               read(File_index,*) temp_character
            end do
            
            !Read node coordinates
            do i=1,node_number
                read(File_index,*) temp_int,(node(i,j),j=1,dim)
                !write(3,*) i,(node(i,j),j=1,dim)
            end do
            
            !Skip two head titles
            skip_lines=2
            do i=1,skip_lines
               read(File_index,*) temp_character
            end do

            !Read the node index in each element
            do i=1,element_number
                read(File_index,*) temp_int,temp_int,temp_int,(element(i,j),j=1,gambit_element_node_number)
                !write(3,*) i,(element(i,j),j=1,gambit_element_node_number)
            end do

            !Calculate the center points for initialization

            do i=1,element_number                             !The actual partcile number is equal to element number
                do j=1,dim
                    position=0.0d0
                    do k=1,gambit_element_node_number
                        position=position+node(element(i,k),j)
                    end do
                    particle_position(i,j)=position/3.0
                end do    
            end do

            !Assign the density and press
            do k=1,element_number

                !particle_press(k)=water_rho_0*g*(domain_size_y-particle_position(k,2))           
                !particle_rho(k)=water_rho_0*((1+7*g*(domain_size_y-particle_position(k,2))/square_c_0)**(0.142857142)) 

                particle_rho(k)=water_rho_0*(-6.0*g*(particle_position(k,2)-domain_size_y)/square_c_0+1.0)**(1.0/6.0)
                particle_press(k)=square_c_0*water_rho_0*((particle_rho(k)/water_rho_0)**7.0-1.0)/7.0

            end do

            !Calculate the fluid particle number and the initial particle number and mass
            particle_ture_number=k                                                        !actual partcile number
            
        endif
        
        ! Synchronize all processors calculation           
        call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

        !close the file index
        close(File_index)
        !******************************************************************************************************

    end subroutine Initial_Fluid_Particle_From_Gambit
    !==========================================================================================================








    !==========================================================================================================
    ! Subroutine Initial_Fluid_Particle_From_Regular_Distribution
    ! Initialze the fluid particle from Regular Distribution
    subroutine Initial_Fluid_Particle_From_Regular_Distribution()
    
        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                                             ! Variables for looping
        
        real(kind=8)::temp_press                                                       ! Partcile temp press
        real(kind=8)::C_k                                                              ! Coefficients
        real(kind=8)::temp_C_k                                                         ! Temp coefficients
        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Initial_Fluid_Particle_From_Regular_Distribution
        
        !******************************************************************************************************
        k=0

        do i=1,interior_particle_number_y                ! 1-m row
            do j=1,interior_particle_number_x            ! 1-m column
                
                k=k+1
                particle_position(k,1)=(j-1)*interior_dx+interior_dx/2.0             
                particle_position(k,2)=(i-1)*interior_dy+interior_dy/2.0             
                
                ! Particle velocity
                particle_velocity(k,:)=0.0D0

                ! Particle pressure
                particle_press(k)=0.0d0+Background_Pressure

                ! Particle density
                particle_rho(k)=(particle_press(k)-Background_Pressure)/square_c_0+water_rho_0

                ! Delete the particle in the backward step
                if (particle_position(k,1)>=(Block_Center_X-0.5*Block_Length) .and. particle_position(k,1)<=(Block_Center_X+0.5*Block_Length) ) then

                    if (particle_position(k,2)>=(Block_Center_Y-0.5*Block_Length) .and. particle_position(k,2)<=(Block_Center_Y+0.5*Block_Length)) then
                        k=k-1   
                    endif

                endif

                ! Particle division degree
                particle_division_degree(k)=0                                         ! Mother particle
                  
            end do
        end do
        
        particle_ture_number=k                                                        ! Actual partcile number
        !------------------------------------------------------------------------------------------------------

        !******************************************************************************************************

    end subroutine Initial_Fluid_Particle_From_Regular_Distribution
    !==========================================================================================================







    !==========================================================================================================
    subroutine Point_Project_Into_Gobal_Background_Grid(Projected_Particle_Position,Grid_Position_Index,Grid_Name)

        implicit none
        
        !------------------------------------------------------------------------------------------------------
        ! Variables for superior-subroutine
        real(kind=8),dimension(dim),intent(in)::Projected_Particle_Position            ! Position of the projected particle
        
        integer,dimension(dim),intent(out)::Grid_Position_Index                        ! The position index of grid which point in 
        integer,intent(out)::Grid_Name                                                 ! The grid name which point in 

        !------------------------------------------------------------------------------------------------------

        ! Body of subroutine Point_Project_Into_Gobal_Background_Grid


        !------------------------------------------------------------------------------------------------------
        
        ! Calculate the coordinates of the particles in the list
        ! Attention: The projected particle should minus the background gird origin and then move forward 1 column 

        if (dim==2) then

            Grid_Position_Index(1) = Floor( (Projected_Particle_Position(1)-chain_origin_x)/chain_dx+1 )
            Grid_Position_Index(2) = Floor( (Projected_Particle_Position(2)-chain_origin_y)/chain_dy+1 )

            ! Calculate the grid number of the particle in the background grid
            Grid_Name = (Grid_Position_Index(2)-1)*chain_x_number+Grid_Position_Index(1)       ! Row by row

        elseif (dim==3) then

            Grid_Position_Index(1) = Floor( (Projected_Particle_Position(1)-chain_origin_x)/chain_dx+1 )
            Grid_Position_Index(2) = Floor( (Projected_Particle_Position(2)-chain_origin_y)/chain_dy+1 )
            Grid_Position_Index(3) = Floor( (Projected_Particle_Position(3)-chain_origin_z)/chain_dz+1 )

            ! Calculate the grid number of the particle in the background grid
            Grid_Name = (Grid_Position_Index(3)-1)*chain_x_number*chain_y_number+(Grid_Position_Index(2)-1)*chain_x_number+Grid_Position_Index(1)       ! Row by row
            
        endif
        !------------------------------------------------------------------------------------------------------

    end subroutine Point_Project_Into_Gobal_Background_Grid
    !==========================================================================================================   





    !==========================================================================================================
    subroutine Point_Project_Into_Subdomain_Background_Grid(Projected_Particle_Position,Grid_Position_Index,Grid_Name)

        implicit none
        
        !------------------------------------------------------------------------------------------------------
        ! Variables for superior-subroutine
        real(kind=8),dimension(dim),intent(in)::Projected_Particle_Position            ! Position of the projected particle

        integer,dimension(dim),intent(out)::Grid_Position_Index                        ! The position index of grid which point in 
        integer,intent(out)::Grid_Name                                                 ! The grid name which point in 

        !------------------------------------------------------------------------------------------------------

        ! Body of subroutine Point_Project_Into_Subdomain_Background_Grid


        !------------------------------------------------------------------------------------------------------

        ! Calculate the coordinates of the particles in the list
        ! Attention: The projected particle should minus the background gird origin and then move forward 2 column 

        if (dim==2) then

            Grid_Position_Index(1) = Floor( (Projected_Particle_Position(1)-subdomain_chain_origin_x)/chain_dx+1 )+1
            Grid_Position_Index(2) = Floor( (Projected_Particle_Position(2)-subdomain_chain_origin_y)/chain_dy+1 )+1

            ! Calculate the grid number of the particle in the background grid
            Grid_Name = (Grid_Position_Index(2)-1)*subdomain_chain_x_number+Grid_Position_Index(1)     ! Row by row

        elseif (dim==3) then

            Grid_Position_Index(1) = Floor( (Projected_Particle_Position(1)-subdomain_chain_origin_x)/chain_dx+1 )+1
            Grid_Position_Index(2) = Floor( (Projected_Particle_Position(2)-subdomain_chain_origin_y)/chain_dy+1 )+1
            Grid_Position_Index(3) = Floor( (Projected_Particle_Position(3)-subdomain_chain_origin_z)/chain_dz+1 )+1

            ! Calculate the grid number of the particle in the background grid
            Grid_Name = (Grid_Position_Index(3)-1)*subdomain_chain_x_number*subdomain_chain_y_number+(Grid_Position_Index(2)-1)*subdomain_chain_x_number+Grid_Position_Index(1)       ! Row by row
            
        endif
        !------------------------------------------------------------------------------------------------------


    end subroutine Point_Project_Into_Subdomain_Background_Grid
    !========================================================================================================== 






    !==========================================================================================================
    ! Initialized the variables for pressure Sampling
    ! Define the sampling points and wave probe position
    ! For 2D simulation, we just need (x,y) position for pressure Sampling and (x) for wave probe 
    ! For 3D simulation, we just need (x,y,z) position for pressure Sampling and (x,y) for wave probe 
    subroutine Sampling_and_Wave_Probe_Initialization()

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                                  ! Loop Variables

        character(len=20)::Char_File_Index                                  ! Character file index
        integer::File_index                                                 ! File index
        character(len=100)::File_name                                       ! File name ( The file name character lenghth should be 100 )
        
        real(kind=8)::Gap_Between_Reflect_boundary                          ! Gap Between Samping Point and Reflect boundary
        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Sampling_and_Wave_Probe_Initialization

        !******************************************************************************************************
        if (dim==2) then

            !--------------------------------------------------------------------------------------------------
            if ( Sampling_Point_number>0 ) then
                
                ! Pressure Sampling position ( For solitary wave simulation )
                ! Sampling_Point_Position(1,2)=0.05
                ! Sampling_Point_Position(1,1)=incline_x+Sampling_Point_Position(1,2)/tan(theta_rad)

                ! Sampling_Point_Position(2,2)=0.15
                ! Sampling_Point_Position(2,1)=incline_x+Sampling_Point_Position(2,2)/tan(theta_rad)
                
                ! Sampling_Point_Position(3,2)=0.25
                ! Sampling_Point_Position(3,1)=incline_x+Sampling_Point_Position(3,2)/tan(theta_rad)
                
                ! Sampling_Point_Position(4,2)=0.35
                ! Sampling_Point_Position(4,1)=incline_x+Sampling_Point_Position(4,2)/tan(theta_rad)

                ! Gap between the reflect boundary
                Gap_Between_Reflect_boundary=0.5d0

                ! ( 0.5, 2.5)
                Sampling_Point_Position(1,2)=domain_size_y/2.0
                Sampling_Point_Position(1,1)=0.0d0+Gap_Between_Reflect_boundary
                
                ! ( 9.5, 2.5)
                Sampling_Point_Position(2,2)=domain_size_y/2.0
                Sampling_Point_Position(2,1)=domain_size_x-Gap_Between_Reflect_boundary

                ! ( 1.0, 2.5)
                Sampling_Point_Position(3,2)=domain_size_y/2.0
                Sampling_Point_Position(3,1)=Block_Center_X-0.50d0
                
                ! ( 4.0, 2.5)
                Sampling_Point_Position(4,2)=domain_size_y/2.0
                Sampling_Point_Position(4,1)=Block_Center_X+2.50d0
                
                ! ( 1.0, 0.5)
                Sampling_Point_Position(5,2)=0.0d0+Gap_Between_Reflect_boundary
                Sampling_Point_Position(5,1)=Block_Center_X-0.50d0
                
                ! ( 1.0, 4.5)
                Sampling_Point_Position(6,2)=domain_size_y-Gap_Between_Reflect_boundary
                Sampling_Point_Position(6,1)=Block_Center_X-0.50d0
                
                ! ( 4.0, 0.5)
                Sampling_Point_Position(7,2)=0.0d0+Gap_Between_Reflect_boundary
                Sampling_Point_Position(7,1)=Block_Center_X+2.50d0
                
                ! ( 4.0, 4.5)
                Sampling_Point_Position(8,2)=domain_size_y-Gap_Between_Reflect_boundary
                Sampling_Point_Position(8,1)=Block_Center_X+2.50d0

                ! ( 7.5, 0.5)
                Sampling_Point_Position(9,2)=0.0d0+Gap_Between_Reflect_boundary
                Sampling_Point_Position(9,1)=Block_Center_X+6.0d0
                
                ! ( 7.5, 4.0)
                Sampling_Point_Position(10,2)=domain_size_y/2.0
                Sampling_Point_Position(10,1)=Block_Center_X+6.0d0

                ! ( 7.5, 4.5)
                Sampling_Point_Position(11,2)=domain_size_y-Gap_Between_Reflect_boundary
                Sampling_Point_Position(11,1)=Block_Center_X+6.0d0

            endif
            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            ! Wave probe position
            if ( Wave_Probe_number>0 ) then

                Wave_Probe_Position(1,1)=0.3d0
                Wave_Probe_Position(2,1)=0.865d0
                Wave_Probe_Position(3,1)=1.114d0
                Wave_Probe_Position(4,1)=1.3625d0
                
            endif
            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            ! Sampling Line position
            SamplingLine:if ( SamplingLineNumber>0 ) then

                !----First Line----
                SamplingLineInitialPoint(1,1) = 0.25d0
                SamplingLineInitialPoint(1,2) = 0.0d0
                
                SamplingLineEndPoint(1,1) = 0.25d0
                SamplingLineEndPoint(1,2) = 1.0d0

                !----Second Line----
                SamplingLineInitialPoint(2,1) = 0.50d0
                SamplingLineInitialPoint(2,2) = 0.00d0

                SamplingLineEndPoint(2,1) = 0.50d0
                SamplingLineEndPoint(2,2) = 1.0d0

                !----Third Line----
                SamplingLineInitialPoint(3,1) = 0.75d0
                SamplingLineInitialPoint(3,2) = 0.00d0

                SamplingLineEndPoint(3,1) = 0.75d0
                SamplingLineEndPoint(3,2) = 1.00d0

                !----Fourth Line----
                SamplingLineInitialPoint(4,1) = 1.0d0
                SamplingLineInitialPoint(4,2) = 0.00d0

                SamplingLineEndPoint(4,1) = 1.0d0
                SamplingLineEndPoint(4,2) = 1.0d0

                !----Fifth Line----
                SamplingLineInitialPoint(5,1) = 2.5d0
                SamplingLineInitialPoint(5,2) = 0.00d0

                SamplingLineEndPoint(5,1) = 2.5d0
                SamplingLineEndPoint(5,2) = 1.0d0

                !Calculate the intervals
                do i=1,SamplingLineNumber
                    do j=1,dim
                       SamplingLinePointInterval(i,j)=(SamplingLineEndPoint(i,j)-SamplingLineInitialPoint(i,j))/(SamplingLinePointNumber-1)
                    enddo
                enddo

                !Calculate the sampling points on sampling line
                do k=1,SamplingLineNumber
                    do i=1,SamplingLinePointNumber
                       do j=1,dim
                          SamplingLinePointPosition(k,i,j)=SamplingLineInitialPoint(k,j)+(i-1)*SamplingLinePointInterval(k,j)
                       enddo
                    enddo
                enddo

            endif SamplingLine
            !--------------------------------------------------------------------------------------------------

        !******************************************************************************************************
        elseif(dim==3) then


        endif
        !******************************************************************************************************

        !******************************************************************************************************
        !Output the Sampling points on Sampling line
        if (Current_Processor_ID==Main_Processor) then
            
            do k=1,SamplingLineNumber

                !----------------------------------------------------------------------------------------------
                ! Transfer the sampling line index from integer to character
                write(char_file_index,'(I2)') k
                File_name  = "./Sampling/"//trim(adjustl(char_file_index))//"th_SamplingLine_Position.dat"

                call Initialziting_Writing_File( General_File_Port,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
          
                if(IOERROR==0) then

                    !Tecplot header
                    write(General_File_Port,*) "TITLE='DISTRIBUTION'"
                    write(General_File_Port,*) "VARIABLES= 'X' 'Y'"

                    do i=1,SamplingLinePointNumber
                        write(General_File_Port,'(6F12.6)') (SamplingLinePointPosition(k,i,j),j=1,dim)
                    end do

                endif
                close(General_File_Port)
                !----------------------------------------------------------------------------------------------

            enddo

        endif
        !******************************************************************************************************
 

    end subroutine Sampling_and_Wave_Probe_Initialization
    !==========================================================================================================











    !==========================================================================================================
    ! Check the set-up for SPH calculation which contains:(1) CFL checking for SPH time iteration 
    !                                                     (2) Local Reynolds number checking for inlet/outlet boundary condition
    !                                                     (3) Wave model checking for wave making
    subroutine Simulation_Setup_Checking()

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                                  ! Loop Variables

        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Simulation_Setup_Checking
        
        !******************************************************************************************************
        ! (1) CFL condition checking

        !calculate the CFL conditions
        !Morris (1997)
        CFL_Morris_dt=0.125*smooth_length**2/(1.01E-6)
        Min_CFL_dt(1)=CFL_Morris_dt

        !Monaghan (1989,1992)
        CFL_Monaghan_dt=0.25*smooth_length/(1.1*sqrt(square_c_0))
        Min_CFL_dt(2)=CFL_Monaghan_dt

        !Minimum CFL dt (0.1 is the coefficient)
        if (Iteration_Method>=3) then
            CFL_dt=2*0.1*minval(Min_CFL_dt)                               ! For Rugger-kuta ,time step can be a little larger                            
        else
            CFL_dt=0.1*minval(Min_CFL_dt)                                 ! For others ,keep same
        endif
        !******************************************************************************************************

        !******************************************************************************************************
        ! (2) Local Reynolds number checking for inlet/outlet boundary condition
        !Calculate the local Reynolds number (Re_h, Re_dx) for inlet/outlet boundary conditions
        Reynolds_Number_Dx = NonDimensionalization_Reference_Velocity*interior_dx  /water_kinematic_viscosity    
        Reynolds_Number_h  = NonDimensionalization_Reference_Velocity*smooth_length/water_kinematic_viscosity

        !******************************************************************************************************

        !******************************************************************************************************
        ! (3) Wave model checking for wave making 
        ! Assign the initial particle position (At t=0)
        Initial_particle_Position=particle_position
        
        ! Check the input wave information and search the suitable wave model
        ! All the processors should finished this work as wave information are changed in this subroutine
        ! The master processor responses for the output
        Wavemaker_Initialization: if ( trim(adjustl(wave_maker_type))=='WaveTheory' .and. make_wave_or_not==1 ) then
            
            ! Check the wave model is right or not, then recommand the best wave model for the wave generation
            call Wave_Model_Verification()

            !--------------------------------------------------------------------------------------------------
            ! Define the particles in wavemaker damping zone
            if (trim(adjustl(WavemakerDamping_Switch))=='on') then

                !Check the  wavemaker damping zone have changed or not
                if (abs(WavemakerDampingLengh-WavemakerDampingLenghScale*wave_length)<=1.0E-7) then
                    
                    if (Current_Processor_ID==1) then
                        write(*,*) "The length of wavemaker damping (",WavemakerDampingLengh,") zone is right"
                    endif

                    !Initialized the data for wavemaker damping 
                    WaveMaker_Analytical_position=particle_position      ! Analytical position in wavemaker damping zone
                    WaveMaker_Analytical_velocity=particle_velocity      ! Analytical velocity in wavemaker damping zone
                    WaveMaker_Analytical_rho     =particle_rho           ! Analytical density in wavemaker damping zone
                    WaveMaker_Analytical_press   =particle_press         ! Analytical pressure in wavemaker damping zone

                    !Define the initial particles in wavemaker damping zone and jeep these same in the calculation
                    In_WaveMakerDampingZone_OrNot=0                      ! Initialization 
                    FluidParticleNumberInWaveMakerDamping=0
                    position_soft_coefficient=0
                    do i=1,particle_ture_number

                        if (particle_position(i,1)<=WavemakerDampingLengh) then

                            !****This value doesn't need be refreshed during the calculation, it is constant during the iteration*****
                            In_WaveMakerDampingZone_OrNot(i)=1                   ! 1 is in; 0 is not in

                            FluidParticleNumberInWaveMakerDamping=FluidParticleNumberInWaveMakerDamping+1

                            position_soft_coefficient(i)=(1-(exp((particle_position(i,1)/WavemakerDampingLengh)**3.5)-1.0)/(exp(1.0)-1.0))

                        endif

                    enddo

                else
                    if (Current_Processor_ID==1) then
                       write(*,*) "The length of wavemaker damping (",WavemakerDampingLengh,") zone is not right"
                       write(*,*) "Actual wave damping length is :",WavemakerDampingLengh
                       write(*,*) "Desired wave damping length is :",WavemakerDampingLenghScale*wave_length
                    endif
                endif
                
            endif
            !--------------------------------------------------------------------------------------------------

        endif Wavemaker_Initialization
        !******************************************************************************************************


    end subroutine Simulation_Setup_Checking
    !==========================================================================================================





    !==========================================================================================================
    ! Subroutine Additional_Initialization
    ! Initialized the rest variables such as : (1) Smoke value
    !                                          (2) W(dx) (This part has been move to the initialization subroutine)
    !                                          (3) Identity_matrix
    !                                          (4) Initial total kinetic energy
    !                                          (5) analytical results such as for Green vortex case
    subroutine Additional_Initialization()

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                                  ! Loop Variables

        real(kind=8)::current_time                                          ! Current time
        character(len=20)::char_time                                        ! Time character type
        real(kind=8)::Velocity_Module                                       ! Velocity Module
        
        integer::File_index                                                 ! File index
        character(len=100)::File_name                                       ! File name ( The file name character lenghth should be 100 )
        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Additional_Initialization

        !******************************************************************************************************
        ! (1) Smoke value
        Smoke_Initialization: if ( DistributeInletOutletBoundaryParticleOrNot==1 .and. StartFrom==0 ) then

            !--------------------------------------------------------------------------------------------------
            ! Define the smoke value 
            particle_smoke_value=0.0d0
            
            do i=1,ture_total_particle_number

                ! All the particle except the fixed boundary particle
                if ( particle_type(i)==Fix_Ghost_Particle_Label ) cycle

                !----------------------------------------------------------------------------------------------
                ! Calculate the smoke value 
                particle_smoke_value(i)=(particle_position(i,dim)-domain_size_y*0.50)/(0.5*Smoke_Width)
                if ( particle_smoke_value(i)>1.0d0 ) then
                    particle_smoke_value(i)=1.0d0
                elseif ( particle_smoke_value(i)<-1.0d0 ) then
                    particle_smoke_value(i)=-1.0d0
                endif
                !----------------------------------------------------------------------------------------------

            enddo
            !--------------------------------------------------------------------------------------------------


            !--------------------------------------------------------------------------------------------------
            if (Current_Processor_ID==Main_Processor) then
                
                ! Open the Smoke Line Output file
                File_name  = "./Initial_Data/Initial_Smoke_Line.dat"
                call Initialziting_Writing_File( General_File_Port,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
          
                if(IOERROR==0) then

                    ! Tecplot header
                    write(General_File_Port,*) "TITLE='DISTRIBUTION'"
                    write(General_File_Port,*) "VARIABLES= 'X' 'Y' 'Smoke'"

                    ! Output fluid particle
                    write(General_File_Port,*) "ZONE I=",particle_ture_number," F=POINT"
                    do i=1,particle_ture_number
                         write(General_File_Port,100) (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),particle_smoke_value(i)
100                      format(12F10.4) 
                    end do
                    
                    ! Output wave maker particle
                    if (make_wave_or_not==1) then
                        write(General_File_Port,*) "ZONE I=",wave_maker_particle_number," F=POINT"
                        do k=1,wave_maker_particle_number
                            i=k+particle_ture_number
                            write(General_File_Port,100) (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),particle_smoke_value(i)
                        end do
                    endif

                    ! Output fix ghost boundary particle
                    if (DistributeFixedBoundaryParticleOrNot==1) then
                        write(General_File_Port,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
                        do k=1,fix_ghost_particle_number
                            i=k+particle_ture_number+wave_maker_particle_number
                            write(General_File_Port,100) (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),particle_smoke_value(i)
                        end do
                    endif

                    ! Output inlet and outlet boundary particles
                    if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then
                        write(General_File_Port,*) "ZONE I=",InletOutlet_Boundary_Particle_Number," F=POINT"
                        do k=1,InletOutlet_Boundary_Particle_Number
                            i=k+particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
                            write(General_File_Port,100) (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),particle_smoke_value(i)
                        end do
                    endif

                endif

                close(General_File_Port)

            endif
            !--------------------------------------------------------------------------------------------------
        
        endif Smoke_Initialization
        !******************************************************************************************************


        !******************************************************************************************************
        ! (2) W(dx) (This part has been move to the initialization subroutine)
        !Method 1: P.K. Stanby ----This is better for no free surface simulation

        !------------------------------------------------------------------------------------------------------
        ! Calculate w(dx) (This part has been move to the initialization subroutine)
        position_difference=0.0d0
        distance=interior_dx
        average_smooth_length=smooth_length
        
        call compute_kernel(Kernel_function_model,&          ! Kernel function model(in)   
                            dim,&                            ! Simulation dimension(in)   
                            distance,&                       ! Distance(in)                
                            position_difference,&            ! Position difference(in)
                            average_smooth_length,&          ! Smooth length(in)
                            WdxForShifting,&                 ! Kernel function(out)
                            temp_dwdx&                       ! Derivates of kernel function(out)
                            )

        !write(*,*) WdxForShifting,w(2),smooth_length,distance
        !------------------------------------------------------------------------------------------------------

        !******************************************************************************************************



        !******************************************************************************************************
        ! (3) Identity_matrix

        !------------------------------------------------------------------------------------------------------
        ! Initialize identity matrix
        identity_matrix=0.0d0
        do k=1,dim
            identity_matrix(k,k)=1.0d0
        end do
        !------------------------------------------------------------------------------------------------------

        !******************************************************************************************************


        !******************************************************************************************************
        ! (4) Initial total kinetic energy

        !------------------------------------------------------------------------------------------------------
        ! Calculate the initial kinetic energy calculation

        Initial_total_kinetic_energy=0.0d0

        do i=1,particle_ture_number

           Initial_total_kinetic_energy=Initial_total_kinetic_energy+0.5*particle_mass(i)*DOT_PRODUCT(particle_velocity(i,:),particle_velocity(i,:))
        
        enddo
        !------------------------------------------------------------------------------------------------------

        !******************************************************************************************************


        ! !******************************************************************************************************
        ! ! (5) Analytical results such as for Green vortex case

        ! !------------------------------------------------------------------------------------------------------
        ! ! Output the analytical results of Green vortex analytical results
        ! if (Current_Processor_ID==Main_Processor) then

        !     do k=0,5

        !         !----------------------------------------------------------------------------------------------
        !         current_time=k*0.1d0

        !         ! Calculate the analytical results
        !         do i=1,total_element_node_number

        !             element_node_velocity(i,2)=-U_Module*exp(-8.0*PI**2/Reynolds_Number*current_time)*cos(2*PI*element_node_position(i,1))*sin(2*PI*element_node_position(i,2))
        !             element_node_velocity(i,1)= U_Module*exp(-8.0*PI**2/Reynolds_Number*current_time)*sin(2*PI*element_node_position(i,1))*cos(2*PI*element_node_position(i,2))

        !             element_node_press(i)=0.25*water_rho_0*U_Module**2*exp(-16.0*PI**2/Reynolds_Number*current_time)*(cos(4*PI*element_node_position(i,1))+cos(4*PI*element_node_position(i,2)))
                
        !         end do
        !         !----------------------------------------------------------------------------------------------

        !         !----------------------------------------------------------------------------------------------
        !         ! Open file and output the analytical results
        !         ! Transfer the Current time from real to character
        !         write(char_time,'(f6.3)') current_time
        !         File_name = "./Initial_Data/"//trim(adjustl(char_time))//"_Analytical_Results.dat"

        !         call Initialziting_Writing_File( General_File_Port,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR

        !         if(IOERROR==0) then
                   
        !             !Tecplot header
        !             write(General_File_Port,*) "TITLE='DISTRIBUTION'"
        !             write(General_File_Port,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'U'"

        !             write(General_File_Port,*) 'ZONE T=" '//trim(adjustl(char_time))//' ",'//"N=",total_element_node_number,",E=",total_element_number,",DATAPACKING=POINT,ZONETYPE=FEQUADRILATERAL "   
                    
        !             !Domain name, number of nodes, number of cells, whether the data is on points or blocks (whether the data is points or blocks), cell type
        !             !Cell type: FETRIANGLE,FEQUADRILATERAL,FETETRAHEDRON,FEBRICK
        !             !For the hybrid grid, we should output the grid in different ZONES 
                    
        !             do i=1,total_element_node_number
        !                 Velocity_Module=sqrt( DOT_PRODUCT( element_node_velocity(i,:),element_node_velocity(i,:) ) )
        !                 write(General_File_Port,100) (element_node_position(i,j)/domain_size_y,j=1,dim),(element_node_velocity(i,j)/velocity_0,j=1,dim),element_node_press(i)/pressure_0,Velocity_Module/velocity_0
        !             end do
                    
        !             do i=1,total_element_number
        !                 write(General_File_Port,'(8I8)') (element_node_index(i,j),j=1,NodeNumberInOneElement)
        !             end do

        !         endif
               
        !         close( General_File_Port )
        !         !----------------------------------------------------------------------------------------------

        !     enddo

        ! endif
        ! !******************************************************************************************************

    end subroutine Additional_Initialization
    !==========================================================================================================









    !==========================================================================================================
    ! Subroutine Initialization_Result_Outputting
    ! Output SPH simulation base information and initial particle results
    subroutine Initialization_Result_Outputting()
        
        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                                  ! Loop Variables
        
        integer::File_index                                                 ! File index
        character(len=100)::File_name                                       ! File name ( The file name character lenghth should be 100 )
        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Initialization_Result_Outputting

        !******************************************************************************************************
        ! (1) Output base information of SPH simulation

        ! Output the end notice on screen 
        write(*,'(A)') " ====================================================================================="
        
        !------------------------------------------------------------------------------------------------------
        !This one is very important
        !As the intel complier cannot note the array out of the boundary in some conditions, so we should pay attention to the array memory
        !The n_total is very important, some conditions the boundary particles will be larger than the interior nodes, this may make errors
        !So we check the memory of the array
        write(*,'(A,F8.4,A)') " Running memory checking                         :"
        write(*,'(A,F6.3,A)') " Recommand value of 'n_total_Factor' is          :",1.1*(particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number)/particle_ture_number

        if (particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number>n_total) then

            write(*,'(A,I8,A)')   " The allocated memory size in partcile number is :",n_total
            write(*,'(A,I8,A)')   " The safety memory size in partcile number is    :",particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number  
            
            write(*,'(A)')        " *Errors may happen during the calculation, please increase the memory for this case*"

        else

            write(*,'(A,I8,A)')   " The allocated memory size in partcile number is :",n_total
            write(*,'(A,I8,A)')   " The current memory size in partcile number is   :",particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
            write(*,'(A,F5.2,A,F5.2,A)')   " The n_total_Factor(",n_total_Factor,") larger than recommanded value (",1.1*(particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number)/particle_ture_number,")"
            write(*,'(A)')        " The memory for the calculation is enough !"

        endif
        write(*,*)
        !------------------------------------------------------------------------------------------------------


        write(*,'(A,F8.4,A)') " Recommand 'dt' for iteration  :",CFL_dt,"s"
        write(*,*)
        write(*,'(A,F9.2,A)') " Simulation Reynolds Number    : ",Reynolds_Number
        write(*,*)
        write(*,'(A,F8.2,A)') " Reynolds-Cell Number          : ",Reynolds_Number_dx
        write(*,'(A,F8.2,A)') " Reynolds-Smooth Length Number : ",Reynolds_Number_h
        write(*,'(A,F8.4,A)') " Both Reynolds-Cell and Reynolds-Smooth Length Number should be less than 10"

        write(*,*)
        write(*,'(A,F8.4,A)') " Basic information of current SPH simulation : "
        write(*,'(2A)')            " Kernel Function    : ",trim(adjustl(Kernel_function_name))
        write(*,'(A,F9.6,A,F9.6)') " Particle resolution:",intial_volume_0**(1.0d0/dim)," Particle smooth length:",smooth_length
        write(*,'(A,F9.6,A,F9.6)') " Particle size      :",intial_volume_0,             " Particle mass         :",water_particle_mass_0
        write(*,*)

        write(*,'(A,F8.4,A)') " Post proceeding information     :"
        write(*,'(A,F8.4,A)') " Desired particle intervel dx    :",interior_dx
        write(*,'(A,F8.4,A)') " Post proceeding grid intervel dx:",element_dx
        write(*,*)

        write(*,'(A,I10)')     " Initial particle number summarization:"
        write(*,'(A,I10)')     " Water particle number                :",Particle_ture_number
        
        if ( make_wave_or_not==1 ) then
            write(*,'(A,I10)') " Wave maker particle number           :",Wave_maker_particle_number    
        endif

        if ( DistributeFixedBoundaryParticleOrNot==1 ) then
            write(*,'(A,I10)') " Fixed ghost particle number          :",Fix_ghost_particle_number    
        endif

        if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then
            write(*,'(A,I10)') " Inlet/Outlet boundary particle number:",InletOutlet_Boundary_Particle_Number    
        endif

        write(*,'(A)') " ====================================================================================="

        ! Output the end notice to log file
        write(Log_File_Port,'(A)') " ====================================================================================="
        
        !------------------------------------------------------------------------------------------------------
        write(Log_File_Port,'(A,F8.4,A)') " Running memory checking                         :"
        write(Log_File_Port,'(A,F6.3,A)') " Recommand value of 'n_total_Factor' is          :",1.1*(particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number)/particle_ture_number

        if ( particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number >= n_total) then

            write(Log_File_Port,'(A,I8,A)')   " The allocated memory size in partcile number is :",n_total
            write(Log_File_Port,'(A,I8,A)')   " The safety memory size in partcile number is    :",particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number  
            
            write(Log_File_Port,'(A)')        " *Errors may happen during the calculation, please increase the memory for this case*"

        else

            write(Log_File_Port,'(A,I8,A)')   " The allocated memory size in partcile number is :",n_total
            write(Log_File_Port,'(A,I8,A)')   " The current memory size in partcile number is   :",particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
            write(Log_File_Port,'(A,F5.2,A,F5.2,A)')   " The n_total_Factor(",n_total_Factor,") larger than recommanded value (",1.1*(particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number)/particle_ture_number,")"
            write(Log_File_Port,'(A)')        " The memory for the calculation is enough !"

        endif
        write(*,*)
        !------------------------------------------------------------------------------------------------------


        write(Log_File_Port,'(A,F8.4,A)') " Recommand 'dt' for iteration  :",CFL_dt,"s"
        write(Log_File_Port,*)
        write(Log_File_Port,'(A,F9.2,A)') " Simulation Reynolds Number    : ",Reynolds_Number
        write(Log_File_Port,*)
        write(Log_File_Port,'(A,F8.2,A)') " Reynolds-Cell Number          : ",Reynolds_Number_Dx
        write(Log_File_Port,'(A,F8.2,A)') " Reynolds-Smooth Length Number : ",Reynolds_Number_h
        write(Log_File_Port,'(A,F8.4,A)') " Recommandation: Both Reynolds-Cell and Reynolds-Smooth Length Number should be less than 10"
        write(Log_File_Port,*)
        write(Log_File_Port,'(A,F8.4,A)') " Basic information of current SPH simulation : "
        write(Log_File_Port,'(2A)')            " Kernel Function    : ",trim(adjustl(Kernel_function_name))
        write(Log_File_Port,'(A,F9.6,A,F9.6)') " Particle resolution:",intial_volume_0**(1.0d0/dim)," Particle smooth length:",smooth_length
        write(Log_File_Port,'(A,F9.6,A,F9.6)') " Particle size      :",intial_volume_0,             " Particle mass         :",water_particle_mass_0
        write(Log_File_Port,*)
            
        write(Log_File_Port,'(A,I10)') " Initial particle number summarization:"
        write(Log_File_Port,'(A,I10)')     " Water particle number                :",Particle_ture_number
        
        if ( make_wave_or_not==1 ) then
            write(Log_File_Port,'(A,I10)') " Wave maker particle number           :",Wave_maker_particle_number    
        endif

        if ( DistributeFixedBoundaryParticleOrNot==1 ) then
            write(Log_File_Port,'(A,I10)') " Fixed ghost particle number          :",Fix_ghost_particle_number    
        endif

        if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then
            write(Log_File_Port,'(A,I10)') " Inlet/Outlet boundary particle number:",InletOutlet_Boundary_Particle_Number    
        endif

        write(Log_File_Port,'(A)') " ====================================================================================="
        
        !******************************************************************************************************




        !******************************************************************************************************
        ! (2) Output initial particle results
        ! Open the Output file
        
        File_name  = "./Initial_Data/Particle_Information.dat"
        call Initialziting_Writing_File( General_File_Port,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
  
        if(IOERROR==0) then

            ! Output tecplot header
            write(General_File_Port,*) "TITLE='DISTRIBUTION'"
            write(General_File_Port,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
            ! Output fluid particle
            write(General_File_Port,*) "ZONE I=",particle_ture_number," F=POINT"
            do i=1,particle_ture_number
                 write(General_File_Port,100) (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
    100          format(8F10.4) 
            end do
            
            ! Output wave maker particle
            if (make_wave_or_not==1) then
                write(General_File_Port,*) "ZONE I=",wave_maker_particle_number," F=POINT"
                do k=1,wave_maker_particle_number
                    i=k+particle_ture_number
                    write(General_File_Port,100) (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                end do
            endif

            ! Output fix ghost boundary particle
            if (DistributeFixedBoundaryParticleOrNot==1) then
                write(General_File_Port,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
                do k=1,fix_ghost_particle_number
                    i=k+particle_ture_number+wave_maker_particle_number
                    write(General_File_Port,100) (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                end do
            endif

            ! Output inlet and outlet boundary particles
            if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then

                write(General_File_Port,*) "ZONE I=",InletOutlet_Boundary_Particle_Number," F=POINT"
                do k=1,InletOutlet_Boundary_Particle_Number
                    i=k+particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
                    write(General_File_Port,100) (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                end do
                
            endif

        end if

        close(General_File_Port)
        !******************************************************************************************************

    end subroutine Initialization_Result_Outputting
    !==========================================================================================================









    !==========================================================================================================
    ! Subroutine Display_Running_Information_On_Screen
    ! Output SPH simulation base information on screen
    subroutine Display_Running_Information_On_Screen()
        
        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                                  ! Loop Variables
        
        integer::File_index                                                 ! File index
        character(len=100)::File_name                                       ! File name ( The file name character lenghth should be 100 )
        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Display_Running_Information_On_Screen

        !******************************************************************************************************
        ! (1) Output base information of SPH simulation

        ! Output the end notice on screen 
        write(*,'(A)') " ======================================================================================="
        write(*,'(A,F9.2,A)') " Simulation Reynolds Number    : ",Reynolds_Number
        write(*,'(A,F8.2,A)') " Reynolds-Cell Number          : ",Reynolds_Number_dx
        write(*,'(A,F8.2,A)') " Reynolds-Smooth Length Number : ",Reynolds_Number_h

        write(*,*)

        write(*,'(A,I10)')     " Present particle number summarization      :"
        write(*,'(A,I10)')     " Water particle number                      :",Particle_ture_number
        
        if ( make_wave_or_not==1 ) then
            write(*,'(A,I10)') " Wave maker particle number                 :",Wave_maker_particle_number    
        endif

        if ( DistributeFixedBoundaryParticleOrNot==1 ) then
            write(*,'(A,I10)') " Fixed ghost particle number                :",Fix_ghost_particle_number    
        endif

        if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then
            write(*,'(A,I10)') " Inlet/Outlet boundary particle number      :",InletOutlet_Boundary_Particle_Number    
        endif

        if ( DistributePeriodBoundaryParticleOrNot==1 ) then
            write(*,'(A,I10)') " Period bounary particle number             :",Period_Boundary_Particle_Number    
        endif

        if ( DistributeTraditionalGhostBoundaryParticleOrNot==1 ) then
            write(*,'(A,I10)') " Traditional ghost boundary particle number :",Traditional_ghost_particle_number    
        endif


        if ( Calculate_Body_Force_or_Not==1 ) then

            write(*,*)
            write(*,'(A,I10)')    " Present force summarization :"
            write(*,'(A,3F12.6)') " Pressure force              :", (Body_Pressure_Force(i)/force_0,i=1,dim)
            write(*,'(A,3F12.6)') " Viscous force               :", (Body_Viscous_Force(i) /force_0,i=1,dim)
            write(*,'(A,3F12.6)') " Total force                 :", (Body_Total_Force(i)   /force_0,i=1,dim)
            
        endif
        write(*,'(A)') " ======================================================================================="

        !******************************************************************************************************


    end subroutine Display_Running_Information_On_Screen
    !==========================================================================================================









    !**********************************************************************************************************
    !
    !  SUBROUTINE : Transfer_External_Boundary(i_time_step) 
    !
    !  PURPOSE    : (1) Transfer the external boundary particles during the running time
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
    !**********************************************************************************************************
    subroutine  Transfer_External_Boundary(i_time_step) 

        implicit none

        !======================================================================================================
        ! Variables from the superior-subroutien
        integer,intent(in)::i_time_step                                     ! Current time step 

        ! Variables in local subroutine
        integer::i,j,k,L,m                                                  ! Variables for looping

        integer::m_virtual,n_virtual                                        ! Virtual particle order
        
        ! Read file varibales
        integer::File_index                                                 ! File index
        character(len=100)::File_name                                       ! File name ( The file name character lenghth should be 100 )
            
        ! Particle packing varibales
        real(kind=8)::x_0,y_0,x_1,y_1                                       ! Starting and edn point on angular bisector
        real(kind=8)::fix_particle_distance                                 ! Real fix particle distance
        real(kind=8)::position_x,position_y,position_z                      ! Partcile temp coordiantes
        real(kind=8)::line_distance                                         ! Line distance
        real(kind=8)::particle_x,particle_y


        !======================================================================================================


        ! Body of subroutine Transfer_External_Boundary


        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! Calculate the present particle number, only need to save the in/outlet boundary partilces 
        ! If no in/outlet boundary stop running, delete all the InOutlet boundary particles
        if ( DistributeInletOutletBoundaryParticleOrNot==0 ) then

            InletOutlet_Boundary_Particle_Number =0

        endif

        Ture_total_particle_number=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number
        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! Distribute the fixed boundary particles at the external boundary
        n_virtual=0
        !------------------------------------------------------------------------------------------------------
        ! Bottom fixed boundary particles
        External_Bottom_Boundary: do i=1,fix_ghost_layer                                   ! Fixed ghost layer
            
            x_0 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_0 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            x_1 = boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_1 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            
            ! Define the initial number in first row
            if (i==1) then
               j=0
               do 
                  j=j+1                                                                    ! Virtual particle number
                  position_x=(j-1)*fix_ghost_dx+x_0                                        ! Virtual particle x coordinate
                  position_y=y_0                                                           ! Virtual particle y coordinate
                
                  ! Exceed the angular bisector yes or no
                  if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit
                         
               end do
               k=j                       ! The virtual fixed boundary particle number of first row 
               
            end if
            
            ! Generate all Bottom particles for each row
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                    ! Has K points and (k-1) interval
            
            do j=1,k-1                   ! From 1 to k-1 points,as the last points will be repeated in next interval
                
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Ture_total_particle_number                             ! Virtual particle index
                particle_position(m_virtual,1)=(j-1)*fix_particle_distance+x_0             ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0                                         ! Virtual particle y coordinate
                free_surface_type(m_virtual)=0
                
                particle_type(m_virtual)=Fix_Ghost_Particle_Label                          ! Fixed Ghost Particle Label
                particle_initial_type(m_virtual)=Fix_Ghost_Particle_Label                  ! Fixed Ghost Particle Label

                ! Form 0 to boundary_size_x is Horizontial mirror rule. others are Center mirror
                if ( 0.0d0 <= particle_position(m_virtual,1) .and. particle_position(m_virtual,1) <= boundary_size_x ) then
                    Boundary_particle_type(m_virtual)=-6                                   ! Horizontial type
                else
                    Boundary_particle_type(m_virtual)=-7                                   ! Center type
                endif
                
            end do
            
            k=k+2                        ! The virtual particle number of next row              
            
        end do External_Bottom_Boundary
        !-------------------------------------------------------------------------------------------------------


        !-------------------------------------------------------------------------------------------------------
        ! Bottom fixed boundary particles
        External_Right_Boundary: do i=1,fix_ghost_layer                                    ! Fixed ghost layer
            
            x_0 = boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_0 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            x_1 = boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_1 = boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            
            ! Define the initial number in first row
            if (i==1) then
               j=0
               do 
                  j=j+1                                                                    ! Virtual particle number
                  position_x=x_0                                                           ! Virtual particle x coordinate
                  position_y=y_0+(j-1)*fix_ghost_dx                                        ! Virtual particle y coordinate
                
                  ! Exceed the angular bisector yes or no
                  if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit
                         
               end do
               k=j                       ! The virtual fixed boundary particle number of first row 
               
            end if
            
            ! Generate all Bottom particles for each row
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                    ! Has K points and (k-1) interval
            
            do j=1,k-1                   ! From 1 to k-1 points,as the last points will be repeated in next interval
                
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Ture_total_particle_number                             ! Virtual particle index
                particle_position(m_virtual,1)=x_0                                         ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0+(j-1)*fix_particle_distance             ! Virtual particle y coordinate
                free_surface_type(m_virtual)=0
                
                particle_type(m_virtual)=Fix_Ghost_Particle_Label                          ! Fixed Ghost Particle Label
                particle_initial_type(m_virtual)=Fix_Ghost_Particle_Label                  ! Fixed Ghost Particle Label

                ! Form 0 to boundary_size_x is Vertical mirror rule. others are Center mirror
                if ( 0.0d0 <= particle_position(m_virtual,2) .and. particle_position(m_virtual,2) <= boundary_size_y ) then
                    Boundary_particle_type(m_virtual)=-8                                   ! Vertical type
                else
                    Boundary_particle_type(m_virtual)=-7                                   ! Center type
                endif
                
            end do
            
            k=k+2                        ! The virtual particle number of next row              
            
        end do External_Right_Boundary
        !-------------------------------------------------------------------------------------------------------


        !-------------------------------------------------------------------------------------------------------
        ! Top fixed boundary particles (From right to left)
        External_Top_Boundary: do i=1,fix_ghost_layer                                      ! Fixed ghost layer
            
            x_0 = boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_0 = boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            x_1 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_1 = boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0

            ! write(*,*) x_1,y_1
            ! Define the initial number in first row
            if (i==1) then
               j=0
               do 
                   j=j+1                                                                   ! Virtual particle index
                   position_x=x_0-(j-1)*fix_ghost_dx                                       ! Virtual particle x coordinate  
                   position_y=y_0                                                          ! Virtual particle y coordinate
                   
                   ! Define the initial number in first row
                   if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit

                end do
                k=j                        ! The virtual particle number of first row 
                
            end if
            
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                    ! k points and k-1 intervals
            
            do j=1,k-1                     ! 1 to k-1 points (Delte the end point)
                
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Ture_total_particle_number                             ! Virtual particle index

                particle_position(m_virtual,1)=x_0-(j-1)*fix_particle_distance             ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0                                         ! Virtual particle y coordinate
                free_surface_type(m_virtual)=0

                particle_type(m_virtual)=Fix_Ghost_Particle_Label                          ! Fixed ghost particle Label
                particle_initial_type(m_virtual)=Fix_Ghost_Particle_Label                  ! Fixed ghost pParticle Label

                ! Form 0 to boundary_size_x is Horizontial mirror rule. others are Center mirror
                if ( 0.0d0 <= particle_position(m_virtual,1) .and. particle_position(m_virtual,1) <= boundary_size_x ) then
                    Boundary_particle_type(m_virtual)=-6                                   ! Horizontial type
                else
                    Boundary_particle_type(m_virtual)=-7                                   ! Center type
                endif

            end do
            
            k=k+2                         ! The virtual particle number of first row                

        end do External_Top_Boundary
        !-------------------------------------------------------------------------------------------------------
        


        !-------------------------------------------------------------------------------------------------------
        ! Top fixed boundary particles (From right to left)
        External_Left_Boundary: do i=1,fix_ghost_layer                                      ! Fixed ghost layer
            
            x_0 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_0 = boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            x_1 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_1 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0

            ! write(*,*) x_1,y_1
            ! Define the initial number in first row
            if (i==1) then

               j=0
               do 
                   j=j+1                                                                   ! Virtual particle index
                   position_x=x_0                                                          ! Virtual particle x coordinate  
                   position_y=y_0-(j-1)*fix_ghost_dx                                       ! Virtual particle y coordinate
                   
                   ! Define the initial number in first row
                   if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit

                end do
                k=j                        ! The virtual particle number of first row 
                
            end if
            
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                    ! k points and k-1 intervals
            
            do j=1,k-1                     ! 1 to k-1 points (Delte the end point)
                
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Ture_total_particle_number                             ! Virtual particle index

                particle_position(m_virtual,1)=x_0                                         ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0-(j-1)*fix_particle_distance             ! Virtual particle y coordinate
                
                free_surface_type(m_virtual)=0

                particle_type(m_virtual)=Fix_Ghost_Particle_Label                          ! Fixed ghost particle Label
                particle_initial_type(m_virtual)=Fix_Ghost_Particle_Label                  ! Fixed ghost pParticle Label

                ! Form 0 to boundary_size_x is Horizontial mirror rule. others are Center mirror
                if ( 0.0d0 <= particle_position(m_virtual,2) .and. particle_position(m_virtual,2) <= boundary_size_y ) then
                    Boundary_particle_type(m_virtual)=-8                                   ! Horizontial type
                else
                    Boundary_particle_type(m_virtual)=-7                                   ! Center type
                endif

            end do
            
            k=k+2                         ! The virtual particle number of first row                

        end do External_Left_Boundary
        !-------------------------------------------------------------------------------------------------------



        New_Fix_ghost_particle_number = n_virtual

        ! Initialize other information for fixed ghost boundary particles
        do i=1,New_Fix_ghost_particle_number
            
            m_virtual=i+Ture_total_particle_number                                           ! Virtual particle index
             
            particle_velocity(m_virtual,:)=0.0                                               ! Particle velocity
            particle_rho(m_virtual)=water_rho_0                                              ! Particle rho
            particle_mass(m_virtual)=water_particle_mass_0                                   ! Particle mass
            particle_smooth_lengh(m_virtual)=smooth_length                                   ! Particle smooth lengh
            particle_press(m_virtual)=0.0d0+Background_Pressure                              ! Particle press
            particle_c(m_virtual)=sqrt(square_c_0)                                           ! Particle c
            particle_energy(m_virtual)=e_0                                                   ! Particle energy
            particle_volume(m_virtual)=intial_volume_0                                       ! Particle volume

        end do
        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! Checking the boundary transfer results

        if (Current_Processor_ID==Main_Processor) then

            ! Open the Output file
            File_index = General_File_Port
            File_name  = "./Initial_Data/Boundary_Transfer_Results_Checking.dat"
            call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
      
            if(IOERROR==0) then

                ! Input tecplot header
                write(File_index,*) "TITLE='DISTRIBUTION'"
                write(File_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
                
                ! Output fluid particle
                write(File_index,*) "ZONE I=",Particle_ture_number," F=POINT"
                do i=1,Particle_ture_number
                    write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
100                 format(8F10.4) 
                enddo
                
                ! Output wave maker particle
                if (make_wave_or_not==1) then
                    write(File_index,*) "ZONE I=",Wave_maker_particle_number," F=POINT"
                    do k=1,Wave_maker_particle_number
                        i=k+Particle_ture_number
                        write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                    end do
                endif

                ! Output fix ghost boundary particle
                if (DistributeFixedBoundaryParticleOrNot==1) then
                    write(File_index,*) "ZONE I=",Fix_ghost_particle_number," F=POINT"
                    do k=1,Fix_ghost_particle_number
                        i=k+Particle_ture_number+Wave_maker_particle_number
                        write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                    end do
                endif

                ! Inlet and outlet boundary particles
                if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then

                    write(File_index,*) "ZONE I=",InletOutlet_Boundary_Particle_Number," F=POINT"
                    do k=1,InletOutlet_Boundary_Particle_Number
                        i=k+Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number
                        write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                    end do
                    
                endif

                ! New fixed boundary particles
                write(File_index,*) "ZONE I=",New_Fix_ghost_particle_number," F=POINT"
                do k=1,New_Fix_ghost_particle_number
                    i=k+Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number
                    write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                end do
                      
            endif
            
            close(File_index)

        endif
        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



        !------------------------------------------------------------------------------------------------------



    end subroutine Transfer_External_Boundary
    !==========================================================================================================








    !**************************************************************************************************************
    !  SUBROUTINE: MLS_Interpolation_In_Gobal_Domain
    !
    !  Purpose: interpolation calculation for particle in full domian 
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
    !  Note: MPI version: mpich-3.2
    !        Fortran version: Inter Fortran (ifort) 18.0.0
    !**************************************************************************************************************

    !==============================================================================================================
    subroutine MLS_Interpolation_In_Gobal_Domain(MLSOrder,Interpolate_particle_position,Ture_effect_number,Ture_effect_name,W_kernel,Fai)


        implicit none

        !==========================================================================================================
        ! Variables
        ! Variables from superior subroutine
        integer,intent(in)::MLSOrder                                                                      ! MLSOrder
        real(kind=8),dimension(dim),intent(in)::Interpolate_particle_position                             ! Interpolation particle position
        integer,intent(in)::Ture_effect_number                                                            ! Ture effect particle number
        integer,dimension(PredictionNumberInSupportDomain),intent(in)::Ture_effect_name                   ! Ture effect particle index
        real(kind=8),dimension(PredictionNumberInSupportDomain),intent(in)::W_kernel                      ! Weight function value
        real(kind=8),dimension(PredictionNumberInSupportDomain),intent(out)::Fai                          ! Output fai value
        
        ! Variables in subroutine
        integer::i,j,k,l,m                                                                                ! Variables for looping
        
        ! MLS Variables
        real(kind=8),dimension(PredictionNumberInSupportDomain,PredictionNumberInSupportDomain)::w_matrix ! Weight function value matrix 
        real(kind=8),dimension(MLSOrder,PredictionNumberInSupportDomain)::P                               ! Basis P  
        real(kind=8),dimension(MLSOrder,MLSOrder)::Modified_matrix                                        ! Modified identity matrix 
        real(kind=8),dimension(MLSOrder,1)::P_i                                                           ! P matrix of i (Dimension of P is : MLSOrder*1)
        real(kind=8),dimension(MLSOrder,1)::Q_i                                                           ! Q matrix of i (Dimension of P is : MLSOrder*1)
        real(kind=8),dimension(MLSOrder,PredictionNumberInSupportDomain)::B                               ! Basis B (Dimension of B is : MLSOrder*k)  
        real(kind=8),dimension(MLSOrder,MLSOrder)::A,A_inversion                                          ! A matrix (Dimension of A and A_inversion is : MLSOrder*MLSOrder)
        real(kind=8),dimension(MLSOrder,PredictionNumberInSupportDomain)::A_inversion_times_B
        real(kind=8),dimension(1,1)::Temp_Fai
        real(kind=8)::sum_w_kernel

        integer::check_error
        !==========================================================================================================


        ! Body of subroutine MLS_Interpolation


        !==========================================================================================================
        ! Weight function value matrix
        w_matrix=0.0d0                                       ! (Dimension of W is : k*k)
        do k=1,Ture_effect_number 
            w_matrix(k,k)=W_kernel(k)
        end do
            
        ! Assign the Modified_matrix
        Modified_matrix=0.0d0
        do k=1,MLSOrder
            Modified_matrix(k,k)=slighted_modified_alpha 
        end do

        ! Initialization
        P=0.0d0
        B=0.0d0
        A=0.0d0

        check_error=1

        !----------------------------------------------------------------------------------------------------------        
        ! MLS interpolation
        ! Evaluate basis P, B Matrix and their derivates  (m is the order of MLS)
        do k=1,Ture_effect_number
    
            ! effect particle index
            j=Ture_effect_name(k)
        
            if(MLSOrder==1) then
    
                P(1,k)=1                                       ! P_i=(1)            (Dimension of P is : k*MLSOrder)
    
            else if(MLSOrder==3) then
                                
                P(1,k)=1                                       ! P_i=(1,x_j/h,y_j/h)
                P(2,k)=(particle_position(j,1)-Interpolate_particle_position(1))/smooth_length
                P(3,k)=(particle_position(j,2)-Interpolate_particle_position(2))/smooth_length
        
            else if(MLSOrder==6) then
    
                P(1,k)=1                                       ! P_i=(1,x_j/h,y_j/h,(x_j/h)^2,x_j/h*y_j/h,(y_j/h)^2)
                P(2,k)=(particle_position(j,1)-Interpolate_particle_position(1))/smooth_length
                P(3,k)=(particle_position(j,2)-Interpolate_particle_position(2))/smooth_length
                P(4,k)=P(2,k)**2
                P(5,k)=P(2,k)*P(3,k)
                P(6,k)=P(3,k)**2
    
            else
                write(*,'(A)') " The order of MLS input is not right , Please ckeck 'MLSOrder' value ! (Interpolation Gobal)" 
            end if
    
        enddo
        !----------------------------------------------------------------------------------------------------------
         
        !----------------------------------------------------------------------------------------------------------
        ! Calculate P_i and its derivates
        if(MLSOrder==1) then
        
            P_i(1,1)=1                                       ! P_i=(1)            (Dimension of P is : MLSOrder*1)
        
        else if(MLSOrder==3) then
                                
            P_i(1,1)=1                                       ! P_i=(1,x_i,y_i)
            P_i(2,1)=0.0d0
            P_i(3,1)=0.0d0
        
        else if(MLSOrder==6) then
        
            P_i(1,1)=1                                       ! P_i=(1,x_i,y_i,x_i^2,x_i*y_i,y_i^2)
            P_i(2,1)=0.0d0
            P_i(3,1)=0.0d0
            P_i(4,1)=0.0d0
            P_i(5,1)=0.0d0
            P_i(6,1)=0.0d0
        
        else
            write(*,'(A)') " The order of MLS input is not right , Please ckeck 'MLSOrder' value ! (Interpolation Gobal)" 
        end if
        !----------------------------------------------------------------------------------------------------------
       
        !----------------------------------------------------------------------------------------------------------
        ! Modified the P and W matrix (Makesure non-singal happens)
        if(MLSOrder/=1) then
        
            ! Modified P: New_P=[P Modified_matrix]   (dimension: m*(k+m) )
            ! Modified W: New_W=[ W               0           ]
            !                   [ 0       Modified_matrix     ]                      (dimension: (k+m)*(k+m) )
            do m=1,MLSOrder
                do L=1,MLSOrder
                    k=Ture_effect_number+L
                    P(m,k)=Modified_matrix(m,L)                     ! Comble P and Modified_matrix( as the add matrix is constant so its derivated keep same)
                end do
            end do 
        
            do L=1,MLSOrder
                 k=Ture_effect_number+L
                 W_matrix(k,k)=Modified_matrix(L,L)                 ! Comble P and Modified_matrix( as the add matrix is constant so its derivated keep same)
            end do
        
        end if
        !----------------------------------------------------------------------------------------------------------
    
        !----------------------------------------------------------------------------------------------------------
        ! Calculate B,  A_inversion
        ! Assign the B=P*W                
        B=matmul(P,W_matrix)                          ! (Dimension of B is : MLSOrder*k)
            
        ! Evaluate matrices A=P*W*transpose(P) and its derivatives
        A=matmul(B,transpose(P))                      ! (Dimension of A is : MLSOrder*MLSOrder)
    
        check_error=1

        if(MLSOrder==1) then
            
            A_inversion(1,1)=1.0/A(1,1)
    
        else if(MLSOrder==3) then
            
            call matrix_inversion_less_three(MLSOrder,A,A_inversion,check_error)

            if(check_error==0) then
                write(*,*) "The A_inversion compute is wrong! (Interpolation Gobal)"
            end if
        
        else if(MLSOrder==6) then 
                    
            A_inversion=A                                          !Initialized the A inversion matrix
            !Calculate the inverse  matrix of A
            !call matrix_inversion_less_six(MLSOrder,A,A_inversion,check_error)
            call BRINV(A_inversion,MLSOrder,check_error)
            
            if(check_error==0) then
                write(*,*) "The A_inversion compute is wrong! (Interpolation Gobal)"
            end if
    
        end if
        !----------------------------------------------------------------------------------------------------------

        !----------------------------------------------------------------------------------------------------------
        ! Check inversion calculation is right or not
        if(check_error/=0) then

            !******************************************************************************************************
            ! Calculate A_inversion times B  and its derivates (dimension: MLSOrder*k)
            A_inversion_times_B=matmul(A_inversion,B)
     
            Q_i=P_i
            do k=1,Ture_effect_number
    
                ! Calculate P_i and B_i   
                Temp_Fai=0.0d0
                do m=1,MLSOrder
                
                    Temp_Fai=Temp_Fai+Q_i(m,1)*A_inversion_times_B(m,k)
              
                end do
        
                Fai(k)=Temp_Fai(1,1)

            end do
            !******************************************************************************************************

        ! For the case, the inverse oricedure is not right
        elseif(check_error==0) then

            sum_w_kernel=0.0d0
            do k=1,Ture_effect_number
                sum_w_kernel=sum_w_kernel+W_kernel(k)  
            end do
 
            do k=1,Ture_effect_number
                Fai(k)=W_kernel(k)/sum_w_kernel
            end do
        
        end if
        !----------------------------------------------------------------------------------------------------------

        !----------------------------------------------------------------------------------------------------------
        !For high order calculation, when the near particle is not enough, decrease the MLS Order
        if(MLSOrder==6 .and. Ture_effect_number<=8) then

            sum_w_kernel=0.0d0
            do k=1,Ture_effect_number
                sum_w_kernel=sum_w_kernel+W_kernel(k)  
            end do
 
            do k=1,Ture_effect_number
                Fai(k)=W_kernel(k)/sum_w_kernel
            end do
        
        end if
        !----------------------------------------------------------------------------------------------------------
        
    end subroutine MLS_Interpolation_In_Gobal_Domain
    !==============================================================================================================





    !**************************************************************************************************************
    !  SUBROUTINE: MLS_Interpolation_In_Subdomain_Domain
    !
    !  Purpose: interpolation calculation for particle in subdomain domian 
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
    !  Note: MPI version: mpich-3.2
    !        Fortran version: Inter Fortran (ifort) 18.0.0
    !**************************************************************************************************************

    !==============================================================================================================
    subroutine MLS_Interpolation_In_Subdomain_Domain(MLSOrder,Interpolate_subdomain_particle_position,Ture_effect_number,Ture_effect_name,W_kernel,Fai)

        implicit none

        !==========================================================================================================
        ! Variables
        ! Variables from superior subroutine
        integer,intent(in)::MLSOrder                                                                      ! MLSOrder
        real(kind=8),dimension(dim),intent(in)::Interpolate_subdomain_particle_position                   ! Interpolation particle position
        integer,intent(in)::Ture_effect_number                                                            ! Ture effect particle number
        integer,dimension(PredictionNumberInSupportDomain),intent(in)::Ture_effect_name                   ! Ture effect particle index
        real(kind=8),dimension(PredictionNumberInSupportDomain),intent(in)::W_kernel                      ! Weight function value
        real(kind=8),dimension(PredictionNumberInSupportDomain),intent(out)::Fai                          ! Output fai value
        
        ! Variables in subroutine
        integer::i,j,k,l,m                                                                                ! Variables for looping
        
        ! MLS Variables
        real(kind=8),dimension(PredictionNumberInSupportDomain,PredictionNumberInSupportDomain)::w_matrix ! Weight function value matrix 
        real(kind=8),dimension(MLSOrder,PredictionNumberInSupportDomain)::P                               ! Basis P  
        real(kind=8),dimension(MLSOrder,MLSOrder)::Modified_matrix                                        ! Modified identity matrix 
        real(kind=8),dimension(MLSOrder,1)::P_i                                                           ! P matrix of i (Dimension of P is : MLSOrder*1)
        real(kind=8),dimension(MLSOrder,1)::Q_i                                                           ! Q matrix of i (Dimension of P is : MLSOrder*1)
        real(kind=8),dimension(MLSOrder,PredictionNumberInSupportDomain)::B                               ! Basis B (Dimension of B is : MLSOrder*k)  
        real(kind=8),dimension(MLSOrder,MLSOrder)::A,A_inversion                                          ! A matrix (Dimension of A and A_inversion is : MLSOrder*MLSOrder)
        real(kind=8),dimension(MLSOrder,PredictionNumberInSupportDomain)::A_inversion_times_B
        real(kind=8),dimension(1,1)::Temp_Fai
        real(kind=8)::sum_w_kernel
                      
        integer::check_error
        !==========================================================================================================

        ! Body of subroutine MLS_Interpolation


        !==========================================================================================================
        ! Weight function value matrix
        w_matrix=0.0d0                                       ! (Dimension of W is : k*k)
        do k=1,Ture_effect_number 
            w_matrix(k,k)=W_kernel(k)
        end do
            
        ! Assign the Modified_matrix
        Modified_matrix=0.0d0
        do k=1,MLSOrder
            Modified_matrix(k,k)=slighted_modified_alpha 
        end do

        !Initialization
        P=0.0d0
        B=0.0d0
        A=0.0d0
        
        check_error=1
            
        !----------------------------------------------------------------------------------------------------------
        ! MLS interpolation
        ! Evaluate basis P, B Matrix and their derivates  (m is the order of MLS)
        do k=1,Ture_effect_number
    
            ! effect particle index
            j=Ture_effect_name(k)
        
            if(MLSOrder==1) then
    
                P(1,k)=1                                       ! P_i=(1)            (Dimension of P is : k*MLSOrder)
    
            else if(MLSOrder==3) then
                                
                P(1,k)=1                                       ! P_i=(1,x_j/h,y_j/h)
                P(2,k)=(subdomain_particle_position(j,1)-Interpolate_subdomain_particle_position(1))/smooth_length
                P(3,k)=(subdomain_particle_position(j,2)-Interpolate_subdomain_particle_position(2))/smooth_length
        
            else if(MLSOrder==6) then
    
                P(1,k)=1                                       ! P_i=(1,x_j/h,y_j/h,(x_j/h)^2,x_j/h*y_j/h,(y_j/h)^2)
                P(2,k)=(subdomain_particle_position(j,1)-Interpolate_subdomain_particle_position(1))/smooth_length
                P(3,k)=(subdomain_particle_position(j,2)-Interpolate_subdomain_particle_position(2))/smooth_length
                P(4,k)=P(2,k)**2
                P(5,k)=P(2,k)*P(3,k)
                P(6,k)=P(3,k)**2
    
            else
                write(*,'(A)') " The order of MLS input is not right , Please ckeck 'MLSOrder' value ! (Interpolation subdomain)" 
            end if
    
        end do
        !----------------------------------------------------------------------------------------------------------
         
        !----------------------------------------------------------------------------------------------------------
        ! Calculate P_i and its derivates
        if(MLSOrder==1) then
        
            P_i(1,1)=1                                       ! P_i=(1)            (Dimension of P is : MLSOrder*1)
        
        else if(MLSOrder==3) then
                                
            P_i(1,1)=1                                       ! P_i=(1,x_i,y_i)
            P_i(2,1)=0.0d0
            P_i(3,1)=0.0d0
        
        else if(MLSOrder==6) then
        
            P_i(1,1)=1                                       ! P_i=(1,x_i,y_i,x_i^2,x_i*y_i,y_i^2)
            P_i(2,1)=0.0d0
            P_i(3,1)=0.0d0
            P_i(4,1)=0.0d0
            P_i(5,1)=0.0d0
            P_i(6,1)=0.0d0
        
        else
            write(*,'(A)') " The order of MLS input is not right , Please ckeck 'MLSOrder' value ! (Interpolation subdomain)" 
        end if
        !----------------------------------------------------------------------------------------------------------
       
        !----------------------------------------------------------------------------------------------------------
        ! Modified the P and W matrix (Makesure non-singal happens)
        if(MLSOrder/=1) then
        
            ! Modified P: New_P=[P Modified_matrix]   (dimension: m*(k+m) )
            ! Modified W: New_W=[ W               0           ]
            !                   [ 0       Modified_matrix     ]                      (dimension: (k+m)*(k+m) )
            do m=1,MLSOrder
                do L=1,MLSOrder
                    k=Ture_effect_number+L
                    P(m,k)=Modified_matrix(m,L)                     ! Comble P and Modified_matrix( as the add matrix is constant so its derivated keep same)
                end do
            end do 
        
            do L=1,MLSOrder
                 k=Ture_effect_number+L
                 W_matrix(k,k)=Modified_matrix(L,L)                 ! Comble P and Modified_matrix( as the add matrix is constant so its derivated keep same)
            end do
        
        end if
        !----------------------------------------------------------------------------------------------------------
    
        !----------------------------------------------------------------------------------------------------------
        ! Calculate B,  A_inversion
        ! Assign the B=P*W                
        B=matmul(P,W_matrix)                          ! (Dimension of B is : MLSOrder*k)
            
        ! Evaluate matrices A=P*W*transpose(P) and its derivatives
        A=matmul(B,transpose(P))                      ! (Dimension of A is : MLSOrder*MLSOrder)
    
        check_error=1
        if(MLSOrder==1) then
            
            A_inversion(1,1)=1.0/A(1,1)
    
        else if(MLSOrder==3) then
            
            call matrix_inversion_less_three(MLSOrder,A,A_inversion,check_error)

            if(check_error==0) then
                write(*,*) "The A_inversion compute is wrong! (Interpolation Subdomain)"
            end if
        
        else if(MLSOrder==6) then 
                    
            A_inversion=A                                          !Initialized the A inversion matrix
            !Calculate the inverse  matrix of A
            !call matrix_inversion_less_six(MLSOrder,A,A_inversion,check_error)
            call BRINV(A_inversion,MLSOrder,check_error)
            
            if(check_error==0) then
                write(*,*) "The A_inversion compute is wrong! (Interpolation Subdomain)"
            end if
    
        end if
        !----------------------------------------------------------------------------------------------------------

        !----------------------------------------------------------------------------------------------------------
        ! Check inversion calculation is right or not
        if(check_error/=0) then

            !******************************************************************************************************
            ! Calculate A_inversion times B  and its derivates (dimension: MLSOrder*k)
            A_inversion_times_B=matmul(A_inversion,B)
     
            Q_i=P_i
            do k=1,Ture_effect_number
    
                ! Calculate P_i and B_i   
                Temp_Fai=0.0d0
                do m=1,MLSOrder
                
                    Temp_Fai=Temp_Fai+Q_i(m,1)*A_inversion_times_B(m,k)
              
                end do
        
                Fai(k)=Temp_Fai(1,1)

            end do
            !******************************************************************************************************

        ! For the case, the inverse oricedure is not right
        elseif(check_error==0) then

            sum_w_kernel=0.0d0
            do k=1,Ture_effect_number
                sum_w_kernel=sum_w_kernel+W_kernel(k)  
            end do
 
            do k=1,Ture_effect_number
                Fai(k)=W_kernel(k)/sum_w_kernel
            end do
        
        end if
        !----------------------------------------------------------------------------------------------------------

        !----------------------------------------------------------------------------------------------------------
        !For high order calculation, when the near particle is not enough, decrease the MLS Order
        if(MLSOrder==6 .and. Ture_effect_number<=8) then

            sum_w_kernel=0.0d0
            do k=1,Ture_effect_number
                sum_w_kernel=sum_w_kernel+W_kernel(k)  
            end do
 
            do k=1,Ture_effect_number
                Fai(k)=W_kernel(k)/sum_w_kernel
            end do
        
        end if
        !----------------------------------------------------------------------------------------------------------

    end subroutine MLS_Interpolation_In_Subdomain_Domain
    !==============================================================================================================







    !**************************************************************************************************************
    !  SUBROUTINE: Compatibility_Interpolation_In_Subdomain_Domain
    !
    !  Purpose   : Compatibility interpolation calculation for particle in subdomain domian 
    !
    !  Programer : Shanqin Jin (E-mail:sjin@mun.ca)
    !
    !  Location  : MUN
    !
    !  Time      : 2017.3.18
    !
    !  Copyright : Memorial University
    !
    !  Version   : 1.0
    !
    !  Note      : MPI version: mpich-3.2
    !              Fortran version: Inter Fortran (ifort) 18.0.0
    !**************************************************************************************************************

    !==============================================================================================================
    subroutine Compatibility_Interpolation_In_Subdomain_Domain(Interpolate_particle_position,Ture_effect_number,Ture_effect_name,W_kernel,DwDx_kernel,Compatibility_Result_Matrix,check_error)

        implicit none

        !==========================================================================================================
        ! Variables
        ! Variables from superior subroutine
        real(kind=8),dimension(dim),intent(in)::Interpolate_particle_position                              ! Interpolation particle position
        integer,intent(in)::Ture_effect_number                                                             ! Ture effect particle number
        integer,dimension(PredictionNumberInSupportDomain),intent(in)::Ture_effect_name                    ! Ture effect particle index
        real(kind=8),dimension(PredictionNumberInSupportDomain),intent(in)::W_kernel                       ! Weight function value
        real(kind=8),dimension(PredictionNumberInSupportDomain,dim),intent(in)::DwDx_kernel                ! Weight function derivates value
        
        real(kind=8),dimension(Compatibility_Equation_dim,4),intent(out)::Compatibility_Result_Matrix      ! Result matrix
        integer,intent(out)::check_error                                                                   ! Interpolation successful---1; Failed---0 

        ! Variables in local subroutine
        integer::i,j,k,L,m                                                                                 ! Variables for looping
        real(kind=8),dimension(Compatibility_Equation_dim,Compatibility_Equation_dim)::Compatibility_Coefficient_Matrix
        real(kind=8),dimension(Compatibility_Equation_dim,4)::Compatibility_RHS_Matrix                     ! Right hand side matrix
        real(kind=8),dimension(Compatibility_Equation_dim,Compatibility_Equation_dim)::A_inversion         ! Inversion matrix of compatibility coefficient matrix   

        !==========================================================================================================

        ! Body of subroutine Compatibility_Interpolation_In_Subdomain_Domain

        !----------------------------------------------------------------------------------------------------------
        ! Compatibility interpolation

        Compatibility_Coefficient_Matrix=0.0d0
        Compatibility_RHS_Matrix=0.0d0

        do k=1,Ture_effect_number
    
            ! Effect particle index
            j=Ture_effect_name(k)

            ! Position difference
            position_difference(:)=Interpolate_particle_position(:)-subdomain_particle_position(j,:)

            !------------------------------------------------------------------------------------------------------
            ! Coefficient matrix for compatibility interpolation
            if ( dim==2 ) then

                Compatibility_Coefficient_Matrix(1,1)=Compatibility_Coefficient_Matrix(1,1)+W_kernel(k)
                Compatibility_Coefficient_Matrix(1,2)=Compatibility_Coefficient_Matrix(1,2)-W_kernel(k)*position_difference(1)
                Compatibility_Coefficient_Matrix(1,3)=Compatibility_Coefficient_Matrix(1,3)-W_kernel(k)*position_difference(2)
                
                Compatibility_Coefficient_Matrix(2,1)=Compatibility_Coefficient_Matrix(2,1)+DwDx_kernel(k,1)
                Compatibility_Coefficient_Matrix(2,2)=Compatibility_Coefficient_Matrix(2,2)-DwDx_kernel(k,1)*position_difference(1)
                Compatibility_Coefficient_Matrix(2,3)=Compatibility_Coefficient_Matrix(2,3)-DwDx_kernel(k,1)*position_difference(2)
                
                Compatibility_Coefficient_Matrix(3,1)=Compatibility_Coefficient_Matrix(3,1)+DwDx_kernel(k,2)
                Compatibility_Coefficient_Matrix(3,2)=Compatibility_Coefficient_Matrix(3,2)-DwDx_kernel(k,2)*position_difference(1)
                Compatibility_Coefficient_Matrix(3,3)=Compatibility_Coefficient_Matrix(3,3)-DwDx_kernel(k,2)*position_difference(2)
            

                !--------------------------------------------------------------------------------------------------
                ! Right hand side for pressure interpolation
                Compatibility_RHS_Matrix(1,1)=Compatibility_RHS_Matrix(1,1)+subdomain_particle_press(j)*W_kernel(k)
                Compatibility_RHS_Matrix(2,1)=Compatibility_RHS_Matrix(2,1)+subdomain_particle_press(j)*DwDx_kernel(k,1)
                Compatibility_RHS_Matrix(3,1)=Compatibility_RHS_Matrix(3,1)+subdomain_particle_press(j)*DwDx_kernel(k,2)

                ! Right hand side for velocity interpolation
                do L=1,dim
                    
                    m=L+1

                    Compatibility_RHS_Matrix(1,m)=Compatibility_RHS_Matrix(1,m)+subdomain_particle_velocity(j,L)*W_kernel(k)
                    Compatibility_RHS_Matrix(2,m)=Compatibility_RHS_Matrix(2,m)+subdomain_particle_velocity(j,L)*DwDx_kernel(k,1)
                    Compatibility_RHS_Matrix(3,m)=Compatibility_RHS_Matrix(3,m)+subdomain_particle_velocity(j,L)*DwDx_kernel(k,2)

                enddo
                !--------------------------------------------------------------------------------------------------


            elseif( dim==3 ) then

                Compatibility_Coefficient_Matrix(1,1)=Compatibility_Coefficient_Matrix(1,1)+W_kernel(k)
                Compatibility_Coefficient_Matrix(1,2)=Compatibility_Coefficient_Matrix(1,2)-W_kernel(k)*position_difference(1)
                Compatibility_Coefficient_Matrix(1,3)=Compatibility_Coefficient_Matrix(1,3)-W_kernel(k)*position_difference(2)
                Compatibility_Coefficient_Matrix(1,4)=Compatibility_Coefficient_Matrix(1,4)-W_kernel(k)*position_difference(3)

                Compatibility_Coefficient_Matrix(2,1)=Compatibility_Coefficient_Matrix(2,1)+DwDx_kernel(k,1)
                Compatibility_Coefficient_Matrix(2,2)=Compatibility_Coefficient_Matrix(2,2)-DwDx_kernel(k,1)*position_difference(1)
                Compatibility_Coefficient_Matrix(2,3)=Compatibility_Coefficient_Matrix(2,3)-DwDx_kernel(k,1)*position_difference(2)
                Compatibility_Coefficient_Matrix(2,4)=Compatibility_Coefficient_Matrix(2,4)-DwDx_kernel(k,1)*position_difference(3)

                Compatibility_Coefficient_Matrix(3,1)=Compatibility_Coefficient_Matrix(3,1)+DwDx_kernel(k,2)
                Compatibility_Coefficient_Matrix(3,2)=Compatibility_Coefficient_Matrix(3,2)-DwDx_kernel(k,2)*position_difference(1)
                Compatibility_Coefficient_Matrix(3,3)=Compatibility_Coefficient_Matrix(3,3)-DwDx_kernel(k,2)*position_difference(2)
                Compatibility_Coefficient_Matrix(3,4)=Compatibility_Coefficient_Matrix(3,4)-DwDx_kernel(k,2)*position_difference(3)

                Compatibility_Coefficient_Matrix(4,1)=Compatibility_Coefficient_Matrix(4,1)+DwDx_kernel(k,3)
                Compatibility_Coefficient_Matrix(4,2)=Compatibility_Coefficient_Matrix(4,2)-DwDx_kernel(k,3)*position_difference(1)
                Compatibility_Coefficient_Matrix(4,3)=Compatibility_Coefficient_Matrix(4,3)-DwDx_kernel(k,3)*position_difference(2)
                Compatibility_Coefficient_Matrix(4,4)=Compatibility_Coefficient_Matrix(4,4)-DwDx_kernel(k,3)*position_difference(3)


                !--------------------------------------------------------------------------------------------------
                ! Right hand side for pressure interpolation
                Compatibility_RHS_Matrix(1,1)=Compatibility_RHS_Matrix(1,1)+subdomain_particle_press(j)*W_kernel(k)
                Compatibility_RHS_Matrix(2,1)=Compatibility_RHS_Matrix(2,1)+subdomain_particle_press(j)*DwDx_kernel(k,1)
                Compatibility_RHS_Matrix(3,1)=Compatibility_RHS_Matrix(3,1)+subdomain_particle_press(j)*DwDx_kernel(k,2)
                Compatibility_RHS_Matrix(4,1)=Compatibility_RHS_Matrix(4,1)+subdomain_particle_press(j)*DwDx_kernel(k,3)

                ! Right hand side for velocity interpolation
                do L=1,dim
                    
                    m=L+1

                    Compatibility_RHS_Matrix(1,m)=Compatibility_RHS_Matrix(1,m)+subdomain_particle_velocity(j,L)*W_kernel(k)
                    Compatibility_RHS_Matrix(2,m)=Compatibility_RHS_Matrix(2,m)+subdomain_particle_velocity(j,L)*DwDx_kernel(k,1)
                    Compatibility_RHS_Matrix(3,m)=Compatibility_RHS_Matrix(3,m)+subdomain_particle_velocity(j,L)*DwDx_kernel(k,2)
                    Compatibility_RHS_Matrix(4,m)=Compatibility_RHS_Matrix(4,m)+subdomain_particle_velocity(j,L)*DwDx_kernel(k,3)

                enddo
                !--------------------------------------------------------------------------------------------------


            endif
            !------------------------------------------------------------------------------------------------------

        end do


        ! Solve the linear equations get the pressure, velocity and their derivates
        check_error=1
        call matrix_inversion_less_six(Compatibility_Equation_dim,Compatibility_Coefficient_Matrix,A_inversion,check_error)

        if (check_error/=0) then
            
            ! First column is the pressure and its derivates
            ! Second column is the Ux and its derivates
            ! Third  column is the Uy and its derivates
            ! Fourth column is the Uz and its derivates (if has)
            Compatibility_Result_Matrix=matmul(A_inversion,Compatibility_RHS_Matrix)
            
        endif
        !----------------------------------------------------------------------------------------------------------
         

    end subroutine Compatibility_Interpolation_In_Subdomain_Domain
    !==============================================================================================================







    !==============================================================================================================
    ! Sampling the pressure at the pressure probes
    Subroutine Pressure_Sampling(SamplingPosition,AveragePressure,MLSPressure,NumericalVelocity)
        
        implicit none

        !==========================================================================================================
        
        !----------------------------------------------------------------------------------------------------------
        ! Variables from superior subroutine
        real(kind=8),dimension(dim),intent(in)::SamplingPosition            ! Current current Sampling Point Position
        real(kind=8),intent(out)::AveragePressure,MLSPressure               ! Average and Normalized Pressure 
        real(kind=8),dimension(dim),intent(out)::NumericalVelocity          ! Sampling Numerical Velocity
        !----------------------------------------------------------------------------------------------------------


        !----------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m,v                                                ! Loop Variables                                             
        
        ! Grid information
        integer::x_n,y_n,mesh_name
        integer::near_mesh_number,near_mesh_name,particle_in_mesh
        integer::Sampling_near_particle_number                              ! Sampling near particle number
        integer,dimension(PredictionNumberInSupportDomain)::Sampling_near_particle_name
             
        ! Weight function variables
        real(kind=8),dimension(PredictionNumberInSupportDomain)::W_kernel   ! Kernel function value
        real(kind=8),dimension(PredictionNumberInSupportDomain)::Fai        ! Weight Value 
        integer,dimension(dim)::Grid_Position_Index                         ! The position index of grid which point in   
        !----------------------------------------------------------------------------------------------------------

        !==========================================================================================================

        ! Body of subroutine Pressure_Sampling

        !----------------------------------------------------------------------------------------------------------
        ! Get the location of the probe based on the background grid (***This is very important***)
        ! The pressure probe grid index should follow the same format of the Gobal background grid calculation
        ! Calculate the grid name and its near grids number

        call Point_Project_Into_Gobal_Background_Grid( SamplingPosition,Grid_Position_Index,mesh_name )

        near_mesh_number=chain_near_mesh_number(mesh_name)         ! Near mesh number
        !----------------------------------------------------------------------------------------------------------
         

        !----------------------------------------------------------------------------------------------------------
        ! Get the effect particle in the near grids
        ! Only for the particle in the subdomain of the Sampling points
        Sampling_near_particle_number=0
        Sampling_near_particle_name=0                                  
         
        do k=1,near_mesh_number                         
             
             near_mesh_name=chain_near_mesh_name(mesh_name,k)      ! Current near grid name
             particle_in_mesh=number(near_mesh_name)               ! Particle number in current grid
           
             do v=1,particle_in_mesh                               ! Search for each particle

                i=particle_in_chain(near_mesh_name,v)               

                if( i<=particle_ture_number+wave_maker_particle_number ) then

                    position_difference(:)=particle_position(i,:)-SamplingPosition(:)
                    distance=sqrt(DOT_PRODUCT(position_difference,position_difference))

                    !make sure the Memory for the effect_particle_number is enough
                    if( distance<=kernel_scale*smooth_length ) then

                        if( Sampling_near_particle_number<PredictionNumberInSupportDomain ) then

                            Sampling_near_particle_number=Sampling_near_particle_number+1
                            Sampling_near_particle_name(Sampling_near_particle_number)=i
                            
                        else
                            
                            write(*,*) "Memory for effect_particle is not engough, it should be larger than 300! (Pressure Sampling)"

                        endif

                    endif

                endif
               
            enddo

         enddo
        !----------------------------------------------------------------------------------------------------------
        

        !----------------------------------------------------------------------------------------------------------
        !initialized the Sampling Variables
        AveragePressure=0.0d0
        MLSPressure=0.0d0
        NumericalVelocity=0.0d0

        !calculate the pressure
        if(Sampling_near_particle_number>3) then
                     
            do L=1,Sampling_near_particle_number

                j=Sampling_near_particle_name(L)
                     
                ! position difference
                position_difference(:)=SamplingPosition(:)-particle_position(j,:)
                distance=sqrt(DOT_PRODUCT(position_difference,position_difference))
            
                call compute_kernel(Kernel_function_model,&                ! Kernel function model(in)   
                                    dim,&                                  ! dimension of the simulation
                                    distance,&                             ! distance(in)                
                                    position_difference,&                  ! position_difference(in)
                                    smooth_length,&                        ! smooth length(in)
                                    temp_w,&                               ! kernel value(out)
                                    temp_dwdx&                             ! derivates of the kernel function(out)
                                    )

                ! kernel function value
                W_kernel(L)=temp_w*particle_mass(j)/particle_rho(j)
                
                ! Sampling_Point_Average_Pressure
                AveragePressure=AveragePressure+particle_press(j)

            end do

            ! Get the weight value of the effect particle
            call MLS_Interpolation_In_Gobal_Domain(MLS_order_Sampling,SamplingPosition,Sampling_near_particle_number,Sampling_near_particle_name,W_kernel,Fai)

            ! Sampling Point MLS Pressure
            do L=1,Sampling_near_particle_number
        
                ! Get the fluid particle index
                j=Sampling_near_particle_name(L)

                ! Samping point - fluid particle
                position_difference(:)=SamplingPosition(:)-particle_position(j,:)

                ! Interpolation
                MLSPressure=MLSPressure+particle_press(j)*Fai(L)+particle_rho(j)*Fai(L)*DOT_PRODUCT(Gravity_grad,position_difference)

                NumericalVelocity(:)=NumericalVelocity(:)+particle_velocity(j,:)*Fai(L)

            enddo
         
            ! Sampling Point Normalized Pressure and velocity
            AveragePressure=AveragePressure/Sampling_near_particle_number

        endif
        !----------------------------------------------------------------------------------------------------------

    end subroutine Pressure_Sampling
    !==============================================================================================================






    !==============================================================================================================
    subroutine Cubic_Spline_Interpolation( Data_row_number,Data_column_number,Import_Data,Coeffcient_Result )

        implicit none
        
        !----------------------------------------------------------------------------------------------------------
        ! Variables from superior subroutine
        integer,intent(in)::Data_row_number                                                         ! Data row number
        integer,intent(in)::Data_column_number                                                      ! Data column number
        real(kind=8),dimension(Data_row_number,Data_column_number),intent(in)::Import_Data          ! Import data
        real(kind=8),dimension(Data_row_number,Data_column_number,4),intent(out)::Coeffcient_Result ! Coeffcient_Result


        ! Variables in local subroutine
        integer::i,j,k,l,v,o,m                                                                      ! Variables for looping

        real(kind=8),allocatable,dimension(:)::h                                                    ! x interval
        real(kind=8),allocatable,dimension(:)::miu                                                 ! miu for coefficient matrix calculation
        real(kind=8),allocatable,dimension(:)::lamada                                               ! lamada for coefficient matrix calculation

        real(kind=8),allocatable,dimension(:,:)::Cubic_Spline_Coefficient_Matrix                    ! Coefficient matrix for cubic spline interpolation
        real(kind=8),allocatable,dimension(:,:)::Coefficient_Inversion_Matrix                       ! Coefficient inversion matrix for cubic spline interpolation
        

        real(kind=8),allocatable,dimension(:,:)::Cubic_Spline_RHS                                   ! Right hand side for cubic spline interpolation
        real(kind=8),allocatable,dimension(:,:)::Cubic_Spline_M_Result                              ! M result for cubic spline interpolation
        
        integer::check_error                                                                        ! Inverse operation return value
        
        !----------------------------------------------------------------------------------------------------------


        ! Body of subroutine Initialziting_Reading_File

        !----------------------------------------------------------------------------------------------------------

        Allocate( h      ( Data_row_number ) )
        Allocate( miu    ( Data_row_number ) )
        Allocate( lamada ( Data_row_number ) )

        Allocate( Cubic_Spline_Coefficient_Matrix( Data_row_number,Data_row_number ) )
        Allocate( Coefficient_Inversion_Matrix   ( Data_row_number,Data_row_number ) )
        
        Allocate( Cubic_Spline_RHS      ( Data_row_number,1 ) )
        Allocate( Cubic_Spline_M_Result ( Data_row_number,1 ) )
        

        Coeffcient_Result=0.0d0

        ! Calculate the interval h=x_(i+1)-x_(i)
        h=0.0d0
        do i=1, Data_row_number-1                                                 ! The calculate interval is i----i+1
            
            h(i)=Import_Data(i+1,1)-Import_Data(i,1)                              ! Distance between i and i+1

        enddo

        ! Calculate the coefficient for cubic spline interpolation
        miu=0.0d0
        lamada=0.0d0
        do i=2, Data_row_number-1                                                
            
            miu(i)    = h(i-1)/(h(i-1)+h(i))
            lamada(i) = h(i)  /(h(i-1)+h(i))
        
        enddo

        ! Combine the coefficient array
        Cubic_Spline_Coefficient_Matrix=0.0d0
        do i=1, Data_row_number  

            if( i==1 ) then
            
                Cubic_Spline_Coefficient_Matrix(i,i)  =2.0d0    ! Diagonal value
                Cubic_Spline_Coefficient_Matrix(i,i+1)=0.0d0    ! For free boundary condition (the second order derivatity is zero)
                
            elseif( i==Data_row_number ) then
                
                Cubic_Spline_Coefficient_Matrix(i,i-1)=0.0d0    ! For free boundary condition (the second order derivatity is zero)
                Cubic_Spline_Coefficient_Matrix(i,i)  =2.0d0    ! Diagonal value

            else
                
                Cubic_Spline_Coefficient_Matrix(i,i-1)=miu(i)
                Cubic_Spline_Coefficient_Matrix(i,i)  =2.0d0    ! Diagonal value
                Cubic_Spline_Coefficient_Matrix(i,i+1)=lamada(i)
                
            endif
            
        enddo


        ! Synchronize all processors calculation 
        call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)


        ! Calculate the inversion martix
        Coefficient_Inversion_Matrix=Cubic_Spline_Coefficient_Matrix
        call BRINV( Coefficient_Inversion_Matrix,Data_row_number,check_error )
        

        ! Synchronize all processors calculation 
        call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)


        ! Calculate the right hand side of function for various y data cubic spline interpolation
        ! As the first column data is the time series, it doesn't have the cubic spline coefficient
        Coeffcient_Result = 0.0d0
        do j=2,Data_column_number

            Cubic_Spline_RHS=0.0
            do i=1, Data_row_number  
                
                if( i==1 ) then
                    Cubic_Spline_RHS(i,1)=0.0d0
                elseif( i==Data_row_number ) then
                    Cubic_Spline_RHS(i,1)=0.0d0
                else
                    Cubic_Spline_RHS(i,1)=6*((Import_Data(i+1,j)-Import_Data(i,j))/(Import_Data(i+1,1)-Import_Data(i,1))-(Import_Data(i,j)-Import_Data(i-1,j))/(Import_Data(i,1)-Import_Data(i-1,1)))/(h(i-1)+h(i))
                endif
                
            enddo

            Cubic_Spline_M_Result = matmul( Coefficient_Inversion_Matrix,Cubic_Spline_RHS )

            ! Calculate the cubic spline interpolation coefficient
            ! S(x)=a*[ x(i+1)-x ]^3 + b*[ x-x(i) ]^3+c*[ x(i+1)-x ]+d*[ x-x(i) ]   x in [x(i),x(i+1) ]

            do i=1,Data_row_number-1

                Coeffcient_Result(i,j,1) = Cubic_Spline_M_Result(i,1)  /( 6*h(i) )   ! a
                Coeffcient_Result(i,j,2) = Cubic_Spline_M_Result(i+1,1)/( 6*h(i) )   ! b
                Coeffcient_Result(i,j,3) = ( Import_Data(i,j)  -Cubic_Spline_M_Result(i,1)  *h(i)**2/6.0 )/h(i)   ! c
                Coeffcient_Result(i,j,4) = ( Import_Data(i+1,j)-Cubic_Spline_M_Result(i+1,1)*h(i)**2/6.0 )/h(i)   ! d

            enddo

        enddo

        !----------------------------------------------------------------------------------------------------------



    end subroutine Cubic_Spline_Interpolation
    !==============================================================================================================







    !==============================================================================================================
    subroutine Output_particle_identification()

        implicit none
        
        !----------------------------------------------------------------------------------------------------------

        ! Variables in local subroutine
        integer::i,j,k,L                                                                    ! Variables for looping

        !----------------------------------------------------------------------------------------------------------


        ! Body of subroutine Output_particle_identification


        !----------------------------------------------------------------------------------------------------------
        ! Calculate the particle number and index for outputting
        Output_particle_number = 0
        Output_particle_index  = 0

        do i=1,ture_total_particle_number

            !------------------------------------------------------------------------------------------------------
            ! You can control the particle type you want output
            
            if ( particle_type(i) /= Water_Particle_Label ) cycle   
             
            if (  OutputDomainXStart <= particle_position(i,1) .and. particle_position(i,1) <= OutputDomainXEnd ) then

                if ( OutputDomainYStart <= particle_position(i,2) .and. particle_position(i,2) <= OutputDomainYEnd ) then

                    Output_particle_number=Output_particle_number+1

                    Output_particle_index(Output_particle_number)=i
                    
                endif

            endif
            !------------------------------------------------------------------------------------------------------
            


            ! !------------------------------------------------------------------------------------------------------
            ! ! For debugging (Output all the particle)

            ! Output_particle_number=Output_particle_number+1

            ! Output_particle_index(Output_particle_number)=i

            ! !------------------------------------------------------------------------------------------------------


        enddo
        !----------------------------------------------------------------------------------------------------------


    end subroutine Output_particle_identification
    !==============================================================================================================

























































    !==========================================================================================================
    !                                       MPI Subroutine (Data Transfer)
    !                                                                _(\_/) 
    !                                                              ,((((^`\
    !                                                             ((((  (6 \ 
    !                                                           ,((((( ,    \
    !                                       ,,,_              ,(((((  /"._  ,`,
    !                                      ((((\\ ,...       ,((((   /    `-.-'
    !                                      )))  ;'    `"'"'""((((   (      
    !                                     (((  /            (((      \
    !                                      )) |                      |
    !                                     ((  |        .       '     |
    !                                     ))  \     _ '      `t   ,.')
    !                                     (   |   y;- -,-""'"-.\   \/  
    !                                     )   / ./  ) /         `\  \
    !                                        |./   ( (           / /'
    !                                        ||     \\          //'|
    !                                        ||      \\       _//'||
    !                                        ||       ))     |_/  ||
    !                                        \_\     |_/          ||
    !                                        `'"                  \_\
    !                                                             `'"
    !==========================================================================================================


    !==========================================================================================================
    !MPI Message exchange subroutines

    !==========================================================================================================
    subroutine MPI_RealScalar_MessageExchange(ExchangeData)

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                           !Variables for looping
        real(kind=8),dimension(subdomain_ntotal)::ExchangeData       !Input Exchange Data
        !------------------------------------------------------------------------------------------------------

        ! Body of MPI_RealScalar_MessageExchange


        !------------------------------------------------------------------------------------------------------
        if (MPI_Dim==1) then

            !--------------------------------------------------------------------------------------------------
            !particle_press
            left_Message_ID=1
            right_Message_ID=2

            !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
            !Data form left processor to right processor
            call MPI_SENDRECV(ExchangeData(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                              ExchangeData(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Data form right processor to left processor
            call MPI_SENDRECV(ExchangeData(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                              ExchangeData(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

        elseif (MPI_Dim==2) then 


            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !Send and receive the four corners

            !--------------------------------------------------------------------------------------------------
            left_Message_ID=1
            right_Message_ID=2

            !Front Corner Data from left processor to right processor
            call MPI_SENDRECV(ExchangeData(RFC_message_FPI),RFC_message_length,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                              ExchangeData(Buffer_LFC_message_FPI_from_left),Buffer_LFC_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Back Corner Data from left processor to right processor
            call MPI_SENDRECV(ExchangeData(RBC_message_FPI),RBC_message_length,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                              ExchangeData(Buffer_LBC_message_FPI_from_left),Buffer_LBC_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            

            !Front Corner Data from right processor to left processor
            call MPI_SENDRECV(ExchangeData(LFC_message_FPI),LFC_message_length,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                              ExchangeData(Buffer_RFC_message_FPI_from_right),Buffer_RFC_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !Back Corner Data from right processor to left processor
            call MPI_SENDRECV(ExchangeData(LBC_message_FPI),LBC_message_length,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                              ExchangeData(Buffer_RBC_message_FPI_from_right),Buffer_RBC_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            Front_Message_ID=3
            Back_Message_ID=4 

            !Left Corner Data from front processor to back processor
            call MPI_SENDRECV(ExchangeData(LBC_message_FPI),LBC_message_length,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                              ExchangeData(Buffer_LFC_message_FPI_from_front),Buffer_LFC_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Right Corner Data from front processor to back processor
            call MPI_SENDRECV(ExchangeData(RBC_message_FPI),RBC_message_length,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                              ExchangeData(Buffer_RFC_message_FPI_from_front),Buffer_RFC_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !Left Corner Data from back processor to front processor
            call MPI_SENDRECV(ExchangeData(LFC_message_FPI),LFC_message_length,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                              ExchangeData(Buffer_LBC_message_FPI_from_back),Buffer_LBC_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !Right Corner Data from back processor to front processor
            call MPI_SENDRECV(ExchangeData(RFC_message_FPI),RFC_message_length,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                              ExchangeData(Buffer_RBC_message_FPI_from_back),Buffer_RBC_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            LeftFront_Message_ID=5
            RightBack_Message_ID=6

            !Corner Data from LeftFront processor to RightBack processor (LeftFront Message)
            call MPI_SENDRECV(ExchangeData(RBC_message_FPI),RBC_message_length,MPI_DOUBLE_PRECISION,RightBack_processor_ID,LeftFront_Message_ID,&
                              ExchangeData(Buffer_LFC_message_FPI_from_LeftFront),Buffer_LFC_message_length_from_LeftFront,MPI_DOUBLE_PRECISION,LeftFront_processor_ID,LeftFront_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Corner Data from RightBack processor to LeftFront processor (RightBack Message)
            call MPI_SENDRECV(ExchangeData(LFC_message_FPI),LFC_message_length,MPI_DOUBLE_PRECISION,LeftFront_processor_ID,RightBack_Message_ID,&
                              ExchangeData(Buffer_RBC_message_FPI_from_RightBack),Buffer_RBC_message_length_from_RightBack,MPI_DOUBLE_PRECISION,RightBack_processor_ID,RightBack_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------

            LeftBack_Message_ID=7
            RightFront_Message_ID=8

            !Corner Data from LeftBack processor to RightFront processor (LeftBack Message)
            call MPI_SENDRECV(ExchangeData(RFC_message_FPI),RFC_message_length,MPI_DOUBLE_PRECISION,RightFront_processor_ID,LeftBack_Message_ID,&
                              ExchangeData(Buffer_LBC_message_FPI_from_LeftBack),Buffer_LBC_message_length_from_LeftBack,MPI_DOUBLE_PRECISION,LeftBack_processor_ID,LeftBack_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Corner Data from RightFront processor to LeftBack processor (RightFront Message)
            call MPI_SENDRECV(ExchangeData(LBC_message_FPI),LBC_message_length,MPI_DOUBLE_PRECISION,LeftBack_processor_ID,RightFront_Message_ID,&
                              ExchangeData(Buffer_RFC_message_FPI_from_RightFront),Buffer_RFC_message_length_from_RightFront,MPI_DOUBLE_PRECISION,RightFront_processor_ID,RightFront_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !Send and receive the four sides message

            !--------------------------------------------------------------------------------------------------
            left_Message_ID=9
            right_Message_ID=10

            !Side Data form left processor to right processor
            call MPI_SENDRECV(ExchangeData(side_message_FPI_to_right),side_message_length_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                              ExchangeData(side_message_FPI_from_left),side_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Side Data from right processor to left processor
            call MPI_SENDRECV(ExchangeData(side_message_FPI_to_left),side_message_length_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                              ExchangeData(side_message_FPI_from_right),side_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            Front_Message_ID=11
            Back_Message_ID=12

            !Side Data from front processor to back processor
            call MPI_SENDRECV(ExchangeData(side_message_FPI_to_back),side_message_length_to_back,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                              ExchangeData(side_message_FPI_from_front),side_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Side Data from back processor to front processor
            call MPI_SENDRECV(ExchangeData(side_message_FPI_to_front),side_message_length_to_front,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                              ExchangeData(side_message_FPI_from_back),side_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          
        endif

    end subroutine MPI_RealScalar_MessageExchange
    !==========================================================================================================



    !==========================================================================================================
    subroutine MPI_IntScalar_MessageExchange(ExchangeData)

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in program
        integer::i,j,k,L,m                                           !Variables for looping
        integer,dimension(subdomain_ntotal)::ExchangeData            !Input Exchange Data
        !------------------------------------------------------------------------------------------------------

        ! Body of MPI_IntScalar_MessageExchange



        !------------------------------------------------------------------------------------------------------
        if (MPI_Dim==1) then

            !--------------------------------------------------------------------------------------------------
            !particle_press
            left_Message_ID=1
            right_Message_ID=2

            !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
            !Data form left processor to right processor
            call MPI_SENDRECV(ExchangeData(message_first_particle_index_to_right),message_length_send_to_right,MPI_INTEGER,right_processor_ID,left_Message_ID,&
                              ExchangeData(message_first_particle_index_from_left),message_length_from_left,MPI_INTEGER,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Data form right processor to left processor
            call MPI_SENDRECV(ExchangeData(message_first_particle_index_to_left),message_length_send_to_left,MPI_INTEGER,left_processor_ID,right_Message_ID,&
                              ExchangeData(message_first_particle_index_from_right),message_length_from_right,MPI_INTEGER,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

        elseif (MPI_Dim==2) then

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !Send and receive the four corners

            !--------------------------------------------------------------------------------------------------
            left_Message_ID=1
            right_Message_ID=2

            !Front Corner Data from left processor to right processor
            call MPI_SENDRECV(ExchangeData(RFC_message_FPI),RFC_message_length,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                              ExchangeData(Buffer_LFC_message_FPI_from_left),Buffer_LFC_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Back Corner Data from left processor to right processor
            call MPI_SENDRECV(ExchangeData(RBC_message_FPI),RBC_message_length,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                              ExchangeData(Buffer_LBC_message_FPI_from_left),Buffer_LBC_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            

            !Front Corner Data from right processor to left processor
            call MPI_SENDRECV(ExchangeData(LFC_message_FPI),LFC_message_length,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                              ExchangeData(Buffer_RFC_message_FPI_from_right),Buffer_RFC_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !Back Corner Data from right processor to left processor
            call MPI_SENDRECV(ExchangeData(LBC_message_FPI),LBC_message_length,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                              ExchangeData(Buffer_RBC_message_FPI_from_right),Buffer_RBC_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            Front_Message_ID=3
            Back_Message_ID=4 

            !Left Corner Data from front processor to back processor
            call MPI_SENDRECV(ExchangeData(LBC_message_FPI),LBC_message_length,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                              ExchangeData(Buffer_LFC_message_FPI_from_front),Buffer_LFC_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Right Corner Data from front processor to back processor
            call MPI_SENDRECV(ExchangeData(RBC_message_FPI),RBC_message_length,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                              ExchangeData(Buffer_RFC_message_FPI_from_front),Buffer_RFC_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !Left Corner Data from back processor to front processor
            call MPI_SENDRECV(ExchangeData(LFC_message_FPI),LFC_message_length,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                              ExchangeData(Buffer_LBC_message_FPI_from_back),Buffer_LBC_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !Right Corner Data from back processor to front processor
            call MPI_SENDRECV(ExchangeData(RFC_message_FPI),RFC_message_length,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                              ExchangeData(Buffer_RBC_message_FPI_from_back),Buffer_RBC_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            LeftFront_Message_ID=5
            RightBack_Message_ID=6

            !Corner Data from LeftFront processor to RightBack processor (LeftFront Message)
            call MPI_SENDRECV(ExchangeData(RBC_message_FPI),RBC_message_length,MPI_DOUBLE_PRECISION,RightBack_processor_ID,LeftFront_Message_ID,&
                              ExchangeData(Buffer_LFC_message_FPI_from_LeftFront),Buffer_LFC_message_length_from_LeftFront,MPI_DOUBLE_PRECISION,LeftFront_processor_ID,LeftFront_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Corner Data from RightBack processor to LeftFront processor (RightBack Message)
            call MPI_SENDRECV(ExchangeData(LFC_message_FPI),LFC_message_length,MPI_DOUBLE_PRECISION,LeftFront_processor_ID,RightBack_Message_ID,&
                              ExchangeData(Buffer_RBC_message_FPI_from_RightBack),Buffer_RBC_message_length_from_RightBack,MPI_DOUBLE_PRECISION,RightBack_processor_ID,RightBack_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------

            LeftBack_Message_ID=7
            RightFront_Message_ID=8

            !Corner Data from LeftBack processor to RightFront processor (LeftBack Message)
            call MPI_SENDRECV(ExchangeData(RFC_message_FPI),RFC_message_length,MPI_DOUBLE_PRECISION,RightFront_processor_ID,LeftBack_Message_ID,&
                              ExchangeData(Buffer_LBC_message_FPI_from_LeftBack),Buffer_LBC_message_length_from_LeftBack,MPI_DOUBLE_PRECISION,LeftBack_processor_ID,LeftBack_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Corner Data from RightFront processor to LeftBack processor (RightFront Message)
            call MPI_SENDRECV(ExchangeData(LBC_message_FPI),LBC_message_length,MPI_DOUBLE_PRECISION,LeftBack_processor_ID,RightFront_Message_ID,&
                              ExchangeData(Buffer_RFC_message_FPI_from_RightFront),Buffer_RFC_message_length_from_RightFront,MPI_DOUBLE_PRECISION,RightFront_processor_ID,RightFront_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !Send and receive the four sides message

            !--------------------------------------------------------------------------------------------------
            left_Message_ID=9
            right_Message_ID=10

            !Side Data form left processor to right processor
            call MPI_SENDRECV(ExchangeData(side_message_FPI_to_right),side_message_length_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                              ExchangeData(side_message_FPI_from_left),side_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Side Data from right processor to left processor
            call MPI_SENDRECV(ExchangeData(side_message_FPI_to_left),side_message_length_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                              ExchangeData(side_message_FPI_from_right),side_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            Front_Message_ID=11
            Back_Message_ID=12

            !Side Data from front processor to back processor
            call MPI_SENDRECV(ExchangeData(side_message_FPI_to_back),side_message_length_to_back,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                              ExchangeData(side_message_FPI_from_front),side_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
            
            !Side Data from back processor to front processor
            call MPI_SENDRECV(ExchangeData(side_message_FPI_to_front),side_message_length_to_front,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                              ExchangeData(side_message_FPI_from_back),side_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            !--------------------------------------------------------------------------------------------------

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        endif

    end subroutine MPI_IntScalar_MessageExchange
    !==========================================================================================================




    !==========================================================================================================
    subroutine MPI_RealVector_MessageExchange(ExchangeData)

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in program
        integer::i,j,k,L,m                                               !Variables for looping
        real(kind=8),dimension(subdomain_ntotal,dim)::ExchangeData       !Input Exchange Data
        !------------------------------------------------------------------------------------------------------

        ! Body of MPI_RealVector_MessageExchange


        !------------------------------------------------------------------------------------------------------

        if ( MPI_Dim==1 ) then

            !--------------------------------------------------------------------------------------------------
            ! particle_velocity
            do j=1,dim

                left_Message_ID=1
                right_Message_ID=2

                !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
                !Data form left processor to right processor
                call MPI_SENDRECV(ExchangeData(message_first_particle_index_to_right,j),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                                  ExchangeData(message_first_particle_index_from_left,j),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                
                !Data form right processor to left processor
                call MPI_SENDRECV(ExchangeData(message_first_particle_index_to_left,j),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                                  ExchangeData(message_first_particle_index_from_right,j),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

            end do
            !--------------------------------------------------------------------------------------------------

        elseif ( MPI_Dim==2 ) then

            do j=1,dim

                !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                !Send and receive the four corners
                !----------------------------------------------------------------------------------------------
                left_Message_ID=1
                right_Message_ID=2

                !Front Corner Data from left processor to right processor
                call MPI_SENDRECV(ExchangeData(RFC_message_FPI,j),RFC_message_length,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                                  ExchangeData(Buffer_LFC_message_FPI_from_left,j),Buffer_LFC_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                
                !Back Corner Data from left processor to right processor
                call MPI_SENDRECV(ExchangeData(RBC_message_FPI,j),RBC_message_length,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                                  ExchangeData(Buffer_LBC_message_FPI_from_left,j),Buffer_LBC_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                

                !Front Corner Data from right processor to left processor
                call MPI_SENDRECV(ExchangeData(LFC_message_FPI,j),LFC_message_length,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                                  ExchangeData(Buffer_RFC_message_FPI_from_right,j),Buffer_RFC_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                !Back Corner Data from right processor to left processor
                call MPI_SENDRECV(ExchangeData(LBC_message_FPI,j),LBC_message_length,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                                  ExchangeData(Buffer_RBC_message_FPI_from_right,j),Buffer_RBC_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                !----------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------
                Front_Message_ID=3
                Back_Message_ID=4 

                !Left Corner Data from front processor to back processor (Front message)
                call MPI_SENDRECV(ExchangeData(LBC_message_FPI,j),LBC_message_length,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                                  ExchangeData(Buffer_LFC_message_FPI_from_front,j),Buffer_LFC_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                
                !Right Corner Data from front processor to back processor (Front message)
                call MPI_SENDRECV(ExchangeData(RBC_message_FPI,j),RBC_message_length,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                                  ExchangeData(Buffer_RFC_message_FPI_from_front,j),Buffer_RFC_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                !Left Corner Data from back processor to front processor (Back message)
                call MPI_SENDRECV(ExchangeData(LFC_message_FPI,j),LFC_message_length,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                                  ExchangeData(Buffer_LBC_message_FPI_from_back,j),Buffer_LBC_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                !Right Corner Data from back processor to front processor (Back message)
                call MPI_SENDRECV(ExchangeData(RFC_message_FPI,j),RFC_message_length,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                                  ExchangeData(Buffer_RBC_message_FPI_from_back,j),Buffer_RBC_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                !----------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------
                LeftFront_Message_ID=5
                RightBack_Message_ID=6

                !Corner Data from LeftFront processor to RightBack processor (LeftFront Message)
                call MPI_SENDRECV(ExchangeData(RBC_message_FPI,j),RBC_message_length,MPI_DOUBLE_PRECISION,RightBack_processor_ID,LeftFront_Message_ID,&
                                  ExchangeData(Buffer_LFC_message_FPI_from_LeftFront,j),Buffer_LFC_message_length_from_LeftFront,MPI_DOUBLE_PRECISION,LeftFront_processor_ID,LeftFront_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                
                !Corner Data from RightBack processor to LeftFront processor (RightBack Message)
                call MPI_SENDRECV(ExchangeData(LFC_message_FPI,j),LFC_message_length,MPI_DOUBLE_PRECISION,LeftFront_processor_ID,RightBack_Message_ID,&
                                  ExchangeData(Buffer_RBC_message_FPI_from_RightBack,j),Buffer_RBC_message_length_from_RightBack,MPI_DOUBLE_PRECISION,RightBack_processor_ID,RightBack_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                !----------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------
                LeftBack_Message_ID=7
                RightFront_Message_ID=8

                !Corner Data from LeftBack processor to RightFront processor (LeftBack Message)
                call MPI_SENDRECV(ExchangeData(RFC_message_FPI,j),RFC_message_length,MPI_DOUBLE_PRECISION,RightFront_processor_ID,LeftBack_Message_ID,&
                                  ExchangeData(Buffer_LBC_message_FPI_from_LeftBack,j),Buffer_LBC_message_length_from_LeftBack,MPI_DOUBLE_PRECISION,LeftBack_processor_ID,LeftBack_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                
                !Corner Data from RightFront processor to LeftBack processor (RightFront Message)
                call MPI_SENDRECV(ExchangeData(LBC_message_FPI,j),LBC_message_length,MPI_DOUBLE_PRECISION,LeftBack_processor_ID,RightFront_Message_ID,&
                                  ExchangeData(Buffer_RFC_message_FPI_from_RightFront,j),Buffer_RFC_message_length_from_RightFront,MPI_DOUBLE_PRECISION,RightFront_processor_ID,RightFront_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                !----------------------------------------------------------------------------------------------

                !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


                !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                !Send and receive the four sides message

                !----------------------------------------------------------------------------------------------
                left_Message_ID=9
                right_Message_ID=10

                !Side Data form left processor to right processor
                call MPI_SENDRECV(ExchangeData(side_message_FPI_to_right,j),side_message_length_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                                  ExchangeData(side_message_FPI_from_left,j),side_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                
                !Side Data from right processor to left processor
                call MPI_SENDRECV(ExchangeData(side_message_FPI_to_left,j),side_message_length_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                                  ExchangeData(side_message_FPI_from_right,j),side_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                !----------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------
                Front_Message_ID=11
                Back_Message_ID=12

                !Side Data from front processor to back processor
                call MPI_SENDRECV(ExchangeData(side_message_FPI_to_back,j),side_message_length_to_back,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                                  ExchangeData(side_message_FPI_from_front,j),side_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                
                !Side Data from back processor to front processor
                call MPI_SENDRECV(ExchangeData(side_message_FPI_to_front,j),side_message_length_to_front,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                                  ExchangeData(side_message_FPI_from_back,j),side_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                !----------------------------------------------------------------------------------------------

                !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            enddo
          
        endif

    end subroutine MPI_RealVector_MessageExchange
    !==========================================================================================================


    !==========================================================================================================
    subroutine MPI_RealTensor_MessageExchange(ExchangeData)

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in program
        integer::i,j,k,L,m                                               !Variables for looping
        real(kind=8),dimension(subdomain_ntotal,dim,dim)::ExchangeData   !Input Exchange Data
        !------------------------------------------------------------------------------------------------------

        ! Body of MPI_RealTensor_MessageExchange


        !------------------------------------------------------------------------------------------------------

        if ( MPI_Dim==1 ) then

            !--------------------------------------------------------------------------------------------------
            !Exchange the Tensor_dyadic_inverse_matrix form the neighbour processor
            do j=1,dim
               do k=1,dim
                
                    left_Message_ID=1
                    right_Message_ID=2

                    !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
                    !Data form left processor to right processor
                    call MPI_SENDRECV(ExchangeData(message_first_particle_index_to_right,k,j),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                                      ExchangeData(message_first_particle_index_from_left,k,j),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                    
                    !Data form right processor to left processor
                    call MPI_SENDRECV(ExchangeData(message_first_particle_index_to_left,k,j),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                                      ExchangeData(message_first_particle_index_from_right,k,j),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                    
                enddo
            enddo
            !--------------------------------------------------------------------------------------------------

        elseif ( MPI_Dim==2 ) then
            
            do j=1,dim
               do k=1,dim

                    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    !Send and receive the four corners
                    !------------------------------------------------------------------------------------------
                    left_Message_ID=1
                    right_Message_ID=2

                    !Front Corner Data from left processor to right processor
                    call MPI_SENDRECV(ExchangeData(RFC_message_FPI,k,j),RFC_message_length,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                                      ExchangeData(Buffer_LFC_message_FPI_from_left,k,j),Buffer_LFC_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                    
                    !Back Corner Data from left processor to right processor
                    call MPI_SENDRECV(ExchangeData(RBC_message_FPI,k,j),RBC_message_length,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                                      ExchangeData(Buffer_LBC_message_FPI_from_left,k,j),Buffer_LBC_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                    

                    !Front Corner Data from right processor to left processor
                    call MPI_SENDRECV(ExchangeData(LFC_message_FPI,k,j),LFC_message_length,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                                      ExchangeData(Buffer_RFC_message_FPI_from_right,k,j),Buffer_RFC_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                    !Back Corner Data from right processor to left processor
                    call MPI_SENDRECV(ExchangeData(LBC_message_FPI,k,j),LBC_message_length,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                                      ExchangeData(Buffer_RBC_message_FPI_from_right,k,j),Buffer_RBC_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                    !------------------------------------------------------------------------------------------

                    !------------------------------------------------------------------------------------------
                    Front_Message_ID=3
                    Back_Message_ID=4 

                    !Left Corner Data from front processor to back processor
                    call MPI_SENDRECV(ExchangeData(LBC_message_FPI,k,j),LBC_message_length,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                                      ExchangeData(Buffer_LFC_message_FPI_from_front,k,j),Buffer_LFC_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                    
                    !Right Corner Data from front processor to back processor
                    call MPI_SENDRECV(ExchangeData(RBC_message_FPI,k,j),RBC_message_length,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                                      ExchangeData(Buffer_RFC_message_FPI_from_front,k,j),Buffer_RFC_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                    !Left Corner Data from back processor to front processor
                    call MPI_SENDRECV(ExchangeData(LFC_message_FPI,k,j),LFC_message_length,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                                      ExchangeData(Buffer_LBC_message_FPI_from_back,k,j),Buffer_LBC_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                    !Right Corner Data from back processor to front processor
                    call MPI_SENDRECV(ExchangeData(RFC_message_FPI,k,j),RFC_message_length,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                                      ExchangeData(Buffer_RBC_message_FPI_from_back,k,j),Buffer_RBC_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                    !------------------------------------------------------------------------------------------

                    !------------------------------------------------------------------------------------------
                    LeftFront_Message_ID=5
                    RightBack_Message_ID=6

                    !Corner Data from LeftFront processor to RightBack processor (LeftFront Message)
                    call MPI_SENDRECV(ExchangeData(RBC_message_FPI,k,j),RBC_message_length,MPI_DOUBLE_PRECISION,RightBack_processor_ID,LeftFront_Message_ID,&
                                      ExchangeData(Buffer_LFC_message_FPI_from_LeftFront,k,j),Buffer_LFC_message_length_from_LeftFront,MPI_DOUBLE_PRECISION,LeftFront_processor_ID,LeftFront_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                    
                    !Corner Data from RightBack processor to LeftFront processor (RightBack Message)
                    call MPI_SENDRECV(ExchangeData(LFC_message_FPI,k,j),LFC_message_length,MPI_DOUBLE_PRECISION,LeftFront_processor_ID,RightBack_Message_ID,&
                                      ExchangeData(Buffer_RBC_message_FPI_from_RightBack,k,j),Buffer_RBC_message_length_from_RightBack,MPI_DOUBLE_PRECISION,RightBack_processor_ID,RightBack_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                    !------------------------------------------------------------------------------------------

                    !------------------------------------------------------------------------------------------
                    LeftBack_Message_ID=7
                    RightFront_Message_ID=8

                    !Corner Data from LeftBack processor to RightFront processor (LeftBack Message)
                    call MPI_SENDRECV(ExchangeData(RFC_message_FPI,k,j),RFC_message_length,MPI_DOUBLE_PRECISION,RightFront_processor_ID,LeftBack_Message_ID,&
                                      ExchangeData(Buffer_LBC_message_FPI_from_LeftBack,k,j),Buffer_LBC_message_length_from_LeftBack,MPI_DOUBLE_PRECISION,LeftBack_processor_ID,LeftBack_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                    
                    !Corner Data from RightFront processor to LeftBack processor (RightFront Message)
                    call MPI_SENDRECV(ExchangeData(LBC_message_FPI,k,j),LBC_message_length,MPI_DOUBLE_PRECISION,LeftBack_processor_ID,RightFront_Message_ID,&
                                      ExchangeData(Buffer_RFC_message_FPI_from_RightFront,k,j),Buffer_RFC_message_length_from_RightFront,MPI_DOUBLE_PRECISION,RightFront_processor_ID,RightFront_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                    !------------------------------------------------------------------------------------------

                    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


                    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    !Send and receive the four sides message

                    !------------------------------------------------------------------------------------------
                    left_Message_ID=9
                    right_Message_ID=10

                    !Side Data form left processor to right processor
                    call MPI_SENDRECV(ExchangeData(side_message_FPI_to_right,k,j),side_message_length_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                                      ExchangeData(side_message_FPI_from_left,k,j),side_message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                    
                    !Side Data from right processor to left processor
                    call MPI_SENDRECV(ExchangeData(side_message_FPI_to_left,k,j),side_message_length_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                                      ExchangeData(side_message_FPI_from_right,k,j),side_message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                    !------------------------------------------------------------------------------------------

                    !------------------------------------------------------------------------------------------
                    Front_Message_ID=11
                    Back_Message_ID=12

                    !Side Data from front processor to back processor
                    call MPI_SENDRECV(ExchangeData(side_message_FPI_to_back,k,j),side_message_length_to_back,MPI_DOUBLE_PRECISION,back_processor_ID,Front_Message_ID,&
                                      ExchangeData(side_message_FPI_from_front,k,j),side_message_length_from_front,MPI_DOUBLE_PRECISION,front_processor_ID,Front_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
                    
                    !Side Data from back processor to front processor
                    call MPI_SENDRECV(ExchangeData(side_message_FPI_to_front,k,j),side_message_length_to_front,MPI_DOUBLE_PRECISION,front_processor_ID,Back_Message_ID,&
                                      ExchangeData(side_message_FPI_from_back,k,j),side_message_length_from_back,MPI_DOUBLE_PRECISION,back_processor_ID,Back_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

                    !------------------------------------------------------------------------------------------

                    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

                enddo
            enddo
          
        endif

    end subroutine MPI_RealTensor_MessageExchange
    !==========================================================================================================





    !==========================================================================================================
    !subroutine for Reduce And Bcast
    subroutine MPI_Real_Message_Refresh(SubdomainReduceBcastData,ReduceBcastData)

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in program
        integer::i,j,k,L,m                                                  !Variables for looping
        real(kind=8),dimension(subdomain_ntotal)::SubdomainReduceBcastData  !Subdomain Reduce BcastData
        real(kind=8),dimension(n_total)::ReduceBcastData                    !Input Reduce and Bcast Data 
        !------------------------------------------------------------------------------------------------------

        ! Body of MPI_Real_Message_Refresh

        !------------------------------------------------------------------------------------------------------
        ReduceBcastData=0.0d0

        Temp_Reduce_Real_Array=0.0d0

        do i=1,actual_particle_number_in_subdomain                      
          
           Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=SubdomainReduceBcastData(i)                            

        end do


        call MPI_REDUCE( Temp_Reduce_Real_Array,ReduceBcastData,n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
        call MPI_BCAST( ReduceBcastData,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

        !------------------------------------------------------------------------------------------------------

    end subroutine MPI_Real_Message_Refresh
    !==========================================================================================================


    !==========================================================================================================
    subroutine MPI_Int_Message_Refresh(SubdomainReduceBcastData,ReduceBcastData)

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in program
        integer::i,j,k,L,m                                                  !Variables for looping
        integer,dimension(subdomain_ntotal)::SubdomainReduceBcastData       !Subdomain Reduce BcastData
        integer,dimension(n_total)::ReduceBcastData                         !Input Reduce and Bcast Data 
        !------------------------------------------------------------------------------------------------------

        ! Body of MPI_Int_Message_Refresh


        !------------------------------------------------------------------------------------------------------
        ReduceBcastData=0

        Temp_Reduce_Int_Array=0

        do i=1,actual_particle_number_in_subdomain                      
          
           Temp_Reduce_Int_Array(actual_particle_name_in_subdomain(i))=SubdomainReduceBcastData(i)                            

        end do


        call MPI_REDUCE( Temp_Reduce_Int_Array,ReduceBcastData,n_total,MPI_INTEGER,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
        call MPI_BCAST( ReduceBcastData,n_total,MPI_INTEGER,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

        !------------------------------------------------------------------------------------------------------


    end subroutine MPI_Int_Message_Refresh
    !==========================================================================================================


    !==========================================================================================================
    subroutine MPI_RealVector_Message_Refresh(SubdomainReduceBcastData,ReduceBcastData)

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in program
        integer::i,j,k,L,m                                                       !Variables for looping
        real(kind=8),dimension(subdomain_ntotal,dim)::SubdomainReduceBcastData   !Subdomain Reduce BcastData
        real(kind=8),dimension(n_total,dim)::ReduceBcastData                     !Input Reduce and Bcast Data 
        !------------------------------------------------------------------------------------------------------

        ! Body of MPI_RealVector_Message_Refresh


        !------------------------------------------------------------------------------------------------------
        ReduceBcastData=0.0d0
        
        do j=1,dim

          Temp_Reduce_Real_Array=0.0d0

          do i=1,actual_particle_number_in_subdomain                      
            
             Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=SubdomainReduceBcastData(i,j)                            

          end do


          call MPI_REDUCE( Temp_Reduce_Real_Array,ReduceBcastData(:,j),n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
          call MPI_BCAST( ReduceBcastData(:,j),n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

        end do 
        !------------------------------------------------------------------------------------------------------


    end subroutine MPI_RealVector_Message_Refresh
    !==========================================================================================================



    !==========================================================================================================
    subroutine MPI_RealTensor_Message_Refresh(SubdomainReduceBcastData,ReduceBcastData)

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in program
        integer::i,j,k,L,m                                                                    !Variables for looping
        real(kind=8),dimension(subdomain_ntotal,dim,dim)::SubdomainReduceBcastData            !Input Reduce and Bcast Data 
        real(kind=8),dimension(n_total,dim,dim)::ReduceBcastData                              !Input Reduce and Bcast Data 
        !------------------------------------------------------------------------------------------------------

        ! Body of MPI_RealTensor_Message_Refresh

        !------------------------------------------------------------------------------------------------------
        ReduceBcastData=0.0d0

        do k=1,dim
            do j=1,dim

              Temp_Reduce_Real_Array=0.0d0

              do i=1,actual_particle_number_in_subdomain                      
                
                 Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=SubdomainReduceBcastData(i,k,j)                            

              end do


              call MPI_REDUCE( Temp_Reduce_Real_Array,ReduceBcastData(:,k,j),n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
              call MPI_BCAST( ReduceBcastData(:,k,j),n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

            end do
        end do
        !------------------------------------------------------------------------------------------------------

    end subroutine MPI_RealTensor_Message_Refresh
    !==========================================================================================================






end module SPH_Subroutine_Module
