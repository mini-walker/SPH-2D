!**************************************************************************************************************
!
!  PROGRAM   : Parallel_Moving_Square_Tradictionary_Ghost
!
!  PURPOSE   : (1) SPHERIC benchmark test six--- Moving square inside a rectangular box tradictionary (Parallel - Release Version)
!              (2) Use tradictionary ghost boundary treatment for the reference value
!
!  Programer : Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location  : MUN
!
!  Time      : 2018.12.30
!
!  Copyright : Memorial University
!
!  Version   : 1.1
!
!  Note      : MPI version: mpich-3.2
!              Fortran version: Inter Fortran (ifort) 18.0.0
!            (1) #SBATCH --nodes=2
!            (2) #SBATCH --ntasks-per-node=8        --- Task number is the MPI CPU number   2*8
!            (3) #SBATCH --cpus-per-task=1          --- Openmp threads number for each MPI  1
!            (4) #SBATCH --mem-per-cpu=2048M        --- Memory for each threads             2G per thred
!            (5) ntasks-per-node * cpus-per-task should be the cpus number of one node
!**************************************************************************************************************


Program Parallel_Moving_Square_Tradictionary_Ghost

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables in program
    integer::i,j,k,L                                                        ! Variables for looping

    !Variables for file opeartion
    integer::File_index                                                     ! File index
    character(len=100)::File_name                                           ! File name ( the file name character lenghth should be 100 )
    !==========================================================================================================

    ! Body of Parallel_Moving_Square_Tradictionary_Ghost

    !**********************************************************************************************************
    ! Initialize the MPI runing system
    call MPI_INIT( ierror_MPI )
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    call MPI_GET_PROCESSOR_NAME(Current_Processor_Group_Name,Name_Length,ierror_MPI)

    !----------------------------------------------------------------------------------------------------------
    Program_Beginning: if (ierror_MPI/=0) then

        write(*,*) "There are some errors in MPI Initialization!"

        goto 666                                                               ! Go to the end.

    else

        !******************************************************************************************************
        ! Generate the header of the program
        if ( Current_Processor_ID==Main_Processor ) then

            !--------------------------------------------------------------------------------------------------
            ! Initializing the log file to save the output information 
            File_name="Log_file.dat"
            call Initialziting_Writing_File( Log_File_Port,File_name,IOERROR ) ! Input variables : File_index,File_name,IOERROR
            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            ! Get the start date and time of the 'Main Processor'
            call DATE_AND_TIME( PRDA, TIME )
            PRETTY_TIME = TIME(1:2) // ':' // TIME(3:4) // ':' // TIME(5:10)  
            !--------------------------------------------------------------------------------------------------
            ! Output the running information on screen
            write(*,'(A)') " ====================================================================================="
            write(*,'(A)') " !                                                                                   !"
            write(*,'(A)') " !                               MSPH Program Version 1.1                            !"
            write(*,'(A)') " !                                Programer: Shanqin Jin                             !"
            write(*,'(4A)')' !                           START TIME: ',PRDA,PRETTY_TIME,'                        !'
            write(*,'(A)') " !                                                                                   !"
            write(*,'(A)') " ====================================================================================="

            ! Output the running information to log file
            write(Log_File_Port,'(A)') " ====================================================================================="
            write(Log_File_Port,'(A)') " !                                                                                   !"
            write(Log_File_Port,'(A)') " !                               MSPH Program Version 1.1                            !"
            write(Log_File_Port,'(A)') " !                                Programer: Shanqin Jin                             !"
            write(Log_File_Port,'(4A)')' !                           START TIME: ',PRDA,PRETTY_TIME,'                        !'
            write(Log_File_Port,'(A)') " !                                                                                   !"
            write(Log_File_Port,'(A)') " ====================================================================================="

        endif

        ! Synchronize all processors calculation
        call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
        !******************************************************************************************************

        !******************************************************************************************************
        ! The Main_Processor Output the error information
        if (Total_Processors_Number/=Desired_Processor_Number .and. Current_Processor_ID==Main_Processor) then

            write(*,20) Total_Processors_Number,Desired_Processor_Number
20          Format(" The input processor number(",I5," ) is not equal to desired processor number (",I5," )")
            
            write(*,30) Desired_Processor_Number
30          Format(" Please change the Variable, 'Desired_Processor_Number', to (",I5," )")
            
            call MPI_FINALIZE(ierror_MPI)
            
            goto 666                                                       ! Go to the end.

        endif

        ! Note the user the current active Processor ID on screen
        write(*,10) Current_Processor_ID,Total_Processors_Number,trim(adjustl(Current_Processor_Group_Name))
10      Format(" Processor",I5," of",I5," on ",A," is alive !")
        !******************************************************************************************************
        
    endif Program_Beginning
   
    call sleep( 1 )                                                        ! Sleep 1 second
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)                            ! Synchronize all processors calculation
    !----------------------------------------------------------------------------------------------------------





    !----------------------------------------------------------------------------------------------------------
    ! All Processors are running for the setting initial SPH Field
    call Initialize_SPH_Calculation()                                      ! Set the initial SPH Field

    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)                            ! Synchronize all processors calculation

    ! Get the start time of SPH iteration on the 'Main Processor'
    if ( Current_Processor_ID==Main_Processor ) then
        call CPU_TIME(Main_Processor_start_time)                           ! Call time Function to get the start time                      
    endif


    call SPH_Time_Iteration()                                              ! SPH time iteration            


    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)                            ! Synchronize all processors calculation

    ! Get the end time of SPH iteration on the 'Main Processor'
    if ( Current_Processor_ID==Main_Processor ) then
        call CPU_TIME(Main_Processor_end_time)                             ! Call time Function to get the end time                      
    endif
    !----------------------------------------------------------------------------------------------------------






    !----------------------------------------------------------------------------------------------------------
    ! Notice the user end of the program, calculate and save the running time on the 'Main Processor'
    Program_Ending: if (Current_Processor_ID==Main_Processor) then

        !--------------------Calculate the running time and show the time on screen----------------------------
        Main_Processor_cost_time = Main_Processor_end_time-Main_Processor_start_time
        Every_step_time          = Main_Processor_cost_time/max_time_step

        ! Get end date for the Main Processor
        call DATE_AND_TIME( FIDA, TIME )                                          ! The end time and date
        FINISH_TIME = TIME(1:2) // ':' // TIME(3:4) // ':' // TIME(5:10)

        ! Output the end notice on screen 
        write(*,'(A)') " ====================================================================================="
        write(*,'(3A)') ' Program start time: ', PRDA, PRETTY_TIME
        write(*,'(3A)') ' Program end time  : ', FIDA, FINISH_TIME
        write(*,'(A,F8.4,A)') " The time of the program for Main Processor cost is:",Main_Processor_cost_time,"s"
        write(*,'(A,F8.4,A)') " The average time cosuming of one step cost is     :",every_step_time,"s"
        write(*,'(3A)') ' The program run to the end, Please press "enter" to finish it!'
        write(*,'(A)') " ====================================================================================="

        ! Output the end notice to log file
        write(Log_File_Port,'(A)')        " ====================================================================================="
        write(Log_File_Port,'(3A)')       ' Program start time: ', PRDA, PRETTY_TIME
        write(Log_File_Port,'(3A)')       ' Program end time  : ', FIDA, FINISH_TIME
        write(Log_File_Port,'(A,F8.4,A)') " The time of the program for Main Processor cost is:",Main_Processor_cost_time,"s"
        write(Log_File_Port,'(A,F8.4,A)') " The average time cosuming of one step cost is     :",every_step_time,"s"
        write(Log_File_Port,'(3A)')       ' The program run to the end, Please press "enter" to finish it!'
        write(Log_File_Port,'(A)')        " ====================================================================================="
       
        close( Log_File_Port )                                                    ! Close the log file port 
        !------------------------------------------------------------------------------------------------------


        !--------------------------------Open the file to save the running time--------------------------------
        File_name  = "Compute_Time.dat"
        call Initialziting_Writing_File( General_File_Port,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
  
        if(IOERROR==0) then
            write(General_File_Port,*) "The time of the program for Main Processor cost is:",Main_Processor_cost_time
            write(General_File_Port,*) "The average time cosuming of one step cost is     :",every_step_time
        endif

        close(General_File_Port)
        !------------------------------------------------------------------------------------------------------

    endif Program_Ending
    !----------------------------------------------------------------------------------------------------------

    !**********************************************************************************************************



    ! Finish the MPI running and the program
666 call MPI_FINALIZE(ierror_MPI)
    
    stop


end program Parallel_Moving_Square_Tradictionary_Ghost