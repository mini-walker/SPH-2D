!****************************************************************************
!
!  PROGRAM: Parallel_Taylor_Green_Vortices_2D_Release
!
!  PURPOSE: Taylor Green Vortices 2D Case Simulation (Parallel -Version)
!           Update : Exchange two column and change the boundary points calculation
!
!  Programer: Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location: MUN
!
!  Time£º2017.3.18
!
!  Copyright: Memorial University
!
!  Version: 1.0
!
!  Note: MPI version: mpich-3.2
!        Fortran version: Inter Fortran (ifort) 18.0.0
!****************************************************************************

  program Parallel_Taylor_Green_Vortices_2D_Release

    use information_module
    use Public_variable_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables in program

    integer::i,j,k,l                                                           !Variables for loop
    
    !Variables for file opeartion
    integer::ioerror=0                                                         !open return value
    integer::stat                                                              !read return value
    integer::status                                                            !allocation return value
  
    !==========================================================================================================

    ! Body of CSPH_MPI_Parallel_wave_maker_2D

    !********************************************************************************
    ! Initial the MPI runing system

    call MPI_INIT( ierror_MPI )
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    call MPI_GET_PROCESSOR_NAME(Current_Processor_Group_Name,Name_Length,ierror_MPI)

    if (ierror_MPI/=0) then

        write(*,*) "There are some errors in MPI Initial!"

    else

      !Note the usr the current active Processor ID
      write(*,10) Current_Processor_ID,Total_Processors_Number,Current_Processor_Group_Name
10    Format(" Processor:",I3," of",I3," on ",A10," is alive !")
        
    endif

    ! The Main_Processor Output the error information
    if (Total_Processors_Number/=Desired_Processor_Number .and. Current_Processor_ID==Main_Processor) then

        write(*,20) Total_Processors_Number,Desired_Processor_Number
20      Format(" The input processor number(",I3,") is not equal to desired processor number (",I3,")")
        
        write(*,30) Desired_Processor_Number
30      Format(" Please change the Variable, 'Desired_Processor_Number', to (',I3,')")
        
        goto 666                        !Go to the end.
        
    endif

    ! Synchronize all processors calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

    !********************************************************************************
    ! Call Date and time Function
    if (Current_Processor_ID==Main_Processor) then
        ! Get the start date and time of the Main Processor
        call DATE_AND_TIME( DATE, TIME )
        PRETTY_TIME = TIME(1:2) // ':' // TIME(3:4) // ':' // TIME(5:10)  
        write(*,'(3A)') ' START TIME: ', DATE, PRETTY_TIME
    endif
    !********************************************************************************
    

    !********************************************************************************
    ! Only the main Processor are running for the setting initial SPH Field
    ! Set the initial SPH Field
    call Initialize_SPH_Calculation()

    !Assign the initial particle position (At t=0)
    Initial_particle_Position=particle_position

    !********************************************************************************
    ! Check the input wave information and search the suitable wave model
    ! All the processors should finished this work as wave information are changed in this subroutine
    ! The master processor responses for the output
    if (trim(adjustl(wave_maker_type))=='WaveTheory') then
        
        !Check the wave model is right or not, then recommand the best wave model for the wave generation
        call Wave_Model_Verification()

        !----------------------------------------------------------------------------
        !Define the particles in wavemaker damping zone
        if (trim(adjustl(WavemakerDamping_Switch))=='on') then

            !Check the  wavemaker damping zone have changed or not
            if (abs(WavemakerDampingLengh-WavemakerDampingLenghScale*wave_length)<=1.0E-7) then
                
                if (Current_Processor_ID==1) then
                    write(*,*) "The length of wavemaker damping (",WavemakerDampingLengh,") zone is right"
                endif

                !Initialized the data for wavemaker damping 
                WaveMaker_Analytical_position=particle_position      !Analytical position in wavemaker damping zone
                WaveMaker_Analytical_velocity=particle_velocity      !Analytical velocity in wavemaker damping zone
                WaveMaker_Analytical_rho=particle_rho                !Analytical density in wavemaker damping zone
                WaveMaker_Analytical_press=particle_press            !Analytical pressure in wavemaker damping zone

                !Define the initial particles in wavemaker damping zone and jeep these same in the calculation
                In_WaveMakerDampingZone_OrNot=0                              !Initialize 
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
        !----------------------------------------------------------------------------

    endif

    !********************************************************************************
    ! Synchronize all processors calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
    !********************************************************************************

    !********************************************************************************
    ! Call time Function
    if (Current_Processor_ID==Main_Processor) then
        ! Get the start time of the Main Processor
        call CPU_TIME(Main_Processor_start_time)
    endif

    !********************************************************************************
    ! Start time iteration
    call time_iteration()              

    !********************************************************************************


    !********************************************************************************
    ! Synchronize all processors calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

    ! Main Processor running to the end
    if (Current_Processor_ID==Main_Processor) then

        !-------------Get the final field results of the compute domain--------------
        ! Open the saving file
        open(unit=1,file="Final_Doamin_Results.dat",status="replace",position="rewind",action="write",iostat=ioerror)    
        
        if(ioerror==0) then
            write(*,*) "The open action of the Final_Doamin_Results.dat is successful!"
            
            !Output the header of tecplot
            write(1,*) "TITLE='DISTRIBUTION'"
            write(1,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
                              
           !Output the results
            write(1,*) "ZONE I=",particle_ture_number," F=POINT"
            do i=1,particle_ture_number
               write(1,110) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)
110            format(8F20.10)            
            end do
      
        else
            write(*,*) "The open action of the Final_Doamin_Results.dat is fail!"
        end if
       
        close(1)
        !----------------------------------------------------------------------------


        !-------------------Open the file to save the running time-------------------
        open(unit=1,file="compute_time.dat",status="replace",position="rewind",action="write",iostat=ioerror)

        !Get the start time of the Main Processor
        call CPU_TIME(Main_Processor_end_time)
        Main_Processor_cost_time=Main_Processor_end_time-Main_Processor_start_time
        every_step_time=Main_Processor_cost_time/max_time_step

        !----------------------------------------------------------------------------
        ! Output  time of the start and end date for the Main Processor
        write(*,*) "==============================================================="
        write(*,'(3A)') ' Program start time: ', DATE, PRETTY_TIME
        
        call DATE_AND_TIME( DATE, TIME )                                             ! The end time and date
        PRETTY_TIME = TIME(1:2) // ':' // TIME(3:4) // ':' // TIME(5:10)  
        write(*,'(3A)') ' Program end time: ', DATE, PRETTY_TIME
        !----------------------------------------------------------------------------

        write(*,*) "The time of the program for Main Processor cost is:",Main_Processor_cost_time
        write(1,*) "The time of the program for Main Processor cost is:",Main_Processor_cost_time
        write(*,*) "The average time cosuming of one step cost is:",every_step_time
        write(1,*) "The average time cosuming of one step cost is:",every_step_time
        write(*,*) 'The program have run to the end,Please press "enter" to finish it!'
        write(*,*) "==============================================================="

        close(1)
        !----------------------------------------------------------------------------

    endif
    !********************************************************************************

    ! Finish the MPI running
666 call MPI_FINALIZE(ierror_MPI)

    stop

  end program Parallel_Taylor_Green_Vortices_2D_Release