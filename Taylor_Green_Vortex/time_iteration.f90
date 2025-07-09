!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE:: Time iteration
!
!  PURPOSE: Leapforg, Eular or Runge–Kutta Iteration
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
!****************************************************************************


subroutine time_iteration()    
                            
    use Public_variable_module
    use information_module
    use MPI
   
    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    
    integer::i,j,k,L,m                                           !loop Variables
    integer::i_time_step=0                                       !current time step
    character(len=20)::filename                                  !filename
  
    !Variables for file opeartion
    integer::ioerror=0                                           !open return value
    integer::stat                                                !read return value
    integer::status                                              !allocation return value

    integer::file_index
    character(len=100)::file_name
    character(len=4)::char_Current_probe_order
    character(len=4)::char_Current_Processor_ID

    !Variables for folder opeartion
    Logical::dirExists

    !==========================================================================================================
        
    ! Body of subroutine time_integration

    ! Call MPI functions
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )

    !Initialize the arraies for iteration
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



    !Only the main processor output the data
    if (Current_Processor_ID==Main_Processor) then

      !***************************************************************************************************

      !This is for linux, on Windows the foler opeartion is not same
      if(trim(adjustl(startFrom))=='InitialTime') then
          
          !----------------------------------------------------------------------------------------------------------------
          !Wave sampling files
          do i=1,Wave_Probe_number

              Wave_Probe_File_index(i)=Wave_Probe_File_Initial_index+i     !Wave_Probe_File_index=Wave_Probe_File_Initial_index+wave_probe_order
              
              !transfer the Current_Processor_ID from integer to character
              write(char_Current_probe_order,'(I4)') i

              file_name="./Sampling/WaveHeight/wave_gauge_"//trim(adjustl(char_Current_probe_order))//"_result.dat"

              open(unit=Wave_Probe_File_index(i),file=file_name,status="replace",position="rewind",action="write",iostat=ioerror)

          end do
          !----------------------------------------------------------------------------------------------------------------

          !----------------------------------------------------------------------------------------------------------------
          !Pressure sampling files
          do i=1,Sampling_Point_number

              Normalized_Pressure_Sampling_File_index(i)=Normalized_Pressure_Sampling_File_Initial_index+i     !Normalized_Pressure_Sampling_File_index=Normalized_Pressure_Sampling_File_Initial_index+sampling_probe_order
              Average_Pressure_Sampling_File_index(i)=Average_Pressure_Sampling_File_Initial_index+i           !Average_Pressure_Sampling_File_index=Average_Pressure_Sampling_File_Initial_index+sampling_probe_order
              Kalman_Pressure_Sampling_File_index(i)=Kalman_Pressure_Sampling_File_Initial_index+i             !Kalman_Pressure_Sampling_File_index=Kalman_Pressure_Sampling_File_Initial_index+sampling_probe_order

              !transfer the Current_Processor_ID from integer to character
              write(char_Current_probe_order,'(I4)') i

              file_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Normalized_Pressure.dat"
              open(unit=Normalized_Pressure_Sampling_File_index(i),file=file_name,status="replace",position="rewind",action="write",iostat=ioerror)

              file_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Average_Pressure.dat"
              open(unit=Average_Pressure_Sampling_File_index(i),file=file_name,status="replace",position="rewind",action="write",iostat=ioerror)

              file_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Kalman_Pressure.dat"
              open(unit=Kalman_Pressure_Sampling_File_index(i),file=file_name,status="replace",position="rewind",action="write",iostat=ioerror)

          end do
          !----------------------------------------------------------------------------------------------------------------

          open(unit=119,file="Latest_Time_Results.dat",status="replace",position="rewind",action="write",iostat=ioerror)

          open(unit=10,file="Free_Surface.dat",status="replace",position="rewind",action="write",iostat=ioerror)

          write(10,*) "TITLE='DISTRIBUTION'"
          write(10,*) "VARIABLES= 'X' 'Y' 'Nor_X' 'Nor_Y' "
                    
          !Tecplot header
          write(119,*) "TITLE='DISTRIBUTION'"
          write(119,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'Rho' 'Mass' 'Ene' 'C' 'Type' 'Wmak_L' 'Step' "

          !Results files
          open(unit=2,file="video_result.dat",status="replace",position="rewind",action="write",iostat=ioerror) 
    
          !Tecplot header
          write(2,*) "TITLE='DISTRIBUTION'"
          write(2,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'U' 'Vor' 'EVX' 'EVY' 'EP' 'EU' "

          open(unit=121,file="Post_Proceeding_Grid.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          !Tecplot header
          write(121,*) "TITLE='mesh'"
          write(121,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'U' 'Vor' 'EVX' 'EVY' 'EP' 'EU' 'Fai'"

          !Wave front results
          open(unit=127,file="Wave_front_results.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          !kinetic energy
          open(unit=128,file="Kinetic_energy.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          !maximum velocity magnitude
          open(unit=129,file="Maximum_velocity_magnitude.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          !norm results
          open(unit=130,file="Norm_Horizontal_velocity.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          open(unit=131,file="Norm_Vertical_velocity.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          open(unit=132,file="Norm_Pressure.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          open(unit=133,file="Norm_Velocity_magnitude.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          !----------------------------------------------------------------------------------------------------------------
      
      else if(trim(adjustl(startFrom))=='LatestTime') then
          
          !----------------------------------------------------------------------------------------------------------------
          !Wave sampling files
          do i=1,Wave_Probe_number

              Wave_Probe_File_index(i)=Wave_Probe_File_Initial_index+i     !Wave_Probe_File_index=Wave_Probe_File_Initial_index+wave_probe_order
              
              !transfer the Current_Processor_ID from integer to character
              write(char_Current_probe_order,'(I4)') i

              file_name="./Sampling/WaveHeight/wave_gauge_"//trim(adjustl(char_Current_probe_order))//"_result.dat"

              open(unit=Wave_Probe_File_index(i),file=file_name,status="replace",position="rewind",action="write",iostat=ioerror)

          end do
          !----------------------------------------------------------------------------------------------------------------

          !----------------------------------------------------------------------------------------------------------------
          !Pressure sampling files
          do i=1,Sampling_Point_number

              Normalized_Pressure_Sampling_File_index(i)=Normalized_Pressure_Sampling_File_Initial_index+i     !Normalized_Pressure_Sampling_File_index=Normalized_Pressure_Sampling_File_Initial_index+sampling_probe_order
              Average_Pressure_Sampling_File_index(i)=Average_Pressure_Sampling_File_Initial_index+i           !Average_Pressure_Sampling_File_index=Average_Pressure_Sampling_File_Initial_index+sampling_probe_order
              Kalman_Pressure_Sampling_File_index(i)=Kalman_Pressure_Sampling_File_Initial_index+i             !Kalman_Pressure_Sampling_File_index=Kalman_Pressure_Sampling_File_Initial_index+sampling_probe_order

              !transfer the Current_Processor_ID from integer to character
              write(char_Current_probe_order,'(I4)') i

              file_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Normalized_Pressure.dat"
              open(unit=Normalized_Pressure_Sampling_File_index(i),file=file_name,status="replace",position="rewind",action="write",iostat=ioerror)

              file_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Average_Pressure.dat"
              open(unit=Average_Pressure_Sampling_File_index(i),file=file_name,status="replace",position="rewind",action="write",iostat=ioerror)

              file_name="./Sampling/Press/Probe_"//trim(adjustl(char_Current_probe_order))//"_Kalman_Pressure.dat"
              open(unit=Kalman_Pressure_Sampling_File_index(i),file=file_name,status="replace",position="rewind",action="write",iostat=ioerror)

          end do
          !----------------------------------------------------------------------------------------------------------------

          open(unit=119,file="Latest_Time_Results_latest.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          open(unit=10,file="Free_Surface.dat",status="replace",position="rewind",action="write",iostat=ioerror)

          write(10,*) "TITLE='DISTRIBUTION'"
          write(10,*) "VARIABLES= 'X' 'Y' 'Nor_X' 'Nor_Y' "
          
          !Tecplot header
          write(119,*) "TITLE='DISTRIBUTION'"
          write(119,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'Rho' 'Mass' 'Ene' 'C' 'Type' 'Wmak_L' 'Step' "


          !Results files
          open(unit=2,file="video_result.dat",status="replace",position="rewind",action="write",iostat=ioerror) 
    
          !Tecplot header
          write(2,*) "TITLE='DISTRIBUTION'"
          write(2,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'U' 'Vor' 'EVX' 'EVY' 'EP' 'EU' "

          open(unit=121,file="Post_Proceeding_Grid.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          !Tecplot header
          write(121,*) "TITLE='mesh'"
          write(121,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'U' 'Vor' 'EVX' 'EVY' 'EP' 'EU' 'Fai'"

          !Wave front results
          open(unit=127,file="Wave_front_results.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          !kinetic energy
          open(unit=128,file="Kinetic_energy.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          !maximum velocity magnitude
          open(unit=129,file="Maximum_velocity_magnitude.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          !norm results
          open(unit=130,file="Norm_Horizontal_velocity.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          open(unit=131,file="Norm_Vertical_velocity.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          open(unit=132,file="Norm_Pressure.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          open(unit=133,file="Norm_Velocity_magnitude.dat",status="replace",position="rewind",action="write",iostat=ioerror)
          
          !----------------------------------------------------------------------------------------------------------------

      else

          write(*,*) "The startFrom variable is not right. ('InitialTime' or 'latestTime')"
          
      end if

    end if
    !***************************************************************************************************
    !Synchronize all processors to start calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

    !------------------------------------------------------------------------------------------
    !open the output tecplot data files
    file_index=400+Current_Processor_ID

    !transfer the Current_Processor_ID from integer to character
    write(char_Current_Processor_ID,'(I4)') Current_Processor_ID
    !write(*,*) char_Current_Processor_ID

    file_name="./Subdomain/Results_in_subdomain_"//trim(adjustl(char_Current_Processor_ID))//".dat"
    !write(*,*) file_name

    open(unit=file_index,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    

    !tecplot header
    write(file_index,*) "TITLE='DISTRIBUTION'"
    write(file_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "

    !***************************************************************************************************
    ! Synchronize all processors to start calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

    !***************************************************************************************************
    !Leapforg or Eular iteration
    if (leapforg_or_eular==1) then
   
       do i_time_step=actual_start_time_step,max_time_step

          !Synchronize all processors calculation
          !call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

          Real_time=i_time_step*dt
          Nondimensionalized_real_time=Real_time/time_0     
    
          !output current time step
          if (Current_Processor_ID==Main_Processor) then
             write(*,10) i_time_step,Real_time,Nondimensionalized_real_time
10           Format("Current step:",I8," Real Time:",F8.4," Nondimensionalized Time:",F8.4)
          endif

          !---------------------------------------------------------------------------------------------------------
          !Wave maker
          if ( make_wave_or_not==1 ) then
             call Wavemaker_Subroutine(i_time_step)
          endif

          !Distribute Period Boundary Particle
          if (DistributePeriodBoundaryParticleOrNot==1) then
              call Distribute_Period_Boundary_Particle(i_time_step)
          endif

          !Devide blocks and send tasks
          call Devide_block_and_send_tasks(i_time_step)
          !---------------------------------------------------------------------------------------------------------
    
          !---------------------------------------------------------------------------------------------------------
          ! Use the information from last step and moving forward 0.5dt
          if (i_time_step/=1) then

             do i=1,actual_particle_number_in_subdomain

                if(subdomain_particle_type(i)/=2) cycle
           
                !Save the density and energy from last step
                before_rho(i)=subdomain_particle_rho(i)
                before_energy(i)=subdomain_particle_energy(i)
                
                !Renew density and energy (half step) 
                subdomain_particle_rho(i)=subdomain_particle_rho(i)+(dt/2.0)*drhodt(i)
                subdomain_particle_energy(i)=subdomain_particle_energy(i)+(dt/2.0)*dedt(i)
               
                !Renew velocity (half step) 
                do j=1,dim
                   before_velocity(i,j)=subdomain_particle_velocity(i,j)
                   subdomain_particle_velocity(i,j)=subdomain_particle_velocity(i,j)+(dt/2.0)*particle_acceleration(i,j)
                end do

             end do

          end if
          !---------------------------------------------------------------------------------------------------------

          !---------------------------------------------------------------------------------------------------------
          !Exchange buffer domain information
          !***As the exchange will add errors, so it may cause unsteady for higher-order calculation
          !***We can ignore the changes from the Prediction step

          call exchange_buffer_domain_information(i_time_step)

          !calculate iteration information
          call single_step_compute(i_time_step)
          !---------------------------------------------------------------------------------------------------------

          !---------------------------------------------------------------------------------------------------------
          if (i_time_step==1) then


             do i=1,actual_particle_number_in_subdomain                 

                if(subdomain_particle_type(i)/=2) cycle

                !Renew density and energy (half step) 
                subdomain_particle_rho(i)=subdomain_particle_rho(i)+(dt/2.0)*drhodt(i)
                subdomain_particle_energy(i)=subdomain_particle_energy(i)+(dt/2.0)*dedt(i)
                
                
                !Renew velocity and position (half step) 
                do j=1,dim
                   subdomain_particle_velocity(i,j)=subdomain_particle_velocity(i,j)+(dt/2.0)*particle_acceleration(i,j)+average_velocity(i,j)
                   subdomain_particle_position(i,j)=subdomain_particle_position(i,j)+dt*subdomain_particle_velocity(i,j)+0.5*particle_acceleration(i,j)*dt*dt
                end do
          
              end do
              
          else
        
              do i=1,actual_particle_number_in_subdomain

                 if(subdomain_particle_type(i)/=2) cycle
          
                 !Renew density and energy (Full step) 
                 subdomain_particle_rho(i)=before_rho(i)+dt*drhodt(i)
                 subdomain_particle_energy(i)=before_energy(i)+dt*dedt(i)
                 
                 !Renew velocity and position (Full step) 
                 do j=1,dim
                    subdomain_particle_velocity(i,j)=before_velocity(i,j)+dt*particle_acceleration(i,j)+average_velocity(i,j)
                    subdomain_particle_position(i,j)=subdomain_particle_position(i,j)+dt*subdomain_particle_velocity(i,j)+0.5*particle_acceleration(i,j)*dt*dt
                 end do
              end do
       
          end if
          !---------------------------------------------------------------------------------------------------------

          !---------------------------------------------------------------------------------------------------------
          !Fliter the reults by using the Wavemaker damping ZONE
          !Parts of the Variables are replaced by using the analytical results
          if (trim(adjustl(WavemakerDamping_Switch))=='on' .and. trim(adjustl(wave_maker_type))=='WaveTheory' .and. make_wave_or_not==1) then

            call wave_maker_damping_zone(i_time_step)

          end if

          !---------------------------------------------------------------------------------------------------------

          !---------------------------------------------------------------------------------------------------------
          !Reinitialize particle density by using MLS
          if (trim(adjustl(RefreshDensityByMLS_Switch))=='on' .and. mod(i_time_step,MLS_time_step)==0) then
            call Reinitialize_particle_density(i_time_step)
          endif
          !---------------------------------------------------------------------------------------------------------

          !---------------------------------------------------------------------------------------------------------
          !Shift particle position
          if (trim(adjustl(ShiftParticle_Switch))=='on') then
            call Shift_particle_position(i_time_step)
          endif
          !---------------------------------------------------------------------------------------------------------

          !Refresh the domain
          call Refresh_domain(i_time_step)

          !Calculate the data for Post-Proceeding
          if(mod(i_time_step,save_time_step)==0 .and. trim(adjustl(OutputForPostProcedingGrid_Switch))=='on') then
            call Calculation_For_PostProceeding(i_time_step)
          endif

          ! Synchronize all processors calculation
          ! MPI_BCAST function has the character "Synchronize", so we don't need Synchronize
          !call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
          !---------------------------------------------------------------------------------------------------------

          ! !**************************************************************************************************************
          ! !Save fluid information of subdomain for debugging
          ! if(mod(i_time_step,save_time_step)==0) then
             
          !   write(file_index,*) "ZONE I=",total_particle_number_in_subdomain," F=POINT"
          !   do i=1,total_particle_number_in_subdomain
          !      write(file_index,100) (subdomain_particle_position(i,j),j=1,dim),(subdomain_particle_velocity(i,j),j=1,dim),subdomain_particle_rho(i),subdomain_particle_press(i)/pressure_0
          !   end do
              
          ! end if
          ! !**************************************************************************************************************


          if (Current_Processor_ID==Main_Processor) then
   
            !**************************************************************************************************************
            !Save fluid information
            if(mod(i_time_step,save_time_step)==0) then

                !----------------------------------------------------------------------------------------------------------
                !For tecplot
                if (trim(adjustl(OutputForTecplot_Switch))=='on') then
                   call Output_for_tecplot(i_time_step)
                endif
                !----------------------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------------------
                !Save the data for Paraview
                if (trim(adjustl(OutputForParaview_Switch))=='on') then
                   call Output_for_paraview(i_time_step)
                endif
                !----------------------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------------------
                !Save the grid data for Post-Proceeding
                if (trim(adjustl(OutputForPostProcedingGrid_Switch))=='on') then
                   call Output_for_Grid_PostProceding(i_time_step)
                endif
                !----------------------------------------------------------------------------------------------------------

            end if
            !**************************************************************************************************************
            
            !**************************************************************************************************************
            !Sampling
            if (mod(i_time_step,sampling_time_step)==0) then
                call Sampling_subroutine(i_time_step)
            end if
            !**************************************************************************************************************

            !**************************************************************************************************************
            !Save Results for Debugging
            if (mod(i_time_step,save_time_step_for_debug)==0) then

                write(119,*) "ZONE I=",ture_total_particle_number," F=POINT"
                do i=1,ture_total_particle_number
                    write(119,500) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_press(i),particle_rho(i),particle_mass(i),particle_energy(i),particle_c(i),particle_initial_type(i),wave_maker_particle_layer(i),i_time_step          
500                 format(9F20.8,3I) 
                end do
                
            end if
            !**************************************************************************************************************

          end if
  
      end do

  !*******************************************************************************************************************
  elseif (leapforg_or_eular==2) then
          
     !Eular iteration
     do i_time_step=actual_start_time_step,max_time_step

  
        Real_time=i_time_step*dt
        Nondimensionalized_real_time=Real_time/time_0     
  
        !output current time step
        if (Current_Processor_ID==Main_Processor) then
           write(*,10) i_time_step,Real_time,Nondimensionalized_real_time
        endif

        ! !Synchronize all processors calculation
        ! call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

        !-----------------------------------------------------------------------------------------
        !Wave maker
        if ( make_wave_or_not==1 ) then
           call Wavemaker_Subroutine(i_time_step)
        endif

        !Distribute Period Boundary Particle
        if (DistributePeriodBoundaryParticleOrNot==1) then
            call Distribute_Period_Boundary_Particle(i_time_step)
        endif

        !Devide blocks and send tasks
        call Devide_block_and_send_tasks(i_time_step)
        !-----------------------------------------------------------------------------------------

        if (i_time_step==1) then

            !*************************************************************************************
            !Save last step information
            before_rho=subdomain_particle_rho                           
            before_velocity=subdomain_particle_velocity                 
            before_position=subdomain_particle_position                 
            before_energy=subdomain_particle_energy                     
            !*************************************************************************************

           !Prediction step
           call single_step_compute(i_time_step)
           
           do i=1,actual_particle_number_in_subdomain                 

              if(subdomain_particle_type(i)/=2) cycle
              
              !*************************************************************************************
              subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt      
              subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt  
              
              !Renew velocity and position
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+dt*particle_acceleration(i,j)+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+dt*subdomain_particle_velocity(i,j)
              end do
              !*************************************************************************************

           end do

          !*************************************************************************************
          !Save Prediction step information
          temp_particle_velocity=subdomain_particle_velocity
          temp_drhodt=drhodt
          temp_dedt=dedt
          temp_particle_acceleration=particle_acceleration
          !*************************************************************************************

          !Exchange buffer domain information
          !***As the exchange will add errors, so it may cause unsteady for higher-order calculation
          !***We can ignore the changes from the Prediction step

          call exchange_buffer_domain_information(i_time_step)

          !Correction step
          call single_step_compute(i_time_step)
          
          do i=1,actual_particle_number_in_subdomain             

              if(subdomain_particle_type(i)/=2) cycle
              
              subdomain_particle_rho(i)=before_rho(i)+0.5*(drhodt(i)+temp_drhodt(i))*dt   
              subdomain_particle_energy(i)=before_energy(i)+0.5*(dedt(i)+temp_dedt(i))*dt  
              
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+0.5*dt*(particle_acceleration(i,j)+temp_particle_acceleration(i,j))+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+0.5*dt*(subdomain_particle_velocity(i,j)+temp_particle_velocity(i,j))+0.5*dt*dt*particle_acceleration(i,j)
              end do
          end do

        else

            !**************************************************************************************
            !Save last step information
            before_rho=subdomain_particle_rho                           
            before_velocity=subdomain_particle_velocity                 
            before_position=subdomain_particle_position                 
            before_energy=subdomain_particle_energy                     
            !**************************************************************************************

           do i=1,actual_particle_number_in_subdomain                 

              if(subdomain_particle_type(i)/=2) cycle
              
              !*************************************************************************************
              subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt     
              subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt  
              
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+dt*particle_acceleration(i,j)+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+dt*subdomain_particle_velocity(i,j)
              end do
              !*************************************************************************************

           end do

          !*************************************************************************************
          !Save Prediction step information
          temp_particle_velocity=subdomain_particle_velocity
          temp_drhodt=drhodt
          temp_dedt=dedt
          temp_particle_acceleration=particle_acceleration
          !*************************************************************************************

          !Exchange buffer domain information
          !***As the exchange will add errors, so it may cause unsteady for the higher-order calculation
          !***We can ignore the changes from the Prediction step

          call exchange_buffer_domain_information(i_time_step)

          !Correction step
          call single_step_compute(i_time_step)
          
          do i=1,actual_particle_number_in_subdomain             

              if(subdomain_particle_type(i)/=2) cycle
              
              subdomain_particle_rho(i)=before_rho(i)+0.5*(drhodt(i)+temp_drhodt(i))*dt    
              subdomain_particle_energy(i)=before_energy(i)+0.5*(dedt(i)+temp_dedt(i))*dt  
              
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+0.5*dt*(particle_acceleration(i,j)+temp_particle_acceleration(i,j))+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+0.5*dt*(subdomain_particle_velocity(i,j)+temp_particle_velocity(i,j))+0.5*dt*dt*particle_acceleration(i,j)
              end do
          end do
 
        end if
        !---------------------------------------------------------------------------------------------------------

        !---------------------------------------------------------------------------------------------------------
        !Fliter the reults by using the Wavemaker damping ZONE
        !Parts of the Variables are replaced by using the analytical results
        if (trim(adjustl(WavemakerDamping_Switch))=='on' .and. trim(adjustl(wave_maker_type))=='WaveTheory' .and. make_wave_or_not==1 ) then

          call wave_maker_damping_zone(i_time_step)

        end if


        !---------------------------------------------------------------------------------------------------------

        !---------------------------------------------------------------------------------------------------------
        !Reinitialize particle density by using MLS
        if (trim(adjustl(RefreshDensityByMLS_Switch))=='on' .and. mod(i_time_step,MLS_time_step)==0) then
          call Reinitialize_particle_density(i_time_step)
        endif
        !---------------------------------------------------------------------------------------------------------

        !---------------------------------------------------------------------------------------------------------
        !Shift particle position
        if (trim(adjustl(ShiftParticle_Switch))=='on') then
          call Shift_particle_position(i_time_step)
        endif
        !---------------------------------------------------------------------------------------------------------
            
        !---------------------------------------------------------------------------------------------------------
        ! Refresh the domain
        call Refresh_domain(i_time_step)

        !Calculate the data for Post-Proceeding
        if(mod(i_time_step,save_time_step)==0 .and. trim(adjustl(OutputForPostProcedingGrid_Switch))=='on') then
          call Calculation_For_PostProceeding(i_time_step)
        endif

        ! Synchronize all processors calculation
        ! MPI_BCAST function has the character "Synchronize", so we don't need Synchronize
        !call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
        !---------------------------------------------------------------------------------------------------------
      
  !     !Note the user the current active Processor ID have finished the iterator information calculation
  !      write(*,310) Current_Processor_ID
  ! 310  Format(" Processor:",I3," has finished Eular interation!") 
        !---------------------------------------------------------------------------------------------------------  


        ! !**************************************************************************************************************
        ! !Save fluid information of subdomain for debugging
        ! if(mod(i_time_step,save_time_step)==0) then
           
        !   write(file_index,*) "ZONE I=",total_particle_number_in_subdomain," F=POINT"
        !   do i=1,total_particle_number_in_subdomain
        !      write(file_index,100) (subdomain_particle_position(i,j)/domain_size_y,j=1,dim),(subdomain_particle_velocity(i,j)/velocity_0,j=1,dim),subdomain_particle_rho(i),subdomain_particle_press(i)/pressure_0
        !   end do
            
        ! end if
        ! !**************************************************************************************************************


       !The main processor output the data
       if (Current_Processor_ID==Main_Processor) then
       
          !*************************************************************************************************************
          !Save fluid information
          if(mod(i_time_step,save_time_step)==0) then

            !----------------------------------------------------------------------------------------------------------
            !For tecplot
            if (trim(adjustl(OutputForTecplot_Switch))=='on') then
               call Output_for_tecplot(i_time_step)
            endif
            !----------------------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------------------
            !Save the data for Paraview
            if (trim(adjustl(OutputForParaview_Switch))=='on') then
                call Output_for_paraview(i_time_step)
            endif
            !----------------------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------------------
            !Save the grid data for Post-Proceeding
            if (trim(adjustl(OutputForPostProcedingGrid_Switch))=='on') then
                call Output_for_Grid_PostProceding(i_time_step)
            endif
            !----------------------------------------------------------------------------------------------------------

          end if
          !*************************************************************************************************************

          !*************************************************************************************************************
          !Sampling
          if(mod(i_time_step,sampling_time_step)==0) then
              call Sampling_subroutine(i_time_step)
          end if
          !*************************************************************************************************************
          
          !*************************************************************************************************************
          !Save Results for Debugging
          if(mod(i_time_step,save_time_step_for_debug)==0) then

              write(119,*) "ZONE I=",ture_total_particle_number," F=POINT"
              do i=1,ture_total_particle_number
                  write(119,500) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_press(i),particle_rho(i),particle_mass(i),particle_energy(i),particle_c(i),particle_initial_type(i),wave_maker_particle_layer(i),i_time_step          
              end do
              
          end if
          !*************************************************************************************************************

        end if

     end do

     !******************************************************************************************************************
  
  else 

     !******************************************************************************************************************
          
     !Runge–Kutta iteration

     do i_time_step=actual_start_time_step,max_time_step


        Real_time=i_time_step*dt
        Nondimensionalized_real_time=Real_time/time_0     
  
        !output current time step
        if (Current_Processor_ID==Main_Processor) then
           write(*,10) i_time_step,Real_time,Nondimensionalized_real_time
        endif

        ! !Synchronize all processors calculation
        ! call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

        !-----------------------------------------------------------------------------------------
        !Wave maker
        if ( make_wave_or_not==1 ) then
            call Wavemaker_Subroutine(i_time_step)
        endif

        !Distribute Period Boundary Particle
        if (DistributePeriodBoundaryParticleOrNot==1) then
            call Distribute_Period_Boundary_Particle(i_time_step)
        endif

        !Devide blocks and send tasks
        call Devide_block_and_send_tasks(i_time_step)
        !-----------------------------------------------------------------------------------------

        if (i_time_step==1) then

            !*************************************************************************************
            !Save last step information (K_1)
            before_rho=subdomain_particle_rho                           
            before_velocity=subdomain_particle_velocity                 
            before_position=subdomain_particle_position                 
            before_energy=subdomain_particle_energy                     
            !*************************************************************************************

           !Prediction step (Calculate K_1 and update postion with K_1 acceleration and 0.5dt for K_2)
           call single_step_compute(i_time_step)
           
           do i=1,actual_particle_number_in_subdomain                 

              if(subdomain_particle_type(i)/=2) cycle
              
              !*************************************************************************************
              subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt/2      
              subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt/2  
              
              !Renew velocity and position
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+0.5*dt*particle_acceleration(i,j)+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+0.5*dt*subdomain_particle_velocity(i,j)
              end do
              !*************************************************************************************

           end do

          !*************************************************************************************
          !Save Prediction step information (Save K_1)
          K_1_temp_particle_velocity=subdomain_particle_velocity
          K_1_temp_drhodt=drhodt
          K_1_temp_dedt=dedt
          K_1_temp_particle_acceleration=particle_acceleration
          !*************************************************************************************

          !Exchange buffer domain information
          !***As the exchange will add errors, so it may cause unsteady for higher-order calculation
          !***We can ignore the changes from the Prediction step

          call exchange_buffer_domain_information(i_time_step)

          !Prediction step (Calculate K_2 and update postion with K_2 acceleration and 0.5dt for K_3)
           call single_step_compute(i_time_step)
           
           do i=1,actual_particle_number_in_subdomain                 

              if(subdomain_particle_type(i)/=2) cycle
              
              !*************************************************************************************
              subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt/2      
              subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt/2  
              
              !Renew velocity and position
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+0.5*dt*particle_acceleration(i,j)+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+0.5*dt*subdomain_particle_velocity(i,j)
              end do
              !*************************************************************************************

           end do

          !*************************************************************************************
          !Save Prediction step information (Save K_2)
          K_2_temp_particle_velocity=subdomain_particle_velocity
          K_2_temp_drhodt=drhodt
          K_2_temp_dedt=dedt
          K_2_temp_particle_acceleration=particle_acceleration
          !*************************************************************************************

          call exchange_buffer_domain_information(i_time_step)

          !Prediction step (Calculate K_3 and update postion with K_3 acceleration and dt for K_4)
          call single_step_compute(i_time_step)

           do i=1,actual_particle_number_in_subdomain                 

              if(subdomain_particle_type(i)/=2) cycle
              
              !*************************************************************************************
              subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt     
              subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt  
              
              !Renew velocity and position
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+dt*particle_acceleration(i,j)+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+dt*subdomain_particle_velocity(i,j)
              end do
              !*************************************************************************************

           end do

          !*************************************************************************************
          !Save Prediction step information (Save K_3)
          K_3_temp_particle_velocity=subdomain_particle_velocity
          K_3_temp_drhodt=drhodt
          K_3_temp_dedt=dedt
          K_3_temp_particle_acceleration=particle_acceleration
          !*************************************************************************************

          call exchange_buffer_domain_information(i_time_step)

          !Prediction step (Calculate K_4 and update final postion)
          call single_step_compute(i_time_step)
          
           do i=1,actual_particle_number_in_subdomain             

              if(subdomain_particle_type(i)/=2) cycle
              
              subdomain_particle_rho(i)=before_rho(i)+(K_1_temp_drhodt(i)+2*K_2_temp_drhodt(i)+2*K_3_temp_drhodt(i)+drhodt(i))*dt/6.0d0   
              subdomain_particle_energy(i)=before_energy(i)+(K_1_temp_dedt(i)+2*K_2_temp_dedt(i)+2*K_3_temp_dedt(i)+dedt(i))*dt/6.0d0   
              
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+(K_1_temp_particle_acceleration(i,j)+2*K_2_temp_particle_acceleration(i,j)+2*K_3_temp_particle_acceleration(i,j)+particle_acceleration(i,j))*dt/6.0d0+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+(K_1_temp_particle_velocity(i,j)+2*K_2_temp_particle_velocity(i,j)+2*K_3_temp_particle_velocity(i,j)+subdomain_particle_velocity(i,j))*dt/6.0d0+0.5*dt*dt*particle_acceleration(i,j)
              end do
            end do

        else

            !*************************************************************************************
            !Save last step information (K_1)
            before_rho=subdomain_particle_rho                           
            before_velocity=subdomain_particle_velocity                 
            before_position=subdomain_particle_position                 
            before_energy=subdomain_particle_energy                     
            !*************************************************************************************

           !Prediction step (update postion with K_1 acceleration in last step and 0.5dt for K_2)
           
           do i=1,actual_particle_number_in_subdomain                 

              if(subdomain_particle_type(i)/=2) cycle
              
              !*************************************************************************************
              subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt/2      
              subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt/2  
              
              !Renew velocity and position
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+0.5*dt*particle_acceleration(i,j)+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+0.5*dt*subdomain_particle_velocity(i,j)
              end do
              !*************************************************************************************

           end do

          !*************************************************************************************
          !Save Prediction step information (Save K_1)
          K_1_temp_particle_velocity=subdomain_particle_velocity
          K_1_temp_drhodt=drhodt
          K_1_temp_dedt=dedt
          K_1_temp_particle_acceleration=particle_acceleration
          !*************************************************************************************

          !Exchange buffer domain information
          !***As the exchange will add errors, so it may cause unsteady for higher-order calculation
          !***We can ignore the changes from the Prediction step

          call exchange_buffer_domain_information(i_time_step)

          !Prediction step (Calculate K_2 and update postion with K_2 acceleration and 0.5dt for K_3)
           call single_step_compute(i_time_step)
           
           do i=1,actual_particle_number_in_subdomain                 

              if(subdomain_particle_type(i)/=2) cycle
              
              !*************************************************************************************
              subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt/2      
              subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt/2  
              
              !Renew velocity and position
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+0.5*dt*particle_acceleration(i,j)+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+0.5*dt*subdomain_particle_velocity(i,j)
              end do
              !*************************************************************************************

           end do

          !*************************************************************************************
          !Save Prediction step information (Save K_2)
          K_2_temp_particle_velocity=subdomain_particle_velocity
          K_2_temp_drhodt=drhodt
          K_2_temp_dedt=dedt
          K_2_temp_particle_acceleration=particle_acceleration
          !*************************************************************************************

          call exchange_buffer_domain_information(i_time_step)

          !Prediction step (Calculate K_3 and update postion with K_3 acceleration and dt for K_4)
          call single_step_compute(i_time_step)

           do i=1,actual_particle_number_in_subdomain                 

              if(subdomain_particle_type(i)/=2) cycle
              
              !*************************************************************************************
              subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt     
              subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt  
              
              !Renew velocity and position
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+dt*particle_acceleration(i,j)+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+dt*subdomain_particle_velocity(i,j)
              end do
              !*************************************************************************************

           end do

          !*************************************************************************************
          !Save Prediction step information (Save K_3)
          K_3_temp_particle_velocity=subdomain_particle_velocity
          K_3_temp_drhodt=drhodt
          K_3_temp_dedt=dedt
          K_3_temp_particle_acceleration=particle_acceleration
          !*************************************************************************************

          call exchange_buffer_domain_information(i_time_step)

          !Prediction step (Calculate K_4 and update final postion)
          call single_step_compute(i_time_step)
          
          do i=1,actual_particle_number_in_subdomain             

              if(subdomain_particle_type(i)/=2) cycle
              
              subdomain_particle_rho(i)=before_rho(i)+(K_1_temp_drhodt(i)+2*K_2_temp_drhodt(i)+2*K_3_temp_drhodt(i)+drhodt(i))*dt/6.0d0   
              subdomain_particle_energy(i)=before_energy(i)+(K_1_temp_dedt(i)+2*K_2_temp_dedt(i)+2*K_3_temp_dedt(i)+dedt(i))*dt/6.0d0   
              
              do j=1,dim
                 subdomain_particle_velocity(i,j)=before_velocity(i,j)+(K_1_temp_particle_acceleration(i,j)+2*K_2_temp_particle_acceleration(i,j)+2*K_3_temp_particle_acceleration(i,j)+particle_acceleration(i,j))*dt/6.0d0+average_velocity(i,j)
                 subdomain_particle_position(i,j)=before_position(i,j)+(K_1_temp_particle_velocity(i,j)+2*K_2_temp_particle_velocity(i,j)+2*K_3_temp_particle_velocity(i,j)+subdomain_particle_velocity(i,j))*dt/6.0d0+0.5*dt*dt*particle_acceleration(i,j)
              end do
          end do

      end if
      !---------------------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------
      !Fliter the reults by using the Wavemaker damping ZONE
      !Parts of the Variables are replaced by using the analytical results
      if (trim(adjustl(WavemakerDamping_Switch))=='on' .and. trim(adjustl(wave_maker_type))=='WaveTheory' .and. make_wave_or_not==1 ) then

        call wave_maker_damping_zone(i_time_step)

      end if
      !---------------------------------------------------------------------------------------------------------
      
      !---------------------------------------------------------------------------------------------------------
      !Reinitialize particle density by using MLS
      if (trim(adjustl(RefreshDensityByMLS_Switch))=='on' .and. mod(i_time_step,MLS_time_step)==0) then
        call Reinitialize_particle_density(i_time_step)
      endif
      !---------------------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------
      !Shift particle position
      if (trim(adjustl(ShiftParticle_Switch))=='on') then
        call Shift_particle_position(i_time_step)
      endif
      !---------------------------------------------------------------------------------------------------------
          
      !---------------------------------------------------------------------------------------------------------
      ! Refresh the domain
      call Refresh_domain(i_time_step)

      !Calculate the data for Post-Proceeding
      if(mod(i_time_step,save_time_step)==0 .and. trim(adjustl(OutputForPostProcedingGrid_Switch))=='on') then
        call Calculation_For_PostProceeding(i_time_step)
      endif

      ! Synchronize all processors calculation
      ! MPI_BCAST function has the character "Synchronize", so we don't need Synchronize
      !call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
    
!     !Note the user the current active Processor ID have finished the iterator information calculation
!      write(*,310) Current_Processor_ID
! 310  Format(" Processor:",I3," has finished Eular interation!") 
      !---------------------------------------------------------------------------------------------------------  


      ! !**************************************************************************************************************
      ! !Save fluid information of subdomain for debugging
      ! if(mod(i_time_step,save_time_step)==0) then
         
      !   write(file_index,*) "ZONE I=",total_particle_number_in_subdomain," F=POINT"
      !   do i=1,total_particle_number_in_subdomain
      !      write(file_index,100) (subdomain_particle_position(i,j)/domain_size_y,j=1,dim),(subdomain_particle_velocity(i,j)/velocity_0,j=1,dim),subdomain_particle_rho(i),subdomain_particle_press(i)/pressure_0
      !   end do
          
      ! end if
      ! !**************************************************************************************************************

       !The main processor output the data
       if (Current_Processor_ID==Main_Processor) then
       
          !************************************************************************************************************
          !Save fluid information
          if(mod(i_time_step,save_time_step)==0) then

              !----------------------------------------------------------------------------------------------------------
              !For tecplot
              if (trim(adjustl(OutputForTecplot_Switch))=='on') then
                 call Output_for_tecplot(i_time_step)
              endif
              !----------------------------------------------------------------------------------------------------------

              !----------------------------------------------------------------------------------------------------------
              !Save the data for Paraview
              if (trim(adjustl(OutputForParaview_Switch))=='on') then
                  call Output_for_paraview(i_time_step)
              endif
              !----------------------------------------------------------------------------------------------------------

              !----------------------------------------------------------------------------------------------------------
              !Save the grid data for Post-Proceeding
              if (trim(adjustl(OutputForPostProcedingGrid_Switch))=='on') then
                  call Output_for_Grid_PostProceding(i_time_step)
              endif
              !----------------------------------------------------------------------------------------------------------

          end if
          !************************************************************************************************************

          !************************************************************************************************************
          !Sampling
          if(mod(i_time_step,sampling_time_step)==0) then
              call Sampling_subroutine(i_time_step)
          end if
          !************************************************************************************************************
          
          !************************************************************************************************************
          !Save Results for Debugging
          if(mod(i_time_step,save_time_step_for_debug)==0) then

              write(119,*) "ZONE I=",ture_total_particle_number," F=POINT"
              do i=1,ture_total_particle_number
                  write(119,500) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_press(i),particle_rho(i),particle_mass(i),particle_energy(i),particle_c(i),particle_initial_type(i),wave_maker_particle_layer(i),i_time_step          
              end do
              
          end if
          !************************************************************************************************************

        end if

     end do
         
  end if 
       

end subroutine time_iteration