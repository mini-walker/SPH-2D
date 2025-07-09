!**************************************************************************************************************
!  SUBROUTINE: SPH Time Iteration
!
!  PURPOSE   : Leapforg, Eular or Runge–Kutta Iteration
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

subroutine SPH_Time_Iteration()    
                            
    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
   
    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,L,m                                                                     ! Variables for looping
    integer::i_time_step                                                                   ! Present time step
  
    !==========================================================================================================
        
    ! Body of subroutine time_integration

    !**********************************************************************************************************
    ! Synchronize all processors to start calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
    


    !**********************************************************************************************************
    ! Iteration method
    ! Time iteration type: 1 is leapforg; 2 is Eular; 3 is Predictor-Corrector; 4 is Beeman; other is RUGGER KUTA
    ! First, use current velocity update the position, S(n+1)=V(n)*dt+S(n)
    ! Then use current accerlation update the velocity , V(n+1)=a(n)*dt+V(n)*dt
    ! As in this time iteration, it start from 1, so we save the data after SPH renew
    ! Define the real time
    if     ( StartFrom == 0 ) then                             
        Real_time=0.0d0
    elseif ( StartFrom == 1 ) then                 
        Real_time=Restart_Time                 
    else
        write(*,'(A)') " The 'startFrom' variable is not right. ('0---InitialTime' or '1---LatestTime')"
    end if



    TimeIterationDo: do i_time_step = Actual_start_time_step , Max_time_step

        !------------------------------------------------------------------------------------------------------
        ! Calculate the time step length
        dt=Fixed_dt 


        ! Calculate the real physical time and non-dimensional real physical time
        Real_time=Real_time+dt
        Real_physical_time=Real_time-Standing_time

        Nondimensionalized_real_time=Real_time/time_0  


        ! Define th in/outlet boundary ramping coefficient
        if ( 0.0d0 < Real_physical_time .and. Real_physical_time <= Boundary_ramping_time .and. Boundary_ramping_time > 1.0E-9 ) then
            InOutlet_Ramping_coefficient             = 0.5*( sin( ( Real_physical_time/Boundary_ramping_time-0.5 )*PI )+1 )
            InOutlet_Ramping_accleration_coefficient = 0.5*  cos( ( Real_physical_time/Boundary_ramping_time-0.5 )*PI )*PI/Boundary_ramping_time
        else
            InOutlet_Ramping_coefficient             = 1.0d0
            InOutlet_Ramping_accleration_coefficient = 0.0d0  
        endif
        !------------------------------------------------------------------------------------------------------



        !------------------------------------------------------------------------------------------------------
        ! Output current time step
        if ( Current_Processor_ID==Main_Processor ) then

            write(*,10) i_time_step, dt, Real_time, Nondimensionalized_real_time
10          Format(" Present step:",I8," Dt:",F10.6," Real time:",F10.6," Nondimensional time:",F10.6)

        endif 
        !------------------------------------------------------------------------------------------------------



        !------------------------------------------------------------------------------------------------------

        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! Attention: The variable, 'Ture_total_particle_number', is the sum of (1) Particle_ture_number;
        !                                                                      (2) Wave_maker_particle_number;
        !                                                                      (3) Fix_ghost_particle_number;
        !                                                                      (4) InletOutlet_Boundary_Particle_Number.
        ! (1), (2), (3) and (4) are nearly fixed in 'Body_Motion_Solver', 'Wavemaker_Subroutine' and 'Updating_InletOutlet_Boundary_Particle' 
        ! We just need to modify the number
        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        ! Body motion solver
        if ( Body_Motion_Solver_On_or_Off==1 ) then
            call Body_Motion_Solver(i_time_step)
        endif

        ! Wave-maker subroutine
        if ( Make_wave_or_not==1 ) then
            call Wavemaker_Subroutine(i_time_step)
        endif


        !------------------------------------------------------------------------------------------------------
        ! Dynamic external boundary condition distributed
        ! It seems that changing the boundary type in the running period, it will make unstable, shockwave will generated
        ! You'd better not changing the external boundary condition during the simulation, even the permeable boundary shut down suddenly
        ! A shockwave will generated in the domain
        ! May be shut down with a ramping function will decrease the shockwave, but I don't test
        if ( Real_physical_time>=InletOutletBoundaryRunningTime ) then
            DistributeInletOutletBoundaryParticleOrNot=0
        endif

        if ( Real_physical_time>=PeriodBoundaryRunningTime ) then
            DistributePeriodBoundaryParticleOrNot=0
        endif 
        
        if ( Real_physical_time>=TraditionalGhostBoundaryRunningTime ) then
            DistributeTraditionalGhostBoundaryParticleOrNot=0
        endif 

        ! Transfer the external boundary particles during the running time
        New_Fix_ghost_particle_number = 0
        
        if ( Real_physical_time>=ExternalBoundaryTransferTime .and. BoundaryTransferOpeartion==0 ) then

            call Transfer_External_Boundary(i_time_step) 

            ! The boundary transfer only do one time
            BoundaryTransferOpeartion=1

        endif


        ! Distribute InletOutlet Boundary Particle 
        ! Attentation: The period boundary condition is distributed after the InletOutlet boundary
        if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then
            
            ! Ramping velocity of the inlet boundary
            Ramping_Inlet_Velocity = InOutlet_Ramping_coefficient*Inflow_velocity

            call Updating_InletOutlet_Boundary_Particle(i_time_step)

        endif

        ! Reorder the particle index : Delete the ill particle and the oulet particle out the computational domain
        ! Renew particle number of each particle type
        call Reorder_Particle_Index(i_time_step)


        ! In 'Distribute_Period_Boundary_Particle' and 'Distribute_Traditional_Ghost_Boundary_Particle'
        ! The 'Period_Boundary_Particle_Number' and 'Traditional_ghost_particle_number' are dynamic changing
        ! They should be modified at every time step

        ! Distribute Period Boundary Particle
        if ( DistributePeriodBoundaryParticleOrNot==1 ) then

            call Distribute_Period_Boundary_Particle(i_time_step)

        endif


        ! Distribute Tradictional Ghost Boundary Particle
        if ( DistributeTraditionalGhostBoundaryParticleOrNot==1 ) then

            call Distribute_Traditional_Ghost_Boundary_Particle(i_time_step)

        endif

        !------------------------------------------------------------------------------------------------------



        ! Divide blocks and allot tasks
        call Divide_Block_and_Allot_Tasks(i_time_step)












        !======================================================================================================
      
        ! Using SPH for iteration, moving forward with dt (Update the information in the next step)
        IterationMethodIf: if ( Iteration_Method==1 ) then

            !**************************************************************************************************
            ! Leapforg Method

            !==================================================================================================
            ! Use the information from last step and moving forward 0.5dt to agree with the position which has been moving forward with dt
            if ( i_time_step /=actual_start_time_step ) then

                !----------------------------------------------------------------------------------------------
                ! SPH iteration 

                ! Save the density and energy from last step
                before_rho=subdomain_particle_rho
                before_energy=subdomain_particle_energy
                before_velocity=subdomain_particle_velocity

                do i=1,actual_particle_number_in_subdomain

                    if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then

                        ! Renew density and energy (half step) 
                        subdomain_particle_rho(i)=subdomain_particle_rho(i)+0.5*dt*drhodt(i)
                        subdomain_particle_energy(i)=subdomain_particle_energy(i)+0.5*dt*dedt(i)
                       
                        ! Renew velocity (half step) 
                        subdomain_particle_velocity(i,:)=subdomain_particle_velocity(i,:)+0.5*dt*particle_acceleration(i,:)

                    endif

                end do
                !----------------------------------------------------------------------------------------------

            end if
            !==================================================================================================


            !--------------------------------------------------------------------------------------------------
            ! Exchange buffer domain information
            ! As the exchange will add errors, so it may cause unsteady for higher-order calculation
            ! We can ignore the changes from the Prediction step

            call Exchange_buffer_domain_information(i_time_step)

            ! Calculate iteration information
            call Single_step_compute(i_time_step)
            !--------------------------------------------------------------------------------------------------


            !==================================================================================================
            ! For first step, the position moving forward full step but the others moving forward with half step.
            if ( i_time_step == actual_start_time_step ) then

                !----------------------------------------------------------------------------------------------
                ! SPH iteration 
                do i=1,actual_particle_number_in_subdomain                 

                    if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then

                        ! Renew density and energy (Half step) 
                        subdomain_particle_rho(i)=subdomain_particle_rho(i)+0.5*dt*drhodt(i)
                        subdomain_particle_energy(i)=subdomain_particle_energy(i)+0.5*dt*dedt(i)
                        
                        ! Renew position (Full step) 
                        subdomain_particle_position(i,:)=subdomain_particle_position(i,:)+dt*( subdomain_particle_velocity(i,:)+average_velocity(i,:) )+0.5*particle_acceleration(i,:)*dt**2

                        ! Renew velocity (Half step) 
                        subdomain_particle_velocity(i,:)=subdomain_particle_velocity(i,:)+0.5*dt*particle_acceleration(i,:)
                        
                    endif

                end do
                !----------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            ! For last step in one step iteration, all data moving forward with dt, but the data except for position moving start form last half step    
            else

                !----------------------------------------------------------------------------------------------
                ! SPH iteration 
                do i=1,actual_particle_number_in_subdomain

                    if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then

                        ! Renew density and energy (Full step start from last step) 
                        subdomain_particle_rho(i)=before_rho(i)+dt*drhodt(i)
                        subdomain_particle_energy(i)=before_energy(i)+dt*dedt(i)
                       
                        ! Renew velocity and position (Full step start from current half step) 
                        subdomain_particle_position(i,:)=subdomain_particle_position(i,:)+dt*( subdomain_particle_velocity(i,:)+average_velocity(i,:) )+0.5*particle_acceleration(i,:)*dt**2

                        subdomain_particle_velocity(i,:)=before_velocity(i,:)+dt*particle_acceleration(i,:)

                    endif

                enddo
                !----------------------------------------------------------------------------------------------

            end if
            !==================================================================================================

            !**************************************************************************************************

        elseif ( Iteration_Method==2 ) then
              
            !**************************************************************************************************
            ! Euler iteration

            ! For the first time step, we should run 2 SPH single step calculation

            !==================================================================================================
            ! Prediction step (All data moving forward with dt)
            ! Save last step information
            before_rho      = subdomain_particle_rho                           
            before_velocity = subdomain_particle_velocity                 
            before_position = subdomain_particle_position                 
            before_energy   = subdomain_particle_energy 
            
            ! For the first time step, we need calculate the SPH iteration information
            ! For the oter time step, we use the iteration information in last step as the prediction uodating information
            if ( i_time_step == actual_start_time_step ) then
                call Single_step_compute(i_time_step)                     ! Calculate SPH iteration information 
            endif

            ! Save Prediction step information
            temp_drhodt=drhodt
            temp_dedt=dedt
            temp_particle_acceleration=particle_acceleration 
            temp_particle_velocity=subdomain_particle_velocity

            !--------------------------------------------------------------------------------------------------
            do i=1,actual_particle_number_in_subdomain                 

                if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then
              
                    ! Renew density and energy 
                    subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt      
                    subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt  
                
                    ! Renew velocity and position
                    subdomain_particle_position(i,:)=before_position(i,:)+dt*(subdomain_particle_velocity(i,:)+average_velocity(i,:))+0.5*particle_acceleration(i,:)*dt**2
                    subdomain_particle_velocity(i,:)=before_velocity(i,:)+dt* particle_acceleration(i,:)

                endif

            enddo
            !--------------------------------------------------------------------------------------------------

            !==================================================================================================

            ! Exchange buffer domain information
            ! As the exchange will add errors, so it may cause unsteady for higher-order calculation
            ! We can ignore the changes from the Prediction step

            call Exchange_buffer_domain_information(i_time_step)

            !==================================================================================================
            ! Correction step
            call Single_step_compute(i_time_step)

            !--------------------------------------------------------------------------------------------------
            ! SPH iteration 

            do i=1,actual_particle_number_in_subdomain             

                if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then
                    
                    ! Renew density and energy 
                    subdomain_particle_rho(i)=before_rho(i)+0.5*(drhodt(i)+temp_drhodt(i))*dt   
                    subdomain_particle_energy(i)=before_energy(i)+0.5*(dedt(i)+temp_dedt(i))*dt  
                    
                    ! Renew velocity and position
                    subdomain_particle_position(i,:)=before_position(i,:)+0.5*dt*(subdomain_particle_velocity(i,:)+temp_particle_velocity(i,:))+average_velocity(i,:)*dt+0.5*particle_acceleration(i,:)*dt**2
                    subdomain_particle_velocity(i,:)=before_velocity(i,:)+0.5*dt*(particle_acceleration(i,:)+temp_particle_acceleration(i,:))

                endif

            end do
            !--------------------------------------------------------------------------------------------------

            !==================================================================================================

            !**************************************************************************************************


        elseif ( Iteration_Method==3 ) then


            !**************************************************************************************************
            ! Predictor-Corrector iteration (Reference: SPHysics theory 2nd-order)

            !==================================================================================================

            !--------------------------------------------------------------------------------------------------
            ! Prediction step (All data moving forward with dt)

            ! Save last step information
            before_rho=subdomain_particle_rho                           
            before_velocity=subdomain_particle_velocity                 
            before_position=subdomain_particle_position                 
            before_energy=subdomain_particle_energy 

            ! SPH calculation
            call Single_step_compute(i_time_step)

            !--------------------------------------------------------------------------------------------------
            ! Moving forward with 0.5dt
            do i=1,actual_particle_number_in_subdomain                 

                if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then
                   
                    ! Renew density and energy 
                    subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt*0.5d0      
                    subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt*0.5d0   
                
                    ! Renew velocity and position
                    subdomain_particle_position(i,:)=before_position(i,:)+0.5d0*dt*( subdomain_particle_velocity(i,:)+average_velocity(i,:) )+0.5*particle_acceleration(i,:)*(0.5d0*dt)**2
                    subdomain_particle_velocity(i,:)=before_velocity(i,:)+0.5d0*dt*particle_acceleration(i,:) 

                endif

            end do
            !--------------------------------------------------------------------------------------------------

            !==================================================================================================

            ! Exchange buffer domain information
            ! As the exchange will add errors, so it may cause unsteady for higher-order calculation
            ! We can ignore the changes from the Prediction step

            call Exchange_buffer_domain_information(i_time_step)

            ! Correction step
            call Single_step_compute(i_time_step)

            !==================================================================================================

            !--------------------------------------------------------------------------------------------------
            ! SPH iteration 
            do i=1,actual_particle_number_in_subdomain             

                if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then

                    !------------------------------------------------------------------------------------------
                    ! Moving forward with 0.5dt
                    subdomain_particle_rho(i)=before_rho(i)+0.5*drhodt(i)*dt   
                    subdomain_particle_energy(i)=before_energy(i)+0.5*dedt(i)*dt  
                  
                    subdomain_particle_position(i,:)=before_position(i,:)+0.5*dt*(subdomain_particle_velocity(i,:)+average_velocity(i,:))+0.5*particle_acceleration(i,:)*(0.5d0*dt)**2
                    subdomain_particle_velocity(i,:)=before_velocity(i,:)+0.5*dt*particle_acceleration(i,:)
                    !------------------------------------------------------------------------------------------

                    !------------------------------------------------------------------------------------------
                    !Final Results
                    subdomain_particle_rho(i)=2*subdomain_particle_rho(i)-before_rho(i)
                    subdomain_particle_energy(i)=2*subdomain_particle_energy(i)-before_energy(i) 
                  

                    subdomain_particle_velocity(i,:)=2*subdomain_particle_velocity(i,:)-before_velocity(i,:)
                    subdomain_particle_position(i,:)=2*subdomain_particle_position(i,:)-before_position(i,:)
                    !------------------------------------------------------------------------------------------

                endif

            end do
            !--------------------------------------------------------------------------------------------------

            !==================================================================================================

            !**************************************************************************************************

        elseif ( Iteration_Method==4 ) then 

            !**************************************************************************************************
            ! Beeman iteration: Beeman predictor step and Adams-Bashforth-Moulton corrector step. 
            ! (Reference: SPHysics theory 4th-order)

            !==================================================================================================
            ! Save last step information
            before_rho=subdomain_particle_rho                           
            before_velocity=subdomain_particle_velocity                 
            before_position=subdomain_particle_position                 
            before_energy=subdomain_particle_energy

            !--------------------------------------------------------------------------------------------------
            ! For first step moving with dt
            if ( i_time_step==actual_start_time_step ) then

                ! Prediction step (All data moving forward with dt)
                call Single_step_compute(i_time_step)

                !----------------------------------------------------------------------------------------------
                ! Moving forward with dt
                do i=1,actual_particle_number_in_subdomain                 

                    if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then
                  
                        subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt      
                        subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt  
                    
                        ! Renew velocity and position
                        subdomain_particle_position(i,:)=before_position(i,:)+dt*( subdomain_particle_velocity(i,:)+average_velocity(i,:) )+0.5*particle_acceleration(i,:)*dt**2
                        subdomain_particle_velocity(i,:)=before_velocity(i,:)+dt*particle_acceleration(i,:) 

                    endif

                enddo
                !----------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            !For other step moving with 1.5dt-0.5dt
            else

                !----------------------------------------------------------------------------------------------
                ! Save the information in last step
                K_1_temp_drhodt=drhodt
                K_1_temp_dedt=dedt
                K_1_temp_particle_acceleration=particle_acceleration
                K_1_temp_particle_velocity=subdomain_particle_velocity  

                ! Prediction step (All data moving forward with dt)
                call Single_step_compute(i_time_step)

                !----------------------------------------------------------------------------------------------


                !----------------------------------------------------------------------------------------------
                ! Moving forward with dt
                do i=1,actual_particle_number_in_subdomain                 

                    if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then
                  
                        subdomain_particle_rho(i)=before_rho(i)+(1.5*drhodt(i)-0.5*K_1_temp_drhodt(i))*dt       ! This K_1 is the accleration from last step   
                        subdomain_particle_energy(i)=before_energy(i)+(1.5*dedt(i)-0.5*K_1_temp_dedt(i))*dt  
                        
                        ! Renew velocity and position
                        subdomain_particle_velocity(i,:)=before_velocity(i,:)+( 1.5*particle_acceleration(i,:)-0.5*K_1_temp_particle_acceleration(i,:) )*dt
                        subdomain_particle_position(i,:)=before_position(i,:)+( 1.5*subdomain_particle_velocity(i,:)-0.5*K_1_temp_particle_velocity(i,:)+average_velocity(i,:) )*dt+2.0/3.0*particle_acceleration(i,:)*dt**2-1.0/6.0*K_1_temp_particle_acceleration(i,:)*dt**2

                    endif

                enddo
                !----------------------------------------------------------------------------------------------
                
                ! Save current step information V(n), a(n)
                K_2_temp_drhodt=drhodt
                K_2_temp_dedt=dedt
                K_2_temp_particle_acceleration=particle_acceleration
                K_2_temp_particle_velocity=subdomain_particle_velocity

                !==============================================================================================

                ! Exchange buffer domain information
                ! As the exchange will add errors, so it may cause unsteady for higher-order calculation
                ! We can ignore the changes from the Prediction step

                call Exchange_buffer_domain_information(i_time_step)

                ! Correction step
                call Single_step_compute(i_time_step)

                !==============================================================================================

                !----------------------------------------------------------------------------------------------
                ! SPH iteration 
                do i=1,actual_particle_number_in_subdomain             

                    if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then

                        !--------------------------------------------------------------------------------------
                        ! Moving forward with 0.5dt
                        subdomain_particle_rho(i)=before_rho(i)+( 5.0/12.0*drhodt(i)+8.0/12.0*K_2_temp_drhodt(i)-1.0/12.0*K_1_temp_drhodt(i) )*dt   
                        subdomain_particle_energy(i)=before_energy(i)+( 5.0/12.0*dedt(i)+8.0/12.0*K_2_temp_dedt(i)-1.0/12.0*K_1_temp_dedt(i) )*dt    
                      
                        subdomain_particle_velocity(i,:)=before_velocity(i,:)+( 5.0/12.0*particle_acceleration(i,:)+8.0/12.0*K_2_temp_particle_acceleration(i,:)-1.0/12.0*K_1_temp_particle_acceleration(i,:) )*dt  
                        subdomain_particle_position(i,:)=before_position(i,:)+( 5.0/12.0*subdomain_particle_velocity(i,:)+8.0/12.0*K_2_temp_particle_velocity(i,:)-1.0/12.0*K_1_temp_particle_velocity(i,:)+average_velocity(i,:) )*dt+1.0/6.0*particle_acceleration(i,:)*dt**2+1.0/3.0*K_2_temp_particle_acceleration(i,:)*dt**2
                        !--------------------------------------------------------------------------------------

                    endif
                    
                enddo
                !----------------------------------------------------------------------------------------------

                !==============================================================================================

            endif
            !--------------------------------------------------------------------------------------------------

            !**************************************************************************************************

        else

            !**************************************************************************************************
            ! Runge–Kutta iteration (4th-order)

            !==================================================================================================
            ! Save last step information (K_1)
            before_rho=subdomain_particle_rho                           
            before_velocity=subdomain_particle_velocity                 
            before_position=subdomain_particle_position                 
            before_energy=subdomain_particle_energy 

            ! Prediction step (Calculate K_1 and update postion with K_1 acceleration and 0.5dt for K_2)
            call Single_step_compute(i_time_step)

            !--------------------------------------------------------------------------------------------------
            ! SPH iteration 
            ! Save Prediction step information (Save K_1)
            K_1_temp_drhodt=drhodt
            K_1_temp_dedt=dedt
            K_1_temp_particle_acceleration=particle_acceleration
            K_1_temp_particle_velocity=subdomain_particle_velocity

            do i=1,actual_particle_number_in_subdomain                 

                if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then

                    subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt*0.5d0     
                    subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt*0.5d0  
                    
                    ! Renew velocity and position
                    subdomain_particle_position(i,:)=before_position(i,:)+0.5d0*dt*(subdomain_particle_velocity(i,:)+average_velocity(i,:))+0.5*0.25*particle_acceleration(i,:)*dt**2
                    subdomain_particle_velocity(i,:)=before_velocity(i,:)+0.5d0*dt*particle_acceleration(i,:)

                endif

            enddo
            !--------------------------------------------------------------------------------------------------

            !==================================================================================================

            ! Exchange buffer domain information
            ! As the exchange will add errors, so it may cause unsteady for higher-order calculation
            ! We can ignore the changes from the Prediction step

            call Exchange_buffer_domain_information(i_time_step)

            ! Prediction step (Calculate K_2 and update postion with K_2 acceleration and 0.5dt for K_3)
            call Single_step_compute(i_time_step)

            !==================================================================================================

            !--------------------------------------------------------------------------------------------------
            ! SPH iteration 
            ! Save Prediction step information (Save K_2)
            K_2_temp_particle_velocity=subdomain_particle_velocity
            K_2_temp_drhodt=drhodt
            K_2_temp_dedt=dedt
            K_2_temp_particle_acceleration=particle_acceleration

            do i=1,actual_particle_number_in_subdomain                 

                if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then
                
                    subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt/2      
                    subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt/2  
                    
                    ! Renew velocity and position
                    subdomain_particle_position(i,:)=before_position(i,:)+0.5*dt*(subdomain_particle_velocity(i,:)+average_velocity(i,:))+0.5*0.25*particle_acceleration(i,:)*dt**2
                    subdomain_particle_velocity(i,:)=before_velocity(i,:)+0.5*dt*particle_acceleration(i,:)

                endif

            enddo
            !--------------------------------------------------------------------------------------------------

            !==================================================================================================

            call Exchange_buffer_domain_information(i_time_step)

            ! Prediction step (Calculate K_3 and update postion with K_3 acceleration and dt for K_4)
            call Single_step_compute(i_time_step)

            !==================================================================================================


            !--------------------------------------------------------------------------------------------------
            ! SPH iteration 
            ! Save Prediction step information (Save K_3)
            K_3_temp_particle_velocity=subdomain_particle_velocity
            K_3_temp_drhodt=drhodt
            K_3_temp_dedt=dedt
            K_3_temp_particle_acceleration=particle_acceleration

            do i=1,actual_particle_number_in_subdomain                 

                if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then
                
                    subdomain_particle_rho(i)=before_rho(i)+drhodt(i)*dt      
                    subdomain_particle_energy(i)=before_energy(i)+dedt(i)*dt 
                    
                    ! Renew velocity and position
                    subdomain_particle_position(i,:)=before_position(i,:)+dt*(subdomain_particle_velocity(i,:)+average_velocity(i,:))+0.5*particle_acceleration(i,:)*dt**2
                    subdomain_particle_velocity(i,:)=before_velocity(i,:)+dt*particle_acceleration(i,:)

                endif

            enddo
            !--------------------------------------------------------------------------------------------------

            !==================================================================================================

            call Exchange_buffer_domain_information(i_time_step)

            ! Prediction step (Calculate K_4 and update final postion)
            call Single_step_compute(i_time_step)
          
            !==================================================================================================

            !--------------------------------------------------------------------------------------------------
            ! SPH iteration
            do i=1,actual_particle_number_in_subdomain             

                if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label ) then
                  
                    subdomain_particle_rho(i)=before_rho(i)+(K_1_temp_drhodt(i)+2*K_2_temp_drhodt(i)+2*K_3_temp_drhodt(i)+drhodt(i))*dt/6.0d0   
                    subdomain_particle_energy(i)=before_energy(i)+(K_1_temp_dedt(i)+2*K_2_temp_dedt(i)+2*K_3_temp_dedt(i)+dedt(i))*dt/6.0d0   
                      
                    subdomain_particle_position(i,:)=before_position(i,:)+(K_1_temp_particle_velocity(i,:)+2*K_2_temp_particle_velocity(i,:)+2*K_3_temp_particle_velocity(i,:)+subdomain_particle_velocity(i,:))*dt/6.0d0+dt*average_velocity(i,:)+0.5*particle_acceleration(i,:)*dt**2
                    subdomain_particle_velocity(i,:)=before_velocity(i,:)+(K_1_temp_particle_acceleration(i,:)+2*K_2_temp_particle_acceleration(i,:)+2*K_3_temp_particle_acceleration(i,:)+particle_acceleration(i,:))*dt/6.0d0

                endif

            enddo
            !--------------------------------------------------------------------------------------------------

            !==================================================================================================

            !**************************************************************************************************

        endif IterationMethodIf

        !======================================================================================================

















        !------------------------------------------------------------------------------------------------------
        ! Fliter the reults by using the Wavemaker damping ZONE
        ! Parts of the Variables are replaced by using the analytical results
        if (trim(adjustl(WavemakerDamping_Switch))=='on' .and. trim(adjustl(wave_maker_type))=='WaveTheory' .and. make_wave_or_not==1 ) then

            call Wave_Maker_Damping_Zone(i_time_step)

        endif
        !------------------------------------------------------------------------------------------------------
      
        !------------------------------------------------------------------------------------------------------ 
        ! Reinitialize particle density by using MLS
        if (trim(adjustl(RefreshDensityByMLS_Switch))=='on' .and. mod(i_time_step,RefreshDensity_time_step)==0) then
            call Reinitialize_particle_density(i_time_step)
        endif
        !------------------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------------------
        ! Shift particle position
        if (trim(adjustl(ShiftParticle_Switch))=='on') then
            call Shift_particle_position(i_time_step)
        endif
        !------------------------------------------------------------------------------------------------------
          
        !------------------------------------------------------------------------------------------------------
        ! Refresh the domain
        call Refresh_domain(i_time_step)
        !------------------------------------------------------------------------------------------------------




        !======================================================================================================
        ! Sampling data in current data

        ! Calculate the data for Post-Proceeding
        if(mod(i_time_step,save_time_step)==0 .and. trim(adjustl(OutputForPostProceedingGrid_Switch))=='on') then
            call Calculation_For_PostProceeding(i_time_step)
        endif

        ! Synchronize all processors calculation
        ! MPI_BCAST function has the character "Synchronize", so we don't need Synchronize
        ! call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
        !------------------------------------------------------------------------------------------------------


        ! !******************************************************************************************************
        ! !Save fluid information of subdomain for debugging
        ! if(mod(i_time_step,save_time_step)==0) then
         
        !     write(Current_Processor_File_Index,*) "ZONE I=",total_particle_number_in_subdomain," F=POINT"
        !     do i=1,total_particle_number_in_subdomain
        !        write(Current_Processor_File_Index,'(8F10.4)') (subdomain_particle_position(i,j),j=1,dim),(subdomain_particle_velocity(i,j),j=1,dim),subdomain_particle_rho(i),subdomain_particle_press(i)/pressure_0
        !     end do
          
        ! end if
        ! !******************************************************************************************************


        !======================================================================================================
        ! Saving Results
        if ( Current_Processor_ID==Main_Processor ) then

            !--------------------------------------------------------------------------------------------------
            ! Output running information on screen
            if ( Mod(i_time_step,Display_time_step)==0 ) then

                call Display_Running_Information_On_Screen()

            endif 
            !--------------------------------------------------------------------------------------------------


            !**************************************************************************************************
            ! Save fluid information
            if(mod(i_time_step,save_time_step)==0) then

                ! This SUBROUTINE is in 'SPH_SUBROUTINE_MODULE'
                call Output_particle_identification()

                !----------------------------------------------------------------------------------------------
                ! For tecplot
                if ( trim(adjustl(OutputForTecplot_Switch))=='on' ) then
                     call Output_for_tecplot(i_time_step)
                endif
                !----------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------
                ! Save the data for Paraview
                if ( trim(adjustl(OutputForParaview_Switch))=='on' ) then
                    call Output_for_paraview(i_time_step)
                endif
                !----------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------
                ! Save the grid data for Post-Proceeding
                if ( trim(adjustl(OutputForPostProceedingGrid_Switch))=='on' ) then
                    call Output_for_Grid_PostProceding(i_time_step)
                endif
                !----------------------------------------------------------------------------------------------

            endif
            !**************************************************************************************************
             
            !**************************************************************************************************
            ! Sampling
            if ( mod(i_time_step,sampling_time_step)==0 ) then
                call Sampling_subroutine(i_time_step)
            endif
            !**************************************************************************************************

            !**************************************************************************************************
            ! Save Results for Debugging
            if ( mod(i_time_step,save_time_step_for_debug)==0 .and. trim(adjustl(OutputForRestartSimulation_Switch))=='on' ) then

                call Output_for_Restart_Calculation(i_time_step)
     
            endif
            !**************************************************************************************************

        endif
        !======================================================================================================




    enddo TimeIterationDo


end subroutine SPH_Time_Iteration

