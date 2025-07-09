!******************************************************************************************************************
!  PROGRAM: Sampling_subroutine
!
!  PURPOSE: Sampling the physical attributes at the sampling points
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
!******************************************************************************************************************

subroutine Sampling_subroutine(i_time_step)   
                            
    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
   
    implicit none

    !==============================================================================================================
    
    !--------------------------------------------------------------------------------------------------------------
    ! Variables from superior subroutine
    integer,intent(in)::i_time_step                                                           ! Current time step
    !--------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------
    ! Variables in subroutine
    integer::i,j,k,L,m,v                                                                      ! Loop Variables
    integer::current_time_number
    
    ! File operation Variables
    integer::File_index                                                                       ! File index
    character(len=20)::char_File_index                                                        ! Character file index
    character(len=50)::File_name                                                              ! File name                                        
    
         
    ! Time
    character(len=20)::char_time                                                              ! Time character type             
        
    real(kind=8)::Q,R,x_0,current_K

    ! Wave front
    real(kind=8)::wave_front
    real(kind=8),dimension(dim)::NumericalVelocity                                            ! Sampling Numerical Velocity
    real(kind=8)::velocity_magnitude,Maximum_velocity_magnitude                               ! velocity_magnitude
    !--------------------------------------------------------------------------------------------------------------

    !==============================================================================================================


    ! Body of subroutine Sampling_subroutine
    
 
    !--------------------------------------------------------------------------------------------------------------
    ! Initialized the Variables for pressure Sampling
    Sampling_Point_Normalized_Pressure=0.0d0                             
    Sampling_Point_Average_Pressure=0.0d0                                
    Sampling_Point_Kalman_filter_Pressure=0.0d0                          
    Sampling_Position_near_particle_number=0
    Sampling_Position_near_particle_name=0
    !--------------------------------------------------------------------------------------------------------------


    !*********************************************Get the wave elevation*******************************************
    ! For the case with free surface, you can get the free surface information
    if ( Has_FreeSurface_Or_Not==1 ) then

        ! Initialized the wave height
        Wave_Probe_elevation=0.0d0

        ! Initialized the wave front data
        wave_front=0.0d0                          
                    
        do i=1,particle_ture_number         ! All the fluid particles

            !----------------------------------------------------------------------------------------------------------
            !----wave height----
            !horizontial distance
            do j=1,Wave_Probe_number

                !Horizontial distance between wave probe and fluid particles
                distanceBetween_WP_and_P(j)=abs(Wave_Probe_Position(j,1)-particle_position(i,1))
                
                !For the particles close to the wave probe
                if(distanceBetween_WP_and_P(j)<=0.55*interior_dx .and. particle_position(i,2)>Wave_Probe_elevation(j) ) then
                    Wave_Probe_elevation(j)=particle_position(i,2)
                end if
                
            end do
            !----------------------------------------------------------------------------------------------------------

            !----------------------------------------------------------------------------------------------------------
            !----wave front data----
            if(particle_position(i,1)>=wave_front ) then
                wave_front=particle_position(i,1)
            end if
            !----------------------------------------------------------------------------------------------------------
        
        end do
        
        ! Save the wave information
        do i=1,Wave_Probe_number
            write(Wave_Probe_File_index(i),*) Nondimensionalized_real_time,(Wave_Probe_elevation(i)+interior_dx/2.0)/domain_size_y
        end do

        ! Save the wave front
        write(Wave_Front_File_Port,*) Nondimensionalized_real_time,wave_front/domain_size_y
        
    endif
    !**************************************************************************************************************



    !*******************************************Sampling the pressure data*****************************************
    ! Get the pressure at the pressure probe
    do m=1,Sampling_Point_number
    
        !----------------------------------------------------------------------------------------------------------
        ! Calculate the pressure
        ! Call Sampling pressure subroutine
        call Pressure_Sampling(Sampling_Point_Position(m,:),Sampling_Point_Average_Pressure(m),Sampling_Point_Normalized_Pressure(m),NumericalVelocity)
        
        ! Simplify kalman filter
        current_time_number=int(i_time_step/sampling_time_step)
                        
        ! Kalman_filter(data,Q,R,x0,P0) 
        Q=1.0e-6
        R=1.0e-5
        x_0=1.0d0
        
        if(current_time_number==1) then
            Sampling_Point_Kalman_filter_Pressure(m)=Sampling_Point_Normalized_Pressure(m)
            last_step_X(m)=x_0
            last_step_P(m)=Sampling_Point_Normalized_Pressure(m)
        else
           
            !kalman_filter
            current_K=last_step_P(m)/(last_step_P(m)+R)
            Sampling_Point_Kalman_filter_Pressure(m)=last_step_X(m)+current_K*(Sampling_Point_Normalized_Pressure(m)-last_step_X(m))
            
            ! save the result in the last kalman filter
            last_step_X(m)=Sampling_Point_Kalman_filter_Pressure(m)
            last_step_P(m)=last_step_P(m)-current_K*last_step_P(m)+Q
            
        end if
        !----------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------
        ! Save the pressure data
        write( Normalized_Pressure_Sampling_File_index(m),'(8F20.6)') Nondimensionalized_real_time,Sampling_Point_Normalized_Pressure(m)   /pressure_0
        write( Average_Pressure_Sampling_File_index(m)   ,'(8F20.6)') Nondimensionalized_real_time,Sampling_Point_Average_Pressure(m)      /pressure_0
        write( Kalman_Pressure_Sampling_File_index(m)    ,'(8F20.6)') Nondimensionalized_real_time,Sampling_Point_Kalman_filter_Pressure(m)/pressure_0
        !----------------------------------------------------------------------------------------------------------
                 
    end do
    !**************************************************************************************************************


    ! !**************************************************************************************************************
    ! ! Output the kinetic energy, Maximum velocity magnitude, norms for the Taylor-Green Vortex Case
    ! Total_kinetic_energy=0.0d0
    ! Maximum_velocity_magnitude=0.0d0

    ! particle_velocity_error=0.0d0
    ! particle_press_error=0.0d0

    ! !--------------------------------------------------------------------------------------------------------------
    ! do i=1,particle_ture_number

    !    ! Total kinetic energy
    !    Total_kinetic_energy=Total_kinetic_energy+0.5*particle_mass(i)*DOT_PRODUCT(particle_velocity(i,:),particle_velocity(i,:))
    
    !    ! Maximum velocity magnitude
    !    velocity_magnitude=sqrt(DOT_PRODUCT(particle_velocity(i,:),particle_velocity(i,:)))
    !    if (velocity_magnitude>=Maximum_velocity_magnitude) then
    !        Maximum_velocity_magnitude=velocity_magnitude
    !    endif

    !    ! Errors
    !    particle_velocity_error(i,:)= (particle_velocity(i,:)-Analytical_particle_velocity(i,:))/velocity_0            !Particle velocity error 
    !    particle_press_error(i)= (particle_press(i)- Analytical_particle_press(i))/pressure_0                          !Particle pressure error
    !    particle_velocity_magnitude_error(i)=(velocity_magnitude-Analytical_particle_velocity_magnitude(i))/velocity_0 !Particle velocity magnitude error

    ! enddo

    ! write(Kinetic_Energy_File_Port,            '(8F20.6)') Nondimensionalized_real_time,Total_kinetic_energy/Initial_total_kinetic_energy
    ! write(Maximum_velocity_magnitude_File_Port,'(8F20.6)') Nondimensionalized_real_time,Maximum_velocity_magnitude/velocity_0
    ! !--------------------------------------------------------------------------------------------------------------

    ! !--------------------------------------------------------------------------------------------------------------
    ! Norm_1=0.0d0
    ! Norm_2=0.0d0
    ! Norm_Infinity=0.0d0

    ! !Horizontal velocity Norm Result
    ! do i=1,particle_ture_number

    !     !Norm Infinity
    !     if (abs(particle_velocity_error(i,1))>=Norm_Infinity) then
    !         Norm_Infinity=abs(particle_velocity_error(i,1))
    !     endif

    !     Norm_1=Norm_1+abs(particle_velocity_error(i,1))
    !     Norm_2=Norm_2+(abs(particle_velocity_error(i,1)))**2

    ! enddo

    ! Norm_1=Norm_1/particle_ture_number

    ! Norm_2=sqrt(Norm_2/particle_ture_number)

    ! write(Norm_Horizontal_velocity_File_Port,'(8F20.6)') Nondimensionalized_real_time,Norm_Infinity,Norm_1,Norm_2
    ! !--------------------------------------------------------------------------------------------------------------
    

    ! !--------------------------------------------------------------------------------------------------------------
    ! Norm_1=0.0d0
    ! Norm_2=0.0d0
    ! Norm_Infinity=0.0d0

    ! !Vertical velocity Norm Result
    ! do i=1,particle_ture_number

    !     !Norm Infinity
    !     if (abs(particle_velocity_error(i,2))>=Norm_Infinity) then
    !         Norm_Infinity=abs(particle_velocity_error(i,2))
    !     endif

    !     Norm_1=Norm_1+abs(particle_velocity_error(i,2))
    !     Norm_2=Norm_2+(abs(particle_velocity_error(i,2)))**2

    ! enddo

    ! Norm_1=Norm_1/particle_ture_number

    ! Norm_2=sqrt(Norm_2/particle_ture_number)

    ! write(Norm_Vertical_velocity_File_Port,'(8F20.6)') Nondimensionalized_real_time,Norm_Infinity,Norm_1,Norm_2
    ! !--------------------------------------------------------------------------------------------------------------


    ! !--------------------------------------------------------------------------------------------------------------
    ! Norm_1=0.0d0
    ! Norm_2=0.0d0
    ! Norm_Infinity=0.0d0

    ! !Pressure Norm Result
    ! do i=1,particle_ture_number

    !     !Norm Infinity
    !     if (abs(particle_press_error(i))>=Norm_Infinity) then
    !         Norm_Infinity=abs(particle_press_error(i))
    !     endif

    !     Norm_1=Norm_1+abs(particle_press_error(i))
    !     Norm_2=Norm_2+(abs(particle_press_error(i)))**2

    ! enddo

    ! Norm_1=Norm_1/particle_ture_number

    ! Norm_2=sqrt(Norm_2/particle_ture_number)

    ! write(Norm_Pressure_File_Port,'(8F20.6)') Nondimensionalized_real_time,Norm_Infinity,Norm_1,Norm_2
    ! !--------------------------------------------------------------------------------------------------------------


    ! !--------------------------------------------------------------------------------------------------------------
    ! Norm_1=0.0d0
    ! Norm_2=0.0d0
    ! Norm_Infinity=0.0d0
    ! !Velocity magnitude norm Result

    ! do i=1,particle_ture_number

    !     !Norm Infinity
    !     if (abs(particle_velocity_magnitude_error(i))>=Norm_Infinity) then
    !         Norm_Infinity=abs(particle_velocity_magnitude_error(i))
    !     endif

    !     Norm_1=Norm_1+abs(particle_velocity_magnitude_error(i))
    !     Norm_2=Norm_2+(abs(particle_velocity_magnitude_error(i)))**2

    ! enddo

    ! Norm_1=Norm_1/particle_ture_number

    ! Norm_2=sqrt(Norm_2/particle_ture_number)

    ! write(Norm_Velocity_magnitude_File_Port,'(8F20.6)') Nondimensionalized_real_time,Norm_Infinity,Norm_1,Norm_2
    ! !--------------------------------------------------------------------------------------------------------------

    ! !**************************************************************************************************************




    !**************************************************************************************************************
    !Output the velocity and pressure at the desired line

    !--------------------------------------------------------------------------------------------------------------
    !Calculate the errors and norm
    do k=1,SamplingLineNumber
        
        do i=1,SamplingLinePointNumber

           !Call Sampling Pressure subroutine
           call Pressure_Sampling(SamplingLinePointPosition(k,i,:),SamplingLinePointAveragePressure(k,i),SamplingLinePointMLSPressure(k,i),SamplingLinePointVelocity(k,i,:))
           
           ! !Calculate the analytical results of Taylor-Green Vortex Case

           ! SamplingLinePointAnalyticalVelocity(k,i,1)= U_Module*sin(2*PI*SamplingLinePointPosition(k,i,1))*cos(2*PI*SamplingLinePointPosition(k,i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time) 
           
           ! SamplingLinePointAnalyticalVelocity(k,i,2)=-U_Module*cos(2*PI*SamplingLinePointPosition(k,i,1))*sin(2*PI*SamplingLinePointPosition(k,i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time)  
            
           ! SamplingLinePointAnalyticalPressure(k,i)=0.25*water_rho_0*U_Module**2*(cos(4*PI*SamplingLinePointPosition(k,i,1))+cos(4*PI*SamplingLinePointPosition(k,i,2)))*exp(-16*PI**2/Reynolds_Number*Nondimensionalized_real_time) 

           ! !Calculate the errors and norm


        enddo

        !Output the results
        
        !transfer the sampling line index from integer to character
        write(char_File_index,'(I2)') k
           
        !transfer the Current time from real to character
        write(char_time,'(f7.4)') Nondimensionalized_real_time

        File_index = General_File_Port
        File_name  = "./Sampling/"//trim(adjustl(char_File_index))//"th_"//trim(adjustl(char_time))//"_SamplingLine_Results.dat"

        !Open the saving file
        open(unit=General_File_Port,file=trim(adjustl(File_name)),status="replace",position="rewind",action="write",iostat=ioerror)    
        
        if(ioerror==0) then

            do i=1,SamplingLinePointNumber
                write(General_File_Port,'(20f12.6)') (SamplingLinePointPosition(k,i,j),j=1,dim),(SamplingLinePointVelocity(k,i,j)/U_Module,j=1,dim),(SamplingLinePointMLSPressure(k,i)-Background_Pressure)/pressure_0,(SamplingLinePointAveragePressure(k,i)-Background_Pressure)/pressure_0
            end do

        end if
       
        close(General_File_Port)

    enddo
    !--------------------------------------------------------------------------------------------------------------


    !--------------------------------------------------------------------------------------------------------------
    !Output the body force
    if ( Calculate_Body_Force_or_Not==1 ) then

        write(Body_Force_File_Port,'(20f12.6)') Nondimensionalized_real_time,(Body_Pressure_Force(i)/force_0,i=1,dim),(Body_Viscous_Force(i)/force_0,i=1,dim),(Body_Total_Force(i)/force_0,i=1,dim)

    endif
    !--------------------------------------------------------------------------------------------------------------


    !**************************************************************************************************************


end subroutine Sampling_subroutine