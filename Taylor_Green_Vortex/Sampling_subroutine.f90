!****************************************************************************
!
!  PROGRAM: Sampling_subroutine
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


!==================================================================================================================
subroutine Sampling(SamplingPosition,AveragePressure,MLSPressure,NumericalVelocity)
    
    use information_module
    use Public_variable_module
    use function_module

    implicit none

    !Variables from mother subroutine

    real(kind=8),dimension(dim),intent(in)::SamplingPosition       !current current Sampling Point Position
    real(kind=8),intent(out)::AveragePressure,MLSPressure          !Average and Normalized Pressure 
    real(kind=8),dimension(dim),intent(out)::NumericalVelocity     !Sampling Numerical Velocity

    !Variables in subroutine
    integer::i,j,k,L,m,v                                           !Loop Variables                                             
    
    !Grid information
    integer::x_n,y_n,mesh_name
    integer::near_mesh_number,near_mesh_name,particle_in_mesh,effect_particle_number
    integer,dimension(NumberInEachGridPrediction)::effect_particle
    real(kind=8)::distance                                         !distance  
    real(kind=8),dimension(dim)::position_difference               !position difference
    integer::Sampling_near_particle_number                         !Sampling near particle number
    integer,dimension(NumberInEachGridPrediction)::Sampling_near_particle_name
         
    !Weight function variables
    real(kind=8)::sum_w,temp_w
    real(kind=8),dimension(dim)::temp_dwdx  
    real(kind=8),dimension(NumberInEachGridPrediction)::W_kernel   !kernel function value
    real(kind=8),dimension(NumberInEachGridPrediction)::Fai        !weight Value           
        
    ! Body of subroutine Sampling Pressure subroutine
        
    !************************************************************************************
    ! Get the location of the probe based on the background grid
    x_n=floor((SamplingPosition(1)-chain_origin_x)/chain_dx+1)
    y_n=floor((SamplingPosition(2)-chain_origin_y)/chain_dy+1)

    !calculate the grid name and its near grids number
    mesh_name=(y_n-1)*chain_x_number+x_n                !grid name
    near_mesh_number=chain_near_mesh_number(mesh_name)  !near mesh number
    !************************************************************************************
     
     !***********************************************************************************
     !Get the effect particle in the near grids
     
     !initialized the effect near index and number
     effect_particle_number=0                                 
     effect_particle=0                                         
     
     do k=1,near_mesh_number                         
         
         near_mesh_name=chain_near_mesh_name(mesh_name,k)      !current near grid name
         particle_in_mesh=number(near_mesh_name)               !particle number in current grid
       
         do v=1,particle_in_mesh                               !search for each particle

             if(particle_in_chain(near_mesh_name,v)<=particle_ture_number) then
                 
                 effect_particle_number=effect_particle_number+1 
                 
                !make sure the Memory for the effect_particle_number is enough
                if(effect_particle_number<=NumberInEachGridPrediction) then
                    effect_particle(effect_particle_number)=particle_in_chain(near_mesh_name,v)
                else
                   write(*,*) "Memory for effect_particle is not engough, it should be larger than 300! (Sampling)"
                   !write(*,*) i,particle_position(i,1),particle_position(i,2)
                end if
            end if
           
        end do
     end do
     !***********************************************************************************
    
     !***********************************************************************************
     ! Only for the particle in the subdomain of the Sampling points
     Sampling_near_particle_number=0
     Sampling_near_particle_name=0
     
     do k=1,effect_particle_number
                
        i=effect_particle(k)

        position_difference(:)=particle_position(i,:)-SamplingPosition(:)
        distance=sqrt(DOT_PRODUCT(position_difference,position_difference))

        if (distance<=kernel_scale*smooth_length .and. distance>=1.0E-7) then                                                                                                                                
            Sampling_near_particle_number=Sampling_near_particle_number+1
            Sampling_near_particle_name(Sampling_near_particle_number)=i
        end if
           
     end do
             
    !initialized the Sampling Variables
    sum_w=0.0d0
    AveragePressure=0.0d0
    MLSPressure=0.0d0
    NumericalVelocity=0.0d0

    !calculate the pressure
    if(Sampling_near_particle_number>3) then
                 
        do L=1,Sampling_near_particle_number

            j=Sampling_near_particle_name(L)
                 
            !position difference
            position_difference(:)=SamplingPosition(:)-particle_position(j,:)
            distance=sqrt(DOT_PRODUCT(position_difference,position_difference))
        
            !write(*,*) distance,smooth_length,temp_w
            call compute_kernel(dim,&                                  !dimension of the simulation
                                distance,&                             !distance(in)                
                                position_difference,&                  !position_difference(in)
                                smooth_length,&                        !smooth length(in)
                                temp_w,&                               !kernel value(out)
                                temp_dwdx&                             !derivates of the kernel function(out)
                                )

            !sum of the weight function
            sum_w=sum_w+temp_w*particle_mass(j)/particle_rho(j)

            !kernel function value
            W_kernel(L)=temp_w
            
            !Sampling_Point_Average_Pressure
            AveragePressure=AveragePressure+particle_press(j)

        end do

        !Get the weight value of the effect particle
        call MLS_Interpolation_In_Gobal_Domain(MLS_order_Sampling,SamplingPosition,Sampling_near_particle_number,Sampling_near_particle_name,W_kernel,Fai)

        !Sampling Point MLS Pressure
        do L=1,Sampling_near_particle_number
    
            !Get the fluid particle index
            j=Sampling_near_particle_name(L)

            !interpolation
            if (trim(adjustl(Gravity_Switch))=='on') then
                MLSPressure=MLSPressure+particle_press(j)*Fai(L)-g*particle_rho(j)*(SamplingPosition(dim)-particle_position(j,dim))*Fai(L)
            else
                MLSPressure=MLSPressure+particle_press(j)*Fai(L)
            endif

            NumericalVelocity(:)=NumericalVelocity(:)+particle_velocity(j,:)*Fai(L)

        enddo
     
        !Sampling Point Normalized Pressure and velocity
        AveragePressure=AveragePressure/Sampling_near_particle_number

    endif

end subroutine Sampling
!==================================================================================================================










!==================================================================================================================
subroutine Sampling_subroutine(i_time_step)   
                            

    use information_module
    use Public_variable_module
    use function_module
   
    implicit none

    !Variables from mother subroutine

    integer,intent(in)::i_time_step                                !current time step
    
    !Variables in subroutine
    integer::i,j,k,L,m,v                                           !Loop Variables
    integer::current_time_number
    
    !File operation Variables
    integer::ioerror=0                                             !Open file return value
    integer::stat                                                  !Read file data return value
    integer::status                                                !Allocate memory status
    integer::file_index                                            !File index
    character(len=20)::char_file_index                             !Character file index
    character(len=50)::file_name                                   !File name                                        
    
    !Grid information
    real(kind=8)::distance                                         !Distance  
    real(kind=8),dimension(dim)::position_difference               !Position difference
         
    !Time
    real(kind=8)::current_time 
    character(len=20)::char_time                                   !Time character type             
        
    real(kind=8)::Q,R,x_0,current_K

    !Wave front
    real(kind=8)::wave_front
    real(kind=8),dimension(dim)::NumericalVelocity                 !Sampling Numerical Velocity
    real(kind=8)::velocity_magnitude,Maximum_velocity_magnitude    !velocity_magnitude
  

    ! Body of subroutine Sampling_subroutine
    
    !*********************************************************************************************
    !current time
    current_time=(i_time_step-standing_time_step)*dt                              

    !*********************************************************************************************
    !initialized the Variables for pressure Sampling
    Sampling_Point_Normalized_Pressure=0.0d0                             
    Sampling_Point_Average_Pressure=0.0d0                                
    Sampling_Point_Kalman_filter_Pressure=0.0d0                          
    Sampling_Position_near_particle_number=0
    Sampling_Position_near_particle_name=0
    
    !-----------------------------------------------------------------------------------------------------

    !************************************Get the wave elevation*****************************
    !initialized the wave height
    Wave_Probe_elevation=0.0d0

    !initialized the wave front data
    wave_front=0.0d0                          
                
    do i=1,particle_ture_number         !all the fluid particles
        
        !-----------------------------------------------------------------------------------
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
        !-----------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------
        !----wave front data----
        if(particle_position(i,1)>=wave_front ) then
            wave_front=particle_position(i,1)
        end if
        !-----------------------------------------------------------------------------------
    
    end do
    
    !save the wave information
    do i=1,Wave_Probe_number
        write(Wave_Probe_File_index(i),*) current_time*time_0,(Wave_Probe_elevation(i)+interior_dx/2.0)/domain_size_y
    end do

    !save the wave front
    write(127,*) current_time*time_0,wave_front/domain_size_y
    !***************************************************************************************

    !-----------------------------------------------------------------------------------------------------


    !-----------------------------------------------------------------------------------------------------

    !********************************Sampling the pressure data*****************************
    do m=1,Sampling_Point_number
    
        !Call Sampling pressure subroutine
        call Sampling(Sampling_Point_Position(m,:),Sampling_Point_Average_Pressure(m),Sampling_Point_Normalized_Pressure(m),NumericalVelocity)
        
        ! simplify kalman filter
        current_time_number=int(i_time_step/sampling_time_step)
                        
        ! kalman_filter(data,Q,R,x0,P0) 
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

        !***********************************************************************************
        
        !***********************************save the pressure data**************************
        write(Normalized_Pressure_Sampling_File_index(m),310) current_time*time_0,Sampling_Point_Normalized_Pressure(m)/pressure_0
310     format(20F20.8) 
        write(Average_Pressure_Sampling_File_index(m),310) current_time*time_0,Sampling_Point_Average_Pressure(m)/pressure_0
        write(Kalman_Pressure_Sampling_File_index(m),310) current_time*time_0,Sampling_Point_Kalman_filter_Pressure(m)/pressure_0
        !*************************************************************************************************
                 
    end do
    !-----------------------------------------------------------------------------------------------------


    !-----------------------------------------------------------------------------------------------------
    !Output the kinetic energy, Maximum velocity magnitude, norms
    Total_kinetic_energy=0.0d0
    Maximum_velocity_magnitude=0.0d0

    particle_velocity_error=0.0d0
    particle_press_error=0.0d0

    do i=1,particle_ture_number

       !Total kinetic energy
       Total_kinetic_energy=Total_kinetic_energy+0.5*particle_mass(i)*DOT_PRODUCT(particle_velocity(i,:),particle_velocity(i,:))
    
       !Maximum velocity magnitude
       velocity_magnitude=sqrt(DOT_PRODUCT(particle_velocity(i,:),particle_velocity(i,:)))
       if (velocity_magnitude>=Maximum_velocity_magnitude) then
           Maximum_velocity_magnitude=velocity_magnitude
       endif

       !errors
       particle_velocity_error(i,:)= (particle_velocity(i,:)-Analytical_particle_velocity(i,:))/velocity_0            !Particle velocity error 
       particle_press_error(i)= (particle_press(i)- Analytical_particle_press(i))/pressure_0                          !Particle pressure error
       particle_velocity_magnitude_error(i)=(velocity_magnitude-Analytical_particle_velocity_magnitude(i))/velocity_0 !Particle velocity magnitude error

    enddo

    write(128,310) Nondimensionalized_real_time,Total_kinetic_energy/Initial_total_kinetic_energy
    write(129,310) Nondimensionalized_real_time,Maximum_velocity_magnitude/velocity_0
    !-----------------------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------------------
    Norm_1=0.0d0
    Norm_2=0.0d0
    Norm_Infinity=0.0d0
    !Horizontal velocity Norm Result

    do i=1,particle_ture_number

        !Norm Infinity
        if (abs(particle_velocity_error(i,1))>=Norm_Infinity) then
            Norm_Infinity=abs(particle_velocity_error(i,1))
        endif

        Norm_1=Norm_1+abs(particle_velocity_error(i,1))
        Norm_2=Norm_2+(abs(particle_velocity_error(i,1)))**2

    enddo

    Norm_1=Norm_1/particle_ture_number

    Norm_2=sqrt(Norm_2/particle_ture_number)

    write(130,310) Nondimensionalized_real_time,Norm_Infinity,Norm_1,Norm_2

    
    Norm_1=0.0d0
    Norm_2=0.0d0
    Norm_Infinity=0.0d0
    !Vertical velocity Norm Result

    do i=1,particle_ture_number

        !Norm Infinity
        if (abs(particle_velocity_error(i,2))>=Norm_Infinity) then
            Norm_Infinity=abs(particle_velocity_error(i,2))
        endif

        Norm_1=Norm_1+abs(particle_velocity_error(i,2))
        Norm_2=Norm_2+(abs(particle_velocity_error(i,2)))**2

    enddo

    Norm_1=Norm_1/particle_ture_number

    Norm_2=sqrt(Norm_2/particle_ture_number)

    write(131,310) Nondimensionalized_real_time,Norm_Infinity,Norm_1,Norm_2


    Norm_1=0.0d0
    Norm_2=0.0d0
    Norm_Infinity=0.0d0
    !Pressure Norm Result

    do i=1,particle_ture_number

        !Norm Infinity
        if (abs(particle_press_error(i))>=Norm_Infinity) then
            Norm_Infinity=abs(particle_press_error(i))
        endif

        Norm_1=Norm_1+abs(particle_press_error(i))
        Norm_2=Norm_2+(abs(particle_press_error(i)))**2

    enddo

    Norm_1=Norm_1/particle_ture_number

    Norm_2=sqrt(Norm_2/particle_ture_number)

    write(132,310) Nondimensionalized_real_time,Norm_Infinity,Norm_1,Norm_2

    Norm_1=0.0d0
    Norm_2=0.0d0
    Norm_Infinity=0.0d0
    !Velocity magnitude norm Result

    do i=1,particle_ture_number

        !Norm Infinity
        if (abs(particle_velocity_magnitude_error(i))>=Norm_Infinity) then
            Norm_Infinity=abs(particle_velocity_magnitude_error(i))
        endif

        Norm_1=Norm_1+abs(particle_velocity_magnitude_error(i))
        Norm_2=Norm_2+(abs(particle_velocity_magnitude_error(i)))**2

    enddo

    Norm_1=Norm_1/particle_ture_number

    Norm_2=sqrt(Norm_2/particle_ture_number)

    write(133,310) Nondimensionalized_real_time,Norm_Infinity,Norm_1,Norm_2

    !-----------------------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------------------
    !Output the velocity and pressure at the desired position

    !-----------------------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------------------
    !Calculate the errors and norm
    do k=1,SamplingLineNumber
        
        do i=1,SamplingLinePointNumber

           !Call Sampling Pressure subroutine
           call Sampling(SamplingLinePointPosition(k,i,:),SamplingLinePointAveragePressure(k,i),SamplingLinePointMLSPressure(k,i),SamplingLinePointVelocity(k,i,:))
           
           !Calculate the analytical results

           SamplingLinePointAnalyticalVelocity(k,i,1)= U_Module*sin(2*PI*SamplingLinePointPosition(k,i,1))*cos(2*PI*SamplingLinePointPosition(k,i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time) 
           
           SamplingLinePointAnalyticalVelocity(k,i,2)=-U_Module*cos(2*PI*SamplingLinePointPosition(k,i,1))*sin(2*PI*SamplingLinePointPosition(k,i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time)  
            
           SamplingLinePointAnalyticalPressure(k,i)=0.25*water_rho_0*U_Module**2*(cos(4*PI*SamplingLinePointPosition(k,i,1))+cos(4*PI*SamplingLinePointPosition(k,i,2)))*exp(-16*PI**2/Reynolds_Number*Nondimensionalized_real_time) 

           !Calculate the errors and norm


        enddo

        !Output the results
        
        !transfer the sampling line index from integer to character
        write(char_file_index,'(I2)') k
           
        !transfer the Current time from real to character
        write(char_time,'(f7.4)') Nondimensionalized_real_time

        file_name="./Sampling/"//trim(adjustl(char_file_index))//"th_"//trim(adjustl(char_time))//"_SamplingLine_Results.dat"

        !Open the saving file
        open(unit=1,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    
        
        if(ioerror==0) then

            do i=1,SamplingLinePointNumber
                write(1,'(20f12.6)') (SamplingLinePointPosition(k,i,j),j=1,dim),(SamplingLinePointAnalyticalVelocity(k,i,j)/U_Module,j=1,dim),SamplingLinePointAnalyticalPressure(k,i)/pressure_0,(SamplingLinePointVelocity(k,i,j)/U_Module,j=1,dim),SamplingLinePointMLSPressure(k,i)/pressure_0,SamplingLinePointAveragePressure(k,i)/pressure_0
            end do

        end if
       
        close(1)

    enddo

    !-----------------------------------------------------------------------------------------------------

end subroutine Sampling_subroutine
!==================================================================================================================