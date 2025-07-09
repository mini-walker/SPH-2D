!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: Wavemaker_Subroutine
!
!  PURPOSE: Assign wavemaker particles information
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

subroutine Wavemaker_Subroutine(i_time_step)

    use function_module
    use information_module
    use Public_variable_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables

    ! Variables from the mother subroutine
    integer,intent(in)::i_time_step                               !current time step

    !----------------------------------------------------------------------------------------------------------
    !wave Variables
    integer::i,j,k,l,m,o,v,r                                      !loop variables
    real(kind=8)::distance                                        !distance
    real(kind=8)::along_cline,vertical_cline,angle
    real(kind=8)::damping_coefficient                             !damping coefficient
    real(kind=8)::current_time                                    !current time
    real(kind=8)::body_g                                          !gravity acceleration

    real(kind=8)::current_wave_angle                              !current wave angle
    real(kind=8)::temp,temp_1,temp_2,temp_3,temp_4                !temp value
    real(kind=8),dimension(dim)::temp_analytical_particle_position!current particle position in the loop
    real(kind=8)::time_soft_coefficient                           !soft coefficient for time
    real(kind=8)::current_wave_elevation                          !wave elevation at current particle z direction
    real(kind=8)::m_1                                             !m value
    real(kind=8)::piston_stroke                                   !piston stroke
    real(kind=8),dimension(dim)::current_particle_position        !current particle position in the loop

    !----------------------------------------------------------------------------------------------------------
    !==========================================================================================================
    ! Body of subroutine Assign_boundary_particles

    ! Call MPI functions
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )

    !**********************************Assign the wave maker particles information******************************
    !Piston Type wavemaker for solitary waves or others
    WaveMakerType:if(trim(adjustl(wave_maker_type))=='PistonType') then

        !Define the time soft coefficient
        if (i_time_step<=wave_maker_soft_time+standing_time_step) then
            time_soft_coefficient=0.5*(sin((-0.5+(i_time_step-standing_time_step)/wave_maker_soft_time)*PI)+1)
        else
            time_soft_coefficient=1.0d0  
        endif

        body_g=g
        wave_maker_acceleration=0.0d0
        wave_maker_velocity=0.0d0
        wave_maker_position=0.0d0
        
        if(trim(adjustl(wave_type))=='Solitarywave') then

            !-------------------------------------------------------------------------
            !The solitary wave wavemaker position

            if (i_time_step>standing_time_step) then
                current_time=dt*(i_time_step-standing_time_step)
                call wavemaker_plate_motion(wave_maker_acceleration,wave_maker_velocity,wave_maker_position,current_time,wave_height,water_depth,dt,body_g)
            endif
            !-------------------------------------------------------------------------


        elseif (trim(adjustl(wave_type))=='Strokes_2nd') then

            !-------------------------------------------------------------------------
            !The Strokes_2nd wave wavemaker position
            !Reference: Long-crested wave generation and absorption for SPH-based DualSPHysics model
            if (i_time_step>standing_time_step) then
                
                current_time=dt*(i_time_step-standing_time_step)

                m_1=2*(sinh(wave_number*water_depth))**2/(sinh(wave_number*water_depth)*cosh(wave_number*water_depth)+wave_number*water_depth)

                current_wave_angle=wave_frequency*current_time+Wave_initial_phase_angle

                piston_stroke=wave_height/m_1

                wave_maker_position=piston_stroke/2.0*sin(current_wave_angle)+((wave_height**2/(32*water_depth))*(3*cosh(water_depth*wave_number)/(sinh(water_depth*wave_number))**3)-2.0/m_1)*sin(2*current_wave_angle)
                wave_maker_velocity=piston_stroke/2.0*wave_frequency*cos(current_wave_angle)+((wave_height**2/(32*water_depth))*(3*cosh(water_depth*wave_number)/(sinh(water_depth*wave_number))**3)-2.0/m_1)*2*wave_frequency*cos(2*current_wave_angle)
                wave_maker_acceleration=-piston_stroke/2.0*wave_frequency**2*sin(current_wave_angle)-((wave_height**2/(32*water_depth))*(3*cosh(water_depth*wave_number)/(sinh(water_depth*wave_number))**3)-2.0/m_1)*4*wave_frequency**2*sin(2*current_wave_angle)

                wave_maker_position=time_soft_coefficient*wave_maker_position
                wave_maker_velocity=time_soft_coefficient*wave_maker_velocity
                wave_maker_acceleration=time_soft_coefficient*wave_maker_acceleration

            endif
            !-------------------------------------------------------------------------

        elseif (trim(adjustl(wave_type))=='Airy') then

            !-------------------------------------------------------------------------
            !The linear wave wavemaker position
            !Reference: Long-crested wave generation and absorption for SPH-based DualSPHysics model
            if (i_time_step>standing_time_step) then
                current_time=dt*(i_time_step-standing_time_step)

                m_1=2*(sinh(wave_number*water_depth))**2/(sinh(wave_number*water_depth)*cosh(wave_number*water_depth)+wave_number*water_depth)

                current_wave_angle=wave_frequency*current_time+Wave_initial_phase_angle

                piston_stroke=wave_height/m_1

                wave_maker_position=piston_stroke/2.0*sin(current_wave_angle)
                wave_maker_velocity=piston_stroke/2.0*wave_frequency*cos(current_wave_angle)
                wave_maker_acceleration=-piston_stroke/2.0*wave_frequency**2*sin(current_wave_angle)

                wave_maker_position=time_soft_coefficient*wave_maker_position
                wave_maker_velocity=time_soft_coefficient*wave_maker_velocity
                wave_maker_acceleration=time_soft_coefficient*wave_maker_acceleration

            endif
            !-------------------------------------------------------------------------

        end if

        !-------------------------------------------------------------------------
        !检测造波端和水池末端波高
        left_wave_height=0.0d0
        right_wave_height=0.0d0

        do k=1,particle_ture_number
                                    
            !计算横坐标位置
            distance=abs(wave_maker_position-particle_position(k,1))
            
            !将位于其中的数据进行对比寻找最大值
            if(distance<=kernel_scale*interior_dx .and. particle_position(k,2)>left_wave_height ) then
                left_wave_height=particle_position(k,2)
            end if

            !计算点到直线的距离
            distance=abs(particle_position(k,2)-tan(theta_rad)*(particle_position(k,1)-incline_x))/sqrt(1+tan(theta_rad)**2)
            
            !将位于其中的数据进行对比寻找最大值
            if(distance<=kernel_scale*interior_dx .and. particle_position(k,2)>right_wave_height ) then
                right_wave_height=particle_position(k,2)
            end if

        end do

        left_wave_height=left_wave_height+0.5*interior_dx
        right_wave_height=right_wave_height+0.5*interior_dx
        !-------------------------------------------------------------------------
            
        !write(*,*) left_wave_height
        
        !-------------------------------------------------------------------------
        do j=1,wave_maker_particle_number
           
            i=j+particle_ture_number
            particle_position(i,1)=wave_maker_position-(wave_maker_particle_layer(i)-1)*interior_dx-interior_dx/2.0                      !虚粒子编号，横坐标
            particle_velocity(i,1)=wave_maker_velocity                                                                                                                  !虚粒子编号，x方向速度(初始时刻的位移、速度为0)
                
            if(particle_position(i,2)>left_wave_height) then
                
                particle_rho(i)=water_rho_0                                                                                              !密度赋值不能出现0
                particle_press(i)=0.0d0
                particle_mass(i)=0.0d0                                                       
            
            else

                particle_rho(i)=water_rho_0*(-6.0*g*(particle_position(i,2)-left_wave_height)/square_c_0+1.0)**(1.0/6.0)
                particle_press(i)=square_c_0*water_rho_0*((particle_rho(i)/water_rho_0)**7.0-1.0)/7.0
                particle_mass(i)=particle_rho(i)*intial_volume_0                             
                
            end if

        end do
        
        !write(*,*) wave_maker_position,wave_maker_velocity
        !-------------------------------------------------------------------------
        
        !-------------------------------------------------------------------------
        !将造波版左侧的粒子重新拉回流域内
        do i=1,particle_ture_number
            
            if(particle_type(i)==120) cycle        !排除病态粒子
            
            if(particle_position(i,1)<wave_maker_position) then
                particle_position(i,1)=2*wave_maker_position-particle_position(i,1)
            end if
                
            !加入数值项
            if((particle_position(i,1)-wave_maker_position)<interior_dx) then
                particle_velocity(i,1)=wave_maker_velocity
            end if

            ! if(particle_position(i,1)>9.5) then
               
            !     !计算点到直线的距离
            !     distance=abs(particle_position(i,2)-tan(theta_rad)*(particle_position(i,1)-10))/sqrt(1+tan(theta_rad)**2)    !点到右部倾斜板的距离
           
            !     !对于靠近边界粒子只保留切向速度
            !     if(distance<=0.2*smooth_length ) then
            !         along_cline=particle_velocity(i,1)*cos(theta_rad)+particle_velocity(i,2)*sin(theta_rad)
            !         vertical_cline=-particle_velocity(i,1)*sin(theta_rad)+particle_velocity(i,2)*cos(theta_rad)
        
            !         particle_velocity(i,1)=along_cline*cos(theta_rad)
            !         particle_velocity(i,2)=along_cline*sin(theta_rad)
        
            !     end if
                  
            ! end if
            
            ! !检测底部是否有漏
            ! if(particle_position(i,2)<=0.5*smooth_length ) then         !底部边界
            !    particle_velocity(i,2)=0.5*particle_velocity(i,2)
            ! end if
           
            ! !检测是否越界
            ! if(particle_position(i,2)<0) then         !底部边界拉回底部边界
            !     if(particle_position(i,2)>=-1.5*smooth_length) then
            !         particle_position(i,2)=-particle_position(i,2)
            !         particle_press(i)=particle_press(i)-2*9.81*1000*abs(particle_position(i,2))
            !     else
            !         particle_type(i)=120     
            !     end if   
                
            ! end if
            
            ! if(particle_position(i,1)>10 .and. particle_position(i,2)>0) then         !底部边界
            !     angle=atan2(particle_position(i,2),(particle_position(i,1)-10))
            !     if(angle<=theta_rad) then
            !         particle_type(i)=120
            !     end if
            ! end if
     
        end do
        !-------------------------------------------------------------------------

    !=================================================================================================================
    !Assign wavemaker particles directly basedon wave theory
    elseif(trim(adjustl(wave_maker_type))=='WaveTheory') then

        wave_maker_position=0.0d0

        !=============================================================================================================
        !Check the fluid particles in the domain or not
        do i=1,particle_ture_number
            
            if(particle_type(i)==120) cycle        !No error particles
            
            if(particle_position(i,2)<0.0d0) then
                particle_type(i)=120
            end if
                
        end do
        !=============================================================================================================

        !=============================================================================================================
        current_time=dt*(i_time_step-standing_time_step)

        !Define the time soft coefficient
        if (i_time_step<=wave_maker_soft_time+standing_time_step) then
            time_soft_coefficient=0.5*(sin((-0.5+(i_time_step-standing_time_step)/wave_maker_soft_time)*PI)+1)
        else
            time_soft_coefficient=1.0d0  
        endif

        !Define the right_wave_height for regular wave
        right_wave_height=water_depth

        !=============================================================================================================
        !For the time start wave making
        Skip_standing_time:if (current_time>=0.0d0) then

            Wavemaker_zone:if(trim(adjustl(wave_type))=='Stokes_2nd') then
            
                !----------------------------------------------------------------------
                !Wave_elevation=H/2*cos(kx-wt+Fai0)+(PI*H)/8*(H/L)*cosh(kh)[cos(2kh)+2]*cos(2(kx-wt+Fai0))/(sinh(kh))**3
                !Fai0: wave initial phase angle;
                !H: wave height;
                !k: wave number;
                !w: wave frequency;
                !t: time;
                !L: wave length;
                !h: water depth.

                !Define the wavemaker particles velocity 
                do j=1,wave_maker_particle_number
               
                    i=j+particle_ture_number

                    current_wave_angle=wave_number*particle_position(i,1)-wave_frequency*current_time+Wave_initial_phase_angle
                    
                    current_wave_elevation=time_soft_coefficient*(0.5*wave_height*cos(current_wave_angle)+PI*wave_height*wave_height/(8*wave_length)*cosh(wave_number*water_depth)/(sinh(wave_number*water_depth))**3*(2+cosh(2*wave_number*water_depth))*cos(2*current_wave_angle))
             
                    !Define the particle base on wave theory
                    particle_velocity(i,1)=PI*wave_height/wave_period*cosh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)+3.0/4.0*(PI**2)*wave_height/wave_period*(wave_height/wave_length)*cosh(2*wave_number*particle_position(i,dim))*cos(2*current_wave_angle)/(sinh(wave_number*water_depth))**4                  
                    particle_velocity(i,dim)=PI*wave_height/wave_period*sinh(wave_number*particle_position(i,dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)+3.0/4.0*(PI**2)*wave_height/wave_period*(wave_height/wave_length)*sinh(2*wave_number*particle_position(i,dim))*sin(2*current_wave_angle)/(sinh(wave_number*water_depth))**4
                    
                    !------------------------------------------------------------------
                    !Define the particle base on soft_coefficient
                    particle_velocity(i,1)=time_soft_coefficient*particle_velocity(i,1)
                    particle_velocity(i,dim)=time_soft_coefficient*particle_velocity(i,dim)

                    !Update the wave maker particles position
                    particle_position(i,1)=Initial_particle_Position(i,1)-time_soft_coefficient*(Wave_amplitue*cosh(wave_number*particle_position(i,dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)+3.0*PI*wave_height**2/(16*wave_length)*cosh(2*wave_number*particle_position(i,dim))*sin(2*current_wave_angle)/(sinh(wave_number*water_depth))**4)
                    particle_position(i,dim)=Initial_particle_Position(i,dim)+time_soft_coefficient*(Wave_amplitue*sinh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)+3.0*PI*wave_height**2/(16*wave_length)*sinh(2*wave_number*particle_position(i,dim))*cos(2*current_wave_angle)/(sinh(wave_number*water_depth))**4)

                    
                    if (particle_position(i,dim)<=current_wave_elevation+water_depth) then

                        !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                        particle_press(i)=water_rho_0*g*(current_wave_elevation+water_depth-particle_position(i,dim))+time_soft_coefficient*(water_rho_0*g*wave_height/2.0*cosh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/cosh(wave_number*water_depth)+3.0/4.0*water_rho_0*g*wave_height*(PI*wave_height/wave_length)*(1.0/sinh(2*wave_number*water_depth))*(cosh(2*wave_number*particle_position(i,dim))/(sinh(wave_number*water_depth))**2-1.0/3.0)*cos(2*current_wave_angle)-1.0/4.0*water_rho_0*g*wave_height*(PI*wave_height/wave_length)*cosh(2*wave_number*particle_position(i,dim))/sinh(2*wave_number*water_depth))
                        particle_rho(i)=water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                        particle_mass(i)=water_particle_mass_0
                        particle_c(i)=sqrt(square_c_0)

                    else
                       
                        particle_press(i)=0.0d0
                        particle_rho(i)=water_rho_0
                        particle_mass(i)=0.0d0
                        particle_c(i)=sqrt(square_c_0)                        
                        
                    endif
                    !------------------------------------------------------------------

                end do
                !----------------------------------------------------------------------

            elseif(trim(adjustl(wave_type))=='Airy') then
                
                !----------------------------------------------------------------------
                !Wave_elevation=H/2*cos(kx-wt+Fai0)
                
                !Define the wavemaker particles velocity 
                do j=1,wave_maker_particle_number
               
                    i=j+particle_ture_number

                    current_wave_angle=wave_number*particle_position(i,1)-wave_frequency*current_time+Wave_initial_phase_angle
                    current_wave_elevation=time_soft_coefficient*0.5*wave_height*cos(current_wave_angle)

                    !Define the particle base on wave theory
                    particle_velocity(i,1)=PI*wave_height/wave_period*cosh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)              
                    particle_velocity(i,dim)=PI*wave_height/wave_period*sinh(wave_number*particle_position(i,dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)
        
                    !-------------------------------------------------------------------
                    !Define the particle base on soft_coefficient
                    particle_velocity(i,1)=time_soft_coefficient*particle_velocity(i,1)
                    particle_velocity(i,dim)=time_soft_coefficient*particle_velocity(i,dim)

                    !Update the wave maker particles position
                    particle_position(i,1)=Initial_particle_Position(i,1)-time_soft_coefficient*Wave_amplitue*cosh(wave_number*particle_position(i,dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)
                    particle_position(i,dim)=Initial_particle_Position(i,dim)+time_soft_coefficient*Wave_amplitue*sinh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)

                    
                    if (particle_position(i,dim)<=current_wave_elevation+water_depth) then

                        !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                        particle_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-particle_position(i,dim))+time_soft_coefficient*(water_rho_0*g*wave_height/2.0*cosh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/cosh(wave_number*water_depth))
                        particle_rho(i)=water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                        particle_mass(i)=water_particle_mass_0
                        particle_c(i)=sqrt(square_c_0)

                    else
                       
                        particle_press(i)=0.0d0
                        particle_rho(i)=water_rho_0
                        particle_mass(i)=0.0d0
                        particle_c(i)=sqrt(square_c_0)                        
                        
                    endif
                    !-------------------------------------------------------------------

                end do
                !-----------------------------------------------------------------------

            elseif(trim(adjustl(wave_type))=='Stokes_5th') then
                
                !Define the wavemaker particles velocity 
                do j=1,wave_maker_particle_number
               
                    i=j+particle_ture_number

                    current_wave_angle=wave_number*particle_position(i,1)-wave_frequency*current_time+Wave_initial_phase_angle


                    !-------------------------------------------------------------------
                    ! Coefficients (Deep water wave is the initial condition)
                    ! "Hydrodynamics of Offshore Structures.P59"
                    current_wave_elevation=0.0d0
                    do m=1,5
                        current_wave_elevation=current_wave_elevation+Fifth_stokes_E(m)*cos(m*current_wave_angle)
                    enddo
                    current_wave_elevation=time_soft_coefficient*current_wave_elevation/wave_number

                    !Define the particle base on wave theory
                    particle_velocity(i,:)=0.0d0
                    temp_analytical_particle_position(:)=0.0d0
                    do m=1,5

                        particle_velocity(i,1)=particle_velocity(i,1)+m*Fifth_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)
                        particle_velocity(i,dim)=particle_velocity(i,dim)+m*Fifth_stokes_D(m)*sinh(m*wave_number*particle_position(i,dim))*sin(m*current_wave_angle)

                        temp_analytical_particle_position(1)=temp_analytical_particle_position(1)+Fifth_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*sin(m*current_wave_angle)
                        temp_analytical_particle_position(dim)=temp_analytical_particle_position(dim)+Fifth_stokes_D(m)*sinh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)
                    
                    enddo

                    particle_velocity(i,1)=time_soft_coefficient*Wave_celerity*particle_velocity(i,1)
                    particle_velocity(i,dim)=time_soft_coefficient*Wave_celerity*particle_velocity(i,dim)
                    
                    temp_analytical_particle_position(1)=time_soft_coefficient*(-Wave_celerity/wave_frequency)*temp_analytical_particle_position(1)
                    temp_analytical_particle_position(dim)=time_soft_coefficient*(Wave_celerity/wave_frequency)*temp_analytical_particle_position(dim)

                    !Update the wave maker particles position
                    particle_position(i,1)=Initial_particle_Position(i,1)+temp_analytical_particle_position(1)
                    particle_position(i,dim)=Initial_particle_Position(i,dim)+temp_analytical_particle_position(dim)

                    ! Stokes_DFai_Dt=0.0
                    ! do m=1,5

                    !     Stokes_DFai_Dt=Stokes_DFai_Dt+m*Fifth_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)

                    ! enddo
                    ! Stokes_DFai_Dt=-Wave_celerity**2*Stokes_DFai_Dt*time_soft_coefficient

                    if (particle_position(i,dim)<=current_wave_elevation+water_depth) then

                        !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                        !particle_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-particle_position(i,dim))-water_rho_0*0.5*(particle_velocity(i,1)**2+particle_velocity(i,dim)**2)+water_rho_0*Stokes_DFai_Dt
                        
                        particle_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-particle_position(i,dim))-water_rho_0*0.5*((particle_velocity(i,1)-time_soft_coefficient*Wave_celerity)**2+particle_velocity(i,dim)**2)-water_rho_0*time_soft_coefficient*Bernoulli_constant_R
                        
                        particle_rho(i)=water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                        particle_mass(i)=water_particle_mass_0
                        particle_c(i)=sqrt(square_c_0)

                    else
                       
                        particle_press(i)=0.0d0
                        particle_rho(i)=water_rho_0
                        particle_mass(i)=0.0d0
                        particle_c(i)=sqrt(square_c_0)                        
                        
                    endif
                    !-------------------------------------------------------------------

                end do
                !-----------------------------------------------------------------------

            elseif(trim(adjustl(wave_type))=='Stokes_3rd') then
                
                !Define the wavemaker particles velocity 
                do j=1,wave_maker_particle_number
               
                    i=j+particle_ture_number

                    current_wave_angle=wave_number*particle_position(i,1)-wave_frequency*current_time+Wave_initial_phase_angle

                    !-------------------------------------------------------------------
                    !Reference: Song Zhiyao "On the universal third order stokes wave solution"
                    current_wave_elevation=0.0d0

                    do m=1,3
                        current_wave_elevation=current_wave_elevation+Third_stokes_E(m)*cos(m*current_wave_angle)
                    enddo
                    current_wave_elevation=time_soft_coefficient*current_wave_elevation/wave_number

                    !Define the particle base on wave theory
                    particle_velocity(i,:)=0.0d0
                    temp_analytical_particle_position(:)=0.0d0
                    do m=1,3

                        particle_velocity(i,1)=particle_velocity(i,1)+m*Third_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)
                        particle_velocity(i,dim)=particle_velocity(i,dim)+m*Third_stokes_D(m)*sinh(m*wave_number*particle_position(i,dim))*sin(m*current_wave_angle)

                        temp_analytical_particle_position(1)=temp_analytical_particle_position(1)+Third_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*sin(m*current_wave_angle)
                        temp_analytical_particle_position(dim)=temp_analytical_particle_position(dim)+Third_stokes_D(m)*sinh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)
                    
                    enddo

                    particle_velocity(i,1)=time_soft_coefficient*Wave_celerity*particle_velocity(i,1)
                    particle_velocity(i,dim)=time_soft_coefficient*Wave_celerity*particle_velocity(i,dim)
                    
                    temp_analytical_particle_position(1)=time_soft_coefficient*(-Wave_celerity/wave_frequency)*temp_analytical_particle_position(1)
                    temp_analytical_particle_position(dim)=time_soft_coefficient*(Wave_celerity/wave_frequency)*temp_analytical_particle_position(dim)

                    !Update the wave maker particles position
                    particle_position(i,1)=Initial_particle_Position(i,1)+temp_analytical_particle_position(1)
                    particle_position(i,dim)=Initial_particle_Position(i,dim)+temp_analytical_particle_position(dim)

                    Stokes_DFai_Dt=0.0
                    do m=1,3

                        Stokes_DFai_Dt=Stokes_DFai_Dt+m*Third_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)

                    enddo
                    Stokes_DFai_Dt=-Wave_celerity**2*Stokes_DFai_Dt*time_soft_coefficient

                    if (particle_position(i,dim)<=current_wave_elevation+water_depth) then

                        !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                        particle_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-particle_position(i,dim))-water_rho_0*0.5*(particle_velocity(i,1)**2+particle_velocity(i,dim)**2)+water_rho_0*Stokes_DFai_Dt

                        particle_rho(i)=water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                        particle_mass(i)=water_particle_mass_0
                        particle_c(i)=sqrt(square_c_0)

                    else
                       
                        particle_press(i)=0.0d0
                        particle_rho(i)=water_rho_0
                        particle_mass(i)=0.0d0
                        particle_c(i)=sqrt(square_c_0)                        
                        
                    endif
                    !-------------------------------------------------------------------

                end do
                !-----------------------------------------------------------------------

            elseif(trim(adjustl(wave_type))=='Solitarywave') then
                
                !Define the wavemaker particles velocity 
                do j=1,wave_maker_particle_number
               
                    i=j+particle_ture_number

                    Solitary_s=1.0/cosh(Solitary_alpha*particle_position(i,1))
                    Solitary_t=tanh(Solitary_alpha*particle_position(i,1))
                    
                    current_wave_angle=wave_number*particle_position(i,1)-wave_frequency*current_time+Wave_initial_phase_angle

                    !-------------------------------------------------------------------
                    !Reference: John Fenton "A ninth-order solution for the solitary wave"
                    current_wave_elevation=1+Solitary_Epsilon*Solitary_s**2-3.0/4.0*Solitary_Epsilon**2*Solitary_s**2*Solitary_t**2+Solitary_Epsilon**3*(5.0/8.0*Solitary_s**2*Solitary_t**2-101/80*Solitary_s**4*Solitary_t**2)


                    !Define the particle base on wave theory
                    particle_velocity(i,1)=1+0.5*Solitary_Epsilon-3.0/20.0*Solitary_Epsilon**2+3/56.0*Solitary_Epsilon**3-Solitary_Epsilon*Solitary_s**2 &
                                            +Solitary_Epsilon**2*(-0.25*Solitary_s**2+Solitary_s**4+particle_position(i,dim)**2*(1.5*Solitary_s**2-2.25*Solitary_s**4))&
                                            +Solitary_Epsilon**3*(19.0/40.0*Solitary_s**2+0.2*Solitary_s**4-1.2*Solitary_s**6+particle_position(i,dim)**2*(-1.5*Solitary_s**2-3.75*Solitary_s**4+7.5*Solitary_s**6)+particle_position(i,dim)**4*(-0.375*Solitary_s**2+45.0/16.0*Solitary_s**4-45.0/16.0*Solitary_s**6))
                    
                    
                    particle_velocity(i,dim)=-Solitary_Epsilon*Solitary_s**2 &
                                             +Solitary_Epsilon**2*(0.375*Solitary_s**2+2*Solitary_s**4+particle_position(i,dim)**2*(0.5*Solitary_s**2-1.5*Solitary_s**4)) &
                                             +Solitary_Epsilon**3*(49/640.0*Solitary_s**2-17.0/20.0*Solitary_s**4-18/5.0*Solitary_s**6+particle_position(i,dim)**2*(-13.0/16.0*Solitary_s**2-25.0/16.0*Solitary_s**4+7.5*Solitary_s**6)+particle_position(i,dim)**4*(-3/40.0*Solitary_s**2+9.0/8.0*Solitary_s**4-27/16.0*Solitary_s**6))

                    
                    particle_velocity(i,1)=sqrt(g*water_depth)*particle_velocity(i,1)
                    particle_velocity(i,dim)=sqrt(3*Solitary_Epsilon)*particle_position(i,dim)*Solitary_t*sqrt(g*water_depth)*particle_velocity(i,dim)

                    if (particle_position(i,dim)<=current_wave_elevation+water_depth) then

                        !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                        particle_press(i)=water_rho_0*g*water_depth*(1-particle_position(i,dim)+Solitary_Epsilon*Solitary_s**2+Solitary_Epsilon**2*(0.75*Solitary_s**2-1.5*Solitary_s**4+particle_position(i,dim)**2*(-1.5*Solitary_s**2+2.25*Solitary_s**4))&
                                                                    +Solitary_Epsilon**3*(-0.5*Solitary_s**2-19.0/20.0*Solitary_s**4+11.0/5*Solitary_s**6+particle_position(i,dim)**2*(0.75*Solitary_s**2+39.0/8.0*Solitary_s**4-33.0/4.0*Solitary_s**6)+particle_position(i,dim)**4*(0.375*Solitary_s**2-45.0/16.0*Solitary_s**4+45/16.0*Solitary_s**6))&
                                                                     )

                        particle_rho(i)=water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                        particle_mass(i)=water_particle_mass_0
                        particle_c(i)=sqrt(square_c_0)

                    else
                       
                        particle_press(i)=0.0d0
                        particle_rho(i)=water_rho_0
                        particle_mass(i)=0.0d0
                        particle_c(i)=sqrt(square_c_0)                        
                        
                    endif

                    !Update the wave maker particles position
                    particle_position(i,1)=particle_position(i,1)+particle_velocity(i,1)*dt
                    particle_position(i,dim)=particle_position(i,dim)+particle_velocity(i,dim)*dt
                    !-------------------------------------------------------------------

                end do
                !-----------------------------------------------------------------------

            else
                
                !The main processor output the error information
                if (Current_Processor_ID==Main_Processor) then

                    write(*,*) "The wave type is not include in the code, please review the type."

                endif

            endif Wavemaker_zone
            !=============================================================================================================


            !=============================================================================================================
            !Update the particles in For the wavemaker damping zone
            !This is an independent data, not included in the SPH calculation
            WavemakerDamping:if (trim(adjustl(WavemakerDamping_Switch))=='on') then

                do i=1,particle_ture_number

                    if(particle_type(i)==120) cycle        !No error particles

                    if (In_WaveMakerDampingZone_OrNot(i)/=1) cycle                        !Only for the particles in Wavemaker Damping Zone
                        
                    !-----------------------------------------------------------------------------------------------------
                    !Wave_elevation=H/2*cos(kx-wt+Fai0)+(PI*H)/8*(H/L)*cosh(kh)[cos(2kh)+2]*cos(2(kx-wt+Fai0))/(sinh(kh))**3
                    if(trim(adjustl(wave_type))=='Stokes_2nd') then

                        current_particle_position(:)=WaveMaker_Analytical_position(i,:)

                        current_wave_angle=wave_number*current_particle_position(1)-wave_frequency*current_time+Wave_initial_phase_angle
                        current_wave_elevation=time_soft_coefficient*(0.5*wave_height*cos(current_wave_angle)+PI*wave_height*wave_height*(1.0/tanh(wave_number*water_depth))*(1+3.0/(2*(sinh(wave_number*water_depth))**3)*cos(2*current_wave_angle)))

                        !Define the particle base on wave theory
                        WaveMaker_Analytical_velocity(i,1)=PI*wave_height/wave_period*cosh(wave_number*current_particle_position(dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)-3.0/4.0*(PI**2)*wave_height/wave_period*(wave_height/wave_length)*cosh(2*wave_number*current_particle_position(dim))*cos(2*current_wave_angle)/(sinh(wave_number*water_depth))**4                  
                        WaveMaker_Analytical_velocity(i,dim)=PI*wave_height/wave_period*sinh(wave_number*current_particle_position(dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)-3.0/4.0*(PI**2)*wave_height/wave_period*(wave_height/wave_length)*sinh(2*wave_number*current_particle_position(dim))*sin(2*current_wave_angle)/(sinh(wave_number*water_depth))**4
                        
                        !------------------------------------------------------------------
                        !Define the particle base on soft_coefficient
                        WaveMaker_Analytical_velocity(i,1)=time_soft_coefficient*WaveMaker_Analytical_velocity(i,1)
                        WaveMaker_Analytical_velocity(i,dim)=time_soft_coefficient*WaveMaker_Analytical_velocity(i,dim)

                        !Update the wave maker particles position  
                        WaveMaker_Analytical_position(i,1)=Initial_particle_Position(i,1)-time_soft_coefficient*(Wave_amplitue*cosh(wave_number*current_particle_position(dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)+3.0*PI*wave_height**2/(16*wave_length)*cosh(2*wave_number*current_particle_position(dim))*sin(2*current_wave_angle)/(sinh(wave_number*water_depth))**4)
                        WaveMaker_Analytical_position(i,dim)=Initial_particle_Position(i,dim)+time_soft_coefficient*(Wave_amplitue*sinh(wave_number*current_particle_position(dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)+3.0*PI*wave_height**2/(16*wave_length)*sinh(2*wave_number*current_particle_position(dim))*cos(2*current_wave_angle)/(sinh(wave_number*water_depth))**4)

                        if (WaveMaker_Analytical_position(i,dim)<=current_wave_elevation+water_depth) then

                            !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                            WaveMaker_Analytical_press(i)=water_rho_0*g*(current_wave_elevation+water_depth-WaveMaker_Analytical_position(i,dim))+time_soft_coefficient*(water_rho_0*g*wave_height/2.0*cosh(wave_number*WaveMaker_Analytical_position(i,dim))*cos(current_wave_angle)/cosh(wave_number*water_depth)+3.0/4.0*water_rho_0*g*wave_height*(PI*wave_height/wave_length)*(1.0/sinh(2*wave_number*water_depth))*(cosh(2*wave_number*WaveMaker_Analytical_position(i,dim))/(sinh(wave_number*water_depth))**2-1.0/3.0)*cos(2*current_wave_angle)-1.0/4.0*water_rho_0*g*wave_height*(PI*wave_height/wave_length)*cosh(2*wave_number*WaveMaker_Analytical_position(i,dim))/sinh(2*wave_number*water_depth))
                            WaveMaker_Analytical_rho(i)=water_rho_0*(1+7*WaveMaker_Analytical_press(i)/(square_c_0*water_rho_0))**(0.142857142)

                        else
                           
                            WaveMaker_Analytical_press(i)=0.0d0
                            WaveMaker_Analytical_rho(i)=water_rho_0                      
                            
                        endif

                        !------------------------------------------------------------------

                    !-----------------------------------------------------------------------------------------------------
                    !Wave_elevation=H/2*cos(kx-wt+Fai0)
                    elseif(trim(adjustl(wave_type))=='Airy') then

                        !-------------------------------------------------------------------
                        current_particle_position(:)=WaveMaker_Analytical_position(i,:)

                        current_wave_angle=wave_number*current_particle_position(1)-wave_frequency*current_time+Wave_initial_phase_angle
                        current_wave_elevation=time_soft_coefficient*0.5*wave_height*cos(current_wave_angle)

                        !Define the particle base on wave theory
                        WaveMaker_Analytical_velocity(i,1)=PI*wave_height/wave_period*cosh(wave_number*current_particle_position(dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)              
                        WaveMaker_Analytical_velocity(i,dim)=PI*wave_height/wave_period*sinh(wave_number*current_particle_position(dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)
            
                        !-------------------------------------------------------------------
                        !Define the particle base on soft_coefficient
                        WaveMaker_Analytical_velocity(i,1)=time_soft_coefficient*WaveMaker_Analytical_velocity(i,1)
                        WaveMaker_Analytical_velocity(i,dim)=time_soft_coefficient*WaveMaker_Analytical_velocity(i,dim)

                        !Update the wave maker particles position
                        !position=initial_position+the ellpside
                        
                        WaveMaker_Analytical_position(i,1)=Initial_particle_Position(i,1)-time_soft_coefficient*Wave_amplitue*cosh(wave_number*current_particle_position(dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)              
                        WaveMaker_Analytical_position(i,dim)=Initial_particle_Position(i,dim)+time_soft_coefficient*Wave_amplitue*sinh(wave_number*current_particle_position(dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)              
                       

                        if (WaveMaker_Analytical_position(i,dim)<=current_wave_elevation+water_depth) then

                            !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                            WaveMaker_Analytical_press(i)=water_rho_0*g*(current_wave_elevation+water_depth-WaveMaker_Analytical_position(i,dim))+time_soft_coefficient*(water_rho_0*g*wave_height/2.0*cosh(wave_number*WaveMaker_Analytical_position(i,dim))*cos(current_wave_angle)/cosh(wave_number*water_depth))
                            WaveMaker_Analytical_rho(i)=water_rho_0*(1+7*WaveMaker_Analytical_press(i)/(square_c_0*water_rho_0))**(0.142857142)

                        else
                           
                            WaveMaker_Analytical_press(i)=0.0d0
                            WaveMaker_Analytical_rho(i)=water_rho_0                        
                            
                        endif
                        !-------------------------------------------------------------------

                    elseif(trim(adjustl(wave_type))=='Stokes_5th') then

                        !-------------------------------------------------------------------
                        current_particle_position(:)=WaveMaker_Analytical_position(i,:)

                        current_wave_angle=wave_number*current_particle_position(1)-wave_frequency*current_time+Wave_initial_phase_angle
                        

                        !-------------------------------------------------------------------
                        ! Coefficients (Deep water wave is the initial condition)
                        ! "Hydrodynamics of Offshore Structures.P59"
                        current_wave_elevation=0.0d0
                        do m=1,5
                            current_wave_elevation=current_wave_elevation+Fifth_stokes_E(m)*cos(m*current_wave_angle)
                        enddo
                        current_wave_elevation=time_soft_coefficient*current_wave_elevation/wave_number


                        !Define the particle base on wave theory
                        WaveMaker_Analytical_velocity(i,:)=0.0d0
                        temp_analytical_particle_position(:)=0.0d0
                        do m=1,5

                            WaveMaker_Analytical_velocity(i,1)=WaveMaker_Analytical_velocity(i,1)+m*Fifth_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)
                            WaveMaker_Analytical_velocity(i,dim)=WaveMaker_Analytical_velocity(i,dim)+m*Fifth_stokes_D(m)*sinh(m*wave_number*current_particle_position(dim))*sin(m*current_wave_angle)

                            temp_analytical_particle_position(1)=temp_analytical_particle_position(1)+Fifth_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*sin(m*current_wave_angle)
                            temp_analytical_particle_position(dim)=temp_analytical_particle_position(dim)+Fifth_stokes_D(m)*sinh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)
                        
                        enddo

                        WaveMaker_Analytical_velocity(i,1)=time_soft_coefficient*Wave_celerity*WaveMaker_Analytical_velocity(i,1)
                        WaveMaker_Analytical_velocity(i,dim)=time_soft_coefficient*Wave_celerity*WaveMaker_Analytical_velocity(i,dim)
                        
                        temp_analytical_particle_position(1)=time_soft_coefficient*(-Wave_celerity/wave_frequency)*temp_analytical_particle_position(1)
                        temp_analytical_particle_position(dim)=time_soft_coefficient*(Wave_celerity/wave_frequency)*temp_analytical_particle_position(dim)

                        !Update the wave maker particles position
                        WaveMaker_Analytical_position(i,1)=Initial_particle_Position(i,1)+temp_analytical_particle_position(1)
                        WaveMaker_Analytical_position(i,dim)=Initial_particle_Position(i,dim)+temp_analytical_particle_position(dim)


                        ! Stokes_DFai_Dt=0.0
                        ! do m=1,5

                        !     Stokes_DFai_Dt=Stokes_DFai_Dt+m*Fifth_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)

                        ! enddo
                        ! Stokes_DFai_Dt=-Wave_celerity**2*Stokes_DFai_Dt*time_soft_coefficient

                        if (WaveMaker_Analytical_position(i,dim)<=current_wave_elevation+water_depth) then

                            !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                            !WaveMaker_Analytical_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-current_particle_position(dim))-water_rho_0*0.5*(WaveMaker_Analytical_velocity(i,1)**2+WaveMaker_Analytical_velocity(i,dim)**2)+water_rho_0*Stokes_DFai_Dt
                            
                            WaveMaker_Analytical_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-WaveMaker_Analytical_position(i,dim))-water_rho_0*0.5*((WaveMaker_Analytical_velocity(i,1)-time_soft_coefficient*Wave_celerity)**2+WaveMaker_Analytical_velocity(i,dim)**2)-time_soft_coefficient*water_rho_0*Bernoulli_constant_R
                            
                            WaveMaker_Analytical_rho(i)=water_rho_0*(1+7*WaveMaker_Analytical_press(i)/(square_c_0*water_rho_0))**(0.142857142)

                        else
                           
                            WaveMaker_Analytical_press(i)=0.0d0
                            WaveMaker_Analytical_rho(i)=water_rho_0                        
                            
                        endif
                        !-------------------------------------------------------------------

                    elseif(trim(adjustl(wave_type))=='Stokes_3rd') then

                        !-------------------------------------------------------------------
                        current_particle_position(:)=WaveMaker_Analytical_position(i,:)

                        current_wave_angle=wave_number*current_particle_position(1)-wave_frequency*current_time+Wave_initial_phase_angle
                        
                        !-------------------------------------------------------------------
                        ! Coefficients (Deep water wave is the initial condition)
                        ! "Hydrodynamics of Offshore Structures.P59"
                        current_wave_elevation=0.0d0
                        do m=1,3
                            current_wave_elevation=current_wave_elevation+Third_stokes_E(m)*cos(m*current_wave_angle)
                        enddo
                        current_wave_elevation=time_soft_coefficient*current_wave_elevation/wave_number


                        !Define the particle base on wave theory
                        WaveMaker_Analytical_velocity(i,:)=0.0d0
                        temp_analytical_particle_position(:)=0.0d0
                        do m=1,3

                            WaveMaker_Analytical_velocity(i,1)=WaveMaker_Analytical_velocity(i,1)+m*Third_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)
                            WaveMaker_Analytical_velocity(i,dim)=WaveMaker_Analytical_velocity(i,dim)+m*Third_stokes_D(m)*sinh(m*wave_number*current_particle_position(dim))*sin(m*current_wave_angle)

                            temp_analytical_particle_position(1)=temp_analytical_particle_position(1)+Third_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*sin(m*current_wave_angle)
                            temp_analytical_particle_position(dim)=temp_analytical_particle_position(dim)+Third_stokes_D(m)*sinh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)
                        
                        enddo

                        WaveMaker_Analytical_velocity(i,1)=time_soft_coefficient*Wave_celerity*WaveMaker_Analytical_velocity(i,1)
                        WaveMaker_Analytical_velocity(i,dim)=time_soft_coefficient*Wave_celerity*WaveMaker_Analytical_velocity(i,dim)
                        
                        temp_analytical_particle_position(1)=time_soft_coefficient*(-Wave_celerity/wave_frequency)*temp_analytical_particle_position(1)
                        temp_analytical_particle_position(dim)=time_soft_coefficient*(Wave_celerity/wave_frequency)*temp_analytical_particle_position(dim)

                        !Update the wave maker particles position
                        WaveMaker_Analytical_position(i,1)=Initial_particle_Position(i,1)+temp_analytical_particle_position(1)
                        WaveMaker_Analytical_position(i,dim)=Initial_particle_Position(i,dim)+temp_analytical_particle_position(dim)

                        Stokes_DFai_Dt=0.0
                        do m=1,3

                            Stokes_DFai_Dt=Stokes_DFai_Dt+m*Third_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)

                        enddo
                        Stokes_DFai_Dt=-Wave_celerity**2*Stokes_DFai_Dt*time_soft_coefficient

                        if (WaveMaker_Analytical_position(i,dim)<=current_wave_elevation+water_depth) then

                            !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                            WaveMaker_Analytical_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-current_particle_position(dim))-water_rho_0*0.5*(WaveMaker_Analytical_velocity(i,1)**2+WaveMaker_Analytical_velocity(i,dim)**2)+water_rho_0*Stokes_DFai_Dt
                            
                            !WaveMaker_Analytical_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-WaveMaker_Analytical_position(i,dim))+water_rho_0*0.5*((WaveMaker_Analytical_velocity(i,1)-time_soft_coefficient*Wave_celerity)**2+WaveMaker_Analytical_velocity(i,dim)**2)-time_soft_coefficient*water_rho_0*Bernoulli_constant_R
                            
                            WaveMaker_Analytical_rho(i)=water_rho_0*(1+7*WaveMaker_Analytical_press(i)/(square_c_0*water_rho_0))**(0.142857142)

                        else
                           
                            WaveMaker_Analytical_press(i)=0.0d0
                            WaveMaker_Analytical_rho(i)=water_rho_0                        
                            
                        endif
                        !-------------------------------------------------------------------


                    !-------------------------------------------------------------------------------------------------
                    else
                        !The main processor output the error information
                        if (Current_Processor_ID==Main_Processor) then

                          write(*,*) "The wave type is not include in the code, please review the type.(wavemaker subroutine) "

                        endif

                    endif
                    !-------------------------------------------------------------------------------------------------

                end do

            endif WavemakerDamping
        !============================================================================================================
        
        else

        !============================================================================================================
        !For the standing time All the results for the hydrostatic 

            !----------------------------------------------------------------------
            do j=1,wave_maker_particle_number
           
                i=j+particle_ture_number

                !Define the particle base on wave theory
                particle_velocity(i,:)=0.0d0

                !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                particle_press(i)=water_rho_0*g*(water_depth-particle_position(i,dim))+water_rho_0*g*0.25*interior_dx  
                particle_rho(i)=water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                particle_mass(i)=water_particle_mass_0
                particle_c(i)=sqrt(square_c_0)

            end do
            !----------------------------------------------------------------------


            !----------------------------------------------------------------------
            do i=1,particle_ture_number

                if(particle_type(i)==120) cycle        !No error particles

                if (In_WaveMakerDampingZone_OrNot(i)/=1) cycle                        !Only for the particles in Wavemaker Damping Zone
                        
                !Define the particle base on wave theory
                WaveMaker_Analytical_velocity(i,:)=0.0d0

                !Update the wave maker particles position  
                WaveMaker_Analytical_position(i,:)=Initial_particle_Position(i,:)

                WaveMaker_Analytical_press(i)=water_rho_0*g*(water_depth-WaveMaker_Analytical_position(i,dim))
                WaveMaker_Analytical_rho(i)=water_rho_0*(1+7*WaveMaker_Analytical_press(i)/(square_c_0*water_rho_0))**(0.142857142)

            enddo
            !----------------------------------------------------------------------

        endif Skip_standing_time
        !=============================================================================================================
            

    !=================================================================================================================
    !The main processor output the error information
    else

       if (Current_Processor_ID==Main_Processor) then

          write(*,*) "The wave maker type is not include in the code, please review the type."

       endif

    endif WaveMakerType

!     !***********************************************************************************************************
!     !For Debug
!     !The main Processor Output the initial particle information
!     if (Current_Processor_ID==Main_Processor) then

!         !Open the Output file
!         open(unit=5,file="Check_Wavemaker.dat")    

!         !input tecplot header
!         write(5,*) "TITLE='DISTRIBUTION'"
!         write(5,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
        
!         !output fluid particle
!         write(5,*) "ZONE I=",particle_ture_number," F=POINT"
!         do i=1,particle_ture_number
!              write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
! 100          format(8F20.10) 
!         end do
        
!         !output wave maker particle
!         write(5,*) "ZONE I=",wave_maker_particle_number," F=POINT"
!         do k=1,wave_maker_particle_number
!             i=k+particle_ture_number
!             write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
!         end do
            
        
!         close(5)

!     endif

!     ! Synchronize all processors calculation
!     ! MPI_BCAST function has the character "Synchronize", so we don't need Synchronize
!     call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
!     !***********************************************************************************************************


end subroutine Wavemaker_Subroutine