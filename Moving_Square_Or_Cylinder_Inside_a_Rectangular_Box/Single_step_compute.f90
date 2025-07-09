!**************************************************************************************************************
!  PROGRAM   : Single_step_compute
!
!  PURPOSE   : Calculate the iteration information for single SPH step
!
!  Programer : Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location  : MUN
!
!  Time£     : 2017.3.18
!
!  Copyright : Memorial University
!
!  Version   : 1.0
!
!  Note      : MPI version: mpich-3.2
!              Fortran version: Inter Fortran (ifort) 18.0.0
!**************************************************************************************************************

subroutine Single_step_compute(i_time_step)              

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables form the mother subroutine
    integer,intent(in)::i_time_step                                     ! Current time step

    ! Variables in subroutine
    integer::i,j,k,L                                                    ! Loop Variables
    real(kind=8)::damping_coefficient                                   ! Damping coefficient
    real(kind=4)::current_time                                          ! Current time

    real(kind=8)::current_wave_angle                                    ! Current wave angle
    real(kind=8)::current_wave_elevation                                ! Wave elevation at current particle z direction

    ! Variables for file opeartion
    integer::File_index                                                 ! File index
    character(len=100)::File_name                                       ! File name
    character(len=4)::Char_Current_Processor_ID                         ! Current Processor ID in character

    !==========================================================================================================

    ! Body of subroutine Single_step_compute


    !**********************************************************************************************************
    ! Call the EOS(Equation of Status)
    ! Get the particle pressure
    do i=1,total_particle_number_in_subdomain
       
        if(subdomain_particle_type(i)==Ill_Particle_Label) cycle         ! No error particles

        call status_equation(subdomain_particle_rho(i),&                 ! Subdomain Particle Density(in)
                             subdomain_particle_energy(i),&              ! Subdomain Particle Energy(in)
                             subdomain_particle_type(i),&                ! Subdomain Particle Type(in)
                             subdomain_particle_press(i),&               ! Subdomain Particle Pressure(out)
                             subdomain_particle_c(i)&                    ! Subdomain Particle Sound Velocity(out)
                             ) 
    
    end do
    !**********************************************************************************************************


    !**********************************************************************************************************
    ! Assign boundary particles
    call Assign_boundary_particles(i_time_step)

    !**********************************************************************************************************


    !**********************************************************************************************************
    ! Initailized the Variables
    pair_i=0                                    ! Particle i index in kth pair (out) 
    pair_j=0                                    ! Particle j index in kth pair (out)
    w=0.0d0                                     ! Weight function value in kth pair(out)
    dwdx=0.0d0                                  ! Derivats of weight function in kth pair(out)
    pair_number=0                               ! Pair number(out)
    
    ! Get the particle pairs information
    call Find_particle_pair()                   
    
    !!Check w,dwdx
    !write(*,*) Current_Processor_ID,all_particle
    !do i=1,pair_number
    !    write(3,*) i,pair_i(i),pair_j(i),w(i),dwdx(i,1),dwdx(i,2)
    !end do
    !**********************************************************************************************************


    !**********************************************************************************************************
    ! Define the free surface and normalized the derivates of the kernel function
    Subdomain_fra_v=0.0d0
    Tensor_dyadic_matrix=0.0d0
    Tensor_dyadic_inverse_matrix=0.0d0
    
    call Define_freesurface_and_Normailized_derivates()

    !**********************************************************************************************************
    

    !**********************************************************************************************************
    ! Initailized the acceleration and energy change rates

    ! SPH iteration information
    average_velocity=0.0d0

    drhodt=0.0d0
    artifical_rho_correct=0.0d0

    dedt=0.0d0
    artifical_viscosity_dedt=0.0d0
    artifical_heat_dedt=0.0d0

    internal_acceleration=0.0d0  
    artifical_internal_acceleration=0.0d0     
    artifical_viscosity_acceleration=0.0d0 
    external_acceleration=0.0d0      
    laplacian_viscous_acceleration=0.0d0
    SPS_acceleration=0.0d0

    !-------------------------------------------------------------------------------------------------------
    !SPH iteration

    !call Solver_Standard_SPH()                                         !Standard SPH Solver
    !call Solver_Delta_SPH()                                            !Delta SPH Solver
    !call Solver_Modified_Delta_SPH()                                   !Modified Delta SPH Solver
    !call Solver_Standard_SPH_SPS()                                     !Standard SPH+SPS Solver
    !call Solver_Normalized_Standard_SPH()                              !Normalized Standard SPH Solver
    !call Solver_Normalized_Standard_SPH_SPS()                          !Normalized Standard SPH+SPS Solver
    !call Solver_Normalized_Delta_SPH_SPS()                             !Normalized Delta SPH+SPS Solver
    !call Solver_Normalized_Modified_Delta_SPH_SPS()                    !Normalized Delta SPH+SPS Solver

    !Modified Riemann is better for pressure
    !call Solver_Normalized_Modified_Linear_Riemann()                   !Modified Linear Riemann (Low dissipation) Normalized SPH 
    !call Solver_Normalized_Modified_Linear_Riemann_Delta_SPS()         !Modified Linear Riemann (Low dissipation) Normalized Delta SPH SPS Solver
    call Solver_Normalized_Modified_Linear_Riemann_Delta()              !Modified Linear Riemann (Low dissipation) Normalized Delta SPH Solver
    
    !call Solver_Normalized_Linear_Riemann_Delta_SPS()                  !Linear Riemann Normalized Delta SPH SPS Solver
    
    !HLLC is High dissipation, but more stable
    !call Solver_Normalized_HLLC_Riemann_Delta()                        !HLLC Riemann (High dissipation) Normalized Delta SPH Solver
    !call Solver_Normalized_HLLC_Riemann_Delta_SPS()                    !HLLC Riemann (High dissipation) Normalized Delta SPH SPS Solver
    
    !Higer order MUSCL approach for Riemann solver 
    !call Solver_Normalized_MUSCL_Linear_Riemann_Delta()                !MUSCL Linear Riemann (Low dissipation) Normalized Delta SPH Solver
    !call Solver_Normalized_MUSCL_HLLC_Riemann_Delta()                  !MUSCL HLLC Riemann (High dissipation) Normalized Delta SPH Solver
    
    !call Solver_Normalized_MUSCL_Upwind_Linear_Riemann_Delta()          !Upwind MUSCL HLLC Riemann (Low dissipation) Normalized Delta SPH Solver
    
    !call Solver_Normalized_MUSCL_Linear_Riemann_Delta_SPS()            !MUSCL Linear Riemann (Low dissipation) Normalized Delta SPH SPS Solver

    !**********************************************************************************************************



    !**********************************************************************************************************
    ! Non reflection buffer zone
    ! Reference: Nonreflecting outlet boundary conditions for incompressible flow using SPH

    if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then


        !------------------------------------------------------------------------------------------------------
        ! Re-assigned the inlet/outlet boundary particles
        do i=1,actual_particle_number_in_subdomain                
           
            ! As we update the in/outlet particle in the time iteration subroutine
            ! we should remove the iteration information of those particles
            if ( subdomain_particle_type(i) <= Inlet_Particle_Label ) then
                
                average_velocity(i,:)=0.0d0

                dedt(i)=0.0d0
                artifical_viscosity_dedt(i)=0.0d0
                artifical_heat_dedt(i)=0.0d0

                drhodt(i)=0.0d0
                artifical_rho_correct(i)=0.0d0

                internal_acceleration(i,:)=0.0d0       
                artifical_viscosity_acceleration(i,:)=0.0d0   
                external_acceleration(i,:)=0.0d0  
                artifical_internal_acceleration(i,:)=0.0d0  
                SPS_acceleration(i,:)=0.0d0  
                laplacian_viscous_acceleration(i,:)=0.0d0

            endif

        enddo
        !------------------------------------------------------------------------------------------------------

        ! Call the non-reflection method by Carlos
        ! And this method only works on the outlet particles
        if( Calos_Non_Reflection_Boundary==1 .and. 0.0d0 < Real_physical_time .and. Real_physical_time <= Shock_wave_release_time ) then

            call Carlos_Non_Reflection_Medthod()     ! Only Calculate the moment equation of one direction wave progration

        endif

    endif
    !**********************************************************************************************************



    !**********************************************************************************************************
    ! Calculate the force on the Body
    if ( Calculate_Body_Force_or_Not==1 .and. mod(i_time_step,sampling_time_step)==0 ) then

        ! Initailized the acceleration for body force calculation
        Body_Force_Pressure_Acceleration=0.0d0
        Body_Force_Viscous_Acceleration =0.0d0
        Body_Force_External_Acceleration=0.0d0

        call Calculate_Body_Force()                                               ! Calculate Body Force

        ! Assemble all the rates
        subdomain_Body_Pressure_Force=0.0d0
        subdomain_Body_Viscous_Force =0.0d0
        subdomain_Body_Total_Force   =0.0d0

        do i=1,actual_particle_number_in_subdomain                
           
            ! Only for the fluid particles near the body
            if(subdomain_Near_Body_Or_Not(i)==1) then               

                ! Renew the particles acceleration
                Body_Force_External_Acceleration(i,:)=Body_Force_External_Acceleration(i,:)+subdomain_particle_rho(i)*Gravity_grad(:)

                ! Calculate the force (As we calculete the force actel on the water, so the loading on body is the negative)
                subdomain_Body_Pressure_Force(:)=subdomain_Body_Pressure_Force(:)+subdomain_particle_volume(i)*Body_Force_Pressure_Acceleration(i,:)
                subdomain_Body_Viscous_Force(:)=subdomain_Body_Viscous_Force(:)+subdomain_particle_volume(i)*Body_Force_Viscous_Acceleration(i,:)
                subdomain_Body_Total_Force(:)=subdomain_Body_Total_Force(:)+subdomain_particle_volume(i)*( Body_Force_Pressure_Acceleration(i,:)+Body_Force_Viscous_Acceleration(i,:)+Body_Force_External_Acceleration(i,:) )
            
            endif
               
        enddo

    endif
    !**********************************************************************************************************




    !**********************************************************************************************************
    ! Assemble all the rates
    particle_acceleration=0.0d0

    do i=1,actual_particle_number_in_subdomain                
       
        ! Only for the fluid particles and inoutlet boundary particles
        if(subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)<=Inlet_Particle_Label) then               

            ! Define the gravity
            external_acceleration(i,:)=external_acceleration(i,:)+Gravity_grad(:)

            ! Accelerate the fluid with the ramping function
            ! V=0.5*U*( sin( (t/ramping_time-0.5)*PI ) + 1.0 )
            ! a=0.5*U*( cos( (t/ramping_time-0.5)*PI ) * PI/ramping_time )
            
            ! **Attention**: the ramping time is the real physical time, not the time step, it should be time_step*dt
            if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then
                external_acceleration(i,:)=external_acceleration(i,:)+Inflow_velocity(:)*InOutlet_Ramping_accleration_coefficient
            endif
            
            ! Renew the particles acceleration
            particle_acceleration(i,:)=internal_acceleration(i,:)+external_acceleration(i,:)+artifical_viscosity_acceleration(i,:)+artifical_internal_acceleration(i,:)+SPS_acceleration(i,:)+laplacian_viscous_acceleration(i,:) 
      
            ! Renew the particles energy change rate
            dedt(i)=dedt(i)+artifical_viscosity_dedt(i)+artifical_heat_dedt(i)

            ! Renew the particles density change rate
            drhodt(i)=drhodt(i)+artifical_rho_correct(i)

        endif
           
    enddo
    !**********************************************************************************************************











    !==========================================================================================================
    ! For the end damping zone (current velocity minus the damping coefficient times velocity directly)
    ! The value of Damping coefficient from 0 to 1 (From the start point to the end of the damping zone)
    ! Attention: we just need damping the velocity, the other parameters are not included in the damping zone
    EndDamping_if:if (trim(adjustl(EndDampingZone_Switch))=='on' .and. make_wave_or_not==1) then


        if (trim(adjustl(EndDampingZoneType))=='Cosin_Type') then

            !--------------------------------------------------------------------------------------------------
            !In Dr. Zheng Thesis
            do i=1,actual_particle_number_in_subdomain

               if(subdomain_particle_type(i)/=2) cycle                !Only for the fluid particles
              
               !In the damping zone or not
               if((domain_size_x-subdomain_particle_position(i,1))<=EndZoneDampingLengh .and. subdomain_particle_position(i,1)<=domain_size_x) then                  
                   
                   !calculate the damping coefficient
                   damping_coefficient=0.5*EndZoneDampingAmplitude*(1-cos(PI*(subdomain_particle_position(i,1)+EndZoneDampingLengh-domain_size_x)/EndZoneDampingLengh))
                           
                   !Damping the velocity
                   do j=1,dim
                       particle_acceleration(i,j)=particle_acceleration(i,j)*(1-damping_coefficient)
                       average_velocity(i,j)=average_velocity(i,j)*(1-damping_coefficient)
                       subdomain_particle_velocity(i,j)=subdomain_particle_velocity(i,j)*(1-damping_coefficient)
                   end do
               end if 
            end do
            !--------------------------------------------------------------------------------------------------

        elseif (trim(adjustl(EndDampingZoneType))=='OpenFoam_Type') then

            !--------------------------------------------------------------------------------------------------
            !In Damping zone in OpenFoam
            do i=1,actual_particle_number_in_subdomain

               if(subdomain_particle_type(i)/=2) cycle                !Only for the fluid particles
              
               !In the damping zone or not
               if((domain_size_x-subdomain_particle_position(i,1))<=EndZoneDampingLengh .and. subdomain_particle_position(i,1)<=domain_size_x) then         
                   
                   !calculate the damping coefficient
                   !value: (1-0)
                   damping_coefficient=(exp(((domain_size_x-subdomain_particle_position(i,1))/EndZoneDampingLengh)**3.5)-1.0)/(exp(1.0)-1.0)*(subdomain_particle_position(i,2)/water_depth)
                           
                   !Damping the velocity
                   do j=1,dim
                       ! particle_acceleration(i,j)=particle_acceleration(i,j)*(1-damping_coefficient)
                       ! average_velocity(i,j)=average_velocity(i,j)*(1-damping_coefficient)
                       ! subdomain_particle_velocity(i,j)=subdomain_particle_velocity(i,j)*(1-damping_coefficient)
                       particle_acceleration(i,j)=particle_acceleration(i,j)*damping_coefficient
                       average_velocity(i,j)=average_velocity(i,j)*damping_coefficient
                       subdomain_particle_velocity(i,j)=subdomain_particle_velocity(i,j)*damping_coefficient
                   end do
               end if 
            end do
            !--------------------------------------------------------------------------------------------------

        elseif (trim(adjustl(EndDampingZoneType))=='DualSPH_Type') then

            !--------------------------------------------------------------------------------------------------
            !Damping zone in S. Shao (ISPH wave simulation by using an internal wave maker)
            do i=1,actual_particle_number_in_subdomain

               if(subdomain_particle_type(i)/=2) cycle                !Only for the fluid particles
              
               !In the damping zone or not
               if((domain_size_x-subdomain_particle_position(i,1))<=EndZoneDampingLengh .and. subdomain_particle_position(i,1)<=domain_size_x) then                  
                   
                   !calculate the damping coefficient
                   !Value: 0-1
                   damping_coefficient=6*((subdomain_particle_position(i,1)+EndZoneDampingLengh-domain_size_x)/EndZoneDampingLengh)**2*((water_depth-subdomain_particle_position(i,2))/water_depth)
                           
                   !Damping the velocity
                   do j=1,dim
                       particle_acceleration(i,j)=particle_acceleration(i,j)*(1-damping_coefficient)
                       average_velocity(i,j)=average_velocity(i,j)*(1-damping_coefficient)
                       subdomain_particle_velocity(i,j)=subdomain_particle_velocity(i,j)*(1-damping_coefficient)
                   end do
               end if 
            end do
            !--------------------------------------------------------------------------------------------------
            
        endif
        
    endif EndDamping_if
    !==========================================================================================================





    
    
end subroutine Single_step_compute
