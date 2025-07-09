!**************************************************************************************************************
!  SUBROUTINE: Assign_boundary_particles
!
!  PURPOSE   : Assign boundary particles information
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

subroutine Assign_boundary_particles(i_time_step)

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables

    ! Variables from the mother subroutine
    integer,intent(in)::i_time_step                                                  ! Current time step
    !-----------------------------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------------------------
    ! Variables in subroutine
    integer::i,j,k,l,o,v,r                                                           ! Variables for looping
    integer::fix_ghost_name                                                          ! Fix ghost particle index
    integer::mesh_name                                                               ! Mesh index
    integer::Ture_effect_number                                                      ! Ture effect particle number
    integer,dimension(PredictionNumberInSupportDomain)::Ture_effect_name             ! Ture effect particle index

    integer::near_mesh_number                                                        ! Near mesh number
    integer::near_mesh_name                                                          ! Near mesh index
    integer::particle_in_mesh                                                        ! Particle number in current mesh
    real(kind=8),dimension(dim)::Temp_Velocity                                       ! Temp velocity
    real(kind=8),dimension(dim)::Mirror_Particle_Position                            ! Mirror Particle Position
    real(kind=8),dimension(dim)::Gradience                                           ! Gradience
    integer::x_n,y_n,z_n                                                             ! Grid index in x, y and z
    integer,dimension(dim)::Grid_Position_Index                                      ! The position index of grid which point in 

    ! Variables for Interpolation
    real(kind=8),dimension(PredictionNumberInSupportDomain)::W_kernel                ! Weight kernel function value
    real(kind=8),dimension(PredictionNumberInSupportDomain,dim)::DwDx_kernel         ! Weight kernel function derivate value
    real(kind=8),dimension(PredictionNumberInSupportDomain)::Fai                     ! Output MLS weight function value

    real(kind=8),dimension(dim)::Mirror_Particle_Interpolation_Velocity              ! Interpolation velocity of the mirror particle
    real(kind=8)::Mirror_Particle_Interpolation_Pressure                             ! Interpolation pressure of the mirror particle


    ! Variables for non-reflection boundary (Martin Lastiwka, Mihai Basa paper)
    real(kind=8)::rho_reference,pressure_reference,u_reference                       ! Reference value of the inlet
    real(kind=8)::Temp_c
    real(kind=8)::J_1,J_2,J_3

    ! Variables for Compatibility Interpolation
    real(kind=8),dimension(Compatibility_Equation_dim,Compatibility_Equation_dim)::Compatibility_Coefficient_Matrix
    real(kind=8),dimension(Compatibility_Equation_dim,4)::Compatibility_Result_Matrix! Right hand side matrix

    integer::File_index                                                              ! Output file index
    character(len=100)::File_name                                                    ! Output file name
    character(len=4)::Char_Current_Processor_ID                                      ! Character of the pressure processor ID

    integer::check_error                                                             ! Interpolation successful---1; Failed---0 

    !----------------------------------------------------------------------------------------------------------

    !==========================================================================================================

    ! Body of subroutine Assign_boundary_particles

    !----------------------------------------------------------------------------------------------------------
    ! Initialized the variables
    subdomain_number=0
    subdomain_particle_in_chain=0
    subdomain_particle_chain_number=0

    ! Project all the particles (Total_particle_number_in_subdomain) to the background grid
    do i=1,Total_particle_number_in_subdomain
        
        if( subdomain_particle_type(i)==Ill_Particle_Label .or. subdomain_particle_type(i)==0 ) cycle         ! No error particles

        ! Calculate the coordinates of the particles in the list

        ! Input value: particle_position
        ! Return value: grid position index, grid name
        call Point_Project_Into_Subdomain_Background_Grid( subdomain_particle_position(i,:),Grid_Position_Index,k )

        x_n = Grid_Position_Index(1);
        y_n = Grid_Position_Index(2);

        if ( x_n>=1 .and. y_n>=1 .and. x_n<=subdomain_chain_x_number .and. y_n<=subdomain_chain_y_number) then
            
            subdomain_particle_chain_number(i)=k                               ! Save the particle background grid index

            !write(*,*) k,chain_max_number

            if( subdomain_number(k) < PredictionNumberInSingleGrid ) then
                subdomain_number(k)=subdomain_number(k)+1                      ! Particle number in current grid index plus 1
                subdomain_particle_in_chain(k,subdomain_number(k))=i           ! The particle index in kth grid is i
            else
                write(*,'(A)') " The subdomain chain's containablity is not enough! (Assign)"
            endif

        else

            subdomain_particle_type(i) = Ill_Particle_Label

        endif

        ! if (i==1) then
        !     write(*,*) i,x_n,y_n,k,subdomain_chain_x_number,subdomain_chain_y_number
        ! endif

        ! if (i==Actual_particle_number_in_subdomain+1) then
        !     write(*,*) i,x_n,y_n,k,subdomain_chain_x_number,subdomain_chain_y_number
        ! endif

    enddo
    !----------------------------------------------------------------------------------------------------------


    !----------------------------------------------------------------------------------------------------------
    ! Find all the vicinity particles of each particle for total particles(actual+buffer zone)
    
    ! Initialized the subdomain effect particle and number
    subdomain_effect_particle=0
    subdomain_effect_particle_number=0

    do i=1,Total_particle_number_in_subdomain

        ! The background grid index of current particle
        mesh_name=subdomain_particle_chain_number(i)                  ! Mesh index
        near_mesh_number=subdomain_chain_near_mesh_number(mesh_name)  ! Vicinity grid index of background grid index of current particle
        
        do k=1,near_mesh_number                         
           
           near_mesh_name=subdomain_chain_near_mesh_name(mesh_name,k) ! Current background grid index
           particle_in_mesh=subdomain_number(near_mesh_name)          ! Particle number in current background grid index
           
           do v=1,particle_in_mesh                                    ! Search all the particles in current background grid
    
                j=subdomain_particle_in_chain(near_mesh_name,v)
               
                ! Position difference and distance i and j particle
                position_difference(:)=subdomain_particle_position(i,:)-subdomain_particle_position(j,:)
                distance=sqrt(DOT_PRODUCT(position_difference,position_difference))
                
                ! For vicinity particle searching, it should be a little large)
                average_smooth_length = 1.3*smooth_length
                
                ! In the support domain or not
                if( distance <= kernel_scale*average_smooth_length ) then

                    ! Save the particle index
                    if( subdomain_effect_particle_number(i) < PredictionNumberInAllVicinityGrids ) then

                        subdomain_effect_particle_number(i)=subdomain_effect_particle_number(i)+1 
                        subdomain_effect_particle(i,subdomain_effect_particle_number(i)) = j

                    else

                        write(*,'(A)') " 'PredictionNumberInAllVicinityGrids' is not engough, it should be larger than 2000! (Assign)"
                        
                    endif

                endif
               
           enddo

        enddo

    enddo
    !----------------------------------------------------------------------------------------------------------


    !**********************************************************************************************************




    !********************************Assign the Outlet boundary particles velocity*****************************

    subdomain_Particle_Interpolation_Velocity=0.0d0
    subdomain_Particle_Interpolation_Pressure=0.0d0
    subdomain_Particle_Interpolation_Rho=0.0d0

    subdomain_in_non_reflection_zone_or_not=0

    ! For the inlet (***Very important***)
    ! The velocity are defined by the inlet definition; 
    ! And the pressure and velocity are exterpolation from the fluid domain;
    if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then

        !======================================================================================================
        do i=1,Actual_particle_number_in_subdomain

            !--------------------------------------------------------------------------------------------------
            ! Interpolate the density, pressure and velocity at the inlet/outlet boundary particles
            if( subdomain_particle_type(i)<=Inlet_Particle_Label ) then

                !----------------------------------------------------------------------------------------------
                ! Get the mirror particle of the inlet/outlet boundary particle
                ! Calculate the coordinates of the mirror particle 
                if ( subdomain_particle_type(i)==Inlet_Particle_Label ) then        ! (For inlet particle)

                    Mirror_Particle_Position(1) = 2*Inlet_Boundary_X - subdomain_particle_position(i,1)
                    Mirror_Particle_Position(2) = subdomain_particle_position(i,2)

                elseif ( subdomain_particle_type(i)==Outlet_Particle_Label ) then   ! (For outlet particle)

                    Mirror_Particle_Position(1) = 2*Outlet_Boundary_X - subdomain_particle_position(i,1)
                    Mirror_Particle_Position(2) = subdomain_particle_position(i,2)

                    ! If we use Calos method, it is only applied for the outlet particle
                    if( Calos_Non_Reflection_Boundary==1 ) then
                        subdomain_in_non_reflection_zone_or_not(i)=1
                    endif

                endif

                ! Calculate the grid name of the particle in the subdomain background grid
                call Point_Project_Into_Subdomain_Background_Grid( Mirror_Particle_Position,Grid_Position_Index,mesh_name )

                near_mesh_number=subdomain_chain_near_mesh_number(mesh_name)      ! Vicinity grid index of background grid index of current particle

                !----------------------------------------------------------------------------------------------


                !----------------------------------------------------------------------------------------------
                ! Initialized the Variables
                W_kernel    = 0.0d0
                DwDx_kernel = 0.0d0

                Ture_effect_number = 0                                            ! Ture effect number
                Ture_effect_name   = 0                                            ! Effect particle index
                !----------------------------------------------------------------------------------------------

                ! Search all the vicinity particles of the mirror particles
                do k=1,near_mesh_number                         
                   
                    near_mesh_name=subdomain_chain_near_mesh_name(mesh_name,k)    ! Current background grid index
                    particle_in_mesh=subdomain_number(near_mesh_name)             ! Particle number in current background grid index
                   
                    do v=1,particle_in_mesh                                       ! Search all the particles in current background grid
            
                        j=subdomain_particle_in_chain(near_mesh_name,v)

                        if( subdomain_particle_type(j)==-3 ) cycle                ! No other boundary particles

                        ! Position difference and distance mirror particle and j particle
                        position_difference(:)=Mirror_Particle_Position(:)-subdomain_particle_position(j,:)
                        distance=sqrt(DOT_PRODUCT(position_difference,position_difference))
                        
                        ! Get the average smooth length (For boundary interpolation, it should be a little large)
                        average_smooth_length = 1.2*smooth_length
                        
                        ! In the support domain or not
                        if( distance <= kernel_scale*average_smooth_length ) then
                            
                            Ture_effect_number=Ture_effect_number+1               ! Ture effect particle number
                            Ture_effect_name(Ture_effect_number)=j                ! Ture effect particle index
                            
                            ! Get the kernel function
                            call compute_kernel(Kernel_function_model,&           ! kernel function model(in)  
                                                dim,&                             ! dimension(in) 
                                                distance,&                        ! distance(in)                
                                                position_difference,&             ! position difference(in)
                                                average_smooth_length,&           ! average smooth length(in)
                                                temp_w,&                          ! kernel function value(out)
                                                temp_dwdx&                        ! derivates of kernel function(out)
                                                )
                        
                            ! Kernel function value
                            W_kernel(Ture_effect_number)      = temp_w   *subdomain_particle_mass(j)/subdomain_particle_rho(j)
                            DwDx_kernel(Ture_effect_number,:) = temp_dwdx*subdomain_particle_mass(j)/subdomain_particle_rho(j)

                        endif   

                    enddo   
                    
                enddo
                !----------------------------------------------------------------------------------------------


                !----------------------------------------------------------------------------------------------
                ! Interpolate in/outlet particle quality
                
                subdomain_Particle_Interpolation_Velocity(i,:)=0.0d0
                subdomain_Particle_Interpolation_Pressure(i)=0.0d0

                if( Ture_effect_number>=3 ) then
                

                    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                    ! Attention: Two method result are nearly same, as compatibility interpolation is faster than MLS so we choose Method 1
                    ! Method 1: Get physical variable value and their derivates by compatibility interpolation
                    call Compatibility_Interpolation_In_Subdomain_Domain(Mirror_Particle_Position,Ture_effect_number,Ture_effect_name,W_kernel,DwDx_kernel,Compatibility_Result_Matrix,check_error)
               
                    ! Compatibility interpolation successful
                    if ( check_error/=0 ) then

                        !--------------------------------------------------------------------------------------
                        ! (Interpolation point - Mirror point )*Gradience
                        position_difference(:)=subdomain_particle_position(i,:)-Mirror_Particle_Position(:)
                        
                        ! Pressure gradience
                        do L=1,dim
                            Gradience(L)=Compatibility_Result_Matrix(L+1,1)
                        enddo

                        ! Pressure interpolation result from mirror particle
                        subdomain_Particle_Interpolation_Pressure(i)=Compatibility_Result_Matrix(1,1)+DOT_PRODUCT(position_difference,Gradience)


                        ! Velocity interpolation result from mirror particle

                        ! For inlet, the velocity are defined by the inlet definition; 
                        ! And the pressure and velocity are exterpolation from the fluid domain;
                        ! If you need to interpolate the inlet particle velocity, such as in SPHERIC Benchmark 6.
                        ! The velocity gradience at inlet is from upstream to downstream, then the inlet velocity will be small, so a shock wave will generated at inlet
                        ! For the outlet has the same problem, if the gradience at outlet from upstream to down stream, it will be larger at outlet
                        ! So you'd better assign the inlet particle velocity, or the velocity should not inclued in the interpolation
                        ! When the inlet or outlet is far from the interior body, which means the flow is stable,
                        ! Then the velocity gradience can be inclued
                        if ( subdomain_particle_type(i)==Inlet_Particle_Label ) then
                        

                            !----------------------------------------------------------------------------------
                            ! Assign the inlet velocity directly
                            ! subdomain_Particle_Interpolation_Velocity(i,:) = Ramping_Inlet_Velocity
                            !----------------------------------------------------------------------------------


                            !----------------------------------------------------------------------------------
                            ! Interpolation without velocity gradience
                            do k=1,dim

                                subdomain_Particle_Interpolation_Velocity(i,k)=Compatibility_Result_Matrix(1,k+1)

                            enddo
                            !----------------------------------------------------------------------------------


                            !----------------------------------------------------------------------------------
                            ! ! Interpolation with velocity gradience
                            ! do k=1,dim

                            !     ! Velocity gradience
                            !     do L=1,dim
                            !         Gradience(L)=Compatibility_Result_Matrix(L+1,k+1)
                            !     enddo

                            !     subdomain_Particle_Interpolation_Velocity(i,k)=Compatibility_Result_Matrix(1,k+1)+DOT_PRODUCT(position_difference,Gradience)

                            ! enddo
                            !----------------------------------------------------------------------------------


                        elseif( subdomain_particle_type(i)==Outlet_Particle_Label ) then
                            

                            !----------------------------------------------------------------------------------
                            ! ! Interpolation without velocity gradience
                            ! do k=1,dim

                            !     subdomain_Particle_Interpolation_Velocity(i,k)=Compatibility_Result_Matrix(1,k+1)

                            ! enddo
                            !----------------------------------------------------------------------------------


                            !----------------------------------------------------------------------------------
                            ! Interpolation with velocity gradience
                            do k=1,dim

                                ! Velocity gradience
                                do L=1,dim
                                    Gradience(L)=Compatibility_Result_Matrix(L+1,k+1)
                                enddo

                                subdomain_Particle_Interpolation_Velocity(i,k)=Compatibility_Result_Matrix(1,k+1)+DOT_PRODUCT(position_difference,Gradience)

                            enddo
                            !----------------------------------------------------------------------------------


                        endif
                        !--------------------------------------------------------------------------------------

                    else

                        !--------------------------------------------------------------------------------------
                        ! If the compatibility interpolation is not successful, return back to MLS interpolation
                        subdomain_Particle_Interpolation_Velocity(i,:) = subdomain_particle_velocity(i,:)
                        subdomain_Particle_Interpolation_Pressure(i)   = subdomain_particle_press(i)

                    endif
                    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


                else

                    subdomain_Particle_Interpolation_Velocity(i,:) = subdomain_particle_velocity(i,:)
                    subdomain_Particle_Interpolation_Pressure(i)   = subdomain_particle_press(i)
                    
                endif
                !----------------------------------------------------------------------------------------------
                ! Calculate the density interpolation results
                subdomain_Particle_Interpolation_Rho(i)=(subdomain_Particle_Interpolation_Pressure(i)-Background_Pressure)/square_c_0+water_rho_0

                ! if ( Current_Processor_ID==0 ) then

                !     write(*,*) subdomain_Particle_Interpolation_Pressure(i),subdomain_Particle_Interpolation_Rho(i),subdomain_Particle_Interpolation_Velocity(i,1)
                    
                ! endif

            endif

        enddo
        !======================================================================================================




        ! !======================================================================================================
        ! ! Only assign inlet boundary particles and outlet keep same updaing by Calos method
        ! ! Reference: Nonrefecting outlet boundary conditions for incompressible fows using SPH
        ! ! (Carlos E. Alvarado-Rodrgueza, Jaime Klapp)
        ! do i=1,Actual_particle_number_in_subdomain

        !     !--------------------------------------------------------------------------------------------------
        !     ! Interpolate the density, pressure and velocity at the inlet/outlet boundary particles
        !     if( subdomain_particle_type(i)==Inlet_Particle_Label ) then

        !         subdomain_particle_velocity(i,:) = subdomain_Particle_Interpolation_Velocity(i,:)
        !         subdomain_particle_press(i)      = subdomain_Particle_Interpolation_Pressure(i)
        !         subdomain_particle_rho(i)        = subdomain_Particle_Interpolation_Rho(i)

        !     endif
        !     !--------------------------------------------------------------------------------------------------

        ! enddo
        ! !====================================================================================================== 





        ! !======================================================================================================
        ! ! Only assign inlet boundary particles and outlet keep same
        ! ! Reference: Simulating 2D open-channel flows through an SPH model
        ! ! (I. Federico , S. Marroneb, A. Colagrossi, F. Aristodemo, M. Antuonoc)
        ! do i=1,Actual_particle_number_in_subdomain

        !     !--------------------------------------------------------------------------------------------------
        !     ! Interpolate the density, pressure and velocity at the inlet/outlet boundary particles
        !     if( subdomain_particle_type(i)==Inlet_Particle_Label ) then

        !         subdomain_particle_velocity(i,:) = subdomain_Particle_Interpolation_Velocity(i,:)
        !         subdomain_particle_press(i)      = subdomain_Particle_Interpolation_Pressure(i)
        !         subdomain_particle_rho(i)        = subdomain_Particle_Interpolation_Rho(i)

        !     endif
        !     !--------------------------------------------------------------------------------------------------

        ! enddo
        ! !====================================================================================================== 





        ! !======================================================================================================
        ! ! Assign inlet/outlet boundary particles
        ! ! Reference: A versatile alogrithm for the treatment of open boundary conditions in Sommthed particle hydrodynamics GPU models
        ! ! (A. Tafuni, R. Vacondio)
        ! do i=1,Actual_particle_number_in_subdomain

        !     !--------------------------------------------------------------------------------------------------
        !     ! Interpolate the density, pressure and velocity at the inlet/outlet boundary particles
        !     if( subdomain_particle_type(i)<=Inlet_Particle_Label ) then

        !         subdomain_particle_velocity(i,:) = subdomain_Particle_Interpolation_Velocity(i,:)
        !         subdomain_particle_press(i)      = subdomain_Particle_Interpolation_Pressure(i)
        !         subdomain_particle_rho(i)        = subdomain_Particle_Interpolation_Rho(i)

        !     endif
        !     !--------------------------------------------------------------------------------------------------

        ! enddo
        ! !====================================================================================================== 




        !======================================================================================================
        ! Assign Non-reflecting inlet boundary particles
        ! Reference: Permeable and non-reflecting boundary conditions in SPH (Martin Lastiwka, Mihai Basa)
        do i=1,Actual_particle_number_in_subdomain


           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           ! Base on the SPHERIC Benchmark Test 6---Moving square in a box.
           ! The in/outlet boundary should far away enough from the interaction body
           ! And the reference is not the interpolation, it should be the far-filed fluid physical attributes
           ! Which means: (a) u_reference    = Ramping_Inlet_Velocity(1)
           !              (b) rho_reference  = water_rho_0
           !              (c) u_reference    = 0.0d0+Background_Pressure

           ! For the inlet boundary particles
           if( subdomain_particle_type(i)==Inlet_Particle_Label ) then

               !-----------------------------------------------------------------------------------------------
               ! Define the reference value
               ! For the inlet particles, the velocity are defined by the boundary conditions
               ! For the Pressure and density, they were exterpolate from the interior domain
               u_reference        = Ramping_Inlet_Velocity(1)                     ! Horzontial velocity
               rho_reference      = water_rho_0
               pressure_reference = 0.0d0

               subdomain_particle_velocity(i,:) = subdomain_Particle_Interpolation_Velocity(i,:)

               ! When the wave have moving out the domain, we should stop the Non-reflecting boundary
               if ( 0.0d0 < Real_physical_time .and. Real_physical_time <= Shock_wave_release_time ) then

                    Temp_c=sqrt(square_c_0)*(subdomain_Particle_Interpolation_Rho(i)/water_rho_0)**3
                   
                    J_1= 0.0d0;
                    J_2= 0.0d0;
                    J_3= -subdomain_Particle_Interpolation_Rho(i)*Temp_c*(subdomain_Particle_Interpolation_Velocity(i,1)-u_reference)+(subdomain_Particle_Interpolation_Pressure(i)-pressure_reference) 

                    subdomain_particle_velocity(i,1)=u_reference+(J_2-J_3)/(2*subdomain_Particle_Interpolation_Rho(i)*Temp_c)
                    subdomain_particle_rho(i)=rho_reference+(-J_1+0.5d0*J_2+0.5d0*J_3)/(Temp_c**2)
                    subdomain_particle_press(i)=pressure_reference+0.5d0*(J_2+J_3)

               else

                    subdomain_particle_press(i)      = subdomain_Particle_Interpolation_Pressure(i)
                    subdomain_particle_rho(i)        = subdomain_Particle_Interpolation_Rho(i)

               endif
               !-----------------------------------------------------------------------------------------------

            ! For the outlet boundary particles
            elseif( subdomain_particle_type(i)==Outlet_Particle_Label ) then

               !-----------------------------------------------------------------------------------------------
               ! Define the reference value (As we dont know the exact physical value of the outlet, so we use the current value as the reference)
            
               ! Reference value is the present outlet particle value
               u_reference       = Ramping_Inlet_Velocity(1)
               rho_reference     = water_rho_0
               pressure_reference= 0.0d0

               subdomain_particle_velocity(i,:) = subdomain_Particle_Interpolation_Velocity(i,:)

               ! When the wave have moving out the domain, we should stop the Non-reflecting boundary
               if ( 0.0d0 < Real_physical_time .and. Real_physical_time <= Shock_wave_release_time ) then

                    Temp_c=sqrt(square_c_0)*(subdomain_Particle_Interpolation_Rho(i)/water_rho_0)**3
                   
                    J_1=-Temp_c**2*(subdomain_Particle_Interpolation_Rho(i)-rho_reference)+(subdomain_Particle_Interpolation_Pressure(i)-pressure_reference)
                    J_2= Temp_c*subdomain_Particle_Interpolation_Rho(i)*(subdomain_Particle_Interpolation_Velocity(i,1)-u_reference)+(subdomain_Particle_Interpolation_Pressure(i)-pressure_reference)
                    J_3= 0.0d0;


                    subdomain_particle_velocity(i,1)=u_reference+(J_2-J_3)/(2*subdomain_Particle_Interpolation_Rho(i)*Temp_c)
                    subdomain_particle_rho(i)=rho_reference+(-J_1+0.5d0*J_2+0.5d0*J_3)/(Temp_c**2)
                    subdomain_particle_press(i)=pressure_reference+0.5d0*(J_2+J_3)

               else

                    subdomain_particle_press(i)      = subdomain_Particle_Interpolation_Pressure(i)
                    subdomain_particle_rho(i)        = subdomain_Particle_Interpolation_Rho(i)

                endif
               !-----------------------------------------------------------------------------------------------


            endif
           !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


        enddo   
        !======================================================================================================
        

    endif
    !**********************************************************************************************************




    !*******************************Assign the Fixed boundary particles information****************************

    ! Assign the boundary particles directly
    FixedBoundary:if (DistributeFixedBoundaryParticleOrNot==1) then
        
        do i=1,Actual_particle_number_in_subdomain

            if( subdomain_particle_type(i)/=Fix_Ghost_Particle_Label ) cycle

            !Initialized the Variables
            subdomain_particle_velocity(i,:)=0.0d0
            subdomain_particle_press(i)=0.0d0
            subdomain_particle_volume(i)=0.0d0
            subdomain_particle_mass(i)=0.0d0
            
            !-------------------------------------------------------------------------------------------------
            W_kernel=0.0d0
            Ture_effect_number=0                                      
            Ture_effect_name=0                                        
            
            ! Search all the vicinity particles
            do L=1,subdomain_effect_particle_number(i)

                ! Effect particle index
                j=subdomain_effect_particle(i,L)

                if( subdomain_particle_type(j)==Fix_Ghost_Particle_Label ) cycle      ! No other boundary particles
                
                ! Position difference and distance between i and j particle
                position_difference(:)=subdomain_particle_position(i,:)-subdomain_particle_position(j,:)
                distance=sqrt(DOT_PRODUCT(position_difference,position_difference))
                
                ! Average smooth length
                average_smooth_length=1.2*smooth_length
                
                ! In the support domain or not
                if(distance<= kernel_scale*average_smooth_length ) then
                    
                    Ture_effect_number=Ture_effect_number+1              ! Ture effect particle number
                    Ture_effect_name(Ture_effect_number)=j               ! Ture effect particle index
                    
                    ! Get the kernel function
                    call compute_kernel(Kernel_function_model,&          ! Kernel function model(in)  
                                        dim,&                            ! dimension
                                        distance,&                       ! distance(in)                
                                        position_difference,&            ! position difference(in)
                                        average_smooth_length,&          ! average smooth length(in)
                                        temp_w,&                         ! kernel function value(out)
                                        temp_dwdx&                       ! derivates of kernel function(out)
                                        )
                
                    ! Kernel function value
                    W_kernel(Ture_effect_number)=temp_w*subdomain_particle_mass(j)/subdomain_particle_rho(j)
      
                endif      
                
            enddo
            !--------------------------------------------------------------------------------------------------



            !--------------------------------------------------------------------------------------------------
            ! If body particle
            if ( subdomain_particle_initial_type(i)==Body_Particle_Label ) then

                External_force_grad=Gravity_grad+Custom_accleration

            else
                
                External_force_grad=Gravity_grad
                
            endif
            !--------------------------------------------------------------------------------------------------


            !--------------------------------------------------------------------------------------------------
            ! Interpolate fixed particle quality
            if(Ture_effect_number>=3) then
                     
                ! Get the weight value of the effect particle
                call MLS_Interpolation_In_Subdomain_Domain(MLS_order_Boundary,subdomain_particle_position(i,:),Ture_effect_number,Ture_effect_name,W_kernel,Fai)
           
                ! Interpolation
                do L=1,Ture_effect_number
            
                    ! Effect particle index
                    j=Ture_effect_name(L)

                    ! Wall particle - fluid particle
                    position_difference(:)=subdomain_particle_position(i,:)-subdomain_particle_position(j,:)

                    subdomain_Particle_Interpolation_Pressure(i)   = subdomain_Particle_Interpolation_Pressure(i)+subdomain_particle_press(j)*Fai(L)+subdomain_particle_rho(j)*Fai(L)*DOT_PRODUCT(External_force_grad,position_difference)

                    subdomain_Particle_Interpolation_Velocity(i,:) = subdomain_Particle_Interpolation_Velocity(i,:)+subdomain_particle_velocity(j,:)*Fai(L)

                    subdomain_particle_mass(i)=subdomain_particle_mass(i)+subdomain_particle_mass(j)*Fai(L)

                enddo
                
            endif
            !--------------------------------------------------------------------------------------------------








            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Method 1: The reference value is the interpolation value (This method is better)
            
            !--------------------------------------------------------------------------------------------------
            ! Permeable and non-reflecting boundary
            ! Reference: Permeable and non-reflecting boundary conditions in SPH (Martin Lastiwka, Mihai Basa)

            ! For the external boundary particles
            if( subdomain_particle_initial_type(i)/=Body_Particle_Label ) then

                !----------------------------------------------------------------------------------------------
                ! Define the reference value
                ! For the inlet particles, the velocity are defined by the boundary conditions
                ! For the Pressure and density, they were exterpolate from the interior domain
                u_reference        = 0.0d0                                  ! Horzontial velocity
                rho_reference      = water_rho_0
                pressure_reference = 0.0d0

                ! When the wave have moving out the domain, we should stop the Non-reflecting boundary
                if ( 0.0d0 < Real_physical_time .and. Real_physical_time <= Shock_wave_release_time ) then

                    subdomain_Particle_Interpolation_Rho(i)=(subdomain_Particle_Interpolation_Pressure(i)-Background_Pressure)/square_c_0+water_rho_0
                    Temp_c=sqrt(square_c_0)*(subdomain_Particle_Interpolation_Rho(i)/water_rho_0)**3

                    subdomain_particle_velocity(i,:)=subdomain_Particle_Interpolation_Velocity(i,:) 

                    ! Horzontial wall
                    if ( subdomain_Boundary_particle_type(i)==-6 ) then
                        
                        J_1= -Temp_c**2*(subdomain_Particle_Interpolation_Rho(i)-rho_reference)+(subdomain_Particle_Interpolation_Pressure(i)-pressure_reference);
                        J_2=  subdomain_Particle_Interpolation_Rho(i)*Temp_c*(subdomain_Particle_Interpolation_Velocity(i,2)-u_reference)+(subdomain_Particle_Interpolation_Pressure(i)-pressure_reference) 
                        J_3= -subdomain_Particle_Interpolation_Rho(i)*Temp_c*(subdomain_Particle_Interpolation_Velocity(i,2)-u_reference)+(subdomain_Particle_Interpolation_Pressure(i)-pressure_reference) 

                        subdomain_particle_velocity(i,2)=u_reference+(J_2-J_3)/(2*subdomain_Particle_Interpolation_Rho(i)*Temp_c)
                        subdomain_particle_press(i)=pressure_reference+0.5d0*(J_2+J_3)

                    ! Vertical wall
                    elseif ( subdomain_Boundary_particle_type(i)==-7 ) then

                        J_1= -Temp_c**2*(subdomain_Particle_Interpolation_Rho(i)-rho_reference)+(subdomain_Particle_Interpolation_Pressure(i)-pressure_reference);
                        J_2=  subdomain_Particle_Interpolation_Rho(i)*Temp_c*(subdomain_Particle_Interpolation_Velocity(i,1)-u_reference)+(subdomain_Particle_Interpolation_Pressure(i)-pressure_reference) 
                        J_3= -subdomain_Particle_Interpolation_Rho(i)*Temp_c*(subdomain_Particle_Interpolation_Velocity(i,1)-u_reference)+(subdomain_Particle_Interpolation_Pressure(i)-pressure_reference) 

                        subdomain_particle_velocity(i,1)=u_reference+(J_2-J_3)/(2*subdomain_Particle_Interpolation_Rho(i)*Temp_c)
                        subdomain_particle_press(i)=pressure_reference+0.5d0*(J_2+J_3)

                    endif
                    
                else

                    subdomain_particle_velocity(i,:) = subdomain_Particle_Interpolation_Velocity(i,:)
                    subdomain_particle_press(i)      = subdomain_Particle_Interpolation_Pressure(i)

                endif
                !----------------------------------------------------------------------------------------------

            else

                subdomain_particle_velocity(i,:) = subdomain_Particle_Interpolation_Velocity(i,:)
                subdomain_particle_press(i)      = subdomain_Particle_Interpolation_Pressure(i)

            endif
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^









            !--------------------------------------------------------------------------------------------------
            ! Assign the boundary velocity
            Body_Particle_if: if ( subdomain_particle_initial_type(i)==Body_Particle_Label ) then

                ! Body particle is assigned as custom velocity
                subdomain_particle_velocity(i,:) = Custom_velocity(:)

                !----------------------------------------------------------------------------------------------
                ! ! For the moving body fixed boundary particle
                ! ! Assgin the velocity 

                ! ! Horzontial mirror rule
                ! if(subdomain_Boundary_particle_type(i)==-6) then

                !     if ( BodyBoundaryCondition==1 ) then

                !         subdomain_particle_velocity(i,1) =  subdomain_particle_velocity(i,1)
                !         subdomain_particle_velocity(i,2) = -subdomain_particle_velocity(i,2)

                !     elseif ( BodyBoundaryCondition==0 ) then

                !         subdomain_particle_velocity(i,1) = 2*Custom_velocity(1)-subdomain_particle_velocity(i,1)
                !         subdomain_particle_velocity(i,2) = 2*Custom_velocity(2)-subdomain_particle_velocity(i,2)

                !     else

                !         write(*,'(A)') " The 'BodyBoundaryCondition' variable is not defined successful. (Assign Boundary)"
                        
                !     endif

                ! ! Center mirror rule
                ! elseif(subdomain_Boundary_particle_type(i)==-7) then
            
                !     subdomain_particle_velocity(i,1) = 2*Custom_velocity(1)-subdomain_particle_velocity(i,1)
                !     subdomain_particle_velocity(i,2) = 2*Custom_velocity(2)-subdomain_particle_velocity(i,2)

                ! ! Vertical mirror rule
                ! elseif(subdomain_Boundary_particle_type(i)==-8) then         
                    
                !     if ( BodyBoundaryCondition==1 ) then

                !         subdomain_particle_velocity(i,1) = -subdomain_particle_velocity(i,1)
                !         subdomain_particle_velocity(i,2) =  subdomain_particle_velocity(i,2)

                !     elseif ( BodyBoundaryCondition==0 ) then

                !         subdomain_particle_velocity(i,1) = 2*Custom_velocity(1)-subdomain_particle_velocity(i,1)
                !         subdomain_particle_velocity(i,2) = 2*Custom_velocity(2)-subdomain_particle_velocity(i,2)

                !     else

                !         write(*,'(A)') " The 'BodyBoundaryCondition' variable is not defined successful. (Assign Boundary)"
                        
                !     endif
                
                ! elseif(subdomain_Boundary_particle_type(i)==-9) then
                    
                !     ! Fixed particle along the circle
                !     call radius_free_slip(particle_position(i,1),particle_position(i,2),particle_velocity(i,1),particle_velocity(i,2),particle_velocity(fix_ghost_name,1),particle_velocity(fix_ghost_name,2))

                ! endif
                !----------------------------------------------------------------------------------------------

            else

                !----------------------------------------------------------------------------------------------
                ! For external boundary particle
                ! Assgin the velocity 

                ! Horzontial mirror rule
                if(subdomain_Boundary_particle_type(i)==-6) then

                    if ( ExternalBoundaryCondition==1 ) then

                        subdomain_particle_velocity(i,1) =  subdomain_particle_velocity(i,1)
                        subdomain_particle_velocity(i,2) = -subdomain_particle_velocity(i,2)

                    elseif ( ExternalBoundaryCondition==0 ) then

                        subdomain_particle_velocity(i,1) = -subdomain_particle_velocity(i,1)
                        subdomain_particle_velocity(i,2) = -subdomain_particle_velocity(i,2)

                    else

                        write(*,'(A)') " The 'ExternalBoundaryCondition' variable is not defined successful. (Assign Boundary)" 
                        
                    endif

                ! Center mirror rule
                elseif(subdomain_Boundary_particle_type(i)==-7) then
            
                    subdomain_particle_velocity(i,1) = -subdomain_particle_velocity(i,1)
                    subdomain_particle_velocity(i,2) = -subdomain_particle_velocity(i,2)

                ! Vertical mirror rule
                elseif(subdomain_Boundary_particle_type(i)==-8) then         
                    
                    if ( ExternalBoundaryCondition==1 ) then

                        subdomain_particle_velocity(i,1) = -subdomain_particle_velocity(i,1)
                        subdomain_particle_velocity(i,2) =  subdomain_particle_velocity(i,2)

                    elseif ( ExternalBoundaryCondition==0 ) then

                        subdomain_particle_velocity(i,1) = -subdomain_particle_velocity(i,1)
                        subdomain_particle_velocity(i,2) = -subdomain_particle_velocity(i,2)

                    else

                        write(*,'(A)') " The 'ExternalBoundaryCondition' variable is not defined successful. (Assign Boundary)" 
                        
                    endif
                
                elseif(subdomain_Boundary_particle_type(i)==-9) then
                    
                    ! Fixed particle along the circle
                    call radius_free_slip(particle_position(i,1),particle_position(i,2),particle_velocity(i,1),particle_velocity(i,2),particle_velocity(fix_ghost_name,1),particle_velocity(fix_ghost_name,2))

                end if
                !----------------------------------------------------------------------------------------------
                
            endif Body_Particle_if
            !--------------------------------------------------------------------------------------------------



            !--------------------------------------------------------------------------------------------------
            ! Assgin the other qualities of fixed particle 

            if (subdomain_particle_mass(i)/=0.0d0) then

                subdomain_particle_energy(i)=e_0
                !subdomain_particle_rho(i)=water_rho_0*(1+7*subdomain_particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                subdomain_particle_rho(i)=(subdomain_particle_press(i)-Background_Pressure)/square_c_0+water_rho_0
                !subdomain_particle_mass(i)=subdomain_particle_rho(i)*subdomain_particle_volume(i)
                !subdomain_particle_mass(i)=subdomain_particle_rho(i)*intial_volume_0
                subdomain_particle_c(i)=sqrt(square_c_0)*(subdomain_particle_rho(i)/water_rho_0)**3
            
            else

                subdomain_particle_velocity(i,:)=0.0d0
                subdomain_particle_press(i)=0.0d0
                subdomain_particle_energy(i)=e_0
                subdomain_particle_rho(i)=water_rho_0
                subdomain_particle_mass(i)=0.0d0
                subdomain_particle_c(i)=sqrt(square_c_0)

            endif
            !--------------------------------------------------------------------------------------------------

        enddo

    endif FixedBoundary






    ! Exchange the fixed boundary particles in buffer domain
    call Exchange_buffer_domain_information(i_time_step)






    !----------------------------------------------------------------------------------------------------------
    ! Don't forget calculate the particle volume

    do k=1,Total_particle_number_in_subdomain  

       !if (subdomain_particle_rho(k)<1.0e-8) cycle
        
       ! Particle volume                          
       subdomain_particle_volume(k)=subdomain_particle_mass(k)/subdomain_particle_rho(k)

       ! if (isnan(subdomain_particle_volume(k))) then
       !      write(*,*) subdomain_particle_mass(k),subdomain_particle_rho(k)
       !     subdomain_particle_volume(k)=intial_volume_0
       ! end if

       ! if (subdomain_particle_volume(k)<1.0e-8) then
       !    subdomain_particle_volume(k)=intial_volume_0
       ! endif

       !write(*,*) subdomain_particle_volume(k)

    end do
    !----------------------------------------------------------------------------------------------------------


    !**********************************************************************************************************











    ! !**********************************************************************************************************
    ! ! Check the boundary interpolation results for each processors
    ! if( i_time_step==Actual_start_time_step ) then

    !     !------------------------------------------------------------------------------------------------------
    !     ! Open the output tecplot data files
    !     Standby_Current_Processor_File_Index = Standby_Initial_File_Port_For_Each_Processor+Current_Processor_ID

    !     File_index = Standby_Current_Processor_File_Index

    !     ! Transfer the Current_Processor_ID from integer to character
    !     write(Char_Current_Processor_ID,'(I4)') Current_Processor_ID

    !     File_name="./Subdomain/Boundary_in_subdomain_"//trim(adjustl(Char_Current_Processor_ID))//".dat"

    !     call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
   
    !     if ( IOERROR==0 ) then

    !         !--------------------------------------------------------------------------------------------------
    !         ! Tecplot header
    !         write(File_index,*) "TITLE='DISTRIBUTION'"
    !         write(File_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
    !         !--------------------------------------------------------------------------------------------------

    !         !--------------------------------------------------------------------------------------------------
    !         ! Output the result
    !         write(File_index,*) "ZONE I=",Total_particle_number_in_subdomain," F=POINT"
    !         do i=1,Total_particle_number_in_subdomain
    !             write( File_index,'(8F16.8)' ) (subdomain_particle_position(i,j),j=1,dim),(subdomain_particle_velocity(i,j),j=1,dim),subdomain_particle_rho(i),(subdomain_particle_press(i)-Background_Pressure)/pressure_0
    !         enddo
    !         !--------------------------------------------------------------------------------------------------
   
    !     endif

    !     close(File_index)
    
    ! endif
    ! !**********************************************************************************************************



        
    
end subroutine Assign_boundary_particles