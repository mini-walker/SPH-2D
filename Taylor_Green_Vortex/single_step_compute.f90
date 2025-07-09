!****************************************************************************
!
!  PROGRAM: single_step_compute
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

subroutine single_step_compute(i_time_step)              

    use information_module
    use Public_variable_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables form the mother subroutine
    integer,intent(in)::i_time_step                               !current time step

    !Variables in subroutine
    integer::i,j,k,L                                              !loop Variables
    real(kind=8)::damping_coefficient                             !damping coefficient
    real(kind=4)::current_time                                    !current time

    real(kind=8)::current_wave_angle                              !current wave angle
    real(kind=8)::time_soft_coefficient                           !soft coefficient for time
    real(kind=8)::current_wave_elevation                          !wave elevation at current particle z direction

    !Variables for file opeartion
    integer::ioerror=0                                            !open return value
    integer::stat                                                 !read return value
    integer::status                                               !allocation return value

    integer::file_index
    character(len=40)::file_name
    character(len=4)::char_Current_Processor_ID

    !==========================================================================================================

    ! Body of subroutine single_step_compute

    ! !*********************************************************************************************************
    ! !Get the MPI runing information for dugging, you can annotate this part when runnging
    ! call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    ! call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    ! !*********************************************************************************************************

    !***********************************************************************************************************
    !Call the EOS(Equation of Status)
    !Get the particle pressure
    do i=1,total_particle_number_in_subdomain
       
        if(subdomain_particle_type(i)==120) cycle                        !No error particles

        call status_equation(subdomain_particle_rho(i),&                 !Subdomain Particle Density(in)
                             subdomain_particle_energy(i),&              !Subdomain Particle Energy(in)
                             subdomain_particle_type(i),&                !Subdomain Particle Type(in)
                             subdomain_particle_press(i),&               !Subdomain Particle Pressure(out)
                             subdomain_particle_c(i)&                    !Subdomain Particle Sound Velocity(out)
                             )
    
    end do
    !***********************************************************************************************************

    !***********************************************************************************************************
    !Assign boundary particles
    subdomain_particle_volume=0.0d0             !particle volume(out)
    
    call Assign_boundary_particles(i_time_step)

    !***********************************************************************************************************

    !***********************************************************************************************************
    !Initailized the Variables
    pair_i=0                                    !particle i index in kth pair (out) 
    pair_j=0                                    !particle j index in kth pair (out)
    w=0.0d0                                     !weight function value in kth pair(out)
    dwdx=0.0d0                                  !derivats of weight function in kth pair(out)
    pair_number=0                               !pair number(out)
    
    !Get the particle pairs information
    call find_particle_pair()                   
    
    !!Check w,dwdx
    !write(*,*) Current_Processor_ID,all_particle
    !do i=1,pair_number
    !    write(3,*) i,pair_i(i),pair_j(i),w(i),dwdx(i,1),dwdx(i,2)
    !end do
    !***********************************************************************************************************

    !***********************************************************************************************************
    !Define the free surface and normalized the derivates of the kernel function
    subdomain_fra_v=0.0d0
    Tensor_dyadic_matrix=0.0d0
    Tensor_dyadic_inverse_matrix=0.0d0
    
    call Define_freesurface_and_Normailized_derivates()

    !***********************************************************************************************************
    
    !***********************************************************************************************************
    !Initailized the acceleration and energy change rates
    drhodt=0.0d0
    dedt=0.0d0
    internal_acceleration=0.0d0  
    artifical_internal_acceleration=0.0d0     
    artifical_viscosity_acceleration=0.0d0 
    artifical_viscosity_dedt=0.0d0
    artifical_heat_dedt=0.0d0
    external_acceleration=0.0d0      
    average_velocity=0.0d0
    artifical_rho_correct=0.0d0 
    laplacian_viscous_acceleration=0.0d0
    SPS_acceleration=0.0d0

    !call Solver_Standard_SPH()                                         !Standard SPH Solver
    !call Solver_Delta_SPH()                                            !Delta SPH Solver
    !call Solver_Modified_Delta_SPH()                                   !Modified Delta SPH Solver
    !call Solver_Standard_SPH_SPS()                                     !Standard SPH+SPS Solver
    !call Solver_Normalized_Standard_SPH()                              !Normalized Standard SPH Solver
    !call Solver_Normalized_Standard_SPH_SPS()                          !Normalized Standard SPH+SPS Solver
    !call Solver_Normalized_Delta_SPH_SPS()                             !Normalized Delta SPH+SPS Solver
    !call Solver_Normalized_Modified_Delta_SPH_SPS()                    !Normalized Delta SPH+SPS Solver

    !Modified Riemann is better for pressure
    call Solver_Normalized_Modified_Linear_Riemann_Delta_SPS()         !Modified Linear Riemann (Low dissipation) Normalized Delta SPH SPS Solver
    !call Solver_Normalized_Modified_Linear_Riemann_Delta()             !Modified Linear Riemann (Low dissipation) Normalized Delta SPH Solver
    
    !call Solver_Normalized_Linear_Riemann_Delta_SPS()                  !Linear Riemann Normalized Delta SPH SPS Solver
    
    !HLLC is not stable
    !call Solver_Normalized_HLLC_Linear_Riemann_Delta_SPS()             !HLLC Linear Riemann (Low dissipation) Normalized Delta SPH SPS Solver
    
    !Higer order MUSCL approach for Riemann solver
    !call Solver_Normalized_MUSCL_Linear_Riemann_Delta_SPS()             !MUSCL Linear Riemann (Low dissipation) Normalized Delta SPH SPS Solver
    !***********************************************************************************************************

     
    !***********************************************************************************************************
    !Assemble all the rates
    particle_acceleration=0.0d0

    do i=1,actual_particle_number_in_subdomain                
       
        if(subdomain_particle_type(i)/=2) cycle                !Only for the fluid particles

        !renew the particles acceleration
        do j=1,dim
            
            if(j==dim .and. trim(adjustl(Gravity_Switch))=='on') then
                external_acceleration(i,j)=external_acceleration(i,j)-g
            end if
    
            particle_acceleration(i,j)=internal_acceleration(i,j)+external_acceleration(i,j)+artifical_viscosity_acceleration(i,j)+artifical_internal_acceleration(i,j)+SPS_acceleration(i,j)+laplacian_viscous_acceleration(i,j) 
        end do
        
        !renew the particles energy
        dedt(i)=dedt(i)+artifical_viscosity_dedt(i)+artifical_heat_dedt(i)
           
    end do
    !***********************************************************************************************************

    !==========================================================================================
    !For the end damping zone (current velocity minus the damping coefficient times velocity directly)
    !The value of Damping coefficient from 0 to 1 (From the start point to the end of the damping zone)
    !Attention: we just need damping the velocity, the other parameters are not included in the damping zone
    EndDamping_if:if (trim(adjustl(EndDampingZone_Switch))=='on' .and. make_wave_or_not==1) then


        if (trim(adjustl(EndDampingZoneType))=='Cosin_Type') then

            !----------------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------------

        elseif (trim(adjustl(EndDampingZoneType))=='OpenFoam_Type') then

            !----------------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------------

        elseif (trim(adjustl(EndDampingZoneType))=='DualSPH_Type') then

            !----------------------------------------------------------------------------------
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
            !----------------------------------------------------------------------------------
            
        endif
        
    endif EndDamping_if
    !==========================================================================================
    



!     !------------------------------------------------------------------------------------------
!     !open the output tecplot data files
!     file_index=1+Current_Processor_ID

!     !transfer the Current_Processor_ID from integer to character
!     write(char_Current_Processor_ID,'(I4)') Current_Processor_ID
!     !write(*,*) char_Current_Processor_ID

!     file_name="Check_Iteration"//trim(adjustl(char_Current_Processor_ID))//".dat"
!     !write(*,*) file_name

!     open(unit=file_index,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    

!     !output tecplot header
!     write(file_index,*) "TITLE='DISTRIBUTION'"
!     write(file_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "

!     !output the results
!     write(file_index,*) "ZONE I=",total_particle_number_in_subdomain," F=POINT"
!     do i=1,total_particle_number_in_subdomain

!         write(file_index,100) (subdomain_particle_position(i,j),j=1,dim),(particle_acceleration(i,j),j=1,dim),subdomain_particle_rho(i),subdomain_particle_press(i)/pressure_0
! 100      format(8F20.10) 
!     end do

!     close(file_index)
!     !------------------------------------------------------------------------------------------
    
    
!     !Note the user the current active Processor ID have finished the iterator information calculation
!     write(*,10) Current_Processor_ID
! 10   Format(" Processor:",I3," has finished iterator information calculation!")
    
    
    
end subroutine single_step_compute  