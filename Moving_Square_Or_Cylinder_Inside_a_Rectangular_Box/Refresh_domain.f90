!**************************************************************************************************************
!
!  SUBROUTINE : Refresh_domain
!
!  PURPOSE    : Refresh domain
!
!  Programer  : Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location   : MUN
!
!  Time       : 2017.3.18
!
!  Copyright  : Memorial University
!
!  Version    : 1.0
!
!  Note       : MPI version: mpich-3.2
!               Fortran version: Inter Fortran (ifort) 18.0.0
!**************************************************************************************************************

subroutine Refresh_domain(i_time_step)

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the superior subroutine
    integer,intent(in)::i_time_step                                     ! Current time step 

    ! Variables in program
    integer::i,j,k,L,m                                                  ! Variables for loop

    !==========================================================================================================

    ! Body of Refresh_domain

    !**********************************************************************************************************
    ! Some Variables donot need to be Refreshed when they are not updated during the calculation
    !----------------------------------------------------------------------------------------------------------
    ! Refresh domain

    ! As we have distribute all particles into all subdomain, then we can define all domain information to be zero


    call MPI_RealVector_Message_Refresh(subdomain_particle_position,particle_position)
    call MPI_RealVector_Message_Refresh(subdomain_particle_velocity,particle_velocity)

    call MPI_Real_Message_Refresh(subdomain_particle_press,particle_press)
    call MPI_Real_Message_Refresh(subdomain_particle_mass,particle_mass)
    call MPI_Real_Message_Refresh(subdomain_particle_rho,particle_rho)
    call MPI_Real_Message_Refresh(subdomain_particle_volume,particle_volume)
    call MPI_Real_Message_Refresh(subdomain_particle_energy,particle_energy)


    call MPI_Int_Message_Refresh(subdomain_particle_initial_type,particle_initial_type)
    call MPI_Int_Message_Refresh(subdomain_particle_type,particle_type)
    call MPI_Int_Message_Refresh(subdomain_particle_division_degree,particle_division_degree)
    

    ! Only when we have free surface and the step we need for output
    if(mod(i_time_step,save_time_step)==0 .and. Has_FreeSurface_Or_Not==1) then
        call MPI_Int_Message_Refresh(subdomain_free_surface_type,free_surface_type)
        call MPI_RealVector_Message_Refresh(subdomain_Normal_vector,Normal_vector)
    endif


    ! As we donot modify the smooth length, so we donot refresh this values
    ! call MPI_Real_Message_Refresh(subdomain_particle_smooth_lengh,particle_smooth_lengh)

    !----------------------------------------------------------------------------------------------------------
    ! As the vorticity is only used for post-proceeding, so only refresh when doing output
    if(mod(i_time_step,save_time_step)==0) then

        !------------------------------------------------------------------------------------------------------
        ! Calculate the vorticity

        subdomain_particle_vorticity=0.0d0

        do k=1,pair_number
            
           i=pair_i(k)                        !particle i index in kth pair
           j=pair_j(k)                        !particle j index in kth pair

           velocity_difference(:)=velocity_difference_ij(k,:)

           subdomain_particle_vorticity(i,:)=subdomain_particle_vorticity(i,:)+subdomain_particle_volume(j)*Cross_Product(velocity_difference,Normalized_dwdx_i(k,:),dim)
           
           subdomain_particle_vorticity(j,:)=subdomain_particle_vorticity(j,:)+subdomain_particle_volume(i)*Cross_Product(velocity_difference,Normalized_dwdx_j(k,:),dim)

        end do

        
        !particle vorticity
        particle_vorticity=0.0d0
        Temp_Reduce_Real_Array=0.0d0
        
        do j=1,3

          do i=1,actual_particle_number_in_subdomain                      
            
             Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_vorticity(i,j)                            

          end do

          call MPI_REDUCE( Temp_Reduce_Real_Array,particle_vorticity(:,j),n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
          call MPI_BCAST( particle_vorticity(:,j),n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
     
        end do
        !------------------------------------------------------------------------------------------------------

        ! As the vorticity has 3 dimension for 2D case so we cannot use the normal refresh subroutine
        
        ! call MPI_RealVector_Message_Refresh(subdomain_particle_vorticity,particle_vorticity)

    endif
    !----------------------------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------------------------
    ! Reduce the body force
    if ( Calculate_Body_Force_or_Not==1 .and. mod(i_time_step,sampling_time_step)==0 ) then

        ! We should Initialzed the variable before reducing

        Body_Pressure_Force=0.0d0                               

        call MPI_REDUCE( subdomain_Body_Pressure_Force,Body_Pressure_Force,dim,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
        call MPI_BCAST( Body_Pressure_Force,dim,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

        Body_Viscous_Force=0.0d0                               

        call MPI_REDUCE( subdomain_Body_Viscous_Force,Body_Viscous_Force,dim,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
        call MPI_BCAST( Body_Viscous_Force,dim,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

        Body_Total_Force=0.0d0                               

        call MPI_REDUCE( subdomain_Body_Total_Force,Body_Total_Force,dim,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
        call MPI_BCAST( Body_Total_Force,dim,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    endif
    !----------------------------------------------------------------------------------------------------------



  end subroutine Refresh_domain

