!**************************************************************************************************************
!  SUBROUTINE: Carlos_Non_Reflection_Medthod
!
!  PURPOSE   : Calculate the iteration Information based on Carlos Non-Reflection Medthod
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

 subroutine Carlos_Non_Reflection_Medthod()
 
    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI

    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,L,v,o,m                                          ! Variables for loop

    !artifical viscosity term
    real(kind=8),dimension(dim)::x_ij,v_ij                          ! Position and velocity difference
    real(kind=8)::v_ij_multiply_x_ij                                ! x_ij times v_ij
    real(kind=8)::temp,temp_i,temp_j                                ! Temp value
    real(kind=8)::square_distance                                   ! Square distance
    real(kind=8)::average_rho                                       ! Average density

    !==========================================================================================================

    ! Body of subroutine Carlos_Non_Reflection_Medthod

    !**********************************************************************************************************
    ! Reference: Nonreflecting outlet boundary conditions for incompressible flow using SPH
    ! Calculate the outlet boundary particles

    do k=1,pair_number
        
        i=pair_i(k)                                                ! particle i index in kth pair
        j=pair_j(k)                                                ! particle j index in kth pair

        if ( subdomain_in_non_reflection_zone_or_not(i) /= 1 .and. subdomain_in_non_reflection_zone_or_not(j) /= 1 ) cycle

        !------------------------------------------------------------------------------------------------------
        ! Velocity_difference * Position_difference_ij
        x_ij(:)=position_difference_ij(k,:)                        ! i-j
        v_ij(:)=velocity_difference_ij(k,:)                        ! i-j
        average_smooth_length=average_smooth_length_ij(k)
        average_rho=average_rho_ij(k)

        square_distance=DOT_PRODUCT(x_ij,x_ij)
        !------------------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------------------
        ! One dimension wave progration equation in SPH format
        if ( subdomain_in_non_reflection_zone_or_not(i)==1 ) then
            internal_acceleration(i,:)=internal_acceleration(i,:)+subdomain_Particle_Interpolation_Velocity(i,1)*subdomain_particle_mass(j)/average_rho*Normalized_dwdx_i(k,1)*v_ij(:)
        endif

        if ( subdomain_in_non_reflection_zone_or_not(j)==1 ) then
            internal_acceleration(j,:)=internal_acceleration(j,:)+subdomain_Particle_Interpolation_Velocity(j,1)*subdomain_particle_mass(i)/average_rho*Normalized_dwdx_j(k,1)*v_ij(:)
        endif
        !------------------------------------------------------------------------------------------------------


        !------------------------------------------------------------------------------------------------------
        ! Artifical viscousity acceleration

        if ( laplacianTerm==1 ) then

            !--------------------------------------------------------------------------------------------------
            ! It is useful for non-free surface flow
            ! Reference:  [Morris] 
            ! v_ij_multiply_x_ij in Equ.(29) in Liu
            temp_i=0.0d0;  temp_j=0.0d0;    
            do L=2,dim
                temp_i=temp_i+x_ij(L)*Normalized_dwdx_i(k,L)/(square_distance+(fai_coefficient*average_smooth_length)**2)  
                temp_j=temp_j+x_ij(L)*Normalized_dwdx_j(k,L)/(square_distance+(fai_coefficient*average_smooth_length)**2)           
            end do    

            ! The laplacian term is Morris type
            if ( subdomain_in_non_reflection_zone_or_not(i)==1 ) then
                artifical_viscosity_acceleration(i,:)=artifical_viscosity_acceleration(i,:)-2*water_kinematic_viscosity*subdomain_particle_volume(j)*temp_i*v_ij(:)
            endif

            if ( subdomain_in_non_reflection_zone_or_not(j)==1 ) then
                artifical_viscosity_acceleration(j,:)=artifical_viscosity_acceleration(j,:)+2*water_kinematic_viscosity*subdomain_particle_volume(i)*temp_j*v_ij(:)
            endif
            !--------------------------------------------------------------------------------------------------

        elseif ( laplacianTerm==0 ) then

            !--------------------------------------------------------------------------------------------------
            ! It is useful for free surface flow 
            ! Reference: [Theoretical analysis and numerical verification of the consistency of viscous smoothed-particle-hydrodynamics formulations in simulating free-surface flows]
            !            [Adaptive particle refinement and derefinement applied to the smoothed particle hydrodynamics methods]

            v_ij_multiply_x_ij=0.0d0
            do L=2,dim
                v_ij_multiply_x_ij=v_ij_multiply_x_ij+x_ij(L)*v_ij(L)
            end do 

            ! The laplacian term is Morris type
            if ( subdomain_in_non_reflection_zone_or_not(i)==1 ) then
                artifical_viscosity_acceleration(i,:)=artifical_viscosity_acceleration(i,:)-subdomain_particle_volume(j)*2*(dim+2)*water_kinematic_viscosity*v_ij_multiply_x_ij/square_distance*Normalized_dwdx_i(k,:)
            endif

            if ( subdomain_in_non_reflection_zone_or_not(j)==1 ) then
                artifical_viscosity_acceleration(j,:)=artifical_viscosity_acceleration(j,:)+subdomain_particle_volume(i)*2*(dim+2)*water_kinematic_viscosity*v_ij_multiply_x_ij/square_distance*Normalized_dwdx_j(k,:)
            endif
            !--------------------------------------------------------------------------------------------------

        else

            write(*,*) 'The defination of laplacianTerm is not right! (Solver)'
            
        endif
        !------------------------------------------------------------------------------------------------------


    enddo
    !**********************************************************************************************************

 
 end subroutine Carlos_Non_Reflection_Medthod