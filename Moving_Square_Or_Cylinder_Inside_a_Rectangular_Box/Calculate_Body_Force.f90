!**************************************************************************************************************
!
!  SUBROUTINE : Calculate_Body_Force
!
!  PURPOSE    : Calculate the body force
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


subroutine Calculate_Body_Force()
 
    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI

    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,l,v,o,m                                          ! Variables for looping

    
    ! Artifical viscosity term
    real(kind=8),dimension(dim)::x_ij,v_ij                          ! Position and velocity difference
    real(kind=8)::v_ij_multiply_x_ij                                ! x_ij times v_ij
    real(kind=8)::Temp                                              ! Temp value
    real(kind=8)::square_distance                                   ! Square distance

    ! Laminar laplacian term
    real(kind=8)::RijDotProductDerivates                            ! Rij*D Wij
    !==========================================================================================================

    ! Body of subroutine Calculate_Body_Force

    !*******************************************Calculate Body Force*******************************************
    subdomain_Near_Body_Or_Not=0

    do k=1,pair_number
        
        i=pair_i(k)                                                 ! Particle i index in kth pair
        j=pair_j(k)                                                 ! Particle j index in kth pair

        ! For the pair do not contain the body particles
        ! Particle initial type is used for body force calculation and body motion solver                       
        if ( subdomain_particle_initial_type(i)/=Body_Particle_Label .and. subdomain_particle_initial_type(j)/=Body_Particle_Label ) cycle

        ! Velocity_difference*Position_difference_ij
        v_ij(:)=velocity_difference_ij(k,:)
        x_ij(:)=position_difference_ij(k,:)

        square_distance=DOT_PRODUCT(x_ij,x_ij)
        v_ij_multiply_x_ij=DOT_PRODUCT(v_ij,x_ij)

        ! For i is the fluid particles
        if ( subdomain_particle_initial_type(j)==Body_Particle_Label ) then

            ! I th particle near the body
            subdomain_Near_Body_Or_Not(i)=1

            ! Pressure term
            Body_Force_Pressure_Acceleration(i,:)=Body_Force_Pressure_Acceleration(i,:)-subdomain_particle_volume(j)*(subdomain_particle_press(i)+subdomain_particle_press(j))*Normalized_dwdx_i(k,:)

            ! Viscous term (These two viscous force are nearly same)
            if ( laplacianTerm==1 ) then

                !----------------------------------------------------------------------------------------------
                !It is useful for non-free surface flow
                !Reference:  [Morris]  

                RijDotProductDerivates=DOT_PRODUCT(x_ij,Normalized_dwdx_i(k,:))
                Body_Force_Viscous_Acceleration(i,:)=Body_Force_Viscous_Acceleration(i,:)+subdomain_particle_rho(i)*subdomain_particle_volume(j)*2*water_kinematic_viscosity*RijDotProductDerivates/square_distance*v_ij(:)
              
                !----------------------------------------------------------------------------------------------

            elseif ( laplacianTerm==0 ) then

                !----------------------------------------------------------------------------------------------
                !It is useful for free surface flow 
                !Reference: [Theoretical analysis and numerical verification of the consistency of viscous smoothed-particle-hydrodynamics formulations in simulating free-surface flows]
                !           [Adaptive particle refinement and derefinement applied to the smoothed particle hydrodynamics methods]

                Body_Force_Viscous_Acceleration(i,:)=Body_Force_Viscous_Acceleration(i,:)+subdomain_particle_rho(i)*subdomain_particle_volume(j)*2*(dim+2)*water_kinematic_viscosity*v_ij_multiply_x_ij/square_distance*Normalized_dwdx_i(k,:)
                
                !----------------------------------------------------------------------------------------------

            endif

        endif


        ! For j is the fluid particles
        if ( subdomain_particle_initial_type(i)==Body_Particle_Label ) then

            ! I th particle near the body
            subdomain_Near_Body_Or_Not(j)=1

            ! Pressure term
            Body_Force_Pressure_Acceleration(j,:)=Body_Force_Pressure_Acceleration(j,:)+subdomain_particle_volume(i)*(subdomain_particle_press(i)+subdomain_particle_press(j))*Normalized_dwdx_j(k,:)

            ! Viscous term (These two viscous force are nearly same)
            if ( laplacianTerm==1 ) then

                !----------------------------------------------------------------------------------------------
                !It is useful for non-free surface flow
                !Reference:  [Morris]  

                RijDotProductDerivates=DOT_PRODUCT(x_ij,Normalized_dwdx_j(k,:))
                Body_Force_Viscous_Acceleration(j,:)=Body_Force_Viscous_Acceleration(j,:)-subdomain_particle_rho(j)*subdomain_particle_volume(i)*2*water_kinematic_viscosity*RijDotProductDerivates/square_distance*v_ij(:)
                
                !----------------------------------------------------------------------------------------------

            elseif ( laplacianTerm==0 ) then

                !----------------------------------------------------------------------------------------------
                !It is useful for free surface flow 
                !Reference: [Theoretical analysis and numerical verification of the consistency of viscous smoothed-particle-hydrodynamics formulations in simulating free-surface flows]
                !           [Adaptive particle refinement and derefinement applied to the smoothed particle hydrodynamics methods]
                Body_Force_Viscous_Acceleration(j,:)=Body_Force_Viscous_Acceleration(j,:)-subdomain_particle_rho(j)*subdomain_particle_volume(i)*2*(dim+2)*water_kinematic_viscosity*v_ij_multiply_x_ij/square_distance*Normalized_dwdx_j(k,:)
                
                !----------------------------------------------------------------------------------------------

            endif

        endif


    enddo
    !***********************************************************************************************************

 end subroutine Calculate_Body_Force