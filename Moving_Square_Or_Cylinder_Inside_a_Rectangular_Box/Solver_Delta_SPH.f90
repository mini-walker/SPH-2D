!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: compute_iterator_information_Delta
!
!  PURPOSE: Compute iterator information for time integration
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

 subroutine Solver_Delta_SPH()
 
    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI

    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,l,v,o,m                                        !Variables for loop

    !SPS Variables
    real(kind=8)::sps_blin,sps_k,sps_niu_t,SPS_Delta_L
    real(kind=8)::master_element_square,last_element_square
    real(kind=8)::ABS_SPS_Strain_Tensor_Square

    real(kind=8)::hxx,hyy,hzz,hxy,hxz,hyz,hvcc,temp_e,rho_ij,temp_p !Temp value
    
    real(kind=8)::average_c                                         !average c(4.68)
    real(kind=8)::average_rho                                       !average rho(4.69)
    real(kind=8),dimension(dim)::x_ij,v_ij                          !position and velocity difference
    real(kind=8)::v_ij_multiply_x_ij                                !x_ij times v_ij
    real(kind=8)::capital_pai_ij                                    !capital_pai
    real(kind=8)::capital_fai_ij_1                                  !4.67
    real(kind=8),dimension(dim)::capital_fai_ij                     !4.67
    real(kind=8)::temp                                              !temp
    real(kind=8)::square_distance                                   !square distance

    !Laminar laplacian term
    real(kind=8)::RijDotProductDerivates                            !Rij*D Wij
    
    !Delta SPH
    real(kind=8)::rho_difference                                                   !rho difference
    real(kind=8),dimension(dim)::psai                                              !psai in Delta-SPH
    real(kind=8),dimension(dim,dim)::i_input_matrix                                !input matrix i
    real(kind=8),dimension(dim,dim)::j_input_matrix                                !input matrix j
    real(kind=8),dimension(dim,dim)::i_matrix_inversion                            !inversion matrix i
    real(kind=8),dimension(dim,dim)::j_matrix_inversion                            !inversion matrix j
    real(kind=8),dimension(dim,dim)::right_hand_side
    real(kind=8),dimension(dim,dim)::A,A_inversion

    real(kind=8),dimension(dim)::dwdx_i,dwdx_j                                     !Normalized derivates of kernel function
    real(kind=8),dimension(dim)::kth_pair_dwdx,kth_pair_dwdx_i,kth_pair_dwdx_j     !Normalized derivates of kernel function

    !External Boundary force
    real(kind=8)::boundary_force_chi
    real(kind=8)::boundary_force_alpha
    real(kind=8)::boundary_force_F_alpha
        
    integer::check_error=0

    !==========================================================================================================
    
    ! Body of subroutine compute_iterator_information

    !**********************************************************************************************************
    ! Get the MPI runing information
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    !**********************************************************************************************************
    

    !**************************************Stand SPH calculation***********************************
    do k=1,pair_number
        
        i=pair_i(k)                        !i index in pair k
        j=pair_j(k)                        !j index in pair k
        
        kth_pair_dwdx(:)=dwdx(k,:)

        !velocity difference
        velocity_difference(:)=velocity_difference_ij(k,:)+average_velocity(i,:)-average_velocity(j,:)  !Attention: i-j
        
        !velocity difference times weight derivates
        temp=DOT_PRODUCT(velocity_difference,kth_pair_dwdx)

        !------------------------------------------------------------------------------------------
        !Stand SPH continum equation
        drhodt(i)=drhodt(i)+subdomain_particle_rho(i)*subdomain_particle_volume(j)*temp
        drhodt(j)=drhodt(j)+subdomain_particle_rho(j)*subdomain_particle_volume(i)*temp
        !------------------------------------------------------------------------------------------
        
        !************************************************************************************
        !Stand SPH moment equation
        rho_ij=1.0/(subdomain_particle_rho(i)*subdomain_particle_rho(j))
        temp_e=0.0d0
            
        !velocity difference
        velocity_difference(:)=-velocity_difference_ij(k,:)    !Attention: j-i
        do L=1,dim
                
           !(pi+pj)/(rho_i*rho_j)*dwdx(i)
           temp_p=-(subdomain_particle_press(i)+subdomain_particle_press(j))*dwdx(k,L)*rho_ij
           temp_e=temp_e+velocity_difference(l)*temp_p
 
          !-------------------------------------------------------------------------------------------------------------
          !Stand SPH acceleration
          internal_acceleration(i,l)=internal_acceleration(i,l)+subdomain_particle_mass(j)*temp_p
          internal_acceleration(j,l)=internal_acceleration(j,l)-subdomain_particle_mass(i)*temp_p
          !-------------------------------------------------------------------------------------------------------------
                
        end do

        dedt(i)=dedt(i)+0.5*subdomain_particle_mass(j)*temp_e
        dedt(j)=dedt(j)+0.5*subdomain_particle_mass(i)*temp_e
        !************************************************************************************
        !artifical viscosity acceleration
        average_smooth_length=average_smooth_length_ij(k)
        
        !velocity_difference*position_difference_ij
        v_ij(:)=velocity_difference_ij(k,:)
        x_ij(:)=position_difference_ij(k,:)

        square_distance=DOT_PRODUCT(x_ij,x_ij)
        v_ij_multiply_x_ij=DOT_PRODUCT(v_ij,x_ij)
        distance=sqrt(square_distance)

        average_c=average_c_ij(k)
        average_rho=average_rho_ij(k)

        !v_ij_multiply_x_ij in Equ.(4.66) in Liu
        if(v_ij_multiply_x_ij<0.0) then
            
            !Equ.(4.67) in Liu
            capital_fai_ij_1=average_smooth_length*v_ij_multiply_x_ij/(square_distance+(fai_coefficient*average_smooth_length)**2)
            
            !Equ.(4.66) in Liu
            capital_pai_ij=(-alpha*average_c*capital_fai_ij_1+beta*(capital_fai_ij_1**2))/average_rho
            
            !Equ.(4.66) in Liu
            do L=1,dim
                temp=-capital_pai_ij*dwdx(k,L)            
                artifical_viscosity_acceleration(i,L)=artifical_viscosity_acceleration(i,L)+subdomain_particle_mass(j)*temp
                artifical_viscosity_acceleration(j,L)=artifical_viscosity_acceleration(j,L)-subdomain_particle_mass(i)*temp
            end do    
        end if
        !************************************************************************************

        !************************************************************************************
        !laplacian viscos term
        
        !------------------------------------------------------------------------------------
        ! !laplacian viscos term
        ! !Reference: Simulation of near-shore solitary wave mechanics by an incompressible SPH method 
        
        ! RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_i)
        ! laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_mass(j)*4*water_kinematic_viscosity*RijDotProductDerivates/(square_distance*(subdomain_particle_rho(i)+subdomain_particle_rho(j)))*v_ij(:)
        
        ! RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_j)
        ! laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_mass(i)*4*water_kinematic_viscosity*RijDotProductDerivates/(square_distance*(subdomain_particle_rho(i)+subdomain_particle_rho(j)))*v_ij(:)
        !------------------------------------------------------------------------------------

        if ( laplacianTerm==1 ) then

            !------------------------------------------------------------------------------------
            !It is useful for non-free surface flow
            !Reference:  [Morris]  
            RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_i)
            laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_volume(j)*2*water_kinematic_viscosity*RijDotProductDerivates/square_distance*v_ij(:)
          
            RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_j)
            laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_volume(i)*2*water_kinematic_viscosity*RijDotProductDerivates/square_distance*v_ij(:)
            !------------------------------------------------------------------------------------

        elseif ( laplacianTerm==0 ) then

            !------------------------------------------------------------------------------------
            !It is useful for free surface flow 
            !Reference: [Theoretical analysis and numerical verification of the consistency of viscous smoothed-particle-hydrodynamics formulations in simulating free-surface flows]
            !           [Adaptive particle refinement and derefinement applied to the smoothed particle hydrodynamics methods]

            laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_volume(j)*2*(dim+2)*water_kinematic_viscosity*v_ij_multiply_x_ij/square_distance*kth_pair_dwdx_i(:)

            laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_volume(i)*2*(dim+2)*water_kinematic_viscosity*v_ij_multiply_x_ij/square_distance*kth_pair_dwdx_j(:)
            !------------------------------------------------------------------------------------

        else

            write(*,*) 'The defination of laplacianTerm is not right!(Solver)'
            
        endif
        
        !************************************************************************************

        !************************************************************************************
        !artifical viscosity energy
        temp=-DOT_PRODUCT(v_ij,kth_pair_dwdx)

        artifical_viscosity_dedt(i)=artifical_viscosity_dedt(i)+0.5*Delta_alpha*average_smooth_length*average_c*water_rho_0/subdomain_particle_rho(i)*v_ij_multiply_x_ij/square_distance*subdomain_particle_volume(j)*temp
        artifical_viscosity_dedt(j)=artifical_viscosity_dedt(j)+0.5*Delta_alpha*average_smooth_length*average_c*water_rho_0/subdomain_particle_rho(j)*v_ij_multiply_x_ij/square_distance*subdomain_particle_volume(i)*temp

        !********************************************************************************
        !artifical internal acceleration in Delta-SPH(out)
        temp=Delta_alpha*average_smooth_length*average_c*v_ij_multiply_x_ij/square_distance
        do L=1,dim
           artifical_internal_acceleration(i,L)=artifical_internal_acceleration(i,L)+temp*dwdx(k,L)*water_rho_0/subdomain_particle_rho(i)*subdomain_particle_volume(j)
           artifical_internal_acceleration(j,L)=artifical_internal_acceleration(j,L)-temp*dwdx(k,L)*water_rho_0/subdomain_particle_rho(j)*subdomain_particle_volume(i)
        end do
        
        !capital_fai_ij

        x_ij=-position_difference_ij(k,:)                             !Attention: j-i
        
        capital_fai_ij=0.0d0
        do l=1,dim
            capital_fai_ij(l)=2*(subdomain_particle_energy(j)-subdomain_particle_energy(i))*x_ij(l)/square_distance
        end do
        
        temp=DOT_PRODUCT(capital_fai_ij,kth_pair_dwdx)

        temp=chi*temp*average_smooth_length*average_c*water_rho_0/average_rho
        
        artifical_heat_dedt(i)=artifical_heat_dedt(i)+temp*subdomain_particle_volume(j)
        artifical_heat_dedt(j)=artifical_heat_dedt(j)-temp*subdomain_particle_volume(i)

        !************************************************************************************

        !************************************************************************************

        velocity_difference(:)=velocity_difference_ij(k,:) ! particle i-j
        
        do L=1,dim
            average_velocity(i,L)=average_velocity(i,L)-epsilon*subdomain_particle_mass(j)*velocity_difference(L)/(subdomain_particle_rho(i)+subdomain_particle_rho(j))*w(k)
            average_velocity(j,L)=average_velocity(j,L)+epsilon*subdomain_particle_mass(i)*velocity_difference(L)/(subdomain_particle_rho(i)+subdomain_particle_rho(j))*w(k)                   
        end do
        !************************************************************************************

        ! !************************************************************************************
        ! !Calculate the force from the boundary particles
        ! if (subdomain_particle_type(j)==-3 .or. subdomain_particle_type(i)==-3) then

        !     position_difference(:)=position_difference_ij(k,:) !Attention i-j

        !     if (distance<average_smooth_length) then

        !         boundary_force_chi=1-distance/average_smooth_length
        !         boundary_force_alpha=distance/(0.5*average_smooth_length)

        !         if (boundary_force_alpha<=2.0/3.0) then

        !             boundary_force_F_alpha=2.0/3.0

        !         elseif(boundary_force_alpha>2.0/3.0 .and. boundary_force_alpha<=1.0) then

        !             boundary_force_F_alpha=2*boundary_force_alpha-1.5*boundary_force_alpha**2

        !         elseif(boundary_force_alpha>1.0 .and. boundary_force_alpha<=2.0) then

        !             boundary_force_F_alpha=0.5*(2-boundary_force_alpha)**2

        !         else

        !             boundary_force_F_alpha=0.0d0
                    
        !         endif

        !         !Only for the fluid particles
        !         if (subdomain_particle_type(i)/=-3) then
        !             external_acceleration(i,:)=external_acceleration(i,:)+0.01*square_c_0*boundary_force_chi*boundary_force_F_alpha/square_distance*position_difference(:)*subdomain_particle_mass(j)/subdomain_particle_mass(i)
        !         endif

        !         if (subdomain_particle_type(j)/=-3) then
        !             external_acceleration(j,:)=external_acceleration(j,:)-0.01*square_c_0*boundary_force_chi*boundary_force_F_alpha/square_distance*position_difference(:)*subdomain_particle_mass(i)/subdomain_particle_mass(j)
        !         endif
                
        !     endif
            
        ! endif
        ! !************************************************************************************

         
    end do
    !******************************************************************************************

    
    !*******************************Delta-SPH Density coorection*******************************
    !-----------------------------------------------------------------------------------------
    grad_rho=0.0d0
    do k=1,pair_number
       
       i=pair_i(k)                        !i index in pair k
       j=pair_j(k)                        !j index in pair k

       !Density difference
       rho_difference=-rho_difference_ij(k)  !Attention: j-i

       ! inversion matrix
       i_matrix_inversion=Tensor_dyadic_inverse_matrix(i,:,:)
       j_matrix_inversion=Tensor_dyadic_inverse_matrix(j,:,:)

       !Density gradience
       grad_rho(i,:)=grad_rho(i,:)+rho_difference*subdomain_particle_volume(j)*matmul(i_matrix_inversion,dwdx(k,:))
       grad_rho(j,:)=grad_rho(j,:)+rho_difference*subdomain_particle_volume(i)*matmul(j_matrix_inversion,dwdx(k,:))

    end do  
    !-----------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------
    !Exchange the grad_rho form the neighbour processor

    call MPI_RealVector_MessageExchange(grad_rho)           ! Exchange the gradience of the density
    
    !-----------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------
    do k=1,pair_number
      
       i=pair_i(k)                        !particle i index in pair k
       j=pair_j(k)                        !particle j index in pair k

       kth_pair_dwdx(:)=dwdx(k,:)

       !Density difference
       rho_difference=-rho_difference_ij(k)  !Attention: j-i
       
       average_smooth_length=average_smooth_length_ij(k)
       average_c=average_c_ij(k)
       
       !distance and position difference
       position_difference(:)=-position_difference_ij(k,:)                  !Attention: j-i
       
       square_distance=DOT_PRODUCT(position_difference,position_difference)

       !psai 
       psai(:)=2*rho_difference/square_distance*position_difference(:)-(grad_rho(i,:)+grad_rho(j,:))
       
    
       !Temp Value
       temp=DOT_PRODUCT(psai,kth_pair_dwdx)
       
       artifical_rho_correct(i)=artifical_rho_correct(i)+delta*average_smooth_length*average_c*temp*subdomain_particle_volume(j)
       artifical_rho_correct(j)=artifical_rho_correct(j)-delta*average_smooth_length*average_c*temp*subdomain_particle_volume(i)
       
    end do
    !-----------------------------------------------------------------------------------------
    
    !***********************************************************************************************
    
    ! !******************************A simplify Delta SPH correction*********************************
    ! !This useful for the No free surface flow
    ! do k=1,pair_number
        
    !    i=pair_i(k)                        !i index in pair k
    !    j=pair_j(k)                        !j index in pair k

    !    kth_pair_dwdx(:)=dwdx(k,:)
              
    !    !Density difference
    !    rho_difference=-rho_difference_ij(k)  !Attention: j-i
       
    !    average_smooth_length=average_smooth_length_ij(k)
    !    average_c=average_c_ij(k)
         
    !    !distance and position difference
    !    position_difference(:)=-position_difference_ij(k,:)                  !Attention: j-i
       
    !    square_distance=DOT_PRODUCT(position_difference,position_difference)
    
    !     !psai
    !     do L=1,dim
    !         psai(L)=2*rho_difference*position_difference(L)/(square_distance+(fai_coefficient*average_smooth_length)**2)
    !     end do
    
    !    !Temp Value
    !    temp=DOT_PRODUCT(psai,kth_pair_dwdx)
    
    !     artifical_rho_correct(i)=artifical_rho_correct(i)+delta*average_smooth_length*average_c*temp*subdomain_particle_volume(j)
    !     artifical_rho_correct(j)=artifical_rho_correct(j)-delta*average_smooth_length*average_c*temp*subdomain_particle_volume(i)
    
    ! end do
    ! !********************************************************************************************
    
    !Final Density change rate
    do i=1,actual_particle_number_in_subdomain
        
        if(subdomain_particle_type(i)/=2) cycle                !Only for the fluid particles
        drhodt(i)=drhodt(i)+artifical_rho_correct(i)

    end do
    
    
 
 end subroutine Solver_Delta_SPH