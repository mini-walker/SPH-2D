!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: Normalized_Linear_Riemann_Delta_Solver
!
!  PURPOSE: Compute iterator information for time integration with Normalized kernel function Delta SPH and SPS Riemann solver
!
!  Programer: Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location: MUN
!
!  Time: 2017.3.18
!
!  Copyright: Memorial University
!
!  Version: 1.0
!
!  Note: MPI version: mpich-3.2
!        Fortran version: Inter Fortran (ifort) 18.0.0
!****************************************************************************

 subroutine Solver_Normalized_Modified_Linear_Riemann_Delta()
 
    use information_module
    use function_module
    use Public_variable_module

    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,l,v,o,m                                        !Variables for loop
    real(kind=8)::distance                                        !distance between i and j
    real(kind=8)::average_smooth_length                           !average smooth length
    real(kind=8),dimension(dim)::position_difference              !position difference

    !SPS Variables
    real(kind=8)::sps_blin,sps_k,sps_niu_t,SPS_Delta_L
    real(kind=8)::master_element_square,last_element_square
    real(kind=8)::ABS_SPS_Strain_Tensor_Square

    real(kind=8),dimension(dim)::velocity_difference                !velocity difference
    real(kind=8)::hxx,hyy,hzz,hxy,hxz,hyz,hvcc,rho_ij               !
    
    !artifical viscosity term
    real(kind=8)::average_c                                         !average c(4.68)
    real(kind=8)::average_rho                                       !average rho(4.69)
    real(kind=8),dimension(dim)::x_ij,v_ij                          !position and velocity difference
    real(kind=8)::v_ij_multiply_x_ij                                !x_ij times v_ij
    real(kind=8)::capital_pai_ij                                    !capital_pai
    real(kind=8)::capital_fai_ij_1                                  !4.67
    real(kind=8),dimension(dim)::capital_fai_ij                     !4.67
    real(kind=8)::temp,temp_i,temp_j                                !temp value
    real(kind=8)::square_distance                                   !square distance
    real(kind=8)::temp_e,temp_e_i,temp_e_j
    real(kind=8)::temp_p,temp_p_i,temp_p_j

    !Laminar laplacian term
    real(kind=8)::RijDotProductDerivates                            !Rij*D Wij
    
    !Delta SPH
    real(kind=8)::rho_difference                                                   !rho difference
    real(kind=8),dimension(dim)::psai                                              !psai in Delta-SPH
    real(kind=8)::psai_ij                                                          !psai in modified Delta-SPH
    real(kind=8),dimension(dim)::sum_grad_rho_ij                                   !grad_rho_i+grad_rho_j
    real(kind=8),dimension(dim,dim)::i_input_matrix                                !input matrix i
    real(kind=8),dimension(dim,dim)::j_input_matrix                                !input matrix j
    real(kind=8),dimension(dim,dim)::i_matrix_inversion                            !inversion matrix i
    real(kind=8),dimension(dim,dim)::j_matrix_inversion                            !inversion matrix j
    real(kind=8),dimension(dim)::dwdx_right_hand_side

    real(kind=8),dimension(dim)::dwdx_i,dwdx_j                                     !Normalized derivates of kernel function
    real(kind=8),dimension(dim)::kth_pair_dwdx,kth_pair_dwdx_i,kth_pair_dwdx_j     !Normalized derivates of kernel function

    !External Boundary force
    real(kind=8)::boundary_force_chi
    real(kind=8)::boundary_force_alpha
    real(kind=8)::boundary_force_F_alpha

    !Riemann slover 
    real(kind=8)::P_Left,P_right,rho_Left,rho_right,U_Left,U_right
    real(kind=8)::P_Bar,U_Bar                                                      !Average pressure and velocity
    real(kind=8),dimension(dim)::I_Velocity,J_Velocity,IJ_Velocity_Bar
    real(kind=8),dimension(dim)::IJ_Direction_Vector
    real(kind=8)::U_asterisk,P_asterisk,I_P_asterisk,J_P_asterisk
    real(kind=8),dimension(dim)::velocity_asterisk,I_velocity_asterisk,J_velocity_asterisk
    real(kind=8)::Belta_Riemman

        
    integer::check_error=0

    !==========================================================================================================
    !DOT_PRODUCT(vector_a,vector_b) two vector DOT_PRODUCT
    !MATMUL(matrix_a,matrix_b)two matrix DOT_PRODUCT

    ! Body of subroutine compute_iterator_information

    !**********************************************************************************************************
    ! Get the MPI runing information
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    !**********************************************************************************************************
    

    !**************************************************************************************** 
    !Calculate the normalized derivates
    Normalized_dwdx_i=0.0d0
    Normalized_dwdx_j=0.0d0

    do k=1,pair_number
    
       i=pair_i(k)                        !particle i index in kth pair
       j=pair_j(k)                        !particle j index in kth pair

       i_matrix_inversion=Tensor_dyadic_inverse_matrix(i,:,:)
       j_matrix_inversion=Tensor_dyadic_inverse_matrix(j,:,:)

       dwdx_right_hand_side(:)=dwdx(k,:)

       Normalized_dwdx_i(k,:)=matmul(i_matrix_inversion,dwdx_right_hand_side)
       Normalized_dwdx_j(k,:)=matmul(j_matrix_inversion,dwdx_right_hand_side)

    end do
    !**************************************************************************************** 


    !*******************************Stand SPH calculation*********************************
    do k=1,pair_number
        
        i=pair_i(k)                                     !particle i index in kth pair
        j=pair_j(k)                                     !particle j index in kth pair

        !The pressure and velocity in Riemann problem
        P_Left=subdomain_particle_press(i)
        P_Right=subdomain_particle_press(j)
        P_Bar=(P_Left+P_Right)/2.0

        I_Velocity=subdomain_particle_velocity(i,:)
        J_Velocity=subdomain_particle_velocity(j,:)
        IJ_Velocity_Bar=(I_Velocity+J_Velocity)/2.0


        !Get the direction vector for from i to j
        x_ij(:)=position_difference_ij(k,:)
        distance=sqrt(DOT_PRODUCT(x_ij,x_ij))
        IJ_Direction_Vector=-x_ij/distance

        U_Left=DOT_PRODUCT(I_Velocity,IJ_Direction_Vector)
        U_Right=DOT_PRODUCT(J_Velocity,IJ_Direction_Vector)
        U_Bar=(U_Left+U_Right)/2.0

        average_c=average_c_ij(k)
        average_rho=average_rho_ij(k)

        U_asterisk=U_Bar+0.5*(P_Left-P_Right)/(average_rho*sqrt(square_c_0))

        velocity_asterisk=(U_asterisk-U_Bar)*IJ_Direction_Vector+IJ_Velocity_Bar

        Belta_Riemman=min(3.0*(max((U_Left-U_Right),0.0d0)),average_c)
        P_asterisk=P_Bar+0.5*Belta_Riemman*average_rho*(U_Left-U_Right)


        !The normalized derivates of weight function
        kth_pair_dwdx_i(:)=Normalized_dwdx_i(k,:)
        kth_pair_dwdx_j(:)=Normalized_dwdx_j(k,:)

        !velocity difference times weight derivates
        !-------------------------------------------------------------------------------------------------------------
        !For the particles near the free surface
        if (subdomain_near_free_surface_or_not(i)==1) then
            I_velocity_asterisk=IJ_Velocity_Bar
            I_P_asterisk=P_Bar
        else
            I_velocity_asterisk=velocity_asterisk
            I_P_asterisk=P_asterisk
        end if

        if (subdomain_near_free_surface_or_not(j)==1) then
            J_velocity_asterisk=IJ_Velocity_Bar
            J_P_asterisk=P_Bar
        else
            J_velocity_asterisk=velocity_asterisk
            J_P_asterisk=P_asterisk
        end if

        !write(*,*) I_velocity_asterisk,J_velocity_asterisk

        velocity_difference(:)=I_Velocity-I_velocity_asterisk     !Attention: i-j
        temp_i=DOT_PRODUCT(velocity_difference,kth_pair_dwdx_i)
        velocity_difference(:)=J_Velocity-J_velocity_asterisk     !Attention: i-j
        temp_j=-DOT_PRODUCT(velocity_difference,kth_pair_dwdx_j)
        !-------------------------------------------------------------------------------------------------------------

        !-------------------------------------------------------------------------------------------------------------
        !Stand SPH continum equation
        drhodt(i)=drhodt(i)+2*subdomain_particle_rho(i)*subdomain_particle_volume(j)*temp_i
        drhodt(j)=drhodt(j)+2*subdomain_particle_rho(j)*subdomain_particle_volume(i)*temp_j
        !-------------------------------------------------------------------------------------------------------------
        
        !************************************************************************************
        !Stand SPH moment equation
        rho_ij=1.0/(subdomain_particle_rho(i)*subdomain_particle_rho(j))

        do L=1,dim
                
          !(pi+pj)/(rho_i*rho_j)*dwdx(i)
          temp_p_i=-2*I_P_asterisk*Normalized_dwdx_i(k,L)*rho_ij
          temp_p_j=-2*J_P_asterisk*Normalized_dwdx_j(k,L)*rho_ij
    
          !-------------------------------------------------------------------------------------------------------------
          !Stand SPH acceleration
          internal_acceleration(i,L)=internal_acceleration(i,L)+subdomain_particle_mass(j)*temp_p_i
          internal_acceleration(j,L)=internal_acceleration(j,L)-subdomain_particle_mass(i)*temp_p_j
          !-------------------------------------------------------------------------------------------------------------
                
        end do

        !************************************************************************************
        !artifical viscosity acceleration
        average_smooth_length=average_smooth_length_ij(k)
        
        !velocity_difference*position_difference_ij
        v_ij(:)=velocity_difference_ij(k,:)
        x_ij(:)=position_difference_ij(k,:)

        square_distance=DOT_PRODUCT(x_ij,x_ij)
        v_ij_multiply_x_ij=DOT_PRODUCT(v_ij,x_ij)
        distance=sqrt(square_distance)

        !v_ij_multiply_x_ij in Equ.(4.66) in Liu
        if(v_ij_multiply_x_ij<0.0) then
            
            !Equ.(4.67) in Liu
            capital_fai_ij_1=average_smooth_length*v_ij_multiply_x_ij/(square_distance+(fai_coefficient*average_smooth_length)**2)
            
            !Equ.(4.66) in Liu
            capital_pai_ij=(-alpha*average_c*capital_fai_ij_1+beta*(capital_fai_ij_1**2))/average_rho
            
            do L=1,dim
                temp_i=-capital_pai_ij*Normalized_dwdx_i(k,L)  
                temp_j=-capital_pai_ij*Normalized_dwdx_j(k,L)           
                artifical_viscosity_acceleration(i,L)=artifical_viscosity_acceleration(i,L)+subdomain_particle_mass(j)*temp_i
                artifical_viscosity_acceleration(j,L)=artifical_viscosity_acceleration(j,L)-subdomain_particle_mass(i)*temp_j
            end do    
        end if
        !************************************************************************************

        !************************************************************************************
        ! !laplacian viscos term
        ! !Reference: Simulation of near-shore solitary wave mechanics by an incompressible SPH method []
        
        ! RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_i)
        ! laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_mass(j)*4*water_kinematic_viscosity*RijDotProductDerivates/(square_distance*(subdomain_particle_rho(i)+subdomain_particle_rho(j)))*v_ij(:)
        
        ! RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_j)
        ! laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_mass(i)*4*water_kinematic_viscosity*RijDotProductDerivates/(square_distance*(subdomain_particle_rho(i)+subdomain_particle_rho(j)))*v_ij(:)


        !Reference:  [Morris]  
        RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_i)
        laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_volume(j)*2*water_kinematic_viscosity*RijDotProductDerivates/square_distance*v_ij(:)
      
        RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_j)
        laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_volume(i)*2*water_kinematic_viscosity*RijDotProductDerivates/square_distance*v_ij(:)

        !************************************************************************************

        !************************************************************************************
        !artifical internal acceleration in Delta-SPH(out)
        temp=Delta_alpha*average_smooth_length*average_c*v_ij_multiply_x_ij/square_distance
        do L=1,dim
           artifical_internal_acceleration(i,L)=artifical_internal_acceleration(i,L)+temp*Normalized_dwdx_i(k,L)*water_rho_0/subdomain_particle_rho(i)*subdomain_particle_volume(j)
           artifical_internal_acceleration(j,L)=artifical_internal_acceleration(j,L)-temp*Normalized_dwdx_j(k,L)*water_rho_0/subdomain_particle_rho(j)*subdomain_particle_volume(i)
        end do
        
        !************************************************************************************
        velocity_difference(:)=velocity_difference_ij(k,:) ! particle i-j
        
        do L=1,dim
            average_velocity(i,L)=average_velocity(i,L)-epsilon*subdomain_particle_mass(j)*velocity_difference(L)/(subdomain_particle_rho(i)+subdomain_particle_rho(j))*w(k)
            average_velocity(j,L)=average_velocity(j,L)+epsilon*subdomain_particle_mass(i)*velocity_difference(L)/(subdomain_particle_rho(i)+subdomain_particle_rho(j))*w(k)                   
        end do
        !************************************************************************************
   
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

    do j=1,dim

        left_Message_ID=300
        right_Message_ID=301

        !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
        !Data form left processor to right processor
        call MPI_SENDRECV(grad_rho(message_first_particle_index_to_right,j),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                          grad_rho(message_first_particle_index_from_left,j),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
        
        !Data form right processor to left processor
        call MPI_SENDRECV(grad_rho(message_first_particle_index_to_left,j),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                          grad_rho(message_first_particle_index_from_right,j),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    end do
    !-----------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------
    !psai=(rho_j-rho_i)-0.5*(grad_rho_i+grad_rho_j)*(r_j-r_i)
    !artifical_rho_correct=2*delta*h*c*psai*(r_j-r_i)*dwdx*V_j/square_distance
    do k=1,pair_number
      
       i=pair_i(k)                        !i index in pair k
       j=pair_j(k)                        !j index in pair k

       !The normalized derivates of weight function
       kth_pair_dwdx_i(:)=Normalized_dwdx_i(k,:)
       kth_pair_dwdx_j(:)=Normalized_dwdx_j(k,:)
           
       !Density difference
       rho_difference=-rho_difference_ij(k)  !Attention: j-i
       
       average_smooth_length=average_smooth_length_ij(k)
       average_c=average_c_ij(k)
       
       !distance and position difference
       position_difference(:)=-position_difference_ij(k,:)                  !Attention: j-i
       
       square_distance=DOT_PRODUCT(position_difference,position_difference)  
       
       !psai 
       sum_grad_rho_ij(:)=grad_rho(i,:)+grad_rho(j,:)
       psai_ij=rho_difference-0.5*DOT_PRODUCT(sum_grad_rho_ij,position_difference)

       !Temp Value
       temp_i=DOT_PRODUCT(position_difference,kth_pair_dwdx_i)/square_distance
       temp_j=DOT_PRODUCT(position_difference,kth_pair_dwdx_j)/square_distance
       
       artifical_rho_correct(i)=artifical_rho_correct(i)+2*delta*average_smooth_length*average_c*psai_ij*temp_i*subdomain_particle_volume(j)
       artifical_rho_correct(j)=artifical_rho_correct(j)-2*delta*average_smooth_length*average_c*psai_ij*temp_j*subdomain_particle_volume(i)
       
    end do
    !-----------------------------------------------------------------------------------------
    
    !**********************************************************************************************
    
    ! !******************************A simplify Delta SPH correction*********************************
    ! !This useful for the No free surface flow
    ! do k=1,pair_number
        
    !    i=pair_i(k)                        !i index in pair k
    !    j=pair_j(k)                        !j index in pair k

    !    !The normalized derivates of weight function
    !    kth_pair_dwdx_i(:)=Normalized_dwdx_i(k,:)
    !    kth_pair_dwdx_j(:)=Normalized_dwdx_j(k,:)
              
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
    !    temp_i=DOT_PRODUCT(psai,kth_pair_dwdx_i)
    !    temp_j=DOT_PRODUCT(psai,kth_pair_dwdx_j)
       
    !    artifical_rho_correct(i)=artifical_rho_correct(i)+delta*average_smooth_length*average_c*temp_i*subdomain_particle_volume(j)
    !    artifical_rho_correct(j)=artifical_rho_correct(j)-delta*average_smooth_length*average_c*temp_j*subdomain_particle_volume(i)
       

    ! end do
    ! !********************************************************************************************
    
    !Final Density change rate
    do i=1,actual_particle_number_in_subdomain
        
        if(subdomain_particle_type(i)/=2) cycle                !Only for the fluid particles

        !if (subdomain_near_free_surface_or_not(i)==1) cycle    !Delta-SPH only used for the interior particles 
        
        drhodt(i)=drhodt(i)+artifical_rho_correct(i)

    end do
    
    
 
 end subroutine Solver_Normalized_Modified_Linear_Riemann_Delta