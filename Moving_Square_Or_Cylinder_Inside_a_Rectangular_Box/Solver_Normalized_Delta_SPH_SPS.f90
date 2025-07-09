!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: Solver_Normalized_Delta_SPH_SPS
!
!  PURPOSE: Compute iterator information for time integration with Normalized kernel function Delta SPH and SPS
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

 subroutine Solver_Normalized_Delta_SPH_SPS()
 
    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI

    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,l,v,o,m                                          !Variables for loop

    !SPS Variables
    real(kind=8)::sps_blin,sps_k,sps_niu_t,SPS_Delta_L
    real(kind=8)::master_element_square,last_element_square
    real(kind=8)::ABS_SPS_Strain_Tensor_Square

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
    
    !initial the shear viscos terms
    tao_xx=0.0d0
    tao_yy=0.0d0
    tao_zz=0.0d0
    tao_xy=0.0d0
    tao_xz=0.0d0
    tao_yz=0.0d0

    SPS_Strain_Tensor_xx=0.0d0  !SPS Strain Tensor
    SPS_Strain_Tensor_yy=0.0d0
    SPS_Strain_Tensor_zz=0.0d0
    SPS_Strain_Tensor_xy=0.0d0
    SPS_Strain_Tensor_xz=0.0d0
    SPS_Strain_Tensor_yz=0.0d0 

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

    !**************************************************************************************** 
    ! Calculate the SPS Strain Tensor
    do k=1,pair_number
    
       i=pair_i(k)                        !particle i index in kth pair
       j=pair_j(k)                        !particle j index in kth pair 
       
       !velocity difference
        velocity_difference(:)=-velocity_difference_ij(k,:)     !Attention: j-i
       
       !SPS Strain Tensor
       if(dim==1) then
           SPS_Strain_Tensor_xx(i)=SPS_Strain_Tensor_xx(i)+velocity_difference(1)*Normalized_dwdx_i(k,1)*subdomain_particle_volume(j)
           SPS_Strain_Tensor_xx(j)=SPS_Strain_Tensor_xx(j)+velocity_difference(1)*Normalized_dwdx_j(k,1)*subdomain_particle_volume(i)
       else if(dim==2) then
           SPS_Strain_Tensor_xx(i)=SPS_Strain_Tensor_xx(i)+velocity_difference(1)*Normalized_dwdx_i(k,1)*subdomain_particle_volume(j)
           SPS_Strain_Tensor_xy(i)=SPS_Strain_Tensor_xy(i)+(velocity_difference(1)*Normalized_dwdx_i(k,2)+velocity_difference(2)*Normalized_dwdx_i(k,1))*subdomain_particle_volume(j)
           SPS_Strain_Tensor_yy(i)=SPS_Strain_Tensor_yy(i)+velocity_difference(2)*Normalized_dwdx_i(k,2)*subdomain_particle_volume(j)

           SPS_Strain_Tensor_xx(j)=SPS_Strain_Tensor_xx(j)+velocity_difference(1)*Normalized_dwdx_j(k,1)*subdomain_particle_volume(i)
           SPS_Strain_Tensor_xy(j)=SPS_Strain_Tensor_xy(j)+(velocity_difference(1)*Normalized_dwdx_j(k,2)+velocity_difference(2)*Normalized_dwdx_j(k,1))*subdomain_particle_volume(i)
           SPS_Strain_Tensor_yy(j)=SPS_Strain_Tensor_yy(j)+velocity_difference(2)*Normalized_dwdx_j(k,2)*subdomain_particle_volume(i)
       else if(dim==3) then
           SPS_Strain_Tensor_xx(i)=SPS_Strain_Tensor_xx(i)+velocity_difference(1)*Normalized_dwdx_i(k,1)*subdomain_particle_volume(j)
           SPS_Strain_Tensor_xy(i)=SPS_Strain_Tensor_xy(i)+(velocity_difference(1)*Normalized_dwdx_i(k,2)+velocity_difference(2)*Normalized_dwdx_i(k,1))*subdomain_particle_volume(j)
           SPS_Strain_Tensor_xz(i)=SPS_Strain_Tensor_xz(i)+(velocity_difference(1)*Normalized_dwdx_i(k,3)+velocity_difference(3)*Normalized_dwdx_i(k,1))*subdomain_particle_volume(j)
           SPS_Strain_Tensor_yy(i)=SPS_Strain_Tensor_yy(i)+velocity_difference(2)*Normalized_dwdx_i(k,2)*subdomain_particle_volume(j)
           SPS_Strain_Tensor_yz(i)=SPS_Strain_Tensor_yz(i)+(velocity_difference(2)*Normalized_dwdx_i(k,3)+velocity_difference(3)*Normalized_dwdx_i(k,2))*subdomain_particle_volume(j)
           SPS_Strain_Tensor_zz(i)=SPS_Strain_Tensor_zz(i)+velocity_difference(3)*Normalized_dwdx_i(k,3)*subdomain_particle_volume(j)

           SPS_Strain_Tensor_xx(j)=SPS_Strain_Tensor_xx(j)+velocity_difference(1)*Normalized_dwdx_j(k,1)*subdomain_particle_volume(i)
           SPS_Strain_Tensor_xy(j)=SPS_Strain_Tensor_xy(j)+(velocity_difference(1)*Normalized_dwdx_j(k,2)+velocity_difference(2)*Normalized_dwdx_j(k,1))*subdomain_particle_volume(i)
           SPS_Strain_Tensor_xz(j)=SPS_Strain_Tensor_xz(j)+(velocity_difference(1)*Normalized_dwdx_j(k,3)+velocity_difference(3)*Normalized_dwdx_j(k,1))*subdomain_particle_volume(i)
           SPS_Strain_Tensor_yy(j)=SPS_Strain_Tensor_yy(j)+velocity_difference(2)*Normalized_dwdx_j(k,2)*subdomain_particle_volume(i)
           SPS_Strain_Tensor_yz(j)=SPS_Strain_Tensor_yz(j)+(velocity_difference(2)*Normalized_dwdx_j(k,3)+velocity_difference(3)*Normalized_dwdx_j(k,2))*subdomain_particle_volume(i)
           SPS_Strain_Tensor_zz(j)=SPS_Strain_Tensor_zz(j)+velocity_difference(3)*Normalized_dwdx_j(k,3)*subdomain_particle_volume(i)
       else
           write(*,*) "The dim is not right!"
       end if
    
    end do
    !******************************************************************************************

    !-----------------------------------------------------------------------------------------
    !Exchange the SPS_Strain_Tensor form the neighbour processor
    left_Message_ID=300
    right_Message_ID=301

    !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
    !Data form left processor to right processor
    call MPI_SENDRECV(SPS_Strain_Tensor_xx(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      SPS_Strain_Tensor_xx(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(SPS_Strain_Tensor_xx(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      SPS_Strain_Tensor_xx(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    left_Message_ID=302
    right_Message_ID=303

    !Data form left processor to right processor
    call MPI_SENDRECV(SPS_Strain_Tensor_xy(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      SPS_Strain_Tensor_xy(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(SPS_Strain_Tensor_xy(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      SPS_Strain_Tensor_xy(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    left_Message_ID=304
    right_Message_ID=305

    !Data form left processor to right processor
    call MPI_SENDRECV(SPS_Strain_Tensor_xz(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      SPS_Strain_Tensor_xz(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(SPS_Strain_Tensor_xz(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      SPS_Strain_Tensor_xz(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    left_Message_ID=306
    right_Message_ID=307

    !Data form left processor to right processor
    call MPI_SENDRECV(SPS_Strain_Tensor_yy(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      SPS_Strain_Tensor_yy(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(SPS_Strain_Tensor_yy(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      SPS_Strain_Tensor_yy(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    left_Message_ID=308
    right_Message_ID=309

    !Data form right processor to left processor
    call MPI_SENDRECV(SPS_Strain_Tensor_yz(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      SPS_Strain_Tensor_yz(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    !Data form left processor to right processor
    call MPI_SENDRECV(SPS_Strain_Tensor_yz(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      SPS_Strain_Tensor_yz(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    left_Message_ID=310
    right_Message_ID=311

    !Data form right processor to left processor
    call MPI_SENDRECV(SPS_Strain_Tensor_zz(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      SPS_Strain_Tensor_zz(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    !Data form left processor to right processor
    call MPI_SENDRECV(SPS_Strain_Tensor_zz(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      SPS_Strain_Tensor_zz(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !-----------------------------------------------------------------------------------------

     do i=1,total_particle_number_in_subdomain

        if(subdomain_particle_type(i)/=2) cycle                !Only for the fluid particles
         
        master_element_square=SPS_Strain_Tensor_xx(i)**2+SPS_Strain_Tensor_yy(i)**2+SPS_Strain_Tensor_zz(i)**2
        last_element_square=SPS_Strain_Tensor_xy(i)**2+SPS_Strain_Tensor_xz(i)**2+SPS_Strain_Tensor_yz(i)**2
        ABS_SPS_Strain_Tensor_Square=2*master_element_square+last_element_square

        if (dim==2) then

            SPS_Delta_L=sqrt(2.0)/2.0*interior_dx

        elseif (dim==3) then

            SPS_Delta_L=sqrt(3.0)/3.0*interior_dx
            
        endif

        sps_niu_t=(0.12*SPS_Delta_L)**2*sqrt(ABS_SPS_Strain_Tensor_Square)

        if (dim==2) then
 
            sps_k=2.0/3.0*sps_niu_t*(SPS_Strain_Tensor_xx(i)+SPS_Strain_Tensor_yy(i))
            
        elseif (dim==3) then

            sps_k=2.0/3.0*sps_niu_t*(SPS_Strain_Tensor_xx(i)+SPS_Strain_Tensor_yy(i)+SPS_Strain_Tensor_zz(i))
            
        endif

        sps_blin=2.0/3.0*0.0066*(SPS_Delta_L)**2*ABS_SPS_Strain_Tensor_Square

         
        tao_xx(i)=(2*sps_niu_t*SPS_Strain_Tensor_xx(i)-sps_k-sps_blin)*subdomain_particle_rho(i)
        tao_xy(i)=(sps_niu_t*SPS_Strain_Tensor_xy(i))*subdomain_particle_rho(i)
        tao_xz(i)=(sps_niu_t*SPS_Strain_Tensor_xz(i))*subdomain_particle_rho(i)
        tao_yy(i)=(2*sps_niu_t*SPS_Strain_Tensor_yy(i)-sps_k-sps_blin)*subdomain_particle_rho(i)
        tao_yz(i)=(sps_niu_t*SPS_Strain_Tensor_yz(i))*subdomain_particle_rho(i)
        tao_zz(i)=(2*sps_niu_t*SPS_Strain_Tensor_zz(i)-sps_k-sps_blin)*subdomain_particle_rho(i)

     enddo
     !********************************************************************************************

    ! !*****************************************XSPH correction***********************************
    ! ! Calculate XSPH term (To make the simulation stable)
    ! do k=1,pair_number
        
    !     i=pair_i(k)                        !particle i index in pair k
    !     j=pair_j(k)                        !particle j index in pair k

    !     velocity_difference(:)=velocity_difference_ij(k,:) ! particle i-j
        
    !     do L=1,dim
    !         average_velocity(i,L)=average_velocity(i,L)-epsilon*subdomain_particle_mass(j)*velocity_difference(L)/(subdomain_particle_rho(i)+subdomain_particle_rho(j))*w(k)
    !         average_velocity(j,L)=average_velocity(j,L)+epsilon*subdomain_particle_mass(i)*velocity_difference(L)/(subdomain_particle_rho(i)+subdomain_particle_rho(j))*w(k)                   
    !     end do

    ! end do
    ! !----------------------------------------------------------------------------------------------

    ! !----------------------------------------------------------------------------------------------
    ! !Exchange the XSPH term
    ! do j=1,dim
      
    !       left_Message_ID=650+j
    !       right_Message_ID=750+j

    !       !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
    !       !Data form left processor to right processor
    !       call MPI_SENDRECV(average_velocity(message_first_particle_index_to_right,j),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
    !                         average_velocity(message_first_particle_index_from_left,j),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
          
    !       !Data form right processor to left processor
    !       call MPI_SENDRECV(average_velocity(message_first_particle_index_to_left,j),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
    !                         average_velocity(message_first_particle_index_from_right,j),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
          
    ! end do
    ! !---------------------------------------------------------------------------------------------- 

    ! !**********************************************************************************************

    !*******************************Stand SPH calculation*********************************
    do k=1,pair_number
        
        i=pair_i(k)                        !particle i index in kth pair
        j=pair_j(k)                        !particle j index in kth pair

        !The normalized derivates of weight function
        kth_pair_dwdx_i(:)=Normalized_dwdx_i(k,:)
        kth_pair_dwdx_j(:)=Normalized_dwdx_j(k,:)
        
        !velocity difference
        velocity_difference(:)=velocity_difference_ij(k,:)+average_velocity(i,:)-average_velocity(j,:)     !Attention: i-j
        
        !velocity difference times weight derivates
        temp_i=DOT_PRODUCT(velocity_difference,kth_pair_dwdx_i)
        temp_j=DOT_PRODUCT(velocity_difference,kth_pair_dwdx_j)

        !-------------------------------------------------------------------------------------------------------------
        !Stand SPH continum equation
        drhodt(i)=drhodt(i)+subdomain_particle_rho(i)*subdomain_particle_volume(j)*temp_i
        drhodt(j)=drhodt(j)+subdomain_particle_rho(j)*subdomain_particle_volume(i)*temp_j
        !-------------------------------------------------------------------------------------------------------------
        
        !************************************************************************************
        !Stand SPH moment equation
        rho_ij=1.0/(subdomain_particle_rho(i)*subdomain_particle_rho(j))

        temp_p_i=0.0d0
        temp_p_j=0.0d0
        temp_e_i=0.0d0
        temp_e_j=0.0d0
            
        !velocity difference
        velocity_difference(:)=-velocity_difference_ij(k,:)    !Attention: j-i
        do L=1,dim
                
           !(pi+pj)/(rho_i*rho_j)*dwdx(i)
           temp_p_i=-(subdomain_particle_press(i)+subdomain_particle_press(j))*Normalized_dwdx_i(k,L)*rho_ij
           temp_p_j=-(subdomain_particle_press(i)+subdomain_particle_press(j))*Normalized_dwdx_j(k,L)*rho_ij

           ! !(pj-pi)/(rho_i*rho_j)*dwdx(i)
           ! temp_p_i=-(-subdomain_particle_press(i)+subdomain_particle_press(j))*Normalized_dwdx_i(k,L)*rho_ij
           ! temp_p_j=-(subdomain_particle_press(i)-subdomain_particle_press(j))*Normalized_dwdx_j(k,L)*rho_ij
           
           temp_e_i=temp_e_i+velocity_difference(l)*temp_p_i
           temp_e_j=temp_e_j+velocity_difference(l)*temp_p_j
                
          !-------------------------------------------------------------------------------------------------------------
          !Stand SPH acceleration
          internal_acceleration(i,l)=internal_acceleration(i,l)+subdomain_particle_mass(j)*temp_p_i
          internal_acceleration(j,l)=internal_acceleration(j,l)-subdomain_particle_mass(i)*temp_p_j
          !-------------------------------------------------------------------------------------------------------------
                
        end do

        dedt(i)=dedt(i)+0.5*subdomain_particle_mass(j)*temp_e_i
        dedt(j)=dedt(j)+0.5*subdomain_particle_mass(i)*temp_e_j
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
            
            do L=1,dim
                temp_i=-capital_pai_ij*Normalized_dwdx_i(k,L)  
                temp_j=-capital_pai_ij*Normalized_dwdx_j(k,L)           
                artifical_viscosity_acceleration(i,L)=artifical_viscosity_acceleration(i,L)+subdomain_particle_mass(j)*temp_i
                artifical_viscosity_acceleration(j,L)=artifical_viscosity_acceleration(j,L)-subdomain_particle_mass(i)*temp_j
            end do    
        end if
        !************************************************************************************

        !************************************************************************************
        
        !------------------------------------------------------------------------------------
        ! !laplacian viscos term
        ! !Reference: Simulation of near-shore solitary wave mechanics by an incompressible SPH method []
        
        ! RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_i)
        ! laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_mass(j)*4*water_kinematic_viscosity*RijDotProductDerivates/(square_distance*(subdomain_particle_rho(i)+subdomain_particle_rho(j)))*v_ij(:)
        
        ! RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_j)
        ! laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_mass(i)*4*water_kinematic_viscosity*RijDotProductDerivates/(square_distance*(subdomain_particle_rho(i)+subdomain_particle_rho(j)))*v_ij(:)
        !------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------
        !It is useful for non-free surface flow
        !Reference:  [Morris]  
        RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_i)
        laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_volume(j)*2*water_kinematic_viscosity*RijDotProductDerivates/square_distance*v_ij(:)
      
        RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_j)
        laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_volume(i)*2*water_kinematic_viscosity*RijDotProductDerivates/square_distance*v_ij(:)
        !------------------------------------------------------------------------------------

        ! !------------------------------------------------------------------------------------
        ! !It is useful for free surface flow 
        ! !Reference: [Theoretical analysis and numerical verification of the consistency of viscous smoothed-particle-hydrodynamics formulations in simulating free-surface flows]
        ! !           [Adaptive particle refinement and derefinement applied to the smoothed particle hydrodynamics methods]

        ! laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_volume(j)*2*(dim+2)*water_kinematic_viscosity*v_ij_multiply_x_ij/square_distance*kth_pair_dwdx_i(:)

        ! laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_volume(i)*2*(dim+2)*water_kinematic_viscosity*v_ij_multiply_x_ij/square_distance*kth_pair_dwdx_j(:)
        ! !------------------------------------------------------------------------------------
        
        !************************************************************************************

        !************************************************************************************
        !artifical viscosity energy
        temp_i=-DOT_PRODUCT(v_ij,kth_pair_dwdx_i)
        temp_j=-DOT_PRODUCT(v_ij,kth_pair_dwdx_j)


        artifical_viscosity_dedt(i)=artifical_viscosity_dedt(i)+0.5*Delta_alpha*average_smooth_length*average_c*water_rho_0/subdomain_particle_rho(i)*v_ij_multiply_x_ij/square_distance*subdomain_particle_volume(j)*temp_i
        artifical_viscosity_dedt(j)=artifical_viscosity_dedt(j)+0.5*Delta_alpha*average_smooth_length*average_c*water_rho_0/subdomain_particle_rho(j)*v_ij_multiply_x_ij/square_distance*subdomain_particle_volume(i)*temp_j

        !********************************************************************************
        !artifical internal acceleration in Delta-SPH(out)
        temp=Delta_alpha*average_smooth_length*average_c*v_ij_multiply_x_ij/square_distance
        do L=1,dim
           artifical_internal_acceleration(i,L)=artifical_internal_acceleration(i,L)+temp*Normalized_dwdx_i(k,L)*water_rho_0/subdomain_particle_rho(i)*subdomain_particle_volume(j)
           artifical_internal_acceleration(j,L)=artifical_internal_acceleration(j,L)-temp*Normalized_dwdx_j(k,L)*water_rho_0/subdomain_particle_rho(j)*subdomain_particle_volume(i)
        end do
        
        !capital_fai_ij

        x_ij=-position_difference_ij(k,:)                              !Attention: j-i
        
        capital_fai_ij=0.0d0
        do l=1,dim
            capital_fai_ij(l)=2*(subdomain_particle_energy(j)-subdomain_particle_energy(i))*x_ij(l)/square_distance
        end do
        

        temp_i=DOT_PRODUCT(capital_fai_ij,kth_pair_dwdx_i)
        temp_j=DOT_PRODUCT(capital_fai_ij,kth_pair_dwdx_j)

        temp_i=chi*temp_i*average_smooth_length*average_c*water_rho_0/average_rho
        temp_j=chi*temp_j*average_smooth_length*average_c*water_rho_0/average_rho
        
        artifical_heat_dedt(i)=artifical_heat_dedt(i)+temp_i*subdomain_particle_volume(j)
        artifical_heat_dedt(j)=artifical_heat_dedt(j)-temp_j*subdomain_particle_volume(i)

        !************************************************************************************
        !SPS viscous term
        temp_p=0.0d0

        !x方向加速度
        do L=1,dim
            
            if(L==1) then
              temp_p=temp_p+(tao_xx(i)+tao_xx(j))*dwdx(k,1) !一维
              !二维，三维情况
              if(dim .ge. 2) then
               
                temp_p=temp_p+(tao_xy(i)+tao_xy(j))*dwdx(k,2)
                
                if(dim==3) then
                    temp_p=temp_p+(tao_xz(i)+tao_xz(j))*dwdx(k,3)
                end if

              end if
                   
            !y方向加速度
            else if(L==2) then
               
               !二维，三维情况
               temp_p=temp_p+(tao_yy(i)+tao_yy(j))*dwdx(k,2)+(tao_xy(i)+tao_xy(j))*dwdx(k,1)
               
               if(dim==3) then
                   temp_p=temp_p+(tao_yz(i)+tao_yz(j))*dwdx(k,3)
               end if
                   
            !z方向速度（说明是三维情况）
            else if(L==3) then
               
               temp_p=temp_p+(tao_zz(i)+tao_zz(j))*dwdx(k,3)+(tao_yz(i)+tao_yz(j))*dwdx(k,2)+(tao_xz(i)+tao_xz(j))*dwdx(k,1)
               
            end if

            !中间值再乘以密度乘积的倒数
            temp_p=temp_p*rho_ij
                
            !存储i，j粒子点加速度,所有粒子对中i,j粒子的加速度的和即为，i,j的最后加速度
            SPS_acceleration(i,L)=SPS_acceleration(i,L)+particle_mass(j)*temp_p
            SPS_acceleration(j,L)=SPS_acceleration(j,L)-particle_mass(i)*temp_p

        enddo
        !**************************************************************************************

        ! !************************************************************************************
        ! !Calculate the force from the boundary particles
        ! if (subdomain_particle_type(j)==-3 .or. subdomain_particle_type(i)==-3) then

        !     position_difference(:)=position_difference_ij(k,:) !×¢ÒâÊÇi-j

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
       psai(:)=2*rho_difference/square_distance*position_difference(:)-(grad_rho(i,:)+grad_rho(j,:))

       !Temp Value
       temp_i=DOT_PRODUCT(psai,kth_pair_dwdx_i)
       temp_j=DOT_PRODUCT(psai,kth_pair_dwdx_j)
       
       artifical_rho_correct(i)=artifical_rho_correct(i)+delta*average_smooth_length*average_c*temp_i*subdomain_particle_volume(j)
       artifical_rho_correct(j)=artifical_rho_correct(j)-delta*average_smooth_length*average_c*temp_j*subdomain_particle_volume(i)
       
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
        drhodt(i)=drhodt(i)+artifical_rho_correct(i)

    end do
    
    
 
 end subroutine Solver_Normalized_Delta_SPH_SPS