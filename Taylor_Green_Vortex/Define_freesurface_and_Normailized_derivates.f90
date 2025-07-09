!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE:: Define_freesurface_and_Normailized_derivates
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

subroutine Define_freesurface_and_Normailized_derivates()             

    use information_module
    use Public_variable_module
    use function_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,L,v,o                                          !Loop Variables
    real(kind=8),dimension(dim)::velocity_difference              !velocity difference
    real(kind=8),dimension(dim)::position_difference              !position difference
    real(kind=8)::support_volume                                  !support doamin volume(PI*R^2)
    real(kind=8)::Normal_vector_Length                            !Normal vector Length
    real(kind=8),dimension(dim)::Virtual_Particle_Position        !Position of virtual particle for free surface searching 

    real(kind=8),dimension(dim,dim)::i_matrix_inversion           !i并矢的逆
    real(kind=8),dimension(dim,dim)::j_matrix_inversion           !i并矢的逆
    real(kind=8),dimension(dim)::right_hand_side
    real(kind=8),dimension(dim,dim)::A,A_inversion
    real(kind=8),dimension(dim,dim)::identity_matrix
    integer::check_error=0

    !Variables for searching near particles 
    real(kind=8)::distance,JI_distance,JT_distance                !i，j粒子点的距离
    real(kind=8)::Normal_distance,Tangential_distance             !i，j粒子点的距离
    real(kind=8),dimension(dim)::I_Tangential_vector              !Tangential vector of ith particle
    real(kind=8),dimension(dim)::I_Normal_vector                  !Normal vector of ith particle

    !Variables for file opeartion
    integer::ioerror=0                                            !open return value
    integer::stat                                                 !read return value
    integer::status                                               !allocation return value

    integer::file_index
    character(len=40)::file_name
    character(len=4)::char_Current_Processor_ID

    !==========================================================================================================

    ! Body of subroutine Define_freesurface_and_Normailized_derivates

    support_volume=PI*(kernel_scale*smooth_length)**dim             !support domain volume

    identity_matrix=0.0d0
    do k=1,dim
        identity_matrix(k,k)=1.0d0
    end do
            
    !**********************************************************************************************************
    !Only calculate the actual particle in the subdomain, donot contain the buffer zone

    subdomain_fra_v=0.0d0                                            !volume fraction
    subdomain_sum_w=0.0d0                                            !sum of the weight function
    
    !calculate the volume
    do k=1,pair_number
        
       i=pair_i(k)                        !particle i index in kth pair
       j=pair_j(k)                        !particle j index in kth pair
         
        subdomain_fra_v(i)=subdomain_fra_v(i)+subdomain_particle_volume(j)/support_volume
        subdomain_fra_v(j)=subdomain_fra_v(j)+subdomain_particle_volume(i)/support_volume

        !Only for the same phase fluid particles
        if (subdomain_particle_type(i)==subdomain_particle_type(j)) then
           subdomain_sum_w(i)=subdomain_sum_w(i)+w(k)*subdomain_particle_volume(j)
           subdomain_sum_w(j)=subdomain_sum_w(j)+w(k)*subdomain_particle_volume(i)
        end if

    end do
    !-----------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------
    !Exchange the sum weight fucntion value from the neighbour processor

    left_Message_ID=4000
    right_Message_ID=4001

    !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
    !Data form left processor to right processor
    call MPI_SENDRECV(subdomain_sum_w(message_first_particle_index_to_right),message_length_send_to_right,MPI_INTEGER,right_processor_ID,left_Message_ID,&
                      subdomain_sum_w(message_first_particle_index_from_left),message_length_from_left,MPI_INTEGER,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(subdomain_sum_w(message_first_particle_index_to_left),message_length_send_to_left,MPI_INTEGER,left_processor_ID,right_Message_ID,&
                      subdomain_sum_w(message_first_particle_index_from_right),message_length_from_right,MPI_INTEGER,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    !-----------------------------------------------------------------------------------------


    !Check all the interior particles (except the buffer zone)
    subdomain_surface_particle_number=0                                 !Initialzed the Variables
    subdomain_surface_particle_name=0
    subdomain_free_surface_type=0

    do i=1,actual_particle_number_in_subdomain

        if(subdomain_particle_type(i)/=2) cycle                         !Only the fluid particle is considered
 
        if(subdomain_fra_v(i)<=0.65) then                                                !这里设置0.7为自由表面粒子，相当于扩宽了自由表面，增加稳定性
           
            subdomain_surface_particle_number=subdomain_surface_particle_number+1
            subdomain_surface_particle_name(subdomain_surface_particle_number)=i         !save the particle index
            subdomain_free_surface_type(i)=1

        end if
               
    end do

    !-----------------------------------------------------------------------------------------


    !*************************************计算CSPH的张量矩阵************************************
    Tensor_dyadic_matrix=0.0d0            !Tensor dyadic matrix
    subdomain_sum_dwdx=0.0d0              !sum dwdx
         
    !Tensor dyadic matrix
    do k=1,pair_number  
         
       i=pair_i(k)                        !particle i index in kth pair
       j=pair_j(k)                        !particle j index in kth pair
       
       ! position difference
       position_difference(:)=-position_difference_ij(k,:)          !Attentation j-i
       
       do L=1,dim
           do o=1,dim
               Tensor_dyadic_matrix(i,L,o)=Tensor_dyadic_matrix(i,L,o)+position_difference(L)*dwdx(k,o)*subdomain_particle_volume(j)
               Tensor_dyadic_matrix(j,L,o)=Tensor_dyadic_matrix(j,L,o)+position_difference(L)*dwdx(k,o)*subdomain_particle_volume(i)
           end do

            subdomain_sum_dwdx(i,L)=subdomain_sum_dwdx(i,L)+dwdx(k,L)*subdomain_particle_volume(j)
            subdomain_sum_dwdx(j,L)=subdomain_sum_dwdx(j,L)-dwdx(k,L)*subdomain_particle_volume(i)
       end do
         
    end do
    !******************************************************************************************
     
    !******************************************************************************************
    !Calculate Tensor dyadic inverse matrix
   
    Tensor_dyadic_inverse_matrix=0.0d0
    !
    !只更新内部(2,k-1)的，边界的要靠邻近粒子传递过来。
    do i=1,actual_particle_number_in_subdomain 

        Tensor_dyadic_inverse_matrix(i,:,:)=identity_matrix                                                                                                     !只对内部粒子点进行更新
        
       if(subdomain_particle_type(i)/=2) cycle          
                
       A=Tensor_dyadic_matrix(i,:,:)                                      !先传到给a的逆阵

       !!调用子程序求解矩阵a的逆(a矩阵会被逆矩阵覆盖)
       ! A_inversion=Tensor_dyadic_matrix(i,:,:)
       ! call BRINV(A_inversion,dim,check_error)

       call matrix_inversion_less_three(dim,A,A_inversion,check_error)
       if(check_error==0) then
           write(*,*) "The tensor inversion compute is wrong! (Normalized) "
       else
           A_inversion=identity_matrix
       end if
       Tensor_dyadic_inverse_matrix(i,:,:)=A_inversion
             
    end do
    !******************************************************************************************
    
    !-----------------------------------------------------------------------------------------
    !Exchange the Tensor_dyadic_inverse_matrix form the neighbour processor
    do j=1,dim
      
       do k=1,dim
        
          left_Message_ID=600+(j-1)*dim+k
          right_Message_ID=700+(j-1)*dim+k

          !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
          !Data form left processor to right processor
          call MPI_SENDRECV(Tensor_dyadic_inverse_matrix(message_first_particle_index_to_right,k,j),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                            Tensor_dyadic_inverse_matrix(message_first_particle_index_from_left,k,j),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
          
          !Data form right processor to left processor
          call MPI_SENDRECV(Tensor_dyadic_inverse_matrix(message_first_particle_index_to_left,k,j),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                            Tensor_dyadic_inverse_matrix(message_first_particle_index_from_right,k,j),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
          
        end do

    end do
    !------------------------------------------------------------------------------------------

    !******************************************************************************************
    !Calculate the normal vector for the free surface particle
    !Out the fluid domain
    New_subdomain_surface_particle_number=subdomain_surface_particle_number

    ! !===========================================================================================
    ! !Method 1:
    ! do L=1,subdomain_surface_particle_number

    !     i=subdomain_surface_particle_name(L)                        !The free surface particle index in subdomain
 
    !     right_hand_side(:)=subdomain_sum_dwdx(i,:)

    !     i_matrix_inversion=Tensor_dyadic_inverse_matrix(i,:,:)

    !     subdomain_Normal_vector(i,:)=matmul(i_matrix_inversion,right_hand_side)

    !     Normal_vector_Length=0.0d0

    !     do k=1,dim
    !         Normal_vector_Length=Normal_vector_Length+subdomain_Normal_vector(i,k)**2
    !     end do
    !     Normal_vector_Length=sqrt(Normal_vector_Length)

    !     subdomain_Normal_vector(i,:)=-subdomain_Normal_vector(i,:)/Normal_vector_Length

    !     !--------------------------------------------------------------------------------------
    !     !Check the free surface particle is right or not.

    !     !Normal and Tangential vector
    !     I_Normal_vector(:)=subdomain_Normal_vector(i,:)
    !     I_Tangential_vector(1)=-subdomain_Normal_vector(i,2)
    !     I_Tangential_vector(2)=subdomain_Normal_vector(i,1)

    !     !Generate a virtual particle in normal vector direction with the length of smooth length
    !     Virtual_Particle_Position(:)=subdomain_particle_position(i,:)+smooth_length*I_Normal_vector(:)

    !     !Search all the near particles of current free surface
    !     do k=1,subdomain_effect_particle_number(i)

    !         j=subdomain_effect_particle(i,k)

    !         if (subdomain_particle_type(j)/=2) cycle          !Only search the fluid particle in vicinity background grid

    !         !calculate the distance and position difference
    !         position_difference(:)=subdomain_particle_position(j,:)-subdomain_particle_position(i,:)
    !         JI_distance=sqrt(DOT_PRODUCT(position_difference,position_difference))

    !         position_difference(:)=subdomain_particle_position(j,:)-Virtual_Particle_Position(:)
    !         JT_distance=sqrt(DOT_PRODUCT(position_difference,position_difference))

    !         Normal_distance=abs(DOT_PRODUCT(position_difference,I_Normal_vector))

    !         Tangential_distance=abs(DOT_PRODUCT(position_difference,I_Tangential_vector))

    !         !Once one criteria is ok, this particle is not free surface particle and we can stop the loop
    !         if (JI_distance>= 1.414*smooth_length .and. JT_distance<smooth_length) then

    !            subdomain_free_surface_type(i)=0
    !            New_subdomain_surface_particle_number=New_subdomain_surface_particle_number-1 

    !            cycle                                           

    !         elseif (JI_distance < 1.414*smooth_length .and. (Normal_distance+Tangential_distance)<smooth_length) then

    !            subdomain_free_surface_type(i)=0
    !            New_subdomain_surface_particle_number=New_subdomain_surface_particle_number-1 

    !            cycle    
              
    !         endif

    !     enddo

    !     !--------------------------------------------------------------------------------------
             
    ! end do
    ! !===========================================================================================

    !===========================================================================================
    !Method 2:
    do L=1,subdomain_surface_particle_number

        i=subdomain_surface_particle_name(L)                        !The free surface particle index in subdomain
 
        right_hand_side(:)=subdomain_sum_dwdx(i,:)

        i_matrix_inversion=Tensor_dyadic_inverse_matrix(i,:,:)

        subdomain_Normal_vector(i,:)=matmul(i_matrix_inversion,right_hand_side)

        Normal_vector_Length=0.0d0

        do k=1,dim
            Normal_vector_Length=Normal_vector_Length+subdomain_Normal_vector(i,k)**2
        end do
        Normal_vector_Length=sqrt(Normal_vector_Length)

        !For the surface particle not in the fluid domain
        if (isnan(Normal_vector_Length)) then
            subdomain_free_surface_type(i)=0
            New_subdomain_surface_particle_number=New_subdomain_surface_particle_number-1 
            cycle
        endif

        !Unit normal vector
        if (Normal_vector_Length/=0.0d0) then
            subdomain_Normal_vector(i,:)=-subdomain_Normal_vector(i,:)/Normal_vector_Length
        else
            subdomain_free_surface_type(i)=0
            New_subdomain_surface_particle_number=New_subdomain_surface_particle_number-1 
            cycle 
        endif
        
        !--------------------------------------------------------------------------------------
        !Check the free surface particle is right or not.

        !Normal and Tangential vector
        I_Normal_vector(:)=subdomain_Normal_vector(i,:)
        I_Tangential_vector(1)=-subdomain_Normal_vector(i,2)
        I_Tangential_vector(2)=subdomain_Normal_vector(i,1)

        !Generate a virtual particle in normal vector direction with the length of smooth length
        Virtual_Particle_Position(:)=subdomain_particle_position(i,:)+smooth_length*I_Normal_vector(:)

        !Search all the near particles of current free surface
        do k=1,subdomain_effect_particle_number(i)

            j=subdomain_effect_particle(i,k)

            if (subdomain_particle_type(j)/=2) cycle          !Only search the fluid particle in vicinity background grid

            !calculate the distance and position difference
            position_difference(:)=subdomain_particle_position(j,:)-subdomain_particle_position(i,:)
            JI_distance=sqrt(DOT_PRODUCT(position_difference,position_difference))

            Normal_distance=DOT_PRODUCT(I_Normal_vector,position_difference)

            position_difference(:)=subdomain_particle_position(j,:)-Virtual_Particle_Position(:)
            JT_distance=sqrt(DOT_PRODUCT(position_difference,position_difference))

            !Once one criteria is ok, this particle is not free surface particle and we can stop the loop
            if (JI_distance>= 1.414*smooth_length .and. JT_distance<smooth_length) then

               subdomain_free_surface_type(i)=0
               New_subdomain_surface_particle_number=New_subdomain_surface_particle_number-1 

               cycle                                           

            elseif (JI_distance < 1.414*smooth_length .and. acos(Normal_distance/JI_distance)<=(PI/4.0)) then

               subdomain_free_surface_type(i)=0
               New_subdomain_surface_particle_number=New_subdomain_surface_particle_number-1 

               cycle    
              
            endif

        enddo

        !--------------------------------------------------------------------------------------
             
    end do
    !===========================================================================================

    !Renew the surface particle number in subdomain
    subdomain_surface_particle_number=New_subdomain_surface_particle_number

    ! !-----------------------------------------------------------------------------------------
    ! !Reduce free surface particle number

    ! !We should Initialzed the variable before reducing
    ! surface_particle_number=0                               
    ! Temp_Reduce_Int_Variable=subdomain_surface_particle_number

    ! call MPI_REDUCE( Temp_Reduce_Int_Variable,surface_particle_number,1,MPI_INTEGER,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    ! call MPI_BCAST( surface_particle_number,1,MPI_INTEGER,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    ! !MPI_BCAST has the Synchronize function so we donot add this here

    ! !-----------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------
    !Exchange the free surface type from the neighbour processor

    left_Message_ID=400
    right_Message_ID=401

    !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
    !Data form left processor to right processor
    call MPI_SENDRECV(subdomain_free_surface_type(message_first_particle_index_to_right),message_length_send_to_right,MPI_INTEGER,right_processor_ID,left_Message_ID,&
                      subdomain_free_surface_type(message_first_particle_index_from_left),message_length_from_left,MPI_INTEGER,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(subdomain_free_surface_type(message_first_particle_index_to_left),message_length_send_to_left,MPI_INTEGER,left_processor_ID,right_Message_ID,&
                      subdomain_free_surface_type(message_first_particle_index_from_right),message_length_from_right,MPI_INTEGER,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    !-----------------------------------------------------------------------------------------

    !*****************************************************************************************

    !*****************************************************************************************
    !Check the particle near free surface or not
    !Initialzed the Variables for near surface checking
    subdomain_near_free_surface_or_not=0
    
    !check all particle pair
    do k=1,pair_number
        
       i=pair_i(k)                        !particle i index in kth pair
       j=pair_j(k)                        !particle j index in kth pair
         
        !if particle i is located at free surface, j should be assgined as the particle near free surface
        if (subdomain_free_surface_type(i)==1) then

            subdomain_near_free_surface_or_not(j)=1

        end if

        !if particle j is located at free surface, i should be assgined as the particle near free surface
        if (subdomain_free_surface_type(j)==1) then

            subdomain_near_free_surface_or_not(i)=1

        end if

    end do
    !***************************************************************************************


end subroutine Define_freesurface_and_Normailized_derivates 