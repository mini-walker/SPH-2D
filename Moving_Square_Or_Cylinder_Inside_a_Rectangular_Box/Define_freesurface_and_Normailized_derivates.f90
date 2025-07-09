!**************************************************************************************************************
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
!**************************************************************************************************************

subroutine Define_freesurface_and_Normailized_derivates()             

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,L,v,o                                          !Loop Variables
    real(kind=8)::support_volume                                  !support doamin volume(PI*R^2)
    real(kind=8)::Normal_vector_Length                            !Normal vector Length
    real(kind=8),dimension(dim)::Virtual_Particle_Position        !Position of virtual particle for free surface searching 

    real(kind=8),dimension(dim,dim)::i_matrix_inversion           !Tensor dyadic matrix of particle i 
    real(kind=8),dimension(dim,dim)::j_matrix_inversion           !Tensor dyadic matrix of particle j
    real(kind=8),dimension(dim)::right_hand_side
    real(kind=8),dimension(dim,dim)::A,A_inversion
    integer::check_error=0

    !Variables for searching near particles 
    real(kind=8)::JI_distance,JT_distance                         !particle distance
    real(kind=8)::Normal_distance,Tangential_distance             !particle distance
    real(kind=8),dimension(dim)::I_Tangential_vector              !Tangential vector of ith particle
    real(kind=8),dimension(dim)::I_Normal_vector                  !Normal vector of ith particle

    !Variables for file
    integer::file_index
    character(len=40)::file_name
    character(len=4)::char_Current_Processor_ID

    !==========================================================================================================

    ! Body of subroutine Define_freesurface_and_Normailized_derivates

    support_volume=PI*(kernel_scale*smooth_length)**dim             !support domain volume

            
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
    !----------------------------------------------------------------------------------------------------------



    !----------------------------------------------------------------------------------------------------------
    !Exchange the sum weight fucntion value from the neighbour processor

    call MPI_RealScalar_MessageExchange(subdomain_sum_w)

    !----------------------------------------------------------------------------------------------------------
    !Check all the interior particles (except the buffer zone)
    subdomain_surface_particle_number=0                                 !Initialzed the Variables
    subdomain_surface_particle_name=0
    subdomain_free_surface_type=0
    
    if (Has_FreeSurface_Or_Not==1 ) then

        do i=1,actual_particle_number_in_subdomain

            if(subdomain_particle_type(i)/=2) cycle                         !Only the fluid particle is considered
     
            if(subdomain_fra_v(i)<=0.65) then                                                !这里设置0.7为自由表面粒子，相当于扩宽了自由表面，增加稳定性
               
                subdomain_surface_particle_number=subdomain_surface_particle_number+1
                subdomain_surface_particle_name(subdomain_surface_particle_number)=i         !save the particle index
                subdomain_free_surface_type(i)=1

            end if
                   
        end do
        
    endif
    !----------------------------------------------------------------------------------------------------------


    !***********************************************计算CSPH的张量矩阵*******************************************
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
    !**********************************************************************************************************


    !**********************************************************************************************************
    !Calculate Tensor dyadic inverse matrix
   
    Tensor_dyadic_inverse_matrix=0.0d0

    !Only the fluid particles and the outlet boundary particles are calculated
    do i=1,actual_particle_number_in_subdomain 

        Tensor_dyadic_inverse_matrix(i,:,:)=identity_matrix                                                                                                     !只对内部粒子点进行更新
        
       if( subdomain_particle_type(i)==Water_Particle_Label .or. subdomain_particle_type(i)==Outlet_Particle_Label ) then          
                
           A=Tensor_dyadic_matrix(i,:,:)                                      !先传到给a的逆阵

           !!调用子程序求解矩阵a的逆(a矩阵会被逆矩阵覆盖)
           ! A_inversion=Tensor_dyadic_matrix(i,:,:)
           ! call BRINV(A_inversion,dim,check_error)

           call matrix_inversion_less_three(dim,A,A_inversion,check_error)
           
           if(check_error==0) then
               write(*,*) "The tensor inversion compute is wrong! (Define Freesurface) "
           else
               A_inversion=identity_matrix
           end if
           Tensor_dyadic_inverse_matrix(i,:,:)=A_inversion

       endif
             
    end do
    !**********************************************************************************************************



    !----------------------------------------------------------------------------------------------------------
    !Exchange the Tensor_dyadic_inverse_matrix form the neighbour processor
    
    call MPI_RealTensor_MessageExchange(Tensor_dyadic_inverse_matrix)

    !----------------------------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------------------------
    !Calculate the normalized derivates
    
    Normalized_dwdx_i=0.0d0
    Normalized_dwdx_j=0.0d0

    do k=1,pair_number
    
       i=pair_i(k)                        !particle i index in kth pair
       j=pair_j(k)                        !particle j index in kth pair

       Normalized_dwdx_i(k,:)=matmul(Tensor_dyadic_inverse_matrix(i,:,:),dwdx(k,:))
       Normalized_dwdx_j(k,:)=matmul(Tensor_dyadic_inverse_matrix(j,:,:),dwdx(k,:))

    end do
    !----------------------------------------------------------------------------------------------------------


    !**********************************************************************************************************
    subdomain_near_free_surface_or_not=0

    if (Has_FreeSurface_Or_Not==1) then
    

        !******************************************************************************************************
        !Calculate the normal vector for the free surface particle
        !Out the fluid domain
        New_subdomain_surface_particle_number=subdomain_surface_particle_number

        ! !=====================================================================================================
        ! !Method 1:
        ! do L=1,subdomain_surface_particle_number

        !     i=subdomain_surface_particle_name(L)                        !The free surface particle index in subdomain
     
        !     subdomain_Normal_vector(i,:)=matmul(Tensor_dyadic_inverse_matrix(i,:,:),subdomain_sum_dwdx(i,:))

        !     Normal_vector_Length=sqrt(DOT_PRODUCT(subdomain_Normal_vector(i,:),subdomain_Normal_vector(i,:)))

        !     subdomain_Normal_vector(i,:)=-subdomain_Normal_vector(i,:)/Normal_vector_Length

        !     !--------------------------------------------------------------------------------------------------
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

        !     !--------------------------------------------------------------------------------------------------
                 
        ! end do
        ! !======================================================================================================


        !======================================================================================================
        !Method 2:
        do L=1,subdomain_surface_particle_number

            i=subdomain_surface_particle_name(L)                        !The free surface particle index in subdomain
     
            subdomain_Normal_vector(i,:)=matmul(Tensor_dyadic_inverse_matrix(i,:,:),subdomain_sum_dwdx(i,:))

            Normal_vector_Length=sqrt(DOT_PRODUCT(subdomain_Normal_vector(i,:),subdomain_Normal_vector(i,:)))

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
            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
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

            !--------------------------------------------------------------------------------------------------
                 
        end do
        !======================================================================================================

        !Renew the surface particle number in subdomain
        subdomain_surface_particle_number=New_subdomain_surface_particle_number

        ! !------------------------------------------------------------------------------------------------------
        ! !Reduce free surface particle number

        ! !We should Initialzed the variable before reducing
        ! surface_particle_number=0                               
        ! Temp_Reduce_Int_Variable=subdomain_surface_particle_number

        ! call MPI_REDUCE( Temp_Reduce_Int_Variable,surface_particle_number,1,MPI_INTEGER,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
        ! call MPI_BCAST( surface_particle_number,1,MPI_INTEGER,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

        ! !MPI_BCAST has the Synchronize function so we donot add this here

        ! !------------------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------------------
        !Exchange the free surface type from the neighbour processor

        call MPI_IntScalar_MessageExchange(subdomain_free_surface_type)

        !------------------------------------------------------------------------------------------------------

        !******************************************************************************************************

        !******************************************************************************************************
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
        !******************************************************************************************************

    endif
    !**********************************************************************************************************


end subroutine Define_freesurface_and_Normailized_derivates 