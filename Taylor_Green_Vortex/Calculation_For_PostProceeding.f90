!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE:: Calculation_For_PostProceeding
!
!  PURPOSE: Calculate the data for Post-Proceeding Grid
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

subroutine Calculation_For_PostProceeding(i_time_step)

    use Public_variable_module
    use function_module
    use information_module
    use MPI

    implicit none

    !==========================================================================================================
    
    !----------------------------------------------------------------------------------------------------------
    !Variables from mother subroutine
    integer,intent(in)::i_time_step                                                  !current time step
    !---------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------
    !Variables for interpolation
    integer::i,j,k,L,v,o,m                                                           !Loop Variables

    !---------------------------------------------------------------------------------------------------------
    !Variables for deviding to each processor
    real(kind=8)::Average_Processor_Length
    
    !---------------------------------------------------------------------------------------------------------
    integer::mesh_name                                                               !Background grid index
    real(kind=8)::temp_w                                                             !Kernel function value
    real(kind=8),dimension(dim)::temp_dwdx                                           !Derivates of kernel function 
    integer,dimension(NumberInEachGridPrediction)::effect_particle                   !Effect particles array
    integer::effect_particle_number                                                  !Effect particle number
    integer::ture_effect_number                                                      !Ture effect number
    integer,dimension(NumberInEachGridPrediction)::ture_effect_name                  !Ture effect particle index
    integer::near_mesh_number                                                        !Near mesh number
    integer::near_mesh_name                                                          !Near mesh name
    integer::particle_in_mesh                                                        !Particle number in mesh
    real(kind=8)::distance                                                           !Distance between i and j
    real(kind=8)::average_smooth_length                                              !Average smooth length
    real(kind=8),dimension(dim)::position_difference                                 !Position difference
    real(kind=8),dimension(dim)::EN_FS_Position_Difference                           !Position difference between element node and free surface particle
    real(kind=8),dimension(dim)::Temp_Velocity                                       !Temp_Velocity
    real(kind=8),dimension(dim)::Nearest_FSP_Normal_Vector                           !Nearest Free Surface Particel Normal Vector
    real(kind=8)::distance_FSP_EN                                                    !Distance between free surface particle and element node
    real(kind=8)::Detected_FSP_EN                                                    !Detected value for free surface particle and element node
    integer::Nearest_Particle_Index                                                  !Nearest Particle Index

    real(kind=8),dimension(NumberInEachGridPrediction)::W_kernel                     !kernel function value
    real(kind=8),dimension(NumberInEachGridPrediction)::Fai                          !weight Value
    !---------------------------------------------------------------------------------------------------------

    !==========================================================================================================

    !**********************************************************************************************************
    ! Get the MPI runing information
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    !**********************************************************************************************************

    !**********************************************************************************************************
    Average_Processor_Length=1.01*PostProceeding_Grid_size_x/Total_Processors_Number

    ! !Devide the element node for each processor
    ! do i=1,total_element_node_number

    !     if ( element_node_position(i,1)>=Current_Processor_ID*Average_Processor_Length .and. element_node_position(i,1)<(Current_Processor_ID+1)*Average_Processor_Length ) then
    !         element_node_processor_index(i)=Current_Processor_ID
    !     endif

    ! enddo
       
    !==========================================================================================================
    ! Assign the element node physical character

    !Get all free surface particles
    surface_particle_number=0
    do i=1,particle_ture_number
        if (free_surface_type(i)==1) then
           surface_particle_number=surface_particle_number+1
           surface_particle_name(surface_particle_number)=i
        end if
    end do
        
    !initialized the Interpolated variables
    Temp_element_node_velocity=0.0d0
    Temp_element_node_press=0.0d0
    Temp_element_node_Fai=0.0d0
    Temp_element_node_vorticity=0.0d0

    interpolation:do i=1,total_element_node_number

        ! !Only for the current processor ID
        ! if ( element_node_processor_index(i)/=Current_Processor_ID ) cycle

        CurrentProcessor:if ( element_node_position(i,1)>=Current_Processor_ID*Average_Processor_Length .and. element_node_position(i,1)<(Current_Processor_ID+1)*Average_Processor_Length ) then
        
            !**************************************************************************************************
            !For the grid not in the dynamic domain
            if ( element_node_position(i,1) < wave_maker_position .and. make_wave_or_not==1 ) then

                temp_element_node_Fai(i)=-1.0d0                            !Not in fluid domain

                cycle                                                      !End current element node
                
            endif
            !**************************************************************************************************    

            !**************************************************************************************************
            !Calculate the leve-set value
            
            if (surface_particle_number>0) then
                
                Distance_FSP_EN=1000000000.0d0                           ! distance between free surface particle and element node
                Nearest_Particle_Index=0                                 ! Nearest Particle Index

                do k=1,surface_particle_number

                   j=surface_particle_name(k)                            ! Surface particle index

                    !----------------------------------------------------------------------------------------------
                    !Get the nearest free surface particle and the Fai value for level-set
                    !Search the nearest the free surface particle
                    !Get the position difference and distance between current node and j_th particle
                    
                    position_difference(:)=element_node_position(i,:)-particle_position(j,:)
                    distance=sqrt(DOT_PRODUCT(position_difference,position_difference))

                    if (distance<=Distance_FSP_EN) then
                        Distance_FSP_EN=distance
                        Nearest_Particle_Index=j
                        EN_FS_Position_Difference= -position_difference        !Element node minus free surface particle
                    endif
                    !----------------------------------------------------------------------------------------------

                enddo

                !Level set value
                !Reference: Fast free-surface detection and level-set function definition in SPH solvers
                Nearest_FSP_Normal_Vector(:)=Normal_vector(Nearest_Particle_Index,:)

                Detected_FSP_EN=DOT_PRODUCT(EN_FS_Position_Difference,Nearest_FSP_Normal_Vector)

                if (Detected_FSP_EN <= -2*smooth_length) then                  

                    temp_element_node_Fai(i)=-1.0d0                            !Not in fluid domain

                elseif (Detected_FSP_EN >= 2*smooth_length) then

                    temp_element_node_Fai(i)=1.0d0                             !In fluid domain
                
                else

                    temp_element_node_Fai(i)=Detected_FSP_EN/(2*smooth_length) !In vicinity of free surface
                    
                endif

            else

                temp_element_node_Fai(i)=1.0d0
                        
            endif

            !**************************************************************************************************

            if (temp_element_node_Fai(i)<=-0.9d0) cycle   !No interpolation for the particle not in fluid domain

            !**************************************************************************************************
            !Interpolation 

            !Get the background grid index and the vicinity grid index
            mesh_name=element_node_chain_number(i)                    !background grid index for current element node
            near_mesh_number=chain_near_mesh_number(mesh_name)        !vicinity grid index for current element node
            
            !**************************************************************************************************
            !Get the effect particle in the near grids
             
            !initialized the effect near index and number
            effect_particle_number=0                                 
            effect_particle=0                                         
             
            do k=1,near_mesh_number                         
                 
                near_mesh_name=chain_near_mesh_name(mesh_name,k)      !current near grid name
                particle_in_mesh=number(near_mesh_name)               !particle number in current grid
               
                do v=1,particle_in_mesh                               !search for each particle
        
                    j=particle_in_chain(near_mesh_name,v)             !Particle index

                    !Only for the fluid particles and not error particle
                    !As the period boundary particles will be deleted after refreshing
                    if(j<=particle_ture_number+wave_maker_particle_number) then
                         
                        effect_particle_number=effect_particle_number+1 
                         
                        !make sure the Memory for the effect_particle_number is enough
                        if(effect_particle_number<=NumberInEachGridPrediction) then
                            effect_particle(effect_particle_number)=j
                        else
                            write(*,'(A)') " Memory for effect_particle is not engough, it should be larger than 300! (Post-Proceding)"
                        endif

                    endif
                   
                enddo
            enddo
            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            W_kernel=0.0d0
            ture_effect_number=0                                     ! ture effect number for i_th node
            ture_effect_name=0                                       ! ture effect name for i_th node
            
            do L=1,effect_particle_number
            
                !fluid particle index
                j=effect_particle(L)
                
                !Get the position difference and distance between current node and j_th particle

                position_difference(:)=element_node_position(i,:)-particle_position(j,:)
                distance=sqrt(DOT_PRODUCT(position_difference,position_difference))

                !Average smooth length
                average_smooth_length=smooth_length
                !write(*,*) particle_smooth_lengh(i)
                
                !Get the particles in support domain (except for the distance==0.0)
                if(distance<=(kernel_scale*average_smooth_length) .and. distance>=1.0e-7) then
                    
                    ture_effect_number=ture_effect_number+1              !ture effect particle number
                    ture_effect_name(ture_effect_number)=j               !ture effect particle index
                    
                    !Get the kernel function
                    call compute_kernel(dim,&                            !dimension
                                        distance,&                       !distance(in)                
                                        position_difference,&            !position difference(in)
                                        average_smooth_length,&          !average smooth length(in)
                                        temp_w,&                         !kernel function value(out)
                                        temp_dwdx&                       !derivates of kernel function(out)
                                        )
                    
                    !kernel function value
                    W_kernel(ture_effect_number)=temp_w
      
                endif      
                
            enddo
            
            !--------------------------------------------------------------------------------------------------

            !Interpolate fixed particle quality
            if(ture_effect_number>=3) then
                     
                !Get the weight value of the effect particle
                call MLS_Interpolation_In_Gobal_Domain(MLS_order_Grid,element_node_position(i,:),ture_effect_number,ture_effect_name,W_kernel,Fai)
           
                do L=1,ture_effect_number
            
                    !Get the fluid particle index
                    j=ture_effect_name(L)

                    !interpolation
                    if (trim(adjustl(Gravity_Switch))=='on') then
                        Temp_element_node_press(i)=Temp_element_node_press(i)+particle_press(j)*Fai(L)-g*particle_rho(j)*(element_node_position(i,dim)-particle_position(j,dim))*Fai(L)
                    else
                        Temp_element_node_press(i)=Temp_element_node_press(i)+particle_press(j)*Fai(L)
                    endif

                    Temp_element_node_velocity(i,:)=Temp_element_node_velocity(i,:)+particle_velocity(j,:)*Fai(L)

                    Temp_element_node_vorticity(i,:)=Temp_element_node_vorticity(i,:)+particle_vorticity(j,:)*Fai(L)

                enddo
                
            endif
            
            !**************************************************************************************************

        endif CurrentProcessor

    enddo interpolation
    !==========================================================================================================

    !==========================================================================================================
    !Reduce all the data in processor 

    Reduce_Real_Array_For_Grid=Temp_element_node_press

    call MPI_REDUCE( Reduce_Real_Array_For_Grid,element_node_press,total_element_node_number,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( element_node_press,total_element_node_number,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    Reduce_Real_Array_For_Grid=temp_element_node_Fai

    call MPI_REDUCE( Reduce_Real_Array_For_Grid,element_node_Fai,total_element_node_number,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( element_node_Fai,total_element_node_number,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    do j=1,dim

        Reduce_Real_Array_For_Grid=Temp_element_node_velocity(:,j)

        call MPI_REDUCE( Reduce_Real_Array_For_Grid,element_node_velocity(:,j),total_element_node_number,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
        call MPI_BCAST( element_node_velocity(:,j),total_element_node_number,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    enddo

    do j=1,3

        Reduce_Real_Array_For_Grid=temp_element_node_vorticity(:,j)

        call MPI_REDUCE( Reduce_Real_Array_For_Grid,element_node_vorticity(:,j),total_element_node_number,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
        call MPI_BCAST( element_node_vorticity(:,j),total_element_node_number,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    enddo
    !==========================================================================================================

end subroutine Calculation_For_PostProceeding