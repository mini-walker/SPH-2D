!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: Devide_block_and_send_tasks
!
!  PURPOSE: Dynamic balance the particle amount in each process:amount_procs 
!           Update: Exchange two column
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

  subroutine Devide_block_and_send_tasks(i_time_step)

    use information_module
    use Public_variable_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the mother subroutien
    integer,intent(in)::i_time_step                                     !current time step 

    ! Variables in program
    integer::i,j,k,L,m                                                  !Variables for loop
    
    !Variables for file opeartion
    integer::ioerror=0                                                  !open return value
    integer::stat                                                       !read return value
    integer::status                                                     !allocation return value
    
    integer::file_index                                                 !file index
    character(len=40)::file_name                                        !file name
    character(len=4)::char_Current_Processor_ID                         !Current Processor ID in character

    !Variables of the MPI Processor
    integer,dimension(Desired_Processor_Number)::buffer_colum_name      !processor interface colum number
        
    integer::sum_to_last_column=0
    integer::sum_to_current_column=0
    integer::start_next_processor=0
    integer::subdomain_index

    integer::chain_name_in_full_view,chain_name_in_block

    integer::particle_number_in_current_chain,Current_Processor_column,current_column_name
    integer::current_full_view_chain_name_in_current_processor
    integer::x_n,y_n,z_n                                                 

    !==========================================================================================================

    ! Body of CSPH_MPI_Parallel_wave_maker_2D

    !**********************************************************************************************************
    ! Get the MPI runing information
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    !**********************************************************************************************************

    !**********************************************************************************************************
    !Check the fluid particles in the domain or not

    do i=1,particle_ture_number
        
        if(particle_type(i)==120) cycle        !No error particles
        
        if(particle_position(i,dim)<0.0d0) then
            particle_type(i)=120
        end if
            
    end do
    !**********************************************************************************************************

    !**********************************************************************************************************
    !Calculate the analytical results
    do i=1,particle_ture_number
                        
        !particle velocity
        Analytical_particle_velocity(i,1)= U_Module*sin(2*PI*particle_position(i,1))*cos(2*PI*particle_position(i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time)   

        Analytical_particle_velocity(i,2)=-U_Module*cos(2*PI*particle_position(i,1))*sin(2*PI*particle_position(i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time) 

        !particle pressure
        Analytical_particle_press(i)=0.25*water_rho_0*U_Module**2*(cos(4*PI*particle_position(i,1))+cos(4*PI*particle_position(i,2)))*exp(-16*PI**2/Reynolds_Number*Nondimensionalized_real_time) 
    
        !velocity magnitude
        Analytical_particle_velocity_magnitude(i)=sqrt(DOT_PRODUCT(Analytical_particle_velocity(i,:),Analytical_particle_velocity(i,:))) 

    end do

    do i=1,total_element_node_number

        !element node velocity
        Analytical_element_node_velocity(i,1)= U_Module*sin(2*PI*element_node_position(i,1))*cos(2*PI*element_node_position(i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time)  

        Analytical_element_node_velocity(i,2)=-U_Module*cos(2*PI*element_node_position(i,1))*sin(2*PI*element_node_position(i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time)

        !element node pressure
        Analytical_element_node_press(i)=0.25*water_rho_0*U_Module**2*(cos(4*PI*element_node_position(i,1))+cos(4*PI*element_node_position(i,2)))*exp(-16*PI**2/Reynolds_Number*Nondimensionalized_real_time) 

        !velocity magnitude
        Analytical_element_node_velocity_magnitude(i)=sqrt(DOT_PRODUCT(Analytical_element_node_velocity(i,:),Analytical_element_node_velocity(i,:))) 

    enddo

    !**********************************************************************************************************

    !**********************************************************************************************************
    !将所有粒子放入网格中
    !注意，在边界粒子赋值，采样子程序和子块网格划分中都有计算网格位置，一个地方改了，其他都要改
    total_particle_number_for_MPI=particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number+Period_Boundary_Particle_Number
    
    number=0                                                         !Initialized the particle number in grid
    do i=1,total_particle_number_for_MPI

        !'0' Type particle means it is a error particle in last step iteration
        !The error particle will be defined as '0' type in the refresh subroutine as type will be reinitialized as 0 at first
        if(particle_type(i)==120 .or. particle_type(i)==0) cycle     !No error particles

        !Calculate the coordinates of the particles in the list
        !注意，这里计算的块的坐标，要减去子块的起始坐标,然后要右移一列
        x_n=floor((particle_position(i,1)-chain_origin_x)/chain_dx+1)
        y_n=floor((particle_position(i,2)-chain_origin_y)/chain_dy+1)

        !Calculate the grid number of the particle in the background grid
        k=(y_n-1)*chain_x_number+x_n                                !Row by row

        !Make sure the particle is located in the Backgroud grid
        if (k<=chain_max_number) then

            !write(*,*) k,chain_max_number
            number(k)=number(k)+1                                   !Particle number in current grid index plus 1
            particle_chain_number(i)=k                              !Save the particle background grid index
            !write(3,*) k
            if(number(k)<=NumberInEachGridPrediction) then
                particle_in_chain(k,number(k))=i                    !The particle index in kth grid is i
            else
                write(*,*) "The chain's containablity is not enough! (Devide Block)"
            end if

        else
            particle_type(i)=120
        endif

        ! !Check error
        ! if (i==1 .and. Current_Processor_ID==0) then
        !   write(*,*) particle_mass(i)
        ! endif

    enddo
    !***********************************************************************************************************

    !Calculate the average particle number in per process
    average_particle_number_per_process=int(total_particle_number_for_MPI/Total_Processors_Number)
    
    ! if (Total_Processors_Number<=6) then
    !     average_particle_number_per_process=int(total_particle_number_for_MPI/Total_Processors_Number)
    ! else
    !     average_particle_number_per_process=int(0.8*total_particle_number_for_MPI/Total_Processors_Number)
    ! endif
    

    !write(*,*) Current_Processor_ID,average_particle_number_per_process

    !*************************************************************************************************************
    !Calculate the particle number in a column
    particle_number_in_each_column_mesh=0

    do i=1,chain_x_number

       particle_number_in_each_column_mesh(i)=0

       do j=1,chain_y_number                                                 !loop in y direction

          k=i+(j-1)*chain_x_number                                           !chain number

          particle_number_in_each_column_mesh(i)=particle_number_in_each_column_mesh(i)+number(k)

      end do 

      ! if (Current_Processor_ID==0) then
      !     write(*,*) i,particle_number_in_each_column_mesh(i)
      ! endif
          
    enddo
    !*************************************************************************************************************

    !*************************************************************************************************************
    !Start the Processors task deviding
    sum_to_last_column=0
    sum_to_current_column=0
    start_next_processor=0
    subdomain_index=0

    !--------------------------------------------------------------------------------------------------------------
    ! Find the start and end column index for each processor.
    !--------------------------------------------------------------------------------------------------------------
    if (Total_Processors_Number/=1) then  ! More than one processor

        do i=1,chain_x_number                                              ! From 1 to max column

          if (i==1 .or. start_next_processor==1) then
             start_next_processor=0
             sum_to_last_column=0
             sum_to_current_column=particle_number_in_each_column_mesh(i)
          else
             sum_to_last_column=sum_to_last_column+particle_number_in_each_column_mesh(i-1)
             sum_to_current_column=sum_to_current_column+particle_number_in_each_column_mesh(i)   
          endif

          if (sum_to_last_column<=average_particle_number_per_process .and. sum_to_current_column>average_particle_number_per_process .and. subdomain_index+1<Total_Processors_Number) then
              
              subdomain_index=subdomain_index+1
              buffer_colum_name(subdomain_index)=i-1
              start_next_processor=1
              sum_to_last_column=0
              sum_to_current_column=0

          endif

          ! IF the Calculator run to the end column, this is the last processor
          if (i==chain_x_number) then

              subdomain_index=Total_Processors_Number
              
          endif

          ! IF the subdomain index meet with the (Processor ID+1), assign the start and end column 
          if (subdomain_index==Current_Processor_ID+1) then

              if (subdomain_index==1) then

                  ! For the first subdomain, we should delete the empty column in the left
                  do k=1,buffer_colum_name(subdomain_index)
                      if (particle_number_in_each_column_mesh(k)/=0) exit  
                  end do

                  start_column(subdomain_index)=k              
                  end_column(subdomain_index)=buffer_colum_name(subdomain_index)

              elseif(subdomain_index==Total_Processors_Number)then

                  ! For the last subdomain, we should delete the empty column in the right
                  start_column(subdomain_index)=buffer_colum_name(subdomain_index-1)+1
                  end_column(subdomain_index)=chain_x_number

                  do k=start_column(subdomain_index),end_column(subdomain_index)
                      if (particle_number_in_each_column_mesh(k)==0) exit  
                  end do

                  end_column(subdomain_index)=k-1

              else

                  ! For the other subdomains, start column name is (end column name of last subdomain) +1
                  ! For the other subdomains, end column name is (buffer_colum_name)
                  start_column(subdomain_index)=buffer_colum_name(subdomain_index-1)+1
                  end_column(subdomain_index)=buffer_colum_name(subdomain_index)

              endif

              exit
              
          endif
          
       enddo
       !--------------------------------------------------------------------------------------------------------------

    else                          
       
       !--------------------------------------------------------------------------------------------------------------
       ! Only one processor
       subdomain_index=1

       do i=1,chain_x_number                                              ! From 1 to max column (Delete the empty column)
          if (particle_number_in_each_column_mesh(i)/=0) exit  
       end do

       start_column(subdomain_index)=i 

       do i=chain_x_number,1,-1                                            ! From max column to 1 (Delete the empty column)
          if (particle_number_in_each_column_mesh(i)/=0) exit  
       end do

       end_column(subdomain_index)=i
       !--------------------------------------------------------------------------------------------------------------

    endif

    !write(*,*) subdomain_index,Current_Processor_ID,start_column(Current_Processor_ID+1),end_column(Current_Processor_ID+1)
    !*************************************************************************************************************


    !*************************************************************************************************************
    !Calculate the chain number for processor
    actual_subdomain_chain_x_number=end_column(Current_Processor_ID+1)-start_column(Current_Processor_ID+1)+1
    total_subdomain_chain_x_number=actual_subdomain_chain_x_number+2*buffer_column_number

    !--------------------------------------------------------------------------------------------------------------
    !Assign the chain number for subdomain
    subdomain_chain_x_number=total_subdomain_chain_x_number
    subdomain_chain_y_number=chain_y_number
    subdomain_chain_z_number=chain_z_number
    subdomain_chain_max_number=subdomain_chain_x_number*subdomain_chain_y_number*subdomain_chain_z_number

    !Assign the origin coordinates of chain number for subdomain
    subdomain_chain_origin_x=chain_origin_x+(start_column(Current_Processor_ID+1)-buffer_column_number)*chain_dx
    subdomain_chain_origin_y=chain_origin_y
    subdomain_chain_origin_z=chain_origin_z
    !--------------------------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------------
    !DeAllocate the memory for subdomain grid
    if( allocated(subdomain_chain_name) ) deallocate( subdomain_chain_name )
    if( allocated(subdomain_column_number) ) deallocate( subdomain_column_number )
    if( allocated(subdomain_each_column_total_particle_number) ) deallocate( subdomain_each_column_total_particle_number )
    
    if( allocated(subdomain_number) ) deallocate( subdomain_number )
    if( allocated(subdomain_particle_in_chain) ) deallocate( subdomain_particle_in_chain )
    if( allocated(subdomain_chain_near_mesh_number) ) deallocate( subdomain_chain_near_mesh_number )
    if( allocated(subdomain_chain_near_mesh_name) ) deallocate( subdomain_chain_near_mesh_name )
    if( allocated(subdomain_mesh_center_positon) ) deallocate( subdomain_mesh_center_positon )

    !--------------------------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------------
    !Allocate the memory for subdomain grid
    Allocate(subdomain_chain_name(subdomain_chain_y_number,subdomain_chain_x_number))
    Allocate(subdomain_column_number(subdomain_chain_x_number))
    Allocate(subdomain_each_column_total_particle_number(subdomain_chain_x_number))
    
    Allocate(subdomain_number(subdomain_chain_max_number))
    Allocate(subdomain_particle_in_chain(subdomain_chain_max_number,NumberInEachGridPrediction))
    Allocate(subdomain_chain_near_mesh_number(subdomain_chain_max_number))
    Allocate(subdomain_chain_near_mesh_name(subdomain_chain_max_number,27))
    Allocate(subdomain_mesh_center_positon(subdomain_chain_max_number,dim))

    !--------------------------------------------------------------------------------------------------------------
    !Distribute Backgroud Grid for subdomain
    !Input: subdomain_chain_x_number,subdomain_chain_y_number,subdomain_chain_z_number
    !Output: subdomain_chain_near_mesh_name,subdomain_chain_near_mesh_number,subdomain_mesh_center_positon
    call Distribute_Backgroud_Grid(dim,&                                      !Dimension of grid
                                   subdomain_chain_x_number,&                 !Grid number in x direction
                                   subdomain_chain_y_number,&                 !Grid number in y direction
                                   subdomain_chain_z_number,&                 !Grid number in z direction
                                   subdomain_chain_max_number,&               !Max grid number (grid_x_number*grid_y_number*grid_z_number)
                                   subdomain_mesh_center_positon,&            !Grid center positon
                                   subdomain_chain_near_mesh_number,&         !Near mesh number
                                   subdomain_chain_near_mesh_name&            !grid near mesh name
                                   )
    
    !--------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------
    !Actual subdomain column
    subdomain_each_column_total_particle_number=0
    do k=1,actual_subdomain_chain_x_number

        i=k+buffer_column_number                                              !The column index in subdomain
        current_column_name=start_column(Current_Processor_ID+1)+k-1
        subdomain_each_column_total_particle_number(i)=particle_number_in_each_column_mesh(current_column_name)

    enddo
    !*************************************************************************************************************

    !*************************************************************************************************************
    !Processor ID for send and receive the exchange data
    ! left_processor_ID
    if (Current_Processor_ID>0) then
        left_processor_ID=Current_Processor_ID-1
    else
        left_processor_ID=MPI_PROC_NULL
    endif

    ! right_processor_ID
    if (Current_Processor_ID<Total_Processors_Number-1) then
        right_processor_ID=Current_Processor_ID+1
    else
        right_processor_ID=MPI_PROC_NULL
    endif
    !*************************************************************************************************************


    !*************************************************************************************************************
    ! Define the column index in subdomain, both method 1 and method 2 are OK!
    !------------------------------------------------------------------------------------------
    ! Method 1: Assign the buffer column directly
    subdomain_chain_name=0
    do k=1,total_subdomain_chain_x_number


        i=start_column(Current_Processor_ID+1)+k-buffer_column_number-1          ! column name

        do j=1,chain_y_number                                                    ! row name

           chain_name_in_full_view=i+(j-1)*chain_x_number
           subdomain_chain_name(j,k)=chain_name_in_full_view

        enddo
       
    enddo
    !------------------------------------------------------------------------------------------
    !write (*,*) i
    !*************************************************************************************************************


    !*************************************************************************************************************
    ! Method 2: Assign the column except buffer zone and exchange column index and chain name of buffer zone
    !------------------------------------------------------------------------------------------
    ! subdomain_chain_name=0
    ! do k=1,total_subdomain_chain_x_number

    !     if (k>=2 .and. k<=total_subdomain_chain_x_number-1) then

    !         i=start_column(Current_Processor_ID+1)+k-2                               ! column name

    !         do j=1,chain_y_number                                                    ! row name

    !            chain_name_in_full_view=i+(j-1)*chain_x_number
    !            subdomain_chain_name(j,k)=chain_name_in_full_view

    !         enddo
            
    !     endif
       
    ! enddo
    ! !------------------------------------------------------------------------------------------

    ! !------------------------------------------------------------------------------------------
    ! ! Define the left and right Processors index

    ! ! send and receive the exchange data
    ! ! left_processor_ID
    ! if (Current_Processor_ID>0) then
    !     left_processor_ID=Current_Processor_ID-1
    ! else
    !     left_processor_ID=MPI_PROC_NULL
    ! endif

    ! ! right_processor_ID
    ! if (Current_Processor_ID<Total_Processors_Number-1) then
    !     right_processor_ID=Current_Processor_ID+1
    ! else
    !     right_processor_ID=MPI_PROC_NULL
    ! endif
    ! !------------------------------------------------------------------------------------------
    
    ! !------------------------------------------------------------------------------------------
    ! ! Exchange column index and chain name of buffer zone 
    ! left_Message_ID=1
    ! right_Message_ID=2

    ! !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
    ! !Data form left processor to right processor
    ! call MPI_SENDRECV(subdomain_chain_name(1,total_subdomain_chain_x_number-1),chain_y_number,MPI_INTEGER,right_processor_ID,left_Message_ID,&
    !                   subdomain_chain_name(1,1),chain_y_number,MPI_INTEGER,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    ! !Data form right processor to left processor
    ! call MPI_SENDRECV(subdomain_chain_name(1,2),chain_y_number,MPI_INTEGER,left_processor_ID,right_Message_ID,&
    !                   subdomain_chain_name(1,total_subdomain_chain_x_number),chain_y_number,MPI_INTEGER,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)


    ! left_Message_ID=3
    ! right_Message_ID=4
    ! !Data form left processor to right processor
    ! call MPI_SENDRECV(subdomain_each_column_total_particle_number(total_subdomain_chain_x_number-1),1,MPI_INTEGER,right_processor_ID,left_Message_ID,&
    !                   subdomain_each_column_total_particle_number(1),1,MPI_INTEGER,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
   
    ! !Data form right processor to left processor
    ! call MPI_SENDRECV(subdomain_each_column_total_particle_number(2),1,MPI_INTEGER,left_processor_ID,right_Message_ID,&
    !                   subdomain_each_column_total_particle_number(total_subdomain_chain_x_number),1,MPI_INTEGER,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    ! !------------------------------------------------------------------------------------------
    ! !*************************************************************************************************************

    !*************************************************************************************************************
    total_particle_number_in_subdomain=0
    actual_particle_number_in_subdomain=0
    particle_number_in_left_buffer_zone=0
    particle_number_in_right_buffer_zone=0

    particle_name_in_left_buffer_zone=0
    particle_name_in_left_buffer_zone=0
    actual_particle_name_in_subdomain=0
    total_particle_name_in_subdomain=0

    do k=1,total_subdomain_chain_x_number

        ! if (Current_Processor_ID==0) then

        !     !write(*,*) subdomain_index,Current_Processor_ID,start_column(Current_Processor_ID+1),end_column(Current_Processor_ID+1)

        !     write(*,*) Current_Processor_ID,k,subdomain_each_column_total_particle_number(k)
            
        ! endif
        
        !write(*,*) Current_Processor_ID,k

        !-------------------------------------------------------------------------------------------------
        !The actual particle in subdomain (2-k-1 column)
        if (k>=buffer_column_number+1 .and. k<=total_subdomain_chain_x_number-buffer_column_number) then

            do j=1,chain_y_number

                current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                if (subdomain_chain_name(j,k)==0) cycle                                        

                do L=1,number(current_full_view_chain_name_in_current_processor)

                    m=particle_in_chain(current_full_view_chain_name_in_current_processor,L) 

                    actual_particle_number_in_subdomain=actual_particle_number_in_subdomain+1

                    actual_particle_name_in_subdomain(actual_particle_number_in_subdomain)=m

                end do

            enddo
        
        endif
        !-------------------------------------------------------------------------------------------------

        !-------------------------------------------------------------------------------------------------
        !The total particle in subdomain (1-k column)
        do j=1,chain_y_number

            current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

            if (subdomain_chain_name(j,k)==0) cycle                                   

            do L=1,number(current_full_view_chain_name_in_current_processor)

                m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)  
                
                total_particle_number_in_subdomain=total_particle_number_in_subdomain+1

                if (total_particle_number_in_subdomain<=subdomain_ntotal) then

                   total_particle_name_in_subdomain(total_particle_number_in_subdomain)=m
                
                else

                  write(*,*) "The meomory for subdomain is not enough, please increase the varibale 'subdomain_ntotal'."
                  
                endif

            end do

        enddo
        !-------------------------------------------------------------------------------------------------


        !-------------------------------------------------------------------------------------------------
        !left buffer zone
        if (k<=buffer_column_number) then

            do j=1,chain_y_number

                current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                if (subdomain_chain_name(j,k)==0) cycle                                   

                do L=1,number(current_full_view_chain_name_in_current_processor)

                    m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)  
                    
                    particle_number_in_left_buffer_zone=particle_number_in_left_buffer_zone+1

                    if (particle_number_in_left_buffer_zone<=column_ntotal) then

                       particle_name_in_left_buffer_zone(particle_number_in_left_buffer_zone)=m

                    else

                       write(*,*) "The meomory for each column is not enough, please increase the varibale 'column_ntotal'."

                    end if

                end do

            enddo
            
        endif
        !-------------------------------------------------------------------------------------------------


        !-------------------------------------------------------------------------------------------------
        !right buffer zone
        if (k>=total_subdomain_chain_x_number-buffer_column_number+1) then

            do j=1,chain_y_number

                current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                if (subdomain_chain_name(j,k)==0) cycle                                  

                do L=1,number(current_full_view_chain_name_in_current_processor)

                    m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)  
                    
                    particle_number_in_right_buffer_zone=particle_number_in_right_buffer_zone+1

                    if (particle_number_in_right_buffer_zone<=column_ntotal) then

                        particle_name_in_right_buffer_zone(particle_number_in_right_buffer_zone)=m

                    else

                        write(*,*) "The meomory for each column is not enough, please increase the varibale 'column_ntotal'."

                    end if


                end do

            enddo
            
        endif
        !-------------------------------------------------------------------------------------------------

    enddo

    !*************************************************************************************************************
    !Calculate the message length and memory index send to left and right
    !Calculate the message length and memory index received from left and right
    message_first_particle_index_to_left=1
    
    message_length_send_to_left=0
    do i=1,buffer_column_number
       message_length_send_to_left=message_length_send_to_left+particle_number_in_each_column_mesh(start_column(Current_Processor_ID+1)+i-1)
    end do
    
    message_length_send_to_right=0
    do i=1,buffer_column_number
       message_length_send_to_right=message_length_send_to_right+particle_number_in_each_column_mesh(end_column(Current_Processor_ID+1)+i-buffer_column_number)
    end do
    message_first_particle_index_to_right=actual_particle_number_in_subdomain-message_length_send_to_right+1


    message_first_particle_index_from_left=actual_particle_number_in_subdomain+1
    message_length_from_left=particle_number_in_left_buffer_zone


    message_first_particle_index_from_right=actual_particle_number_in_subdomain+particle_number_in_left_buffer_zone+1
    message_length_from_right=particle_number_in_right_buffer_zone
    !*************************************************************************************************************


    !write(*,*) Current_Processor_ID,particle_number_in_left_buffer_zone,actual_particle_number_in_subdomain,particle_number_in_right_buffer_zone,total_particle_number_in_subdomain


    !*************************************************************************************************************
    !Re-distribute the points from middle, left, right
    Full_subdomain_particle_name=0

    k=0

    do j=1,actual_particle_number_in_subdomain

       k=k+1

        i=actual_particle_name_in_subdomain(j)

        Full_subdomain_particle_name(k)=i

        subdomain_particle_position(k,:)=particle_position(i,:)
        subdomain_particle_velocity(k,:)=particle_velocity(i,:)
        subdomain_particle_mass(k)=particle_mass(i)                          
        subdomain_particle_type(k)=particle_type(i)                                 
        subdomain_particle_initial_type(k)=particle_initial_type(i)                         
        subdomain_particle_rho(k)=particle_rho(i)                                                 
        subdomain_particle_smooth_lengh(k)=particle_smooth_lengh(i)                    
        subdomain_particle_c(k)=particle_c(i)                               
        subdomain_particle_press(k)=particle_press(i)                           
        subdomain_particle_energy(k)=particle_energy(i)                          
        subdomain_particle_eta(k)=particle_eta(i)                             
        subdomain_free_surface_type(k)=free_surface_type(i)                             
        subdomain_wave_maker_particle_layer(k)=wave_maker_particle_layer(i)
        subdomain_Boundary_particle_type(k)=Boundary_particle_type(i)
        
        subdomain_In_WaveMakerDampingZone_OrNot(k)=In_WaveMakerDampingZone_OrNot(i)
        subdomain_WaveMaker_Analytical_position(k,:)=WaveMaker_Analytical_position(i,:)      
        subdomain_WaveMaker_Analytical_velocity(k,:)=WaveMaker_Analytical_velocity(i,:)     
        subdomain_WaveMaker_Analytical_rho(k)=WaveMaker_Analytical_rho(i)              
        subdomain_WaveMaker_Analytical_press(k)=WaveMaker_Analytical_press(i)  
        subdomain_position_soft_coefficient(k)=position_soft_coefficient(i)           


    end do

    do j=1,particle_number_in_left_buffer_zone

       k=k+1

        i=particle_name_in_left_buffer_zone(j)

        Full_subdomain_particle_name(k)=i

        subdomain_particle_position(k,:)=particle_position(i,:)
        subdomain_particle_velocity(k,:)=particle_velocity(i,:)
        subdomain_particle_mass(k)=particle_mass(i)                          
        subdomain_particle_type(k)=particle_type(i)                                 
        subdomain_particle_initial_type(k)=particle_initial_type(i)                         
        subdomain_particle_rho(k)=particle_rho(i)                                                 
        subdomain_particle_smooth_lengh(k)=particle_smooth_lengh(i)                    
        subdomain_particle_c(k)=particle_c(i)                               
        subdomain_particle_press(k)=particle_press(i)                           
        subdomain_particle_energy(k)=particle_energy(i)                          
        subdomain_particle_eta(k)=particle_eta(i)                             
        subdomain_free_surface_type(k)=free_surface_type(i)                             
        subdomain_wave_maker_particle_layer(k)=wave_maker_particle_layer(i)
        subdomain_Boundary_particle_type(k)=Boundary_particle_type(i)

        subdomain_In_WaveMakerDampingZone_OrNot(k)=In_WaveMakerDampingZone_OrNot(i)
        subdomain_WaveMaker_Analytical_position(k,:)=WaveMaker_Analytical_position(i,:)      
        subdomain_WaveMaker_Analytical_velocity(k,:)=WaveMaker_Analytical_velocity(i,:)     
        subdomain_WaveMaker_Analytical_rho(k)=WaveMaker_Analytical_rho(i)              
        subdomain_WaveMaker_Analytical_press(k)=WaveMaker_Analytical_press(i) 
        subdomain_position_soft_coefficient(k)=position_soft_coefficient(i)

    end do

    do j=1,particle_number_in_right_buffer_zone

       k=k+1

        i=particle_name_in_right_buffer_zone(j)

        Full_subdomain_particle_name(k)=i

        subdomain_particle_position(k,:)=particle_position(i,:)
        subdomain_particle_velocity(k,:)=particle_velocity(i,:)
        subdomain_particle_mass(k)=particle_mass(i)                          
        subdomain_particle_type(k)=particle_type(i)                                 
        subdomain_particle_initial_type(k)=particle_initial_type(i)                         
        subdomain_particle_rho(k)=particle_rho(i)                                                 
        subdomain_particle_smooth_lengh(k)=particle_smooth_lengh(i)                    
        subdomain_particle_c(k)=particle_c(i)                               
        subdomain_particle_press(k)=particle_press(i)                           
        subdomain_particle_energy(k)=particle_energy(i)                          
        subdomain_particle_eta(k)=particle_eta(i)                             
        subdomain_free_surface_type(k)=free_surface_type(i)                             
        subdomain_wave_maker_particle_layer(k)=wave_maker_particle_layer(i)
        subdomain_Boundary_particle_type(k)=Boundary_particle_type(i)
        
        subdomain_In_WaveMakerDampingZone_OrNot(k)=In_WaveMakerDampingZone_OrNot(i)
        subdomain_WaveMaker_Analytical_position(k,:)=WaveMaker_Analytical_position(i,:)      
        subdomain_WaveMaker_Analytical_velocity(k,:)=WaveMaker_Analytical_velocity(i,:)     
        subdomain_WaveMaker_Analytical_rho(k)=WaveMaker_Analytical_rho(i)              
        subdomain_WaveMaker_Analytical_press(k)=WaveMaker_Analytical_press(i)
        subdomain_position_soft_coefficient(k)=position_soft_coefficient(i) 

    end do


!     !*************************************************************************************************************
!     !For Debug
!     if (i_time_step==1) then

!         !------------------------------------------------------------------------------------------
!         !open the output tecplot data files
!         file_index=1+Current_Processor_ID

!         !transfer the Current_Processor_ID from integer to character
!         write(char_Current_Processor_ID,'(I4)') Current_Processor_ID
!         !write(*,*) char_Current_Processor_ID

!         file_name="Total_Partcile_in_subdomain_"//trim(adjustl(char_Current_Processor_ID))//".dat"
!         !write(*,*) file_name

!         open(unit=file_index,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    

!         !写入tecplot抬头
!         write(file_index,*) "TITLE='DISTRIBUTION'"
!         write(file_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "

!         !输出粒子
!         write(file_index,*) "ZONE I=",total_particle_number_in_subdomain," F=POINT"
!         do k=1,total_particle_number_in_subdomain
!             i=total_particle_name_in_subdomain(k)
!             write(file_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
! 100             format(8F20.10) 
!         end do

!         close(file_index)
!         !------------------------------------------------------------------------------------------

!         !------------------------------------------------------------------------------------------
!         !open the output tecplot data files
!         file_index=1+Current_Processor_ID

!         file_name="Actual_Partcile_in_subdomain_"//trim(adjustl(char_Current_Processor_ID))//".dat"
!         !write(*,*) file_name

!         open(unit=file_index,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    

!         !写入tecplot抬头
!         write(file_index,*) "TITLE='DISTRIBUTION'"
!         write(file_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "

!         !输出粒子
!         write(file_index,*) "ZONE I=",actual_particle_number_in_subdomain," F=POINT"
!         do k=1,actual_particle_number_in_subdomain
!             i=actual_particle_name_in_subdomain(k)
!             write(file_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
!         end do

!         close(file_index)
!         !------------------------------------------------------------------------------------------

!         !------------------------------------------------------------------------------------------
!         !open the output tecplot data files
!         file_index=1+Current_Processor_ID

!         file_name="Re_Total_Partcile_in_subdomain_"//trim(adjustl(char_Current_Processor_ID))//".dat"
!         !write(*,*) file_name

!         open(unit=file_index,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    

!         !写入tecplot抬头
!         write(file_index,*) "TITLE='DISTRIBUTION'"
!         write(file_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "

!         !输出粒子
!         write(file_index,*) "ZONE I=",total_particle_number_in_subdomain," F=POINT"
!         do k=1,total_particle_number_in_subdomain
!             i=Full_subdomain_particle_name(k)
!             write(file_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
!         end do

!         close(file_index)
!         !------------------------------------------------------------------------------------------

    
!     endif
!     !*************************************************************************************************************

  end subroutine Devide_block_and_send_tasks