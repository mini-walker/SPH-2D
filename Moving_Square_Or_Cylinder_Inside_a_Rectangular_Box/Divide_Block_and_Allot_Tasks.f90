!**************************************************************************************************************
!  SUBROUTINE : Divide_Block_and_Allot_Tasks
!
!  PURPOSE    : Dynamic balance the particle amount in each process:amount_procs 
!               Update: Exchange two column
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

subroutine Divide_Block_and_Allot_Tasks(i_time_step)

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the superior subroutien
    integer,intent(in)::i_time_step                                     ! Current time step 

    ! Variables in program
    integer::i,j,k,L,m                                                  ! Variables for loop
    
    ! Variables for file opeartion
    integer::File_index                                                 ! File index
    character(len=100)::File_name                                       ! File name
    character(len=4)::Char_Current_Processor_ID                         ! Current Processor ID in character

    ! Variables of the MPI Processor
    integer,dimension(Desired_Processor_Number)::buffer_colum_name      ! Processor interface colum number
        
    integer::sum_to_last_column=0
    integer::sum_to_current_column=0
    integer::start_next_processor=0
    integer::subdomain_index

    integer::chain_name_in_full_view,chain_name_in_block

    integer::particle_number_in_current_chain,Current_Processor_column,current_column_name
    integer::current_full_view_chain_name_in_current_processor
    integer::x_n,y_n,z_n 
    integer,dimension(dim)::Grid_Position_Index                        ! The position index of grid which point in                                               

    !==========================================================================================================


    ! Body of Divide_Block_and_Allot_Tasks


    ! !**********************************************************************************************************
    ! !Calculate the analytical results
    ! do i=1,Particle_ture_number
                        
    !     !particle velocity
    !     Analytical_particle_velocity(i,1)= U_Module*sin(2*PI*particle_position(i,1))*cos(2*PI*particle_position(i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time)   

    !     Analytical_particle_velocity(i,2)=-U_Module*cos(2*PI*particle_position(i,1))*sin(2*PI*particle_position(i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time) 

    !     !particle pressure
    !     Analytical_particle_press(i)=0.25*water_rho_0*U_Module**2*(cos(4*PI*particle_position(i,1))+cos(4*PI*particle_position(i,2)))*exp(-16*PI**2/Reynolds_Number*Nondimensionalized_real_time) 
    
    !     !velocity magnitude
    !     Analytical_particle_velocity_magnitude(i)=sqrt(DOT_PRODUCT(Analytical_particle_velocity(i,:),Analytical_particle_velocity(i,:))) 

    ! end do

    ! do i=1,total_element_node_number

    !     !element node velocity
    !     Analytical_element_node_velocity(i,1)= U_Module*sin(2*PI*element_node_position(i,1))*cos(2*PI*element_node_position(i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time)  

    !     Analytical_element_node_velocity(i,2)=-U_Module*cos(2*PI*element_node_position(i,1))*sin(2*PI*element_node_position(i,2))*exp(-8*PI**2/Reynolds_Number*Nondimensionalized_real_time)

    !     !element node pressure
    !     Analytical_element_node_press(i)=0.25*water_rho_0*U_Module**2*(cos(4*PI*element_node_position(i,1))+cos(4*PI*element_node_position(i,2)))*exp(-16*PI**2/Reynolds_Number*Nondimensionalized_real_time) 

    !     !velocity magnitude
    !     Analytical_element_node_velocity_magnitude(i)=sqrt(DOT_PRODUCT(Analytical_element_node_velocity(i,:),Analytical_element_node_velocity(i,:))) 

    ! enddo
    ! !**********************************************************************************************************

    !**********************************************************************************************************

    ! All particles need project into background grid
    Total_particle_number_for_MPI=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number+Period_Boundary_Particle_Number+Traditional_ghost_particle_number

    number=0                                                                        ! Initialized the particle number in grid
    
    do i=1,Total_particle_number_for_MPI

        ! '0' Type particle means it is a error particle in last step iteration
        ! The error particle will be defined as '0' type in the refresh subroutine as type will be reinitialized as 0 at first
        if(particle_type(i)==Ill_Particle_Label .or. particle_type(i)==0) cycle     ! No error particles

        ! Input value: particle_position
        ! Return value: grid position index, grid name
        call Point_Project_Into_Gobal_Background_Grid(particle_position(i,:),Grid_Position_Index,k)

        x_n = Grid_Position_Index(1);
        y_n = Grid_Position_Index(2);

        ! Make sure the particle is located in the Backgroud grid
        if ( x_n>=1 .and. y_n>=1 .and. x_n<=chain_x_number .and. y_n<=chain_y_number) then

            !-------------------------------------------------------------------------------------------
            !write(*,*) k,chain_max_number

            particle_chain_number(i)=k                                              ! Save the particle background grid index

            if( number(k) < PredictionNumberInSingleGrid ) then
                number(k)=number(k)+1                                               ! Particle number in current grid index plus 1
                particle_in_chain(k,number(k))=i                                    ! The particle index in kth grid is i
            else
                write(*,*) "The chain's containablity is not enough! (Divide Block)"
            endif
            !-------------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------------
            ! Assign the maximum and minimum value of the chain index in x and y 
            if (i==1) then

                MaximunChainIndex_X=x_n
                MinimunChainIndex_X=x_n
                MaximunChainIndex_Y=y_n
                MinimunChainIndex_Y=y_n

            else

                ! maximum chain index in x
                if ( x_n > MaximunChainIndex_X ) then
                    MaximunChainIndex_X=x_n
                endif

                ! minimum chain index in x
                if ( x_n < MinimunChainIndex_X ) then
                    MinimunChainIndex_X=x_n
                endif

                ! maximum chain index in y
                if ( y_n > MaximunChainIndex_Y ) then
                    MaximunChainIndex_Y=y_n
                endif

                ! minimum chain index in y
                if ( y_n < MinimunChainIndex_Y ) then
                    MinimunChainIndex_Y=y_n
                endif
              
            endif
            !-------------------------------------------------------------------------------------------

        else
            particle_type(i) = Ill_Particle_Label
        endif


    enddo
    !***********************************************************************************************************


    !***********************************************************************************************************
    ! Devide the blocks
    DevideBlock_if: if (MPI_Dim==1) then

        !=======================================================================================================

        !*******************************************************************************************************
        !Calculate the average particle number in per process
        average_particle_number_per_process=int(Total_particle_number_for_MPI/Total_Processors_Number)
        
        ! if (Total_Processors_Number<=6) then
        !     average_particle_number_per_process=int(Total_particle_number_for_MPI/Total_Processors_Number)
        ! else
        !     average_particle_number_per_process=int(0.8*Total_particle_number_for_MPI/Total_Processors_Number)
        ! endif
        

        !write(*,*) Current_Processor_ID,average_particle_number_per_process

        !*********************************************************************************************************
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
        !*********************************************************************************************************

        !*********************************************************************************************************
        !Start the Processors task deviding
        sum_to_last_column=0
        sum_to_current_column=0
        start_next_processor=0
        subdomain_index=0

        !----------------------------------------------------------------------------------------------------------
        ! Find the start and end column index for each processor.
        !----------------------------------------------------------------------------------------------------------
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
           !-----------------------------------------------------------------------------------------------------

        else                          
           
           !-----------------------------------------------------------------------------------------------------
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
           !-----------------------------------------------------------------------------------------------------

        endif

        !write(*,*) subdomain_index,Current_Processor_ID,start_column(Current_Processor_ID+1),end_column(Current_Processor_ID+1)
        !*********************************************************************************************************


        !*********************************************************************************************************
        !Calculate the chain number for processor
        actual_subdomain_chain_x_number=end_column(Current_Processor_ID+1)-start_column(Current_Processor_ID+1)+1
        actual_subdomain_chain_y_number=chain_y_number

        total_subdomain_chain_x_number=actual_subdomain_chain_x_number+2*buffer_column_number
        total_subdomain_chain_y_number=chain_y_number

        !----------------------------------------------------------------------------------------------------------
        !Assign the chain number for subdomain
        subdomain_chain_x_number=total_subdomain_chain_x_number
        subdomain_chain_y_number=total_subdomain_chain_y_number
        subdomain_chain_z_number=chain_z_number
        subdomain_chain_max_number=subdomain_chain_x_number*subdomain_chain_y_number*subdomain_chain_z_number

        !Assign the origin coordinates of chain number for subdomain
        subdomain_chain_origin_x=chain_origin_x+(start_column(Current_Processor_ID+1)-buffer_column_number)*chain_dx
        subdomain_chain_origin_y=chain_origin_y
        subdomain_chain_origin_z=chain_origin_z
        !----------------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------------
        !DeAllocate the memory for subdomain grid
        if( allocated(subdomain_chain_name) ) deallocate( subdomain_chain_name )

        if( allocated(subdomain_number) ) deallocate( subdomain_number )
        if( allocated(subdomain_particle_in_chain) ) deallocate( subdomain_particle_in_chain )
        if( allocated(subdomain_chain_near_mesh_number) ) deallocate( subdomain_chain_near_mesh_number )
        if( allocated(subdomain_chain_near_mesh_name) ) deallocate( subdomain_chain_near_mesh_name )
        if( allocated(subdomain_mesh_center_positon) ) deallocate( subdomain_mesh_center_positon )

        !----------------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------------
        !Allocate the memory for subdomain grid
        Allocate(subdomain_chain_name(subdomain_chain_y_number,subdomain_chain_x_number))
        
        Allocate(subdomain_number(subdomain_chain_max_number))
        Allocate(subdomain_particle_in_chain(subdomain_chain_max_number,PredictionNumberInSingleGrid))
        Allocate(subdomain_chain_near_mesh_number(subdomain_chain_max_number))
        Allocate(subdomain_chain_near_mesh_name(subdomain_chain_max_number,27))
        Allocate(subdomain_mesh_center_positon(subdomain_chain_max_number,dim))

        !----------------------------------------------------------------------------------------------------------
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
        
        !----------------------------------------------------------------------------------------------------------


        !***********************************************************************************************************



        !***********************************************************************************************************
        ! Processor ID for send and receive the exchange data
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
        !***********************************************************************************************************


        !***********************************************************************************************************
        ! Define the column index in subdomain, both method 1 and method 2 are OK!
        !-----------------------------------------------------------------------------------------------------------
        ! Method 1: Assign the buffer column directly
        subdomain_chain_name=0
        do k=1,total_subdomain_chain_x_number

            i=start_column(Current_Processor_ID+1)+k-buffer_column_number-1          ! column name

            do j=1,chain_y_number                                                    ! row name

               chain_name_in_full_view=i+(j-1)*chain_x_number
               subdomain_chain_name(j,k)=chain_name_in_full_view

            enddo
           
        enddo
        !-----------------------------------------------------------------------------------------------------------

        !***********************************************************************************************************


        !***********************************************************************************************************
        ! Calculate the particle number and record the particle index in origin full domain
        ! Devide the actual particle domain into four corners, four boundaries and center part
        ! Example: 
        !         |   |         |   |
        !         | L |         | R |
        !         | S |  Center | S |
        !         |   |         |   |
        !       LS: left side; RS: right side
        
        particle_number_in_LS_Of_ActualSubdomain=0 ; particle_name_in_LS_Of_ActualSubdomain=0
        particle_number_in_RS_Of_ActualSubdomain=0 ; particle_name_in_RS_Of_ActualSubdomain=0

        particle_number_in_Center_Of_ActualSubdomain=0 ; particle_name_in_Center_Of_ActualSubdomain=0
        
        particle_number_in_LS_buffer_zone=0 ;  particle_name_in_LS_buffer_zone=0
        particle_number_in_RS_buffer_zone=0 ;  particle_name_in_RS_buffer_zone=0


        do j=1,total_subdomain_chain_y_number
            do k=1,total_subdomain_chain_x_number

              !----------------------------------------------------------------------------------------------------
              ! The actual particle in subdomain (2---k-1 column)
              if (k>=buffer_column_number+1 .and. k<=total_subdomain_chain_x_number-buffer_column_number ) then
                
                  !================================================================================================
                  ! Left part in actual subdomain
                  if ( k<=buffer_column_number+buffer_column_number ) then

                      !--------------------------------------------------------------------------------------------

                      !--------------------------------------------------------------------------------------------
                      ! LS: left side 

                      ! Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_LS_Of_ActualSubdomain=particle_number_in_LS_Of_ActualSubdomain+1
                          particle_name_in_LS_Of_ActualSubdomain(particle_number_in_LS_Of_ActualSubdomain)=m

                      enddo
                      !-------------------------------------------------------------------------------------------- 

                  !================================================================================================
                  ! Right part in actual subdomain
                  elseif ( k>=total_subdomain_chain_x_number-2*buffer_column_number+1 ) then

                      !--------------------------------------------------------------------------------------------
                      ! RS: right side 

                      ! Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          Particle_number_in_RS_Of_ActualSubdomain=Particle_number_in_RS_Of_ActualSubdomain+1
                          particle_name_in_RS_Of_ActualSubdomain(Particle_number_in_RS_Of_ActualSubdomain)=m

                      enddo

                      !-------------------------------------------------------------------------------------------- 

                  !================================================================================================
                  ! Center part
                  else 

                      !--------------------------------------------------------------------------------------------
                      ! Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          Particle_number_in_Center_Of_ActualSubdomain=Particle_number_in_Center_Of_ActualSubdomain+1
                          particle_name_in_Center_Of_ActualSubdomain(Particle_number_in_Center_Of_ActualSubdomain)=m

                      enddo
                      !-------------------------------------------------------------------------------------------- 

                  endif 
                  !================================================================================================ 

              !----------------------------------------------------------------------------------------------------
              ! Receive data from left
              elseif(k<buffer_column_number+1) then

                  ! Get the full view chain name 
                  current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                  if (subdomain_chain_name(j,k)==0) cycle

                  do L=1,number(current_full_view_chain_name_in_current_processor)

                      m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                      particle_number_in_LS_buffer_zone=particle_number_in_LS_buffer_zone+1
                      particle_name_in_LS_buffer_zone(particle_number_in_LS_buffer_zone)=m

                  enddo

              !----------------------------------------------------------------------------------------------------
              ! Receive data from right
              elseif(k>total_subdomain_chain_x_number-buffer_column_number) then

                  ! Get the full view chain name 
                  current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                  if (subdomain_chain_name(j,k)==0) cycle

                  do L=1,number(current_full_view_chain_name_in_current_processor)

                      m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                      particle_number_in_RS_buffer_zone=particle_number_in_RS_buffer_zone+1
                      particle_name_in_RS_buffer_zone(particle_number_in_RS_buffer_zone)=m

                  enddo

              endif 

            enddo
        enddo
        !***********************************************************************************************************


        !***********************************************************************************************************
        ! Re-order the particle index in subdomain

        Actual_particle_number_in_subdomain=0 ;  actual_particle_name_in_subdomain=0
        Total_particle_number_in_subdomain=0  ;  total_particle_name_in_subdomain=0
        
        Full_subdomain_particle_name=0

        !-----------------------------------------------------------------------------------------------------
        ! Two sides in actual subdomain (anti-clockwise)
        ! LS---RS

        do i=1,Particle_number_in_LS_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_LS_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_LS_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_LS_Of_ActualSubdomain(i)

        enddo

        do i=1,Particle_number_in_RS_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_RS_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_RS_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_RS_Of_ActualSubdomain(i)

        enddo

        !-----------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------
        ! Center part in subdomain
        do i=1,Particle_number_in_Center_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_Center_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_Center_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_Center_Of_ActualSubdomain(i)

        enddo
        !-----------------------------------------------------------------------------------------------------


        !-----------------------------------------------------------------------------------------------------
        ! Two sides in the buffer zone, star from the FS (anti-clockwise)
        do i=1,particle_number_in_LS_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_LS_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_LS_buffer_zone(i)

        enddo

        do i=1,particle_number_in_RS_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_RS_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_RS_buffer_zone(i)

        enddo
        !-----------------------------------------------------------------------------------------------------


        !*************************************************************************************************************


        !*************************************************************************************************************
        ! Calculate the message length and memory index send to left and right
        ! Calculate the message length and memory index received from left and right
        message_first_particle_index_to_left=1
        message_length_send_to_left=Particle_number_in_LS_Of_ActualSubdomain

        message_first_particle_index_to_right=Particle_number_in_LS_Of_ActualSubdomain+1
        message_length_send_to_right=Particle_number_in_RS_Of_ActualSubdomain


        message_first_particle_index_from_left=Actual_particle_number_in_subdomain+1
        message_length_from_left=particle_number_in_LS_buffer_zone


        message_first_particle_index_from_right=Actual_particle_number_in_subdomain+particle_number_in_LS_buffer_zone+1
        message_length_from_right=particle_number_in_RS_buffer_zone
        !*************************************************************************************************************


    !===============================================================================================================
    ! MPI 2D division
    elseif (MPI_Dim==2) then

        !Define the chain index for each processor
        ChainLength_X=MaximunChainIndex_X-MinimunChainIndex_X+1
        ChainLength_Y=MaximunChainIndex_Y-MinimunChainIndex_Y+1
        ChainLength_Z=MaximunChainIndex_Z-MinimunChainIndex_Z+1

        ! !-----------------------------------------------------------------------------------------------------------
        ! ! Method 1: Devide by the average length
        ! AverageChainLength_X=int( DFLOAT(ChainLength_X)/DFLOAT(Processor_Number_X) )    ! Attention: this is integer devided by integer 
        ! AverageChainLength_Y=int( DFLOAT(ChainLength_Y)/DFLOAT(Processor_Number_Y) )
        ! AverageChainLength_Z=int( DFLOAT(ChainLength_Z)/DFLOAT(Processor_Number_Z) )

        ! ! Loop each processor
        ! do i=1,Processor_Number_Y                                     ! Row by row

        !     do j=1,Processor_Number_X                                 ! column by column

        !         CurrentProcessorIndex=(i-1)*Processor_Number_X+j-1

        !         !Current Processor Position
        !         CurrentProcessorPosition(1)=j-1
        !         CurrentProcessorPosition(2)=i-1
        !         CurrentProcessorPosition(3)=0

        !         ! Calculate the current processor information special for the boundary 

        !         if ( CurrentProcessorIndex==Current_Processor_ID ) then

        !             !-----------------------------------------------------------------------------------------------
        !             ! Chain index in x direction
                    
        !             ActualCurrentProcessorOriginChain(1)=MinimunChainIndex_X+CurrentProcessorPosition(1)*AverageChainLength_X

        !             if (j/=Processor_Number_X) then

        !                 ActualCurrentProcessorEndChain(1)=MinimunChainIndex_X+(CurrentProcessorPosition(1)+1)*AverageChainLength_X-1

        !             else
                        
        !                 ActualCurrentProcessorEndChain(1)=MaximunChainIndex_X

        !             endif

        !             TotalCurrentProcessorOriginChain(1)=ActualCurrentProcessorOriginChain(1)-buffer_column_number
        !             TotalCurrentProcessorEndChain(1)=ActualCurrentProcessorEndChain(1)+buffer_column_number
        !             !-----------------------------------------------------------------------------------------------

        !             !-----------------------------------------------------------------------------------------------
        !             ! Chain index in y direction

        !             ActualCurrentProcessorOriginChain(2)=MinimunChainIndex_Y+CurrentProcessorPosition(2)*AverageChainLength_Y

        !             if (i/=Processor_Number_Y) then

        !                 ActualCurrentProcessorEndChain(2)=MinimunChainIndex_Y+(CurrentProcessorPosition(2)+1)*AverageChainLength_Y-1

        !             else

        !                 ActualCurrentProcessorEndChain(2)=MaximunChainIndex_Y
                      
        !             endif

        !             TotalCurrentProcessorOriginChain(2)=ActualCurrentProcessorOriginChain(2)-buffer_column_number
        !             TotalCurrentProcessorEndChain(2)=ActualCurrentProcessorEndChain(2)+buffer_column_number
        !             !-----------------------------------------------------------------------------------------------

        !             !-----------------------------------------------------------------------------------------------
        !             ! Actial and Total (Actial+2*buffer_column_number) Chain length for each processor

        !             ActualCurrentProcessorChainLength(1)=ActualCurrentProcessorEndChain(1)-ActualCurrentProcessorOriginChain(1)+1
        !             ActualCurrentProcessorChainLength(2)=ActualCurrentProcessorEndChain(2)-ActualCurrentProcessorOriginChain(2)+1

        !             TotalCurrentProcessorChainLength(1)=ActualCurrentProcessorChainLength(1)+2*buffer_column_number
        !             TotalCurrentProcessorChainLength(2)=ActualCurrentProcessorChainLength(2)+2*buffer_column_number

        !             !-----------------------------------------------------------------------------------------------

        !             ! ! Check the results
        !             ! if ( Current_Processor_ID==0 ) then
        !             !     write(*,*) MinimunChainIndex_X,MaximunChainIndex_X,MinimunChainIndex_Y,MaximunChainIndex_Y
        !             ! endif
                    
        !             ! write(*,*) Current_Processor_ID,ActualCurrentProcessorOriginChain(1),ActualCurrentProcessorEndChain(1),ActualCurrentProcessorOriginChain(2),ActualCurrentProcessorEndChain(2)
        !             !-----------------------------------------------------------------------------------------------

        !             !-----------------------------------------------------------------------------------------------
        !             !Processor ID for send and receive the exchange data

        !             left_processor_ID=Current_Processor_ID-1
        !             right_processor_ID=Current_Processor_ID+1
        !             Front_processor_ID=Current_Processor_ID-Processor_Number_X
        !             Back_processor_ID=Current_Processor_ID+Processor_Number_X
                    
        !             LeftFront_processor_ID=Current_Processor_ID-1-Processor_Number_X
        !             RightBack_processor_ID=Current_Processor_ID+1+Processor_Number_X
        !             LeftBack_processor_ID=Current_Processor_ID-1+Processor_Number_X
        !             RightFront_processor_ID=Current_Processor_ID+1-Processor_Number_X


        !             !special for the sides and corners
        !             if (i==1) then
        !                 Front_processor_ID=MPI_PROC_NULL
        !                 LeftFront_processor_ID=MPI_PROC_NULL
        !                 RightFront_processor_ID=MPI_PROC_NULL
        !             endif

        !             if (i==Processor_Number_Y) then
        !                 Back_processor_ID=MPI_PROC_NULL
        !                 LeftBack_processor_ID=MPI_PROC_NULL
        !                 RightBack_processor_ID=MPI_PROC_NULL
        !             endif

        !             if (j==1) then
        !                 left_processor_ID=MPI_PROC_NULL
        !                 LeftFront_processor_ID=MPI_PROC_NULL
        !                 LeftBack_processor_ID=MPI_PROC_NULL
        !             endif

        !             if (j==Processor_Number_X) then
        !                 right_processor_ID=MPI_PROC_NULL
        !                 RightFront_processor_ID=MPI_PROC_NULL
        !                 RightBack_processor_ID=MPI_PROC_NULL
        !             endif

        !             !-----------------------------------------------------------------------------------------------

        !         endif

        !     enddo

        ! enddo
        ! !-----------------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------------
        ! Method 2: Devide by the average at each length
       
        ! Calculate the actual length of chain number for every processor

        EachProcessorChainLength_X=0
        SumColumnChainLength_X=0

        do j=1,Processor_Number_X                                     ! column by column

            if ( j==1 ) then                                                                                                     !For first column
                EachProcessorChainLength_X(j)=NINT( DFLOAT(ChainLength_X)/DFLOAT(Processor_Number_X) )
            elseif ( j==Processor_Number_X ) then                                                                                !For last column
                EachProcessorChainLength_X(j)=ChainLength_X-SumColumnChainLength_X
            else                                                                                                                 !Others
                EachProcessorChainLength_X(j)=NINT( DFLOAT(ChainLength_X-SumColumnChainLength_X)/DFLOAT(Processor_Number_X-j+1) )
            endif

            SumColumnChainLength_X=SumColumnChainLength_X+EachProcessorChainLength_X(j)

            ! Check the deviding
            if (EachProcessorChainLength_X(j)<=0) then
               write(*,*) 'Errors in X Deviding!'
            endif

        enddo
        
        EachProcessorChainLength_Y=0
        SumColumnChainLength_Y=0

        do i=1,Processor_Number_Y                                     ! Row by row

            if ( i==1 ) then                                                                                                     !For first column
                EachProcessorChainLength_Y(i)=NINT( DFLOAT(ChainLength_Y)/DFLOAT(Processor_Number_Y) )
            elseif ( i==Processor_Number_Y ) then                                                                                !For last column
                EachProcessorChainLength_Y(i)=ChainLength_Y-SumColumnChainLength_Y
            else                                                                                                                 !Others
                EachProcessorChainLength_Y(i)=NINT( DFLOAT(ChainLength_Y-SumColumnChainLength_Y)/DFLOAT(Processor_Number_Y-i+1) )
            endif

            SumColumnChainLength_Y=SumColumnChainLength_Y+EachProcessorChainLength_Y(i)

            ! Check the deviding
            if (EachProcessorChainLength_Y(i)<=0) then
               write(*,*) 'Errors in Y Deviding!'
            endif

        enddo


        ! Loop each processor

        SumColumnChainLength_Y=0

        do i=1,Processor_Number_Y                                                             ! Row by row
            
            SumColumnChainLength_X=0
            
            do j=1,Processor_Number_X                                                         ! column by column

                CurrentProcessorIndex=(i-1)*Processor_Number_X+j-1

                ! Current Processor Position
                CurrentProcessorPosition(1)=j-1
                CurrentProcessorPosition(2)=i-1
                CurrentProcessorPosition(3)=0

                ! Calculate the current processor information special for the boundary 
                if ( CurrentProcessorIndex==Current_Processor_ID ) then

                    !-----------------------------------------------------------------------------------------------
                    ! Chain index in x direction
                    
                    ActualCurrentProcessorOriginChain(1)=MinimunChainIndex_X+SumColumnChainLength_X
                    ActualCurrentProcessorEndChain(1)=MinimunChainIndex_X+SumColumnChainLength_X+EachProcessorChainLength_X(j)-1

                    TotalCurrentProcessorOriginChain(1)=ActualCurrentProcessorOriginChain(1)-buffer_column_number
                    TotalCurrentProcessorEndChain(1)=ActualCurrentProcessorEndChain(1)+buffer_column_number
                    !-----------------------------------------------------------------------------------------------

                    !-----------------------------------------------------------------------------------------------
                    ! Chain index in y direction

                    ActualCurrentProcessorOriginChain(2)=MinimunChainIndex_Y+SumColumnChainLength_Y
                    ActualCurrentProcessorEndChain(2)=MinimunChainIndex_Y+SumColumnChainLength_Y+EachProcessorChainLength_Y(i)-1

                    TotalCurrentProcessorOriginChain(2)=ActualCurrentProcessorOriginChain(2)-buffer_column_number
                    TotalCurrentProcessorEndChain(2)=ActualCurrentProcessorEndChain(2)+buffer_column_number
                    !-----------------------------------------------------------------------------------------------

                    !-----------------------------------------------------------------------------------------------
                    ! Actial and Total (Actial+2*buffer_column_number) Chain length for each processor

                    ActualCurrentProcessorChainLength(1)=ActualCurrentProcessorEndChain(1)-ActualCurrentProcessorOriginChain(1)+1
                    ActualCurrentProcessorChainLength(2)=ActualCurrentProcessorEndChain(2)-ActualCurrentProcessorOriginChain(2)+1

                    TotalCurrentProcessorChainLength(1)=ActualCurrentProcessorChainLength(1)+2*buffer_column_number
                    TotalCurrentProcessorChainLength(2)=ActualCurrentProcessorChainLength(2)+2*buffer_column_number

                    !-----------------------------------------------------------------------------------------------

                    ! ! Check the results
                    ! if ( Current_Processor_ID==0 ) then
                    !     write(*,*) MinimunChainIndex_X,MaximunChainIndex_X,MinimunChainIndex_Y,MaximunChainIndex_Y
                    ! endif
                    
                    ! write(*,*) Current_Processor_ID,ActualCurrentProcessorOriginChain(1),ActualCurrentProcessorEndChain(1),ActualCurrentProcessorOriginChain(2),ActualCurrentProcessorEndChain(2)
                    !-----------------------------------------------------------------------------------------------

                    !-----------------------------------------------------------------------------------------------
                    ! Processor ID for send and receive the exchange data

                    left_processor_ID=Current_Processor_ID-1
                    right_processor_ID=Current_Processor_ID+1
                    Front_processor_ID=Current_Processor_ID-Processor_Number_X
                    Back_processor_ID=Current_Processor_ID+Processor_Number_X
                    
                    LeftFront_processor_ID=Current_Processor_ID-1-Processor_Number_X
                    RightBack_processor_ID=Current_Processor_ID+1+Processor_Number_X
                    LeftBack_processor_ID=Current_Processor_ID-1+Processor_Number_X
                    RightFront_processor_ID=Current_Processor_ID+1-Processor_Number_X


                    ! Special for the sides and corners
                    if (i==1) then
                        Front_processor_ID=MPI_PROC_NULL
                        LeftFront_processor_ID=MPI_PROC_NULL
                        RightFront_processor_ID=MPI_PROC_NULL
                    endif

                    if (i==Processor_Number_Y) then
                        Back_processor_ID=MPI_PROC_NULL
                        LeftBack_processor_ID=MPI_PROC_NULL
                        RightBack_processor_ID=MPI_PROC_NULL
                    endif

                    if (j==1) then
                        left_processor_ID=MPI_PROC_NULL
                        LeftFront_processor_ID=MPI_PROC_NULL
                        LeftBack_processor_ID=MPI_PROC_NULL
                    endif

                    if (j==Processor_Number_X) then
                        right_processor_ID=MPI_PROC_NULL
                        RightFront_processor_ID=MPI_PROC_NULL
                        RightBack_processor_ID=MPI_PROC_NULL
                    endif

                    !-----------------------------------------------------------------------------------------------

                endif
                
                SumColumnChainLength_X=SumColumnChainLength_X+EachProcessorChainLength_X(j)   !The chain in x direction

            enddo

            SumColumnChainLength_Y=SumColumnChainLength_Y+EachProcessorChainLength_Y(i)       !The chain in y direction

        enddo
        !-----------------------------------------------------------------------------------------------------------

        !***********************************************************************************************************
        ! Calculate the chain number for processor
        actual_subdomain_chain_x_number=ActualCurrentProcessorChainLength(1)
        actual_subdomain_chain_y_number=ActualCurrentProcessorChainLength(2)

        ! if (actual_subdomain_chain_x_number<=0.0d0 .or. actual_subdomain_chain_y_number<=0.0d0) then

        !    write(*,*) Current_Processor_ID
          
        ! endif
        
        total_subdomain_chain_x_number=TotalCurrentProcessorChainLength(1)
        total_subdomain_chain_y_number=TotalCurrentProcessorChainLength(2)
        
        !-----------------------------------------------------------------------------------------------------------
        ! Assign the chain number for subdomain
        subdomain_chain_x_number=total_subdomain_chain_x_number
        subdomain_chain_y_number=total_subdomain_chain_y_number
        subdomain_chain_z_number=chain_z_number
        subdomain_chain_max_number=subdomain_chain_x_number*subdomain_chain_y_number*subdomain_chain_z_number

        ! Assign the origin coordinates of chain number for subdomain
        subdomain_chain_origin_x=chain_origin_x+TotalCurrentProcessorOriginChain(1)*chain_dx
        subdomain_chain_origin_y=chain_origin_y+TotalCurrentProcessorOriginChain(2)*chain_dy
        subdomain_chain_origin_z=chain_origin_z

        !write(*,*) ActualCurrentProcessorOriginChain(1),ActualCurrentProcessorOriginChain(2),TotalCurrentProcessorOriginChain(1),TotalCurrentProcessorOriginChain(2),Current_Processor_ID
        !-----------------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------------
        ! Deallocate the memory for subdomain grid
        if( allocated(subdomain_chain_name) ) deallocate( subdomain_chain_name )

        
        if( allocated(subdomain_number) ) deallocate( subdomain_number )
        if( allocated(subdomain_particle_in_chain) ) deallocate( subdomain_particle_in_chain )
        if( allocated(subdomain_chain_near_mesh_number) ) deallocate( subdomain_chain_near_mesh_number )
        if( allocated(subdomain_chain_near_mesh_name) ) deallocate( subdomain_chain_near_mesh_name )
        if( allocated(subdomain_mesh_center_positon) ) deallocate( subdomain_mesh_center_positon )

        !-----------------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------------
        ! Allocate the memory for subdomain grid
        Allocate(subdomain_chain_name(subdomain_chain_y_number,subdomain_chain_x_number))
        
        Allocate(subdomain_number(subdomain_chain_max_number))
        Allocate(subdomain_particle_in_chain(subdomain_chain_max_number,PredictionNumberInSingleGrid))
        Allocate(subdomain_chain_near_mesh_number(subdomain_chain_max_number))
        Allocate(subdomain_chain_near_mesh_name(subdomain_chain_max_number,27))
        Allocate(subdomain_mesh_center_positon(subdomain_chain_max_number,dim))

        !-----------------------------------------------------------------------------------------------------------
        ! Distribute Backgroud Grid for subdomain
        ! Input: subdomain_chain_x_number,subdomain_chain_y_number,subdomain_chain_z_number
        ! Output: subdomain_chain_near_mesh_name,subdomain_chain_near_mesh_number,subdomain_mesh_center_positon
        call Distribute_Backgroud_Grid(dim,&                                      ! Dimension of grid
                                       subdomain_chain_x_number,&                 ! Grid number in x direction
                                       subdomain_chain_y_number,&                 ! Grid number in y direction
                                       subdomain_chain_z_number,&                 ! Grid number in z direction
                                       subdomain_chain_max_number,&               ! Max grid number (grid_x_number*grid_y_number*grid_z_number)
                                       subdomain_mesh_center_positon,&            ! Grid center positon
                                       subdomain_chain_near_mesh_number,&         ! Near mesh number
                                       subdomain_chain_near_mesh_name&            ! Grid near mesh name
                                       )
        
        !-----------------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------------
        ! Actual subdomain column
        ! Assign the buffer column directly
        subdomain_chain_name=0
        k=0  
        do i=TotalCurrentProcessorOriginChain(2),TotalCurrentProcessorEndChain(2)
            
            k=k+1                  ! Row index
            L=0                    ! Initialized column index

            do j=TotalCurrentProcessorOriginChain(1),TotalCurrentProcessorEndChain(1)

                L=L+1              ! Column index

               chain_name_in_full_view=j+(i-1)*chain_x_number
               subdomain_chain_name(K,L)=chain_name_in_full_view

            enddo
        enddo
        !-----------------------------------------------------------------------------------------------------------

        !***********************************************************************************************************
        ! Calculate the particle number and record the particle index in origin full domain
        ! Devide the actual particle domain into four corners, four boundaries and center part
        ! Example: |LBC|    BS   |RBC|
        !          |   |         |   |
        !          | L |         | R |
        !          | S |  Center | S |
        !          |   |         |   |
        !          |LFC|    Fs   |RFC|
        !        LBC: left back corner; RBC: right back corner; LFC: left front corner; LFC: left front corner;
        !        LS: left side; RS: right side
        

        Particle_number_in_LFC_Of_ActualSubdomain=0 ; particle_name_in_LFC_Of_ActualSubdomain=0
        particle_number_in_RFC_Of_ActualSubdomain=0 ; particle_name_in_RFC_Of_ActualSubdomain=0
        particle_number_in_LBC_Of_ActualSubdomain=0 ; particle_name_in_LBC_Of_ActualSubdomain=0
        particle_number_in_RBC_Of_ActualSubdomain=0 ; particle_name_in_RBC_Of_ActualSubdomain=0
        

        particle_number_in_LS_Of_ActualSubdomain=0 ; particle_name_in_LS_Of_ActualSubdomain=0
        particle_number_in_RS_Of_ActualSubdomain=0 ; particle_name_in_RS_Of_ActualSubdomain=0
        particle_number_in_FS_Of_ActualSubdomain=0 ; particle_name_in_FS_Of_ActualSubdomain=0
        particle_number_in_BS_Of_ActualSubdomain=0 ; particle_name_in_BS_Of_ActualSubdomain=0

        particle_number_in_Center_Of_ActualSubdomain=0 ; particle_name_in_Center_Of_ActualSubdomain=0
        

        do j=1,total_subdomain_chain_y_number
            do k=1,total_subdomain_chain_x_number

              !----------------------------------------------------------------------------------------------------
              ! The actual particle in subdomain (2---k-1 column)
              if (k>=buffer_column_number+1 .and. k<=total_subdomain_chain_x_number-buffer_column_number .and. j>=buffer_column_number+1 .and. j<=total_subdomain_chain_y_number-buffer_column_number) then
                
                  !================================================================================================
                  ! Left part in actual subdomain
                  if ( k<=buffer_column_number+buffer_column_number ) then

                      !--------------------------------------------------------------------------------------------
                      ! LFC: left front corner;
                      if ( j<=buffer_column_number+buffer_column_number ) then
                
                          ! Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              Particle_number_in_LFC_Of_ActualSubdomain=Particle_number_in_LFC_Of_ActualSubdomain+1
                              particle_name_in_LFC_Of_ActualSubdomain(Particle_number_in_LFC_Of_ActualSubdomain)=m

                          enddo
                      !--------------------------------------------------------------------------------------------

                      !-------------------------------------------------------------------------------------------- 
                      ! LBC: left back corner;
                      elseif ( j>=total_subdomain_chain_y_number-2*buffer_column_number+1 ) then
                         
                          ! Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              Particle_number_in_LBC_Of_ActualSubdomain=Particle_number_in_LBC_Of_ActualSubdomain+1
                              particle_name_in_LBC_Of_ActualSubdomain(Particle_number_in_LBC_Of_ActualSubdomain)=m

                          enddo
                      !--------------------------------------------------------------------------------------------
                      ! LS: left side 
                      else

                          ! Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              Particle_number_in_LS_Of_ActualSubdomain=Particle_number_in_LS_Of_ActualSubdomain+1
                              particle_name_in_LS_Of_ActualSubdomain(Particle_number_in_LS_Of_ActualSubdomain)=m

                          enddo
                      !-------------------------------------------------------------------------------------------- 
                      endif

                  !================================================================================================
                  ! Right part in actual subdomain
                  elseif ( k>=total_subdomain_chain_x_number-2*buffer_column_number+1 ) then

                      !--------------------------------------------------------------------------------------------
                      ! RFC: right front corner;
                      if ( j<=buffer_column_number+buffer_column_number ) then
                
                          ! Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              Particle_number_in_RFC_Of_ActualSubdomain=Particle_number_in_RFC_Of_ActualSubdomain+1
                              particle_name_in_RFC_Of_ActualSubdomain(Particle_number_in_RFC_Of_ActualSubdomain)=m

                          enddo
                      !--------------------------------------------------------------------------------------------

                      !-------------------------------------------------------------------------------------------- 
                      ! RBC: right back corner;
                      elseif ( j>=total_subdomain_chain_y_number-2*buffer_column_number+1 ) then
                         
                          ! Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              Particle_number_in_RBC_Of_ActualSubdomain=Particle_number_in_RBC_Of_ActualSubdomain+1
                              particle_name_in_RBC_Of_ActualSubdomain(Particle_number_in_RBC_Of_ActualSubdomain)=m

                          enddo
                      !--------------------------------------------------------------------------------------------
                      ! RS: right side 
                      else

                          ! Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              Particle_number_in_RS_Of_ActualSubdomain=Particle_number_in_RS_Of_ActualSubdomain+1
                              particle_name_in_RS_Of_ActualSubdomain(Particle_number_in_RS_Of_ActualSubdomain)=m

                          enddo
                      !-------------------------------------------------------------------------------------------- 
                      endif

                  !================================================================================================
                  ! Center part
                  else 

                      !--------------------------------------------------------------------------------------------
                      ! FS: front side;
                      if ( j<=buffer_column_number+buffer_column_number ) then
                
                          ! Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              Particle_number_in_FS_Of_ActualSubdomain=Particle_number_in_FS_Of_ActualSubdomain+1
                              particle_name_in_FS_Of_ActualSubdomain(Particle_number_in_FS_Of_ActualSubdomain)=m

                          enddo
                      !--------------------------------------------------------------------------------------------

                      !-------------------------------------------------------------------------------------------- 
                      ! BS: back side;
                      elseif ( j>=total_subdomain_chain_y_number-2*buffer_column_number+1 ) then
                         
                          ! Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              Particle_number_in_BS_Of_ActualSubdomain=Particle_number_in_BS_Of_ActualSubdomain+1
                              particle_name_in_BS_Of_ActualSubdomain(Particle_number_in_BS_Of_ActualSubdomain)=m

                          enddo
                      !--------------------------------------------------------------------------------------------
                      ! RS: right side 
                      else

                          ! Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              Particle_number_in_Center_Of_ActualSubdomain=Particle_number_in_Center_Of_ActualSubdomain+1
                              particle_name_in_Center_Of_ActualSubdomain(Particle_number_in_Center_Of_ActualSubdomain)=m

                          enddo
                      !-------------------------------------------------------------------------------------------- 
                      endif

                  endif 
                  !================================================================================================     

              endif 

            enddo
        enddo
        !***********************************************************************************************************

        !***********************************************************************************************************
        ! Define the message received from the near Processors
        
        particle_number_in_LFC_buffer_zone=0 ;  particle_name_in_LFC_buffer_zone=0
        particle_number_in_Back_LFC_buffer_zone=0 ;  particle_name_in_Back_LFC_buffer_zone=0
        particle_number_in_Right_LFC_buffer_zone=0 ;  particle_name_in_Right_LFC_buffer_zone=0
        

        particle_number_in_RFC_buffer_zone=0 ;  particle_name_in_RFC_buffer_zone=0
        particle_number_in_Back_RFC_buffer_zone=0 ;  particle_name_in_Back_RFC_buffer_zone=0
        particle_number_in_Left_RFC_buffer_zone=0 ;  particle_name_in_Left_RFC_buffer_zone=0
        

        particle_number_in_LBC_buffer_zone=0 ;  particle_name_in_LBC_buffer_zone=0
        particle_number_in_Front_LBC_buffer_zone=0 ;  particle_name_in_Front_LBC_buffer_zone=0
        particle_number_in_Right_LBC_buffer_zone=0 ;  particle_name_in_Right_LBC_buffer_zone=0
        

        particle_number_in_RBC_buffer_zone=0 ;  particle_name_in_RBC_buffer_zone=0
        particle_number_in_Front_RBC_buffer_zone=0 ;  particle_name_in_Front_RBC_buffer_zone=0
        particle_number_in_Left_RBC_buffer_zone=0 ;  particle_name_in_Left_RBC_buffer_zone=0


        particle_number_in_LS_buffer_zone=0 ;  particle_name_in_LS_buffer_zone=0
        particle_number_in_RS_buffer_zone=0 ;  particle_name_in_RS_buffer_zone=0
        particle_number_in_FS_buffer_zone=0 ;  particle_name_in_FS_buffer_zone=0
        particle_number_in_BS_buffer_zone=0 ; particle_name_in_BS_buffer_zone=0
        

        do j=1,total_subdomain_chain_y_number
            do k=1,total_subdomain_chain_x_number

              !Donot contain the actual parts only the buffer domain
              if (k>=buffer_column_number+1 .and. k<=total_subdomain_chain_x_number-buffer_column_number .and. j>=buffer_column_number+1 .and. j<=total_subdomain_chain_y_number-buffer_column_number) cycle
            
              !================================================================================================
              !Left part
              if (k<buffer_column_number+1) then

                  !--------------------------------------------------------------------------------------------
                  !Left Front Corner in buffer zone
                  if (j<buffer_column_number+1) then

                      !Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_LFC_buffer_zone=particle_number_in_LFC_buffer_zone+1
                          particle_name_in_LFC_buffer_zone(particle_number_in_LFC_buffer_zone)=m

                      enddo

                  !--------------------------------------------------------------------------------------------
                  !Back Left Front Corner in buffer zone
                  elseif (j>=buffer_column_number+1 .and. j<=2*buffer_column_number) then

                      !Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_Back_LFC_buffer_zone=particle_number_in_Back_LFC_buffer_zone+1
                          particle_name_in_Back_LFC_buffer_zone(particle_number_in_Back_LFC_buffer_zone)=m

                      enddo

                  !--------------------------------------------------------------------------------------------
                  !Left Back Corner in buffer zone
                  elseif (j>total_subdomain_chain_y_number-buffer_column_number) then

                      !Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_LBC_buffer_zone=particle_number_in_LBC_buffer_zone+1
                          particle_name_in_LBC_buffer_zone(particle_number_in_LBC_buffer_zone)=m

                      enddo

                  !--------------------------------------------------------------------------------------------
                  !Front Left Back Corner in buffer zone
                  elseif (j>total_subdomain_chain_y_number-2*buffer_column_number .and. j<=total_subdomain_chain_y_number-buffer_column_number) then

                      !Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_Front_LBC_buffer_zone=particle_number_in_Front_LBC_buffer_zone+1
                          particle_name_in_Front_LBC_buffer_zone(particle_number_in_Front_LBC_buffer_zone)=m

                      enddo

                  !--------------------------------------------------------------------------------------------
                  !LS in buffer zone
                  else

                      !Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_LS_buffer_zone=particle_number_in_LS_buffer_zone+1
                          particle_name_in_LS_buffer_zone(particle_number_in_LS_buffer_zone)=m

                      enddo

                  !--------------------------------------------------------------------------------------------
                  endif


              !================================================================================================
              !Right part
              elseif (k>total_subdomain_chain_x_number-buffer_column_number) then

                  !--------------------------------------------------------------------------------------------
                  !Right Front Corner in buffer zone
                  if (j<buffer_column_number+1) then

                      !Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_RFC_buffer_zone=particle_number_in_RFC_buffer_zone+1
                          particle_name_in_RFC_buffer_zone(particle_number_in_RFC_buffer_zone)=m

                      enddo

                  !--------------------------------------------------------------------------------------------
                  !Back Right Front Corner in buffer zone
                  elseif (j>=buffer_column_number+1 .and. j<=2*buffer_column_number) then

                      !Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_Back_RFC_buffer_zone=particle_number_in_Back_RFC_buffer_zone+1
                          particle_name_in_Back_RFC_buffer_zone(particle_number_in_Back_RFC_buffer_zone)=m

                      enddo

                  !--------------------------------------------------------------------------------------------
                  !Right Back Corner in buffer zone
                  elseif (j>total_subdomain_chain_y_number-buffer_column_number) then

                      !Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_RBC_buffer_zone=particle_number_in_RBC_buffer_zone+1
                          particle_name_in_RBC_buffer_zone(particle_number_in_RBC_buffer_zone)=m

                      enddo

                  !--------------------------------------------------------------------------------------------
                  !Front Right Back Corner in buffer zone
                  elseif (j>total_subdomain_chain_y_number-2*buffer_column_number .and. j<=total_subdomain_chain_y_number-buffer_column_number) then

                      !Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_Front_RBC_buffer_zone=particle_number_in_Front_RBC_buffer_zone+1
                          particle_name_in_Front_RBC_buffer_zone(particle_number_in_Front_RBC_buffer_zone)=m

                      enddo

                  !--------------------------------------------------------------------------------------------
                  !RS in buffer zone
                  else

                      !Get the full view chain name 
                      current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                      if (subdomain_chain_name(j,k)==0) cycle

                      do L=1,number(current_full_view_chain_name_in_current_processor)

                          m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                          particle_number_in_RS_buffer_zone=particle_number_in_RS_buffer_zone+1
                          particle_name_in_RS_buffer_zone(particle_number_in_RS_buffer_zone)=m

                      enddo

                  !--------------------------------------------------------------------------------------------
                  endif

              !================================================================================================
              !Center part
              else

                  !Right of the left corners
                  if (k<=2*buffer_column_number .and. k>buffer_column_number) then

                      !----------------------------------------------------------------------------------------
                      !Right Left Front Corner
                      if (j<=buffer_column_number) then
                          
                          !Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              particle_number_in_Right_LFC_buffer_zone=particle_number_in_Right_LFC_buffer_zone+1
                              particle_name_in_Right_LFC_buffer_zone(particle_number_in_Right_LFC_buffer_zone)=m

                          enddo

                      !----------------------------------------------------------------------------------------
                      !Right Left Back Corner
                      else

                          !Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              particle_number_in_Right_LBC_buffer_zone=particle_number_in_Right_LBC_buffer_zone+1
                              particle_name_in_Right_LBC_buffer_zone(particle_number_in_Right_LBC_buffer_zone)=m

                          enddo
                        
                      endif
                      !----------------------------------------------------------------------------------------

                  !Left of the right corners
                  elseif (k>total_subdomain_chain_x_number-2*buffer_column_number .and. k<=total_subdomain_chain_x_number-buffer_column_number) then

                      !----------------------------------------------------------------------------------------
                      !Left Right Front Corner
                      if (j<=buffer_column_number) then
                          
                          !Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              particle_number_in_Left_RFC_buffer_zone=particle_number_in_Left_RFC_buffer_zone+1
                              particle_name_in_Left_RFC_buffer_zone(particle_number_in_Left_RFC_buffer_zone)=m

                          enddo

                      !----------------------------------------------------------------------------------------
                      !Left Right Back Corner
                      else

                          !Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              particle_number_in_Left_RBC_buffer_zone=particle_number_in_Left_RBC_buffer_zone+1
                              particle_name_in_Left_RBC_buffer_zone(particle_number_in_Left_RBC_buffer_zone)=m

                          enddo
                        
                      endif
                      !----------------------------------------------------------------------------------------

                  !Center part of the FS and BS
                  else

                      !----------------------------------------------------------------------------------------
                      !FS
                      if (j<=buffer_column_number) then
                          
                          !Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              particle_number_in_FS_buffer_zone=particle_number_in_FS_buffer_zone+1
                              particle_name_in_FS_buffer_zone(particle_number_in_FS_buffer_zone)=m

                          enddo

                      !----------------------------------------------------------------------------------------
                      !BS
                      else

                          !Get the full view chain name 
                          current_full_view_chain_name_in_current_processor=subdomain_chain_name(j,k)

                          if (subdomain_chain_name(j,k)==0) cycle

                          do L=1,number(current_full_view_chain_name_in_current_processor)

                              m=particle_in_chain(current_full_view_chain_name_in_current_processor,L)
                              particle_number_in_BS_buffer_zone=particle_number_in_BS_buffer_zone+1
                              particle_name_in_BS_buffer_zone(particle_number_in_BS_buffer_zone)=m

                          enddo
                        
                      endif
                      !----------------------------------------------------------------------------------------
                    
                  endif

              endif
              !================================================================================================

            enddo
        enddo
        !******************************************************************************************************

    
        !******************************************************************************************************
        !Re-order the particle index in the actual domain

        Actual_particle_number_in_subdomain=0
        actual_particle_name_in_subdomain=0

        Total_particle_name_in_subdomain=0
        Total_particle_number_in_subdomain=0

        Full_subdomain_particle_name=0

        !-----------------------------------------------------------------------------------------------------
        !Four corners in actual subdomain (anti-clockwise)
        !LFC---RFC---RBC---LBC
        do i=1,Particle_number_in_LFC_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_LFC_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_LFC_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_LFC_Of_ActualSubdomain(i)

        enddo

        do i=1,Particle_number_in_RFC_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_RFC_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_RFC_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_RFC_Of_ActualSubdomain(i)

        enddo

        do i=1,Particle_number_in_RBC_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_RBC_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_RBC_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_RBC_Of_ActualSubdomain(i)

        enddo

        do i=1,Particle_number_in_LBC_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_LBC_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_LBC_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_LBC_Of_ActualSubdomain(i)

        enddo
        !-----------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------
        !Four sides in actual subdomain (anti-clockwise)
        !FS---RS---BS---LS
        do i=1,Particle_number_in_FS_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_FS_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_FS_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_FS_Of_ActualSubdomain(i)

        enddo

        do i=1,Particle_number_in_RS_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_RS_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_RS_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_RS_Of_ActualSubdomain(i)

        enddo

        do i=1,Particle_number_in_BS_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_BS_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_BS_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_BS_Of_ActualSubdomain(i)

        enddo

        do i=1,Particle_number_in_LS_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_LS_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_LS_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_LS_Of_ActualSubdomain(i)

        enddo
        !-----------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------
        !Center part in subdomain
        do i=1,Particle_number_in_Center_Of_ActualSubdomain

            Actual_particle_number_in_subdomain=Actual_particle_number_in_subdomain+1
            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=Particle_name_in_Center_Of_ActualSubdomain(i)
            actual_particle_name_in_subdomain(Actual_particle_number_in_subdomain)=Particle_name_in_Center_Of_ActualSubdomain(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=Particle_name_in_Center_Of_ActualSubdomain(i)

        enddo
        !-----------------------------------------------------------------------------------------------------


        !-----------------------------------------------------------------------------------------------------
        !Four corners in the buffer zone ( LFC---RFC---RBC---LBC  anti-clockwise)
        do i=1,particle_number_in_LFC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_LFC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_LFC_buffer_zone(i)

        enddo

        do i=1,particle_number_in_RFC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_RFC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_RFC_buffer_zone(i)

        enddo

        do i=1,particle_number_in_RBC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_RBC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_RBC_buffer_zone(i)

        enddo

        do i=1,particle_number_in_LBC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_LBC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_LBC_buffer_zone(i)

        enddo
        !-----------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------
        !Four zone near corners in the buffer zone, star from the right LFC (anti-clockwise)
        do i=1,particle_number_in_Right_LFC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_Right_LFC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_Right_LFC_buffer_zone(i)

        enddo

        do i=1,particle_number_in_Back_RFC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_Back_RFC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_Back_RFC_buffer_zone(i)

        enddo

        do i=1,particle_number_in_Left_RBC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_Left_RBC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_Left_RBC_buffer_zone(i)

        enddo

        do i=1,particle_number_in_Front_LBC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_Front_LBC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_Front_LBC_buffer_zone(i)

        enddo
        !-----------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------
        !Four zone near corners in the buffer zone, star from the back LFC (anti-clockwise)
        do i=1,particle_number_in_Back_LFC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_Back_LFC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_Back_LFC_buffer_zone(i)

        enddo

        do i=1,particle_number_in_Left_RFC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_Left_RFC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_Left_RFC_buffer_zone(i)

        enddo

        do i=1,particle_number_in_Front_RBC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_Front_RBC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_Front_RBC_buffer_zone(i)

        enddo

        do i=1,particle_number_in_Right_LBC_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_Right_LBC_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_Right_LBC_buffer_zone(i)

        enddo
        !-----------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------
        !Four sides in the buffer zone, star from the FS (anti-clockwise)
        do i=1,particle_number_in_FS_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_FS_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_FS_buffer_zone(i)

        enddo

        do i=1,particle_number_in_RS_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_RS_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_RS_buffer_zone(i)

        enddo

        do i=1,particle_number_in_BS_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_BS_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_BS_buffer_zone(i)

        enddo

        do i=1,particle_number_in_LS_buffer_zone

            Total_particle_number_in_subdomain=Total_particle_number_in_subdomain+1

            Full_subdomain_particle_name(Total_particle_number_in_subdomain)=particle_name_in_LS_buffer_zone(i)
            total_particle_name_in_subdomain(Total_particle_number_in_subdomain)=particle_name_in_LS_buffer_zone(i)

        enddo
        !-----------------------------------------------------------------------------------------------------

        !=====================================================================================================

        !Calculate the message length and memory index send to left and right
        !Calculate the message length and memory index received from left and right

        !-----------------------------------------------------------------------------------------------------
        !Message send to the other processors (in actual domain)
        !Four corners (LFC---RFC---RBC---LBC)
        LFC_message_FPI=1
        LFC_message_length=Particle_number_in_LFC_Of_ActualSubdomain

        RFC_message_FPI=Particle_number_in_LFC_Of_ActualSubdomain+1
        RFC_message_length=Particle_number_in_RFC_Of_ActualSubdomain

        RBC_message_FPI=Particle_number_in_LFC_Of_ActualSubdomain+Particle_number_in_RFC_Of_ActualSubdomain+1
        RBC_message_length=Particle_number_in_RBC_Of_ActualSubdomain

        LBC_message_FPI=Particle_number_in_LFC_Of_ActualSubdomain+Particle_number_in_RFC_Of_ActualSubdomain+Particle_number_in_RBC_Of_ActualSubdomain+1
        LBC_message_length=Particle_number_in_LBC_Of_ActualSubdomain

        total_particle_number_in_FourCorners_Of_ActualSubdomain=Particle_number_in_LFC_Of_ActualSubdomain+Particle_number_in_RFC_Of_ActualSubdomain+Particle_number_in_RBC_Of_ActualSubdomain+Particle_number_in_LBC_Of_ActualSubdomain

        Before_Particle_Number=Particle_number_in_LFC_Of_ActualSubdomain+Particle_number_in_RFC_Of_ActualSubdomain+Particle_number_in_RBC_Of_ActualSubdomain+Particle_number_in_LBC_Of_ActualSubdomain

        !Four sides (FS---RS---BS---LS)
        side_message_FPI_to_front=Before_Particle_Number+1
        side_message_length_to_front=Particle_number_in_FS_Of_ActualSubdomain

        side_message_FPI_to_right=Before_Particle_Number+Particle_number_in_FS_Of_ActualSubdomain+1
        side_message_length_to_right=Particle_number_in_RS_Of_ActualSubdomain

        side_message_FPI_to_back=Before_Particle_Number+Particle_number_in_FS_Of_ActualSubdomain+Particle_number_in_RS_Of_ActualSubdomain+1
        side_message_length_to_back=Particle_number_in_BS_Of_ActualSubdomain

        side_message_FPI_to_left=Before_Particle_Number+Particle_number_in_FS_Of_ActualSubdomain+Particle_number_in_RS_Of_ActualSubdomain+Particle_number_in_BS_Of_ActualSubdomain+1
        side_message_length_to_left=Particle_number_in_LS_Of_ActualSubdomain

        !-----------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------
        !Receive message first particle index and the length (Buffer zone)

        Before_Particle_Number=Actual_particle_number_in_subdomain

        !Four corners in buffer zone (LFC---RFC---RBC---LBC)
        Buffer_LFC_message_FPI_from_leftfront=Before_Particle_Number+1
        Buffer_LFC_message_length_from_leftfront=particle_number_in_LFC_buffer_zone

        Buffer_RFC_message_FPI_from_rightfront=Before_Particle_Number+particle_number_in_LFC_buffer_zone+1
        Buffer_RFC_message_length_from_rightfront=particle_number_in_RFC_buffer_zone

        Buffer_RBC_message_FPI_from_rightback=Before_Particle_Number+particle_number_in_LFC_buffer_zone+particle_number_in_RFC_buffer_zone+1
        Buffer_RBC_message_length_from_rightback=particle_number_in_RBC_buffer_zone

        Buffer_LBC_message_FPI_from_leftback=Before_Particle_Number+particle_number_in_LFC_buffer_zone+particle_number_in_RFC_buffer_zone+particle_number_in_RBC_buffer_zone+1
        Buffer_LBC_message_length_from_leftback=particle_number_in_LBC_buffer_zone

        total_particle_number_in_FourCorners_Of_BufferDomain=particle_number_in_LFC_buffer_zone+particle_number_in_RFC_buffer_zone+particle_number_in_RBC_buffer_zone+particle_number_in_LBC_buffer_zone

        Before_Particle_Number=Before_Particle_Number+particle_number_in_LFC_buffer_zone+particle_number_in_RFC_buffer_zone+particle_number_in_RBC_buffer_zone+particle_number_in_LBC_buffer_zone

        !Four near corners in buffer zone (star from Right LFC )
        Buffer_LFC_message_FPI_from_front=Before_Particle_Number+1
        Buffer_LFC_message_length_from_front=particle_number_in_Right_LFC_buffer_zone

        Buffer_RFC_message_FPI_from_right=Before_Particle_Number+particle_number_in_Right_LFC_buffer_zone+1
        Buffer_RFC_message_length_from_right=particle_number_in_Back_RFC_buffer_zone

        Buffer_RBC_message_FPI_from_back=Before_Particle_Number+particle_number_in_Back_RFC_buffer_zone+particle_number_in_Right_LFC_buffer_zone+1
        Buffer_RBC_message_length_from_back=particle_number_in_Left_RBC_buffer_zone

        Buffer_LBC_message_FPI_from_left=Before_Particle_Number+particle_number_in_Left_RBC_buffer_zone+particle_number_in_Back_RFC_buffer_zone+particle_number_in_Right_LFC_buffer_zone+1
        Buffer_LBC_message_length_from_left=particle_number_in_Front_LBC_buffer_zone


        Before_Particle_Number=Before_Particle_Number+particle_number_in_Right_LFC_buffer_zone+particle_number_in_Back_RFC_buffer_zone+particle_number_in_Left_RBC_buffer_zone+particle_number_in_Front_LBC_buffer_zone

        !Four near corners in buffer zone (star from Back LFC )
        Buffer_LFC_message_FPI_from_left=Before_Particle_Number+1
        Buffer_LFC_message_length_from_left=particle_number_in_Back_LFC_buffer_zone

        Buffer_RFC_message_FPI_from_front=Before_Particle_Number+particle_number_in_Back_LFC_buffer_zone+1
        Buffer_RFC_message_length_from_front=particle_number_in_Left_RFC_buffer_zone
      
        Buffer_RBC_message_FPI_from_right=Before_Particle_Number+particle_number_in_Left_RFC_buffer_zone+particle_number_in_Back_LFC_buffer_zone+1
        Buffer_RBC_message_length_from_right=particle_number_in_Front_RBC_buffer_zone
        
        Buffer_LBC_message_FPI_from_back=Before_Particle_Number+particle_number_in_Front_RBC_buffer_zone+particle_number_in_Left_RFC_buffer_zone+particle_number_in_Back_LFC_buffer_zone+1
        Buffer_LBC_message_length_from_back=particle_number_in_Right_LBC_buffer_zone
    
        
        Before_Particle_Number=Before_Particle_Number+particle_number_in_Back_LFC_buffer_zone+particle_number_in_Left_RFC_buffer_zone+particle_number_in_Front_RBC_buffer_zone+particle_number_in_Right_LBC_buffer_zone

        !Four side message in the buffer zone, star from the FS (anti-clockwise)

        side_message_FPI_from_front=Before_Particle_Number+1
        side_message_length_from_front=particle_number_in_FS_buffer_zone

        side_message_FPI_from_right=Before_Particle_Number+particle_number_in_FS_buffer_zone+1
        side_message_length_from_right=particle_number_in_RS_buffer_zone

        side_message_FPI_from_back=Before_Particle_Number+particle_number_in_RS_buffer_zone+particle_number_in_FS_buffer_zone+1
        side_message_length_from_back=particle_number_in_BS_buffer_zone

        side_message_FPI_from_left=Before_Particle_Number+particle_number_in_BS_buffer_zone+particle_number_in_RS_buffer_zone+particle_number_in_FS_buffer_zone+1
        side_message_length_from_left=particle_number_in_LS_buffer_zone
        !-----------------------------------------------------------------------------------------------------

    endif DevideBlock_if

    !=========================================================================================================

    ! Initialized the variables in subdomain
    !---------------------------------------------------------------------------------------------------------
    do k=1,Total_particle_number_in_subdomain

        i=Full_subdomain_particle_name(k)

        subdomain_particle_position(k,:)=particle_position(i,:)
        subdomain_particle_velocity(k,:)=particle_velocity(i,:)
        subdomain_particle_mass(k)=particle_mass(i)                          
        subdomain_particle_type(k)=particle_type(i)                                 
        subdomain_particle_initial_type(k)=particle_initial_type(i)        ! Particle initial type is used for body force calculation and body motion solver                       
        subdomain_particle_rho(k)=particle_rho(i)                                                 
        subdomain_particle_smooth_lengh(k)=particle_smooth_lengh(i)                    
        subdomain_particle_c(k)=particle_c(i)                               
        subdomain_particle_press(k)=particle_press(i)                           
        subdomain_particle_energy(k)=particle_energy(i)                                                     
        subdomain_free_surface_type(k)=free_surface_type(i)                             
        subdomain_Boundary_particle_type(k)=Boundary_particle_type(i)
        subdomain_particle_division_degree(k)=particle_division_degree(i)

        subdomain_wave_maker_particle_layer(k)=wave_maker_particle_layer(i)
        subdomain_In_WaveMakerDampingZone_OrNot(k)=In_WaveMakerDampingZone_OrNot(i)
        subdomain_WaveMaker_Analytical_position(k,:)=WaveMaker_Analytical_position(i,:)      
        subdomain_WaveMaker_Analytical_velocity(k,:)=WaveMaker_Analytical_velocity(i,:)     
        subdomain_WaveMaker_Analytical_rho(k)=WaveMaker_Analytical_rho(i)              
        subdomain_WaveMaker_Analytical_press(k)=WaveMaker_Analytical_press(i)
        subdomain_position_soft_coefficient(k)=position_soft_coefficient(i) 

    end do
    !-----------------------------------------------------------------------------------------------------


    !-----------------------------------------------------------------------------------------------------
    ! Code For Debugging
    ! do k=1,Actual_particle_number_in_subdomain

    !     i=Full_subdomain_particle_name(k)

    !     subdomain_particle_position(k,:)=particle_position(i,:)
    !     subdomain_particle_velocity(k,:)=particle_velocity(i,:)
    !     subdomain_particle_mass(k)=particle_mass(i)                          
    !     subdomain_particle_type(k)=particle_type(i)                                 
    !     subdomain_particle_initial_type(k)=particle_initial_type(i)                         
    !     subdomain_particle_rho(k)=particle_rho(i)                                                 
    !     subdomain_particle_smooth_lengh(k)=particle_smooth_lengh(i)                    
    !     subdomain_particle_c(k)=particle_c(i)                               
    !     subdomain_particle_press(k)=particle_press(i)                           
    !     subdomain_particle_energy(k)=particle_energy(i)                                                    
    !     subdomain_free_surface_type(k)=free_surface_type(i)                             
    !     subdomain_wave_maker_particle_layer(k)=wave_maker_particle_layer(i)
    !     subdomain_Boundary_particle_type(k)=Boundary_particle_type(i)
        
    !     subdomain_In_WaveMakerDampingZone_OrNot(k)=In_WaveMakerDampingZone_OrNot(i)
    !     subdomain_WaveMaker_Analytical_position(k,:)=WaveMaker_Analytical_position(i,:)      
    !     subdomain_WaveMaker_Analytical_velocity(k,:)=WaveMaker_Analytical_velocity(i,:)     
    !     subdomain_WaveMaker_Analytical_rho(k)=WaveMaker_Analytical_rho(i)              
    !     subdomain_WaveMaker_Analytical_press(k)=WaveMaker_Analytical_press(i)
    !     subdomain_position_soft_coefficient(k)=position_soft_coefficient(i) 

    ! end do


    !write(*,*) Buffer_LFC_message_length_from_front,Buffer_RFC_message_length_from_front,Buffer_LBC_message_length_from_back,Buffer_RBC_message_length_from_back,Current_Processor_ID

    ! call MPI_RealVector_MessageExchange(subdomain_particle_position)    ! Exchange the position
    ! call MPI_RealVector_MessageExchange(subdomain_particle_velocity)    ! Exchange the velocity
   
    ! call MPI_RealScalar_MessageExchange(subdomain_particle_press)       ! Exchange the pressure
    ! call MPI_RealScalar_MessageExchange(subdomain_particle_rho)         ! Exchange the density
    !-----------------------------------------------------------------------------------------------------

    !=====================================================================================================
















!     !*****************************************************************************************************
!     ! For debuging: check the particle dividing results
!     if ( i_time_step==Actual_start_time_step ) then

!         !-------------------------------------------------------------------------------------------------
!         ! Open the output tecplot data files
!         Standby_Current_Processor_File_Index = Standby_Initial_File_Port_For_Each_Processor+Current_Processor_ID

!         ! Transfer the Current_Processor_ID from integer to character
!         write(Char_Current_Processor_ID,'(I4)') Current_Processor_ID
!         !write(*,*) Char_Current_Processor_ID

!         File_index = Standby_Current_Processor_File_Index
!         File_name  = "./Subdomain/Total_Particle_in_subdomain_"//trim(adjustl(Char_Current_Processor_ID))//".dat"

!         call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
   
!         if ( IOERROR==0 ) then

!             !---------------------------------------------------------------------------------------------
!             ! Tecplot header
!             write(File_index,*) "TITLE='DISTRIBUTION'"
!             write(File_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
!             !---------------------------------------------------------------------------------------------

!             !---------------------------------------------------------------------------------------------
!             ! Output the result
!             write(File_index,*) "ZONE I=",Total_particle_number_in_subdomain," F=POINT"
!             do i=1,Total_particle_number_in_subdomain
!                 write(File_index,100) (subdomain_particle_position(i,j),j=1,dim),(subdomain_particle_velocity(i,j),j=1,dim),subdomain_particle_rho(i),subdomain_particle_press(i)/pressure_0
! 100             format(8F20.10) 
!             end do
!             !---------------------------------------------------------------------------------------------

!         endif
!         close(File_index)


!         !-------------------------------------------------------------------------------------------------


!         !-------------------------------------------------------------------------------------------------
!         ! Open the output tecplot data files
!         File_index = Standby_Current_Processor_File_Index
!         File_name  = "./Subdomain/Actual_Partcile_in_subdomain_"//trim(adjustl(Char_Current_Processor_ID))//".dat"

!         call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
        
!         if ( IOERROR==0 ) then

!             !---------------------------------------------------------------------------------------------
!             ! Tecplot header
!             write(File_index,*) "TITLE='DISTRIBUTION'"
!             write(File_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
!             !---------------------------------------------------------------------------------------------

!             !---------------------------------------------------------------------------------------------
!             ! Output the result
!             write(File_index,*) "ZONE I=",Actual_particle_number_in_subdomain," F=POINT"
!             do i=1,Actual_particle_number_in_subdomain
!                 write(File_index,100) (subdomain_particle_position(i,j),j=1,dim),(subdomain_particle_velocity(i,j),j=1,dim),subdomain_particle_rho(i),subdomain_particle_press(i)/pressure_0
!             end do
!             !---------------------------------------------------------------------------------------------

!         endif
!         close(File_index)

!         !-------------------------------------------------------------------------------------------------



!         !-------------------------------------------------------------------------------------------------
!         ! Open the output tecplot data files
!         File_index = Standby_Current_Processor_File_Index
!         File_name  = "./Subdomain/Actual_Partcile_in_subdomain_"//trim(adjustl(Char_Current_Processor_ID))//".dat"

!         call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
        
!         if ( IOERROR==0 ) then

!             !---------------------------------------------------------------------------------------------
!             ! Tecplot header
!             write(File_index,*) "TITLE='DISTRIBUTION'"
!             write(File_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
!             !---------------------------------------------------------------------------------------------

!             !---------------------------------------------------------------------------------------------
!             ! Output the result
!             write(File_index,*) "ZONE I=",Total_particle_number_in_subdomain," F=POINT"
!             do k=1,Total_particle_number_in_subdomain
!                 i=Full_subdomain_particle_name(k)
!                 write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
!             end do
!             !---------------------------------------------------------------------------------------------

!         endif
!         close(File_index)

!         !-------------------------------------------------------------------------------------------------


!         ! Synchronize all processors calculation           
!         call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

!     endif
!     !*****************************************************************************************************



end subroutine Divide_Block_and_Allot_Tasks

