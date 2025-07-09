!**************************************************************************************************************
!  Subroutine: Distribute_Backgroud_Grid (Have checked, this is right)
!
!  PURPOSE: Get the near grid name of each grid.
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
!  Input: chain_x_number,chain_y_number,chain_z_number,grid_max_number=chain_x_number*chain_y_number*chain_z_number
!
!  Output: mesh_center_positon,chain_near_mesh_name,chain_near_mesh_number
!
!  Note: MPI version: mpich-3.2
!        Fortran version: Inter Fortran (ifort) 18.0.0
!**************************************************************************************************************
subroutine Distribute_Backgroud_Grid(dim,&                                       ! Dimension of grid
                                     grid_x_number,&                             ! Grid number in x direction
                                     grid_y_number,&                             ! Grid number in y direction
                                     grid_z_number,&                             ! Grid number in z direction
                                     grid_max_number,&                           ! Max grid number (grid_x_number*grid_y_number*grid_z_number)
                                     grid_center_positon,&                       ! Grid center positon
                                     grid_near_mesh_number,&                     ! Near mesh number
                                     grid_near_mesh_name&                        ! grid near mesh name
                                     )

    USE MPI

    implicit none

    !==========================================================================================================
    
    !----------------------------------------------------------------------------------------------------------
    !Variables from superior subroutine
    integer,intent(in)::dim                                                      ! Dimension of grid
    integer,intent(in)::grid_x_number,grid_y_number,grid_z_number                ! Grid number in x,y,z direction
    integer,intent(in)::grid_max_number                                          ! MAX GRID NUMBER (grid_x_number*grid_y_number*grid_z_number)
    integer,dimension(grid_max_number,3),intent(out)::grid_center_positon        ! Near mesh number
    integer,dimension(grid_max_number),intent(out)::grid_near_mesh_number        ! Near mesh number
    integer,dimension(grid_max_number,27),intent(out)::grid_near_mesh_name       ! Near mesh name
    !----------------------------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------------------------
    !Variables in local subroutine
    integer::i,j,k,L,m                                                           ! Loop Variables
    integer::mesh_x,mesh_y,mesh_z                                                ! Temp varibales for backgroud mesh

    ! Variables of the MPI Processor for checking
    integer::ierror_MPI                                                          ! MPI Function return value for success runing or not
    integer,dimension(MPI_STATUS_SIZE)::Status_MPI                               ! MPI running status
    integer::Current_Processor_ID                                                ! Current Processor ID
    integer::Total_Processors_Number                                             ! Total Processors Number (0,1,2,...,Total_Processors_Number-1)

    !----------------------------------------------------------------------------------------------------------

    !==========================================================================================================
    !For 2-D case
    Dimension_If: if (dim==2) then

        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! Distribute the background grid
        ! write(*,*) chain_dx,chain_dy,chain_x_number,chain_y_number
        ! Calculate the center of background grid cell
        do j=1,grid_y_number
            do i=1,grid_x_number
                k=i+(j-1)*grid_x_number                                          ! Row by row
                grid_center_positon(k,1)=i
                grid_center_positon(k,2)=j 
            end do
        end do
        
        ! Initialized the Grid near mesh name(index) and number
        grid_near_mesh_number=0                    
        grid_near_mesh_name=0                      
        
        do i=1,grid_max_number
            
            ! Get the particle center position
            mesh_x=grid_center_positon(i,1)
            mesh_y=grid_center_positon(i,2)
            

            ! Grid in the center domain
            if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number ) then
                
                grid_near_mesh_number(i)=9                                       ! The cell in the center has NINE vicinity grids

                grid_near_mesh_name(i,1)=i                                       ! The first is itself
                grid_near_mesh_name(i,2)=i-1                                     ! The second is left vicinity grid 
                grid_near_mesh_name(i,3)=i+1                                     ! The third is right vicinity grid 
                grid_near_mesh_name(i,4)=i-1+grid_x_number                       ! The fourth is left top vicinity grid 
                grid_near_mesh_name(i,5)=i+grid_x_number                         ! The fifth is top vicinity grid 
                grid_near_mesh_name(i,6)=i+1+grid_x_number                       ! The sixth is right top vicinity grid 
                grid_near_mesh_name(i,7)=i-1-grid_x_number                       ! The seventh is left bottom vicinity grid
                grid_near_mesh_name(i,8)=i-grid_x_number                         ! The eighth is bottom vicinity grid 
                grid_near_mesh_name(i,9)=i+1-grid_x_number                       ! The ninth is right bottom vicinity grid 
                
            ! Left side grid (Not include the corner)
            else if(mesh_x==1 .and. mesh_y>1 .and. mesh_y<grid_y_number) then
                    
                grid_near_mesh_number(i)=6                                       ! Side cell has SIX vicinity grids
            
                grid_near_mesh_name(i,1)=i                                       ! The first is itself
                grid_near_mesh_name(i,2)=i+1                                     ! The second is right vicinity grid 
                grid_near_mesh_name(i,3)=i+grid_x_number                         ! The third is bottom vicinity grid 
                grid_near_mesh_name(i,4)=i+1+grid_x_number                       ! The fourth is right bottom vicinity grid
                grid_near_mesh_name(i,5)=i-grid_x_number                         ! The fifth is top vicinity grid 
                grid_near_mesh_name(i,6)=i+1-grid_x_number                       ! The sixth is right top vicinity grid 
 
            ! Right side grid (Not include the corner)
            else if(mesh_x==grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number) then
                    
                grid_near_mesh_number(i)=6                                       ! Side cell has SIX vicinity grids
            
                grid_near_mesh_name(i,1)=i                                       ! The first is itself
                grid_near_mesh_name(i,2)=i-1                                     ! The second is left vicinity grid 
                grid_near_mesh_name(i,3)=i+grid_x_number                         ! The third is top vicinity grid
                grid_near_mesh_name(i,4)=i-1+grid_x_number                       ! The fourth is left top vicinity grid 
                grid_near_mesh_name(i,5)=i-grid_x_number                         ! The fifth is bottom vicinity grid 
                grid_near_mesh_name(i,6)=i-1-grid_x_number                       ! The sixth is left bottom vicinity grid 
                
            ! Top side grid (Not include the corner)
            else if(mesh_y==grid_y_number.and. mesh_x>1 .and. mesh_x<grid_x_number) then
                    
                grid_near_mesh_number(i)=6                                       ! Side cell has SIX vicinity grids
            
                grid_near_mesh_name(i,1)=i                                       ! The first is itself
                grid_near_mesh_name(i,2)=i-1                                     ! The second is left vicinity grid 
                grid_near_mesh_name(i,3)=i+1                                     ! The third is right vicinity grid 
                grid_near_mesh_name(i,4)=i-1-grid_x_number                       ! The fourth is left bottom vicinity grid 
                grid_near_mesh_name(i,5)=i-grid_x_number                         ! The fifth is bottom vicinity grid
                grid_near_mesh_name(i,6)=i+1-grid_x_number                       ! The sixth is right bottom vicinity grid
                 
            ! Bottom side grid (Not include the corner)
            else if(mesh_y==1.and. mesh_x>1 .and. mesh_x<grid_x_number) then
                    
                grid_near_mesh_number(i)=6                                       ! Side cell has SIX vicinity grids
            
                grid_near_mesh_name(i,1)=i                                       ! The first is itself
                grid_near_mesh_name(i,2)=i-1                                     ! The second is left vicinity grid 
                grid_near_mesh_name(i,3)=i+1                                     ! The third is right vicinity grid 
                grid_near_mesh_name(i,4)=i-1+grid_x_number                       ! The fourth is left top vicinity grid 
                grid_near_mesh_name(i,5)=i+grid_x_number                         ! The fifth is top vicinity grid 
                grid_near_mesh_name(i,6)=i+1+grid_x_number                       ! The sixth is right top vicinity grid
                
            ! Four corners
            ! Left-top corner
            else if(mesh_x==1 .and. mesh_y==grid_y_number) then    
                    
                grid_near_mesh_number(i)=4                                       ! Corner cell has FOUR vicinity grids
            
                grid_near_mesh_name(i,1)=i                                       ! The first is itself
                grid_near_mesh_name(i,2)=i+1                                     ! The second is right vicinity grid 
                grid_near_mesh_name(i,3)=i-grid_x_number                         ! The third is bottom vicinity grid 
                grid_near_mesh_name(i,4)=i+1-grid_x_number                       ! The fourth is right top vicinity grid 
                
            ! Right-top corner
            else if(mesh_x==grid_x_number .and. mesh_y==grid_y_number) then    
                    
                grid_near_mesh_number(i)=4                                       ! Corner cell has FOUR vicinity grids
            
                grid_near_mesh_name(i,1)=i                                       ! The first is itself
                grid_near_mesh_name(i,2)=i-1                                     ! The second is left vicinity grid 
                grid_near_mesh_name(i,3)=i-grid_x_number                         ! The third is bottom vicinity grid 
                grid_near_mesh_name(i,4)=i-1-grid_x_number                       ! The fourth is left bottom vicinity grid
                
            ! Left-bottom corner
            else if(mesh_x==1 .and. mesh_y==1) then    
                    
                grid_near_mesh_number(i)=4                                       ! Corner cell has FOUR vicinity grids
            
                grid_near_mesh_name(i,1)=i                                       ! The first is itself
                grid_near_mesh_name(i,2)=i+1                                     ! The second is right vicinity grid 
                grid_near_mesh_name(i,3)=i+grid_x_number                         ! The third is top vicinity grid 
                grid_near_mesh_name(i,4)=i+1+grid_x_number                       ! The fourth is right top vicinity grid 

            ! Right-bottom corner    
            else if(mesh_x==grid_x_number .and. mesh_y==1) then    
                    
                grid_near_mesh_number(i)=4                                       ! Corner cell has FOUR vicinity grids
            
                grid_near_mesh_name(i,1)=i                                       ! The first is itself
                grid_near_mesh_name(i,2)=i-1                                     ! The second is left vicinity grid 
                grid_near_mesh_name(i,3)=i+grid_x_number                         ! The third is top vicinity grid 
                grid_near_mesh_name(i,4)=i-1+grid_x_number                       ! The fourth is left top vicinity grid 
                
            end if  
            
        end do 
        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    !For 3-D case
    else if (dim==3) then
        
        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! Distribute background grid in 3D
        ! Calculate the center of background grid cell
        do k=1,grid_z_number
            do j=1,grid_y_number
                do i=1,grid_x_number
                    
                    L=i+(j-1)*grid_x_number+(k-1)*grid_x_number*grid_y_number    ! Row by row
                    grid_center_positon(L,1)=i
                    grid_center_positon(L,2)=j 
                    grid_center_positon(L,3)=k 

                end do
            end do
        end do
        
        ! Initialized the chain near mesh name(index) and number
        grid_near_mesh_number=0                    
        grid_near_mesh_name=0                      
        
        do i=1,grid_max_number
            
            ! Get the particle center position
            mesh_x=grid_center_positon(i,1)
            mesh_y=grid_center_positon(i,2)
            mesh_z=grid_center_positon(i,3)

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Cell in the center of the grid block
            if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then
                
                grid_near_mesh_number(i)=27                                                   ! The cell in the center has 27 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The 1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i-1                                                  ! The 2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i+1                                                  ! The 3rd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i-1+grid_x_number                                    ! The 4th is left  top    vicinity grid in Local layer
                grid_near_mesh_name(i,5)=i  +grid_x_number                                    ! The 5th is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,6)=i+1+grid_x_number                                    ! The 6th is right top    vicinity grid in Local layer
                grid_near_mesh_name(i,7)=i-1-grid_x_number                                    ! The 7th is left  bottom vicinity grid in Local layer
                grid_near_mesh_name(i,8)=i  -grid_x_number                                    ! The 8th is       bottom vicinity grid in Local layer
                grid_near_mesh_name(i,9)=i+1-grid_x_number                                    ! The 9th is right bottom vicinity grid in Local layer

                ! Upper layper 
                grid_near_mesh_name(i,10)=i                +grid_x_number*grid_y_number       ! The 10th is itself                     in Upper layer
                grid_near_mesh_name(i,11)=i-1              +grid_x_number*grid_y_number       ! The 11th is left         vicinity grid in Upper layer
                grid_near_mesh_name(i,12)=i+1              +grid_x_number*grid_y_number       ! The 12th is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,13)=i-1+grid_x_number+grid_x_number*grid_y_number       ! The 13th is left  top    vicinity grid in Upper layer
                grid_near_mesh_name(i,14)=i  +grid_x_number+grid_x_number*grid_y_number       ! The 14th is       top    vicinity grid in Upper layer
                grid_near_mesh_name(i,15)=i+1+grid_x_number+grid_x_number*grid_y_number       ! The 15th is right top    vicinity grid in Upper layer
                grid_near_mesh_name(i,16)=i-1-grid_x_number+grid_x_number*grid_y_number       ! The 16th is left  bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,17)=i  -grid_x_number+grid_x_number*grid_y_number       ! The 17th is       bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,18)=i+1-grid_x_number+grid_x_number*grid_y_number       ! The 18th is right bottom vicinity grid in Upper layer

                ! Lower layer
                grid_near_mesh_name(i,19)=i                -grid_x_number*grid_y_number       ! The 19th is itself                     in Lower layer
                grid_near_mesh_name(i,20)=i-1              -grid_x_number*grid_y_number       ! The 20th is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,21)=i+1              -grid_x_number*grid_y_number       ! The 21th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,22)=i-1+grid_x_number-grid_x_number*grid_y_number       ! The 22th is left  top    vicinity grid in Lower layer
                grid_near_mesh_name(i,23)=i  +grid_x_number-grid_x_number*grid_y_number       ! The 23th is       top    vicinity grid in Lower layer
                grid_near_mesh_name(i,24)=i+1+grid_x_number-grid_x_number*grid_y_number       ! The 24th is right top    vicinity grid in Lower layer
                grid_near_mesh_name(i,25)=i-1-grid_x_number-grid_x_number*grid_y_number       ! The 25th is left  bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,26)=i  -grid_x_number-grid_x_number*grid_y_number       ! The 26th is       bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,27)=i+1-grid_x_number-grid_x_number*grid_y_number       ! The 27th is right bottom vicinity grid in Lower layer


            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Six faces of the grid block (Not contain 12 edges and 8 corners ) 
            
            !--------------------------------------------------------------------------------------------------
            ! Left face (Not contain 4 edges and 4 corners ) 
            else if(mesh_x==1 .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then
                
                grid_near_mesh_number(i)=18                                                   ! Cell on face has has 18 vicinity grids

                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The 1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i+1                                                 ! The 2nd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i  +grid_x_number                                   ! The 3rd is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i+1+grid_x_number                                   ! The 4th is right top    vicinity grid in Local layer
                grid_near_mesh_name(i,5) =i  -grid_x_number                                   ! The 5th is       bottom vicinity grid in Local layer
                grid_near_mesh_name(i,6) =i+1-grid_x_number                                   ! The 6th is right bottom vicinity grid in Local layer

                ! Upper layer 
                grid_near_mesh_name(i,7) =i                +grid_x_number*grid_y_number       ! The 7th  is itself                     in Upper layer
                grid_near_mesh_name(i,8) =i+1              +grid_x_number*grid_y_number       ! The 8th  is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,9) =i  +grid_x_number+grid_x_number*grid_y_number       ! The 9rd  is       top    vicinity grid in Upper layer
                grid_near_mesh_name(i,10)=i+1+grid_x_number+grid_x_number*grid_y_number       ! The 10th is right top    vicinity grid in Upper layer
                grid_near_mesh_name(i,11)=i  -grid_x_number+grid_x_number*grid_y_number       ! The 11th is       bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,12)=i+1-grid_x_number+grid_x_number*grid_y_number       ! The 12th is right bottom vicinity grid in Upper layer

                ! Lower layer
                grid_near_mesh_name(i,13)=i                -grid_x_number*grid_y_number       ! The 13th is itself                     in Lower layer
                grid_near_mesh_name(i,14)=i+1              -grid_x_number*grid_y_number       ! The 14th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,15)=i  +grid_x_number-grid_x_number*grid_y_number       ! The 15th is       top    vicinity grid in Lower layer
                grid_near_mesh_name(i,16)=i+1+grid_x_number-grid_x_number*grid_y_number       ! The 16th is right top    vicinity grid in Lower layer
                grid_near_mesh_name(i,17)=i  -grid_x_number-grid_x_number*grid_y_number       ! The 17th is       bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,18)=i+1-grid_x_number-grid_x_number*grid_y_number       ! The 18th is right bottom vicinity grid in Lower layer

            !--------------------------------------------------------------------------------------------------
            ! Right face (Not contain 4 edges and 4 corners ) 
            else if(mesh_x==grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=18                                                   ! Cell on face has has 18 vicinity grids

                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                    in Local layer
                grid_near_mesh_name(i,2) =i-1                                                 ! The  2nd is left        vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i  +grid_x_number                                   ! The  3rd is      top    vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i-1+grid_x_number                                   ! The  4th is left top    vicinity grid in Local layer
                grid_near_mesh_name(i,5) =i  -grid_x_number                                   ! The  5th is      bottom vicinity grid in Local layer
                grid_near_mesh_name(i,6) =i-1-grid_x_number                                   ! The  6th is left bottom vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,7) =i                +grid_x_number*grid_y_number       ! The  7th is itself                    in Upper layer
                grid_near_mesh_name(i,8) =i-1              +grid_x_number*grid_y_number       ! The  8th is left        vicinity grid in Upper layer
                grid_near_mesh_name(i,9) =i  +grid_x_number+grid_x_number*grid_y_number       ! The  9th is      top    vicinity grid in Upper layer
                grid_near_mesh_name(i,10)=i-1+grid_x_number+grid_x_number*grid_y_number       ! The 10th is left top    vicinity grid in Upper layer
                grid_near_mesh_name(i,11)=i  -grid_x_number+grid_x_number*grid_y_number       ! The 11th is      bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,12)=i-1-grid_x_number+grid_x_number*grid_y_number       ! The 12th is left bottom vicinity grid in Upper layer

                ! Lower layer
                grid_near_mesh_name(i,13)=i                -grid_x_number*grid_y_number       ! The 13th is itself                    in Lower layer
                grid_near_mesh_name(i,14)=i-1              -grid_x_number*grid_y_number       ! The 14th is left        vicinity grid in Lower layer
                grid_near_mesh_name(i,15)=i  +grid_x_number-grid_x_number*grid_y_number       ! The 15th is      top    vicinity grid in Lower layer
                grid_near_mesh_name(i,16)=i-1+grid_x_number-grid_x_number*grid_y_number       ! The 16th is left top    vicinity grid in Lower layer
                grid_near_mesh_name(i,17)=i  -grid_x_number-grid_x_number*grid_y_number       ! The 17th is      bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,18)=i-1-grid_x_number-grid_x_number*grid_y_number       ! The 18th is left bottom vicinity grid in Lower layer

            !--------------------------------------------------------------------------------------------------
            ! Front face (Not contain 4 edges and 4 corners ) 
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==1 .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=18                                                   ! Cell on face has has 18 vicinity grids
            
                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i-1                                                  ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i+1                                                  ! The  3rd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i-1+grid_x_number                                    ! The  4th is left  top    vicinity grid in Local layer
                grid_near_mesh_name(i,5)=i  +grid_x_number                                    ! The  5th is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,6)=i+1+grid_x_number                                    ! The  6th is right top    vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,7) =i                +grid_x_number*grid_y_number       ! The  7th is itself                     in Upper layer
                grid_near_mesh_name(i,8) =i-1              +grid_x_number*grid_y_number       ! The  8th is left         vicinity grid in Upper layer
                grid_near_mesh_name(i,9) =i+1              +grid_x_number*grid_y_number       ! The  9th is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,10)=i-1+grid_x_number+grid_x_number*grid_y_number       ! The 10th is left  top    vicinity grid in Upper layer
                grid_near_mesh_name(i,11)=i  +grid_x_number+grid_x_number*grid_y_number       ! The 11th is       top    vicinity grid in Upper layer
                grid_near_mesh_name(i,12)=i+1+grid_x_number+grid_x_number*grid_y_number       ! The 12th is right top    vicinity grid in Upper layer

                ! Lower layer
                grid_near_mesh_name(i,13)=i                -grid_x_number*grid_y_number       ! The 13th is itself                     in Lower layer
                grid_near_mesh_name(i,14)=i-1              -grid_x_number*grid_y_number       ! The 14th is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,15)=i+1              -grid_x_number*grid_y_number       ! The 15th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,16)=i-1+grid_x_number-grid_x_number*grid_y_number       ! The 16th is left  top    vicinity grid in Lower layer
                grid_near_mesh_name(i,17)=i  +grid_x_number-grid_x_number*grid_y_number       ! The 17th is       top    vicinity grid in Lower layer
                grid_near_mesh_name(i,18)=i+1+grid_x_number-grid_x_number*grid_y_number       ! The 18th is right top    vicinity grid in Lower layer


            !--------------------------------------------------------------------------------------------------
            ! Back face (Not contain 4 edges and 4 corners ) 
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then
                
                grid_near_mesh_number(i)=18                                                   ! Cell on face has has 18 vicinity grids
            
                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i-1                                                 ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i+1                                                 ! The  3rd is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,4) =i-1-grid_x_number                                   ! The  4th is left  bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,5) =i  -grid_x_number                                   ! The  5th is       bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,6) =i+1-grid_x_number                                   ! The  6th is right bottom vicinity grid in Lower layer

                ! Upper layer
                grid_near_mesh_name(i,7) =i                +grid_x_number*grid_y_number       ! The  7th is itself                     in Upper layer
                grid_near_mesh_name(i,8) =i-1              +grid_x_number*grid_y_number       ! The  8th is left         vicinity grid in Upper layer
                grid_near_mesh_name(i,9) =i+1              +grid_x_number*grid_y_number       ! The  9th is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,10)=i-1-grid_x_number+grid_x_number*grid_y_number       ! The 10th is left  bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,11)=i  -grid_x_number+grid_x_number*grid_y_number       ! The 11th is       bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,12)=i+1-grid_x_number+grid_x_number*grid_y_number       ! The 12th is right bottom vicinity grid in Upper layer

                ! Lower layer
                grid_near_mesh_name(i,13)=i                -grid_x_number*grid_y_number       ! The 13th is itself                     in Lower layer
                grid_near_mesh_name(i,14)=i-1              -grid_x_number*grid_y_number       ! The 14th is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,15)=i+1              -grid_x_number*grid_y_number       ! The 15th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,16)=i-1-grid_x_number-grid_x_number*grid_y_number       ! The 16th is left  bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,17)=i  -grid_x_number-grid_x_number*grid_y_number       ! The 17th is       bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,18)=i+1-grid_x_number-grid_x_number*grid_y_number       ! The 18th is right bottom vicinity grid in Lower layer


            !--------------------------------------------------------------------------------------------------
            ! Top face (Not contain 4 edges and 4 corners ) 
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==grid_z_number) then

                grid_near_mesh_number(i)=18                                                   ! Cell on face has has 18 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i-1                                                 ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i+1                                                 ! The  3rd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i-1+grid_x_number                                   ! The  4th is left  top    vicinity grid in Local layer
                grid_near_mesh_name(i,5) =i  +grid_x_number                                   ! The  5th is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,6) =i+1+grid_x_number                                   ! The  6th is right top    vicinity grid in Local layer
                grid_near_mesh_name(i,7) =i-1-grid_x_number                                   ! The  7th is left  bottom vicinity grid in Local layer
                grid_near_mesh_name(i,8) =i  -grid_x_number                                   ! The  8th is       bottom vicinity grid in Local layer
                grid_near_mesh_name(i,9) =i+1-grid_x_number                                   ! The  9th is right bottom vicinity grid in Local layer

                ! Lower layer
                grid_near_mesh_name(i,10)=i                -grid_x_number*grid_y_number       ! The 10th is itself                     in Lower layer
                grid_near_mesh_name(i,11)=i-1              -grid_x_number*grid_y_number       ! The 11th is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,12)=i+1              -grid_x_number*grid_y_number       ! The 12th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,13)=i-1+grid_x_number-grid_x_number*grid_y_number       ! The 13th is left  top    vicinity grid in Lower layer
                grid_near_mesh_name(i,14)=i  +grid_x_number-grid_x_number*grid_y_number       ! The 14th is       top    vicinity grid in Lower layer
                grid_near_mesh_name(i,15)=i+1+grid_x_number-grid_x_number*grid_y_number       ! The 15th is right top    vicinity grid in Lower layer
                grid_near_mesh_name(i,16)=i-1-grid_x_number-grid_x_number*grid_y_number       ! The 16th is left  bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,17)=i  -grid_x_number-grid_x_number*grid_y_number       ! The 17th is       bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,18)=i+1-grid_x_number-grid_x_number*grid_y_number       ! The 18th is right bottom vicinity grid in Lower layer

            !--------------------------------------------------------------------------------------------------
            ! Bottom face (Not contain 4 edges and 4 corners ) 
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==1) then

                grid_near_mesh_number(i)=18                                                   ! Cell on face has has 18 vicinity grids

                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i-1                                                 ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i+1                                                 ! The  3rd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i-1+grid_x_number                                   ! The  4th is left  top    vicinity grid in Local layer
                grid_near_mesh_name(i,5) =i  +grid_x_number                                   ! The  5th is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,6) =i+1+grid_x_number                                   ! The  6th is right top    vicinity grid in Local layer
                grid_near_mesh_name(i,7) =i-1-grid_x_number                                   ! The  7th is left  bottom vicinity grid in Local layer
                grid_near_mesh_name(i,8) =i  -grid_x_number                                   ! The  8th is       bottom vicinity grid in Local layer
                grid_near_mesh_name(i,9) =i+1-grid_x_number                                   ! The  9th is right bottom vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,10)=i                +grid_x_number*grid_y_number       ! The 10th is itself                     in Upper layer
                grid_near_mesh_name(i,11)=i-1              +grid_x_number*grid_y_number       ! The 11th is left         vicinity grid in Upper layer
                grid_near_mesh_name(i,12)=i+1              +grid_x_number*grid_y_number       ! The 12th is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,13)=i-1+grid_x_number+grid_x_number*grid_y_number       ! The 13th is left  top    vicinity grid in Upper layer
                grid_near_mesh_name(i,14)=i  +grid_x_number+grid_x_number*grid_y_number       ! The 14th is       top    vicinity grid in Upper layer
                grid_near_mesh_name(i,15)=i+1+grid_x_number+grid_x_number*grid_y_number       ! The 15th is right top    vicinity grid in Upper layer
                grid_near_mesh_name(i,16)=i-1-grid_x_number+grid_x_number*grid_y_number       ! The 16th is left  bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,17)=i  -grid_x_number+grid_x_number*grid_y_number       ! The 17th is       bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,18)=i+1-grid_x_number+grid_x_number*grid_y_number       ! The 18th is right bottom vicinity grid in Upper layer

            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            ! Twelve edges
            ! Four edges on the bottom face
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Left in Four edges on the bottom face
            else if(mesh_x==1 .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i+1                                                 ! The  2nd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i  +grid_x_number                                   ! The  3rd is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i+1+grid_x_number                                   ! The  4th is right top    vicinity grid in Local layer
                grid_near_mesh_name(i,5) =i  -grid_x_number                                   ! The  5th is       bottom vicinity grid in Local layer
                grid_near_mesh_name(i,6) =i+1-grid_x_number                                   ! The  6th is right bottom vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,7) =i                +grid_x_number*grid_y_number       ! The  7th is itself                     in Upper layer
                grid_near_mesh_name(i,8) =i+1              +grid_x_number*grid_y_number       ! The  8th is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,9) =i  +grid_x_number+grid_x_number*grid_y_number       ! The  9th is       top    vicinity grid in Upper layer
                grid_near_mesh_name(i,10)=i+1+grid_x_number+grid_x_number*grid_y_number       ! The 10th is right top    vicinity grid in Upper layer
                grid_near_mesh_name(i,11)=i  -grid_x_number+grid_x_number*grid_y_number       ! The 11th is       bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,12)=i+1-grid_x_number+grid_x_number*grid_y_number       ! The 12th is right bottom vicinity grid in Upper layer
                
            ! Right in Four edges on the bottom face
            else if(mesh_x==grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids

                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                    in Local layer
                grid_near_mesh_name(i,2)=i-1                                                  ! The  2nd is left        vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i-1+grid_x_number                                    ! The  3rh is left top    vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i  +grid_x_number                                    ! The  4th is      top    vicinity grid in Local layer
                grid_near_mesh_name(i,5)=i-1-grid_x_number                                    ! The  5th is left bottom vicinity grid in Local layer
                grid_near_mesh_name(i,6)=i  -grid_x_number                                    ! The  6th is      bottom vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,7) =i                +grid_x_number*grid_y_number       ! The  7th is itself                    in Upper layer
                grid_near_mesh_name(i,8) =i-1              +grid_x_number*grid_y_number       ! The  8th is left        vicinity grid in Upper layer
                grid_near_mesh_name(i,9) =i-1+grid_x_number+grid_x_number*grid_y_number       ! The  9th is left top    vicinity grid in Upper layer
                grid_near_mesh_name(i,10)=i  +grid_x_number+grid_x_number*grid_y_number       ! The 10th is      top    vicinity grid in Upper layer
                grid_near_mesh_name(i,11)=i-1-grid_x_number+grid_x_number*grid_y_number       ! The 11th is left bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,12)=i  -grid_x_number+grid_x_number*grid_y_number       ! The 12th is      bottom vicinity grid in Upper layer

            ! Front in Four edges on the bottom face
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==1 .and. mesh_z==1) then

                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids
                   
                ! Local layer 
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i-1                                                 ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i+1                                                 ! The  3rd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i-1+grid_x_number                                   ! The  4th is left  top    vicinity grid in Local layer
                grid_near_mesh_name(i,5) =i  +grid_x_number                                   ! The  5th is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,6) =i+1+grid_x_number                                   ! The  6th is right top    vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,7) =i                +grid_x_number*grid_y_number       ! The  7th is itself                     in Upper layer
                grid_near_mesh_name(i,8) =i-1              +grid_x_number*grid_y_number       ! The  8th is left         vicinity grid in Upper layer
                grid_near_mesh_name(i,9) =i+1              +grid_x_number*grid_y_number       ! The  9th is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,10)=i-1+grid_x_number+grid_x_number*grid_y_number       ! The 10th is left  top    vicinity grid in Upper layer
                grid_near_mesh_name(i,11)=i  +grid_x_number+grid_x_number*grid_y_number       ! The 11th is       top    vicinity grid in Upper layer
                grid_near_mesh_name(i,12)=i+1+grid_x_number+grid_x_number*grid_y_number       ! The 12th is right top    vicinity grid in Upper layer


            ! Back in Four edges on the bottom face
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==grid_y_number .and. mesh_z==1) then

                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids

                ! Local layer 
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i-1                                                  ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i+1                                                  ! The  3rd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i-1-grid_x_number                                    ! The  4th is left  bottom vicinity grid in Local layer
                grid_near_mesh_name(i,5)=i  -grid_x_number                                    ! The  5th is       bottom vicinity grid in Local layer
                grid_near_mesh_name(i,6)=i+1-grid_x_number                                    ! The  6th is right bottom vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,7) =i                +grid_x_number*grid_y_number       ! The  7th is itself                     in Upper layer
                grid_near_mesh_name(i,8) =i-1              +grid_x_number*grid_y_number       ! The  8th is left         vicinity grid in Upper layer
                grid_near_mesh_name(i,9) =i+1              +grid_x_number*grid_y_number       ! The  9th is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,10)=i-1-grid_x_number+grid_x_number*grid_y_number       ! The 10th is left  bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,11)=i  -grid_x_number+grid_x_number*grid_y_number       ! The 11th is       bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,12)=i+1-grid_x_number+grid_x_number*grid_y_number       ! The 12th is right bottom vicinity grid in Upper layer

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            ! Four edges on the top face
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Left in four edges on the top face
            else if(mesh_x==1 .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids

                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i+1                                                 ! The  2nd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i  +grid_x_number                                   ! The  3rd is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i+1+grid_x_number                                   ! The  4th is right top    vicinity grid in Local layer
                grid_near_mesh_name(i,5) =i  -grid_x_number                                   ! The  5th is       bottom vicinity grid in Local layer
                grid_near_mesh_name(i,6) =i+1-grid_x_number                                   ! The  6th is right bottom vicinity grid in Local layer

                ! Lower layer
                grid_near_mesh_name(i,7) =i                -grid_x_number*grid_y_number       ! The  7th is itself                     in Lower layer
                grid_near_mesh_name(i,8) =i+1              -grid_x_number*grid_y_number       ! The  8th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,9) =i  +grid_x_number-grid_x_number*grid_y_number       ! The  9th is       top    vicinity grid in Lower layer
                grid_near_mesh_name(i,10)=i+1+grid_x_number-grid_x_number*grid_y_number       ! The 10th is right top    vicinity grid in Lower layer
                grid_near_mesh_name(i,11)=i  -grid_x_number-grid_x_number*grid_y_number       ! The 11th is       bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,12)=i+1-grid_x_number-grid_x_number*grid_y_number       ! The 12th is right bottom vicinity grid in Lower layer
                
            ! Right in four edges on the top face
            else if(mesh_x==grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids

                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i-1                                                 ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i-1+grid_x_number                                   ! The  3rd is left  top    vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i  +grid_x_number                                   ! The  4th is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,5) =i-1-grid_x_number                                   ! The  5th is left  bottom vicinity grid in Local layer
                grid_near_mesh_name(i,6) =i  -grid_x_number                                   ! The  6th is       bottom vicinity grid in Local layer

                ! Lower layer
                grid_near_mesh_name(i,7) =i                -grid_x_number*grid_y_number       ! The  7th is itself                     in Lower layer
                grid_near_mesh_name(i,8) =i-1              -grid_x_number*grid_y_number       ! The  8th is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,9) =i-1+grid_x_number-grid_x_number*grid_y_number       ! The  9th is left  top    vicinity grid in Lower layer
                grid_near_mesh_name(i,10)=i  +grid_x_number-grid_x_number*grid_y_number       ! The 10th is       top    vicinity grid in Lower layer
                grid_near_mesh_name(i,11)=i-1-grid_x_number-grid_x_number*grid_y_number       ! The 11th is left  bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,12)=i  -grid_x_number-grid_x_number*grid_y_number       ! The 12th is       bottom vicinity grid in Lower layer

            ! Front in four edges on the top face
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==1 .and. mesh_z==grid_z_number) then

                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids
                  
                ! Local layer 
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i-1                                                 ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i+1                                                 ! The  3rd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i-1+grid_x_number                                   ! The  4th is left  top    vicinity grid in Local layer
                grid_near_mesh_name(i,5) =i  +grid_x_number                                   ! The  5th is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,6) =i+1+grid_x_number                                   ! The  6th is right top    vicinity grid in Local layer

                ! Lower layer
                grid_near_mesh_name(i,7) =i                -grid_x_number*grid_y_number       ! The  7th is itself                     in Lower layer
                grid_near_mesh_name(i,8) =i-1              -grid_x_number*grid_y_number       ! The  8th is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,9) =i+1              -grid_x_number*grid_y_number       ! The  9th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,10)=i-1+grid_x_number-grid_x_number*grid_y_number       ! The 10th is left  top    vicinity grid in Lower layer
                grid_near_mesh_name(i,11)=i  +grid_x_number-grid_x_number*grid_y_number       ! The 11th is       top    vicinity grid in Lower layer
                grid_near_mesh_name(i,12)=i+1+grid_x_number-grid_x_number*grid_y_number       ! The 12th is right top    vicinity grid in Lower layer


            ! Back in four edges on the top face
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==grid_y_number .and. mesh_z==grid_z_number) then

                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids

                ! Local layer 
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i-1                                                 ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i+1                                                 ! The  3rd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i-1-grid_x_number                                   ! The  4th is left  bottom vicinity grid in Local layer
                grid_near_mesh_name(i,5) =i  -grid_x_number                                   ! The  5th is       bottom vicinity grid in Local layer
                grid_near_mesh_name(i,6) =i+1-grid_x_number                                   ! The  6th is right bottom vicinity grid in Local layer

                ! Lower layer
                grid_near_mesh_name(i,7) =i                -grid_x_number*grid_y_number       ! The  7th is itself                     in Lower layer
                grid_near_mesh_name(i,8) =i-1              -grid_x_number*grid_y_number       ! The  8th is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,9) =i+1              -grid_x_number*grid_y_number       ! The  9th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,10)=i-1-grid_x_number-grid_x_number*grid_y_number       ! The 10th is left  bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,11)=i  -grid_x_number-grid_x_number*grid_y_number       ! The 11th is       bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,12)=i+1-grid_x_number-grid_x_number*grid_y_number       ! The 12th is right bottom vicinity grid in Lower layer
            
            ! Four vertical edges
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Left vertical edge on front face
            else if(mesh_x==1 .and. mesh_y==1 .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i+1                                                 ! The  2nd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i  +grid_x_number                                   ! The  3rd is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i+1+grid_x_number                                   ! The  4th is right top    vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,5) =i                +grid_x_number*grid_y_number       ! The  5th is itself                     in Upper layer
                grid_near_mesh_name(i,6) =i+1              +grid_x_number*grid_y_number       ! The  6th is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,7) =i  +grid_x_number+grid_x_number*grid_y_number       ! The  7th is       top    vicinity grid in Upper layer
                grid_near_mesh_name(i,8) =i+1+grid_x_number+grid_x_number*grid_y_number       ! The  8th is right top    vicinity grid in Upper layer

                ! Lower layer
                grid_near_mesh_name(i,9) =i                -grid_x_number*grid_y_number       ! The  9th is itself                     in Lower layer
                grid_near_mesh_name(i,10)=i+1              -grid_x_number*grid_y_number       ! The 10th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,11)=i  +grid_x_number-grid_x_number*grid_y_number       ! The 11th is       top    vicinity grid in Lower layer
                grid_near_mesh_name(i,12)=i+1+grid_x_number-grid_x_number*grid_y_number       ! The 12th is right top    vicinity grid in Lower layer

            ! Left vertical edge on back face
            else if(mesh_x==1 .and. mesh_y==grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i+1                                                 ! The  2nd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i  -grid_x_number                                   ! The  3rd is       bottom vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i+1-grid_x_number                                   ! The  4th is right bottom vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,5) =i                +grid_x_number*grid_y_number       ! The  5th is itself                     in Upper layer
                grid_near_mesh_name(i,6) =i+1              +grid_x_number*grid_y_number       ! The  6th is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,7) =i  -grid_x_number+grid_x_number*grid_y_number       ! The  7th is       bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,8) =i+1-grid_x_number+grid_x_number*grid_y_number       ! The  8th is right bottom vicinity grid in Upper layer

                ! Lower layer
                grid_near_mesh_name(i,9) =i                -grid_x_number*grid_y_number       ! The  9th is itself                     in Lower layer
                grid_near_mesh_name(i,10)=i+1              -grid_x_number*grid_y_number       ! The 10th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,11)=i  -grid_x_number-grid_x_number*grid_y_number       ! The 11th is       bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,12)=i+1-grid_x_number-grid_x_number*grid_y_number       ! The 12th is right bottom vicinity grid in Lower layer

            ! Right vertical edge on front face
            else if(mesh_x==grid_x_number .and. mesh_y==1 .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i-1                                                 ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i-1+grid_x_number                                   ! The  3rd is left  top    vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i  +grid_x_number                                   ! The  4th is       top    vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,5) =i                +grid_x_number*grid_y_number       ! The  5th is itself                     in Upper layer
                grid_near_mesh_name(i,6) =i-1              +grid_x_number*grid_y_number       ! The  6th is left         vicinity grid in Upper layer
                grid_near_mesh_name(i,7) =i-1+grid_x_number+grid_x_number*grid_y_number       ! The  7th isleft   top    vicinity grid in Upper layer
                grid_near_mesh_name(i,8) =i  +grid_x_number+grid_x_number*grid_y_number       ! The  8th is       top    vicinity grid in Upper layer

                ! Lower layer
                grid_near_mesh_name(i,9) =i                -grid_x_number*grid_y_number       ! The  9th is itself                     in Lower layer
                grid_near_mesh_name(i,10)=i-1              -grid_x_number*grid_y_number       ! The 10nd is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,11)=i-1+grid_x_number-grid_x_number*grid_y_number       ! The 11th is left  top    vicinity grid in Lower layer
                grid_near_mesh_name(i,12)=i  +grid_x_number-grid_x_number*grid_y_number       ! The 12th is       top    vicinity grid in Lower layer

            ! Right vertical edge on back face
            else if(mesh_x==grid_x_number .and. mesh_y==grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=12                                                   ! Cell on dege has 12 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1) =i                                                   ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2) =i-1                                                 ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3) =i-1-grid_x_number                                   ! The  3rd is left  bottom vicinity grid in Local layer
                grid_near_mesh_name(i,4) =i  -grid_x_number                                   ! The  4th is       bottom vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,5) =i                +grid_x_number*grid_y_number       ! The  5th is itself                     in Upper layer
                grid_near_mesh_name(i,6) =i-1              +grid_x_number*grid_y_number       ! The  6th is left         vicinity grid in Upper layer
                grid_near_mesh_name(i,7) =i-1-grid_x_number+grid_x_number*grid_y_number       ! The  7th is left  bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,8) =i  -grid_x_number+grid_x_number*grid_y_number       ! The  8th is       bottom vicinity grid in Upper layer

                ! Lower layer
                grid_near_mesh_name(i,9) =i                -grid_x_number*grid_y_number       ! The  9th is itself                     in Lower layer
                grid_near_mesh_name(i,10)=i-1              -grid_x_number*grid_y_number       ! The 10th is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,11)=i-1-grid_x_number-grid_x_number*grid_y_number       ! The 11th is left  bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,12)=i  -grid_x_number-grid_x_number*grid_y_number       ! The 12th is       bottom vicinity grid in Lower layer
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Eight corners
            !--------------------------------------------------------------------------------------------------
            ! Four corners on block bottom face
            ! Left-Front corner on block bottom face  
            else if(mesh_x==1 .and. mesh_y==1 .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=8                                                    ! Cell on corner has 8 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i+1                                                  ! The  2nd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i  +grid_x_number                                    ! The  3rd is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i+1+grid_x_number                                    ! The  4th is right top    vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,5)=i                +grid_x_number*grid_y_number        ! The  5th is itself                     in Lower layer
                grid_near_mesh_name(i,6)=i+1              +grid_x_number*grid_y_number        ! The  6th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,7)=i  +grid_x_number+grid_x_number*grid_y_number        ! The  7th is       top    vicinity grid in Lower layer
                grid_near_mesh_name(i,8)=i+1+grid_x_number+grid_x_number*grid_y_number        ! The  8th is right top    vicinity grid in Lower layer

            ! Right-Front corner on block bottom face  
            else if(mesh_x==grid_x_number .and. mesh_y==1 .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=8                                                    ! Cell on corner has 8 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i-1                                                  ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i-1+grid_x_number                                    ! The  3rd is left  top    vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i  +grid_x_number                                    ! The  4th is       top    vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,5)=i                +grid_x_number*grid_y_number        ! The  5th is itself                     in Upper layer
                grid_near_mesh_name(i,6)=i-1              +grid_x_number*grid_y_number        ! The  6th is left         vicinity grid in Upper layer
                grid_near_mesh_name(i,7)=i-1+grid_x_number+grid_x_number*grid_y_number        ! The  7th is left  top    vicinity grid in Upper layer
                grid_near_mesh_name(i,8)=i  +grid_x_number+grid_x_number*grid_y_number        ! The  8th is       top    vicinity grid in Upper layer

            ! Left-Back corner on block bottom face  
            else if(mesh_x==1 .and. mesh_y==grid_y_number .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=8                                                    ! Cell on corner has 8 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i+1                                                  ! The  2nd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i  -grid_x_number                                    ! The  3rd is       bottom vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i+1-grid_x_number                                    ! The  4th is right bottom vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,5)=i                +grid_x_number*grid_y_number        ! The  5th is itself                     in Upper layer
                grid_near_mesh_name(i,6)=i+1              +grid_x_number*grid_y_number        ! The  6th is right        vicinity grid in Upper layer
                grid_near_mesh_name(i,7)=i  -grid_x_number+grid_x_number*grid_y_number        ! The  7th is       bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,8)=i+1-grid_x_number+grid_x_number*grid_y_number        ! The  8th is right bottom vicinity grid in Upper layer


            ! Right-Back corner on block bottom face  
            else if(mesh_x==grid_x_number .and. mesh_y==grid_y_number .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=8                                                    ! Cell on corner has 8 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i-1                                                  ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i-1-grid_x_number                                    ! The  3rd is left  bottom vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i  -grid_x_number                                    ! The  4th is       bottom vicinity grid in Local layer

                ! Upper layer
                grid_near_mesh_name(i,5)=i                +grid_x_number*grid_y_number        ! The  5th is itself                     in Upper layer
                grid_near_mesh_name(i,6)=i-1              +grid_x_number*grid_y_number        ! The  6th is left         vicinity grid in Upper layer
                grid_near_mesh_name(i,7)=i-1-grid_x_number+grid_x_number*grid_y_number        ! The  7th is left  bottom vicinity grid in Upper layer
                grid_near_mesh_name(i,8)=i  -grid_x_number+grid_x_number*grid_y_number        ! The  8th is       bottom vicinity grid in Upper layer

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Four corners on block top face
            ! Left-Front corner on block top face
            else if(mesh_x==1 .and. mesh_y==1 .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=8                                                    ! Cell on corner has 8 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i+1                                                  ! The  2nd is right        vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i  +grid_x_number                                    ! The  3rd is       top    vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i+1+grid_x_number                                    ! The  4th is right top    vicinity grid in Local layer

                ! Lower layer
                grid_near_mesh_name(i,5)=i                -grid_x_number*grid_y_number        ! The  5th is itself                     in Lower layer
                grid_near_mesh_name(i,6)=i+1              -grid_x_number*grid_y_number        ! The  6th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,7)=i  +grid_x_number-grid_x_number*grid_y_number        ! The  7th is       top    vicinity grid in Lower layer
                grid_near_mesh_name(i,8)=i+1+grid_x_number-grid_x_number*grid_y_number        ! The  8th is right top    vicinity grid in Lower layer

            ! Left-Back corner on block top face
            else if(mesh_x==1 .and. mesh_y==grid_y_number .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=8                                                    ! Cell on corner has 8 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i+1                                                  ! The  2nd is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,3)=i  -grid_x_number                                    ! The  3rd is       bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,4)=i+1-grid_x_number                                    ! The  4th is right bottom vicinity grid in Lower layer

                ! Lower layer
                grid_near_mesh_name(i,5)=i                -grid_x_number*grid_y_number        ! The  5th is itself                     in Local layer
                grid_near_mesh_name(i,6)=i+1              -grid_x_number*grid_y_number        ! The  6th is right        vicinity grid in Lower layer
                grid_near_mesh_name(i,7)=i  -grid_x_number-grid_x_number*grid_y_number        ! The  7th is       bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,8)=i+1-grid_x_number-grid_x_number*grid_y_number        ! The  8th is right bottom vicinity grid in Lower layer

            ! Right-Front corner on block top face
            else if(mesh_x==grid_x_number .and. mesh_y==1 .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=8                                                    ! Cell on corner has 8 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i-1                                                  ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i-1+grid_x_number                                    ! The  3rd is left  top    vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i  +grid_x_number                                    ! The  4th is       top    vicinity grid in Local layer

                ! Lower layer
                grid_near_mesh_name(i,5)=i                -grid_x_number*grid_y_number        ! The  5th is itself                     in Lower layer
                grid_near_mesh_name(i,6)=i-1              -grid_x_number*grid_y_number        ! The  6th is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,7)=i-1+grid_x_number-grid_x_number*grid_y_number        ! The  7th is left  top    vicinity grid in Lower layer
                grid_near_mesh_name(i,8)=i  +grid_x_number-grid_x_number*grid_y_number        ! The  8th is       top    vicinity grid in Lower layer

            ! Right-Back corner on block top face
            else if(mesh_x==grid_x_number .and. mesh_y==grid_y_number .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=8                                                    ! Cell on corner has 8 vicinity grids
                
                ! Local layer
                grid_near_mesh_name(i,1)=i                                                    ! The  1st is itself                     in Local layer
                grid_near_mesh_name(i,2)=i-1                                                  ! The  2nd is left         vicinity grid in Local layer
                grid_near_mesh_name(i,3)=i-1-grid_x_number                                    ! The  3rd is left  bottom vicinity grid in Local layer
                grid_near_mesh_name(i,4)=i  -grid_x_number                                    ! The  4th is       bottom vicinity grid in Local layer

                ! Lower layer
                grid_near_mesh_name(i,5)=i                -grid_x_number*grid_y_number        ! The  5th is itself                     in Lower layer
                grid_near_mesh_name(i,6)=i-1              -grid_x_number*grid_y_number        ! The  6th is left         vicinity grid in Lower layer
                grid_near_mesh_name(i,7)=i-1-grid_x_number-grid_x_number*grid_y_number        ! The  7th is left  bottom vicinity grid in Lower layer
                grid_near_mesh_name(i,8)=i  -grid_x_number-grid_x_number*grid_y_number        ! The  8th is       bottom vicinity grid in Lower layer

            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            end if  
            
        end do

     else 

        write(*,'(A)') " The Dimension input in 'Distribute_Backgroud_Grid.f90' is not right! "

     endif Dimension_If
    !==========================================================================================================







!     !**********************************************************************************************************
!     !For checking
!     !Call MPI functions
!     call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
!     call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    
!     if (Current_Processor_ID==0) then
       
!        open(unit=5,file="Backgroud_Grid.dat",status="replace",position="rewind",action="write",iostat=ioerror)

!         !Output the vicinity gird name and number
!         do i=1,grid_max_number
           
!            !write(5,100) i,(grid_center_positon(i,j),j=1,dim),grid_near_mesh_number(i),(grid_near_mesh_name(i,j),j=1,grid_near_mesh_number(i))
!            write(5,100) i,grid_near_mesh_number(i),(grid_near_mesh_name(i,j),j=1,grid_near_mesh_number(i))
! 100        format(100I10) 
          
!            ! if(i==165) then
!            !     write(*,*) i,chain_near_mesh_number(i)
!            !     k=chain_near_mesh_number(i)
           
!            !     write(*,*) (chain_near_mesh_name(i,j),j=1,k)
!            ! end if

!         end do

!         close(5)

!     end if
!     !**********************************************************************************************************

end subroutine Distribute_Backgroud_Grid