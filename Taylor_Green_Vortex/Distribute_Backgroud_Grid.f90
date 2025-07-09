!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
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
!****************************************************************************
subroutine Distribute_Backgroud_Grid(dim,&                           !Dimension of grid
                                     grid_x_number,&                 !Grid number in x direction
                                     grid_y_number,&                 !Grid number in y direction
                                     grid_z_number,&                 !Grid number in z direction
                                     grid_max_number,&               !Max grid number (grid_x_number*grid_y_number*grid_z_number)
                                     grid_center_positon,&           !Grid center positon
                                     grid_near_mesh_number,&         !Near mesh number
                                     grid_near_mesh_name&            !grid near mesh name
                                     )
    use MPI

    implicit none

    !==========================================================================================================
    
    !----------------------------------------------------------------------------------------------------------
    !Variables from mother subroutine
    integer,intent(in)::dim                                                      !Dimension of grid
    integer,intent(in)::grid_x_number,grid_y_number,grid_z_number                !Grid number in x,y,z direction
    integer,intent(in)::grid_max_number                                          !MAX GRID NUMBER (grid_x_number*grid_y_number*grid_z_number)
    integer,dimension(grid_max_number,3),intent(out)::grid_center_positon        !Near mesh number
    integer,dimension(grid_max_number),intent(out)::grid_near_mesh_number        !Near mesh number
    integer,dimension(grid_max_number,27),intent(out)::grid_near_mesh_name       !Near mesh name
    !----------------------------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------------------------
    !Variables in subroutine
    integer::i,j,k,L,m                                                  !loop Variables
    
    !Read file varibales
    integer::ioerror=0                                                  !open file return value
    integer::stat                                                       !read file data return value
    integer::status                                                     !allocate memory status

    integer::mesh_x,mesh_y,mesh_z                                       !temp varibales for backgroud mesh

    ! ! Variables of the MPI Processor for checking
    ! integer::ierror_MPI                                                                         ! MPI Function return value for success runing or not
    ! integer,dimension(MPI_STATUS_SIZE)::Status_MPI                                              ! MPI running status
    ! integer::Current_Processor_ID                                                               ! Current Processor ID
    ! integer::Total_Processors_Number                                                            ! Total Processors Number (0,1,2,...,Total_Processors_Number-1)

    !----------------------------------------------------------------------------------------------------------

    !==========================================================================================================
    !For 2-D case
    if (dim==2) then

        !*************************************开始划分网格(网格中心点的计算坐标是对的)******************************
        !write(*,*) chain_dx,chain_dy,chain_x_number,chain_y_number
        !网格是一行一行的数的
        do j=1,grid_y_number
            do i=1,grid_x_number
                k=i+(j-1)*grid_x_number
                grid_center_positon(k,1)=i
                grid_center_positon(k,2)=j 
            end do
        end do
        
        !查找相邻网格编号及相邻网格数量(如果边界除以网格间距是整数，那么这是正确的)
        grid_near_mesh_number=0                    !初始化网格相邻网格数量数组
        grid_near_mesh_name=0                      !初始化网格相邻网格名称数组
        
        do i=1,grid_max_number
            
            !取出网格中心点坐标
            mesh_x=grid_center_positon(i,1)
            mesh_y=grid_center_positon(i,2)
            
            !if(i==163) then
            !    write(*,*) mesh_x,mesh_y
            !end if  
            
            !处理中间网格
            if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number ) then
                
                grid_near_mesh_number(i)=9                         !中间的网格有九个网格
                
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                     !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number      !第4个网格为左上侧网格
                grid_near_mesh_name(i,5)=i+grid_x_number        !第5个网格为上侧网格
                grid_near_mesh_name(i,6)=i+1+grid_x_number      !第6个网格为右上侧网格
                grid_near_mesh_name(i,7)=i-1-grid_x_number      !第7个网格为左下侧网格
                grid_near_mesh_name(i,8)=i-grid_x_number        !第8个网格为下侧网格
                grid_near_mesh_name(i,9)=i+1-grid_x_number      !第9个网格为右下侧网格
                
            !左侧边界网格(不包括上下部分)
            else if(mesh_x==1 .and. mesh_y>1 .and. mesh_y<grid_y_number) then
                    
                grid_near_mesh_number(i)=6                         !中间的网格有九个网格
            
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                     !第2个网格为右侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number        !第3个网格为下侧网格
                grid_near_mesh_name(i,4)=i+1+grid_x_number      !第4个网格为右下侧网格
                grid_near_mesh_name(i,5)=i-grid_x_number        !第5个网格为上侧网格
                grid_near_mesh_name(i,6)=i+1-grid_x_number      !第6个网格为右上侧网格
 
            !右侧边界网格(不包括上下部分)
            else if(mesh_x==grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number) then
                    
                grid_near_mesh_number(i)=6                         !中间的网格有九个网格
            
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number        !第3个网格为上侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number      !第4个网格为左上侧网格
                grid_near_mesh_name(i,5)=i-grid_x_number        !第5个网格为下侧网格
                grid_near_mesh_name(i,6)=i-1-grid_x_number      !第6个网格为左下侧网格
                
            !上部边界网格
            else if(mesh_y==grid_y_number.and. mesh_x>1 .and. mesh_x<grid_x_number) then
                    
                grid_near_mesh_number(i)=6                         !中间的网格有九个网格
            
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                     !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1-grid_x_number      !第4个网格为左下侧网格
                grid_near_mesh_name(i,5)=i-grid_x_number        !第5个网格为下侧网格
                grid_near_mesh_name(i,6)=i+1-grid_x_number      !第6个网格为右下侧网格
                 
            !下部边界网格
            else if(mesh_y==1.and. mesh_x>1 .and. mesh_x<grid_x_number) then
                    
                grid_near_mesh_number(i)=6                         !中间的网格有九个网格
            
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                     !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number      !第4个网格为左上侧网格
                grid_near_mesh_name(i,5)=i+grid_x_number        !第5个网格为上侧网格
                grid_near_mesh_name(i,6)=i+1+grid_x_number      !第6个网格为右上侧网格
                
            !四个顶角
            else if(mesh_x==1 .and. mesh_y==grid_y_number) then    !左上
                    
                grid_near_mesh_number(i)=4                         !中间的网格有九个网格
            
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                     !第2个网格为右侧网格
                grid_near_mesh_name(i,3)=i-grid_x_number        !第3个网格为下侧网格
                grid_near_mesh_name(i,4)=i+1-grid_x_number      !第4个网格为右上侧网格
                
            else if(mesh_x==grid_x_number .and. mesh_y==grid_y_number) then    !右上
                    
                grid_near_mesh_number(i)=4                         !中间的网格有九个网格
            
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i-grid_x_number        !第3个网格为下侧网格
                grid_near_mesh_name(i,4)=i-1-grid_x_number      !第4个网格为左下侧网格
                
            else if(mesh_x==1 .and. mesh_y==1) then    !左下
                    
                grid_near_mesh_number(i)=4                         !中间的网格有九个网格
            
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                     !第2个网格为右侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number        !第3个网格为上侧网格
                grid_near_mesh_name(i,4)=i+1+grid_x_number      !第4个网格为右上侧网格

                
            else if(mesh_x==grid_x_number .and. mesh_y==1) then    !右下
                    
                grid_near_mesh_number(i)=4                         !中间的网格有九个网格
            
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number        !第3个网格为上侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number      !第4个网格为左上侧网格
                
            end if  
            
        end do 
        !***********************************************************************************************************

    !For 3-D case
    else if (dim==3) then
        
        !*************************************开始划分网格(网格中心点的计算坐标是对的)******************************
        !write(*,*) chain_dx,chain_dy,chain_x_number,chain_y_number
        !网格是一行一行的数的
        do k=1,grid_z_number
            do j=1,grid_y_number
                do i=1,grid_x_number
                    
                    L=i+(j-1)*grid_x_number+(k-1)*grid_x_number*grid_y_number
                    grid_center_positon(L,1)=i
                    grid_center_positon(L,2)=j 
                    grid_center_positon(L,3)=k 

                end do
            end do
        end do
        
        !查找相邻网格编号及相邻网格数量(如果边界除以网格间距是整数，那么这是正确的)
        grid_near_mesh_number=0                    !初始化网格相邻网格数量数组
        grid_near_mesh_name=0                      !初始化网格相邻网格名称数组
        
        do i=1,grid_max_number
            
            !取出网格中心点坐标
            mesh_x=grid_center_positon(i,1)
            mesh_y=grid_center_positon(i,2)
            mesh_z=grid_center_positon(i,3)
            
            !if(i==163) then
            !    write(*,*) mesh_x,mesh_y
            !end if  
            
            !网格体内部
            if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then
                
                grid_near_mesh_number(i)=27                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number                                   !第4个网格为左上侧网格
                grid_near_mesh_name(i,5)=i+grid_x_number                                     !第5个网格为上侧网格
                grid_near_mesh_name(i,6)=i+1+grid_x_number                                   !第6个网格为右上侧网格
                grid_near_mesh_name(i,7)=i-1-grid_x_number                                   !第7个网格为左下侧网格
                grid_near_mesh_name(i,8)=i-grid_x_number                                     !第8个网格为下侧网格
                grid_near_mesh_name(i,9)=i+1-grid_x_number                                   !第9个网格为右下侧网格

                !上层网格
                grid_near_mesh_name(i,10)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,11)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,12)=i+1+grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,13)=i-1+grid_x_number+grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,14)=i+grid_x_number+grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,15)=i+1+grid_x_number+grid_x_number*grid_y_number       !第6个网格为右上侧网格
                grid_near_mesh_name(i,16)=i-1-grid_x_number+grid_x_number*grid_y_number       !第7个网格为左下侧网格
                grid_near_mesh_name(i,17)=i-grid_x_number+grid_x_number*grid_y_number         !第8个网格为下侧网格
                grid_near_mesh_name(i,18)=i+1-grid_x_number+grid_x_number*grid_y_number       !第9个网格为右下侧网格

                !下层网格
                grid_near_mesh_name(i,19)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,20)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,21)=i+1-grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,22)=i-1+grid_x_number-grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,23)=i+grid_x_number-grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,24)=i+1+grid_x_number-grid_x_number*grid_y_number       !第6个网格为右上侧网格
                grid_near_mesh_name(i,25)=i-1-grid_x_number-grid_x_number*grid_y_number       !第7个网格为左下侧网格
                grid_near_mesh_name(i,26)=i-grid_x_number-grid_x_number*grid_y_number         !第8个网格为下侧网格
                grid_near_mesh_name(i,27)=i+1-grid_x_number-grid_x_number*grid_y_number       !第9个网格为右下侧网格
                
            !网格六个面(棱和角除外)

            !-------------------------------------------------------------------------------------------------------
            !左侧面(棱和角除外)
            else if(mesh_x==1 .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then
                
                grid_near_mesh_number(i)=18                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                     !第2个网格为右侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number         !第3个网格为下侧网格
                grid_near_mesh_name(i,4)=i+1+grid_x_number       !第4个网格为右下侧网格
                grid_near_mesh_name(i,5)=i-grid_x_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,6)=i+1-grid_x_number       !第6个网格为右上侧网格

                !上层网格
                grid_near_mesh_name(i,7)=i+grid_x_number*grid_y_number                        !第1个网格为自身
                grid_near_mesh_name(i,8)=i+1+grid_x_number*grid_y_number                      !第2个网格为右侧网格
                grid_near_mesh_name(i,9)=i+grid_x_number+grid_x_number*grid_y_number          !第3个网格为下侧网格
                grid_near_mesh_name(i,10)=i+1+grid_x_number+grid_x_number*grid_y_number       !第4个网格为右下侧网格
                grid_near_mesh_name(i,11)=i-grid_x_number+grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,12)=i+1-grid_x_number+grid_x_number*grid_y_number       !第6个网格为右上侧网格

                !下层网格
                grid_near_mesh_name(i,13)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,14)=i+1-grid_x_number*grid_y_number                     !第2个网格为右侧网格
                grid_near_mesh_name(i,15)=i+grid_x_number-grid_x_number*grid_y_number         !第3个网格为下侧网格
                grid_near_mesh_name(i,16)=i+1+grid_x_number-grid_x_number*grid_y_number       !第4个网格为右下侧网格
                grid_near_mesh_name(i,17)=i-grid_x_number-grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,18)=i+1-grid_x_number-grid_x_number*grid_y_number       !第6个网格为右上侧网格
            
            !-------------------------------------------------------------------------------------------------------
            !右侧面(棱和角除外)
            else if(mesh_x==grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=18                      !中间的网格有九个网格

                !本层网格
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number         !第3个网格为上侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,5)=i-grid_x_number         !第5个网格为下侧网格
                grid_near_mesh_name(i,6)=i-1-grid_x_number       !第6个网格为左下侧网格

                !上层网格
                grid_near_mesh_name(i,7)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,8)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,9)=i+grid_x_number+grid_x_number*grid_y_number         !第3个网格为上侧网格
                grid_near_mesh_name(i,10)=i-1+grid_x_number+grid_x_number*grid_y_number      !第4个网格为左上侧网格
                grid_near_mesh_name(i,11)=i-grid_x_number+grid_x_number*grid_y_number        !第5个网格为下侧网格
                grid_near_mesh_name(i,12)=i-1-grid_x_number+grid_x_number*grid_y_number      !第6个网格为左下侧网格

                !下层网格
                grid_near_mesh_name(i,13)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,14)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,15)=i+grid_x_number-grid_x_number*grid_y_number         !第3个网格为上侧网格
                grid_near_mesh_name(i,16)=i-1+grid_x_number-grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,17)=i-grid_x_number-grid_x_number*grid_y_number         !第5个网格为下侧网格
                grid_near_mesh_name(i,18)=i-1-grid_x_number-grid_x_number*grid_y_number       !第6个网格为左下侧网格

            !-------------------------------------------------------------------------------------------------------
            !前侧面(棱和角除外)
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==1 .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=18                         !中间的网格有九个网格
            
                !本层网格
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                     !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number      !第4个网格为左上侧网格
                grid_near_mesh_name(i,5)=i+grid_x_number        !第5个网格为上侧网格
                grid_near_mesh_name(i,6)=i+1+grid_x_number      !第6个网格为右上侧网格

                !上层网格
                grid_near_mesh_name(i,7)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,8)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,9)=i+1+grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,10)=i-1+grid_x_number+grid_x_number*grid_y_number      !第4个网格为左上侧网格
                grid_near_mesh_name(i,11)=i+grid_x_number+grid_x_number*grid_y_number        !第5个网格为上侧网格
                grid_near_mesh_name(i,12)=i+1+grid_x_number+grid_x_number*grid_y_number      !第6个网格为右上侧网格

                !下层网格
                grid_near_mesh_name(i,13)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,14)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,15)=i+1-grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,16)=i-1+grid_x_number-grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,17)=i+grid_x_number-grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,18)=i+1+grid_x_number-grid_x_number*grid_y_number       !第6个网格为右上侧网格


            !-------------------------------------------------------------------------------------------------------
            !后侧面(棱和角除外)
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then
                
                grid_near_mesh_number(i)=18                         !中间的网格有九个网格
            
                !本层网格
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                     !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1-grid_x_number       !第4个网格为左下侧网格
                grid_near_mesh_name(i,5)=i-grid_x_number         !第5个网格为下侧网格
                grid_near_mesh_name(i,6)=i+1-grid_x_number       !第6个网格为右下侧网格

                !上层网格
                grid_near_mesh_name(i,7)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,8)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,9)=i+1+grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,10)=i-1-grid_x_number+grid_x_number*grid_y_number       !第4个网格为左下侧网格
                grid_near_mesh_name(i,11)=i-grid_x_number+grid_x_number*grid_y_number         !第5个网格为下侧网格
                grid_near_mesh_name(i,12)=i+1-grid_x_number+grid_x_number*grid_y_number       !第6个网格为右下侧网格

                !下层网格
                grid_near_mesh_name(i,13)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,14)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,15)=i+1-grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,16)=i-1-grid_x_number-grid_x_number*grid_y_number       !第4个网格为左下侧网格
                grid_near_mesh_name(i,17)=i-grid_x_number-grid_x_number*grid_y_number         !第5个网格为下侧网格
                grid_near_mesh_name(i,18)=i+1-grid_x_number-grid_x_number*grid_y_number       !第6个网格为右下侧网格


            !-------------------------------------------------------------------------------------------------------
            !顶侧面(棱和角除外)
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==grid_z_number) then

                grid_near_mesh_number(i)=18                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number                                   !第4个网格为左上侧网格
                grid_near_mesh_name(i,5)=i+grid_x_number                                     !第5个网格为上侧网格
                grid_near_mesh_name(i,6)=i+1+grid_x_number                                   !第6个网格为右上侧网格
                grid_near_mesh_name(i,7)=i-1-grid_x_number                                   !第7个网格为左下侧网格
                grid_near_mesh_name(i,8)=i-grid_x_number                                     !第8个网格为下侧网格
                grid_near_mesh_name(i,9)=i+1-grid_x_number                                   !第9个网格为右下侧网格

                !下层网格
                grid_near_mesh_name(i,10)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,11)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,12)=i+1-grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,13)=i-1+grid_x_number-grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,14)=i+grid_x_number-grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,15)=i+1+grid_x_number-grid_x_number*grid_y_number       !第6个网格为右上侧网格
                grid_near_mesh_name(i,16)=i-1-grid_x_number-grid_x_number*grid_y_number       !第7个网格为左下侧网格
                grid_near_mesh_name(i,17)=i-grid_x_number-grid_x_number*grid_y_number         !第8个网格为下侧网格
                grid_near_mesh_name(i,18)=i+1-grid_x_number-grid_x_number*grid_y_number       !第9个网格为右下侧网格

            !-------------------------------------------------------------------------------------------------------
            !底侧面(棱和角除外)
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==1) then

                grid_near_mesh_number(i)=18                      !中间的网格有九个网格

                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number                                   !第4个网格为左上侧网格
                grid_near_mesh_name(i,5)=i+grid_x_number                                     !第5个网格为上侧网格
                grid_near_mesh_name(i,6)=i+1+grid_x_number                                   !第6个网格为右上侧网格
                grid_near_mesh_name(i,7)=i-1-grid_x_number                                   !第7个网格为左下侧网格
                grid_near_mesh_name(i,8)=i-grid_x_number                                     !第8个网格为下侧网格
                grid_near_mesh_name(i,9)=i+1-grid_x_number                                   !第9个网格为右下侧网格

                !上层网格
                grid_near_mesh_name(i,10)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,11)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,12)=i+1+grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,13)=i-1+grid_x_number+grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,14)=i+grid_x_number+grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,15)=i+1+grid_x_number+grid_x_number*grid_y_number       !第6个网格为右上侧网格
                grid_near_mesh_name(i,16)=i-1-grid_x_number+grid_x_number*grid_y_number       !第7个网格为左下侧网格
                grid_near_mesh_name(i,17)=i-grid_x_number+grid_x_number*grid_y_number         !第8个网格为下侧网格
                grid_near_mesh_name(i,18)=i+1-grid_x_number+grid_x_number*grid_y_number       !第9个网格为右下侧网格

            !-------------------------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------------------------
            !十二条棱
            !底部四条棱。
            !=======================================================================================================
            !底部左棱
            else if(mesh_x==1 .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=12                      !中间的网格有九个网格

                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number                                     !第5个网格为上侧网格
                grid_near_mesh_name(i,4)=i+1+grid_x_number                                   !第6个网格为右上侧网格
                grid_near_mesh_name(i,5)=i-grid_x_number                                     !第8个网格为下侧网格
                grid_near_mesh_name(i,6)=i+1-grid_x_number                                   !第9个网格为右下侧网格

                !上层网格
                grid_near_mesh_name(i,7)=i+grid_x_number*grid_y_number                        !第1个网格为自身
                grid_near_mesh_name(i,8)=i+1+grid_x_number*grid_y_number                      !第3个网格为右侧网格
                grid_near_mesh_name(i,9)=i+grid_x_number+grid_x_number*grid_y_number          !第5个网格为上侧网格
                grid_near_mesh_name(i,10)=i+1+grid_x_number+grid_x_number*grid_y_number        !第6个网格为右上侧网格
                grid_near_mesh_name(i,11)=i-grid_x_number+grid_x_number*grid_y_number          !第8个网格为下侧网格
                grid_near_mesh_name(i,12)=i+1-grid_x_number+grid_x_number*grid_y_number        !第9个网格为右下侧网格
                
            !底部右棱
            else if(mesh_x==grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=12                      !中间的网格有九个网格

                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i-1+grid_x_number                                   !第4个网格为左上侧网格
                grid_near_mesh_name(i,4)=i+grid_x_number                                     !第5个网格为上侧网格
                grid_near_mesh_name(i,5)=i-1-grid_x_number                                   !第7个网格为左下侧网格
                grid_near_mesh_name(i,6)=i-grid_x_number                                     !第8个网格为下侧网格

                !上层网格
                grid_near_mesh_name(i,7)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,8)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,9)=i-1+grid_x_number+grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,10)=i+grid_x_number+grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,11)=i-1-grid_x_number+grid_x_number*grid_y_number       !第7个网格为左下侧网格
                grid_near_mesh_name(i,12)=i-grid_x_number+grid_x_number*grid_y_number         !第8个网格为下侧网格

            !底部前棱
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==1 .and. mesh_z==1) then

                grid_near_mesh_number(i)=12                      !中间的网格有九个网格
                  
                !本层网格 
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                     !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,5)=i+grid_x_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,6)=i+1+grid_x_number       !第6个网格为右上侧网格

                !上层网格
                grid_near_mesh_name(i,7)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,8)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,9)=i+1+grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,10)=i-1+grid_x_number+grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,11)=i+grid_x_number+grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,12)=i+1+grid_x_number+grid_x_number*grid_y_number       !第6个网格为右上侧网格


            !底部后棱
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==grid_y_number .and. mesh_z==1) then

                grid_near_mesh_number(i)=12                      !中间的网格有九个网格

                !本层网格 
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                     !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1-grid_x_number       !第4个网格为左下侧网格
                grid_near_mesh_name(i,5)=i-grid_x_number         !第5个网格为下侧网格
                grid_near_mesh_name(i,6)=i+1-grid_x_number       !第6个网格为右下侧网格

                !上层网格
                grid_near_mesh_name(i,7)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,8)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,9)=i+1+grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,10)=i-1-grid_x_number+grid_x_number*grid_y_number       !第4个网格为左下侧网格
                grid_near_mesh_name(i,11)=i-grid_x_number+grid_x_number*grid_y_number         !第5个网格为下侧网格
                grid_near_mesh_name(i,12)=i+1-grid_x_number+grid_x_number*grid_y_number       !第6个网格为右下侧网格

            !=======================================================================================================

            !顶部四条棱。
            !=======================================================================================================
            !顶部左棱
            else if(mesh_x==1 .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=12                      !中间的网格有九个网格

                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number                                     !第5个网格为上侧网格
                grid_near_mesh_name(i,4)=i+1+grid_x_number                                   !第6个网格为右上侧网格
                grid_near_mesh_name(i,5)=i-grid_x_number                                     !第8个网格为下侧网格
                grid_near_mesh_name(i,6)=i+1-grid_x_number                                   !第9个网格为右下侧网格

                !下层网格
                grid_near_mesh_name(i,7)=i-grid_x_number*grid_y_number                        !第1个网格为自身
                grid_near_mesh_name(i,8)=i+1-grid_x_number*grid_y_number                      !第3个网格为右侧网格
                grid_near_mesh_name(i,9)=i+grid_x_number-grid_x_number*grid_y_number          !第5个网格为上侧网格
                grid_near_mesh_name(i,10)=i+1+grid_x_number-grid_x_number*grid_y_number        !第6个网格为右上侧网格
                grid_near_mesh_name(i,11)=i-grid_x_number-grid_x_number*grid_y_number          !第8个网格为下侧网格
                grid_near_mesh_name(i,12)=i+1-grid_x_number-grid_x_number*grid_y_number        !第9个网格为右下侧网格
                
            !顶部右棱
            else if(mesh_x==grid_x_number .and. mesh_y>1 .and. mesh_y<grid_y_number .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=12                      !中间的网格有九个网格

                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i-1+grid_x_number                                   !第4个网格为左上侧网格
                grid_near_mesh_name(i,4)=i+grid_x_number                                     !第5个网格为上侧网格
                grid_near_mesh_name(i,5)=i-1-grid_x_number                                   !第7个网格为左下侧网格
                grid_near_mesh_name(i,6)=i-grid_x_number                                     !第8个网格为下侧网格

                !下层网格
                grid_near_mesh_name(i,7)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,8)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,9)=i-1+grid_x_number-grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,10)=i+grid_x_number-grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,11)=i-1-grid_x_number-grid_x_number*grid_y_number       !第7个网格为左下侧网格
                grid_near_mesh_name(i,12)=i-grid_x_number-grid_x_number*grid_y_number         !第8个网格为下侧网格

            !顶部前棱
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==1 .and. mesh_z==grid_z_number) then

                grid_near_mesh_number(i)=12                      !中间的网格有九个网格
                  
                !本层网格 
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                     !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1+grid_x_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,5)=i+grid_x_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,6)=i+1+grid_x_number       !第6个网格为右上侧网格

                !下层网格
                grid_near_mesh_name(i,7)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,8)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,9)=i+1-grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,10)=i-1+grid_x_number-grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,11)=i+grid_x_number-grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,12)=i+1+grid_x_number-grid_x_number*grid_y_number       !第6个网格为右上侧网格


            !顶部后棱
            else if(mesh_x>1 .and. mesh_x<grid_x_number .and. mesh_y==grid_y_number .and. mesh_z==grid_z_number) then

                grid_near_mesh_number(i)=12                      !中间的网格有九个网格

                !本层网格 
                grid_near_mesh_name(i,1)=i                       !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                     !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i+1                     !第3个网格为右侧网格
                grid_near_mesh_name(i,4)=i-1-grid_x_number       !第4个网格为左下侧网格
                grid_near_mesh_name(i,5)=i-grid_x_number         !第5个网格为下侧网格
                grid_near_mesh_name(i,6)=i+1-grid_x_number       !第6个网格为右下侧网格

                !下层网格
                grid_near_mesh_name(i,7)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,8)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,9)=i+1-grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,10)=i-1-grid_x_number-grid_x_number*grid_y_number       !第4个网格为左下侧网格
                grid_near_mesh_name(i,11)=i-grid_x_number-grid_x_number*grid_y_number         !第5个网格为下侧网格
                grid_near_mesh_name(i,12)=i+1-grid_x_number-grid_x_number*grid_y_number       !第6个网格为右下侧网格
            
            !垂直四条棱。
            !=======================================================================================================
            !垂直左前棱
            else if(mesh_x==1 .and. mesh_y==1 .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=12                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number                                     !第5个网格为上侧网格
                grid_near_mesh_name(i,4)=i+1+grid_x_number                                   !第6个网格为右上侧网格

                !上层网格
                grid_near_mesh_name(i,5)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i+1+grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,7)=i+grid_x_number+grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,8)=i+1+grid_x_number+grid_x_number*grid_y_number       !第6个网格为右上侧网格

                !下层网格
                grid_near_mesh_name(i,9)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,10)=i+1-grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,11)=i+grid_x_number-grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,12)=i+1+grid_x_number-grid_x_number*grid_y_number       !第6个网格为右上侧网格

            !垂直左后棱
            else if(mesh_x==1 .and. mesh_y==grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=12                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,3)=i-grid_x_number                                     !第8个网格为下侧网格
                grid_near_mesh_name(i,4)=i+1-grid_x_number                                   !第9个网格为右下侧网格

                !上层网格
                grid_near_mesh_name(i,5)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i+1+grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,7)=i-grid_x_number+grid_x_number*grid_y_number         !第8个网格为下侧网格
                grid_near_mesh_name(i,8)=i+1-grid_x_number+grid_x_number*grid_y_number       !第9个网格为右下侧网格

                !下层网格
                grid_near_mesh_name(i,9)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,10)=i+1-grid_x_number*grid_y_number                    !第3个网格为右侧网格
                grid_near_mesh_name(i,11)=i-grid_x_number-grid_x_number*grid_y_number        !第8个网格为下侧网格
                grid_near_mesh_name(i,12)=i+1-grid_x_number-grid_x_number*grid_y_number      !第9个网格为右下侧网格

            !垂直右前棱
            else if(mesh_x==grid_x_number .and. mesh_y==1 .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=12                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i-1+grid_x_number                                   !第4个网格为左上侧网格
                grid_near_mesh_name(i,4)=i+grid_x_number                                     !第5个网格为上侧网格

                !上层网格
                grid_near_mesh_name(i,5)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,7)=i-1+grid_x_number+grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,8)=i+grid_x_number+grid_x_number*grid_y_number         !第5个网格为上侧网格

                !下层网格
                grid_near_mesh_name(i,9)=i-grid_x_number*grid_y_number                        !第1个网格为自身
                grid_near_mesh_name(i,10)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,11)=i-1+grid_x_number-grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,12)=i+grid_x_number-grid_x_number*grid_y_number         !第5个网格为上侧网格

            !垂直右后棱
            else if(mesh_x==grid_x_number .and. mesh_y==grid_y_number .and. mesh_z>1 .and. mesh_z<grid_z_number) then

                grid_near_mesh_number(i)=12                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i-1-grid_x_number                                   !第7个网格为左下侧网格
                grid_near_mesh_name(i,4)=i-grid_x_number                                     !第8个网格为下侧网格

                !上层网格
                grid_near_mesh_name(i,5)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,7)=i-1-grid_x_number+grid_x_number*grid_y_number       !第7个网格为左下侧网格
                grid_near_mesh_name(i,8)=i-grid_x_number+grid_x_number*grid_y_number         !第8个网格为下侧网格

                !下层网格
                grid_near_mesh_name(i,9)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,10)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,11)=i-1-grid_x_number-grid_x_number*grid_y_number       !第7个网格为左下侧网格
                grid_near_mesh_name(i,12)=i-grid_x_number-grid_x_number*grid_y_number         !第8个网格为下侧网格

            !====================================================================================================

            !----------------------------------------------------------------------------------------------------
            !八个角点
            !====================================================================================================
            !底部平面四个角点
            !底部左前
            else if(mesh_x==1 .and. mesh_y==1 .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=8                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number                                     !第5个网格为上侧网格
                grid_near_mesh_name(i,4)=i+1+grid_x_number                                   !第6个网格为右上侧网格

                !上层网格
                grid_near_mesh_name(i,5)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i+1+grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,7)=i+grid_x_number+grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,8)=i+1+grid_x_number+grid_x_number*grid_y_number       !第6个网格为右上侧网格

            !底部右前
            else if(mesh_x==grid_x_number .and. mesh_y==1 .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=8                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i-1+grid_x_number                                   !第4个网格为左上侧网格
                grid_near_mesh_name(i,4)=i+grid_x_number                                     !第5个网格为上侧网格

                !上层网格
                grid_near_mesh_name(i,5)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,7)=i-1+grid_x_number+grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,8)=i+grid_x_number+grid_x_number*grid_y_number         !第5个网格为上侧网格

            !底部左后
            else if(mesh_x==1 .and. mesh_y==grid_y_number .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=8                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,3)=i-grid_x_number                                     !第8个网格为下侧网格
                grid_near_mesh_name(i,4)=i+1-grid_x_number                                   !第9个网格为右下侧网格

                !上层网格
                grid_near_mesh_name(i,5)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i+1+grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,7)=i-grid_x_number+grid_x_number*grid_y_number         !第8个网格为下侧网格
                grid_near_mesh_name(i,8)=i+1-grid_x_number+grid_x_number*grid_y_number       !第9个网格为右下侧网格


            !底部右后
            else if(mesh_x==grid_x_number .and. mesh_y==grid_y_number .and. mesh_z==1) then
                
                grid_near_mesh_number(i)=8                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i-1-grid_x_number                                   !第7个网格为左下侧网格
                grid_near_mesh_name(i,4)=i-grid_x_number                                     !第8个网格为下侧网格

                !上层网格
                grid_near_mesh_name(i,5)=i+grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i-1+grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,7)=i-1-grid_x_number+grid_x_number*grid_y_number       !第7个网格为左下侧网格
                grid_near_mesh_name(i,8)=i-grid_x_number+grid_x_number*grid_y_number         !第8个网格为下侧网格


            !====================================================================================================
            !顶部平面四个角点
            !顶部左前
            else if(mesh_x==1 .and. mesh_y==1 .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=8                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,3)=i+grid_x_number                                     !第5个网格为上侧网格
                grid_near_mesh_name(i,4)=i+1+grid_x_number                                   !第6个网格为右上侧网格


                !下层网格
                grid_near_mesh_name(i,5)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i+1-grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,7)=i+grid_x_number-grid_x_number*grid_y_number         !第5个网格为上侧网格
                grid_near_mesh_name(i,8)=i+1+grid_x_number-grid_x_number*grid_y_number       !第6个网格为右上侧网格


            !顶部左后
            else if(mesh_x==1 .and. mesh_y==grid_y_number .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=8                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i+1                                                 !第3个网格为右侧网格
                grid_near_mesh_name(i,3)=i-grid_x_number                                     !第8个网格为下侧网格
                grid_near_mesh_name(i,4)=i+1-grid_x_number                                   !第9个网格为右下侧网格


                !下层网格
                grid_near_mesh_name(i,5)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i+1-grid_x_number*grid_y_number                     !第3个网格为右侧网格
                grid_near_mesh_name(i,7)=i-grid_x_number-grid_x_number*grid_y_number         !第8个网格为下侧网格
                grid_near_mesh_name(i,8)=i+1-grid_x_number-grid_x_number*grid_y_number       !第9个网格为右下侧网格

            !顶部右前
            else if(mesh_x==grid_x_number .and. mesh_y==1 .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=8                      !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i-1+grid_x_number                                   !第4个网格为左上侧网格
                grid_near_mesh_name(i,4)=i+grid_x_number                                     !第5个网格为上侧网格


                !下层网格
                grid_near_mesh_name(i,5)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,7)=i-1+grid_x_number-grid_x_number*grid_y_number       !第4个网格为左上侧网格
                grid_near_mesh_name(i,8)=i+grid_x_number-grid_x_number*grid_y_number         !第5个网格为上侧网格


            !顶部右后
            else if(mesh_x==grid_x_number .and. mesh_y==grid_y_number .and. mesh_z==grid_z_number) then
                
                grid_near_mesh_number(i)=8                     !中间的网格有九个网格
                
                !本层网格
                grid_near_mesh_name(i,1)=i                                                   !第1个网格为自身
                grid_near_mesh_name(i,2)=i-1                                                 !第2个网格为左侧网格
                grid_near_mesh_name(i,3)=i-1-grid_x_number                                   !第7个网格为左下侧网格
                grid_near_mesh_name(i,4)=i-grid_x_number                                     !第8个网格为下侧网格


                !下层网格
                grid_near_mesh_name(i,5)=i-grid_x_number*grid_y_number                       !第1个网格为自身
                grid_near_mesh_name(i,6)=i-1-grid_x_number*grid_y_number                     !第2个网格为左侧网格
                grid_near_mesh_name(i,7)=i-1-grid_x_number-grid_x_number*grid_y_number       !第7个网格为左下侧网格
                grid_near_mesh_name(i,8)=i-grid_x_number-grid_x_number*grid_y_number         !第8个网格为下侧网格
            !====================================================================================================

            !----------------------------------------------------------------------------------------------------

            end if  
            
        end do

     else 

        write(*,*) " The Dimension input in Distribute_Backgroud_Grid is not right!"

     endif
    
!     !***********************************************************************************************************
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
!     !***********************************************************************************************************

end subroutine Distribute_Backgroud_Grid