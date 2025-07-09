!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  Subroutine: Distribute_Grid_For_Postproceeding 
!
!  PURPOSE: Distribute Grid For Postproceeding
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
subroutine Distribute_Grid_For_Postproceeding()

    use information_module
    use Public_variable_module
    use MPI
    
    implicit none
    !==========================================================================================================
    
    !----------------------------------------------------------------------------------------------------------
    !Variables in subroutine
    integer::i,j,k,L,m                                                  !loop Variables
    integer::x_n,y_n,z_n                                                !Backgroud Grid Location
    
    !Read file varibales
    integer::ioerror=0                                                  !open file return value
    integer::stat                                                       !read file data return value
    integer::status                                                     !allocate memory status

    !==========================================================================================================

    !Call MPI functions
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )

    !For 2D case
    if (dim==2) then

        !------------------------------------------------------------------------------
        !Generate the grid
        !Distribute the grid for post-proceeding
        !The grid origin is same with the Backgroud Grid

        !Generate the node for element
        k=0
        do i=1,element_node_number_y
            do j=1,element_node_number_x

               k=k+1

               element_node_position(k,1)=(j-1)*element_dx !+0.5*element_dx
               element_node_position(k,2)=(i-1)*element_dx !+0.5*element_dx

               element_node_Fai(k)=0.0d0

            enddo   
        enddo

        !Connect grid for each element
        NodeNumberInOneElement=4
        
        do i=1,element_number_y
            do j=1,element_number_x            !1-m column
                
                k=j+(i-1)*element_number_x     !Element index
                
                !网格左下角第一个节点编号为网格编号加上当前行数-1(Anti-clockwise)
                L=k+(i-1)

                element_node_index(k,1)=L
                element_node_index(k,2)=L+1
                element_node_index(k,3)=L+element_node_number_x+1
                element_node_index(k,4)=L+element_node_number_x

            end do
        
        end do
        !------------------------------------------------------------------------------

        !------------------------------------------------------------------------------
        !Output the initial post-proceeding grid 
        output:if (Current_Processor_ID==0) then

            write(*,*) "Post proceeding grid intervel dx is: ",element_dx

            write(*,*) "Particle intervel dx is:",interior_dx

            !--------------------------------------------------------------------------
            !Grid only
            open(unit=5,file="Initial_Post_Proceeding_Grid.dat",status="replace",position="rewind",action="write",iostat=ioerror)
            
            !Input tecplot header
            write(5,*) "TITLE='mesh'"                      !Tecplot title
            write(5,*) "VARIABLES= 'X' 'Y'"                !VARIABLES number
            write(5,*) "ZONE T= '0' ,N=",total_element_node_number,",E=",total_element_number,",DATAPACKING=POINT,ZONETYPE=FEQUADRILATERAL "   
            
            !域名，节点数，单元数，数据是在点上还是块上（数据是点还是块），单元类型
            !单元类型：FETRIANGLE(三角形网格)、FEQUADRILATERAL(四节点四边形网格)、FETETRAHEDRON(四节点四面体网格)、FEBRICK(8节点六面体)
            !对于混合网格需要不同的ZONES
            
            do i=1,total_element_node_number
                write(5,100) (element_node_position(i,j),j=1,dim)
100             format(8F20.10) 
            end do
            
            do i=1,total_element_number
                write(5,110) (element_node_index(i,j),j=1,NodeNumberInOneElement)
110             format(8I8) 
            end do

            close(5)
            !--------------------------------------------------------------------------


            !--------------------------------------------------------------------------
            !Node only
            open(unit=5,file="Initial_Post_Proceeding_Node.dat",status="replace",position="rewind",action="write",iostat=ioerror)
            
            !Input tecplot header
            write(5,*) "TITLE='Node'"                      !Tecplot title
            write(5,*) "VARIABLES= 'X' 'Y'"                !VARIABLES number

            write(5,*) "ZONE I=",total_element_node_number," F=POINT"
            
            do i=1,total_element_node_number
                write(5,100) (element_node_position(i,j),j=1,dim)
            end do

            close(5)
            !--------------------------------------------------------------------------

        endif output

        !------------------------------------------------------------------------------
        !Put all the grid node in the Backgroud grid
        do i=1,total_element_node_number

            x_n=floor((element_node_position(i,1)-chain_origin_x)/chain_dx+1)
            y_n=floor((element_node_position(i,2)-chain_origin_y)/chain_dy+1)

            !Calculate the Backgroud Grid Index
            k=(y_n-1)*chain_x_number+x_n                       !one row by row

            element_node_chain_number(i)=k                     !Save the element node location in Backgroud grid

        enddo
        !------------------------------------------------------------------------------

    !For 3D case
    elseif (dim==3) then

        !Connect grid for each element
        NodeNumberInOneElement=8
        
    endif


!     !***********************************************************************************************************
!     !For checking
    
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


end subroutine Distribute_Grid_For_Postproceeding