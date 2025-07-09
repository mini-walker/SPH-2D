!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE:: Output_for_Grid_PostProceding
!
!  PURPOSE: Save the Grid data for Post-Proceding
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

subroutine Output_for_Grid_PostProceding(i_time_step)

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


    !---------------------------------------------------------------------------------------------------------
    !Variables for output
    real(kind=8)::current_time                                   !current time
    character(len=20)::char_time                                 !time character type
  
    !Variables for file opeartion
    integer::ioerror=0                                           !open return value
    integer::stat                                                !read return value
    integer::status                                              !allocation return value

    integer::file_index
    character(len=40)::file_name                                 !File name

    !Variables for folder opeartion
    Logical::dirExists                                           !Direction Exists or Not
    real(kind=8)::Velocity_Module                                !Velocity Module

    !Variables for tecplot binary output
    integer::Binary_Tecplot_Port                                 !Binary Tecplot Port
    character(len=400)::tecplot_character                        !Tecplot character 
    real(kind=4):: EOHMARKER = 357.0                             !Tecplot EOHMARKER
    real(kind=4):: Zone_marker = 299.0                           !Tecplot Zone marker
    !----------------------------------------------------------------------------------------------------------

    !==========================================================================================================

    ! Body of subroutine Output_for_Grid_PostProceding


    !==========================================================================================================
    !Output the element data
    !----------------------------------------------------------------------------------------------------------
    !Output for ASCII Tecplot Grid
    !current time
    current_time=(i_time_step-standing_time_step)*dt
      
    !transfer the Current_Processor_ID from integer to character
    write(char_time,'(f12.5)') current_time

    element_node_velocity_error=0.0d0
    element_node_velocity_error=0.0d0
    element_node_velocity_magnitude_error=0.0d0

    !Calculate errors
    do i=1,total_element_node_number

       !Maximum velocity magnitude
       Velocity_Module=sqrt(DOT_PRODUCT(element_node_velocity(i,:),element_node_velocity(i,:)))

       !errors
       element_node_velocity_error(i,:)= log10(abs((element_node_velocity(i,:)-Analytical_element_node_velocity(i,:))/velocity_0))            !Element node velocity error 
       element_node_press_error(i)= log10(abs((element_node_press(i)- Analytical_element_node_press(i))/pressure_0))                          !Element node pressure error
       element_node_velocity_magnitude_error(i)=log10(abs((Velocity_Module-Analytical_element_node_velocity_magnitude(i))/velocity_0))        !Element node velocity magnitude error

    enddo

    write(121,*) 'ZONE T=" '//trim(adjustl(char_time))//' ",'//"N=",total_element_node_number,",E=",total_element_number,",DATAPACKING=POINT,ZONETYPE=FEQUADRILATERAL "   
    
    !域名，节点数，单元数，数据是在点上还是块上（数据是点还是块），单元类型
    !单元类型：FETRIANGLE(三角形网格)、FEQUADRILATERAL(四节点四边形网格)、FETETRAHEDRON(四节点四面体网格)、FEBRICK(8节点六面体)
    !对于混合网格需要不同的ZONES
    
    do i=1,total_element_node_number
        Velocity_Module=sqrt(element_node_velocity(i,1)**2+element_node_velocity(i,2)**2)
        write(121,100) (element_node_position(i,j)/domain_size_y,j=1,dim),(element_node_velocity(i,j)/velocity_0,j=1,dim),element_node_press(i)/pressure_0,Velocity_Module/velocity_0,element_node_vorticity(i,3)*time_0,(element_node_velocity_error(i,j),j=1,dim),element_node_press_error(i),element_node_velocity_magnitude_error(i),element_node_Fai(i)
100     format(20F20.8) 
    end do
    
    do i=1,total_element_number
        write(121,'(8I8)') (element_node_index(i,j),j=1,NodeNumberInOneElement)
    end do
    !----------------------------------------------------------------------------------------------------------


    !----------------------------------------------------------------------------------------------------------
    !Output for Binary Tecplot Grid

    ! Binary Data
    ! Tecplot binary strcuture:
    ! +----------------+
    ! | HEADER SECTION |
    ! +----------------+
    ! +---------+
    ! |FLOAT32 | EOHMARKER, value=357.0
    ! +---------+
    ! +----------------+
    ! | DATA SECTION |
    ! +----------------+
    ! Check the folder exist or not
    inquire(directory='Binary_Tecplot_Grid', exist=dirExists )                                      
    
    if(dirExists) then

        !------------------------------------------------------------------------------------------------------
        !transfer the current time from float to character
        write(char_time,'(f8.4)') (i_time_step-standing_time_step)*dt
        !write(*,*) char_time

        file_name='./Binary_Tecplot_Grid/'//trim(adjustl(char_time))//'.dat'
        !write(*,*) file_name

        open(unit=802,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",form="BINARY",iostat=ioerror)
        !------------------------------------------------------------------------------------------------------

        Binary_Tecplot_Port=802

        !    1. The header section.
        !    The header section contains The version number of the file, a title
        !    of the file, the names of the varialles to be plotted, the
        !    descriptions of all zones to be read in and all text and geometry
        !    definitions.
    
        !------1.1 Magic number, Version number
        write(Binary_Tecplot_Port) "#!TDV101"
        
        !-----1.2. Integer value of 1.----------------------------------------------------------
        !         +-----------+
        !         | INT32     |       This is used to determine the byte order
        !         +-----------+       of the reader relative to the writer.
        write(Binary_Tecplot_Port) 1
    
        !-----1.3. Title and variable names.-------------------------------------------------
        !-----1.3.1. The TITLE. 
        tecplot_character="TITLE='"//trim(adjustl(char_time))//"'"
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        
        !-----1.3.2 Number of variables (NumVar) in the datafile.
        !tecplot_character="VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
        write(Binary_Tecplot_Port) 7                                      !Variable number
        !------1.3.3 Variable names. N = L[1] + L[2] + .... L[NumVar]
        tecplot_character='X' 
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        tecplot_character='Y' 
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        tecplot_character='VX' 
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        tecplot_character='VY'
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        tecplot_character='P'
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        tecplot_character='U'
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        tecplot_character='Fai'
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        

        !-----1.4. Zones-------------------------------------------------------------------
        !-----1.4.1 Zone 1+++++++++++++++++++++++++++++++++++++++++++
        !--------Zone marker. Value = 299.0 (keep it same for all zones)
        write(Binary_Tecplot_Port) Zone_marker
        
        !--------Zone name. (See note 1.) N = length of zone name + 1.
        tecplot_character="ZONE "//trim(adjustl(char_time))
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))

        !---------Zone Color (set to -1 if you want tecplot to determine).
        write(Binary_Tecplot_Port) -1
        
        !---------ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
        write(Binary_Tecplot_Port) 0
        
        !---------DataPacking 0=Block, 1=Point
        write(Binary_Tecplot_Port) 1
        
        !---------Specify Var Location. 0 = Don't specify, all data is located at the nodes. 1 = Specify
        write(Binary_Tecplot_Port) 0                
        
        !---------Number of user defined face neighbor connections (value >= 0)
        write(Binary_Tecplot_Port) 0
        
        !---------IMax,JMax,KMax
        write(Binary_Tecplot_Port) element_node_number_x,element_node_number_y,1   
        !write(Binary_Tecplot_Port) I_Max,J_Max,K_Max                              
        
        ! For all zone types (repeat for each Auxiliary data name/value pair):
        write(Binary_Tecplot_Port) 0

        !----I HEADER OVER--------------------------------------------------------------------------------------------
        
        !        EOHMARKER, value=357.0
        write(Binary_Tecplot_Port) EOHMARKER
        
        !++++++++II. Data section+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !-----2.1 zone 1-----------------------------------------------------------------------
        !--------Zone marker Value = 299.0
        write(Binary_Tecplot_Port) Zone_marker
        
        !--------variable data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit           
        write(Binary_Tecplot_Port) 2                                                   ! very important
        write(Binary_Tecplot_Port) 2
        write(Binary_Tecplot_Port) 2
        write(Binary_Tecplot_Port) 2
        write(Binary_Tecplot_Port) 2
        write(Binary_Tecplot_Port) 2
        write(Binary_Tecplot_Port) 2

        !--------Has variable sharing 0 = no, 1 = yes.      
        write(Binary_Tecplot_Port) 0
                                     
        !----------Zone number to share connectivity list with (-1 = no sharing).
        write(Binary_Tecplot_Port) -1
        
        !----------Zone Data. Each variable is in data format asspecified above.
        I=0
        do  K=1,element_node_number_y
            do  L=1,element_node_number_x

                I=I+1
                
                Velocity_Module=sqrt(DOT_PRODUCT(element_node_velocity(i,:),element_node_velocity(i,:)))
                write(Binary_Tecplot_Port) (element_node_position(i,j)/domain_size_y,j=1,dim),(element_node_velocity(i,j)/velocity_0,j=1,dim),element_node_press(i)/pressure_0,Velocity_Module/velocity_0,element_node_Fai(i)

            enddo
        enddo

        close(Binary_Tecplot_Port)

    else

      write(*,*) ' The Binary_Tecplot folder is not exist, you should creat it!'

    endif
    !----------------------------------------------------------------------------------------------------------

    !==========================================================================================================


end subroutine Output_for_Grid_PostProceding