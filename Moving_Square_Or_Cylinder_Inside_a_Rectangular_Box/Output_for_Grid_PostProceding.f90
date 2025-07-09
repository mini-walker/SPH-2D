!**************************************************************************************************************
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
!**************************************************************************************************************


subroutine Output_for_Grid_PostProceding(i_time_step)

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI

    implicit none

    !==========================================================================================================
    
    !----------------------------------------------------------------------------------------------------------
    ! Variables from superior subroutine
    integer,intent(in)::i_time_step                                                     ! Current time step
    !----------------------------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------------------------
    ! Variables in local subroutine
    integer::i,j,k,L,m                                                                  ! Loop Variables
    
    ! Variables for output
    character(len=20)::char_time                                                        ! Time character type
  
    ! Variables for file opeartion
    integer::File_index                                                                 ! File index
    character(len=100)::File_name                                                       ! File name

    ! Variables for folder opeartion
    Logical::dirExists                                                                  ! Direction Exists or Not

    real(kind=8)::Velocity_Module                                                       ! Velocity Module

    ! Variables for tecplot binary output
    integer::Binary_Tecplot_Port                                                        ! Binary Tecplot Port
    character(len=400)::tecplot_character                                               ! Tecplot character 
    integer::Output_variable_number                                                     ! Output variable number
    real(kind=4):: EOHMARKER   = 357.0                                                  ! Tecplot EOHMARKER
    real(kind=4):: Zone_marker = 299.0                                                  ! Tecplot Zone Marker
    !----------------------------------------------------------------------------------------------------------

    !==========================================================================================================


    ! Body of subroutine Output_for_Grid_PostProceding


    !==========================================================================================================
    ! Transfer the real time from real to character
    write(char_time,'(F10.6)') Real_time
    
    !----------------------------------------------------------------------------------------------------------
    ! ! Calculate the velocity error
    ! element_node_velocity_error=0.0d0
    ! element_node_velocity_error=0.0d0
    ! element_node_velocity_magnitude_error=0.0d0

    ! !Calculate errors
    ! do i=1,total_element_node_number

    !    !Maximum velocity magnitude
    !    Velocity_Module=sqrt(DOT_PRODUCT(element_node_velocity(i,:),element_node_velocity(i,:)))

    !    !errors
    !    element_node_velocity_error(i,:)= log10(abs((element_node_velocity(i,:)-Analytical_element_node_velocity(i,:))/velocity_0))            !Element node velocity error 
    !    element_node_press_error(i)= log10(abs((element_node_press(i)- Analytical_element_node_press(i))/pressure_0))                          !Element node pressure error
    !    element_node_velocity_magnitude_error(i)=log10(abs((Velocity_Module-Analytical_element_node_velocity_magnitude(i))/velocity_0))        !Element node velocity magnitude error

    ! enddo
    !==========================================================================================================



    ! !==========================================================================================================

    ! !----------------------------------------------------------------------------------------------------------
    ! ! Output Grid results in Tecplot ASCII in one file

    ! write(Animation_Grid_Result_File_Port,*) 'ZONE T=" '//trim(adjustl(char_time))//' ",'//"N=",total_element_node_number,",E=",total_element_number,",DATAPACKING=POINT,ZONETYPE=FEQUADRILATERAL "   
    
    ! ! Domain name, number of nodes, number of cells, whether the data is on points or blocks (whether the data is points or blocks), cell type
    ! ! Cell type: FETRIANGLE, FEQUADRILATERAL, FETETRAHEDRON, FEBRICK
    ! ! For the hybrid grid, we should output the grid in different ZONES 
    
    ! do i=1,total_element_node_number
    !     Velocity_Module=sqrt(DOT_PRODUCT(element_node_velocity(i,:),element_node_velocity(i,:)))                    
    !     write(Animation_Grid_Result_File_Port,'(20F10.4)') (element_node_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(element_node_velocity(i,j)/velocity_0,j=1,dim),(element_node_press(i)-Background_Pressure)/pressure_0,Velocity_Module/velocity_0,element_node_vorticity(i,3)*time_0,element_node_Fai(i)
    ! end do
    
    ! do i=1,total_element_number
    !     write(Animation_Grid_Result_File_Port,'(8I8)') (element_node_index(i,j),j=1,NodeNumberInOneElement)
    ! end do
    ! !----------------------------------------------------------------------------------------------------------

    ! !==========================================================================================================



    !==========================================================================================================
    ! Output the Grid results in Tecplot Binary format

    ! You can choose either TDV112 or TDV101 version (Suggestion: TDV112 as it has time variable) 
    ! Binary Data (TDV112 Version)
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
        File_name='./Binary_Tecplot_Grid/'//trim(adjustl(char_time))//'.dat'
        !write(*,*) File_name

        Binary_Tecplot_Port=Tecplot_General_File_Port

        open(unit=Binary_Tecplot_Port,file=trim(adjustl(File_name)),status="replace",position="rewind",action="write",form="BINARY",iostat=ioerror)
        !------------------------------------------------------------------------------------------------------

        !    1. The header section.
        !    The header section contains The version number of the file, a title
        !    of the file, the names of the varialles to be plotted, the
        !    descriptions of all zones to be read in and all text and geometry
        !    definitions.

        !-----1.1 Magic number, Version number-----------------------------------------------------------------
        write(Binary_Tecplot_Port) "#!TDV112"

        !-----1.2. Integer value of 1.-------------------------------------------------------------------------
        !         +-----------+
        !         | INT32     |       This is used to determine the byte order
        !         +-----------+       of the reader relative to the writer.
        
        write(Binary_Tecplot_Port) 1   !This is used to determine the byte order of the reader, relative to the writer.

        write(Binary_Tecplot_Port) 0   !FileType 0 = FULL, 1 = GRID, 2 = SOLUTION (Full contains both grid and solution)

        !-----1.3. Title and variable names.-------------------------------------------------------------------
        !-----1.3.1. The TITLE. 
        tecplot_character="TITLE='"//trim(adjustl(char_time))//"'"
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        
        !-----1.3.2 Number of variables (NumVar) in the datafile.
        !tecplot_character="VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'U' 'Vor' 'EVX' 'EVY' 'EP' 'EU' 'Fai'"

        Output_variable_number=8
        write(Binary_Tecplot_Port) Output_variable_number               !Variable number

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
        tecplot_character='Vor'
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        tecplot_character='Fai'
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        

        !-----1.4. Zones---------------------------------------------------------------------------------------

        !-----1.4.1 Zone 1+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !--------Zone marker. Value = 299.0 (keep it same for all zones)
        write(Binary_Tecplot_Port) Zone_marker
        
        !--------Zone name. (See note 1.) N = length of zone name + 1.
        tecplot_character="ZONE "//trim(adjustl(char_time))
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))

        write(Binary_Tecplot_Port) -1                       ! parent zone of this zone
        write(Binary_Tecplot_Port) -2                       ! StrandID
        write(Binary_Tecplot_Port) real(Real_time,kind=8)   ! solution time, 
        write(Binary_Tecplot_Port) -1                       ! reserved data position, useless
        write(Binary_Tecplot_Port) 3                        ! ZoneType 0=ORDERED, 1=FELINESEG, 2=FETRIANGLE, 3=FEQUADRILATERAL, 4=FETETRAHEDRON, 5=FEBRICK, 6=FEPOLYGON, 7=FEPOLYHEDRON
        write(Binary_Tecplot_Port) 0                        ! Specify Var Location. 0 = Don't specify, all data is located at the nodes
        write(Binary_Tecplot_Port) 0                        ! Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE) ORDERED and FELINESEG zones must specify 0
        write(Binary_Tecplot_Port) 0                        ! Number of miscellaneous user-defined face neighbor connections (value >= 0)

        !----------IMax,JMax,KMax (For ORDERED)
        !write(Binary_Tecplot_Port) ture_total_particle_number,1,1                  !(IMax,JMax,KMax Ordered zone)

        !----------NumPts,NumElements (For Finite element)
        write(Binary_Tecplot_Port) total_element_node_number,total_element_number   !(NumPts,NumElements)

        !ICellDim,JCellDim,KCellDim (for future use; set to zero)
        write(Binary_Tecplot_Port) 0,0,0                                            !(ICellDim,JCellDim,KCellDim)
        !------------------------------------------------------------------------------------------------------

        ! For all zone types (repeat for each Auxiliary data name/value pair):
        write(Binary_Tecplot_Port) 0

        !----------I HEADER OVER-------------------------------------------------------------------------------
        ! EOHMARKER, value=357.0
        write(Binary_Tecplot_Port) EOHMARKER

        !++++++++II. Data section++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !-----2.1 zone 1---------------------------------------------------------------------------------------

        !----------Zone marker Value = 299.0
        write(Binary_Tecplot_Port) Zone_marker

        !----------Variable data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
        ! very important           
        do i=1,Output_variable_number
            write(Binary_Tecplot_Port) 2                                                   ! very important
        enddo

        write(Binary_Tecplot_Port) 0 ! Has passive variables: 0 = no, 1 = yes.

        !----------Has variable sharing 0 = no, 1 = yes. ------------------------------------------------------     
        write(Binary_Tecplot_Port) 0
                                     
        !----------Zone number to share connectivity list with (-1 = no sharing).
        write(Binary_Tecplot_Port) -1


        !----------Attention: The significant difference between Version TV101 and TV112-----------------------
        ! First:  We need to show the maximum and minimum value for each variables and they should be in order
        ! Second: We cannot output the data in one loop the only choice is output one by one 
        ! Third:  The maximum and minimum value is used for defining the level show, so you also can output an estimation value
        ! Fourth: The grid index counting from 0 which means the grid node index should minus 1 in TV112

        !------------------------------------------------------------------------------------------------------
        !Write the maximum and mimum value for each Varibales
        do i=1,Output_variable_number
           write(Binary_Tecplot_Port) 0.0d0             !Minimum value            
           write(Binary_Tecplot_Port) 2.0d0             !maximum value  
        enddo

        !------------------------------------------------------------------------------------------------------

        !----------Zone Data. Each variable is in data format as specified above.

        ! Method 1: Out the grid nodes in order and then output the point index in per grid
        
        ! Data
        ! Position
        do j=1,dim
            do i=1,total_element_node_number
               write(Binary_Tecplot_Port) element_node_position(i,j)/NonDimensionalization_Reference_Length
            enddo
        enddo

        ! Velocity
        do j=1,dim
            do i=1,total_element_node_number
               write(Binary_Tecplot_Port) element_node_velocity(i,j)/velocity_0
            enddo
        enddo

        ! Pressure
        do i=1,total_element_node_number
           write(Binary_Tecplot_Port) (element_node_press(i)-Background_Pressure)/pressure_0
        enddo

        ! Velocity Module
        do i=1,total_element_node_number
            Velocity_Module=sqrt( DOT_PRODUCT( element_node_velocity(i,:),element_node_velocity(i,:) ) )
            write(Binary_Tecplot_Port) Velocity_Module/velocity_0
        enddo

        ! Particle Vorticity Module （Vor*D/U=Vor*time_0）
        do i=1,total_element_node_number
           write(Binary_Tecplot_Port) element_node_vorticity(i,3)*time_0   
        end do

        ! Fai value
        do i=1,total_element_node_number
           write(Binary_Tecplot_Port) element_node_Fai(i)
        end do

        ! iii. specific to fe zones if ZoneType is NOT FEPOLYGON or FEPOLYHEDRON: if "zone number to share connectivity lists with" == -1
        
        ! if ZoneType is NOT FEPOLYGON or FEPOLYHEDRON: INT32*N
        ! Zone Connectivity Data N=L*JMax  (JMax is the mesh number)
        ! This represents JMax sets of adjacency zero based indices (Zero based means the index counting from 0)
        ! where each set contains L values and L is
        ! 2 for LINESEGS
        ! 3 for TRIANGLES
        ! 4 for QUADRILATERALS
        ! 4 for TETRAHEDRONS
        ! 8 for BRICKS

        ! Attention: The grid index counting from 0 which means the grid node index should minus 1 in TV112
        do i=1,total_element_number                      ! Grid number
            write(Binary_Tecplot_Port) ((element_node_index(i,j)-1),j=1,NodeNumberInOneElement) 
        enddo
        !------------------------------------------------------------------------------------------------------

        close(Binary_Tecplot_Port)

    else

      write(*,*) ' The Binary_Tecplot folder is not exist, you should creat it!'

    endif
    !----------------------------------------------------------------------------------------------------------

    !==========================================================================================================











    ! !==========================================================================================================
    ! !Output for Binary Tecplot Grid

    ! ! Binary Data (TDV101 Version)
    ! ! Tecplot binary strcuture:
    ! ! +----------------+
    ! ! | HEADER SECTION |
    ! ! +----------------+
    ! ! +---------+
    ! ! |FLOAT32 | EOHMARKER, value=357.0
    ! ! +---------+
    ! ! +----------------+
    ! ! | DATA SECTION |
    ! ! +----------------+
    ! ! Check the folder exist or not
    ! inquire(directory='Binary_Tecplot_Grid', exist=dirExists )                                      
    
    ! if(dirExists) then

    !     !------------------------------------------------------------------------------------------------------
    !     File_name='./Binary_Tecplot_Grid/'//trim(adjustl(char_time))//'.dat'
    !     !write(*,*) File_name

    !     Binary_Tecplot_Port=Tecplot_General_File_Port

    !     open(unit=Binary_Tecplot_Port,file=trim(adjustl(File_name)),status="replace",position="rewind",action="write",form="BINARY",iostat=ioerror)
    !     !------------------------------------------------------------------------------------------------------

    !     !    1. The header section.
    !     !    The header section contains The version number of the file, a title
    !     !    of the file, the names of the varialles to be plotted, the
    !     !    descriptions of all zones to be read in and all text and geometry
    !     !    definitions.

    !     !-----1.1 Magic number, Version number-----------------------------------------------------------------
    !     write(Binary_Tecplot_Port) "#!TDV101"

    !     !-----1.2. Integer value of 1.-------------------------------------------------------------------------
    !     !         +-----------+
    !     !         | INT32     |       This is used to determine the byte order
    !     !         +-----------+       of the reader relative to the writer.
        
    !     write(Binary_Tecplot_Port) 1   !This is used to determine the byte order of the reader, relative to the writer.

    !     !-----1.3. Title and variable names.-------------------------------------------------------------------
    !     !-----1.3.1. The TITLE. 
    !     tecplot_character="TITLE='"//trim(adjustl(char_time))//"'"
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        
    !     !-----1.3.2 Number of variables (NumVar) in the datafile.
    !     !tecplot_character="VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'U' 'Vor' 'EVX' 'EVY' 'EP' 'EU' 'Fai'"

    !     Output_variable_number=8
    !     write(Binary_Tecplot_Port) Output_variable_number               !Variable number

    !     !------1.3.3 Variable names. N = L[1] + L[2] + .... L[NumVar]
    !     tecplot_character='X' 
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
    !     tecplot_character='Y' 
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
    !     tecplot_character='VX' 
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
    !     tecplot_character='VY'
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
    !     tecplot_character='P'
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
    !     tecplot_character='U'
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
    !     tecplot_character='Vor'
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
    !     tecplot_character='Fai'
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        

    !     !-----1.4. Zones---------------------------------------------------------------------------------------

    !     !-----1.4.1 Zone 1+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !     !---------Zone marker. Value = 299.0 (keep it same for all zones)
    !     write(Binary_Tecplot_Port) Zone_marker
        
    !     !---------Zone name. (See note 1.) N = length of zone name + 1.
    !     tecplot_character="ZONE "//trim(adjustl(char_time))
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))

    !     !---------Zone Color (set to -1 if you want tecplot to determine).
    !     write(Binary_Tecplot_Port) -1
        
    !     !---------ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
    !     write(Binary_Tecplot_Port) 3
        
    !     !---------DataPacking 0=Block, 1=Point
    !     write(Binary_Tecplot_Port) 1
        
    !     !---------Specify Var Location. 0 = Don't specify, all data is located at the nodes. 1 = Specify
    !     write(Binary_Tecplot_Port) 0                
        
    !     !---------Number of user defined face neighbor connections (value >= 0)
    !     write(Binary_Tecplot_Port) 0
        
    !     !---------IMax,JMax,KMax (For points)
    !     !write(Binary_Tecplot_Port) ture_total_particle_number,1,1                  !(IMax,JMax,KMax Ordered zone)
        
    !     !---------NumPts,NumElements (For Finite element)
    !     write(Binary_Tecplot_Port) total_element_node_number,total_element_number   !(NumPts,NumElements)
        
    !     !ICellDim,JCellDim,KCellDim (for future use; set to zero)
    !     write(Binary_Tecplot_Port) 0,0,0                                            !(ICellDim,JCellDim,KCellDim)

    !     ! For all zone types (repeat for each Auxiliary data name/value pair):
    !     write(Binary_Tecplot_Port) 0

    !     !----I HEADER OVER-------------------------------------------------------------------------------------

    !     ! EOHMARKER, value=357.0
    !     write(Binary_Tecplot_Port) EOHMARKER


    !     !++++++++II. Data section++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !     !-----2.1 zone 1---------------------------------------------------------------------------------------
    !     !--------Zone marker Value = 299.0
    !     write(Binary_Tecplot_Port) Zone_marker
        
    !     !-----------Variable data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
    !     ! Very important           
    !     do i=1,Output_variable_number
    !         write(Binary_Tecplot_Port) 2                                                   ! very important
    !     enddo

    !     !----------Has variable sharing 0 = no, 1 = yes. ------------------------------------------------------   
    !     write(Binary_Tecplot_Port) 0
                                     
    !     !----------Zone number to share connectivity list with (-1 = no sharing).
    !     write(Binary_Tecplot_Port) -1


    !     !------------------------------------------------------------------------------------------------------
    !     !----------Zone Data. Each variable is in data format as specified above.

    !     !Method 1: Out the grid nodes in order and then output the point index in per grid
        
    !     ! i Zone Data.
    !     do i=1,total_element_node_number
    !         Velocity_Module=sqrt( DOT_PRODUCT( element_node_velocity(i,:),element_node_velocity(i,:) ) )                    
    !         write(Binary_Tecplot_Port) (element_node_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(element_node_velocity(i,j)/velocity_0,j=1,dim),(element_node_press(i)-Background_Pressure)/pressure_0,Velocity_Module/velocity_0,element_node_vorticity(i,3)*time_0,element_node_Fai(i)
    !     enddo
        
    !     ! iii specific to fe zones (if "zone number to share connectivity lists with"== -1)
    !     do i=1,total_element_number
    !         write(Binary_Tecplot_Port) (element_node_index(i,j),j=1,NodeNumberInOneElement)
    !     enddo


    !     close(Binary_Tecplot_Port)

    ! else

    !   write(*,*) ' The Binary_Tecplot folder is not exist, you should creat it!'

    ! endif
    ! !----------------------------------------------------------------------------------------------------------

    ! !==========================================================================================================



end subroutine Output_for_Grid_PostProceding