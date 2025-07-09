!**************************************************************************************************************
!  SUBROUTINE : Output_for_tecplot
!
!  PURPOSE    : Save the data for Tecplot
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


subroutine Output_for_tecplot(i_time_step)

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
    integer,intent(in)::i_time_step                                                       ! Current time step

    !----------------------------------------------------------------------------------------------------------


    !----------------------------------------------------------------------------------------------------------
    ! Variables in local subroutine
    integer::i,j,k,L,m                                                                    ! Variables for looping

  
    character(len=20)::Char_time                                                          ! Time character type
    real(kind=8)::Velocity_Module                                                         ! Velocity Module
  
    ! Variables for file opeartion
    integer::File_index                                                                   ! File index
    character(len=100)::File_name                                                         ! File name

    ! Variables for folder opeartion
    Logical::dirExists

    ! Variables for tecplot binary output
    integer::Binary_Tecplot_Port                                                          ! Binary Tecplot Port
    character(len=400)::tecplot_character                                                 ! Tecplot character 
    integer::Output_variable_number                                                       ! Output variable number
    real(kind=4)::EOHMARKER   = 357.0                                                     ! Tecplot EOHMARKER
    real(kind=4)::Zone_marker = 299.0                                                     ! Tecplot Zone marker  
    !----------------------------------------------------------------------------------------------------------


    !==========================================================================================================

    ! Body of subroutine Output_for_tecplot

    !transfer the Current_Processor_ID from integer to character
    write(Char_time,'(F10.6)') Real_time

    ! !Calculate the errors
    ! particle_velocity_error=0.0d0
    ! particle_press_error=0.0d0

    ! do i=1,particle_ture_number

    !    !Maximum velocity magnitude
    !    Velocity_Module=sqrt(DOT_PRODUCT(particle_velocity(i,:),particle_velocity(i,:)))

    !    !errors
    !    particle_velocity_error(i,:)= log10(abs((particle_velocity(i,:)-Analytical_particle_velocity(i,:))/velocity_0))            !Particle velocity error 
    !    particle_press_error(i)= log10(abs((particle_press(i)- Analytical_particle_press(i))/pressure_0))                          !Particle pressure error
    !    particle_velocity_magnitude_error(i)=log10(abs(((Velocity_Module-Analytical_particle_velocity_magnitude(i)))/velocity_0))  !Particle velocity magnitude error

    ! enddo

    ! !==========================================================================================================
    ! !Animation Result in ASCII type (Method 1 and method 2)
    ! !You just need to choose one type (This is for multiPhase Flow)
    ! !All the flow results output in one file
    ! !Colour for different zone (Method 1)
    ! !----------------------------------------------------------------------------------------------------------
    ! !Zone for water
    ! write(Animation_Result_File_Port,*) 'ZONE T=" '//trim(adjustl(Char_time))//' ",'//' I=',particle_ture_number," F=POINT"
    ! !Time shoule be in quotation marks ("")
    
    ! do i=1,particle_ture_number
    !    write(Animation_Result_File_Port,'(6F12.6)') (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(particle_velocity(i,j)/velocity_0,j=1,dim),particle_press(i)/pressure_0
    ! end do
    ! !----------------------------------------------------------------------------------------------------------

    ! !----------------------------------------------------------------------------------------------------------
    ! !Zone for wavemaker particles
    ! write(Animation_Result_File_Port,*) 'ZONE T=" '//trim(adjustl(Char_time))//' ",'//' I=',wave_maker_particle_number," F=POINT"
    ! !Time shoule be in quotation marks ("")
    
    ! do k=1,wave_maker_particle_number
    !    i=k+particle_ture_number
    !    write(Animation_Result_File_Port,'(6F12.6)') (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(particle_velocity(i,j)/velocity_0,j=1,dim),particle_press(i)/pressure_0   
    ! end do
    ! !----------------------------------------------------------------------------------------------------------

    ! !----------------------------------------------------------------------------------------------------------
    ! !Zone for Boundary particles
    ! write(Animation_Result_File_Port,*) 'ZONE T=" '//trim(adjustl(Char_time))//' ",'//' I=',fix_ghost_particle_number," F=POINT"
    ! !Time shoule be in quotation marks ("")
    
    ! do k=1,fix_ghost_particle_number
    !    i=k+particle_ture_number+wave_maker_particle_number
    !    write(Animation_Result_File_Port,'(6F12.6)') (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(particle_velocity(i,j)/velocity_0,j=1,dim),particle_press(i)/pressure_0   
    ! end do
    ! !----------------------------------------------------------------------------------------------------------

    ! !==========================================================================================================




    ! !==========================================================================================================
    ! !All the particle in one zone (Method 2)
    ! !----------------------------------------------------------------------------------------------------------
    ! !Zone for all particles
    ! write(Animation_Result_File_Port,*) 'ZONE T=" '//trim(adjustl(Char_time))//' ",'//' I=',ture_total_particle_number," F=POINT"
    ! !Time shoule be in quotation marks ("")
    
    ! do i=1,ture_total_particle_number
    !     Velocity_Module=sqrt(DOT_PRODUCT(particle_velocity(i,:),particle_velocity(i,:)))
    !     write(Animation_Result_File_Port,'(20F10.4)') (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(particle_velocity(i,j)/velocity_0,j=1,dim),(particle_press(i)-Background_Pressure)/pressure_0,Velocity_Module/velocity_0,particle_vorticity(i,3)*time_0 
    ! end do
    ! !----------------------------------------------------------------------------------------------------------
    
    ! !==========================================================================================================



    !==========================================================================================================
    ! If the case has free surface
    if ( Has_FreeSurface_Or_Not==1 ) then

        !------------------------------------------------------------------------------------------------------
        ! Output the free surface particle
        write(Free_Surface_File_Port,*) 'ZONE T=" '//trim(adjustl(Char_time))//' ",'//' I=',surface_particle_number," F=POINT"
        ! Time shoule be in quotation marks ("")

        do k=1,surface_particle_number
        
            i=surface_particle_name(k)
            write(Free_Surface_File_Port,'(20F10.4)') (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(Normal_vector(i,j),j=1,dim)
        
        end do
        !------------------------------------------------------------------------------------------------------
        
    endif
    !==========================================================================================================












    !==========================================================================================================
    !                                            .--,       .--,
    !                                           ( (  \.---./  ) )
    !                                            '.__/o   o\__.'
    !                                               {=  ^  =}
    !                                                >  -  <
    !                                               /       \
    !                                              //       \\
    !                                             //|   .   |\\
    !                                             "'\       /'"_.-~^`'-.
    !                                                \  _  /--'         `
    !                                              ___)( )(___
    !                                             (((__) (__)))    
    !==========================================================================================================











    !==========================================================================================================

    !----------------------------------------------------------------------------------------------------------
    ! Smoke Line
    write(Smoke_Line_File_Port,*) 'ZONE T=" '//trim(adjustl(Char_time))//' ",'//' I=',Output_particle_number," F=POINT"
    ! Time shoule be in quotation marks ("")
    
    do k=1,Output_particle_number
        
        i=Output_particle_index(k)
        write(Smoke_Line_File_Port,'(20F10.4)') (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),particle_smoke_value(i)
    
    enddo
    !----------------------------------------------------------------------------------------------------------


    !==========================================================================================================







    !==========================================================================================================
    ! Animation Result in Binary type
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
    inquire(directory='Binary_Tecplot', exist=dirExists )                                      
    
    if(dirExists) then

        !------------------------------------------------------------------------------------------------------
        File_name='./Binary_Tecplot/'//trim(adjustl(Char_time))//'.dat'
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

        write(Binary_Tecplot_Port) 0   !FileType 0 = FULL, 1 = GRID, 2 = SOLUTION

        !-----1.3. Title and variable names.-------------------------------------------------------------------

        !-----1.3.1. The TITLE.--------------------------------------------------------------------------------
        tecplot_character="TITLE='"//trim(adjustl(Char_time))//"'"
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))
        
        !-----1.3.2 Number of variables (NumVar) in the datafile.
        !tecplot_character="VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'U' 'Vor' "
        Output_variable_number=7
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

        !-----1.4. Zones---------------------------------------------------------------------------------------

        !-----1.4.1 Zone 1-------------------------------------------------------------------------------------
        
        !--------Zone marker. Value = 299.0 (keep it same for all zones)
        write(Binary_Tecplot_Port) Zone_marker
        
        !--------Zone name. (See note 1.) N = length of zone name + 1.
        tecplot_character="ZONE "//trim(adjustl(Char_time))
        call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))

        write(Binary_Tecplot_Port) -1                       ! parent zone of this zone
        write(Binary_Tecplot_Port) -2                       ! StrandID
        write(Binary_Tecplot_Port) real(Real_time,kind=8)   ! solution time, 
        write(Binary_Tecplot_Port) -1                       ! reserved data position, useless
        write(Binary_Tecplot_Port) 0                        ! ZoneType 0=ORDERED, 1=FELINESEG, 2=FETRIANGLE, 3=FEQUADRILATERAL, 4=FETETRAHEDRON, 5=FEBRICK, 6=FEPOLYGON, 7=FEPOLYHEDRON
        write(Binary_Tecplot_Port) 0                        ! Specify Var Location. 0 = Don't specify, all data is located at the nodes
        write(Binary_Tecplot_Port) 0                        ! Are raw local 1-to-1 face neighbors supplied? (0=FALSE 1=TRUE) ORDERED and FELINESEG zones must specify 0
        write(Binary_Tecplot_Port) 0                        ! Number of miscellaneous user-defined face neighbor connections (value >= 0)
  
        !--------IMax,JMax,KMax
        write(Binary_Tecplot_Port) Output_particle_number,1,1                  !(IMax)
        
        ! For all zone types (repeat for each Auxiliary data name/value pair):
        write(Binary_Tecplot_Port) 0

        !--------I HEADER OVER---------------------------------------------------------------------------------
        
        !--------EOHMARKER, value=357.0
        write(Binary_Tecplot_Port) EOHMARKER
        
        !--------II. Data section------------------------------------------------------------------------------

        !--------2.1 zone 1------------------------------------------------------------------------------------

        !--------Zone marker Value = 299.0
        write(Binary_Tecplot_Port) Zone_marker
        
        !--------Variable data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit           
        do i=1,Output_variable_number
            write(Binary_Tecplot_Port) 2                                                   ! very important
        enddo
        
        write(Binary_Tecplot_Port) 0 ! Has passive variables: 0 = no, 1 = yes.

        !--------Has variable sharing 0 = no, 1 = yes. --------------------------------------------------------    
        write(Binary_Tecplot_Port) 0
                                     
        !--------Zone number to share connectivity list with (-1 = no sharing).--------------------------------
        write(Binary_Tecplot_Port) -1
        
        !--------Zone Data. Each variable is in data format asspecified above.---------------------------------

        !--------Attention: The significant difference between Version TV101 and TV 112------------------------
        ! First:  We need to show the maximum and minimum value for each variables and they should be in order
        ! Second: We cannot output the data in one loop the only choice is output one by one 
        ! Third:  The maximum and minimum value is used for defining the level show, so you also can output an estimation value

        !------------------------------------------------------------------------------------------------------
        ! Write the maximum and mimum value for each Varibales
        do i=1,Output_variable_number
           write(Binary_Tecplot_Port) 0.0d0             !Minimum value            
           write(Binary_Tecplot_Port) 2.0d0             !maximum value  
        enddo

        !------------------------------------------------------------------------------------------------------
        ! Data
        ! Position
        do j=1,dim
            do k=1,Output_particle_number
                i=Output_particle_index(k)
                write(Binary_Tecplot_Port) (particle_position(i,j)/NonDimensionalization_Reference_Length)
            enddo
        enddo

        ! Velocity
        do j=1,dim
            do k=1,Output_particle_number
                i=Output_particle_index(k)
                write(Binary_Tecplot_Port) (particle_velocity(i,j)/velocity_0)
            enddo
        enddo

        ! Pressure
        do k=1,Output_particle_number
            i=Output_particle_index(k)
            write(Binary_Tecplot_Port) (particle_press(i)-Background_Pressure)/pressure_0
        enddo

        !Velocity Module
        do k=1,Output_particle_number
            i=Output_particle_index(k)
            Velocity_Module=sqrt( DOT_PRODUCT( particle_velocity(i,:),particle_velocity(i,:) ) )
            write(Binary_Tecplot_Port) Velocity_Module/velocity_0
        enddo

        !Particle Vorticity Module
        do k=1,Output_particle_number
            i=Output_particle_index(k)
           write(Binary_Tecplot_Port) particle_vorticity(i,3)*time_0   
        end do
        !------------------------------------------------------------------------------------------------------

        close(Binary_Tecplot_Port)

    else

      write(*,*) ' The Binary_Tecplot folder is not exist, you should creat it!'

    endif
    !==========================================================================================================





















    ! !==========================================================================================================
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
    ! inquire(directory='Binary_Tecplot', exist=dirExists )                                      
    
    ! if(dirExists) then

    !     !------------------------------------------------------------------------------------------------------
    !     File_name='./Binary_Tecplot/'//trim(adjustl(Char_time))//'.dat'
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
    !     write(Binary_Tecplot_Port) 1

    !     !-----1.3. Title and variable names.-------------------------------------------------------------------

    !     !-----1.3.1. The TITLE.--------------------------------------------------------------------------------
    !     tecplot_character="TITLE='"//trim(adjustl(Char_time))//"'"
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))

    !     !-----1.3.2 Number of variables (NumVar) in the datafile.----------------------------------------------
    !     !tecplot_character="VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'U' 'Vor' "
    !     Output_variable_number=7
    !     write(Binary_Tecplot_Port) Output_variable_number               !Variable number

    !     !-----1.3.3 Variable names. N = L[1] + L[2] + .... L[NumVar]-------------------------------------------
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

    !     !-----1.4. Zones---------------------------------------------------------------------------------------

    !     !-----1.4.1 Zone 1-------------------------------------------------------------------------------------
        
    !     !--------Zone marker. Value = 299.0 (keep it same for all zones)
    !     write(Binary_Tecplot_Port) Zone_marker
        
    !     !--------Zone name. (See note 1.) N = length of zone name + 1.
    !     tecplot_character="ZONE "//trim(adjustl(Char_time))
    !     call Write_ASCII_to_BINARY(Binary_Tecplot_Port,trim(adjustl(tecplot_character)),len(trim(adjustl(tecplot_character))))

    !     !--------Zone Color (set to -1 if you want tecplot to determine).
    !     write(Binary_Tecplot_Port) -1
        
    !     !--------ZoneType 0=ORDERED,1=FELINESEG,2=FETRIANGLE,3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
    !     write(Binary_Tecplot_Port) 0
        
    !     !--------DataPacking 0=Block, 1=Point
    !     write(Binary_Tecplot_Port) 1
        
    !     !--------Specify Var Location. 0 = Don't specify, all data is located at the nodes. 1 = Specify
    !     write(Binary_Tecplot_Port) 0                
        
    !     !--------Number of user defined face neighbor connections (value >= 0)
    !     write(Binary_Tecplot_Port) 0
        
    !     !--------IMax,JMax,KMax
    !     write(Binary_Tecplot_Port) Output_particle_number,1,1                      !(IMax)
    !     !write(Binary_Tecplot_Port) I_Max,J_Max,K_Max                              !(IMax)
        
    !     ! For all zone types (repeat for each Auxiliary data name/value pair):
    !     write(Binary_Tecplot_Port) 0

    !     !--------I HEADER OVER---------------------------------------------------------------------------------

    !     !--------EOHMARKER, value=357.0
    !     write(Binary_Tecplot_Port) EOHMARKER
        
    !     !--------II. Data section------------------------------------------------------------------------------

    !     !-----2.1 zone 1---------------------------------------------------------------------------------------

    !     !--------Zone marker Value = 299.0
    !     write(Binary_Tecplot_Port) Zone_marker
        
    !     !--------Variable data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit           
    !     do i=1,Output_variable_number
    !         write(Binary_Tecplot_Port) 2                                                   ! very important
    !     enddo
        
    !     !--------Has variable sharing 0 = no, 1 = yes.      
    !     write(Binary_Tecplot_Port) 0
                                     
    !     !--------Zone number to share connectivity list with (-1 = no sharing).
    !     write(Binary_Tecplot_Port) -1
        
    !     !--------Zone Data. Each variable is in data format asspecified above.

    !     do k=1,Output_particle_number
    !         i=Output_particle_index(k)
    !         Velocity_Module=sqrt( DOT_PRODUCT( particle_velocity(i,:),particle_velocity(i,:) ) )
    !         write(Binary_Tecplot_Port) (particle_position(i,j)/NonDimensionalization_Reference_Length,j=1,dim),(particle_velocity(i,j)/velocity_0,j=1,dim),(particle_press(i)-Background_Pressure)/pressure_0,Velocity_Module/velocity_0,particle_vorticity(i,3)*time_0   
    !     enddo

    !     close(Binary_Tecplot_Port)

    ! else

    !   write(*,*) ' The Binary_Tecplot folder is not exist, you should creat it!'

    ! endif
    ! !==========================================================================================================



end subroutine Output_for_tecplot