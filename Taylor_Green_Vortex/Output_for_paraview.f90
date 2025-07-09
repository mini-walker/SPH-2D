!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE:: Output_for_paraview
!
!  PURPOSE: Save the data for Paraview
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

subroutine Output_for_paraview(i_time_step)

    use Public_variable_module
    use information_module
    use MPI

    !==========================================================================================================
    !Variables from mother subroutine
    integer,intent(in)::i_time_step                              !current time step

    !Variables in subroutine
    integer::i,j,k,L,m                                           !loop Variables
    character(len=20)::filename                                  !filename
  
    !Variables for file opeartion
    integer::ioerror=0                                           !open return value
    integer::stat                                                !read return value
    integer::status                                              !allocation return value

    integer::file_index
    character(len=40)::file_name                                 !File name
    character(len=10)::char_time                                 !time character type

    !Variables for folder opeartion
    Logical::dirExists

    !Variables for tecplot binary output
    integer::Binary_Paraview_Port                                !Binary Paraview Port
    character(len=400)::Paraview_character                       !Paraview character 
    character(len=1)::line_feed_character=char(10)               !Paraview line feed character
    character(len=12)::offset                                    !Paraview memory offset
    character(len=9)::int_to_string_1,int_to_string_2            !integer to string 1 and 2
    real(kind=8)::double_float=0.0d0                             !double float
    integer::single_int                                          !single integer
    integer::position_data_size,velocity_data_size,pressure_data_size
    integer::connectivity_data_size,offsets_data_size,types_data_size
    integer::position_offset,velocity_offset,pressure_offset
    integer::connectivity_offset,offsets_offset,types_offset

    !==========================================================================================================

    ! Body of subroutine Output_for_paraview
      
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Check the folder exist or not
    inquire(directory='Paraview', exist=dirExists )                                      
    
    ASCII_Paraview:if(dirExists) then

      !--------------------------------------------------------------------------------------------------------
      !transfer the Current_Processor_ID from integer to character
      write(char_time,'(f12.5)') i_time_step*dt
      !write(*,*) char_time

      file_name='./Paraview/flow.vtu.'//trim(adjustl(char_time))
      !write(*,*) file_name

      open(unit=801,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)
      !--------------------------------------------------------------------------------------------------------

      
      ! !========================================================================================================
      ! !Unstructure grid ASCII Format 1 (xml style)
      ! !--------------------------------------------------------------------------------------------------------
      ! !Output the data
      ! !Header of the file
      ! write (801,'(a)') '<?xml version="1.0"?>'
      ! write (801,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">' 
      ! write (801,'(a)') '<UnstructuredGrid>'

      ! !Point number
      ! write (801,*) '<Piece NumberOfPoints="',ture_total_particle_number,'" NumberOfCells="',ture_total_particle_number,'" >'

      ! !--------------------------------------------------------------------------------------------------------
      ! !Velocity
      ! write (801,'(a)') '<PointData>'
      ! write (801,'(a)') '<DataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="ascii">'
      ! if (dim==2) then
          
      !     do i = 1,ture_total_particle_number
      !        write (801,"(101E20.9E3)") (particle_velocity(i,j)/velocity_0,j=1,dim),0.0000
      !     enddo

      ! else if (dim==3) then

      !     do i = 1,ture_total_particle_number
      !        write (801,"(101E20.9E3)") (particle_velocity(i,j)/velocity_0,j=1,dim)
      !     enddo
        
      ! endif
      ! write (801,'(a)') '</DataArray>'
      ! !--------------------------------------------------------------------------------------------------------


      ! !--------------------------------------------------------------------------------------------------------
      ! !Pressure
      ! write (801,'(a)') '<DataArray type="Float64" Name="Pressure" format="ascii">'

      ! do i = 1,ture_total_particle_number
      !    write (801,"(101E20.9E3)") particle_press(i)/pressure_0
      ! enddo
      ! write (801,'(a)') '</DataArray>'
      ! !--------------------------------------------------------------------------------------------------------

      ! !--------------------------------------------------------------------------------------------------------
      ! !Position
      ! write (801,'(a)') '</PointData>'
      ! write (801,'(a)') '<Points>'
      ! write (801,'(a)') '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'

      ! if (dim==2) then
          
      !     do i = 1,ture_total_particle_number
      !        write (801,"(101E20.9E3)") (particle_position(i,j)/domain_size_y,j=1,dim),0.0000
      !     enddo

      ! else if (dim==3) then

      !     do i = 1,ture_total_particle_number
      !        write (801,"(101E20.9E3)") (particle_position(i,j)/domain_size_y,j=1,dim)
      !     enddo
        
      ! endif

      ! write (801,'(a)') '</DataArray>'
      ! write (801,'(a)') '</Points>'
      ! !--------------------------------------------------------------------------------------------------------

      ! !--------------------------------------------------------------------------------------------------------
      ! !End of the file
      ! write (801,'(a)') '<Cells>'
      ! write (801,'(a)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
      ! do i = 1,ture_total_particle_number
      !    write (801,"(I5)") i + 1 
      ! enddo
      ! write (801,'(a)') '</DataArray>'

      ! write (801,'(a)') '<DataArray type="Int32" Name="offsets" format="ascii">'
      ! do i = 1,ture_total_particle_number
      !    write (801,"(I5)") 1 
      ! enddo
      ! write (801,'(a)') '</DataArray>'

      ! write (801,'(a)') '<DataArray type="Int32" Name="types" format="ascii">'
      ! do i = 1,ture_total_particle_number
      !    write (801,"(I5)") 1 
      ! enddo
      ! write (801,'(a)') '</DataArray>'
      ! write (801,'(a)') '</Cells>'

      ! write (801,'(a)') '</Piece>'
      ! write (801,'(a)') '</UnstructuredGrid>'
      ! write (801,'(a)') '</VTKFile>'
      ! !--------------------------------------------------------------------------------------------------------
  
      ! close(801)
      ! !========================================================================================================


      !========================================================================================================
      !Unstructure grid ASCII Format 2 (xml style)
      !--------------------------------------------------------------------------------------------------------
      !Output the data
      !Header of the file
      write (801,'(a)') '<?xml version="1.0"?>'
      write (801,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">' 
      write (801,'(a)') '<UnstructuredGrid>'

      !Point number
      write (801,*) '<Piece NumberOfPoints="',ture_total_particle_number,'" NumberOfCells="',0,'" >'
      
      !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      write (801,'(a)') '<Points>'          !Points position
      !Data array type (Position)
      write (801,'(a)') '<DataArray name="Position" type="Float64" NumberOfComponents="3" format="ascii">'
      
      if (dim==2) then
          
          do i = 1,ture_total_particle_number
             write (801,"(101E20.9E3)") (particle_position(i,j)/domain_size_y,j=1,dim),0.0000
          enddo

      else if (dim==3) then

          do i = 1,ture_total_particle_number
             write (801,"(101E20.9E3)") (particle_position(i,j)/domain_size_y,j=1,dim)
          enddo
        
      endif
      write (801,'(a)') '</DataArray>'
      write (801,'(a)') '</Points>'         !End Points position
      !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     
      !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      !Point data
      !Velocity
      write (801,'(a)') '<PointData Scalars="Pressure" Vectors="Velocity" >'       !Data at points
      
      write (801,'(a)') '<DataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="ascii">'
      if (dim==2) then
          
          do i = 1,ture_total_particle_number
             write (801,"(101E20.9E3)") (particle_velocity(i,j)/velocity_0,j=1,dim),0.0000
          enddo

      else if (dim==3) then

          do i = 1,ture_total_particle_number
             write (801,"(101E20.9E3)") (particle_velocity(i,j)/velocity_0,j=1,dim)
          enddo
        
      endif
      write (801,'(a)') '</DataArray>'
      !--------------------------------------------------------------------------------------------------------
      
      !--------------------------------------------------------------------------------------------------------
      !Pressure
      write (801,'(a)') '<DataArray type="Float64" Name="Pressure" format="ascii">'

      do i = 1,ture_total_particle_number
         write (801,"(101E20.9E3)") particle_press(i)/pressure_0
      enddo
      write (801,'(a)') '</DataArray>'
      !--------------------------------------------------------------------------------------------------------

      write (801,'(a)') '</PointData>'                      !End Points Data type
      !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      !--------------------------------------------------------------------------------------------------------
      !File end
      write (801,'(a)') '<Cells>'
      write (801,'(a)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
      write (801,'(a)') '</DataArray>'
      write (801,'(a)') '<DataArray type="Int32" Name="offsets" format="ascii">'
      write (801,'(a)') '</DataArray>'
      write (801,'(a)') '<DataArray type="UInt8" Name="types" format="ascii">'
      write (801,'(a)') '</DataArray>'
      write (801,'(a)') '</Cells>'
      write (801,'(a)') '</Piece>'
      write (801,'(a)') '</UnstructuredGrid>'
      write (801,'(a)') '</VTKFile>'
      !--------------------------------------------------------------------------------------------------------

      close(801)
      !========================================================================================================

    else

      write(*,*) ' The Paraview folder is not exist, you should creat it!'

    endif ASCII_Paraview
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






    ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! ! Check the folder exist or not
    ! ! Still have some bugers in the binary vtk write
    ! inquire(directory='Binary_Paraview', exist=dirExists )                                      
    
    ! Binary_Paraview:if(dirExists) then

    !   !--------------------------------------------------------------------------------------------------------
    !   !transfer the Current_Processor_ID from integer to character
    !   write(char_time,'(f8.4)') i_time_step*dt
    !   !write(*,*) char_time

    !   file_name='./Binary_Paraview/flow_binary.vtu.'//trim(adjustl(char_time))
    !   !write(*,*) file_name

    !   Binary_Paraview_Port=801

    !   open(unit=Binary_Paraview_Port,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",form="BINARY",iostat=ioerror)
    !   !--------------------------------------------------------------------------------------------------------

    !   position_data_size=3*ture_total_particle_number*sizeof(double_float)
    !   velocity_data_size=3*ture_total_particle_number*sizeof(double_float)
    !   pressure_data_size=ture_total_particle_number*sizeof(double_float)
    !   connectivity_data_size=ture_total_particle_number*sizeof(single_int)
    !   offsets_data_size=ture_total_particle_number*sizeof(single_int)
    !   types_data_size=ture_total_particle_number*sizeof(single_int)

    !   ! position_offset=0
    !   ! velocity_offset=position_offset+sizeof(single_int)+position_data_size
    !   ! pressure_offset=velocity_offset+sizeof(single_int)+velocity_data_size
    !   ! connectivity_offset=pressure_offset+sizeof(single_int)+pressure_data_size
    !   ! offsets_offset=connectivity_offset+sizeof(single_int)+connectivity_data_size
    !   ! types_offset=offsets_offset+sizeof(single_int)+offsets_data_size

    !   velocity_offset=0
    !   pressure_offset=velocity_offset+sizeof(single_int)+velocity_data_size

    !   position_offset=pressure_offset+sizeof(single_int)+pressure_data_size

    !   connectivity_offset=position_offset+sizeof(single_int)+position_data_size

    !   offsets_offset=connectivity_offset+sizeof(single_int)+connectivity_data_size

    !   types_offset=offsets_offset+sizeof(single_int)+offsets_data_size

    !   !========================================================================================================
    !   !Binary Unstructure gird VTK Format (xml style)
    !   !--------------------------------------------------------------------------------------------------------
    !   !Output the data
    !   !Header of the file
    !   Paraview_character='<?xml version="1.0"?>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   Paraview_character='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   Paraview_character='<UnstructuredGrid>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   !transfer particle and cell number to character

    !   write(int_to_string_1(1:9),'(I9)') ture_total_particle_number          !particle number
    !   write(int_to_string_2(1:9),'(I9)') ture_total_particle_number          !cell number

    !   Paraview_character='<Piece NumberOfPoints="'//int_to_string_1//'" NumberOfCells="'//int_to_string_2//'">'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)
      
    !   !--------------------------------------------------------------------------------------------------------
    !   Paraview_character='<PointData>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   write(offset(1:12),'(i12)') velocity_offset          !position offset

    !   Paraview_character='<DataArray name="Velocity" type="Float64" NumberOfComponents="3" format="appended" offset="'//offset//'" />'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   write(offset(1:12),'(i12)') pressure_offset          !pressure offset

    !   Paraview_character='<DataArray name="Pressure" type="Float64" format="appended" offset="'//offset//'" />'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   Paraview_character='</PointData>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)
    !   !--------------------------------------------------------------------------------------------------------

    !   Paraview_character= '<CellData>  </CellData>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   !-------------------------------------------------------------------------------------------------------
    !   !write the data header and the offsets
    !   Paraview_character='<Points>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   write(offset(1:12),'(i12)') position_offset          !position offset

    !   Paraview_character='<DataArray name="Position" type="Float64" NumberOfComponents="3" format="appended" offset="'//offset//'" />'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   Paraview_character='</Points>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)
    !   !--------------------------------------------------------------------------------------------------------

    !   Paraview_character='<Cells>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   write(offset(1:12),'(i12)') connectivity_offset          !connectivity offset

    !   Paraview_character='<DataArray type="Int32" Name="connectivity" format="appended" offset="'//offset//'" />'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   write(offset(1:12),'(i12)') offsets_offset                !offsets offset

    !   Paraview_character='<DataArray type="Int32" Name="offsets" format="appended" offset="'//offset//'" />'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   write(offset(1:12),'(i12)') types_offset                  !types offset

    !   Paraview_character='<DataArray type="Int32" Name="types" format="appended" offset="'//offset//'" />'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   Paraview_character='</Cells>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   Paraview_character='</Piece>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   Paraview_character='</UnstructuredGrid>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   Paraview_character='<AppendedData encoding="raw">'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   Paraview_character='_' 
    !   write(Binary_Paraview_Port) trim(Paraview_character)
    !   !--------------------------------------------------------------------------------------------------------

    !   write(Binary_Paraview_Port) velocity_data_size

    !   if (dim==2) then
          
    !       do i = 1,ture_total_particle_number
    !          write (Binary_Paraview_Port) (particle_velocity(i,j)/velocity_0,j=1,dim),double_float
    !       enddo

    !   else if (dim==3) then

    !       do i = 1,ture_total_particle_number
    !          write (Binary_Paraview_Port) (particle_velocity(i,j)/velocity_0,j=1,dim)
    !       enddo
        
    !   endif

    !   write(Binary_Paraview_Port) pressure_data_size

    !   do i = 1,ture_total_particle_number
    !      write (Binary_Paraview_Port) particle_press(i)/pressure_0
    !   enddo

    !   write(Binary_Paraview_Port) position_data_size

    !   if (dim==2) then
          
    !       do i = 1,ture_total_particle_number
    !          write (Binary_Paraview_Port) (particle_position(i,j)/domain_size_y,j=1,dim),double_float
    !       enddo

    !   else if (dim==3) then

    !       do i = 1,ture_total_particle_number
    !          write (Binary_Paraview_Port) (particle_position(i,j)/domain_size_y,j=1,dim)
    !       enddo
        
    !   endif

    !   write(Binary_Paraview_Port) connectivity_data_size

    !   do i = 1,ture_total_particle_number
    !      write (Binary_Paraview_Port) i + 1 
    !   enddo

    !   write(Binary_Paraview_Port) offsets_data_size
      
    !   do i = 1,ture_total_particle_number
    !      write (Binary_Paraview_Port) 1 
    !   enddo

    !   write(Binary_Paraview_Port) types_data_size
      
    !   do i = 1,ture_total_particle_number
    !      write (Binary_Paraview_Port) 1 
    !   enddo
     
    !   !--------------------------------------------------------------------------------------------------------
     
    !   !--------------------------------------------------------------------------------------------------------
    !   Paraview_character=line_feed_character//'</AppendedData>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)

    !   Paraview_character='</VTKFile>'//line_feed_character
    !   write(Binary_Paraview_Port) trim(Paraview_character)
    !   !--------------------------------------------------------------------------------------------------------


    !   ! !--------------------------------------------------------------------------------------------------------
    !   ! Paraview_character='# vtk DataFile Version 3.0'//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)
      
    !   ! Paraview_character='vtk output'//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)

    !   ! Paraview_character='BINARY'//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)

    !   ! Paraview_character='DATASET UNSTRUCTURED_GRID'//line_feed_character//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)

    !   ! !--------------------------------------------------------------------------------------------------------

    !   ! !--------------------------------------------------------------------------------------------------------
    !   ! ! POINTS SECTION
    !   ! write(int_to_string_1(1:9),'(I9)') ture_total_particle_number          !particle number

    !   ! Paraview_character= 'POINTS'//int_to_string_1//'float'//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)

    !   ! if (dim==2) then
          
    !   !     do i = 1,ture_total_particle_number
    !   !        write (Binary_Paraview_Port) (particle_position(i,j)/domain_size_y,j=1,dim),double_float
    !   !     enddo

    !   ! else if (dim==3) then

    !   !     do i = 1,ture_total_particle_number
    !   !        write (Binary_Paraview_Port) (particle_position(i,j)/domain_size_y,j=1,dim)
    !   !     enddo
        
    !   ! endif

    !   ! !--------------------------------------------------------------------------------------------------------
    !   ! ! CELLS SECTION
    !   ! write(int_to_string_1(1:9),'(I9)') ture_total_particle_number          !cell number
    !   ! write(int_to_string_2(1:9),'(I9)') ture_total_particle_number*(1)    !cell number

    !   ! Paraview_character=line_feed_character//line_feed_character//'CELLS'//int_to_string_1//' '//int_to_string_2//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)

    !   ! write(Binary_Paraview_Port) ture_total_particle_number

    !   ! do i = 1,ture_total_particle_number
    !   !    write (Binary_Paraview_Port) 1,i
    !   ! enddo
    !   ! !--------------------------------------------------------------------------------------------------------

    !   ! !--------------------------------------------------------------------------------------------------------
    !   ! ! CELL_TYPES SECTION
    !   ! write(int_to_string_2(1:9),'(I9)') ture_total_particle_number                                   !cell number

    !   ! Paraview_character=line_feed_character//line_feed_character//'CELL_TYPES'//int_to_string_2//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)
      
    !   ! do i = 1,ture_total_particle_number
    !   !    write (Binary_Paraview_Port) 1 
    !   ! enddo
    !   ! !--------------------------------------------------------------------------------------------------------
      
    !   ! !--------------------------------------------------------------------------------------------------------
    !   ! ! POINT_DATA SECTION
    !   ! write(int_to_string_1(1:9),'(I9)') ture_total_particle_number          !particle number
      
    !   ! Paraview_character=line_feed_character//line_feed_character//'POINT_DATA'//int_to_string_1//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)

    !   ! Paraview_character='SCALARS pressure float'//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)

    !   ! Paraview_character='LOOKUP_TABLE default'//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)

    !   ! do i = 1,ture_total_particle_number
    !   !    write (Binary_Paraview_Port) particle_press(i)/pressure_0
    !   ! enddo
    !   ! !--------------------------------------------------------------------------------------------------------

    !   ! !--------------------------------------------------------------------------------------------------------
    !   ! Paraview_character=line_feed_character//line_feed_character//'VECTORS velocity float'//line_feed_character ; write(Binary_Paraview_Port) trim(Paraview_character)

    !   ! if (dim==2) then
          
    !   !     do i = 1,ture_total_particle_number
    !   !        write (Binary_Paraview_Port) (particle_velocity(i,j)/velocity_0,j=1,dim),double_float
    !   !     enddo

    !   ! else if (dim==3) then

    !   !     do i = 1,ture_total_particle_number
    !   !        write (Binary_Paraview_Port) (particle_velocity(i,j)/velocity_0,j=1,dim)
    !   !     enddo
        
    !   ! endif
    !   ! !--------------------------------------------------------------------------------------------------------

    !   !========================================================================================================
    !   close(Binary_Paraview_Port)
    
    ! else

    !   write(*,*) ' The Binary_Paraview folder is not exist, you should creat it!'

    ! endif Binary_Paraview

end subroutine Output_for_paraview

