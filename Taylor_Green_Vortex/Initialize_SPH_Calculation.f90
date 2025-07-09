!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  Subroutine: Initialize_SPH_Calculation
!
!  PURPOSE: Set the initial information from initial or latest domain data for SPH calculation
!           (contains: particles deviding and backgroud mesh)
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
!  Input: domain size and packing information
!
!  Output: all particles and background mesh information
!
!  Note: MPI version: mpich-3.2
!        Fortran version: Inter Fortran (ifort) 18.0.0
!****************************************************************************

subroutine Initialize_SPH_Calculation()

    use information_module
    use Public_variable_module
    use MPI
    
    implicit none

    !Variables in subroutine

    !==========================================================================================================
    integer::i,j,k,L,m                                                  !loop Variables
    
    !Read file varibales
    integer::ioerror=0                                                  !Open file return value
    integer::stat                                                       !Read file data return value
    integer::status                                                     !Allocate memory status

    real(kind=8)::current_time                                          !Current time
    character(len=20)::char_time                                        !Time character type
    real(kind=8)::Velocity_Module                                       !Velocity Module

    integer::file_index                                                 !File index
    character(len=20)::char_file_index                                  !Character file index
    character(len=50)::file_name                                        !File name
    
    !Variables for folder opeartion
    Logical::dirExists
    !==========================================================================================================

    !Body of Initialize_SPH_Calculation

    !Call MPI functions
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )

    !**********************************************************************************************************
    !Initial particle information array
    particle_position=0.0d0
    particle_mass=0.0d0
    particle_rho=0.0d0
    particle_smooth_lengh=0.0d0
    particle_energy=0.0d0
    particle_c=0.0d0
    particle_press=0.0d0

    particle_type=0
    free_surface_type=0
    particle_initial_type=0
    Boundary_particle_type=0

    !**********************************************************************************************************
    
    !---------------------------------------------------------------------------------
    !Subdomain results directory
    inquire(directory='Subdomain', exist=dirExists )                                      
    
    if(dirExists) then
      
      if (Current_Processor_ID==Main_Processor) then
         call system ('rm -r Subdomain')                 !Remove the Subdomain folder first
      end if

      call system ('mkdir -p Subdomain')                 !Make the new folder

    else

      call system ('mkdir -p Subdomain') 

    endif
    !---------------------------------------------------------------------------------

    !Only the main processor output the data
    if (Current_Processor_ID==Main_Processor) then

      !---------------------------------------------------------------------------------
      !For Linux
      !*Create result folder
      !*fusion is for merging data of different time steps in one file
      !---------------------------------------------------------------------------------
      !Sampling directory
      inquire(directory='Sampling', exist=dirExists )                                      
      
      if(dirExists) then

        call system ('rm -r Sampling')                 !Remove the Sampling folder first
        call system ('mkdir -p Sampling/Press')        !Make the new folder
        call system ('mkdir -p Sampling/WaveHeight')

      else

        call system ('mkdir -p Sampling/Press')
        call system ('mkdir -p Sampling/WaveHeight')

      endif
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      !Paraview directory
      inquire(directory='Paraview', exist=dirExists )                                      
      
      if(dirExists) then

        call system ('rm -r Paraview')                 !Remove the Paraview folder first
        call system ('mkdir -p Paraview')              !Make the new folder

      else

        call system ('mkdir -p Paraview') 

      endif
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      !Binary Tecplot directory
      inquire(directory='Binary_Tecplot', exist=dirExists )                                      
      
      if(dirExists) then

        call system ('rm -r Binary_Tecplot')                 !Remove the Binary_Tecplot folder first
        call system ('mkdir -p Binary_Tecplot')              !Make the new folder

      else

        call system ('mkdir -p Binary_Tecplot') 

      endif
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      !Binary Tecplot directory
      inquire(directory='Binary_Paraview', exist=dirExists )                                      
      
      if(dirExists) then

        call system ('rm -r Binary_Paraview')                 !Remove the Binary_Paraview folder first
        call system ('mkdir -p Binary_Paraview')              !Make the new folder

      else

        call system ('mkdir -p Binary_Paraview') 

      endif
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      !Binary Tecplot directory
      inquire(directory='Binary_Tecplot_Grid', exist=dirExists )                                      
      
      if(dirExists) then

        call system ('rm -r Binary_Tecplot_Grid')                 !Remove the Binary_Tecplot_Grid folder first
        call system ('mkdir -p Binary_Tecplot_Grid')              !Make the new folder

      else

        call system ('mkdir -p Binary_Tecplot_Grid') 

      endif
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      !Binary Tecplot directory
      inquire(directory='Binary_Paraview_Grid', exist=dirExists )                                      
      
      if(dirExists) then

        call system ('rm -r Binary_Paraview_Grid')                 !Remove the Binary_Paraview_Grid folder first
        call system ('mkdir -p Binary_Paraview_Grid')              !Make the new folder

      else

        call system ('mkdir -p Binary_Paraview_Grid') 

      endif
      !---------------------------------------------------------------------------------

      ! call system ('touch output/0.00000/fsInfo.dat')
      ! call system ('touch output/0.00000/wallInfo.dat')
      ! call system ('touch output/0.00000/dummyInfo.dat')

      ! !---------------------------------------------------------------------------------
      ! !For Windows
      ! !*Create result folder
      ! !*fusion is for merging data of different time steps in one file
      ! inquire(directory='Sampling', exist=dirExists )                                      
      ! !---------------------------------------------------------------------------------
      ! !Sampling directory
      ! if(dirExists) then

      !   call system ('rmdir Sampling')                 !Remove the Sampling folder first
      !   call system ('mkdir Sampling/Press')        !Make the new folder
      !   call system ('mkdir Sampling/WaveHeight')

      ! else

      !   call system ('mkdir Sampling/Press')
      !   call system ('mkdir Sampling/WaveHeight')

      ! endif
      ! !---------------------------------------------------------------------------------

      ! !---------------------------------------------------------------------------------
      ! !Paraview directory
      ! inquire(directory='Paraview', exist=dirExists )                                      
      
      ! if(dirExists) then

      !   call system ('rmdir Paraview')                 !Remove the Paraview folder first
      !   call system ('mkdir Paraview')        

      ! else

      !   call system ('mkdir Paraview') 

      ! endif
      ! !---------------------------------------------------------------------------------

      ! !---------------------------------------------------------------------------------
      ! !Subdomain results directory
      ! inquire(directory='Subdomain', exist=dirExists )                                      
      
      ! if(dirExists) then

      !   call system ('rmdir Subdomain')                 !Remove the Subdomain folder first
      !   call system ('mkdir Subdomain')              !Make the new folder

      ! else

      !   call system ('mkdir Subdomain') 

      ! endif
      ! !---------------------------------------------------------------------------------
    endif

    !**********************************************************************************************************

    !**********************************************************************************************************
    !Distribute the partciles
    if(trim(adjustl(startFrom))=='InitialTime') then

        call Distribute_Initial_Particle()                   ! Satrt from initial step
         
    else if(trim(adjustl(startFrom))=='LatestTime') then
        
        call Distribute_Latest_Particle()                    ! Satrt from Latest step
        
    else
        
        write(*,*) " The 'startFrom' variable is not right. ('InitialTime' or 'LatestTime')"
        
    end if
    !**********************************************************************************************************

    !**********************************************************************************************************
    !Distribute the mesh for post-proceeding
    call Distribute_Grid_For_Postproceeding()

    !**********************************************************************************************************

    !**********************************************************************************************************
    !Distribute Backgroud Grid
    !Input: chain_x_number,chain_y_number,chain_z_number
    !Output: chain_near_mesh_name
    call Distribute_Backgroud_Grid(dim,&                            !Dimension of grid
                                   chain_x_number,&                 !Grid number in x direction
                                   chain_y_number,&                 !Grid number in y direction
                                   chain_z_number,&                 !Grid number in z direction
                                   chain_max_number,&               !Max grid number (grid_x_number*grid_y_number*grid_z_number)
                                   mesh_center_positon,&            !Grid center positon
                                   chain_near_mesh_number,&         !Near mesh number
                                   chain_near_mesh_name&            !grid near mesh name
                                   )
    
    !***********************************************************************************************************
    !Distribute the postion of fixed boundary particle
    if (DistributeFixedBoundaryParticleOrNot==1) then
        call Distribute_Fixed_Boundary_Particle()
    endif

    !***********************************************************************************************************

    !***********************************************************************************************************
    !initialized the Variables for pressure Sampling
    !Define the sampling points and wave probe position
    !For 2D simulation, we just need (x,y) position for pressure Sampling and (x) for wave probe 
    !For 3D simulation, we just need (x,y,z) position for pressure Sampling and (x,y) for wave probe 
    
    if (dim==2) then
        
        !------------------------------------------------------------------------------------------------------
        if (Sampling_Point_number>0) then
            
            !Pressure Sampling position
            ! Sampling_Point_Position(1,2)=0.05
            ! Sampling_Point_Position(1,1)=incline_x+Sampling_Point_Position(1,2)/tan(theta_rad)

            ! Sampling_Point_Position(2,2)=0.15
            ! Sampling_Point_Position(2,1)=incline_x+Sampling_Point_Position(2,2)/tan(theta_rad)
            
            ! Sampling_Point_Position(3,2)=0.25
            ! Sampling_Point_Position(3,1)=incline_x+Sampling_Point_Position(3,2)/tan(theta_rad)
            
            ! Sampling_Point_Position(4,2)=0.35
            ! Sampling_Point_Position(4,1)=incline_x+Sampling_Point_Position(4,2)/tan(theta_rad)

            Sampling_Point_Position(1,2)=0.003d0
            Sampling_Point_Position(1,1)=1.61d0

            Sampling_Point_Position(2,2)=0.015d0
            Sampling_Point_Position(2,1)=1.61d0
            
            Sampling_Point_Position(3,2)=0.030d0
            Sampling_Point_Position(3,1)=1.61d0
            
            Sampling_Point_Position(4,2)=0.080d0
            Sampling_Point_Position(4,1)=1.61d0

        endif
        !------------------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------------------
        !Wave probe position
        if (Wave_Probe_number>0) then

            Wave_Probe_Position(1,1)=0.3d0
            Wave_Probe_Position(2,1)=0.865d0
            Wave_Probe_Position(3,1)=1.114d0
            Wave_Probe_Position(4,1)=1.3625d0
            
        endif
        !------------------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------------------
        !Sampling Line position
        SamplingLine:if (SamplingLineNumber>0) then

            !----First Line----
            SamplingLineInitialPoint(1,1) = 0.0d0
            SamplingLineInitialPoint(1,2) = 0.25d0
            
            SamplingLineEndPoint(1,1) = 1.0d0
            SamplingLineEndPoint(1,2) = 0.25d0

            !----Second Line----
            SamplingLineInitialPoint(2,1) = 0.0d0
            SamplingLineInitialPoint(2,2) = 0.50d0

            SamplingLineEndPoint(2,1) = 1.0d0
            SamplingLineEndPoint(2,2) = 0.50d0

            !----Third Line----
            SamplingLineInitialPoint(3,1) = 0.0d0
            SamplingLineInitialPoint(3,2) = 0.75d0

            SamplingLineEndPoint(3,1) = 1.0d0
            SamplingLineEndPoint(3,2) = 0.75d0

            !Calculate the intervals
            do i=1,SamplingLineNumber
                do j=1,dim
                   SamplingLinePointInterval(i,j)=(SamplingLineEndPoint(i,j)-SamplingLineInitialPoint(i,j))/(SamplingLinePointNumber-1)
                enddo
            enddo

            !Calculate the sampling points on sampling line
            do k=1,SamplingLineNumber
                do i=1,SamplingLinePointNumber
                   do j=1,dim
                      SamplingLinePointPosition(k,i,j)=SamplingLineInitialPoint(k,j)+(i-1)*SamplingLinePointInterval(k,j)
                   enddo
                enddo
            enddo
            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            !Output the Sampling points on Sampling line
            if (Current_Processor_ID==Main_Processor) then
                
                do k=1,SamplingLineNumber

                    !transfer the sampling line index from integer to character
                    write(char_file_index,'(I2)') k

                    file_name="./Sampling/"//trim(adjustl(char_file_index))//"th_SamplingLine_Position.dat"

                    !Open the saving file
                    open(unit=1,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    
                    
                    if(ioerror==0) then
                       
                        write(*,'(A)') " The open action of the "// trim(adjustl(file_name))//" is successful!"

                        !Tecplot header
                        write(1,*) "TITLE='DISTRIBUTION'"
                        write(1,*) "VARIABLES= 'X' 'Y'"

                        do i=1,SamplingLinePointNumber
                            write(1,'(6f12.6)') (SamplingLinePointPosition(k,i,j),j=1,dim)
                        end do

                    else
                        write(*,'(A)') " The open action of the "// trim(adjustl(file_name))//" is fail!"
                    end if
                   
                    close(1)

                enddo

            endif
            !--------------------------------------------------------------------------------------------------
            
        endif SamplingLine
        !------------------------------------------------------------------------------------------------------

    elseif(dim==3) then


    endif
    !***********************************************************************************************************

    !***********************************************************************************************************
    !calculate the CFL conditions
    !Morris (1997)
    CFL_Morris_dt=0.125*smooth_length**2/(1.01E-6)
    Min_CFL_dt(1)=CFL_Morris_dt

    !Monaghan (1989,1992)
    CFL_Monaghan_dt=0.25*smooth_length/(1.1*sqrt(square_c_0))
    Min_CFL_dt(2)=CFL_Monaghan_dt

    !minimum CFL dt (0.1 is the coefficient)
    if (leapforg_or_eular==3) then
        CFL_dt=2*0.1*minval(Min_CFL_dt)                               ! For Rugger-kuta ,time step can be a little larger                            
    else
        CFL_dt=0.1*minval(Min_CFL_dt)                                 ! For others ,keep same
    endif
    
    !The main Processor Output the initial particle information
    if (Current_Processor_ID==Main_Processor) then

        write(*,*) "==============================================================="
        write(*,*) "Recommand 'dt' for iteration: ",CFL_dt
        write(*,*) "==============================================================="

    endif
    !***********************************************************************************************************

    !***********************************************************************************************************

    !定义粒子动力粘性系数
    do i=1,ture_total_particle_number
        if(abs(particle_type(i))==1) then
            particle_eta(i)=0.0
        else if(abs(particle_type(i))/=1) then
            particle_eta(i)=1.0e-3
        end if
    end do

    !-----------------------------------------------------------------------------------------------------------
    !Calculate the initial kinetic energy calculation

    Initial_total_kinetic_energy=0.0d0

    do i=1,particle_ture_number

       Initial_total_kinetic_energy=Initial_total_kinetic_energy+0.5*particle_mass(i)*DOT_PRODUCT(particle_velocity(i,:),particle_velocity(i,:))
    
    enddo
    !-----------------------------------------------------------------------------------------------------------

    !***********************************************************************************************************
    !The main Processor Output the initial particle information
    if (Current_Processor_ID==Main_Processor) then

        !Open the Output file
        open(unit=5,file="Particle_Information.dat",status="replace",position="rewind",action="write",iostat=ioerror)    
        if(ioerror==0) then
            write(*,*) "The open action of the Particle_Information.dat is Fine!"
            
            !input tecplot header
            write(5,*) "TITLE='DISTRIBUTION'"
            write(5,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
            !output fluid particle
            write(5,*) "ZONE I=",particle_ture_number," F=POINT"
            do i=1,particle_ture_number
                 write(5,100) (particle_position(i,j)/domain_size_y,j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    100          format(8F20.10) 
            end do
            
            !output wave maker particle
            if (make_wave_or_not==1) then
                write(5,*) "ZONE I=",wave_maker_particle_number," F=POINT"
                do k=1,wave_maker_particle_number
                    i=k+particle_ture_number
                    write(5,100) (particle_position(i,j)/domain_size_y,j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
                end do
            endif

            !output fix ghost boundary particle
            if (DistributeFixedBoundaryParticleOrNot==1) then
                write(5,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
                do k=1,fix_ghost_particle_number
                    i=k+particle_ture_number+wave_maker_particle_number
                    write(5,100) (particle_position(i,j)/domain_size_y,j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
                end do
            endif

        else
            write(*,*) "The open action of the Initial_particle_information.dat is fail!"
        end if
        close(5)

    endif
    !***********************************************************************************************************

    !***********************************************************************************************************
    !Output the analytical results

    if (Current_Processor_ID==Main_Processor) then

        do k=0,5

            current_time=k*0.1d0

            !transfer the Current time from real to character
            write(char_time,'(f6.3)') current_time

            file_name=trim(adjustl(char_time))//"_Analytical_Results.dat"

            !Calculate the analytical results
            do i=1,total_element_node_number

                element_node_velocity(i,2)=-U_Module*exp(-8.0*PI**2/Reynolds_Number*current_time)*cos(2*PI*element_node_position(i,1))*sin(2*PI*element_node_position(i,2))
                element_node_velocity(i,1)= U_Module*exp(-8.0*PI**2/Reynolds_Number*current_time)*sin(2*PI*element_node_position(i,1))*cos(2*PI*element_node_position(i,2))

                element_node_press(i)=0.25*water_rho_0*U_Module**2*exp(-16.0*PI**2/Reynolds_Number*current_time)*(cos(4*PI*element_node_position(i,1))+cos(4*PI*element_node_position(i,2)))
            
            end do

            !Open the saving file
            open(unit=1,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    
            
            if(ioerror==0) then
               
                write(*,*) "The open action of the "// trim(adjustl(file_name))//" is successful!"

                !Tecplot header
                write(1,*) "TITLE='DISTRIBUTION'"
                write(1,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'P' 'U'"

                write(1,*) 'ZONE T=" '//trim(adjustl(char_time))//' ",'//"N=",total_element_node_number,",E=",total_element_number,",DATAPACKING=POINT,ZONETYPE=FEQUADRILATERAL "   
                
                !域名，节点数，单元数，数据是在点上还是块上（数据是点还是块），单元类型
                !单元类型：FETRIANGLE(三角形网格)、FEQUADRILATERAL(四节点四边形网格)、FETETRAHEDRON(四节点四面体网格)、FEBRICK(8节点六面体)
                !对于混合网格需要不同的ZONES
                
                do i=1,total_element_node_number
                    Velocity_Module=sqrt(element_node_velocity(i,1)**2+element_node_velocity(i,2)**2)
                    write(1,100) (element_node_position(i,j)/domain_size_y,j=1,dim),(element_node_velocity(i,j)/velocity_0,j=1,dim),element_node_press(i)/pressure_0,Velocity_Module/velocity_0
                end do
                
                do i=1,total_element_number
                    write(1,'(8I8)') (element_node_index(i,j),j=1,NodeNumberInOneElement)
                end do
                

            else
                write(*,*) "The open action of the "// trim(adjustl(file_name))//" is fail!"
            end if
           
            close(1)
            !-------------------------------------------------------------------------------

        enddo

    endif
    !***********************************************************************************************************
    
end subroutine Initialize_SPH_Calculation