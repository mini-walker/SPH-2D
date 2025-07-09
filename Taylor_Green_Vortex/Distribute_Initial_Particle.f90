!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  Subroutine: Distribute_Initial_Particle
!
!  PURPOSE: Set the initial information from initial domain data for SPH calculation
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
!  Output: all particles information
!
!  Note: MPI version: mpich-3.2
!        Fortran version: Inter Fortran (ifort) 18.0.0
!****************************************************************************

subroutine Distribute_Initial_Particle()

    use information_module
    use Public_variable_module
    use MPI
    
    implicit none

    !Variables in subroutine
    
    !==========================================================================================================
    integer::i,j,k,L,m                                                  !loop Variables
    real(kind=8)::temp_press                                            !partcile temp press
    real(kind=8)::C_k                                                   !coefficients
    real(kind=8)::temp_C_k                                              !temp coefficients
    integer::m_virtual,n_virtual                                        !virtual particle order
    
    !Read file varibales
    integer::ioerror=0                                                  !open file return value
    integer::stat                                                       !read file data return value
    integer::status                                                     !allocate memory status
    integer::skip_lines                                                 !The skip lines number
    character(len=20)::temp_character                                   !temp characters for reading
    integer::file_index                                                 !file index 
    
    !particle packing varibales
    real(kind=4),dimension(dim)::position_difference                    !position difference
    real(kind=8)::x_0,y_0,x_1,y_1                                       !starting and edn point on angular bisector
    real(kind=8)::fix_particle_distance                                 !real fix particle distance
    real(kind=8)::position_x,position_y,position_z                      !partcile temp coordiantes
    real(kind=8)::line_distance                                         !line distance
    real(kind=8)::particle_x,particle_y

    !Variables for gambit initialization 
    real(kind=8),allocatable,dimension(:,:)::node                       !node coordinate(node index, x/y/z)
    integer,allocatable,dimension(:,:)::element                         !nodes in gambit element(node index, 1/2/3)
    integer::node_number,element_number                                 !node and element number
    integer::temp_int                                                   !temp integer value
    real(kind=8)::position                                              !position value

    !==========================================================================================================
    
    !Body of Set_initial_SPH_Field

    !Call MPI functions
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )


    !**********************************************************************************************************
    !Deviding the initial particle position 
    if (deviding_from_input_files_or_not==0) then

        !-------------------------Get the particle position data from Tecplot input file-----------------------
        if (trim(adjustl(initial_file_type))=='Tecplot') then

            !open the input Tecplot data files
            file_index=1+Current_Processor_ID
            open(unit=file_index,file='pretreat_particle.txt',status="old",position="rewind",action="read",iostat=ioerror) 

            !Skip the initial three head titles
            skip_lines=3
            do i=1,skip_lines
               read(file_index,*) temp_character
            end do
            
            k=0
            do
                k=k+1                                                     !calculator plus 1
                read(file_index,*,iostat=status) (particle_position(k,j),j=1,dim)
                if(status/=0) exit
                
                !particle_press(k)=rho_0*g*(y-particle_position(k,2))                                 
                !particle_rho(k)=rho_0*((1+7*g*(y-particle_position(k,2))/square_c_0)**(0.142857142)) 

                particle_rho(k)=water_rho_0*(-6.0*g*(particle_position(k,2)-domain_size_y)/square_c_0+1.0)**(1.0/6.0)
                particle_press(k)=square_c_0*water_rho_0*((particle_rho(k)/water_rho_0)**7.0-1.0)/7.0
  
            end do

            particle_ture_number=k-1                                                      !actual partcile number
            !write(*,*) particle_ture_number              
            intial_volume_0=tank_square/k                                                 !intial volume                          

            water_particle_mass_0=water_rho_0*intial_volume_0

            !close the file index
            close(file_index)
        !------------------------------------------------------------------------------------------------------

        !-------------------------Get the particle position data from Gambit input file------------------------
        else if(trim(adjustl(initial_file_type))=='Gambit') then

            !open the input Gambit data files
            file_index=1+Current_Processor_ID
            open(unit=file_index,file='particle.neu',status="old",position="rewind",action="read",iostat=ioerror)  

            !Skip the initial five head titles
            skip_lines=5
            do i=1,skip_lines
               read(file_index,*) temp_character
            end do

            !Read node and element number
            read(file_index,*) node_number,element_number

            !Allocate memory for node, element, and center points
            allocate(node(node_number,dim))                                                              !node coordinate
            allocate(element(element_number,gambit_element_node_number))                                 !element point index

            !Skip two head titles
            skip_lines=2
            do i=1,skip_lines
               read(file_index,*) temp_character
            end do
            
            !Read node coordinates
            do i=1,node_number
                read(file_index,*) temp_int,(node(i,j),j=1,dim)
                !write(3,*) i,(node(i,j),j=1,dim)
            end do
            
            !Skip two head titles
            skip_lines=2
            do i=1,skip_lines
               read(file_index,*) temp_character
            end do

            !Read the node index in each element
            do i=1,element_number
                read(file_index,*) temp_int,temp_int,temp_int,(element(i,j),j=1,gambit_element_node_number)
                !write(3,*) i,(element(i,j),j=1,gambit_element_node_number)
            end do

            !Calculate the center points for initialization

            do i=1,element_number                             !The actual partcile number is equal to element number
                do j=1,dim
                    position=0.0d0
                    do k=1,gambit_element_node_number
                        position=position+node(element(i,k),j)
                    end do
                    particle_position(i,j)=position/3.0
                end do    
            end do

            !Assign the density and press
            do k=1,element_number

                !particle_press(k)=water_rho_0*g*(domain_size_y-particle_position(k,2))           
                !particle_rho(k)=water_rho_0*((1+7*g*(domain_size_y-particle_position(k,2))/square_c_0)**(0.142857142)) 

                particle_rho(k)=water_rho_0*(-6.0*g*(particle_position(k,2)-domain_size_y)/square_c_0+1.0)**(1.0/6.0)
                particle_press(k)=square_c_0*water_rho_0*((particle_rho(k)/water_rho_0)**7.0-1.0)/7.0

            end do

            particle_ture_number=k                                                        !actual partcile number
            !write(*,*) particle_ture_number              
            intial_volume_0=tank_square/k                                                 !intial volume

            water_particle_mass_0=water_rho_0*intial_volume_0

        !-------------------------------------------------------------------------------------------------------

        else

            write(*,*) " The input file type for field initialization is not right!"
            write(*,*) " The type should be 'Tecplot' or 'Gambit' ."

        endif
        !--------------------------------------------------------------------------------------------------------


    else

        !--------------------------Get the particle position data by regular distribution------------------------
        k=0

        do i=1,interior_particle_number_y                !1-m row
            do j=1,interior_particle_number_x            !1-m column
                
                k=k+1
                particle_position(k,1)=(j-1)*interior_dx+interior_dx/2.0             
                particle_position(k,2)=(i-1)*interior_dy+interior_dy/2.0             
                
                !particle velocity
                particle_velocity(k,1)= U_Module*sin(2*PI*particle_position(k,1))*cos(2*PI*particle_position(k,2))  

                particle_velocity(k,2)=-U_Module*cos(2*PI*particle_position(k,1))*sin(2*PI*particle_position(k,2))

                !particle pressure
                particle_press(k)=0.25*water_rho_0*U_Module**2*(cos(4*PI*particle_position(k,1))+cos(4*PI*particle_position(k,2)))

                !particle density
                particle_rho(k)=water_rho_0*(1+7*particle_press(k)/(square_c_0*water_rho_0))**(0.142857142) 
                
                 
            end do
        end do
        
        particle_ture_number=k                                                        !actual partcile number
        !write(*,*) particle_ture_number    
        !intial_volume_0=tank_square/particle_ture_number           
        intial_volume_0=interior_dx*interior_dy                                       !intial volume

        !-------------------------------------------------------------------------------------------------------
        
    endif

    
    !write(*,*) particle_ture_number
    
    !initialization other information for fluid particles
    do i=1,particle_ture_number
        particle_mass(i)=intial_volume_0*particle_rho(i)                                !particle mass
        particle_type(i)=2                                                              !particle type(1-air,2-water,-1-boudary air,-2-boudary water)
        particle_initial_type(i)=2                                                      !particle initial type
        particle_smooth_lengh(i)=smooth_length                                          !particle smooth lengh
        particle_c(i)=sqrt(square_c_0)                                                  !particle c
        free_surface_type(i)=0                                                          !free surface type(0-no,1-yes)
        particle_energy(i)=e_0                                                          !particle energy
    end do
    !***********************************************************************************************************
    

    !***************************************Generate wave maker particles***************************************
    if (fix_wave_maker_or_not/=0 .and. make_wave_or_not==1) then
        
        !Set the wave height
        if(trim(adjustl(wave_maker_type))=='WaveTheory') then
            wave_maker_height=domain_size_y
        else
            wave_maker_height=boundary_size_y 
        endif 

        !-------------------------------------------------------------------------------------------------------
        n_virtual=0                                                     !virtual particle index
        wave_maker_particle_layer=0
        
        do i=1,interior_particle_number_y+200                           !1 to m row
            do j=1,wave_maker_layer                                     !1 to m column
 
                position_x=-(j-1)*interior_dx-interior_dx/2.0
                position_y=(i-1)*interior_dy+interior_dy/2.0
                
                ! make sure the hight of the wave maker boundary
                if(position_y<=wave_maker_height) then
                
                    n_virtual=n_virtual+1
                    k=n_virtual+particle_ture_number

                    particle_position(k,1)=position_x             
                    particle_position(k,2)=position_y             
            
                    !particle_press(k)=water_rho_0*g*(domain_size_y-particle_position(k,2))           
                    !particle_rho(k)=water_rho_0*((1+7*g*(domain_size_y-particle_position(k,2))/square_c_0)**(0.142857142)) 

                    particle_rho(k)=water_rho_0*(-6.0*g*(particle_position(k,2)-domain_size_y)/square_c_0+1.0)**(1.0/6.0)
                    particle_press(k)=square_c_0*water_rho_0*((particle_rho(k)/water_rho_0)**7.0-1.0)/7.0
                    
                    wave_maker_particle_layer(k)=j
                
                end if

            end do
        end do
        !-------------------------------------------------------------------------------------------------------
        
        
        !wave maker particle number
        wave_maker_particle_number=n_virtual
        
        !initialization other information for wave maker particles
        do i=1,wave_maker_particle_number
            m_virtual=i+particle_ture_number                                               !virtual particle index
            particle_velocity(m_virtual,:)=0.0                                             !particle velocity
            particle_mass(m_virtual)=intial_volume_0*particle_rho(i)                       !particle mass
            particle_type(m_virtual)=130                                                   !particle type(130-wave maker particle)
            particle_initial_type(m_virtual)=130                                           !particle initial type
            particle_smooth_lengh(m_virtual)=smooth_length                                 !particle smooth lengh
            particle_c(m_virtual)=sqrt(square_c_0)                                         !particle c
            particle_energy(m_virtual)=e_0                                                 !particle energy
       end do
    
    end if
    !***********************************************************************************************************

    
    !*********************************Distribute the fixed boundary particles***********************************
    if (DistributeFixedBoundaryParticleOrNot==1) then
       
        n_virtual=wave_maker_particle_number                               !virtual particle index
        
        !-----------------------------------------------------------------------------------------------------------
        !Top fixed boundary particles (From right to left)
        Top: do i=1,fix_ghost_layer                                            !fix ghost layer
            
            x_0=boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_0=boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            x_1=-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_1=boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0

            
            !write(*,*) x_1,y_1
            !Define the initial number in first row
            if (i==1) then
               j=0
               do 
                   j=j+1                                                                        !virtual particle index
                   position_x=x_0-(j-1)*fix_ghost_dx                                            !virtual particle x coordinate  
                   position_y=y_0                                                               !virtual particle y coordinate
                   
                   !Define the initial number in first row
                   if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit

                end do
                k=j                        !the virtual particle number of first row 
                
            end if
            
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)               !k points and k-1 intervals
            
            do j=1,k-1                                                                !1 to k-1 points (Delte the end point)
                n_virtual=n_virtual+1
                m_virtual=n_virtual+particle_ture_number                              !virtual particle index

                particle_position(m_virtual,1)=x_0-(j-1)*fix_particle_distance        !virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0                                    !virtual particle y coordinate
                free_surface_type(m_virtual)=0
            end do
            
            k=k+2                         !the virtual particle number of first row                

        end do Top
        !-----------------------------------------------------------------------------------------------------------
        
        !-----------------------------------------------------------------------------------------------------------
        !Bottom fixed boundary particles
        Bottom: do i=1,fix_ghost_layer                                             !fix ghost layer
            x_0=-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_0=-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            x_1=boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_1=-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            
            !Define the initial number in first row
            if (i==1) then
               j=0
               do 
                  j=j+1                                                                      !virtual particle number
                  position_x=(j-1)*fix_ghost_dx+x_0                                          !virtual particle x coordinate
                  position_y=y_0                                                             !virtual particle y coordinate
                
                  !exceed the angular bisector yes or no
                  if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit
                         
               end do
               k=j                        !the virtual fixed boundary particle number of first row 
               
            end if
            
            !Generate all Bottom particles for each row
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                    !Has K points and (k-1) interval
            
            do j=1,k-1                                     !From 1 to k-1 points,as the last points will be repeated in next interval
                n_virtual=n_virtual+1
                m_virtual=n_virtual+particle_ture_number                                   !virtual particle index
                particle_position(m_virtual,1)=(j-1)*fix_particle_distance+x_0             !virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0                                         !virtual particle y coordinate
                free_surface_type(m_virtual)=0
            end do
            k=k+2                        !the virtual particle number of next row              
            
        end do Bottom
        !-----------------------------------------------------------------------------------------------------------


        !-----------------------------------------------------------------------------------------------------------
        !Right fixed boundary particles
        Right: do i=1,fix_ghost_layer                                            !fix ghost layer
            x_0=boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_0=-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            x_1=boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_1=boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            
            !write(*,*) x_1,y_1
            !Define the initial number in first row
            if (i==1) then
               j=0
               do 
                   j=j+1                                                                        !virtual particle index
                   position_x=x_0                                                               !virtual particle x coordinate  
                   position_y=y_0+(j-1)*fix_ghost_dx                                            !virtual particle y coordinate
                   
                   !Define the initial number in first row
                   if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit

                end do
                k=j                        !the virtual particle number of first row 
                
            end if
            
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)
            do j=1,k-1
                n_virtual=n_virtual+1
                m_virtual=n_virtual+particle_ture_number                                          !virtual particle index
                particle_position(m_virtual,1)=x_0                                                !virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0+(j-1)*fix_particle_distance                    !virtual particle y coordinate
                free_surface_type(m_virtual)=0
            end do
            k=k+2                         !the virtual particle number of first row                

        end do Right
        !-----------------------------------------------------------------------------------------------------------


        !-----------------------------------------------------------------------------------------------------------
        !Left fixed boundary particles (From right to left)
        Left: do i=1,fix_ghost_layer                                            !fix ghost layer

            x_0=-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_0=boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            x_1=-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_1=-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            
            !write(*,*) x_1,y_1
            !Define the initial number in first row
            if (i==1) then
               j=0
               do 
                   j=j+1                                                                        !virtual particle index
                   position_x=x_0                                                               !virtual particle x coordinate  
                   position_y=y_0-(j-1)*fix_ghost_dy                                            !virtual particle y coordinate
                   
                   !Define the initial number in first row
                   if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit

                end do
                k=j                        !the virtual particle number of first row 
                
            end if
            
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)
            do j=1,k-1
                n_virtual=n_virtual+1
                m_virtual=n_virtual+particle_ture_number                              !virtual particle index
                particle_position(m_virtual,1)=x_0                                    !virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0-(j-1)*fix_particle_distance        !virtual particle y coordinate
                free_surface_type(m_virtual)=0
            end do
            k=k+2                         !the virtual particle number of first row                

        end do Left
        !-----------------------------------------------------------------------------------------------------------

        
        fix_ghost_particle_number=n_virtual-wave_maker_particle_number
        
        !initialization other information for fixed ghost boundary particles
        do i=1,fix_ghost_particle_number
            m_virtual=i+particle_ture_number+wave_maker_particle_number                  !virtual particle index
            particle_velocity(m_virtual,:)=0.0                                           !particle velocity
            particle_rho(m_virtual)=water_rho_0                                          !particle rho
            particle_mass(m_virtual)=water_particle_mass_0                               !particle mass
            particle_type(m_virtual)=-3                                                  !particle type(-3-fixed ghost boundary particles)
            particle_initial_type(m_virtual)=-3                                          !particle initial type
            particle_smooth_lengh(m_virtual)=smooth_length                               !particle smooth lengh
            particle_press(m_virtual)=0.0                                                !particle press
            particle_c(m_virtual)=sqrt(square_c_0)                                       !particle c
            particle_energy(m_virtual)=e_0                                               !particle energy
        end do

    endif
    
    !***********************************************************************************************************
    
    !-----------------------------------------------------------------------------------------------------------
    !This one is very important
    !As the intel complier cannot note the array out of the boundary in some conditions, so we should pay attention to the array memory
    !The n_total is very important, some conditions the boundary particles will be larger than the interior nodes, this may make errors
    !So we check the memory of the array
    if (Current_Processor_ID==Main_Processor) then
        if (particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number>n_total) then

            write(*,'(A)') " ====================================================================================="
            write(*,30) n_total,particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number
30          Format(" The total partcile number is (",I8,"),and this value should larger than (",I8,")")
            write(*,'(A)') " ====================================================================================="
            write(*,'(A)') " Errors may happen during the calculation, you should increase the memory for the arraies"
            write(*,'(A)') " ====================================================================================="
            write(*,20) 1.1*(particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number)/particle_ture_number
20          Format(" Recommand value of n_total_Factor is (",F6.3,")")
            write(*,'(A)') " ====================================================================================="
        
        else

            write(*,'(A)') " ====================================================================================="
            write(*,10) n_total,particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
10          Format(" The total partcile number is (",I8,"),and the actual particle number is(",I8,")")
            write(*,'(A)') " ====================================================================================="
            write(*,20) 1.1*(particle_ture_number+wave_maker_particle_number+2*fix_ghost_particle_number)/particle_ture_number
            write(*,'(A)') " ====================================================================================="
            write(*,'(A)') " The memory for the calculation is enough !"
            write(*,'(A)') " ====================================================================================="
            
        endif
    end if
    !-----------------------------------------------------------------------------------------------------------

    !***********************************************************************************************************
    !The main Processor Output the initial particle information
    if (Current_Processor_ID==Main_Processor) then

        !Open the Output file
        open(unit=5,file="Initial_particle_information.dat",status="replace",position="rewind",action="write",iostat=ioerror)    
        if(ioerror==0) then
            write(*,*) "The open action of the Initial_particle_information.dat is Fine!"
            
            !input tecplot header
            write(5,*) "TITLE='DISTRIBUTION'"
            write(5,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
            !output fluid particle
            write(5,*) "ZONE I=",particle_ture_number," F=POINT"
            do i=1,particle_ture_number
                 write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    100          format(8F20.10) 
            end do
            
            !output wave maker particle
            if (make_wave_or_not==1) then
                write(5,*) "ZONE I=",wave_maker_particle_number," F=POINT"
                do k=1,wave_maker_particle_number
                    i=k+particle_ture_number
                    write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
                end do
            endif

            !output fix ghost boundary particle
            if (DistributeFixedBoundaryParticleOrNot==1) then
                write(5,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
                do k=1,fix_ghost_particle_number
                    i=k+particle_ture_number+wave_maker_particle_number
                    write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
                end do
            endif
            
        else
            write(*,*) "The open action of the Initial_particle_information.dat is fail!"
        end if
        close(5)

    endif
    !***********************************************************************************************************

    !***********************************************************************************************************
    !Actual start time step and ture total particle number 
    actual_start_time_step=1
    ture_total_particle_number=particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number+Period_Boundary_Particle_Number
    !***********************************************************************************************************


    
end subroutine Distribute_Initial_Particle