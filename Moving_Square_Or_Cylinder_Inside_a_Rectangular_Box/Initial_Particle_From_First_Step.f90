!**************************************************************************************************************
!  Subroutine : Initial_Particle_From_First_Step
!
!  PURPOSE    : Set the initial information from initial domain data for SPH calculation
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
!  Input      : domain size and packing information
!
!  Output     : all particles information
!
!  Note       : MPI version: mpich-3.2
!               Fortran version: Inter Fortran (ifort) 18.0.0
!**************************************************************************************************************

subroutine Initial_Particle_From_First_Step()

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,L,m                                                  ! Variables for looping
    integer::m_virtual,n_virtual                                        ! Virtual particle order
    
    ! Read file varibales
    integer::File_index                                                 ! File index
    character(len=100)::File_name                                       ! File name ( The file name character lenghth should be 100 )
        
    ! Particle packing varibales
    real(kind=8)::x_0,y_0,x_1,y_1                                       ! Starting and edn point on angular bisector
    real(kind=8)::fix_particle_distance                                 ! Real fix particle distance
    real(kind=8)::position_x,position_y,position_z                      ! Partcile temp coordiantes
    real(kind=8)::line_distance                                         ! Line distance
    real(kind=8)::particle_x,particle_y

    !==========================================================================================================
    


    !==========================================================================================================
    ! Body of Initial_Particle_From_First_Step


    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! Deviding the initial particle position 
    Fluid_Particle_Initialization: if (Deviding_from_input_files_or_not==0) then

        !-------------------------Get the particle position data from Tecplot input file-----------------------
        if ( Initial_file_type==0 ) then

            call Initial_Fluid_Particle_From_Tecplot()

        !-------------------------Get the particle position data from Gambit input file------------------------
        elseif ( Initial_file_type==1 ) then

            call Initial_Fluid_Particle_From_Gambit()

        else

            if (Current_Processor_ID==Main_Processor) then

                write(*,*) " The input file type for field initialization is not right! "
                write(*,*) " The type should be '0' or '1' ."
                
            endif

        endif
        !------------------------------------------------------------------------------------------------------

    else

        !--------------------------Get the particle position data by regular distribution----------------------
        
        call Initial_Fluid_Particle_From_Regular_Distribution()

    endif Fluid_Particle_Initialization
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! Particle volume and mass ( The fluid particle number(Particle_ture_number) is defined in the Initial fluid particle subroutine in Subroutine_Module )
    !write(*,*) Particle_ture_number    
    !intial_volume_0=tank_square/Particle_ture_number           
    intial_volume_0=interior_dx*interior_dy                                             ! Intial volume
    
    ! Attentation: the initial water mass is used during the calculation
    ! I cost one night to find this bug, donot forget!
    water_particle_mass_0=water_rho_0*intial_volume_0

    
    !write(*,*) Particle_ture_number
    
    ! Initializate other information for fluid particles
    do i=1,Particle_ture_number
        particle_mass(i)=water_particle_mass_0                                          ! Particle mass
        particle_type(i)=Water_Particle_Label                                           ! Particle type(1-air,2-water,-1-boudary air,-2-boudary water)
        particle_initial_type(i)=Water_Particle_Label                                   ! Particle initial type
        particle_smooth_lengh(i)=smooth_length                                          ! Particle smooth lengh
        particle_c(i)=sqrt(square_c_0)                                                  ! Particle c
        free_surface_type(i)=0                                                          ! Free surface type(0-no,1-yes)
        particle_energy(i)=e_0                                                          ! Particle energy
        particle_volume(i)=intial_volume_0                                              ! Particle volume
    end do
    !***********************************************************************************************************
    

    !***************************************Generate wave maker particles***************************************
    if ( Fix_wave_maker_or_not/=0 .and. Make_wave_or_not==1 ) then
        
        !Set the wave maker height
        if(trim(adjustl(wave_maker_type))=='WaveTheory') then
            wave_maker_height=domain_size_y
        else
            wave_maker_height=boundary_size_y 
        endif 

        !-------------------------------------------------------------------------------------------------------
        n_virtual=0                                                     ! Virtual particle index
        wave_maker_particle_layer=0
        
        do i=1,interior_particle_number_y+200                           ! 1 to m row
            do j=1,wave_maker_layer                                     ! 1 to m column
 
                position_x=-(j-1)*interior_dx-interior_dx/2.0
                position_y=(i-1)*interior_dy+interior_dy/2.0
                
                ! make sure the hight of the wave maker boundary
                if(position_y<=wave_maker_height) then
                
                    n_virtual=n_virtual+1
                    k=n_virtual+Particle_ture_number

                    particle_position(k,1)=position_x             
                    particle_position(k,2)=position_y             
            
                    !particle_press(k)=water_rho_0*g*(domain_size_y-particle_position(k,2))           
                    !particle_rho(k)=water_rho_0*((1+7*g*(domain_size_y-particle_position(k,2))/square_c_0)**(0.142857142)) 

                    particle_rho(k)   = water_rho_0*(-6.0*g*(particle_position(k,2)-domain_size_y)/square_c_0+1.0)**(1.0/6.0)
                    particle_press(k) = square_c_0*water_rho_0*((particle_rho(k)/water_rho_0)**7.0-1.0)/7.0
                    
                    wave_maker_particle_layer(k)=j
                
                end if

            end do
        end do
        !-------------------------------------------------------------------------------------------------------
        
        
        !wave maker particle number
        Wave_maker_particle_number=n_virtual
        
        !initialization other information for wave maker particles
        do i=1,Wave_maker_particle_number
            m_virtual=i+Particle_ture_number                                               ! Virtual particle index
            particle_velocity(m_virtual,:)=0.0                                             ! Particle velocity
            particle_mass(m_virtual)=intial_volume_0*particle_rho(i)                       ! Particle mass
            particle_type(m_virtual)=Wave_Maker_Particle_Label                             ! Particle type(130-wave maker particle)
            particle_initial_type(m_virtual)=Wave_Maker_Particle_Label                     ! Particle initial type
            particle_smooth_lengh(m_virtual)=smooth_length                                 ! Particle smooth lengh
            particle_c(m_virtual)=sqrt(square_c_0)                                         ! Particle c
            particle_energy(m_virtual)=e_0                                                 ! Particle energy
            particle_volume(m_virtual)=intial_volume_0                                     ! Particle volume
       end do
    
    endif
    !***********************************************************************************************************

    
    !*********************************Distribute the fixed boundary particles***********************************
    Fix_ghost_particle_number=0

    if ( DistributeFixedBoundaryParticleOrNot==1 ) then
       
        n_virtual=Wave_maker_particle_number                                               ! Virtual particle index


        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        !-------------------------------------------------------------------------------------------------------
        ! Bottom fixed boundary particles
        External_Bottom_Boundary: do i=1,fix_ghost_layer                                   ! Fixed ghost layer
            
            x_0 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_0 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            x_1 = boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_1 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            
            ! Define the initial number in first row
            if (i==1) then
               j=0
               do 
                  j=j+1                                                                    ! Virtual particle number
                  position_x=(j-1)*fix_ghost_dx+x_0                                        ! Virtual particle x coordinate
                  position_y=y_0                                                           ! Virtual particle y coordinate
                
                  ! Exceed the angular bisector yes or no
                  if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit
                         
               end do
               k=j                       ! The virtual fixed boundary particle number of first row 
               
            end if
            
            ! Generate all Bottom particles for each row
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                    ! Has K points and (k-1) interval
            
            do j=1,k-1                   ! From 1 to k-1 points,as the last points will be repeated in next interval
                
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Particle_ture_number                                   ! Virtual particle index
                particle_position(m_virtual,1)=(j-1)*fix_particle_distance+x_0             ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0                                         ! Virtual particle y coordinate
                free_surface_type(m_virtual)=0
                
                particle_type(m_virtual)=Fix_Ghost_Particle_Label                          ! Fixed Ghost Particle Label
                particle_initial_type(m_virtual)=Fix_Ghost_Particle_Label                  ! Fixed Ghost Particle Label

                ! Form 0 to boundary_size_x is Horizontial mirror rule. others are Center mirror
                if ( 0.0d0 <= particle_position(m_virtual,1) .and. particle_position(m_virtual,1) <= boundary_size_x ) then
                    Boundary_particle_type(m_virtual)=-6                                   ! Horizontial type
                else
                    Boundary_particle_type(m_virtual)=-7                                   ! Center type
                endif
                
            end do
            
            k=k+2                        ! The virtual particle number of next row              
            
        end do External_Bottom_Boundary
        !-------------------------------------------------------------------------------------------------------


        !-------------------------------------------------------------------------------------------------------
        ! Bottom fixed boundary particles
        External_Right_Boundary: do i=1,fix_ghost_layer                                    ! Fixed ghost layer
            
            x_0 = boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_0 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            x_1 = boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_1 = boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            
            ! Define the initial number in first row
            if (i==1) then
               j=0
               do 
                  j=j+1                                                                    ! Virtual particle number
                  position_x=x_0                                                           ! Virtual particle x coordinate
                  position_y=y_0+(j-1)*fix_ghost_dx                                        ! Virtual particle y coordinate
                
                  ! Exceed the angular bisector yes or no
                  if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit
                         
               end do
               k=j                       ! The virtual fixed boundary particle number of first row 
               
            end if
            
            ! Generate all Bottom particles for each row
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                    ! Has K points and (k-1) interval
            
            do j=1,k-1                   ! From 1 to k-1 points,as the last points will be repeated in next interval
                
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Particle_ture_number                                   ! Virtual particle index
                particle_position(m_virtual,1)=x_0                                         ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0+(j-1)*fix_particle_distance             ! Virtual particle y coordinate
                free_surface_type(m_virtual)=0
                
                particle_type(m_virtual)=Fix_Ghost_Particle_Label                          ! Fixed Ghost Particle Label
                particle_initial_type(m_virtual)=Fix_Ghost_Particle_Label                  ! Fixed Ghost Particle Label

                ! Form 0 to boundary_size_x is Vertical mirror rule. others are Center mirror
                if ( 0.0d0 <= particle_position(m_virtual,2) .and. particle_position(m_virtual,2) <= boundary_size_y ) then
                    Boundary_particle_type(m_virtual)=-8                                   ! Vertical type
                else
                    Boundary_particle_type(m_virtual)=-7                                   ! Center type
                endif
                
            end do
            
            k=k+2                        ! The virtual particle number of next row              
            
        end do External_Right_Boundary
        !-------------------------------------------------------------------------------------------------------


        !-------------------------------------------------------------------------------------------------------
        ! Top fixed boundary particles (From right to left)
        External_Top_Boundary: do i=1,fix_ghost_layer                                      ! Fixed ghost layer
            
            x_0 = boundary_size_x+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_0 = boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            x_1 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_1 = boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0

            ! write(*,*) x_1,y_1
            ! Define the initial number in first row
            if (i==1) then
               j=0
               do 
                   j=j+1                                                                   ! Virtual particle index
                   position_x=x_0-(j-1)*fix_ghost_dx                                       ! Virtual particle x coordinate  
                   position_y=y_0                                                          ! Virtual particle y coordinate
                   
                   ! Define the initial number in first row
                   if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit

                end do
                k=j                        ! The virtual particle number of first row 
                
            end if
            
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                    ! k points and k-1 intervals
            
            do j=1,k-1                     ! 1 to k-1 points (Delte the end point)
                
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Particle_ture_number                                   ! Virtual particle index

                particle_position(m_virtual,1)=x_0-(j-1)*fix_particle_distance             ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0                                         ! Virtual particle y coordinate
                free_surface_type(m_virtual)=0

                particle_type(m_virtual)=Fix_Ghost_Particle_Label                          ! Fixed ghost particle Label
                particle_initial_type(m_virtual)=Fix_Ghost_Particle_Label                  ! Fixed ghost pParticle Label

                ! Form 0 to boundary_size_x is Horizontial mirror rule. others are Center mirror
                if ( 0.0d0 <= particle_position(m_virtual,1) .and. particle_position(m_virtual,1) <= boundary_size_x ) then
                    Boundary_particle_type(m_virtual)=-6                                   ! Horizontial type
                else
                    Boundary_particle_type(m_virtual)=-7                                   ! Center type
                endif

            end do
            
            k=k+2                         ! The virtual particle number of first row                

        end do External_Top_Boundary
        !-------------------------------------------------------------------------------------------------------
        


        !-------------------------------------------------------------------------------------------------------
        ! Top fixed boundary particles (From right to left)
        External_Left_Boundary: do i=1,fix_ghost_layer                                      ! Fixed ghost layer
            
            x_0 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_0 = boundary_size_y+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            x_1 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_1 =                -(i-1)*fix_ghost_dy-fix_ghost_dy/2.0

            ! write(*,*) x_1,y_1
            ! Define the initial number in first row
            if (i==1) then

               j=0
               do 
                   j=j+1                                                                   ! Virtual particle index
                   position_x=x_0                                                          ! Virtual particle x coordinate  
                   position_y=y_0-(j-1)*fix_ghost_dx                                       ! Virtual particle y coordinate
                   
                   ! Define the initial number in first row
                   if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit

                end do
                k=j                        ! The virtual particle number of first row 
                
            end if
            
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                    ! k points and k-1 intervals
            
            do j=1,k-1                     ! 1 to k-1 points (Delte the end point)
                
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Particle_ture_number                                   ! Virtual particle index

                particle_position(m_virtual,1)=x_0                                         ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0-(j-1)*fix_particle_distance             ! Virtual particle y coordinate
                
                free_surface_type(m_virtual)=0

                particle_type(m_virtual)=Fix_Ghost_Particle_Label                          ! Fixed ghost particle Label
                particle_initial_type(m_virtual)=Fix_Ghost_Particle_Label                  ! Fixed ghost pParticle Label

                ! Form 0 to boundary_size_x is Horizontial mirror rule. others are Center mirror
                if ( 0.0d0 <= particle_position(m_virtual,2) .and. particle_position(m_virtual,2) <= boundary_size_y ) then
                    Boundary_particle_type(m_virtual)=-8                                   ! Horizontial type
                else
                    Boundary_particle_type(m_virtual)=-7                                   ! Center type
                endif

            end do
            
            k=k+2                         ! The virtual particle number of first row                

        end do External_Left_Boundary
        !-------------------------------------------------------------------------------------------------------


        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




        !-------------------------------------------------------------------------------------------------------
        ! Block fixed boundary particles
        Block_Bottom: do i=1,fix_ghost_layer                                               ! Fixed ghost layer
            
            x_0=(Block_Center_X-0.5*Block_Length)+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_0=(Block_Center_Y-0.5*Block_Length)+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0

            x_1=(Block_Center_X+0.5*Block_Length)-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_1=(Block_Center_Y-0.5*Block_Length)+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0

            ! Define the initial number in first row
            if (i==1) then
               j=0
               do 
                  j=j+1                                                                    ! Virtual particle number
                  position_x=(j-1)*fix_ghost_dx+x_0                                        ! Virtual particle x coordinate
                  position_y=y_0                                                           ! Virtual particle y coordinate
                
                  ! Exceed the angular bisector yes or no
                  if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit
                         
               end do
               k=j                        !the virtual fixed boundary particle number of first row 
               
            end if
            
            ! Generate all Bottom particles for each row
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                      ! Has K points and (k-1) interval
            
            do j=1,k-1                                     !From 1 to k-1 points,as the last points will be repeated in next interval
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Particle_ture_number                                     ! Virtual particle index
                particle_position(m_virtual,1)=x_0+(j-1)*fix_particle_distance               ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0                                           ! Virtual particle y coordinate
                free_surface_type(m_virtual)=0

                particle_type(m_virtual)=Fix_Ghost_Particle_Label                            ! Fixed Ghost Particle Label
                particle_initial_type(m_virtual)=Body_Particle_Label                         ! Body Particle Label

                if (j==1) then
                    Boundary_particle_type(m_virtual)=-7                                     ! Center type
                else
                    Boundary_particle_type(m_virtual)=-6                                     ! Horizontal type
                endif

            end do
            
            k=k-2                        !the virtual particle number of next row              
            
        enddo Block_Bottom

        Block_Right: do i=1,fix_ghost_layer                                                  ! Fixed ghost layer
            
            x_0=(Block_Center_X+0.5*Block_Length)-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_0=(Block_Center_Y-0.5*Block_Length)+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0

            x_1=(Block_Center_X+0.5*Block_Length)-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_1=(Block_Center_Y+0.5*Block_Length)-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0

            ! Define the initial number in first row
            if (i==1) then
               j=0
               do 
                  j=j+1                                                                      ! Virtual particle number
                  position_x=x_0                                                             ! Virtual particle x coordinate
                  position_y=y_0+(j-1)*fix_ghost_dx                                          ! Virtual particle y coordinate
                
                  ! Exceed the angular bisector yes or no
                  if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit
                         
               end do
               k=j                        !the virtual fixed boundary particle number of first row 
               
            end if
            
            ! Generate all Bottom particles for each row
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                      ! Has K points and (k-1) interval
            
            do j=1,k-1                                     !From 1 to k-1 points,as the last points will be repeated in next interval
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Particle_ture_number                                     ! Virtual particle index
                particle_position(m_virtual,1)=x_0                                           ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0+(j-1)*fix_particle_distance               ! Virtual particle y coordinate
                free_surface_type(m_virtual)=0

                particle_type(m_virtual)=Fix_Ghost_Particle_Label                            ! Fix Ghost Particle Label
                particle_initial_type(m_virtual)=Body_Particle_Label                         ! Body Particle Label

                if (j==1) then
                    Boundary_particle_type(m_virtual)=-7                                     ! Center type
                else
                    Boundary_particle_type(m_virtual)=-8                                     ! Vertical type
                endif

            end do
            
            k=k-2                        !the virtual particle number of next row              
            
        end do Block_Right

        Block_Top: do i=1,fix_ghost_layer                                                     ! Fixed ghost layer
            
            x_0=(Block_Center_X+0.5*Block_Length)-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0
            y_0=(Block_Center_Y+0.5*Block_Length)-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0

            x_1=(Block_Center_X-0.5*Block_Length)+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_1=(Block_Center_Y+0.5*Block_Length)-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0

            ! Define the initial number in first row
            if (i==1) then
               j=0
               do 
                  j=j+1                                                                       ! Virtual particle number
                  position_x=x_0-(j-1)*fix_ghost_dx                                           ! viVirtualrtual particle x coordinate
                  position_y=y_0                                                              ! Virtual particle y coordinate
                
                  ! Exceed the angular bisector yes or no
                  if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit
                         
               end do
               k=j                        ! The virtual fixed boundary particle number of first row 
               
            end if
            
            ! Generate all Bottom particles for each row
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                      ! Has K points and (k-1) interval
            
            do j=1,k-1                                     ! From 1 to k-1 points,as the last points will be repeated in next interval
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Particle_ture_number                                     ! Virtual particle index
                particle_position(m_virtual,1)=x_0-(j-1)*fix_particle_distance               ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0                                           ! Virtual particle y coordinate
                free_surface_type(m_virtual)=0

                particle_type(m_virtual)=Fix_Ghost_Particle_Label                            ! Fixed ghost particle label
                particle_initial_type(m_virtual)=Body_Particle_Label                         ! Body Particle Label

                if (j==1) then
                    Boundary_particle_type(m_virtual)=-7                                     ! Center type
                else
                    Boundary_particle_type(m_virtual)=-6                                     ! Horizontal type
                endif

            end do
            
            k=k-2                        ! The virtual particle number of next row              
            
        enddo Block_Top

        Block_Left: do i=1,fix_ghost_layer                                                   ! Fixed ghost layer
            
            x_0=(Block_Center_X-0.5*Block_Length)+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_0=(Block_Center_Y+0.5*Block_Length)-(i-1)*fix_ghost_dy-fix_ghost_dy/2.0

            x_1=(Block_Center_X-0.5*Block_Length)+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0
            y_1=(Block_Center_Y-0.5*Block_Length)+(i-1)*fix_ghost_dy+fix_ghost_dy/2.0

            ! Define the initial number in first row
            if (i==1) then
               j=0
               do 
                  j=j+1                                                                      ! Virtual particle number
                  position_x=x_0                                                             ! Virtual particle x coordinate
                  position_y=y_0-(j-1)*fix_ghost_dx                                          ! Virtual particle y coordinate
                
                  ! Exceed the angular bisector yes or no
                  if(sqrt((x_1-position_x)**2+(y_1-position_y)**2)<=fix_ghost_dx/2.0) exit
                         
               end do
               k=j                        ! The virtual fixed boundary particle number of first row 
               
            end if
            
            ! Generate all Bottom particles for each row
            fix_particle_distance=sqrt((x_1-x_0)**2+(y_1-y_0)**2)/(k-1)                      ! Has K points and (k-1) interval
            
            do j=1,k-1      ! From 1 to k-1 points,as the last points will be repeated in next interval
                n_virtual=n_virtual+1
                m_virtual=n_virtual+Particle_ture_number                                     ! Virtual particle index
                particle_position(m_virtual,1)=x_0                                           ! Virtual particle x coordinate 
                particle_position(m_virtual,2)=y_0-(j-1)*fix_particle_distance               ! Virtual particle y coordinate
                free_surface_type(m_virtual)=0
                
                particle_type(m_virtual)=Fix_Ghost_Particle_Label                            ! Fixed Ghost Particle Label
                particle_initial_type(m_virtual)=Body_Particle_Label                         ! Body Particle Label

                if (j==1) then
                    Boundary_particle_type(m_virtual)=-7                                     ! Center type
                else
                    Boundary_particle_type(m_virtual)=-8                                     ! Vertical type
                endif

            end do
            
            k=k-2                        ! The virtual particle number of next row              
            
        enddo Block_Left
        !-------------------------------------------------------------------------------------------------------

        Fix_ghost_particle_number=n_virtual-Wave_maker_particle_number
        
        ! Initialize other information for fixed ghost boundary particles
        do i=1,Fix_ghost_particle_number
            m_virtual=i+Particle_ture_number+Wave_maker_particle_number                      ! Virtual particle index
             
            particle_velocity(m_virtual,:)=0.0                                               ! Particle velocity
            particle_rho(m_virtual)=water_rho_0                                              ! Particle rho
            particle_mass(m_virtual)=water_particle_mass_0                                   ! Particle mass
            particle_smooth_lengh(m_virtual)=smooth_length                                   ! Particle smooth lengh
            particle_press(m_virtual)=0.0d0+Background_Pressure                              ! Particle press
            particle_c(m_virtual)=sqrt(square_c_0)                                           ! Particle c
            particle_energy(m_virtual)=e_0                                                   ! Particle energy
            particle_volume(m_virtual)=intial_volume_0                                       ! Particle volume

        end do

    endif
    !***********************************************************************************************************


    !***********************************************************************************************************
    ! Distribute the initial inlet and outlet particle (Just do the ghost opeartion)
    Inlet_Particle_Number=0
    Outlet_Particle_Number=0
    InletOutlet_Boundary_Particle_Number=0

    if (DistributeInletOutletBoundaryParticleOrNot==1) then

        if ( Distribute_Inlet_Boundary_X==1 ) then

            do i=1,Particle_ture_number

                !---------------------------------------------------------------------------------------------
                ! Move right outside to the left inside and left outside to the right inside

                if( particle_position(i,1)-Inlet_Boundary_X <= InletOutlet_Boundary_Zone_Width ) then

                    Inlet_Particle_Number=Inlet_Particle_Number+1
                    InletOutlet_Boundary_Particle_Number=InletOutlet_Boundary_Particle_Number+1

                    k=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number

                    particle_position(k,1)=2*Inlet_Boundary_X-particle_position(i,1)
                    particle_position(k,2)=particle_position(i,2)

                    particle_type(k)=Inlet_Particle_Label                                        ! Particle type(-3-fixed ghost boundary particles)
                    particle_initial_type(k)=Inlet_Particle_Label                                ! Particle initial type

                endif
                !---------------------------------------------------------------------------------------------

            enddo
            
        endif


        if ( Distribute_Outlet_Boundary_X==1 ) then

            do i=1,Particle_ture_number

                !---------------------------------------------------------------------------------------------
                ! Move right outside to the left inside and left outside to the right inside

                if( Outlet_Boundary_X-particle_position(i,1) <= InletOutlet_Boundary_Zone_Width ) then

                    Outlet_Particle_Number=Outlet_Particle_Number+1
                    InletOutlet_Boundary_Particle_Number=InletOutlet_Boundary_Particle_Number+1

                    k=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number

                    particle_position(k,1)=2*Outlet_Boundary_X-particle_position(i,1)
                    particle_position(k,2)=particle_position(i,2)

                    particle_type(k)=Outlet_Particle_Label                                       ! Particle type(-3-fixed ghost boundary particles)
                    particle_initial_type(k)=Outlet_Particle_Label                               ! Particle initial type

                endif
                !---------------------------------------------------------------------------------------------

            enddo
            
        endif


        !initialize other information for fixed ghost boundary particles
        do i=1,InletOutlet_Boundary_Particle_Number

            m_virtual=i+Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number    ! Virtual particle index
            
            particle_velocity(m_virtual,:)=0.0                                                       ! Particle velocity
            particle_rho(m_virtual)=water_rho_0                                                      ! Particle rho
            particle_mass(m_virtual)=water_particle_mass_0                                           ! Particle mass
            particle_smooth_lengh(m_virtual)=smooth_length                                           ! Particle smooth lengh
            particle_press(m_virtual)=0.0d0+Background_Pressure                                      ! Particle press
            particle_c(m_virtual)=sqrt(square_c_0)                                                   ! Particle c
            particle_energy(m_virtual)=e_0                                                           ! Particle energy
            particle_volume(m_virtual)=intial_volume_0                                               ! Particle volume

        enddo

    endif

    !***********************************************************************************************************
    ! Synchronize all processors calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
    

    !***********************************************************************************************************
    ! Actual start time step and ture total particle number 
    Actual_start_time_step=1
    Ture_total_particle_number=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number
    !***********************************************************************************************************
























    !***********************************************************************************************************
    ! Checking the particle initialization results
    if (Current_Processor_ID==Main_Processor) then

        ! Open the Output file
        File_index = General_File_Port
        File_name  = "./Initial_Data/Initialization_Particle_Checking.dat"
        call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
  
        if(IOERROR==0) then

            ! Input tecplot header
            write(File_index,*) "TITLE='DISTRIBUTION'"
            write(File_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
            ! Output fluid particle
            write(File_index,*) "ZONE I=",Particle_ture_number," F=POINT"
            do i=1,Particle_ture_number
                 write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
100              format(8F10.4) 
            end do
            
            ! Output wave maker particle
            if (make_wave_or_not==1) then
                write(File_index,*) "ZONE I=",Wave_maker_particle_number," F=POINT"
                do k=1,Wave_maker_particle_number
                    i=k+Particle_ture_number
                    write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                end do
            endif

            ! Output fix ghost boundary particle
            if (DistributeFixedBoundaryParticleOrNot==1) then
                write(File_index,*) "ZONE I=",Fix_ghost_particle_number," F=POINT"
                do k=1,Fix_ghost_particle_number
                    i=k+Particle_ture_number+Wave_maker_particle_number
                    write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                end do
            endif

            ! Inlet and outlet boundary particles
            if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then

                write(File_index,*) "ZONE I=",InletOutlet_Boundary_Particle_Number," F=POINT"
                do k=1,InletOutlet_Boundary_Particle_Number
                    i=k+Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number
                    write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
                end do
                
            endif
            
            
        endif
        
        close(File_index)

    endif
    !***********************************************************************************************************


    ! Synchronize all processors calculation
    call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)


    
end subroutine Initial_Particle_From_First_Step