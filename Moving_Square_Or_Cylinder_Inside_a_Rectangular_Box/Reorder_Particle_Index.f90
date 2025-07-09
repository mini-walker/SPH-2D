!**************************************************************************************************************
!
!  SUBROUTINE : Reorder_Particle_Index
!
!  PURPOSE    : (1) Reorder the particle index : Remove the ill particle and the oulet particle out the computational domain
!               (2) Renew particle number of each particle type
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

subroutine  Reorder_Particle_Index(i_time_step)

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the superior-subroutien
    integer,intent(in)::i_time_step                                     ! Current time step 

    ! Variables in program
    integer::i,j,k,L,m                                                  ! Variables for looping
    integer::Temp_Total_Particle_Number                                 ! Temp total particle number

    ! Variables for file opeartion
    integer::File_index                                                 ! File index
    character(len=100)::File_name                                       ! File name

    !==========================================================================================================


    ! Body of subroutine Reorder_Particle_Index
   
    !==========================================================================================================
    ! We should identify the ill particle after we update the in/outlet particle
    ! Check the fluid particles in the domain or not

    do i=1,Particle_ture_number
        
        if( particle_type(i)==Water_Particle_Label ) then        ! Only the water particle
        
            if( particle_position(i,1)<=-1.0E-10 .or. particle_position(i,1)>boundary_size_x .or. particle_position(i,dim)<=-1.0E-10 .or. particle_position(i,dim)>boundary_size_y ) then
                
                particle_type(i)=Ill_Particle_Label

                ! if ( Current_Processor_ID==Main_Processor ) then
                !     write(*,'(10F10.4)') ( Particle_position(i,j),j=1,dim )
                ! endif

            endif

        endif
            
    enddo
    !==========================================================================================================
   





    !==========================================================================================================

    !----------------------------------------------------------------------------------------------------------
    ! Attention: The variable, old 'Ture_total_particle_number', is the sum of (1) Particle_ture_number;
    !                                                                          (2) Wave_maker_particle_number;
    !                                                                          (3) Fix_ghost_particle_number;
    !                                                                          (4) InletOutlet_Boundary_Particle_Number.
    !                                                                          (5) New_Fix_ghost_particle_number.
    !----------------------------------------------------------------------------------------------------------

    Ture_total_particle_number=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number+New_Fix_ghost_particle_number
    
    Temp_Total_Particle_Number=Ture_total_particle_number+New_Inlet_Particle_Number

    Particle_ture_number       = 0
    Fix_ghost_particle_number  = 0
    Wave_maker_particle_number = 0
    Inlet_Particle_Number      = 0
    Outlet_Particle_Number     = 0

    do i=1,Temp_Total_Particle_Number

        ! Fluid particle
        if (particle_type(i)==Water_Particle_Label) then

            Particle_ture_number=Particle_ture_number+1
            Fluid_Particle_Index(Particle_ture_number)=i

        ! Fix ghost particle number
        elseif (particle_type(i)==Fix_Ghost_Particle_Label) then

            Fix_ghost_particle_number=Fix_ghost_particle_number+1
            Fix_Ghost_Particle_Index(Fix_ghost_particle_number)=i

        ! Inlet Particle Number
        elseif (particle_type(i)==Inlet_Particle_Label) then

            Inlet_Particle_Number=Inlet_Particle_Number+1
            Inlet_Particle_Index(Inlet_Particle_Number)=i    

        ! Outlet Particle Number
        elseif (particle_type(i)==Outlet_Particle_Label) then

            Outlet_Particle_Number=Outlet_Particle_Number+1
            Outlet_Particle_Index(Outlet_Particle_Number)=i

        ! Wavemaker particle number
        elseif (particle_type(i)==Wave_Maker_Particle_Label) then

            Wave_maker_particle_number=Wave_maker_particle_number+1
            Wave_Maker_Particle_Index(Wave_maker_particle_number)=i

        endif
    
    enddo

    !----------------------------------------------------------------------------------------------------------
    ! Attention: The variable, old 'Ture_total_particle_number', is the sum of (1) Particle_ture_number;
    !                                                                          (2) Wave_maker_particle_number;
    !                                                                          (3) Fix_ghost_particle_number;
    !                                                                          (4) InletOutlet_Boundary_Particle_Number.
    !----------------------------------------------------------------------------------------------------------

    ! Renew the in/outlet particle number and new Ture_total_particle_number
    InletOutlet_Boundary_Particle_Number=Inlet_Particle_Number+Outlet_Particle_Number

    Ture_total_particle_number=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number
    
    !==========================================================================================================





    !==========================================================================================================

    !----------------------------------------------------------------------------------------------------------
    ! Re-order the particle index, Delete the error particle and the particle out the domain
    ! New order: 1, Fluid particle; 2, Wave maker; 3, fixed boundary; 4, Period boundary; 5, In/outlet boundary
    k=0
    do L=1,Particle_ture_number

        k=k+1
        i=Fluid_Particle_Index(L)

        Buffer_particle_position(k,:)=particle_position(i,:)
        Buffer_particle_velocity(k,:)=particle_velocity(i,:)
        Buffer_particle_mass(k)=particle_mass(i)                          
        Buffer_particle_type(k)=particle_type(i)                                 
        Buffer_particle_initial_type(k)=particle_initial_type(i)                         
        Buffer_particle_rho(k)=particle_rho(i)                                                 
        Buffer_particle_smooth_lengh(k)=particle_smooth_lengh(i)                    
        Buffer_particle_c(k)=particle_c(i)                               
        Buffer_particle_press(k)=particle_press(i)                           
        Buffer_particle_energy(k)=particle_energy(i)                                                     
        Buffer_free_surface_type(k)=free_surface_type(i)                        
        Buffer_Boundary_particle_type(k)=Boundary_particle_type(i)

        Buffer_particle_smoke_value(k)=particle_smoke_value(i)
        Buffer_particle_division_degree(k)=particle_division_degree(i)

    enddo

    do L=1,Wave_maker_particle_number

        k=k+1
        i=Wave_Maker_Particle_Index(L)

        Buffer_particle_position(k,:)=particle_position(i,:)
        Buffer_particle_velocity(k,:)=particle_velocity(i,:)
        Buffer_particle_mass(k)=particle_mass(i)                          
        Buffer_particle_type(k)=particle_type(i)                                 
        Buffer_particle_initial_type(k)=particle_initial_type(i)                         
        Buffer_particle_rho(k)=particle_rho(i)                                                 
        Buffer_particle_smooth_lengh(k)=particle_smooth_lengh(i)                    
        Buffer_particle_c(k)=particle_c(i)                               
        Buffer_particle_press(k)=particle_press(i)                           
        Buffer_particle_energy(k)=particle_energy(i)                                                    
        Buffer_free_surface_type(k)=free_surface_type(i)                        
        Buffer_Boundary_particle_type(k)=Boundary_particle_type(i)

        Buffer_particle_smoke_value(k)=particle_smoke_value(i)
        Buffer_particle_division_degree(k)=particle_division_degree(i)

    enddo

    do L=1,Fix_ghost_particle_number

        k=k+1
        i=Fix_Ghost_Particle_Index(L)

        Buffer_particle_position(k,:)=particle_position(i,:)
        Buffer_particle_velocity(k,:)=particle_velocity(i,:)
        Buffer_particle_mass(k)=particle_mass(i)                          
        Buffer_particle_type(k)=particle_type(i)                                 
        Buffer_particle_initial_type(k)=particle_initial_type(i)                         
        Buffer_particle_rho(k)=particle_rho(i)                                                 
        Buffer_particle_smooth_lengh(k)=particle_smooth_lengh(i)                    
        Buffer_particle_c(k)=particle_c(i)                               
        Buffer_particle_press(k)=particle_press(i)                           
        Buffer_particle_energy(k)=particle_energy(i)                                                      
        Buffer_free_surface_type(k)=free_surface_type(i)                        
        Buffer_Boundary_particle_type(k)=Boundary_particle_type(i)

        Buffer_particle_smoke_value(k)=particle_smoke_value(i)
        Buffer_particle_division_degree(k)=particle_division_degree(i)

    enddo

    do L=1,Inlet_Particle_Number

        k=k+1
        i=Inlet_Particle_Index(L)

        Buffer_particle_position(k,:)=particle_position(i,:)
        Buffer_particle_velocity(k,:)=particle_velocity(i,:)
        Buffer_particle_mass(k)=particle_mass(i)                          
        Buffer_particle_type(k)=particle_type(i)                                 
        Buffer_particle_initial_type(k)=particle_initial_type(i)                         
        Buffer_particle_rho(k)=particle_rho(i)                                                 
        Buffer_particle_smooth_lengh(k)=particle_smooth_lengh(i)                    
        Buffer_particle_c(k)=particle_c(i)                               
        Buffer_particle_press(k)=particle_press(i)                           
        Buffer_particle_energy(k)=particle_energy(i)                                                     
        Buffer_free_surface_type(k)=free_surface_type(i)                        
        Buffer_Boundary_particle_type(k)=Boundary_particle_type(i)

        Buffer_particle_smoke_value(k)=particle_smoke_value(i)
        Buffer_particle_division_degree(k)=particle_division_degree(i)

    enddo

    do L=1,Outlet_Particle_Number

        k=k+1
        i=Outlet_Particle_Index(L)

        Buffer_particle_position(k,:)=particle_position(i,:)
        Buffer_particle_velocity(k,:)=particle_velocity(i,:)
        Buffer_particle_mass(k)=particle_mass(i)                          
        Buffer_particle_type(k)=particle_type(i)                                 
        Buffer_particle_initial_type(k)=particle_initial_type(i)                         
        Buffer_particle_rho(k)=particle_rho(i)                                                 
        Buffer_particle_smooth_lengh(k)=particle_smooth_lengh(i)                    
        Buffer_particle_c(k)=particle_c(i)                               
        Buffer_particle_press(k)=particle_press(i)                           
        Buffer_particle_energy(k)=particle_energy(i)                                                       
        Buffer_free_surface_type(k)=free_surface_type(i)                        
        Buffer_Boundary_particle_type(k)=Boundary_particle_type(i)

        Buffer_particle_smoke_value(k)=particle_smoke_value(i)
        Buffer_particle_division_degree(k)=particle_division_degree(i)

    enddo

    particle_position=Buffer_particle_position
    particle_velocity=Buffer_particle_velocity
    particle_mass=Buffer_particle_mass                        
    particle_type=Buffer_particle_type                                 
    particle_initial_type=Buffer_particle_initial_type                      
    particle_rho=Buffer_particle_rho                                                
    particle_smooth_lengh=Buffer_particle_smooth_lengh                   
    particle_c=Buffer_particle_c                             
    particle_press=Buffer_particle_press                          
    particle_energy=Buffer_particle_energy                                                    
    free_surface_type=Buffer_free_surface_type                       
    Boundary_particle_type=Buffer_Boundary_particle_type

    particle_smoke_value=Buffer_particle_smoke_value
    particle_division_degree=Buffer_particle_division_degree
    !----------------------------------------------------------------------------------------------------------

    !==========================================================================================================
    






!     !**********************************************************************************************************
!     if ( i_time_step==Actual_start_time_step .and. Current_Processor_ID==Main_Processor ) then

!         ! Open the Output file
!         File_index = General_File_Port
!         File_name  = './Initial_Data/Check_Particle_Reorder_Result.dat'

!         call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
    
!         if( IOERROR==0 ) then

!             ! Input tecplot header
!             write(File_index,*) 'TITLE="DISTRIBUTION"'   
!             write(File_index,*) 'VARIABLES= "X" "Y" "VX" "VY" "RHO" "P" '
            
!             ! Output fluid particle
!             write(File_index,*) "ZONE I=",Particle_ture_number," F=POINT"
!             do i=1,Particle_ture_number
!                  write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
! 100              format(8F10.4) 
!             enddo
            
!             ! Output wave maker particle
!             if ( Make_wave_or_not==1 ) then
!                 write(File_index,*) "ZONE I=",Wave_maker_particle_number," F=POINT"
!                 do k=1,Wave_maker_particle_number
!                     i=k+Particle_ture_number
!                     write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
!                 enddo
!             endif

!             ! Output fix ghost boundary particle
!             if ( DistributeFixedBoundaryParticleOrNot==1 ) then
!                 write(File_index,*) "ZONE I=",Fix_ghost_particle_number," F=POINT"
!                 do k=1,Fix_ghost_particle_number
!                     i=k+Particle_ture_number+Wave_maker_particle_number
!                     write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
!                 enddo
!             endif

!             ! Inlet and outlet boundary particles
!             if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then

!                 write(File_index,*) "ZONE I=",InletOutlet_Boundary_Particle_Number," F=POINT"
!                 do k=1,InletOutlet_Boundary_Particle_Number
!                     i=k+Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number
!                     write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
!                 enddo
                
!             endif

!         endif

!         close(File_index)

!     endif
!     !**********************************************************************************************************


end subroutine  Reorder_Particle_Index
