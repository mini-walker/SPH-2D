!**************************************************************************************************************
!
!  SUBROUTINE : Distribute_Period_Boundary_Particle
!
!  PURPOSE    : Distribute period boundary particles on the desired boundary direction
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

subroutine  Distribute_Period_Boundary_Particle(i_time_step)

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the superior subroutien
    integer,intent(in)::i_time_step                                     ! Current time step 

    ! Variables in local subroutine
    integer::i,j,k,L,m                                                  ! Variables for loop
    
    ! Variables for file opeartion
    integer::File_index                                                 ! File index
    character(len=100)::File_name                                       ! File name

    integer::Temp_Period_Boundary_Particle_Number                       ! Temp Period Boundary Particle Number
    !==========================================================================================================



    ! Body of subroutine Distribute_Period_Boundary_Particle


    !----------------------------------------------------------------------------------------------------------
    ! Attention: The variable, 'ture_total_particle_number', is the sum of (1) Particle_ture_number;
    !                                                                      (2) Wave_maker_particle_number;
    !                                                                      (3) Fix_ghost_particle_number;
    !                                                                      (4) InletOutlet_Boundary_Particle_Number;
    !----------------------------------------------------------------------------------------------------------
    ture_total_particle_number=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number


    !==========================================================================================================
    ! Move the particle out the boundary
    ! Inside particle period motion
    if ( X_circulation==1 ) then

        do i=1,ture_total_particle_number

            if ( particle_type(i)==Ill_Particle_Label ) cycle            ! No error particles

            !-------------------------------------------------------------------------------------------------
            ! Move right outside to the left inside and left outside to the right inside

            if ( particle_position(i,1) > Period_Boundary_X_Max ) then

                particle_position(i,1)=particle_position(i,1)-Period_Boundary_Length_X
                particle_position(i,2)=particle_position(i,2)

            elseif ( particle_position(i,1) < Period_Boundary_X_Min ) then

                particle_position(i,1)=particle_position(i,1)+Period_Boundary_Length_X
                particle_position(i,2)=particle_position(i,2)

            endif
            !-------------------------------------------------------------------------------------------------

        enddo

    endif

    if ( Y_circulation==1 ) then

        do i=1,ture_total_particle_number

            if ( particle_type(i)==Ill_Particle_Label ) cycle             ! No error particles

            !-------------------------------------------------------------------------------------------------
            ! Move top outside to the bottom inside and bottom outside to the top inside

            if ( particle_position(i,2) > Period_Boundary_Y_Max ) then

                particle_position(i,1)=particle_position(i,1)
                particle_position(i,2)=particle_position(i,2)-Period_Boundary_Length_Y

            elseif (particle_position(i,2) < Period_Boundary_Y_Min) then

                particle_position(i,1)=particle_position(i,1)
                particle_position(i,2)=particle_position(i,2)+Period_Boundary_Length_Y

            endif

            !-------------------------------------------------------------------------------------------------

        enddo

    endif
    !=========================================================================================================




    !=========================================================================================================
    ! Period boundary in x direction

    Period_Boundary_Particle_Number=0

    if ( X_circulation==1 ) then

        ! Left and right Period
        do i=1,ture_total_particle_number

            if ( particle_type(i)==Ill_Particle_Label ) cycle        !No error particles
            
            !-------------------------------------------------------------------------------------------------
            ! Move left inside to the right outside
            if ( particle_position(i,1) <= ( Period_Boundary_X_Min+Period_Boundary_Zone_Width) .and. particle_position(i,1) >= Period_Boundary_X_Min ) then

                Period_Boundary_Particle_Number=Period_Boundary_Particle_Number+1            ! Increase Period_Boundary_Particle_Number
                j=ture_total_particle_number+Period_Boundary_Particle_Number                 ! Particle index

                particle_position(j,1)=particle_position(i,1)+Period_Boundary_Length_X
                particle_position(j,2)=particle_position(i,2)

                particle_velocity(j,:)=particle_velocity(i,:)

                particle_rho(j)=particle_rho(i)                                              ! Particle density
                particle_mass(j)=particle_mass(i)                                            ! Particle mass
                particle_type(j)=Period_Particle_Label                                       ! Particle type(Period_Particle_Label-Period boundary particles)
                particle_initial_type(j)=Period_Particle_Label                               ! Particle initial type
                particle_smooth_lengh(j)=particle_smooth_lengh(i)                            ! Particle smooth lengh
                particle_press(j)=particle_press(i)                                          ! Particle pressure
                particle_c(j)=particle_c(i)                                                  ! Particle c
                particle_energy(j)=particle_energy(i)                                        ! Particle energy
                particle_volume(j)=particle_volume(i)                                        ! Particle volume
                particle_smoke_value(j)=particle_smoke_value(i)                              ! Particle smoke line value

            endif
            !-------------------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------------------
            ! Move right inside to the left outside
            if ( particle_position(i,1) >= ( Period_Boundary_X_Max-Period_Boundary_Zone_Width) .and. particle_position(i,1) <= Period_Boundary_X_Max ) then

                Period_Boundary_Particle_Number=Period_Boundary_Particle_Number+1            ! Increase Period_Boundary_Particle_Number
                j=ture_total_particle_number+Period_Boundary_Particle_Number                 ! Particle index

                particle_position(j,1)=particle_position(i,1)-Period_Boundary_Length_X
                particle_position(j,2)=particle_position(i,2)

                particle_velocity(j,:)=particle_velocity(i,:)
                
                particle_rho(j)=particle_rho(i)                                              ! Particle density
                particle_mass(j)=particle_mass(i)                                            ! Particle mass
                particle_type(j)=Period_Particle_Label                                       ! Particle type(Period_Particle_Label-Period boundary particles)
                particle_initial_type(j)=Period_Particle_Label                               ! Particle initial type
                particle_smooth_lengh(j)=particle_smooth_lengh(i)                            ! Particle smooth lengh
                particle_press(j)=particle_press(i)                                          ! Particle pressure
                particle_c(j)=particle_c(i)                                                  ! Particle c
                particle_energy(j)=particle_energy(i)                                        ! Particle energy
                particle_volume(j)=particle_volume(i)                                        ! Particle volume
                particle_smoke_value(j)=particle_smoke_value(i)                              ! Particle smoke line value

            endif
            !-------------------------------------------------------------------------------------------------

        enddo

    endif
    !=========================================================================================================



    !=========================================================================================================
    if ( Y_circulation==1 ) then

        ! Have circulation in x direction or not
        if ( X_circulation==1 ) then

            Temp_Period_Boundary_Particle_Number=Period_Boundary_Particle_Number

        else

            Temp_Period_Boundary_Particle_Number=0
            
        endif

        ! Top and Bottom Period
        do i=1,ture_total_particle_number+Temp_Period_Boundary_Particle_Number

            if(particle_type(i)==Ill_Particle_Label) cycle        !No error particles
            
            !-------------------------------------------------------------------------------------------------
            ! Move bottom inside to the top outside
            if( particle_position(i,2) <= ( Period_Boundary_Y_Min+Period_Boundary_Zone_Width) .and. particle_position(i,2) >= Period_Boundary_Y_Min ) then

                Period_Boundary_Particle_Number=Period_Boundary_Particle_Number+1            ! Increase Period_Boundary_Particle_Number
                j=ture_total_particle_number+Period_Boundary_Particle_Number                 ! Particle index

                particle_position(j,1)=particle_position(i,1)
                particle_position(j,2)=particle_position(i,2)+Period_Boundary_Length_Y

                particle_velocity(j,:)=particle_velocity(i,:)
                
                particle_rho(j)=particle_rho(i)                                              ! Particle density
                particle_mass(j)=particle_mass(i)                                            ! Particle mass
                particle_type(j)=Period_Particle_Label                                       ! Particle type(Period_Particle_Label-Period boundary particles)
                particle_initial_type(j)=Period_Particle_Label                               ! Particle initial type
                particle_smooth_lengh(j)=particle_smooth_lengh(i)                            ! Particle smooth lengh
                particle_press(j)=particle_press(i)                                          ! Particle pressure
                particle_c(j)=particle_c(i)                                                  ! Particle c
                particle_energy(j)=particle_energy(i)                                        ! Particle energy
                particle_volume(j)=particle_volume(i)                                        ! Particle volume
                particle_smoke_value(j)=particle_smoke_value(i)                              ! Particle smoke line value

            end if
            !-------------------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------------------
            ! Move top inside to the bottom outside
            if( particle_position(i,2) >= ( Period_Boundary_Y_Max-Period_Boundary_Zone_Width) .and. particle_position(i,2) <= Period_Boundary_Y_Max ) then

                Period_Boundary_Particle_Number=Period_Boundary_Particle_Number+1            ! Increase Period_Boundary_Particle_Number
                j=ture_total_particle_number+Period_Boundary_Particle_Number                 ! Particle index

                particle_position(j,1)=particle_position(i,1)
                particle_position(j,2)=particle_position(i,2)-Period_Boundary_Length_Y

                particle_velocity(j,:)=particle_velocity(i,:)
                
                particle_rho(j)=particle_rho(i)                                              ! Particle density
                particle_mass(j)=particle_mass(i)                                            ! Particle mass
                particle_type(j)=Period_Particle_Label                                       ! Particle type(Period_Particle_Label-Period boundary particles)
                particle_initial_type(j)=Period_Particle_Label                               ! Particle initial type
                particle_smooth_lengh(j)=particle_smooth_lengh(i)                            ! Particle smooth lengh
                particle_press(j)=particle_press(i)                                          ! Particle pressure
                particle_c(j)=particle_c(i)                                                  ! Particle c
                particle_energy(j)=particle_energy(i)                                        ! Particle energy
                particle_volume(j)=particle_volume(i)                                        ! Particle volume
                particle_smoke_value(j)=particle_smoke_value(i)                              ! Particle smoke line value
                
            end if
            !-------------------------------------------------------------------------------------------------

        end do

    endif
    !===========================================================================================================
    



    !----------------------------------------------Very important-----------------------------------------------
    ! Renew the 'ture_total_particle_number' which contains : (1) Particle_ture_number;
    !                                                         (2) Wave_maker_particle_number;
    !                                                         (3) Fix_ghost_particle_number;
    !                                                         (4) InletOutlet_Boundary_Particle_Number;
    !                                                         (5) Period_Boundary_Particle_Number.
    ture_total_particle_number=Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number+Period_Boundary_Particle_Number

    !-----------------------------------------------------------------------------------------------------------





!     !***********************************************************************************************************
!     if ( i_time_step==Actual_start_time_step .and. Current_Processor_ID==Main_Processor ) then

!         ! Open the Output file
!         File_index = General_File_Port
!         File_name  = './Initial_Data/Check_Period_Particle_Distribution_Result.dat'

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

!             ! Output period boundary particle
!             if ( DistributePeriodBoundaryParticleOrNot==1 ) then
!                 write(File_index,*) "ZONE I=",Period_Boundary_Particle_Number," F=POINT"
!                 do k=1,Period_Boundary_Particle_Number
!                     i=k+Particle_ture_number+Wave_maker_particle_number+Fix_ghost_particle_number+InletOutlet_Boundary_Particle_Number
!                     write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),(particle_press(i)-Background_Pressure)/pressure_0
!                 enddo
!             endif


!         endif

!         close(File_index)

!     endif
!     !***********************************************************************************************************


end subroutine  Distribute_Period_Boundary_Particle
