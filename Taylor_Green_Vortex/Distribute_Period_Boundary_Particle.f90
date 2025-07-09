!****************************************************************************
!
!  PROGRAM: Parallel_Taylor_Green_Vortices_2D_Release
!
!  SUBROUTINE: Distribute_Period_Boundary_Particle
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

subroutine  Distribute_Period_Boundary_Particle(i_time_step)

    use Public_variable_module
    use information_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the mother subroutien
    integer,intent(in)::i_time_step                                     !current time step 

    ! Variables in program
    integer::i,j,k,L,m                                                  !Variables for loop
    
    !Variables for file opeartion
    integer::ioerror=0                                                  !open return value
    integer::stat                                                       !read return value
    integer::status                                                     !allocation return value
    
    integer::file_index                                                 !file index
    character(len=40)::file_name                                        !file name
    character(len=4)::char_Current_Processor_ID                         !Current Processor ID in character

    integer::Particle_Index                                             !Particle Index
    integer::Temp_Period_Boundary_Particle_Number                       !Temp Period Boundary Particle Number

    !==========================================================================================================

    ! Body of subroutine Distribute_Period_Boundary_Particle
   
    !=========================================================================================================
    !Move the particle out the boundary

    !Inside particle period motion
    do i=1,particle_ture_number

        if(particle_type(i)==120) cycle        !No error particles

        !-----------------------------------------------------------------------------------------------------
        !Move right outside to the left inside and left outside to the right inside
        if( particle_position(i,1) > Period_Boundary_X_Max ) then

            particle_position(i,1)=particle_position(i,1)-Period_Boundary_Length_X
            particle_position(i,2)=particle_position(i,2)

        elseif( particle_position(i,1) < Period_Boundary_X_Min ) then

            particle_position(i,1)=particle_position(i,1)+Period_Boundary_Length_X
            particle_position(i,2)=particle_position(i,2)

        endif
        !-----------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------
        !Move top outside to the bottom inside and bottom outside to the top inside
        if( particle_position(i,2) > Period_Boundary_Y_Max ) then

            particle_position(i,1)=particle_position(i,1)
            particle_position(i,2)=particle_position(i,2)-Period_Boundary_Length_Y

        elseif(particle_position(i,2) < Period_Boundary_Y_Min) then

            particle_position(i,1)=particle_position(i,1)
            particle_position(i,2)=particle_position(i,2)+Period_Boundary_Length_Y

        endif
        !-----------------------------------------------------------------------------------------------------

    end do

    !=========================================================================================================

    !=========================================================================================================
    Period_Boundary_Particle_Number=0

    !Left and right Period
    do i=1,particle_ture_number

        if(particle_type(i)==120) cycle        !No error particles
        
        !-----------------------------------------------------------------------------------------------------
        !Move left inside to the right outside
        if( particle_position(i,1) <= ( Period_Boundary_X_Min+kernel_scale*smooth_length) .and. particle_position(i,1) >= Period_Boundary_X_Min ) then

            Period_Boundary_Particle_Number=Period_Boundary_Particle_Number+1           !Increase Period_Boundary_Particle_Number
            j=particle_ture_number+Period_Boundary_Particle_Number                      !Particle_Index

            particle_position(j,1)=particle_position(i,1)+Period_Boundary_Length_X
            particle_position(j,2)=particle_position(i,2)

            particle_velocity(j,:)=particle_velocity(i,:)

            particle_rho(j)=particle_rho(i)                                              !particle density
            particle_mass(j)=particle_mass(i)                                            !particle mass
            particle_type(j)=-101                                                        !particle type(-101-Period boundary particles)
            particle_smooth_lengh(j)=particle_smooth_lengh(i)                            !particle smooth lengh
            particle_press(j)=particle_press(i)                                          !particle pressure
            particle_c(j)=particle_c(i)                                                  !particle c
            particle_energy(j)=particle_energy(i)                                        !particle energy

        end if
        !-----------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------
        !Move right inside to the left outside
        if( particle_position(i,1) >= ( Period_Boundary_X_Max-kernel_scale*smooth_length) .and. particle_position(i,1) <= Period_Boundary_X_Max ) then

            Period_Boundary_Particle_Number=Period_Boundary_Particle_Number+1           !Increase Period_Boundary_Particle_Number
            j=particle_ture_number+Period_Boundary_Particle_Number                      !Particle_Index

            particle_position(j,1)=particle_position(i,1)-Period_Boundary_Length_X
            particle_position(j,2)=particle_position(i,2)

            particle_velocity(j,:)=particle_velocity(i,:)
            
            particle_rho(j)=particle_rho(i)                                              !particle density
            particle_mass(j)=particle_mass(i)                                            !particle mass
            particle_type(j)=-101                                                        !particle type(-101-Period boundary particles)
            particle_smooth_lengh(j)=particle_smooth_lengh(i)                            !particle smooth lengh
            particle_press(j)=particle_press(i)                                          !particle pressure
            particle_c(j)=particle_c(i)                                                  !particle c
            particle_energy(j)=particle_energy(i)                                        !particle energy
            
        end if
        !-----------------------------------------------------------------------------------------------------

    end do
    !=========================================================================================================


    !=========================================================================================================
    Temp_Period_Boundary_Particle_Number=Period_Boundary_Particle_Number

    !Top and Bottom Period
    do i=1,particle_ture_number+Temp_Period_Boundary_Particle_Number

        if(particle_type(i)==120) cycle        !No error particles
        
        !-----------------------------------------------------------------------------------------------------
        !Move bottom inside to the top outside
        if( particle_position(i,2) <= ( Period_Boundary_Y_Min+kernel_scale*smooth_length) .and. particle_position(i,2) >= Period_Boundary_Y_Min ) then

            Period_Boundary_Particle_Number=Period_Boundary_Particle_Number+1           !Increase Period_Boundary_Particle_Number
            j=particle_ture_number+Period_Boundary_Particle_Number                      !Particle_Index

            particle_position(j,1)=particle_position(i,1)
            particle_position(j,2)=particle_position(i,2)+Period_Boundary_Length_Y

            particle_velocity(j,:)=particle_velocity(i,:)
            
            particle_rho(j)=particle_rho(i)                                              !particle density
            particle_mass(j)=particle_mass(i)                                            !particle mass
            particle_type(j)=-101                                                        !particle type(-101-Period boundary particles)
            particle_smooth_lengh(j)=particle_smooth_lengh(i)                            !particle smooth lengh
            particle_press(j)=particle_press(i)                                          !particle pressure
            particle_c(j)=particle_c(i)                                                  !particle c
            particle_energy(j)=particle_energy(i)                                        !particle energy
            
        end if
        !-----------------------------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------------------------
        !Move top inside to the bottom outside
        if( particle_position(i,2) >= ( Period_Boundary_Y_Max-kernel_scale*smooth_length) .and. particle_position(i,2) <= Period_Boundary_Y_Max ) then

            Period_Boundary_Particle_Number=Period_Boundary_Particle_Number+1           !Increase Period_Boundary_Particle_Number
            j=particle_ture_number+Period_Boundary_Particle_Number                      !Particle_Index

            particle_position(j,1)=particle_position(i,1)
            particle_position(j,2)=particle_position(i,2)-Period_Boundary_Length_Y

            particle_velocity(j,:)=particle_velocity(i,:)
            
            particle_rho(j)=particle_rho(i)                                              !particle density
            particle_mass(j)=particle_mass(i)                                            !particle mass
            particle_type(j)=-101                                                        !particle type(-101-Period boundary particles)
            particle_smooth_lengh(j)=particle_smooth_lengh(i)                            !particle smooth lengh
            particle_press(j)=particle_press(i)                                          !particle pressure
            particle_c(j)=particle_c(i)                                                  !particle c
            particle_energy(j)=particle_energy(i)                                        !particle energy
            
        end if
        !-----------------------------------------------------------------------------------------------------

    end do
    !===========================================================================================================
    
!     !=========================================================================================================

!     !Call MPI functions
!     call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
!     call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )

!     if (Current_Processor_ID==Main_Processor .and. i_time_step==1) then

!         !Open the Output file
!         open(unit=5,file="Period_Particle_information.dat",status="replace",position="rewind",action="write",iostat=ioerror)    
!         if(ioerror==0) then
!             write(*,*) "The open action of the Period_Particle_information.dat is Fine!"
            
!             !input tecplot header
!             write(5,*) "TITLE='DISTRIBUTION'"
!             write(5,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
!             !output fluid particle
!             write(5,*) "ZONE I=",particle_ture_number," F=POINT"
!             do i=1,particle_ture_number
!                  write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
! 100              format(8F20.10) 
!             end do
            
!             !output wave maker particle
!             if (make_wave_or_not==1) then
!                 write(5,*) "ZONE I=",wave_maker_particle_number," F=POINT"
!                 do k=1,wave_maker_particle_number
!                     i=k+particle_ture_number
!                     write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
!                 end do
!             endif

!             !output fix ghost boundary particle
!             if (DistributeFixedBoundaryParticleOrNot==1) then
!                 write(5,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
!                 do k=1,fix_ghost_particle_number
!                     i=k+particle_ture_number+wave_maker_particle_number
!                     write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
!                 end do
!             endif

!             !output period boundary particle
!             if (DistributePeriodBoundaryParticleOrNot==1) then
!                 write(5,*) "ZONE I=",Period_Boundary_Particle_Number," F=POINT"
!                 do k=1,Period_Boundary_Particle_Number
!                     i=k+particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
!                     write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
!                 end do
!             endif 
            
!         else
!             write(*,*) "The open action of the Period_Particle_information.dat is fail!"
!         end if
!         close(5)

!     endif
!     !***********************************************************************************************************


end subroutine  Distribute_Period_Boundary_Particle