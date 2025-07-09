!****************************************************************************
!
!  PROGRAM: wave_maker_damping_zone
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

subroutine wave_maker_damping_zone(i_time_step)              

    use information_module
    use Public_variable_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables form the mother subroutine
    integer,intent(in)::i_time_step                               !current time step

    !Variables in subroutine
    integer::i,j,k,L                                              !loop Variables
    real(kind=8)::damping_coefficient                             !damping coefficient
    real(kind=8)::current_time                                    !current time

    real(kind=8)::current_wave_angle                              !current wave angle
    real(kind=8)::time_soft_coefficient                           !soft coefficient for time
    real(kind=8)::current_wave_elevation                          !wave elevation at current particle z direction

    !Variables for file opeartion
    integer::ioerror=0                                            !open return value
    integer::stat                                                 !read return value
    integer::status                                               !allocation return value

    integer::file_index
    character(len=40)::file_name
    character(len=4)::char_Current_Processor_ID

    !==========================================================================================================

    ! Body of subroutine single_step_compute

    ! !*********************************************************************************************************
    ! !Get the MPI runing information for dugging, you can annotate this part when runnging
    ! call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    ! call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    ! !*********************************************************************************************************

    !==========================================================================================
    !For the wavemaker damping zone
    do i=1,actual_particle_number_in_subdomain

       !if(subdomain_particle_type(i)/=2 .or. subdomain_In_WaveMakerDampingZone_OrNot(i)/=1) cycle                !Only for the fluid particles

       ! !Define the particles in the wavemaker damping zone or not
       if(subdomain_In_WaveMakerDampingZone_OrNot(i)==1) then   

        ! if(subdomain_particle_type(i)/=2 ) cycle
        ! if(subdomain_particle_position(i,1)<=WavemakerDampingLengh) then

            ! write(*,*) subdomain_In_WaveMakerDampingZone_OrNot(i),Current_Processor_ID 

            !calculate the position damping coefficient
            
            ! if (Current_Processor_ID==1) then
            !    write(*,*) subdomain_particle_position(i,1),WavemakerDampingLengh,subdomain_position_soft_coefficient(i)
            ! endif
               
            ! !Renew the wavemaker damping zone 
            ! subdomain_particle_press(i)=position_soft_coefficient*subdomain_WaveMaker_Analytical_press(i)+(1-position_soft_coefficient)*subdomain_particle_press(i)   
            ! subdomain_particle_rho(i)=position_soft_coefficient*subdomain_WaveMaker_Analytical_rho(i)+(1-position_soft_coefficient)*subdomain_particle_rho(i)
            
            subdomain_particle_velocity(i,:)=subdomain_position_soft_coefficient(i)*subdomain_WaveMaker_Analytical_velocity(i,:)+(1-subdomain_position_soft_coefficient(i))*subdomain_particle_velocity(i,:)
            subdomain_particle_position(i,:)=subdomain_position_soft_coefficient(i)*subdomain_WaveMaker_Analytical_position(i,:)+(1-subdomain_position_soft_coefficient(i))*subdomain_particle_position(i,:)

            subdomain_particle_press(i)=subdomain_position_soft_coefficient(i)*subdomain_WaveMaker_Analytical_press(i)+(1-subdomain_position_soft_coefficient(i))*subdomain_particle_press(i)   
            subdomain_particle_rho(i)=water_rho_0*(1+7*subdomain_particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)

            !subdomain_particle_rho(i)=subdomain_position_soft_coefficient(i)*subdomain_WaveMaker_Analytical_rho(i)+(1-subdomain_position_soft_coefficient(i))*subdomain_particle_rho(i)
            
            ! subdomain_particle_velocity(i,:)=subdomain_WaveMaker_Analytical_velocity(i,:)
            ! subdomain_particle_position(i,:)=subdomain_WaveMaker_Analytical_position(i,:)

        end if 
    end do
        
    !==========================================================================================

    ! write(*,*) L,FluidParticleNumberInWaveMakerDamping


    



!     !------------------------------------------------------------------------------------------
!     !open the output tecplot data files
!     file_index=1+Current_Processor_ID

!     !transfer the Current_Processor_ID from integer to character
!     write(char_Current_Processor_ID,'(I4)') Current_Processor_ID
!     !write(*,*) char_Current_Processor_ID

!     file_name="Check_Iteration"//trim(adjustl(char_Current_Processor_ID))//".dat"
!     !write(*,*) file_name

!     open(unit=file_index,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    

!     !写入tecplot抬头
!     write(file_index,*) "TITLE='DISTRIBUTION'"
!     write(file_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "

!     !输出粒子
!     write(file_index,*) "ZONE I=",total_particle_number_in_subdomain," F=POINT"
!     do i=1,total_particle_number_in_subdomain

!         write(file_index,100) (subdomain_particle_position(i,j),j=1,dim),(particle_acceleration(i,j),j=1,dim),subdomain_particle_rho(i),subdomain_particle_press(i)/pressure_0
! 100      format(8F20.10) 
!     end do

!     close(file_index)
!     !------------------------------------------------------------------------------------------
    
    
!     !Note the user the current active Processor ID have finished the iterator information calculation
!     write(*,10) Current_Processor_ID
! 10   Format(" Processor:",I3," has finished iterator information calculation!")
    
    
    
end subroutine wave_maker_damping_zone  