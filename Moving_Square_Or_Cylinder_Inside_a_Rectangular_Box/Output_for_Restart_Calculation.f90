!**************************************************************************************************************
!  SUBROUTINE:: Output_for_Restart_Calculation
!
!  PURPOSE: Save the data for restart calculation from breakpoint
!           Every time step will have an independ folder contain results
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


Subroutine Output_for_Restart_Calculation( i_time_step )
    
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
    
    integer::File_index                                                                 ! File index
    character(len=100)::File_name                                                       ! File name ( The file name character lenghth should be 100 )
    
    ! Variables for folder opeartion
    integer::Reserve_Folder_Number                                                      ! Reserve Folder Number
    real(kind=8)::Delete_time                                                           ! Delete time
    character(len=100)::Delete_folder_name                                              ! Delete folder name
    character(len=100)::Current_folder_name                                             ! Current folder name
    character(len=100)::RemoveCommand,BuildCommand                                      ! Remove and Build Command
    Logical::dirExists                                                                  ! Variables for folder opeartion

    !----------------------------------------------------------------------------------------------------------


    !==========================================================================================================


    ! Body of Subroutine Output_for_Restart_Calculation


    !----------------------------------------------------------------------------------------------------------
    ! Only save ten folders for restarting
    Reserve_Folder_Number=10

    Delete_time=Real_time-Reserve_Folder_Number*save_time_step_for_debug*dt

    write(char_time,'(F10.6)') Delete_time

    ! Delete the foler that not need to be reserved
    Delete_folder_name='./Restart_Data/'//trim(adjustl(char_time))

    ! Windows system
    if ( Linux_or_Windows==0 ) then                                        ! Windows command

        RemoveCommand='rd /s /q '//trim(adjustl(Delete_folder_name))

    ! Linux system
    elseif ( Linux_or_Windows==1 ) then                                    ! Linux command

        RemoveCommand='rm -r '//trim(adjustl(Delete_folder_name))

    else
        write(*,'(A)') " The value of 'Linux_or_Windows' in 'Information_Module' is not right! ( 0 or 1)" 
    endif

    ! Inquire the directory exist or not
    inquire(directory=trim(adjustl(Delete_folder_name)), exist=dirExists )
    
    ! Delete the folder
    if(dirExists) then

        call system (trim(adjustl(RemoveCommand)))                         ! Remove the Subdomain folder first

    endif   
    !----------------------------------------------------------------------------------------------------------



    !----------------------------------------------------------------------------------------------------------
    ! Current time
    write(char_time,'(F10.6)') Real_time

    ! Creat independ folder to contain the restart results of every time step
    current_folder_name='./Restart_Data/'//trim(adjustl(char_time))

    ! Windows system
    if ( Linux_or_Windows==0 ) then                                        ! Windows command

        RemoveCommand='rd /s /q '//trim(adjustl(current_folder_name))
        BuildCommand='md '//trim(adjustl(current_folder_name))

    ! Linux system
    elseif ( Linux_or_Windows==1 ) then                                    ! Linux command

        RemoveCommand='rm -r '//trim(adjustl(current_folder_name))
        BuildCommand ='mkdir -p '//trim(adjustl(current_folder_name))

    else
        write(*,'(A)') " The value of 'Linux_or_Windows' in 'Information_Module' is not right! ( 0 or 1)" 
    endif

    ! Inquire the directory exist or not
    inquire(directory=trim(adjustl(current_folder_name)), exist=dirExists )
    
    ! If the SPH calculation start from the first step, we need remove all the existed directories and the files in them,
    ! We just keep the old folders and results when start from the last time step.
    if(dirExists) then

        call system (trim(adjustl(RemoveCommand)))                         ! Remove the Subdomain folder first
        call system (trim(adjustl(BuildCommand)))                          ! Make the new folder

    else

        call system (trim(adjustl(BuildCommand)))                          ! Make the new folder

    endif                    
    !----------------------------------------------------------------------------------------------------------


    !----------------------------------------------------------------------------------------------------------
    ! Open the new file and save the results for restart from the desired breakpoint

    ! File name and file index
    File_name  = './Restart_Data/'//trim(adjustl(char_time))//'/'//trim(adjustl(char_time))//'.dat'
    File_index = Tecplot_General_File_Port
    call Initialziting_Writing_File( File_index,File_name,IOERROR )        ! Input variables : File_index,File_name,IOERROR

    if(IOERROR==0) then
        
        ! Tecplot header
        ! Initial particle type is used in body force calculation, so we should output it 
        ! The freesurface is re-calculated in the single step, we can ignore it 
        write(File_index,'(A)') 'TITLE="DISTRIBUTION"'
        write(File_index,'(A)') 'VARIABLES= "X" "Y" "VX" "VY" "Rho" "Mass" "E" "Smoke" "Type" "Ini_Type" "BoundaryType" "WaveLayer" "Level" '
 
        ! Output fluid particle
        write(File_index,'(A,I,A)') 'ZONE T=" '//trim(adjustl(char_time))//' ",'//' I=',particle_ture_number," F=POINT"
        do i=1,particle_ture_number
            write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_mass(i),particle_energy(i),particle_smoke_value(i),particle_type(i),particle_initial_type(i),Boundary_particle_type(i),wave_maker_particle_layer(i),particle_division_degree(i)          
100         format(8F16.8,5I) 
        enddo
        
        ! Output wave maker particle
        if (make_wave_or_not==1) then

            write(File_index,'(A,I,A)') 'ZONE T=" '//trim(adjustl(char_time))//' ",'//' I=',wave_maker_particle_number," F=POINT"
            do k=1,wave_maker_particle_number
                i=k+particle_ture_number
                write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_mass(i),particle_energy(i),particle_smoke_value(i),particle_type(i),particle_initial_type(i),Boundary_particle_type(i),wave_maker_particle_layer(i),particle_division_degree(i)          
            enddo

        endif

        ! Output fix ghost boundary particle
        if (DistributeFixedBoundaryParticleOrNot==1) then
            
            write(File_index,'(A,I,A)') 'ZONE T=" '//trim(adjustl(char_time))//' ",'//' I=',fix_ghost_particle_number," F=POINT"
            do k=1,fix_ghost_particle_number
                i=k+particle_ture_number+wave_maker_particle_number
                write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_mass(i),particle_energy(i),particle_smoke_value(i),particle_type(i),particle_initial_type(i),Boundary_particle_type(i),wave_maker_particle_layer(i),particle_division_degree(i)          
            enddo

        endif

        ! Inlet and outlet boundary particles
        if ( DistributeInletOutletBoundaryParticleOrNot==1 ) then
            
            write(File_index,'(A,I,A)') 'ZONE T=" '//trim(adjustl(char_time))//' ",'//' I=',InletOutlet_Boundary_Particle_Number," F=POINT"
            do k=1,InletOutlet_Boundary_Particle_Number
                i=k+particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
                write(File_index,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_mass(i),particle_energy(i),particle_smoke_value(i),particle_type(i),particle_initial_type(i),Boundary_particle_type(i),wave_maker_particle_layer(i),particle_division_degree(i)          
            enddo
            
        endif

        ! Close the file index
        close(File_index)

    endif
    !----------------------------------------------------------------------------------------------------------

end subroutine Output_for_Restart_Calculation