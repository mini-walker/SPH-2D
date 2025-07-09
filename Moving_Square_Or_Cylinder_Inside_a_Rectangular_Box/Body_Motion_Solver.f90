!**************************************************************************************************************
!
!  SUBROUTINE : Body_Motion_Solver
!
!  PURPOSE    : Calculate the body force
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


subroutine Body_Motion_Solver(i_time_step)
 
    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI

    implicit none

    !==========================================================================================================
    ! Variables from superior subroutine
    integer::i_time_step                                              ! Present iteration time step

    ! Variables in local subroutine
    integer::i,j,k,l,v,o,m                                            ! Variables for looping

    real(kind=8)::present_time                                        ! Present time = simulation time - standing time

    integer::Custom_data_column_index                                 ! User-defined data location index

    integer::Custom_data_file_index                                   ! User-defined motion data file index
    character(len=100)::Custom_data_file_name                         ! User-defined motion data file name
    integer::Interpolation_result_file_index                          ! File index for outputting cubic spline interpolation result
    character(len=100)::Interpolation_result_file_name                ! File name for outputting cubic spline interpolation result

    character(len=10)::Char_temp                                      ! Character temp value for reading

    real(kind=8)::Interpolation_dt=0.0001d0                           ! Interpolation dt
    real(kind=8)::Interpolation_t                                     ! Interpolation t
    integer::Present_time_location                                    ! Present time location
    !==========================================================================================================


    ! Body of subroutine Body_Motion_Solver


    !==========================================================================================================
    ! Present time = simulation time - standing time
    Present_time=Real_time-Standing_time

    ! Free motion or custom motion: (1) 0 --- Free motion;
    !                                   1 --- User-defined motion.
    if ( Free_motion_or_Custom_motion==0 ) then

        !------------------------------------------------------------------------------------------------------
        ! Free motion

        !------------------------------------------------------------------------------------------------------

    else

        !------------------------------------------------------------------------------------------------------

        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! Calculate the cubic spline interpolation coefficient at the first iteration time step
        if ( i_time_step==Actual_start_time_step ) then
            
            Custom_data_file_index    = General_File_Port
            Custom_data_file_name     = './Input_Data/Motion_Body.dat'

            ! Open the custom motion data file
            call Initialziting_Reading_File( Custom_data_file_index,Custom_data_file_name,ioerror )
                
            Custom_data_column_number = 4             ! The first column is time and the left is the motion data

            ! Writing file open successful
            if ( ioerror==0 ) then
                
                !----------------------------------------------------------------------------------------------
                ! Get the motion data row number
                Custom_data_row_number=0
                do 
                    read(Custom_data_file_index,*,iostat=stat) Char_temp
                    if(stat/=0) exit

                    Custom_data_row_number=Custom_data_row_number+1

                enddo
                !----------------------------------------------------------------------------------------------


                !----------------------------------------------------------------------------------------------
                ! Allocate the memory for cubic spline interpolation
                Allocate ( User_defined_motion_data(Custom_data_row_number,Custom_data_column_number) )

                ! As the first column data is time, it doesn't has the cubic spline coefficient
                Allocate ( Cubic_Spline_Coeffcient (Custom_data_row_number,Custom_data_column_number,4) )
                Allocate ( Interpolation_Result (Custom_data_column_number) )

                ! Read the the user-defined motion data
                Rewind(Custom_data_file_index)               ! Return back to the file head

                do i=1,Custom_data_row_number

                    read(Custom_data_file_index,*) (User_defined_motion_data(i,j),j=1,Custom_data_column_number)
                        
                    ! if ( Current_Processor_ID==Main_Processor ) then
                    !     write(*,'(10F16.8)') (User_defined_motion_data(i,j),j=1,Custom_data_column_number)
                    ! endif

                enddo

                ! Synchronize all processors calculation
                ! MPI_BCAST function has the character "Synchronize", so we don't need Synchronize
                call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

                ! Calculate the cubic spline interpolation coefficient
                call Cubic_Spline_Interpolation( Custom_data_row_number,Custom_data_column_number,User_defined_motion_data,Cubic_Spline_Coeffcient )
                !----------------------------------------------------------------------------------------------


                !----------------------------------------------------------------------------------------------
                ! Open file for outputting
                Interpolation_result_file_name  = './Initial_Data/Cubic_Spline_Interpolation_Motion_Data.dat'
                Interpolation_result_file_index = Standby_General_File_Port
                call Initialziting_Writing_File( Interpolation_result_file_index,Interpolation_result_file_name,Standby_ioerror )    ! Input variables : File_index,File_name,IOERROR
                
                ! Synchronize all processors to make sure the file port is ready for all processors
                call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

                ! Output result file open successful
                if ( Standby_ioerror==0 ) then
                    
                    ! Calculate and output the cubic spline interpolation result
                    Interpolation_t=0.0d0

                    do i=1,1000000

                        Interpolation_t=Interpolation_t+Interpolation_dt

                        if ( Interpolation_t < User_defined_motion_data(1,1) .or. Interpolation_t > User_defined_motion_data(Custom_data_row_number,1) ) exit  ! Out the simulation time exit
                           
                        ! Check the location of Interpolation_t
                        Search_location: do k=1, Custom_data_row_number-1

                            if( Interpolation_t>User_defined_motion_data(k,1) .and. Interpolation_t<=User_defined_motion_data(k+1,1) ) exit Search_location
                            
                        enddo Search_location

                        ! Calculate and output the result 
                        do j=2,Custom_data_column_number
                            Interpolation_Result(j)=Cubic_Spline_Coeffcient(k,j,1)*( User_defined_motion_data(k+1,1)-Interpolation_t )**3 + Cubic_Spline_Coeffcient(k,j,2)*( Interpolation_t-User_defined_motion_data(k,1) )**3 + Cubic_Spline_Coeffcient(k,j,3)*( User_defined_motion_data(k+1,1)-Interpolation_t ) + Cubic_Spline_Coeffcient(k,j,4)*( Interpolation_t-User_defined_motion_data(k,1) )
                        enddo

                        if ( Current_Processor_ID==Main_Processor ) then
                            !write(*,'(10F16.8)') Interpolation_t,( Interpolation_Result(j),j=2,Custom_data_column_number )
                            write(Interpolation_result_file_index,'(8F16.8)') Interpolation_t,( Interpolation_Result(j),j=2,Custom_data_column_number )
                        
                        endif

                    enddo
                    
                endif


                ! Synchronize all processors calculation 
                ! MPI_BCAST function has the character "Synchronize", so we don't need Synchronize
                call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)

                ! You need to wait all the processors finish the file operation and then close the file prot
                close(Interpolation_result_file_index)
                !----------------------------------------------------------------------------------------------

            endif

            close(Custom_data_file_index)
            
        endif
        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! Initialize the custom value
        Custom_position=0.0d0                      ! User-defined position
        Custom_velocity=0.0d0                      ! User-defined velocity
        Custom_accleration=0.0d0                   ! User-defined accleration

        ! Have the user-defined data
        if ( Present_time<User_defined_motion_data(Custom_data_row_number,1) ) then
        
            !--------------------------------------------------------------------------------------------------
            ! Search the location of Present_time
            do k=1, Custom_data_row_number-1

                if( Present_time>User_defined_motion_data(k,1) .and. Present_time<=User_defined_motion_data(k+1,1) ) exit 
                
            enddo 

            Present_time_location=k

            ! Calculate and output the interpolation result 
            do j=2,Custom_data_column_number
                Interpolation_Result(j)=Cubic_Spline_Coeffcient(k,j,1)*( User_defined_motion_data(k+1,1)-Present_time )**3 + Cubic_Spline_Coeffcient(k,j,2)*( Present_time-User_defined_motion_data(k,1) )**3 + Cubic_Spline_Coeffcient(k,j,3)*( User_defined_motion_data(k+1,1)-Present_time ) + Cubic_Spline_Coeffcient(k,j,4)*( Present_time-User_defined_motion_data(k,1) )
            enddo
            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            ! User-defined accleration
            Custom_data_column_index=2

            Custom_accleration(1)=Interpolation_Result(Custom_data_column_index)
            Custom_accleration(2)=0.0d0

            ! User-defined velocity
            Custom_data_column_index=3

            Custom_velocity(1)=Interpolation_Result(Custom_data_column_index)
            Custom_velocity(2)=0.0d0

            ! User-defined position
            Custom_data_column_index=4

            Custom_position(1)=Interpolation_Result(Custom_data_column_index)
            Custom_position(2)=Block_Center_Y
            !--------------------------------------------------------------------------------------------------

        else

            !--------------------------------------------------------------------------------------------------
            Custom_position(1)=User_defined_motion_data(Custom_data_row_number,4)
            Custom_position(2)=Block_Center_Y

            Custom_velocity=0.0d0
            Custom_accleration=0.0d0
            !--------------------------------------------------------------------------------------------------

        endif


        Center_position_difference(1)=Custom_position(1)-Block_Center_X
        Center_position_difference(2)=Custom_position(2)-Block_Center_Y

        ! Renew the new center position
        Block_Center_X=Custom_position(1)
        Block_Center_Y=Custom_position(2)
        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ! Update the body particles
        do i=1,ture_total_particle_number

            ! Only the body particles are renewed
            if ( particle_initial_type(i) == Body_Particle_Label ) then

                particle_position(i,:)= particle_position(i,:)+Center_position_difference(:)

            endif

        enddo
        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


        !------------------------------------------------------------------------------------------------------

    endif
    !==========================================================================================================






 end subroutine Body_Motion_Solver
