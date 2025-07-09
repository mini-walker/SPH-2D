!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: exchange_buffer_domain_information
!
!  PURPOSE: Test wave maker with CSPH with ghost boundary method (Parallel -Version)
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

  subroutine exchange_buffer_domain_information(i_time_step)

    use information_module
    use Public_variable_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the mother subroutien
    integer,intent(in)::i_time_step                                     !current time step

    ! Variables in program
    integer::i,j,k,L,m                                                  !Variables for loop

    !==========================================================================================================

    ! Body of exchange_buffer_domain_information

    !*************************************************************************************************************
    ! Get the MPI runing information
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    !write(*,*) Current_Processor_ID
    !*************************************************************************************************************

    !-------------------------------------------------------------------------------------------------------------
    ! particle_position
    do j=1,dim

        left_Message_ID=j
        right_Message_ID=j

        !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
        !Data form left processor to right processor
        call MPI_SENDRECV(subdomain_particle_position(message_first_particle_index_to_right,j),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                          subdomain_particle_position(message_first_particle_index_from_left,j),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
        
        !Data form right processor to left processor
        call MPI_SENDRECV(subdomain_particle_position(message_first_particle_index_to_left,j),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                          subdomain_particle_position(message_first_particle_index_from_right,j),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    end do
    !-------------------------------------------------------------------------------------------------------------
    ! particle_velocity
    do j=1,dim

        left_Message_ID=(j+1)**2
        right_Message_ID=(j+1)**2

        !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
        !Data form left processor to right processor
        call MPI_SENDRECV(subdomain_particle_velocity(message_first_particle_index_to_right,j),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                          subdomain_particle_velocity(message_first_particle_index_from_left,j),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
        
        !Data form right processor to left processor
        call MPI_SENDRECV(subdomain_particle_velocity(message_first_particle_index_to_left,j),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                          subdomain_particle_velocity(message_first_particle_index_from_right,j),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    end do
    !-------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------
    !particle_press
    left_Message_ID=100
    right_Message_ID=101

    !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
    !Data form left processor to right processor
    call MPI_SENDRECV(subdomain_particle_press(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      subdomain_particle_press(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(subdomain_particle_press(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      subdomain_particle_press(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    !-------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------
    !particle_rho
    left_Message_ID=102
    right_Message_ID=103

    !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
    !Data form left processor to right processor
    call MPI_SENDRECV(subdomain_particle_rho(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      subdomain_particle_rho(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(subdomain_particle_rho(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      subdomain_particle_rho(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    !-------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------
    !particle_mass
    left_Message_ID=104
    right_Message_ID=105

    !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
    !Data form left processor to right processor
    call MPI_SENDRECV(subdomain_particle_mass(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      subdomain_particle_mass(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(subdomain_particle_mass(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      subdomain_particle_mass(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    !-------------------------------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------------------------------
    !particle_energy
    left_Message_ID=106
    right_Message_ID=107

    !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
    !Data form left processor to right processor
    call MPI_SENDRECV(subdomain_particle_energy(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      subdomain_particle_energy(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(subdomain_particle_energy(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      subdomain_particle_energy(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    !-------------------------------------------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------------------------------------------
    !particle_c
    left_Message_ID=108
    right_Message_ID=109

    !write(*,*) Current_Processor_ID,left_processor_ID,right_processor_ID
    !Data form left processor to right processor
    call MPI_SENDRECV(subdomain_particle_c(message_first_particle_index_to_right),message_length_send_to_right,MPI_DOUBLE_PRECISION,right_processor_ID,left_Message_ID,&
                      subdomain_particle_c(message_first_particle_index_from_left),message_length_from_left,MPI_DOUBLE_PRECISION,left_processor_ID,left_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)
    
    !Data form right processor to left processor
    call MPI_SENDRECV(subdomain_particle_c(message_first_particle_index_to_left),message_length_send_to_left,MPI_DOUBLE_PRECISION,left_processor_ID,right_Message_ID,&
                      subdomain_particle_c(message_first_particle_index_from_right),message_length_from_right,MPI_DOUBLE_PRECISION,right_processor_ID,right_Message_ID,MPI_COMM_WORLD,Status_MPI,ierror_MPI)

    !-------------------------------------------------------------------------------------------------------------

  end subroutine exchange_buffer_domain_information

