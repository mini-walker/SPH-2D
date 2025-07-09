!*****************************************************************************************************************
!
!  SUBROUTINE : Exchange_buffer_domain_information
!
!  PURPOSE    : Exchange the pyhsical attributes of the particles in the buffer domain
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
!*****************************************************************************************************************

  subroutine Exchange_buffer_domain_information(i_time_step)

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !=============================================================================================================
    ! Variables from the mother subroutien
    integer,intent(in)::i_time_step                                     ! Current time step

    ! Variables in program
    integer::i,j,k,L,m                                                  ! Variables for looping

    !=============================================================================================================

    ! Body of Exchange_buffer_domain_information


    !*************************************************************************************************************


    call MPI_RealVector_MessageExchange(subdomain_particle_position)    ! Exchange the position
    call MPI_RealVector_MessageExchange(subdomain_particle_velocity)    ! Exchange the velocity
   
    call MPI_RealScalar_MessageExchange(subdomain_particle_press)       ! Exchange the pressure
    call MPI_RealScalar_MessageExchange(subdomain_particle_rho)         ! Exchange the density
    call MPI_RealScalar_MessageExchange(subdomain_particle_mass)        ! Exchange the mass
    call MPI_RealScalar_MessageExchange(subdomain_particle_volume)      ! Exchange the volume
    call MPI_RealScalar_MessageExchange(subdomain_particle_energy)      ! Exchange the energy
    call MPI_RealScalar_MessageExchange(subdomain_particle_c)           ! Exchange the sound velocity


    !*************************************************************************************************************

  end subroutine Exchange_buffer_domain_information

