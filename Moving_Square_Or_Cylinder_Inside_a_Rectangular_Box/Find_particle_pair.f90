!**************************************************************************************************************
!  SUBROUTINE : Find_particle_pair
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

 subroutine Find_particle_pair()
 
    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI

    implicit none

    !==========================================================================================================
    
    ! Variables in subroutine
    integer::i,j,k,L,v                                                  ! Loop Variables

    ! Variables for file opeartion
    integer::File_index                                                 ! File index
    character(len=100)::File_name                                       ! File name
    character(len=4)::Char_Current_Processor_ID                         ! Current Processor ID in character
    
    !==========================================================================================================

    ! Body of find_particle_pair

    !----------------------------------------------------------------------------------------------------------
    ! Initialization  
    pair_i=0
    pair_j=0
    
    w=0.0
    dwdx=0.0
    
    rho_difference_ij=0.0d0                                  ! Density difference     = (rho_i-rho_j)
    position_difference_ij=0.0d0                             ! Position difference    = (position_difference_i-position_difference_j)
    velocity_difference_ij=0.0d0                             ! Velocity difference    = (velocity_difference_i-velocity_difference_j)
    average_c_ij=0.0d0                                       ! Average Sound velocity = (particle_c(i)+particle_c(j))/2.0
    average_rho_ij=0.0d0                                     ! Average density        = (particle_rho(i)+particle_rho(j))/2.0
    average_smooth_length_ij=0.0d0                           ! Average smooth length  = (particle_smooth_lengh(i)+particle_smooth_lengh(j))/2.0


    ! Initialized pair_number
    pair_number=0                                                       
     
    do i=1,Total_particle_number_in_subdomain-1                          ! Total_particle_number_in_subdomain-1
    
        ! Search all the vicinity particles        
        do L=1,subdomain_effect_particle_number(i)
    
            ! Vicinity particle index
            j=subdomain_effect_particle(i,L)

            if( j<=i .or. j>Total_particle_number_in_subdomain ) cycle   ! Only for the particle index larger than current particle

            ! Delete the pair between boundary particle and inlet boundary particles
            if( subdomain_particle_type(i)==Ghost_Water_Particle_Label .and. subdomain_particle_type(j)==Ghost_Water_Particle_Label ) cycle    
            if( subdomain_particle_type(i)==Fix_Ghost_Particle_Label   .and. subdomain_particle_type(j)==Fix_Ghost_Particle_Label   ) cycle    
            if( subdomain_particle_type(i)==Inlet_Particle_Label       .and. subdomain_particle_type(j)==Inlet_Particle_Label       ) cycle 

            ! Position difference and distance between i and j particle
            position_difference(:)=subdomain_particle_position(i,:)-subdomain_particle_position(j,:)

            distance=sqrt(DOT_PRODUCT(position_difference,position_difference))
         
            ! Average smooth length
            average_smooth_length=(subdomain_particle_smooth_lengh(i)+subdomain_particle_smooth_lengh(j))/2.0
            
            ! Particle in support domain or not and Delete the particle in same location
            if( distance<=(kernel_scale*average_smooth_length) .and. distance>=1.0E-10 ) then
                
                if( pair_number < pair_n ) then
                    
                    pair_number=pair_number+1                            ! pair_number increase 1
                    pair_i(pair_number)=i                                ! particile index in kth pair is i
                    pair_j(pair_number)=j                                ! particile index in kth pair is j
                
                    ! Call the weight function
                    call compute_kernel(Kernel_function_model,&          ! Kernel function model(in)  
                                        dim,&                            ! Dimension
                                        distance,&                       ! Distance(in)                
                                        position_difference,&            ! Position difference(in)
                                        average_smooth_length,&          ! Average smooth length(in)
                                        w(pair_number),&                 ! Kernel function value for kth particle pair(out)
                                        temp_dwdx&                       ! Derivates for kth particle pair(out)
                                        )
                    
                    ! Save the results
                    dwdx(pair_number,:)=temp_dwdx(:)
                    position_difference_ij(pair_number,:)=position_difference(:)                                                     
                    velocity_difference_ij(pair_number,:)=subdomain_particle_velocity(i,:)-subdomain_particle_velocity(j,:)          


                    rho_difference_ij(pair_number)=subdomain_particle_rho(i)-subdomain_particle_rho(j)                                    
                    average_c_ij(pair_number)=0.5d0*(subdomain_particle_c(i)+subdomain_particle_c(j))                                       
                    average_rho_ij(pair_number)=0.5d0*(subdomain_particle_rho(i)+subdomain_particle_rho(j))                               
                    average_smooth_length_ij(pair_number)=average_smooth_length                                                           
 
                else
                    write(*,'(A)') " !**ERROR**! There are too many particle pair,Please add pair_n!"
                endif

            endif
            
        enddo
    
    enddo
    !----------------------------------------------------------------------------------------------------------








!     !----------------------------------------------------------------------------------------------------------
!     ! Output results for debugging
!     if ( Current_Processor_ID==Main_Processor .and. i_time_step==Actual_start_time_step ) then
        
!         ! Open the Output file
!         File_index = General_File_Port
!         File_name  = './Initial_Data/Check_Particle_Pair.dat'

!         call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
    
!         if( IOERROR==0 ) then
 
!             do k=1,pair_number

!                 i=pair_i(k)                        ! Particle i index in kth pair
!                 j=pair_j(k)                        ! Particle j index in kth pair
                
!                 write(File_index,*) k,i,j 
!             end do

!             close(File_index)

!         endif


!         ! Open the Output file
!         File_index = General_File_Port
!         File_name  = './Initial_Data/Check_Particle_Pair_Position.dat'

!         call Initialziting_Writing_File( File_index,File_name,IOERROR )    ! Input variables : File_index,File_name,IOERROR
    
!         if( IOERROR==0 ) then   

!             do i=1,Total_particle_number_in_subdomain
                
!                 write(File_index,100) i,subdomain_particle_position(i,1),subdomain_particle_position(i,2)
! 100             format(I8,8F20.10) 

!             enddo

!             close(File_index)

!         endif

!     endif
!     !----------------------------------------------------------------------------------------------------------
 
 end subroutine Find_particle_pair