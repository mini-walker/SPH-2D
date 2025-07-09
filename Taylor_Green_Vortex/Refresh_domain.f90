!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: Refresh_domain
!
!  PURPOSE: Refresh domain
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

  subroutine Refresh_domain(i_time_step)

    use information_module
    use Public_variable_module
    use function_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the mother subroutine
    integer,intent(in)::i_time_step                                     !current time step 

    ! Variables in program
    integer::i,j,k,L,m                                                  !Variables for loop
    
    !Variables for file opeartion
    integer::ioerror=0                                                  !open return value
    integer::stat                                                       !read return value
    integer::status                                                     !allocation return value

    real(kind=8),dimension(dim)::kth_pair_dwdx,velocity_difference

    integer::file_index
    character(len=40)::file_name
    character(len=4)::char_Current_Processor_ID  

    !==========================================================================================================

    ! Body of CSPH_MPI_Parallel_wave_maker_2D

    !*************************************************************************************************************
    ! Get the MPI runing information
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    !*************************************************************************************************************


    !write(*,*) Current_Processor_ID,average_particle_number_per_process

    !*************************************************************************************************************
    !Some Variables donot need to be Refreshed when they are not updated during the calculation
    !-------------------------------------------------------------------------------------------------------------
    !Refresh domain

    !As we have devide all particles into all subdomain, then we can define all domain information to be zero

    particle_position=0.0d0
    particle_velocity=0.0d0
    particle_mass=0.0d0
    particle_rho=0.0d0
    particle_smooth_lengh=0.0d0
    particle_c=0.0d0
    particle_press=0.0d0
    particle_energy=0.0d0

    particle_type=0
    free_surface_type=0


    ! !***********************************************************************************************************
    ! !The main Processor Output the initial particle information
    ! if (Current_Processor_ID==Main_Processor .and. i_time_step==1) then

    !     !Open the Output file
    !     open(unit=5,file="Before_Refresh.dat",status="replace",position="rewind",action="write",iostat=ioerror)    
    !     if(ioerror==0) then
            
    !         !input tecplot header
    !         write(5,*) "TITLE='DISTRIBUTION'"
    !         write(5,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
    !         !output fluid particle
    !         write(5,*) "ZONE I=",particle_ture_number," F=POINT"
    !         do i=1,particle_ture_number
    !              write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    ! 100          format(8F20.10) 
    !         end do
            
    !         !output wave maker particle
    !         write(5,*) "ZONE I=",wave_maker_particle_number," F=POINT"
    !         do k=1,wave_maker_particle_number
    !             i=k+particle_ture_number
    !             write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    !         end do

    !         !output fix ghost boundary particle
    !         write(5,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
    !         do k=1,fix_ghost_particle_number
    !             i=k+particle_ture_number+wave_maker_particle_number
    !             write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    !         end do

    !         !output ghost of fix ghost boundary particle
    !         write(5,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
    !         do k=1,fix_ghost_particle_number
    !             i=k+particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
    !             write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    !         end do
            
    !     else
    !         write(*,*) "The open action of the Initial_particle_information.dat is fail!"
    !     end if
    !     close(5)

    ! endif

    ! ! ! Synchronize all processors calculation
    ! ! call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
    ! ! !***********************************************************************************************************
    
    !--------------------------------------------------------------------------------------------------------------
    !particle position
    Temp_Reduce_Real_Array=0.0d0
    
    do j=1,dim

      do i=1,actual_particle_number_in_subdomain                      
        
         Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_position(i,j)                            

      end do


      call MPI_REDUCE( Temp_Reduce_Real_Array,particle_position(:,j),n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
      call MPI_BCAST( particle_position(:,j),n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    end do 
    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !particle velocity
    Temp_Reduce_Real_Array=0.0d0
    
    do j=1,dim

      do i=1,actual_particle_number_in_subdomain                       
        
         Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_velocity(i,j)                            

      end do


      call MPI_REDUCE( Temp_Reduce_Real_Array,particle_velocity(:,j),n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
      call MPI_BCAST( particle_velocity(:,j),n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    end do 
    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !particle press
    Temp_Reduce_Real_Array=0.0d0

    do i=1,actual_particle_number_in_subdomain                      
      
       Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_press(i)                            

    end do


    call MPI_REDUCE( Temp_Reduce_Real_Array,particle_press,n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_press,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !particle mass
    Temp_Reduce_Real_Array=0.0d0

    do i=1,actual_particle_number_in_subdomain                          
      
       Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_mass(i)                            

    end do


    call MPI_REDUCE( Temp_Reduce_Real_Array,particle_mass,n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_mass,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !particle density
    Temp_Reduce_Real_Array=0.0d0

    do i=1,actual_particle_number_in_subdomain             

       Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_rho(i)                            

    end do


    call MPI_REDUCE( Temp_Reduce_Real_Array,particle_rho,n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_rho,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !particle smooth lengh
    Temp_Reduce_Real_Array=0.0d0

    do i=1,actual_particle_number_in_subdomain                        
      
       Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_smooth_lengh(i)                            

    end do


    call MPI_REDUCE( Temp_Reduce_Real_Array,particle_smooth_lengh,n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_smooth_lengh,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !particle c
    Temp_Reduce_Real_Array=0.0d0

    do i=1,actual_particle_number_in_subdomain                         
      
       Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_c(i)                            

    end do


    call MPI_REDUCE( Temp_Reduce_Real_Array,particle_c,n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_c,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !particle_energy
    Temp_Reduce_Real_Array=0.0d0

    do i=1,actual_particle_number_in_subdomain                         
      
       Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_energy(i)                            

    end do


    call MPI_REDUCE( Temp_Reduce_Real_Array,particle_energy,n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_energy,n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !particle_type
    Temp_Reduce_Int_Array=0

    do i=1,actual_particle_number_in_subdomain             

       Temp_Reduce_Int_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_type(i)                            

    end do


    call MPI_REDUCE( Temp_Reduce_Int_Array,particle_type,n_total,MPI_INTEGER,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( particle_type,n_total,MPI_INTEGER,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !free_surface_type
    Temp_Reduce_Int_Array=0

    do i=1,actual_particle_number_in_subdomain             

       Temp_Reduce_Int_Array(actual_particle_name_in_subdomain(i))=subdomain_free_surface_type(i)                            

    end do


    call MPI_REDUCE( Temp_Reduce_Int_Array,free_surface_type,n_total,MPI_INTEGER,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
    call MPI_BCAST( free_surface_type,n_total,MPI_INTEGER,Main_Processor,MPI_COMM_WORLD,ierror_MPI )

    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !particle position
    Normal_vector=0.0d0
    Temp_Reduce_Real_Array=0.0d0
    
    do j=1,dim

      do i=1,actual_particle_number_in_subdomain                      
        
         Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_Normal_vector(i,j)                            

      end do

      call MPI_REDUCE( Temp_Reduce_Real_Array,Normal_vector(:,j),n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
      call MPI_BCAST( Normal_vector(:,j),n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
 
    end do 
    !---------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------
    !Calculate the vorticity

    subdomain_particle_vorticity=0.0d0

    do k=1,pair_number
        
       i=pair_i(k)                        !particle i index in kth pair
       j=pair_j(k)                        !particle j index in kth pair

       velocity_difference(:)=velocity_difference_ij(k,:)

       subdomain_particle_vorticity(i,:)=subdomain_particle_vorticity(i,:)+subdomain_particle_volume(j)*Cross_Product(velocity_difference,Normalized_dwdx_i(k,:),dim)
       
       subdomain_particle_vorticity(j,:)=subdomain_particle_vorticity(j,:)+subdomain_particle_volume(i)*Cross_Product(velocity_difference,Normalized_dwdx_j(k,:),dim)

    end do


    !particle vorticity
    particle_vorticity=0.0d0
    Temp_Reduce_Real_Array=0.0d0
    
    do j=1,3

      do i=1,actual_particle_number_in_subdomain                      
        
         Temp_Reduce_Real_Array(actual_particle_name_in_subdomain(i))=subdomain_particle_vorticity(i,j)                            

      end do

      call MPI_REDUCE( Temp_Reduce_Real_Array,particle_vorticity(:,j),n_total,MPI_DOUBLE_PRECISION,MPI_SUM,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
      call MPI_BCAST( particle_vorticity(:,j),n_total,MPI_DOUBLE_PRECISION,Main_Processor,MPI_COMM_WORLD,ierror_MPI )
 
    end do 

    !---------------------------------------------------------------------------------------------------------------
    ! !***********************************************************************************************************


    ! ! !***********************************************************************************************************
    ! ! when only one Processor
    ! if (Total_Processors_Number==1) then

    !     actual_particle_number_in_subdomain=particle_ture_number

    !     k=0

    !     do j=1,ture_total_particle_number

    !        k=k+1

    !         i=Full_subdomain_particle_name(k)

    !         particle_position(k,:)=subdomain_particle_position(i,:)
    !         particle_velocity(k,:)=subdomain_particle_velocity(i,:)
    !         particle_mass(k)=subdomain_particle_mass(i)                          
    !         particle_type(k)=subdomain_particle_type(i)                                 
    !         particle_initial_type(k)=subdomain_particle_initial_type(i)                         
    !         particle_rho(k)=subdomain_particle_rho(i)                                                 
    !         particle_smooth_lengh(k)=subdomain_particle_smooth_lengh(i)                    
    !         particle_c(k)=subdomain_particle_c(i)                               
    !         particle_press(k)=subdomain_particle_press(i)                           
    !         particle_energy(k)=subdomain_particle_energy(i)                          
    !         particle_eta(k)=subdomain_particle_eta(i)                             
    !         free_surface_type(k)=subdomain_free_surface_type(i)                             
    !         wave_maker_particle_layer(k)=subdomain_wave_maker_particle_layer(i)

    !     end do

    ! endif

    ! ! !***********************************************************************************************************


    ! !***********************************************************************************************************
    ! !The main Processor Output the initial particle information
    ! if (Current_Processor_ID==Main_Processor .and. i_time_step==1) then

    !     !Open the Output file
    !     open(unit=5,file="After_Refresh.dat",status="replace",position="rewind",action="write",iostat=ioerror)    
    !     if(ioerror==0) then
            
    !         !input tecplot header
    !         write(5,*) "TITLE='DISTRIBUTION'"
    !         write(5,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
    !         !output fluid particle
    !         write(5,*) "ZONE I=",particle_ture_number," F=POINT"
    !         do i=1,particle_ture_number
    !              write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    !         end do
            
    !         !output wave maker particle
    !         write(5,*) "ZONE I=",wave_maker_particle_number," F=POINT"
    !         do k=1,wave_maker_particle_number
    !             i=k+particle_ture_number
    !             write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    !         end do

    !         !output fix ghost boundary particle
    !         write(5,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
    !         do k=1,fix_ghost_particle_number
    !             i=k+particle_ture_number+wave_maker_particle_number
    !             write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    !         end do

    !         !output ghost of fix ghost boundary particle
    !         write(5,*) "ZONE I=",fix_ghost_particle_number," F=POINT"
    !         do k=1,fix_ghost_particle_number
    !             i=k+particle_ture_number+wave_maker_particle_number+fix_ghost_particle_number
    !             write(5,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    !         end do
            
    !     else
    !         write(*,*) "The open action of the Initial_particle_information.dat is fail!"
    !     end if
    !     close(5)

    ! endif
    ! !***********************************************************************************************************



  end subroutine Refresh_domain

