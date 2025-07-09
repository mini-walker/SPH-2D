!**************************************************************************************************************
!
!  SUBROUTINE : Reinitialize_particle_density
!
!  PURPOSE    : Reinitialize particle density with MLS method
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

subroutine Reinitialize_particle_density(i_time_step)           

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from the superior subroutine
    integer,intent(in)::i_time_step                                      ! Current time step

    ! Variables in subroutine
    integer::i,j,k,L,v,o,m                                               ! Variables for loop

    ! Variables for file opeartion
    integer::File_index
    character(len=40)::File_name
    character(len=4)::Char_Current_Processor_ID

    ! Variables for kernel function
    integer::ture_effect_number                                          ! Number of effect particle
    integer,dimension(PredictionNumberInSupportDomain)::ture_effect_name ! Ture effect particle index
    real(kind=8),dimension(PredictionNumberInSupportDomain)::W_kernel    ! Kernel value
    real(kind=8)::sum_w_kernel                                           ! Sum kernel value

    !Variables for Refresh Density with MLS
    real(kind=8),dimension(PredictionNumberInSupportDomain)::Fai         ! The weight function value of MLS

    !==========================================================================================================

    ! Body of subroutine Reinitialize particle density
    
    !**********************************************************************************************************
    !Only calculate the actual particle in the subdomain, donot contain the buffer zone
    Reinitialize_subdomain_particle_rho=subdomain_particle_rho
    Reinitialize_subdomain_particle_velocity=subdomain_particle_velocity

    ! Here, all the particle density have been assgined as zero
    ! Please remember, I cost one day on this bug (Very important)
    subdomain_particle_rho=0.0d0                                    

    if (RefreshDensity_MLS_Order==1) then

      !--------------------------------------------------------------------------------------------------------
      ! Normalized MLS Interpolation
      do k=1,pair_number
          
          i=pair_i(k)                        ! particle i index in kth pair
          j=pair_j(k)                        ! particle j index in kth pair
          
          !Only for the same phase fluid 
          if (subdomain_particle_type(i)==subdomain_particle_type(j)) then

              if(subdomain_particle_type(i)==2) then
                subdomain_particle_rho(i)=subdomain_particle_rho(i)+subdomain_particle_mass(j)*w(k)/subdomain_sum_w(i)
              end if

              if(subdomain_particle_type(j)==2) then
                subdomain_particle_rho(j)=subdomain_particle_rho(j)+subdomain_particle_mass(i)*w(k)/subdomain_sum_w(j)
              end if
            
          endif

      enddo
      !--------------------------------------------------------------------------------------------------------

    else

      !--------------------------------------------------------------------------------------------------------
      ! Linear or quadric basis MLS Interpolation
      do i=1,actual_particle_number_in_subdomain

          if(subdomain_particle_type(i)==2) then

              !------------------------------------------------------------------------------------------------
              ! !For the particle near the free surface, we should not Reinitialize them
              ! if (subdomain_near_free_surface_or_not(i)==1) then

              !   subdomain_particle_rho(i)=Reinitialize_subdomain_particle_rho(i)  !Keep same

              !   cycle                                                             !Skip this loop
     
              ! endif
              !------------------------------------------------------------------------------------------------

              !------------------------------------------------------------------------------------------------
              ! !For the particle near the boundary, we should not Reinitialize them
              ! if (subdomain_near_boundary_or_not(i)==1) then

              !   subdomain_particle_rho(i)=Reinitialize_subdomain_particle_rho(i)  !Keep same

              !   cycle                                                             !Skip this loop
     
              ! endif
              !------------------------------------------------------------------------------------------------

              !Initialized current particle density
              subdomain_particle_rho(i)=0.0d0

              !Initialized the Variables for kernel function
              W_kernel=0.0d0                                           ! weight function value
              ture_effect_number=0                                     ! ture effect particle number
              ture_effect_name=0                                       ! ture effect particle index
            
              !Search all the vicinity particles
              do L=1,subdomain_effect_particle_number(i)

                  !effect particle index
                  j=subdomain_effect_particle(i,L)
                  
                  !position difference and distance between i and j particle
                  position_difference(:)=subdomain_particle_position(i,:)-subdomain_particle_position(j,:)
                  distance=sqrt(DOT_PRODUCT(position_difference,position_difference))
                  
                  !average smooth length (scaled by 1.3 to make sure the vicinity particles are enough)
                  !--------------------------------------------------------------------------------------------
                  !Important:
                  !If the Reinitialize support radius larger than 1.33dx, you should increase the buffer grid size 
                  !If you use just one buffer grid, it may cause uncertainty during the buffer zone.
                  !--------------------------------------------------------------------------------------------
                  average_smooth_length=(subdomain_particle_smooth_lengh(i)+subdomain_particle_smooth_lengh(j))/2.0
                  !write(*,*) particle_smooth_lengh(i)
                  
                  !particle in support domain or not for interpolation, we should contain the particle self
                  !The particle i is contained in support domain
                  if( distance<= kernel_scale*average_smooth_length ) then
                      
                      ture_effect_number=ture_effect_number+1              !ture effect particle number
                      ture_effect_name(ture_effect_number)=j               !ture effect particle index
                      
                      !call the compute kernel and get the weight function value
                      call compute_kernel(Kernel_function_model,&          ! Kernel function model(in)   
                                          dim,&                            !simulation dimension
                                          distance,&                       !distance(in)                
                                          position_difference,&            !position difference(in)
                                          average_smooth_length,&          !average smooth length(in)
                                          temp_w,&                         !weight function value(out)
                                          temp_dwdx&                       !derivates of weight function(out)
                                          )

                      !----------------------------------------------------------------------------------------
                      !Kernel value matrix, B martrix, P Basis matrix
                      W_kernel(ture_effect_number)=temp_w*subdomain_particle_volume(j)

                  endif
     
              enddo

              !================================================================================================
                              
              !Get the weight value of the effect particle
              call MLS_Interpolation_In_Subdomain_Domain(RefreshDensity_MLS_Order,subdomain_particle_position(i,:),ture_effect_number,ture_effect_name,W_kernel,Fai)
               
              !================================================================================================
              !Interpolate fixed particle quality
              do L=1,ture_effect_number
     
                  j=ture_effect_name(L)

                  subdomain_particle_rho(i)=subdomain_particle_rho(i)+Reinitialize_subdomain_particle_rho(j)*Fai(L)
          
              end do
              !------------------------------------------------------------------------------------------------
              !Check the density Interpolation results
              if ( isnan(subdomain_particle_rho(i)) ) then
                  subdomain_particle_rho(i)=Reinitialize_subdomain_particle_rho(i)
              endif
              !------------------------------------------------------------------------------------------------ 

              !================================================================================================
              !write(*,*) subdomain_particle_rho(i),Reinitialize_subdomain_particle_rho(i)

          else
            !Assgin the other particle density
            subdomain_particle_rho(i)=Reinitialize_subdomain_particle_rho(i)

          end if

      end do
      !-------------------------------------------------------------------------------------------------------

    endif

    ! !=======================================================================================================
    ! !Check the Reinitialize density
    ! do i=1,actual_particle_number_in_subdomain

    !     !if(subdomain_particle_type(i)/=2) cycle

    !     if (subdomain_particle_type(i)<=-201) then
    !       write(*,*) subdomain_particle_rho(i),subdomain_particle_volume(i),subdomain_particle_mass(i),subdomain_particle_type(i)
    !     endif
          
    ! enddo
    ! !=======================================================================================================

!     !*******************************************************************************************************

!     !------------------------------------------------------------------------------------------
!     !open the output tecplot data files
!     file_index=700+Current_Processor_ID

!     !transfer the Current_Processor_ID from integer to character
!     write(char_Current_Processor_ID,'(I4)') Current_Processor_ID
!     !write(*,*) char_Current_Processor_ID

!     file_name="./Subdomain/Density_in_subdomain_"//trim(adjustl(char_Current_Processor_ID))//".dat"
!     !write(*,*) file_name

!     open(unit=file_index,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    

!     !tecplot header
!     write(file_index,*) "TITLE='DISTRIBUTION'"
!     write(file_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
!     !------------------------------------------------------------------------------------------

!     !------------------------------------------------------------------------------------------
!     write(file_index,*) "ZONE I=",total_particle_number_in_subdomain," F=POINT"
!     do i=1,total_particle_number_in_subdomain
!        write(file_index,170) (subdomain_particle_position(i,j),j=1,dim),(subdomain_particle_velocity(i,j),j=1,dim),subdomain_particle_rho(i),subdomain_particle_press(i)/pressure_0
! 170          format(8F20.10) 
!     end do
!     !------------------------------------------------------------------------------------------

!     close(file_index)

!     !*******************************************************************************************************

end subroutine Reinitialize_particle_density
