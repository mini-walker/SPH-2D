!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: find_particle_pair
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

 subroutine find_particle_pair()
 
    use information_module
    use Public_variable_module
    use MPI

    implicit none

    !==========================================================================================================
    
    ! Variables in subroutine
    integer::i,j,k,L,v                                                  !Loop Variables
    real(kind=8),dimension(dim)::temp_dwdx                              !Temp derivates of weight function
    real(kind=8)::distance                                              !Distance between i and j particle
    real(kind=8)::average_smooth_length                                 !Average smooth length
    real(kind=8),dimension(dim)::position_difference                    !Position difference

    integer::file_index
    character(len=40)::file_name
    character(len=4)::char_Current_Processor_ID  
    
    !==========================================================================================================

    ! Body of find_particle_pair

    !----------------------------------------------------------------------------------------------------------
    ! Initialization  
    pair_i=0
    pair_j=0
    w=0.0
    dwdx=0.0
    
    rho_difference_ij=0.0d0                                  !第k个粒子对中密度差(rho_i-rho_j)
    position_difference_ij=0.0d0                             !第k个粒子对中位置差(position_difference_i-position_difference_j)
    velocity_difference_ij=0.0d0                             !第k个粒子对中速度差(velocity_difference_i-velocity_difference_j)
    average_c_ij=0.0d0                                       !第k个粒子对声速平均值(particle_c(i)+particle_c(j))/2.0
    average_rho_ij=0.0d0                                     !第k个粒子对密度平均值(particle_rho(i)+particle_rho(j))/2.0
    average_smooth_length_ij=0.0d0                           !第k个粒子对光滑长度平均值(particle_smooth_lengh(i)+particle_smooth_lengh(j))/2.0


    !Initialized pair_number
    pair_number=0                                                       !**********这个对数一定要记得初始化************
    
    do i=1,total_particle_number_in_subdomain-1                         !从第i个点一直到n-1个点
    
        !Search all the vicinity particles        
        do L=1,subdomain_effect_particle_number(i)
    
            !vicinity particle index
            j=subdomain_effect_particle(i,L)

            if(j<=i .or. j>total_particle_number_in_subdomain) cycle    !Only for the particle index larger than current particle


            !Delete the pair between boundary particle
            if(subdomain_particle_type(i)/=2 .and. subdomain_particle_type(j)/=2) cycle     !没有流体粒子参与的粒子对
            
            !position difference and distance between i and j particle
            position_difference(:)=subdomain_particle_position(i,:)-subdomain_particle_position(j,:)

            distance=sqrt(DOT_PRODUCT(position_difference,position_difference))
         
            !average smooth length
            average_smooth_length=(subdomain_particle_smooth_lengh(i)+subdomain_particle_smooth_lengh(j))/2.0
            
            !particle in support domain or not and Delete the particle in same location
            if(distance<=(kernel_scale*average_smooth_length) .and. distance>=1.0E-7) then
                
                if(pair_number<=pair_n) then
                    
                    pair_number=pair_number+1                           !pair_number increase 1
                    pair_i(pair_number)=i                               !particile index in kth pair is i
                    pair_j(pair_number)=j                               !particile index in kth pair is j
                
                    !call the weight function
                    call compute_kernel(dim,&                            !dimension
                                        distance,&                       !distance(in)                
                                        position_difference,&            !position difference(in)
                                        average_smooth_length,&          !average smooth length(in)
                                        w(pair_number),&                 !kernel function value for kth particle pair(out)
                                        temp_dwdx&                       !derivates for kth particle pair(out)
                                        )
                    
                    !Save the results
                    do k=1,dim
                        dwdx(pair_number,k)=temp_dwdx(k)
                        position_difference_ij(pair_number,k)=position_difference(k)                                                     !第k个粒子对中位置差(position_difference_i-position_difference_j)
                        velocity_difference_ij(pair_number,k)=subdomain_particle_velocity(i,k)-subdomain_particle_velocity(j,k)          !第k个粒子对中速度差(velocity_difference_i-velocity_difference_j)
                    end do

                    rho_difference_ij(pair_number)=subdomain_particle_rho(i)-subdomain_particle_rho(j)                                    !第k个粒子对中密度差(rho_i-rho_j)
                    average_c_ij(pair_number)=(subdomain_particle_c(i)+subdomain_particle_c(j))/2.0                                       !第k个粒子对声速平均值(particle_c(i)+particle_c(j))/2.0
                    average_rho_ij(pair_number)=(subdomain_particle_rho(i)+subdomain_particle_rho(j))/2.0                                 !第k个粒子对密度平均值(particle_rho(i)+particle_rho(j))/2.0
                    average_smooth_length_ij(pair_number)=average_smooth_length                                                           !第k个粒子对光滑长度平均值(particle_smooth_lengh(i)+particle_smooth_lengh(j))/2.0
 
                else
                    write(*,*) "!**ERROR**! There are too many particle pair,Please add pair_n!"
                end if
            end if
            
        end do
    
    end do
    
    !----------------------------------------------------------------------------------------------------------


!     !----------------------------------------------------------------------------------------------------------
!     ! Get the MPI runing information
!     call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
!     call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
!     !----------------------------------------------------------------------------------------------------------
  

!     if (Current_Processor_ID==Main_Processor) then
        
!         !Open the Output file
!         open(unit=5,file="Check_particle_pair.dat")    

!         do k=1,pair_number

!             i=pair_i(k)                        !particle i index in kth pair
!             j=pair_j(k)                        !particle j index in kth pair
            
!             write(5,*) k,i,j 
!         end do

!         close(5)
!         !Open the Output file
!         open(unit=5,file="particle_pair.dat")    

!         do i=1,all_particle
            
!             write(5,100) i,subdomain_particle_position(i,1),subdomain_particle_position(i,2)
! 100         format(I5,8F20.10) 
!         end do

!         close(5)
!     endif

!     !----------------------------------------------------------------------------------------------------------
 
 end subroutine find_particle_pair