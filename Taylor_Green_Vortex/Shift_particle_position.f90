!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE:: Shift_particle_position
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

subroutine Shift_particle_position(i_time_step)         

    use information_module
    use Public_variable_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables from mother subroutine

    integer,intent(in)::i_time_step                               !current time step

    !Variables in subroutine
    integer::i,j,k,L                                              !loop Variables
    real(kind=8),dimension(dim)::temp_dwdx                        !derivates of kernel weight
    real(kind=8)::w_dx,temp_w                                     !kernel function at dx
    real(kind=8),dimension(dim)::position_difference              !position difference
    real(kind=8),dimension(dim)::velocity_difference              !velocity difference
    real(kind=8)::distance                                        !distance between particle i and j
    real(kind=8)::average_smooth_length                           !average smooth length
    real(kind=8)::i_velocity,j_velocity,max_velocity_ij           !max velocity of i and j
    real(kind=8),dimension(dim)::iTH_shift_distance               !iTH shift distance
    real(kind=8)::shift_distance_amplitude                        !shift distance amplitude
    real(kind=8)::Temp_Rho                                        !Temp density
    real(kind=8),dimension(dim)::kth_pair_dwdx_i,kth_pair_dwdx_j  !Normalized derivates of kernel function

    !==========================================================================================================

    ! Body of subroutine Shift_particle_position

    !**********************************************************************************************************
    !Method 1: P.K. Stanby ----This is better for no free surface simulation
    !----------------------------------------------------------------------------------------------------------
    ! calculate w(dx) 
    position_difference=0.0d0
    distance=interior_dx
    average_smooth_length=smooth_length
    
    call compute_kernel(dim,&                            !simulation dimension
                        distance,&                       !distance(in)                
                        position_difference,&            !position difference(in)
                        average_smooth_length,&          !smooth length(in)
                        temp_w,&                         !kernel function(out)
                        temp_dwdx&                       !derivates of kernel function(out)
                        )
                     
    W_dx=temp_w

    !write(*,*) W_dx,w(2),smooth_length,distance,sqrt((subdomain_shift_particle_position(pair_i(2),1)-subdomain_shift_particle_position(pair_j(2),1))**2+(subdomain_shift_particle_position(pair_i(2),2)-subdomain_shift_particle_position(pair_j(2),2))**2)
    !----------------------------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------------------------
    !Calculate the shift position
    
    subdomain_shift_particle_position=0.0d0
    
    do k=1,pair_number
       
       i=pair_i(k)                        !particle i index in kth pair
       j=pair_j(k)                        !particle j index in kth pair
    
       ! !-------------------------------------------------------------------------------------------------------
       ! !Standard derivates
       ! do L=1,dim
       !    subdomain_shift_particle_position(i,L)=subdomain_shift_particle_position(i,L)-shift_belta*(smooth_length**2)*(1+shift_R*(w(k)/w_dx)**shift_n)*dwdx(k,L)*subdomain_particle_volume(j)
       !    subdomain_shift_particle_position(j,L)=subdomain_shift_particle_position(j,L)+shift_belta*(smooth_length**2)*(1+shift_R*(w(k)/w_dx)**shift_n)*dwdx(k,L)*subdomain_particle_volume(i)
       ! end do
       ! !-------------------------------------------------------------------------------------------------------

       !-------------------------------------------------------------------------------------------------------
       !Normalized derivates
       !The normalized derivates of weight function
       kth_pair_dwdx_i(:)=Normalized_dwdx_i(k,:)
       kth_pair_dwdx_j(:)=Normalized_dwdx_j(k,:)

       do L=1,dim
          subdomain_shift_particle_position(i,L)=subdomain_shift_particle_position(i,L)-shift_belta*(smooth_length**2)*(1+shift_R*(w(k)/w_dx)**shift_n)*kth_pair_dwdx_i(L)*subdomain_particle_volume(j)
          subdomain_shift_particle_position(j,L)=subdomain_shift_particle_position(j,L)+shift_belta*(smooth_length**2)*(1+shift_R*(w(k)/w_dx)**shift_n)*kth_pair_dwdx_j(L)*subdomain_particle_volume(i)
       end do
       !-------------------------------------------------------------------------------------------------------

    end do
    !----------------------------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------------------------
    !Shift particle position
    do i=1,actual_particle_number_in_subdomain
        
       if(subdomain_particle_type(i)/=2) cycle

       if(subdomain_near_free_surface_or_not(i)==1) cycle           !Only the interior particles not the free surface

       iTH_shift_distance(:)=subdomain_shift_particle_position(i,:)

       shift_distance_amplitude=sqrt(DOT_PRODUCT(iTH_shift_distance,iTH_shift_distance))
       
       if (shift_distance_amplitude>max_shift) then             !maximum is 0.15

          subdomain_shift_particle_position(i,:)=subdomain_shift_particle_position(i,:)/shift_distance_amplitude*max_shift
         
       endif

       subdomain_particle_position(i,:)=subdomain_particle_position(i,:)+subdomain_shift_particle_position(i,:)

       !Renew the particle density and pressure
       do L=1,dim
          subdomain_particle_rho(i)=subdomain_particle_rho(i)+grad_rho(i,L)*subdomain_shift_particle_position(i,L)
       end do

       subdomain_particle_press(i)=(square_c_0*water_rho_0/7.0)*((subdomain_particle_rho(i)/water_rho_0)**7-1.0)

       ! subdomain_particle_press(i)=subdomain_particle_press(i)-water_rho_0*g*subdomain_shift_particle_position(i,dim)
       ! subdomain_particle_rho(i)=water_rho_0*(1+7*subdomain_particle_press(i)/(square_c_0*water_rho_0))**(1.0/7.0)
    
    end do
    !---------------------------------------------------------------------------------------------------------

    !**********************************************************************************************************


    ! !**********************************************************************************************************
    ! !Method 2: R. Vacondio, B. D. Rogers
    ! !Variable resolution for SPH: A dynamic particle coalesing and splitting scheme
    ! !----------------------------------------------------------------------------------------------------------
    ! !Calculate the shift position
    
    ! subdomain_shift_particle_position=0.0d0
    ! sum_mass=0.0d0
    ! sum_distance=0.0d0
    ! near_particle_max_velocity=0.0d0
    ! shift_effect_particle_number=0

    ! do k=1,pair_number
       
    !    i=pair_i(k)                        !particle i index in kth pair
    !    j=pair_j(k)                        !particle j index in kth pair
    
    !    !--------------------------------------------------------------------
    !    !Sum mass, distance and particle number
    !    sum_mass(i)=sum_mass(i)+subdomain_particle_mass(j)                     ! Sum mass of near particles
    !    sum_mass(j)=sum_mass(j)+subdomain_particle_mass(i)

    !    !position difference
    !    position_difference=0.0d0
    !    distance=0.0d0
    !    do L=1,dim
    !       position_difference(L)=subdomain_particle_position(i,L)-subdomain_particle_position(j,L)
    !       distance=distance+position_difference(L)**2
    !    end do
    !    distance=sqrt(distance)

    !    sum_distance(i)=sum_distance(i)+distance                                ! Sum distance of near particles
    !    sum_distance(j)=sum_distance(j)+distance

    !    shift_effect_particle_number(i)=shift_effect_particle_number(i)+1       ! effect particle number
    !    shift_effect_particle_number(j)=shift_effect_particle_number(j)+1
    !    !--------------------------------------------------------------------

    !    !--------------------------------------------------------------------
    !    !Max velocity
    !    i_velocity=0.0d0
    !    j_velocity=0.0d0
    !    do L=1,dim
    !         i_velocity=i_velocity+(subdomain_particle_velocity(i,L)**2)
    !         j_velocity=j_velocity+(subdomain_particle_velocity(j,L)**2)
    !    end do
    !    i_velocity=sqrt(i_velocity)
    !    j_velocity=sqrt(j_velocity)

    !    if (i_velocity>=j_velocity) then
    !       max_velocity_ij=i_velocity
    !    else
    !       max_velocity_ij=j_velocity
    !    endif

    !    if (near_particle_max_velocity(i)<=max_velocity_ij) then

    !       near_particle_max_velocity(i)=max_velocity_ij
         
    !    endif

    !    if (near_particle_max_velocity(j)<=max_velocity_ij) then

    !       near_particle_max_velocity(j)=max_velocity_ij

    !    endif
    !    !--------------------------------------------------------------------

    !    do L=1,dim
    !       subdomain_shift_particle_position(i,L)=subdomain_shift_particle_position(i,L)+position_difference(L)*subdomain_particle_mass(j)/((distance)**3)
    !       subdomain_shift_particle_position(j,L)=subdomain_shift_particle_position(j,L)-position_difference(L)*subdomain_particle_mass(i)/((distance)**3)
    !    end do

    ! end do
    ! !----------------------------------------------------------------------------------------------------------

    ! !----------------------------------------------------------------------------------------------------------
    ! !Shift particle position
    ! do i=1,actual_particle_number_in_subdomain
        
    !    if(subdomain_particle_type(i)/=2) cycle

    !    if(subdomain_near_free_surface_or_not(i)==1) cycle        !For the particle near the free surface, we should not Reinitialize them

    !    subdomain_shift_particle_position(i,:)=subdomain_shift_particle_position(i,:)*shift_belta*(sum_distance(i)/shift_effect_particle_number(i))**2*near_particle_max_velocity(i)*dt/sum_mass(i)

    !    iTH_shift_distance(:)=subdomain_shift_particle_position(i,:)

    !    shift_distance_amplitude=sqrt(DOT_PRODUCT(iTH_shift_distance,iTH_shift_distance))

    !    if (shift_distance_amplitude>=max_shift) then             !maximum is 0.1dx

    !      do L=1,dim
    !         subdomain_shift_particle_position(i,L)=subdomain_shift_particle_position(i,L)/shift_distance_amplitude*max_shift
    !      end do

    !    endif

    !    !Renew the density and pressure after shifting
    !    do L=1,dim
    !       subdomain_particle_rho(i)=subdomain_particle_rho(i)+grad_rho(i,L)*subdomain_shift_particle_position(i,L)
    !    end do

    !    subdomain_particle_press(i)=(square_c_0*water_rho_0/7.0)*((subdomain_particle_rho(i)/water_rho_0)**7-1.0)

    !    ! subdomain_particle_position(i,:)=subdomain_particle_position(i,:)+subdomain_shift_particle_position(i,:)

    !    ! !对粒子压力和密度进行重新赋值
    !    ! subdomain_particle_press(i)=subdomain_particle_press(i)-water_rho_0*g*subdomain_shift_particle_position(i,dim)
    !    ! !subdomain_particle_rho(i)=water_rho_0*(1+7*subdomain_particle_press(i)/(square_c_0*water_rho_0))**(1.0/7.0)
    
    ! end do
    ! !---------------------------------------------------------------------------------------------------------

    !************************************************************************************************************




end subroutine Shift_particle_position 