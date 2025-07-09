!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: Solver_Normalized_MUSCL_Upwind_Linear_Riemann_Delta
!
!  PURPOSE: Compute iterator information for time integration with Normalized kernel function Delta SPH and Modified MUSCL Riemann solver
!
!  Programer: Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location: MUN
!
!  Time: 2017.3.18
!
!  Copyright: Memorial University
!
!  Version: 1.0
!
!  Note: MPI version: mpich-3.2
!        Fortran version: Inter Fortran (ifort) 18.0.0
!****************************************************************************

 subroutine Solver_Normalized_MUSCL_Upwind_Linear_Riemann_Delta()
 
    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI

    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,l,v,o,m                                          !Variables for loop


    !SPS Variables
    real(kind=8)::sps_blin,sps_k,sps_niu_t,SPS_Delta_L
    real(kind=8)::master_element_square,last_element_square
    real(kind=8)::ABS_SPS_Strain_Tensor_Square


    real(kind=8)::hxx,hyy,hzz,hxy,hxz,hyz,hvcc,rho_ij               !
    real(kind=8)::pressure_difference
    
    !artifical viscosity term
    real(kind=8)::average_c                                         !average c(4.68)
    real(kind=8)::average_rho                                       !average rho(4.69)
    real(kind=8),dimension(dim)::x_ij,v_ij                          !position and velocity difference
    real(kind=8)::v_ij_multiply_x_ij                                !x_ij times v_ij
    real(kind=8)::capital_pai_ij                                    !capital_pai
    real(kind=8)::capital_fai_ij_1                                  !4.67
    real(kind=8),dimension(dim)::capital_fai_ij                     !4.67
    real(kind=8)::temp,temp_i,temp_j                                !temp value
    real(kind=8)::square_distance                                   !square distance
    real(kind=8)::temp_e,temp_e_i,temp_e_j
    real(kind=8)::temp_p,temp_p_i,temp_p_j

    !Laminar laplacian term
    real(kind=8)::RijDotProductDerivates                            !Rij*D Wij
    
    !Delta SPH
    real(kind=8)::rho_difference                                                   !rho difference
    real(kind=8),dimension(dim)::psai                                              !psai in Delta-SPH
    real(kind=8)::psai_ij                                                          !psai in modified Delta-SPH
    real(kind=8),dimension(dim)::sum_grad_rho_ij                                   !grad_rho_i+grad_rho_j


    real(kind=8),dimension(dim)::dwdx_i,dwdx_j                                     !Normalized derivates of kernel function
    real(kind=8),dimension(dim)::kth_pair_dwdx,kth_pair_dwdx_i,kth_pair_dwdx_j     !Normalized derivates of kernel function

    !External Boundary force
    real(kind=8)::boundary_force_chi
    real(kind=8)::boundary_force_alpha
    real(kind=8)::boundary_force_F_alpha

    !Riemann slover 
    real(kind=8)::P_Left,P_right,rho_Left,rho_right,U_Left,U_right
    real(kind=8)::P_Bar,U_Bar                                                      !Average pressure and velocity
    real(kind=8),dimension(dim)::I_Velocity,J_Velocity,IJ_Velocity_Bar
    real(kind=8)::I_Velocity_Module,J_Velocity_Module
    real(kind=8),dimension(dim)::IJ_Direction_Vector
    real(kind=8)::U_asterisk,P_asterisk,I_P_asterisk,J_P_asterisk
    real(kind=8),dimension(dim)::velocity_asterisk,I_velocity_asterisk,J_velocity_asterisk
    real(kind=8)::Belta_Riemman
    real(kind=8)::DeltaScalar_I_Left,DeltaScalar_I_Right,DeltaScalar_J_Left,DeltaScalar_J_Right
    real(kind=8),dimension(dim,dim)::Transfer_Matrix                               !Transfer Matrix from global to local
    real(kind=8),dimension(dim)::Velocity_Left,Velocity_Right
    real(kind=8),dimension(dim)::Velocity_Bar_0

    !MUSCL
    real(kind=8)::DelPhi_i,DelPhi_j,DelPhi_ij,DelPhiLim_i,DelPhiLim_j
    real(kind=8)::Scalar_I_Left,Scalar_I_Right,Scalar_J_Left,Scalar_J_Right         !Scalar value of I or J left and right
    real(kind=8)::DeltaScalar_I_Forward,DeltaScalar_I_Backward,DeltaScalar_J_Forward,DeltaScalar_J_Backward
    
    real(kind=8),dimension(dim)::Vector_I_Left,Vector_I_Right,Vector_J_Left,Vector_J_Right
    real(kind=8),dimension(dim)::DeltaVector_I_Forward,DeltaVector_I_Backward,DeltaVector_J_Forward,DeltaVector_J_Backward
    
    integer::check_error=0

    !==========================================================================================================
    !DOT_PRODUCT(vector_a,vector_b) two vector DOT_PRODUCT
    !MATMUL(matrix_a,matrix_b)two matrix DOT_PRODUCT

    ! Body of subroutine compute_iterator_information

    ! !**********************************************************************************************************
    ! ! Get the MPI runing information
    ! call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
    ! call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )
    ! !**********************************************************************************************************
    
    !***************************************XSPH correction***********************************

    !-----------------------------------------------------------------------------------------
    ! Calculate XSPH term (To make the simulation stable)
    do k=1,pair_number
        
        i=pair_i(k)                                         !particle i index in pair k
        j=pair_j(k)                                         !particle j index in pair k

        velocity_difference(:)=velocity_difference_ij(k,:)  ! particle i-j
        
        average_velocity(i,:)=average_velocity(i,:)-epsilon*subdomain_particle_mass(j)*velocity_difference(:)/(subdomain_particle_rho(i)+subdomain_particle_rho(j))*w(k)
        average_velocity(j,:)=average_velocity(j,:)+epsilon*subdomain_particle_mass(i)*velocity_difference(:)/(subdomain_particle_rho(i)+subdomain_particle_rho(j))*w(k)                   

    end do
    !-----------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------
    !Exchange the XSPH term
    call MPI_RealVector_MessageExchange(average_velocity)  ! Exchange the gradience of the XSPH term
    
    !-----------------------------------------------------------------------------------------

    !*****************************************************************************************


    !***************************************************************************************** 
    
    !-----------------------------------------------------------------------------------------
    !Calculate the gradience of density, pressure and velocity derivates

    grad_rho=0.0d0
    grad_pressure=0.0d0
    grad_velocity=0.0d0

    forward_grad_pressure=0.0d0
    backward_grad_pressure=0.0d0
    forward_grad_velocity=0.0d0
    backward_grad_velocity=0.0d0   

    do k=1,pair_number
       
       i=pair_i(k)                                             !i index in pair k
       j=pair_j(k)                                             !j index in pair k

       !--------------------------------------------------------------------------------------
       !Density difference
       rho_difference=-rho_difference_ij(k)                    !Attention: j-i

       !Density gradience
       grad_rho(i,:)=grad_rho(i,:)+rho_difference*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)
       grad_rho(j,:)=grad_rho(j,:)+rho_difference*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)
       !--------------------------------------------------------------------------------------


       !--------------------------------------------------------------------------------------
       !backward and forward gradience

       !velocity difference
       velocity_difference(:)=-velocity_difference_ij(k,:)                                 !Attention: j-i

       !Pressure difference
       pressure_difference=subdomain_particle_press(j)-subdomain_particle_press(i)         !Attention: j-i


       !======================================================================================
       !Upwind direction is the position direction
       if (Position_Or_Velocity_Direction==0) then

            !---------------------------------------------------------------------------------
           if ( subdomain_particle_position(j,1)>=subdomain_particle_position(i,1) ) then  ! i forward and j backward

                do L=1,dim
                    forward_grad_velocity(i,L,:) =forward_grad_velocity(i,L,:) +2*velocity_difference(L)*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)
                    backward_grad_velocity(j,L,:)=backward_grad_velocity(j,L,:)+2*velocity_difference(L)*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)
                enddo
                
                forward_grad_pressure(i,:) =forward_grad_pressure(i,:) +2*pressure_difference*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)
                backward_grad_pressure(j,:)=backward_grad_pressure(j,:)+2*pressure_difference*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)
           
           else                                                                            !i backward and j forward

                do L=1,dim
                    backward_grad_velocity(i,L,:)=backward_grad_velocity(i,L,:)+2*velocity_difference(L)*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)
                    forward_grad_velocity(j,L,:) =forward_grad_velocity(j,L,:) +2*velocity_difference(L)*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)
                enddo
                
                backward_grad_pressure(i,:)=backward_grad_pressure(i,:)+2*pressure_difference*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)
                forward_grad_pressure(j,:) =forward_grad_pressure(j,:) +2*pressure_difference*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)
           
           endif
            !---------------------------------------------------------------------------------

        !Upwind direction is the velocity direction
       else

            !---------------------------------------------------------------------------------
            !Position difference
            position_difference=-position_difference_ij(k,:)                                !j-i 
            
            !Unit vector From i to j
            IJ_Direction_Vector=position_difference/sqrt(DOT_PRODUCT(position_difference,position_difference)) 

            I_Velocity=subdomain_particle_velocity(i,:)!+average_velocity(i,:)
            J_Velocity=subdomain_particle_velocity(j,:)!+average_velocity(j,:)
            
            I_Velocity_Module=sqrt(DOT_PRODUCT(I_Velocity,I_Velocity))
            J_Velocity_Module=sqrt(DOT_PRODUCT(J_Velocity,J_Velocity))

            I_Upwind_Direction=0.0d0
            if ( I_Velocity_Module>=1.0E-8) then
                I_Upwind_Direction=I_Velocity/I_Velocity_Module
            else
                I_Upwind_Direction(1)=1.0d0
            endif

            J_Upwind_Direction=0.0d0
            if ( J_Velocity_Module>=1.0E-8 ) then
                J_Upwind_Direction=J_Velocity/J_Velocity_Module
            else
                J_Upwind_Direction(1)=1.0d0
            endif


            !Upwind of particle i
            if ( DOT_PRODUCT(IJ_Direction_Vector,I_Upwind_Direction)>=0.0d0 ) then

                do L=1,dim
                    forward_grad_velocity(i,L,:) =forward_grad_velocity(i,L,:) +2*velocity_difference(L)*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)
                enddo
                
                forward_grad_pressure(i,:) =forward_grad_pressure(i,:) +2*pressure_difference*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)

            !Backwind of particle i
            else

                do L=1,dim
                    backward_grad_velocity(i,L,:)=backward_grad_velocity(i,L,:)+2*velocity_difference(L)*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)
                enddo
                
                backward_grad_pressure(i,:)=backward_grad_pressure(i,:)+2*pressure_difference*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)

            endif

            !Upwind of particle j
            if (DOT_PRODUCT(IJ_Direction_Vector,J_Upwind_Direction)<=0.0d0) then

                do L=1,dim
                    forward_grad_velocity(j,L,:) =forward_grad_velocity(j,L,:) +2*velocity_difference(L)*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)
                enddo
                
                forward_grad_pressure(j,:) =forward_grad_pressure(j,:) +2*pressure_difference*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)
           
            else

                do L=1,dim
                    backward_grad_velocity(j,L,:)=backward_grad_velocity(j,L,:)+2*velocity_difference(L)*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)
                enddo
                
                backward_grad_pressure(j,:)=backward_grad_pressure(j,:)+2*pressure_difference*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)

                
            endif

            !---------------------------------------------------------------------------------

       endif
       !======================================================================================


    end do 

    grad_pressure=0.50*(backward_grad_pressure+forward_grad_pressure) 
    grad_velocity=0.50*(backward_grad_velocity+forward_grad_velocity)

    !-----------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------
    !Exchange the grad_rho form the neighbour processor

    call MPI_RealVector_MessageExchange(grad_rho)                        ! Exchange the gradience of the density
    call MPI_RealVector_MessageExchange(grad_pressure)                   ! Exchange the gradience of the pressure
    call MPI_RealVector_MessageExchange(forward_grad_pressure)           ! Exchange the gradience of the density
    call MPI_RealVector_MessageExchange(backward_grad_pressure)          ! Exchange the gradience of the pressure


    call MPI_RealTensor_MessageExchange(grad_velocity)                   ! Exchange the gradience of the horizontal velocity
    call MPI_RealTensor_MessageExchange(forward_grad_velocity)           ! Exchange the gradience of the horizontal velocity
    call MPI_RealTensor_MessageExchange(backward_grad_velocity)          ! Exchange the gradience of the horizontal velocity

    !-----------------------------------------------------------------------------------------


    !*****************************************************************************************




    !***********************************Stand SPH calculation*********************************
    do k=1,pair_number
        
        i=pair_i(k)                                     !particle i index in kth pair
        j=pair_j(k)                                     !particle j index in kth pair

        !Get the direction vector for from i to j
        x_ij(:)=-position_difference_ij(k,:)            !j-i 
        distance=sqrt(DOT_PRODUCT(x_ij,x_ij))
        IJ_Direction_Vector=x_ij/distance               !R_ji

        !=====================================================================================
        !MUSCL Scheme
        !Reconstrute the left and right value in the Riemann solver

        !----------Pressure---------
        DeltaScalar_I_Forward =DOT_PRODUCT(forward_grad_pressure(i,:), x_ij(:))
        DeltaScalar_I_Backward=DOT_PRODUCT(backward_grad_pressure(i,:),x_ij(:))

        ! Input variables order: Delta_i_backward,Delta_i_forward,DelPhi_i
        call Upwind_MUSCL_Limiter(DeltaScalar_I_Backward,DeltaScalar_I_Forward,DelPhiLim_i)

        !The pressure in Riemann problem
        Scalar_I_Left =subdomain_particle_press(i)-0.5*DelPhiLim_i
        Scalar_I_Right=subdomain_particle_press(i)+0.5*DelPhiLim_i



        DeltaScalar_J_Forward =-DOT_PRODUCT(forward_grad_pressure(j,:), x_ij(:))
        DeltaScalar_J_Backward=-DOT_PRODUCT(backward_grad_pressure(j,:),x_ij(:))

        ! Input variables order: Delta_i_backward,Delta_i_forward,DelPhi_i
        call Upwind_MUSCL_Limiter(DeltaScalar_J_Backward,DeltaScalar_J_Forward,DelPhiLim_j)

        !The pressure in Riemann problem
        Scalar_J_Left =subdomain_particle_press(j)-0.5*DelPhiLim_j
        Scalar_J_Right=subdomain_particle_press(j)+0.5*DelPhiLim_j

        !-------------------------------------------------------------------------------------
        !Pressure Status for Riemann Solver
        ! I Left and J left
        P_Left =subdomain_particle_press(i)+0.5*(Scalar_I_Left-Scalar_I_Right)
        P_Right=subdomain_particle_press(j)+0.5*(Scalar_J_Left-Scalar_J_Right)

        ! ! I Right and J Right
        ! P_Left =subdomain_particle_press(i)-0.5*(Scalar_I_Left-Scalar_I_Right)
        ! P_Right=subdomain_particle_press(j)-0.5*(Scalar_J_Left-Scalar_J_Right)

        ! ! I Right and J Left
        ! P_Left =subdomain_particle_press(i)-0.5*(Scalar_I_Left-Scalar_I_Right)
        ! P_Right=subdomain_particle_press(j)+0.5*(Scalar_J_Left-Scalar_J_Right)
        !-------------------------------------------------------------------------------------
        
        P_Bar=(P_Left+P_Right)/2.0 

        !----------Velocity---------
        velocity_difference(:)=-velocity_difference_ij(k,:)               !j-1 
        
        !Loop for 3 direction
        do L=1,dim

            DeltaScalar_I_Forward =DOT_PRODUCT(forward_grad_velocity(i,L,:), x_ij(:))
            DeltaScalar_I_Backward=DOT_PRODUCT(backward_grad_velocity(i,L,:),x_ij(:))

            ! Input variables order: Delta_i_backward,Delta_i_forward,DelPhi_i
            call Upwind_MUSCL_Limiter(DeltaScalar_I_Backward,DeltaScalar_I_Forward,DelPhiLim_i)

            !The pressure in Riemann problem
            Scalar_I_Left =subdomain_particle_velocity(i,L)-0.5*DelPhiLim_i
            Scalar_I_Right=subdomain_particle_velocity(i,L)+0.5*DelPhiLim_i


            DeltaScalar_J_Forward =-DOT_PRODUCT(forward_grad_velocity(j,L,:), x_ij(:))
            DeltaScalar_J_Backward=-DOT_PRODUCT(backward_grad_velocity(j,L,:),x_ij(:))

            ! Input variables order: Delta_i_backward,Delta_i_forward,DelPhi_i
            call Upwind_MUSCL_Limiter(DeltaScalar_J_Backward,DeltaScalar_J_Forward,DelPhiLim_j)

            !The pressure in Riemann problem
            Scalar_J_Left =subdomain_particle_velocity(j,L)-0.5*DelPhiLim_j
            Scalar_J_Right=subdomain_particle_velocity(j,L)+0.5*DelPhiLim_j


            !---------------------------------------------------------------------------------
            !Velocity Status for Riemann Solver
            ! I Left and J left
            I_Velocity(L)=subdomain_particle_velocity(i,L)+0.5*(Scalar_I_Left-Scalar_I_Right)
            J_Velocity(L)=subdomain_particle_velocity(j,L)+0.5*(Scalar_J_Left-Scalar_J_Right)

            ! ! I Right and J Right
            ! I_Velocity(L)=subdomain_particle_velocity(i,L)-0.5*(Scalar_I_Left-Scalar_I_Right)
            ! J_Velocity(L)=subdomain_particle_velocity(j,L)-0.5*(Scalar_J_Left-Scalar_J_Right)

            ! ! I Right and J Left
            ! I_Velocity(L)=subdomain_particle_velocity(i,L)-0.5*(Scalar_I_Left-Scalar_I_Right)
            ! J_Velocity(L)=subdomain_particle_velocity(j,L)+0.5*(Scalar_J_Left-Scalar_J_Right)
            !---------------------------------------------------------------------------------

        enddo
        
        IJ_Velocity_Bar=(I_Velocity+J_Velocity)/2.0
        !=====================================================================================

        !--------------------------------Low dissipation Riemann------------------------------

        U_Left =DOT_PRODUCT(I_Velocity,IJ_Direction_Vector)
        U_Right=DOT_PRODUCT(J_Velocity,IJ_Direction_Vector)
        U_Bar=(U_Left+U_Right)/2.0

        average_c=average_c_ij(k)
        average_rho=average_rho_ij(k)

        U_asterisk=U_Bar+0.5*(P_Left-P_Right)/(average_rho*sqrt(square_c_0))

        velocity_asterisk=(U_asterisk-U_Bar)*IJ_Direction_Vector+IJ_Velocity_Bar

        Belta_Riemman=min(3.0*(max((U_Left-U_Right),0.0d0)),average_c)
        P_asterisk=P_Bar+0.5*Belta_Riemman*average_rho*(U_Left-U_Right)
        !-------------------------------------------------------------------------------------


        !The normalized derivates of weight function
        kth_pair_dwdx_i(:)=Normalized_dwdx_i(k,:)
        kth_pair_dwdx_j(:)=Normalized_dwdx_j(k,:)


        !velocity difference times weight derivates
        !-------------------------------------------------------------------------------------
        !For the particles near the free surface
        if (subdomain_near_free_surface_or_not(i)==1) then
            I_velocity_asterisk=( subdomain_particle_velocity(i,:)+subdomain_particle_velocity(j,:) )/2.0d0
            I_P_asterisk=P_Bar
        else
            I_velocity_asterisk=velocity_asterisk
            I_P_asterisk=P_asterisk
        end if

        if (subdomain_near_free_surface_or_not(j)==1) then
            J_velocity_asterisk=( subdomain_particle_velocity(i,:)+subdomain_particle_velocity(j,:) )/2.0d0
            J_P_asterisk=P_Bar
        else
            J_velocity_asterisk=velocity_asterisk
            J_P_asterisk=P_asterisk
        end if

        !write(*,*) I_velocity_asterisk,J_velocity_asterisk

        !XSPH added in the density ratio or not
        !drhodt=sum(m_j*(V_i-V_j))*dwdx   -----the V_j can be replaced by V_i_Riemann

        velocity_difference(:)=subdomain_particle_velocity(i,:)-I_velocity_asterisk+average_velocity(i,:)-average_velocity(j,:)     !Attention: i-j    
        temp_i=DOT_PRODUCT(velocity_difference,kth_pair_dwdx_i)
        velocity_difference(:)=subdomain_particle_velocity(j,:)-J_velocity_asterisk+average_velocity(j,:)-average_velocity(i,:)     !Attention: j-i     
        temp_j=-DOT_PRODUCT(velocity_difference,kth_pair_dwdx_j)                                              !This negative comes for the dwdx
        
        ! velocity_difference(:)=I_Velocity-I_velocity_asterisk    !Attention: i-j     !Attention: i-j
        ! temp_i=DOT_PRODUCT(velocity_difference,kth_pair_dwdx_i)
        ! velocity_difference(:)=J_Velocity-J_velocity_asterisk    !Attention: i-j     !Attention: j-i
        ! temp_j=-DOT_PRODUCT(velocity_difference,kth_pair_dwdx_j)
        !-------------------------------------------------------------------------------------

        !-------------------------------------------------------------------------------------
        !Stand SPH continum equation
        drhodt(i)=drhodt(i)+2*subdomain_particle_rho(i)*subdomain_particle_volume(j)*temp_i
        drhodt(j)=drhodt(j)+2*subdomain_particle_rho(j)*subdomain_particle_volume(i)*temp_j
        !-------------------------------------------------------------------------------------
        
        !*************************************************************************************
        !Stand SPH moment equation
        rho_ij=1.0/(subdomain_particle_rho(i)*subdomain_particle_rho(j))

        do L=1,dim
                
            !(pi+pj)/(rho_i*rho_j)*dwdx(i)
            temp_p_i=-2*I_P_asterisk*Normalized_dwdx_i(k,L)*rho_ij
            temp_p_j=-2*J_P_asterisk*Normalized_dwdx_j(k,L)*rho_ij
        
            !---------------------------------------------------------------------------------
            !Stand SPH acceleration
            internal_acceleration(i,L)=internal_acceleration(i,L)+subdomain_particle_mass(j)*temp_p_i
            internal_acceleration(j,L)=internal_acceleration(j,L)-subdomain_particle_mass(i)*temp_p_j
            !---------------------------------------------------------------------------------
                
        end do

        !*************************************************************************************
        !artifical viscosity acceleration
        average_smooth_length=average_smooth_length_ij(k)
        
        !velocity_difference*position_difference_ij
        v_ij(:)=velocity_difference_ij(k,:)
        x_ij(:)=position_difference_ij(k,:)

        square_distance=DOT_PRODUCT(x_ij,x_ij)
        v_ij_multiply_x_ij=DOT_PRODUCT(v_ij,x_ij)
        distance=sqrt(square_distance)

        !v_ij_multiply_x_ij in Equ.(4.66) in Liu
        if(v_ij_multiply_x_ij<0.0) then
            
            !Equ.(4.67) in Liu
            capital_fai_ij_1=average_smooth_length*v_ij_multiply_x_ij/(square_distance+(fai_coefficient*average_smooth_length)**2)
            
            !Equ.(4.66) in Liu
            capital_pai_ij=(-alpha*average_c*capital_fai_ij_1+beta*(capital_fai_ij_1**2))/average_rho
            
            do L=1,dim
                temp_i=-capital_pai_ij*Normalized_dwdx_i(k,L)  
                temp_j=-capital_pai_ij*Normalized_dwdx_j(k,L)           
                artifical_viscosity_acceleration(i,L)=artifical_viscosity_acceleration(i,L)+subdomain_particle_mass(j)*temp_i
                artifical_viscosity_acceleration(j,L)=artifical_viscosity_acceleration(j,L)-subdomain_particle_mass(i)*temp_j
            end do    
        end if
        !************************************************************************************

        !************************************************************************************
        !laplacian viscos term
        
        !------------------------------------------------------------------------------------
        ! !laplacian viscos term
        ! !Reference: Simulation of near-shore solitary wave mechanics by an incompressible SPH method 
        
        ! RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_i)
        ! laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_mass(j)*4*water_kinematic_viscosity*RijDotProductDerivates/(square_distance*(subdomain_particle_rho(i)+subdomain_particle_rho(j)))*v_ij(:)
        
        ! RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_j)
        ! laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_mass(i)*4*water_kinematic_viscosity*RijDotProductDerivates/(square_distance*(subdomain_particle_rho(i)+subdomain_particle_rho(j)))*v_ij(:)
        !------------------------------------------------------------------------------------

        if ( laplacianTerm==1 ) then

            !--------------------------------------------------------------------------------
            !It is useful for non-free surface flow
            !Reference:  [Morris]  
            RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_i)
            laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_volume(j)*2*water_kinematic_viscosity*RijDotProductDerivates/square_distance*v_ij(:)
          
            RijDotProductDerivates=DOT_PRODUCT(x_ij,kth_pair_dwdx_j)
            laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_volume(i)*2*water_kinematic_viscosity*RijDotProductDerivates/square_distance*v_ij(:)
            !--------------------------------------------------------------------------------

        elseif ( laplacianTerm==0 ) then

            !--------------------------------------------------------------------------------
            !It is useful for free surface flow 
            !Reference: [Theoretical analysis and numerical verification of the consistency of viscous smoothed-particle-hydrodynamics formulations in simulating free-surface flows]
            !           [Adaptive particle refinement and derefinement applied to the smoothed particle hydrodynamics methods]

            laplacian_viscous_acceleration(i,:)=laplacian_viscous_acceleration(i,:)+subdomain_particle_volume(j)*2*(dim+2)*water_kinematic_viscosity*v_ij_multiply_x_ij/square_distance*kth_pair_dwdx_i(:)

            laplacian_viscous_acceleration(j,:)=laplacian_viscous_acceleration(j,:)-subdomain_particle_volume(i)*2*(dim+2)*water_kinematic_viscosity*v_ij_multiply_x_ij/square_distance*kth_pair_dwdx_j(:)
            !--------------------------------------------------------------------------------

        else

            write(*,*) 'The defination of laplacianTerm is not right!(Solver)'
            
        endif
        
        !************************************************************************************

        !************************************************************************************
        !artifical internal acceleration in Delta-SPH(out)
        temp=Delta_alpha*average_smooth_length*average_c*v_ij_multiply_x_ij/square_distance
        do L=1,dim
           artifical_internal_acceleration(i,L)=artifical_internal_acceleration(i,L)+temp*Normalized_dwdx_i(k,L)*water_rho_0/subdomain_particle_rho(i)*subdomain_particle_volume(j)
           artifical_internal_acceleration(j,L)=artifical_internal_acceleration(j,L)-temp*Normalized_dwdx_j(k,L)*water_rho_0/subdomain_particle_rho(j)*subdomain_particle_volume(i)
        end do
        
        ! !************************************************************************************
        ! velocity_difference(:)=velocity_difference_ij(k,:) ! particle i-j
        
        ! average_velocity(i,:)=average_velocity(i,:)-epsilon*subdomain_particle_mass(j)*velocity_difference(:)/(subdomain_particle_rho(i)+subdomain_particle_rho(j))*w(k)
        ! average_velocity(j,:)=average_velocity(j,:)+epsilon*subdomain_particle_mass(i)*velocity_difference(:)/(subdomain_particle_rho(i)+subdomain_particle_rho(j))*w(k)                   
        ! !************************************************************************************

    end do
    !******************************************************************************************

    
    !*******************************Delta-SPH Density coorection*******************************

    !------------------------------------------------------------------------------------------
    !psai=(rho_j-rho_i)-0.5*(grad_rho_i+grad_rho_j)*(r_j-r_i)
    !artifical_rho_correct=2*delta*h*c*psai*(r_j-r_i)*dwdx*V_j/square_distance
    do k=1,pair_number
      
       i=pair_i(k)                        !i index in pair k
       j=pair_j(k)                        !j index in pair k

       !Density difference
       rho_difference=-rho_difference_ij(k)  !Attention: j-i
       
       average_smooth_length=average_smooth_length_ij(k)
       average_c=average_c_ij(k)
       
       !distance and position difference
       position_difference(:)=-position_difference_ij(k,:)                  !Attention: j-i
       
       square_distance=DOT_PRODUCT(position_difference,position_difference)  
       
       !psai 
       sum_grad_rho_ij(:)=grad_rho(i,:)+grad_rho(j,:)
       psai_ij=rho_difference-0.5*DOT_PRODUCT(sum_grad_rho_ij,position_difference)

       !Temp Value
       temp_i=DOT_PRODUCT(position_difference,Normalized_dwdx_i(k,:))/square_distance
       temp_j=DOT_PRODUCT(position_difference,Normalized_dwdx_j(k,:))/square_distance
       
       artifical_rho_correct(i)=artifical_rho_correct(i)+2*delta*average_smooth_length*average_c*psai_ij*temp_i*subdomain_particle_volume(j)
       artifical_rho_correct(j)=artifical_rho_correct(j)-2*delta*average_smooth_length*average_c*psai_ij*temp_j*subdomain_particle_volume(i)
       
    end do
    !------------------------------------------------------------------------------------------

    !******************************************************************************************


 end subroutine Solver_Normalized_MUSCL_Upwind_Linear_Riemann_Delta