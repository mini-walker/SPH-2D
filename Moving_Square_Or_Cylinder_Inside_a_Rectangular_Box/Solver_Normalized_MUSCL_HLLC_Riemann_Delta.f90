!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: Solver_Normalized_MUSCL_HLLC_Riemann_Delta
!
!  PURPOSE: Compute iterator information for time integration with Normalized kernel function Delta SPH and MUSCL Riemann solver
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

 subroutine Solver_Normalized_MUSCL_HLLC_Riemann_Delta()
 
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
    real(kind=8)::P_Left,P_right,rho_Left,rho_right,U_Left,U_right,c_Left,c_Right
    real(kind=8)::P_Bar,U_Bar                                                      !Average pressure and velocity
    real(kind=8),dimension(dim)::I_Velocity,J_Velocity,IJ_Velocity_Bar
    real(kind=8),dimension(dim)::IJ_Direction_Vector
    real(kind=8)::U_asterisk,P_asterisk,I_P_asterisk,J_P_asterisk,Rho_asterisk,Rho_asterisk_Left,Rho_asterisk_Right
    real(kind=8),dimension(dim)::velocity_asterisk,I_velocity_asterisk,J_velocity_asterisk
    real(kind=8)::Belta_Riemman,Q_L,Q_R,Q_factor_L,Q_factor_R
    real(kind=8)::S_Left,S_Right,S_asterisk,x_over_t
    real(kind=8),dimension(dim)::Velocity_Left,Velocity_Right,Velocity_HLLC,Rho_Velocity_HLLC
    real(kind=8)::sin_theta,cos_theta
    real(kind=8)::Rho_HLLC,P_HLLC
    real(kind=8),dimension(dim)::Velocity_Bar_0
    real(kind=8),dimension(dim,dim)::Transfer_Matrix                               !Transfer Matrix from global to local

    !MUSCL
    real(kind=8)::DelPhi_i,DelPhi_j,DelPhi_ij,DelPhiLim_i,DelPhiLim_j

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
    

    !***************************************************************************************** 
    
    !-----------------------------------------------------------------------------------------
    !Calculate the gradience of density, pressure and velocity derivates

    grad_rho=0.0d0
    grad_pressure=0.0d0
    grad_velocity=0.0d0

    do k=1,pair_number
       
       i=pair_i(k)                                            !i index in pair k
       j=pair_j(k)                                            !j index in pair k

       !Density difference
       rho_difference=-rho_difference_ij(k)                    !Attention: j-i

       !Density gradience
       grad_rho(i,:)=grad_rho(i,:)+rho_difference*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)
       grad_rho(j,:)=grad_rho(j,:)+rho_difference*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)

       !velocity difference
       velocity_difference(:)=-velocity_difference_ij(k,:)     !Attention: j-i

       do L=1,dim

           !velocity gradience
           grad_velocity(i,L,:)=grad_velocity(i,L,:)+velocity_difference(L)*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)
           grad_velocity(j,L,:)=grad_velocity(j,L,:)+velocity_difference(L)*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)

       enddo

       !Pressure difference
       pressure_difference=subdomain_particle_press(j)-subdomain_particle_press(i)  !Attention: j-i

       !Pressure gradience
       grad_pressure(i,:)=grad_pressure(i,:)+pressure_difference*subdomain_particle_volume(j)*Normalized_dwdx_i(k,:)
       grad_pressure(j,:)=grad_pressure(j,:)+pressure_difference*subdomain_particle_volume(i)*Normalized_dwdx_j(k,:)

    end do  
    !-----------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------
    !Exchange the grad_rho form the neighbour processor

    call MPI_RealVector_MessageExchange(grad_rho)           ! Exchange the gradience of the density
    call MPI_RealTensor_MessageExchange(grad_velocity)      ! Exchange the gradience of the horizontal velocity
    call MPI_RealVector_MessageExchange(grad_pressure)      ! Exchange the gradience of the pressure
    !-----------------------------------------------------------------------------------------


    !*****************************************************************************************


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


    !***********************************Stand SPH calculation*********************************
    do k=1,pair_number
        
        i=pair_i(k)                                     !particle i index in kth pair
        j=pair_j(k)                                     !particle j index in kth pair

        !Get the direction vector for from i to j
        x_ij(:)=-position_difference_ij(k,:)            !j-1 
        distance=sqrt(DOT_PRODUCT(x_ij,x_ij))
        IJ_Direction_Vector=x_ij/distance               !R_ji

        !=====================================================================================
        !MUSCL Scheme
        !Reconstrute the left and right value in the Riemann solver

        !----------Pressure---------
        DelPhi_i= 0.5*DOT_PRODUCT(grad_pressure(i,:),x_ij(:))
        DelPhi_j=-0.5*DOT_PRODUCT(grad_pressure(j,:),x_ij(:))
        DelPhi_ij=subdomain_particle_press(j)-subdomain_particle_press(i) !j-1 

        call MUSCL_Limiter(DelPhi_i,DelPhi_j,DelPhi_ij,DelPhiLim_i,DelPhiLim_j)

        !The pressure in Riemann problem
        P_Left =subdomain_particle_press(i)+DelPhiLim_i
        P_Right=subdomain_particle_press(j)-DelPhiLim_j
        P_Bar=(P_Left+P_Right)/2.0

        !----------Velocity---------
        velocity_difference(:)=-velocity_difference_ij(k,:)               !j-1 
        
        !Loop for 3 direction
        do L=1,dim

            DelPhi_i= 0.5*DOT_PRODUCT(grad_velocity(i,L,:),x_ij(:))
            DelPhi_j=-0.5*DOT_PRODUCT(grad_velocity(j,L,:),x_ij(:))
            DelPhi_ij=velocity_difference(L) 

            call MUSCL_Limiter(DelPhi_i,DelPhi_j,DelPhi_ij,DelPhiLim_i,DelPhiLim_j)

            !The horizotial velocity in Riemann problem
            I_Velocity(L)=subdomain_particle_velocity(i,L)+DelPhiLim_i
            J_Velocity(L)=subdomain_particle_velocity(j,L)-DelPhiLim_j

        enddo
        
        IJ_Velocity_Bar=(I_Velocity+J_Velocity)/2.0
        !=====================================================================================

        !=====================================================================================
        !HLLC Riemann solver

        !Assign the transfer matrix (From global to local)
        !Reference: A improved MUSCL treatment for the SPH-ALE method: comparison with the standard SPH method for the jet impingement case
        !                [xn  yn  zn]
        !Transfer_Matrix=[xm  ym  zm]
        !                [xo  yo  zo]
        ![xn  yn  zn] is the unit vector from i to j
        ![xm  ym  zm] and [xo  yo  zo] are the components of mij and oij unit vectors, respectively.
        !These vector are normal to each other and normal to the nij vector
        Transfer_Matrix(1,1)= IJ_Direction_Vector(1) ; Transfer_Matrix(1,2)=IJ_Direction_Vector(2) ;
        Transfer_Matrix(2,1)=-IJ_Direction_Vector(2) ; Transfer_Matrix(2,2)=IJ_Direction_Vector(1) ;

        Velocity_Left =matmul(Transfer_Matrix,I_Velocity)
        Velocity_Right=matmul(Transfer_Matrix,J_Velocity)
        Velocity_Bar_0=matmul(Transfer_Matrix,IJ_Velocity_Bar)

        rho_Left=subdomain_particle_rho(i)
        rho_Right=subdomain_particle_rho(j)

        c_Left=subdomain_particle_c(i)
        c_Right=subdomain_particle_c(j)

        average_c=average_c_ij(k)
        average_rho=average_rho_ij(k)
        !-------------------------------------------------------------------------------------

        !Horizontial velocity star
        U_asterisk=(Velocity_Left(1)*c_Left*rho_Left+Velocity_Right(1)*c_Right*rho_Right+P_Left-P_Right)/(rho_Left*c_Left+rho_Right*c_Right)

        !Pressure star
        P_asterisk=(P_Right*c_Left*rho_Left+P_Left*c_Right*rho_Right-c_Left*rho_Left*c_Right*rho_Right*(Velocity_Right(1)-Velocity_Left(1)))/(rho_Left*c_Left+rho_Right*c_Right)

        Rho_asterisk=(P_asterisk-Background_Pressure)/square_c_0+water_rho_0

        Rho_asterisk_Left=rho_Left*((P_asterisk/P_Left)**(0.142857142))

        Rho_asterisk_Right=rho_Right*((P_asterisk/P_Right)**(0.142857142))

        !-------------------------------------------------------------------------------------
        !HLLC approximation
        !--- Left Shock Speed ---
        if ( P_asterisk <= P_Left ) then          !Rarefaction

          Q_L=1.0d0

        else                                      !Shock

            if (rho_Left>Rho_asterisk) then

                Q_L=sqrt( (Rho_asterisk/rho_Left)*(P_asterisk-P_Left)/(rho_Left-Rho_asterisk) )/c_Left

            else

                Q_L=1.0d0

            endif

        endif

        S_Left=Velocity_Left(1)-c_Left*Q_L

        !--- Right Shock Speed ---
        if ( P_asterisk <= P_Right ) then         !Rarefaction

          Q_R=1.0d0

        else                                      !Shock

            if (rho_Right>Rho_asterisk) then

                Q_R=sqrt( (Rho_asterisk/rho_Right)*(P_asterisk-P_Right)/(rho_Right-Rho_asterisk) )/c_Right

            else

                Q_R=1.0d0

            endif

        endif

        S_Right=Velocity_Right(1)+c_Right*Q_R

        !--- Contact discontinuity Shock Speed ---
        S_asterisk=U_asterisk

        x_over_t=Velocity_Bar_0(1)                       !This is very important

        if ( S_Left > x_over_t ) then

          Rho_HLLC=rho_Left
          Velocity_HLLC=Velocity_Left

        elseif ( S_asterisk>= x_over_t .and. x_over_t >= S_Left ) then

          Q_factor_L = rho_Left*((S_Left - Velocity_Left(1))/(S_Left - S_asterisk))
          
          Rho_HLLC=Q_factor_L

          Rho_Velocity_HLLC(1)=S_asterisk*Q_factor_L
          Rho_Velocity_HLLC(2)=Velocity_Left(2)*Q_factor_L

          Velocity_HLLC=Rho_Velocity_HLLC/Rho_HLLC

        elseif ( S_Right >= x_over_t .and. x_over_t >= S_asterisk ) then

          Q_factor_R = rho_Right*((S_Right - Velocity_Right(1))/(S_Right - S_asterisk))

          Rho_HLLC=Q_factor_R

          Rho_Velocity_HLLC(1)=S_asterisk*Q_factor_R
          Rho_Velocity_HLLC(2)=Velocity_Right(2)*Q_factor_R

          Velocity_HLLC=Rho_Velocity_HLLC/Rho_HLLC

        else

          Rho_HLLC=rho_Right
          Velocity_HLLC=Velocity_Right

        endif

        !-------------------------------------------------------------------------------------
        !Re-assigned the star value
        velocity_asterisk =matmul(TRANSPOSE(Transfer_Matrix),Velocity_HLLC)
        P_asterisk=square_c_0*(Rho_HLLC-water_rho_0)+Background_Pressure

        !=====================================================================================

        !The normalized derivates of weight function
        kth_pair_dwdx_i(:)=Normalized_dwdx_i(k,:)
        kth_pair_dwdx_j(:)=Normalized_dwdx_j(k,:)


        !velocity difference times weight derivates
        !-------------------------------------------------------------------------------------
        !For the particles near the free surface
        if (subdomain_near_free_surface_or_not(i)==1) then
            I_velocity_asterisk=IJ_Velocity_Bar
            I_P_asterisk=P_Bar
        else
            I_velocity_asterisk=velocity_asterisk
            I_P_asterisk=P_asterisk
        end if

        if (subdomain_near_free_surface_or_not(j)==1) then
            J_velocity_asterisk=IJ_Velocity_Bar
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


 end subroutine Solver_Normalized_MUSCL_HLLC_Riemann_Delta






