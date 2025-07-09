!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: MLS_Interpolation_In_Subdomain_Domain
!
!  Purpose: interpolation calculation for particle in subdomain domian 
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
    
subroutine MLS_Interpolation_In_Subdomain_Domain(MLSOrder,Interpolate_subdomain_particle_position,ture_effect_number,ture_effect_name,W_kernel,Fai)

    use function_module
    use information_module
    use Public_variable_module

    implicit none

    ! Variables
    ! Variables from mother subroutine
    integer,intent(in)::MLSOrder                                                           !MLSOrder
    real(kind=8),dimension(dim),intent(in)::Interpolate_subdomain_particle_position                  !Interpolation particle position
    integer,intent(in)::ture_effect_number                                                 !ture effect particle number
    integer,dimension(NumberInEachGridPrediction),intent(in)::ture_effect_name             !ture effect particle index
    real(kind=8),dimension(NumberInEachGridPrediction),intent(in)::W_kernel                !weight function value
    real(kind=8),dimension(NumberInEachGridPrediction),intent(out)::Fai                    !Output fai value
    
    ! Variables in subroutine
    integer::i,j,k,l,m                                                                     !Variables for looping
    
    ! MLS Variables
    real(kind=8),dimension(NumberInEachGridPrediction,NumberInEachGridPrediction)::w_matrix! Weight function value matrix 
    real(kind=8),dimension(MLSOrder,NumberInEachGridPrediction)::P                         ! Basis P  
    real(kind=8),dimension(MLSOrder,MLSOrder)::Modified_matrix                             ! Modified identity matrix 
    real(kind=8),dimension(MLSOrder,1)::P_i                                                ! P matrix of i (Dimension of P is : MLSOrder*1)
    real(kind=8),dimension(MLSOrder,1)::Q_i                                                ! Q matrix of i (Dimension of P is : MLSOrder*1)
    real(kind=8),dimension(MLSOrder,NumberInEachGridPrediction)::B                         ! Basis B (Dimension of B is : MLSOrder*k)  
    real(kind=8),dimension(MLSOrder,MLSOrder)::A,A_inversion                               ! A matrix (Dimension of A and A_inversion is : MLSOrder*MLSOrder)
    real(kind=8),dimension(MLSOrder,NumberInEachGridPrediction)::A_inversion_times_B
    real(kind=8),dimension(1,1)::Temp_Fai
    real(kind=8)::sum_w_kernel
    
    ! Updated_Spherd Variables
    real(kind=8)::miu=4
    real(kind=8)::shperd_R
    real(kind=8)::shperd_Belta
    real(kind=8),dimension(NumberInEachGridPrediction)::shperd_w
    real(kind=8)::distance
    real(kind=8),dimension(dim)::position_difference                                       
    
    integer::check_error
    
    ! Body of subroutine MLS_Interpolation

    ! Weight function value matrix
    w_matrix=0.0d0                                       ! (Dimension of W is : k*k)
    do k=1,ture_effect_number 
        w_matrix(k,k)=W_kernel(k)
    end do
        
    ! Assign the Modified_matrix
    Modified_matrix=0.0d0
    do k=1,MLSOrder
        Modified_matrix(k,k)=slighted_modified_alpha 
    end do

    !Initialization
    P=0.0d0
    B=0.0d0
    A=0.0d0
    
    check_error=1
    
    if(MLSOrder==2) then
                
        !************************************************************************************
        ! Updated Spherd Initialization
        shperd_w=0.0d0
        sum_w_kernel=0.0d0
        
        ! Support radius for spherd interpolation method
        shperd_R=kernel_scale*smooth_length

        do k=1,ture_effect_number
    
            ! effect particle index
            j=ture_effect_name(k)
     
            ! Calculate position difference and distance
            position_difference(:)=Interpolate_subdomain_particle_position(:)-subdomain_particle_position(j,:)
            distance=sqrt(DOT_PRODUCT(position_difference,position_difference))
            
            if (distance<=shperd_R) then
                shperd_Belta=exp(-distance**2/(1.0-(distance/shperd_R)**2))
            else
                shperd_Belta=0.0d0
            end if
            
            shperd_w(k)=(shperd_Belta/distance)**miu

            sum_w_kernel=sum_w_kernel+shperd_w(k)
 
        end do

        do k=1,ture_effect_number
             Fai(k)=shperd_w(k)/sum_w_kernel
        end do
        !************************************************************************************
        
    else
        
        !************************************************************************************
        ! MLS interpolation
        ! Evaluate basis P, B Matrix and their derivates  (m is the order of MLS)
        do k=1,ture_effect_number
    
            ! effect particle index
            j=ture_effect_name(k)
        
            if(MLSOrder==1) then
    
                P(1,k)=1                                       ! P_i=(1)            (Dimension of P is : k*MLSOrder)
    
            else if(MLSOrder==3) then
                                
                P(1,k)=1                                       ! P_i=(1,x_j/h,y_j/h)
                P(2,k)=(subdomain_particle_position(j,1)-Interpolate_subdomain_particle_position(1))/smooth_length
                P(3,k)=(subdomain_particle_position(j,2)-Interpolate_subdomain_particle_position(2))/smooth_length
        
            else if(MLSOrder==6) then
    
                P(1,k)=1                                       ! P_i=(1,x_j/h,y_j/h,(x_j/h)^2,x_j/h*y_j/h,(y_j/h)^2)
                P(2,k)=(subdomain_particle_position(j,1)-Interpolate_subdomain_particle_position(1))/smooth_length
                P(3,k)=(subdomain_particle_position(j,2)-Interpolate_subdomain_particle_position(2))/smooth_length
                P(4,k)=((subdomain_particle_position(j,1)-Interpolate_subdomain_particle_position(1))/smooth_length)**2
                P(5,k)=(subdomain_particle_position(j,1)-Interpolate_subdomain_particle_position(1))/smooth_length*(subdomain_particle_position(j,2)-Interpolate_subdomain_particle_position(2))/smooth_length
                P(6,k)=((subdomain_particle_position(j,2)-Interpolate_subdomain_particle_position(2))/smooth_length)**2
    
            else
                write(*,*) "The order of MLS input is not right , Please ckeck 'MLSOrder' value !" 
            end if
    
        end do
        !************************************************************************************************************
         
        !************************************************************************************************************
        ! Calculate P_i and its derivates
        if(MLSOrder==1) then
        
            P_i(1,1)=1                                       ! P_i=(1)            (Dimension of P is : MLSOrder*1)
        
        else if(MLSOrder==3) then
                                
            P_i(1,1)=1                                       ! P_i=(1,x_i,y_i)
            P_i(2,1)=0.0d0
            P_i(3,1)=0.0d0
        
        else if(MLSOrder==6) then
        
            P_i(1,1)=1                                       ! P_i=(1,x_i,y_i,x_i^2,x_i*y_i,y_i^2)
            P_i(2,1)=0.0d0
            P_i(3,1)=0.0d0
            P_i(4,1)=0.0d0
            P_i(5,1)=0.0d0
            P_i(6,1)=0.0d0
        
        else
            write(*,*) "The order of MLS input is not right , Please ckeck 'MLSOrder' value !" 
        end if
        !************************************************************************************************************
       
        !************************************************************************************************************
        ! Modified the P and W matrix (Makesure non-singal happens)
        if(MLSOrder/=1) then
        
            ! Modified P: New_P=[P Modified_matrix]   (dimension: m*(k+m) )
            ! Modified W: New_W=[ W               0           ]
            !                   [ 0       Modified_matrix     ]                      (dimension: (k+m)*(k+m) )
            do m=1,MLSOrder
                do L=1,MLSOrder
                    k=ture_effect_number+L
                    P(m,k)=Modified_matrix(m,L)                     ! Comble P and Modified_matrix( as the add matrix is constant so its derivated keep same)
                end do
            end do 
        
            do L=1,MLSOrder
                 k=ture_effect_number+L
                 W_matrix(k,k)=Modified_matrix(L,L)                 ! Comble P and Modified_matrix( as the add matrix is constant so its derivated keep same)
            end do
        
        end if
        !************************************************************************************************************
    
        !************************************************************************************************************
        ! Calculate B,  A_inversion
        ! Assign the B=P*W                
        B=matmul(P,W_matrix)                          ! (Dimension of B is : MLSOrder*k)
            
        ! Evaluate matrices A=P*W*transpose(P) and its derivatives
        A=matmul(B,transpose(P))                      ! (Dimension of A is : MLSOrder*MLSOrder)
    
        check_error=1
        if(MLSOrder==1) then
            
            A_inversion(1,1)=1.0/A(1,1)
    
        else if(MLSOrder==3) then
            
            !------------------------------------------------------------------------
            call matrix_inversion_less_three(MLSOrder,A,A_inversion,check_error)
            !------------------------------------------------------------------------

            ! !------------------------------------------------------------------------
            ! A_inversion=A
            ! call BRINV(A_inversion,MLSOrder,check_error)
            ! !------------------------------------------------------------------------

            if(check_error==0) then
                write(*,*) "The A_inversion compute is wrong! (Interpolation)"
            end if
        
        else if(MLSOrder==6) then 
                    
            A_inversion=A                                          !Initialized the A inversion matrix
            !Calculate the inverse  matrix of A
            !call matrix_inversion_less_six(MLSOrder,A,A_inversion,check_error)
            call BRINV(A_inversion,MLSOrder,check_error)
            
            if(check_error==0) then
                write(*,*) "The A_inversion compute is wrong! (Interpolation)"
            end if
    
        end if
    
        ! Check inversion calculation is right or not
        if(check_error/=0) then
            ! ***************************************************************************************************************************
            ! Calculate A_inversion times B  and its derivates (dimension: MLSOrder*k)
            A_inversion_times_B=matmul(A_inversion,B)
     
            Q_i=P_i
            do k=1,ture_effect_number
    
                ! Calculate P_i and B_i   
                Temp_Fai=0.0d0
                do m=1,MLSOrder
                
                    Temp_Fai=Temp_Fai+Q_i(m,1)*A_inversion_times_B(m,k)
              
                end do
        
                Fai(k)=Temp_Fai(1,1)

            end do
            !***********************************************************************************************************************

        ! For the case, the inverse oricedure is not right
        elseif(check_error==0) then

            sum_w_kernel=0.0d0
            do k=1,ture_effect_number
                sum_w_kernel=sum_w_kernel+W_kernel(k)  
            end do
 
            do k=1,ture_effect_number
                Fai(k)=W_kernel(k)/sum_w_kernel
            end do
        
        end if

        !For high order calculation, when the near particle is not enough, decrease the MLS Order
        if(MLSOrder==6 .and. ture_effect_number<=8) then

            sum_w_kernel=0.0d0
            do k=1,ture_effect_number
                sum_w_kernel=sum_w_kernel+W_kernel(k)  
            end do
 
            do k=1,ture_effect_number
                Fai(k)=W_kernel(k)/sum_w_kernel
            end do
        
        end if
    
    end if
    

end subroutine MLS_Interpolation_In_Subdomain_Domain