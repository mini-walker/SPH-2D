!**********************************************************************************************************
!  SUBROUTINE: Compute_kernel
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
!**********************************************************************************************************

subroutine compute_kernel(Kernel_function_model,&                       ! Kernel function model
                          dim,&                                         ! Case Dimension
                          distance,&                                    ! Distance(in)                
                          position_difference,&                         ! Position difference(in)
                          average_smooth_length,&                       ! Average smooth length(in)
                          w,&                                           ! The weight function value for current pair(out)
                          temp_dwdx&                                    ! The derivates of weight function for current pair(out)
                          )

    implicit none

    !======================================================================================================
    ! Variables from superior-subroutine
    integer,intent(in)::Kernel_function_model                           ! Kernel function model
    integer,intent(in)::dim                                             ! Case Dimension
    real(kind=8),dimension(dim),intent(in)::position_difference         ! Position difference(i-j),1-x,2-y,3-z
    real(kind=8),intent(in)::distance                                   ! Distance between particle i and j
    real(kind=8),intent(in)::average_smooth_length                      ! Average smooth length of particle i and j
    
    real(kind=8),intent(out)::w                                         ! The weight function value for current pair(out)
    real(kind=8),dimension(dim),intent(out)::temp_dwdx                  ! The derivates of weight function for current pair(out)

    
    ! Variables in local subroutine
    real(kind=8),parameter::PI=3.14159265358                            ! Constant parameter---PI

    integer::i,j,k,l                                                    ! Variables for looping
    real(kind=8)::factor                                                ! Factor  P90
    real(kind=8)::q                                                     ! Ratio: r/h (Distance/Average smooth length)

    !======================================================================================================

    ! Body of subroutine compute_kernel
    
    !Initialize the value
    q=distance/average_smooth_length

    !Quintic spline function
    Kernel_function_model_if: if ( Kernel_function_model==1 ) then

        !--------------------------------------------------------------------------------------------------
        !Factors
        if (dim==1) then
            factor=1/(120.0*average_smooth_length)
        elseif (dim==2) then
            factor=7.0/(478.0*PI*(average_smooth_length*average_smooth_length))
        elseif (dim==3) then
            factor=3/(359.0*PI*average_smooth_length**3)
        endif

        
        !Calculate the derivates and value of weight function (w and temp_dwdx)
        if(0.0<=q .and. q<1.0) then
           
            w=factor*((3.0-q)**5-6.0*(2.0-q)**5+15.0*(1.0-q)**5)

            temp_dwdx(:)=factor*(-5.0*(3.0-q)**4+30.0*(2.0-q)**4-75*(1.0-q)**4)/average_smooth_length*(position_difference(:)/distance)
           
        else if(1.0<=q .and. q<2.0) then
               
            w=factor*((3.0-q)**5-6.0*(2.0-q)**5)

            temp_dwdx(:)=factor*(-5.0*(3.0-q)**4+30.0*(2.0-q)**4)/average_smooth_length*(position_difference(:)/distance)

        else if(2.0<=q .and. q<=3.0) then
           
            w=factor*((3.0-q)**5)

            temp_dwdx(:)=factor*(-5.0*(3.0-q)**4)/average_smooth_length*(position_difference(:)/distance)

        else
            w=0.0d0
            temp_dwdx(:)=0.0d0
        end if
        !--------------------------------------------------------------------------------------------------

    elseif ( Kernel_function_model==2 ) then

        !--------------------------------------------------------------------------------------------------
        !Cubic spline function
        !Factors
        if (dim==1) then
            factor=1/(average_smooth_length)
        elseif (dim==2) then
            factor=15.0/(7.0*PI*(average_smooth_length*average_smooth_length))
        elseif (dim==3) then
            factor=3.0/(2.0*PI*average_smooth_length**3)
        endif

        !Calculate the derivates and value of weight function (w and temp_dwdx)
        if(0.0<=q .and. q<1.0) then
           
            w=factor*(2.0/3.0-q**2+q**3/2.0)
            temp_dwdx(:)=factor*(-2+3.0/2.0*q)/average_smooth_length**2*position_difference(:)
           
        else if(1.0<=q .and. q<=2.0) then
               
            w=factor*((2-q)**3/6.0)
            temp_dwdx(:)=-factor*((2-q)**2/2.0)/average_smooth_length*(position_difference(:)/distance)
        
        else
            w=0.0d0
            temp_dwdx(:)=0.0d0
        end if
        !--------------------------------------------------------------------------------------------------

    elseif ( Kernel_function_model==3 ) then

        !--------------------------------------------------------------------------------------------------
        !Quadratic smooth function
        !Factors
        if (dim==1) then
            factor=1/(average_smooth_length)
        elseif (dim==2) then
            factor=2.0/(PI*(average_smooth_length*average_smooth_length))
        elseif (dim==3) then
            factor=5.0/(4.0*PI*average_smooth_length**3)
        endif

        !Calculate the derivates and value of weight function (w and temp_dwdx)
        if(0.0<=q .and. q<2.0) then
           
            w=factor*(3.0/16.0*q**2 - 3.0/4.0*q + 3.0/4.0)
            temp_dwdx(:)=factor*(3.0/8.0*q - 3.0/4.0)*position_difference(:)/average_smooth_length**2

        else
            w=0.0d0
            temp_dwdx(:)=0.0d0
        end if
        !--------------------------------------------------------------------------------------------------

    elseif ( Kernel_function_model==4 ) then

        !--------------------------------------------------------------------------------------------------
        !Improved Gauss Kernel Fcuntion

        !factor
        factor=1.0/(PI**(dim/2.0)*(average_smooth_length**dim)*(1-10*exp(-9.0)))

        !w and dwdx
        if(distance<=3.0*average_smooth_length) then
            
            !Kernel function value 
            w=(exp(-q**2)-exp(-9.0))*factor;

            !Derivates of kernel function
            temp_dwdx(:)=(exp(-q**2)*(-2*q))*position_difference(:)/(distance*average_smooth_length)*factor;

        else
            w=0.0d0
            temp_dwdx=0.0d0
        end if
        !--------------------------------------------------------------------------------------------------

    elseif ( Kernel_function_model==5 ) then

        !--------------------------------------------------------------------------------------------------
        !For wendland kernel function,we recommand for C2 kernel function;
        !C4 and C6 require more near particles and the viscous is C6>C4>C2; SO C2 is better for wave CSPH_MPI_Parallel_wave_maker_2D
        !h=1.25dx is more stable, however h=2.0dx will make the prussure more smooth 
        !Wendland C2 kernel function
        if (dim==2) then
            factor=7.0/(4*PI*average_smooth_length**2)
        elseif (dim==3) then
            factor=21.0/(16*PI*average_smooth_length**3)
        endif
          
             
        !w and dwdx
        if(distance<=2.0*average_smooth_length) then
            
            !Kernel function value 
            w=factor*(1+2*q)*(1-q/2.0)**4
            
            !Derivates of kernel function
            temp_dwdx(:)=factor*(-5*q*(1-q/2)**3)*position_difference(:)/(distance*average_smooth_length)

        else
            w=0.0d0
            temp_dwdx=0.0d0
        end if
        !--------------------------------------------------------------------------------------------------

    elseif ( Kernel_function_model==6 ) then

        !--------------------------------------------------------------------------------------------------
        !Wendland C4 kernel function
        if (dim==2) then
            factor=9.0/(4*PI*average_smooth_length**2)
        elseif (dim==3) then
            factor=495/(256*PI*average_smooth_length**3)
        endif
          
        !w and dwdx
        if(distance<=2.0*average_smooth_length) then
            
            !Kernel function value 
            w=factor*(1+3*q+35.0/12.0*q**2)*(1-q/2.0)**6

            !Derivates of kernel function
            temp_dwdx(:)=factor*(-3*(1-q/2.0)**5*(35/12.0*q**2+3*q+1)+(1-q/2.0)**6*(35/6.0*q+3))*position_difference(:)/(distance*average_smooth_length)

        else
            w=0.0d0
            temp_dwdx=0.0d0
        end if
        !--------------------------------------------------------------------------------------------------

    elseif ( Kernel_function_model==7 ) then

        !--------------------------------------------------------------------------------------------------
        !Wendland C6 kernel function
        if (dim==2) then
            factor=78/(28*PI*average_smooth_length**2)
        elseif (dim==3) then
            factor=1365/(512*PI*average_smooth_length**3)
        endif
             
        !w and dwdx
        if(distance<=2.0*average_smooth_length) then
            
            !Kernel function value 
            w=factor*(1+4*q+6.25*q**2+4*q**3)*(1-q/2.0)**8

            !Derivates of kernel function
            temp_dwdx(:)=factor*(-4*(1-q/2.0)**7*(4*q**3+6.25*q**2+4*q+1)+(1-q/2.0)**8*(12*q**2+12.5*q+4))*position_difference(:)/(distance*average_smooth_length)

        else
            w=0.0d0
            temp_dwdx=0.0d0
        end if
        !--------------------------------------------------------------------------------------------------

    else

        write(*,'(A)') " Value of 'Kernel_function_model' is not right! ( Compute_Kernel: value 1--7 ) "

    endif Kernel_function_model_if

    ! For the particle it own
    if ( q<=1.0E-8 ) then
        temp_dwdx = 0.0d0
    endif
   
end subroutine compute_kernel 
    