!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: compute_kernel
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

subroutine compute_kernel(dim,&                            !维度
                          distance,&                       !距离(in)                
                          position_difference,&            !坐标差(in)
                          average_smooth_length,&          !平均光滑长度(in)
                          w,&                              !第pair_number对粒子对的核函数值(out)
                          temp_dwdx&                       !第pair_number对粒子对的核函数对x/y/z的偏导值(out)
                          )

implicit none

    ! Variables
    real(kind=8),parameter::PI=3.14159265358                            !常数PI
    
    !传递变量
    integer,intent(in)::dim                                             !维度
    real(kind=8),intent(out)::w                                         !第k个粒子对中核函数值
    real(kind=8),dimension(dim),intent(out)::temp_dwdx                  !核函数的偏导数值中间存储值
    real(kind=8),dimension(dim),intent(in)::position_difference         !坐标的差(i-j),1-x,2-y,3-z
    real(kind=8),intent(in)::distance                                   !i，j粒子点的距离
    real(kind=8),intent(in)::average_smooth_length                      !i，j粒子点的平均光滑长度
    
    !子程序内部使用变量
    integer::i,j,k,l                                                !循环控制变量
    real(kind=8)::factor                                            !比例因子  P90
    real(kind=8)::q                                                 !r/h 距离和光滑长度的比值
    
    ! Body of subroutine compute_kernel
    
    !初始化数据
    q=distance/average_smooth_length

    w=0.0d0
    temp_dwdx=0.0d0

    
    !!五次分段
    !!计算比例因子
    !!factor=1/(120*average_smooth_length)
    !factor=7.0/(478.0*PI*(average_smooth_length*average_smooth_length))
    !!factor=3/(359*PI*average_smooth_length**3)
    !
    !!计算w和temp_dwdt
    !if(0.0<=q .and. q<1.0) then
    !    
    !    w=factor*((3.0-q)**5-6.0*(2.0-q)**5+15.0*(1.0-q)**5)
    !    do i=1,dim
    !        temp_dwdx(i)=factor*(-5.0*(3.0-q)**4+30.0*(2.0-q)**4-75*(1.0-q)**4)/average_smooth_length*(position_difference(i)/distance)
    !    end do
    !    
    !else if(1.0<=q .and. q<2.0) then
    !        
    !    w=factor*((3.0-q)**5-6.0*(2.0-q)**5)
    !    do i=1,dim
    !        temp_dwdx(i)=factor*(-5.0*(3.0-q)**4+30.0*(2.0-q)**4)/average_smooth_length*(position_difference(i)/distance)
    !    end do
    ! 
    !else if(2.0<=q .and. q<=3.0) then
    !    
    !    w=factor*((3.0-q)**5)
    !    do i=1,dim
    !        temp_dwdx(i)=factor*(-5.0*(3.0-q)**4)/average_smooth_length*(position_difference(i)/distance)
    !    end do
    !    
    !else
    !    w=0.0
    !    do i=1,dim
    !        temp_dwdx(i)=0.0
    !    end do
    !end if
    !
    !!分段三次样条
    !!计算比例因子
    !!factor=1/(average_smooth_length)
    !factor=15.0/(7.0*PI*(average_smooth_length*average_smooth_length))
    !!factor=3.0/(2.0*PI*average_smooth_length**3)
    !
    !!计算w和temp_dwdt
    !if(0.0<=q .and. q<1.0) then
    !    
    !    w=factor*(2.0/3.0-q**2+q**3/2.0)
    !    do i=1,dim
    !        temp_dwdx(i)=factor*(-2+3.0/2.0*q)/average_smooth_length**2*position_difference(i)
    !    end do
    !    
    !else if(1.0<=q .and. q<=2.0) then
    !        
    !    w=factor*((2-q)**3/6.0)
    !    do i=1,dim
    !        temp_dwdx(i)=-factor*((2-q)**2/2.0)/average_smooth_length*(position_difference(i)/distance)
    !    end do
    !
    !else
    !    w=0.0
    !    do i=1,dim
    !        temp_dwdx(i)=0.0
    !    end do
    !end if
        
    !--------------------------------------------------------------------------------------------------
    !Improved Gauss Kernel Fcuntion
    
    !factor
    factor=1.0/(PI**(dim/2.0)*(average_smooth_length**dim)*(1-10*exp(-9.0)))

    !w and dwdx
    if(distance<=3.0*average_smooth_length) then
        
        !Kernel function value 
        w=(exp(-q**2)-exp(-9.0))*factor;
    
        !Derivates of kernel function
        do i=1,dim
            temp_dwdx(i)=(exp(-q**2)*(-2*q))*position_difference(i)/(distance*average_smooth_length)*factor;
        end do

    else

        w=0.0
        temp_dwdx=0.0
        
    end if
    !--------------------------------------------------------------------------------------------------

    ! !--------------------------------------------------------------------------------------------------
    ! !For wendland kernel function,we recommand for C2 kernel function;
    ! !C4 and C6 require more near particles and the viscous is C6>C4>C2; SO C2 is better for wave CSPH_MPI_Parallel_wave_maker_2D
    ! !h=1.25dx is more stable, however h=2.0dx will make the prussure more smooth 
    ! !Wendland C2 kernel function
    ! if (dim==2) then
    !     factor=7.0/(4*PI*average_smooth_length**2)
    ! elseif (dim==3) then
    !     factor=21.0/(16*PI*average_smooth_length**3)
    ! endif
      
         
    ! !w and dwdx
    ! if(distance<=2.0*average_smooth_length) then
        
    !     !Kernel function value 
    !     w=factor*(1+2*q)*(1-q/2.0)**4
        
    !     !Derivates of kernel function
    !     do i=1,dim
    !         temp_dwdx(i)=factor*(-5*q*(1-q/2)**3)*position_difference(i)/(distance*average_smooth_length)
    !     end do
    
    ! else
    !     w=0.0
    !     temp_dwdx=0.0
    ! end if
    ! !--------------------------------------------------------------------------------------------------

    ! !--------------------------------------------------------------------------------------------------
    ! !Wendland C4 kernel function
    ! if (dim==2) then
    !     factor=9.0/(4*PI*average_smooth_length**2)
    ! elseif (dim==3) then
    !     factor=495/(256*PI*average_smooth_length**3)
    ! endif
      
         
    ! !w and dwdx
    ! if(distance<=2.0*average_smooth_length) then
        
    !     !Kernel function value 
    !     w=factor*(1+3*q+35.0/12.0*q**2)*(1-q/2.0)**6

    !     !Derivates of kernel function
    !     do i=1,dim
    !         temp_dwdx(i)=factor*(-3*(1-q/2.0)**5*(35/12.0*q**2+3*q+1)+(1-q/2.0)**6*(35/6.0*q+3))*position_difference(i)/(distance*average_smooth_length)
    !     end do
    
    ! else
    !     w=0.0
    !     temp_dwdx=0.0
    ! end if
    ! !--------------------------------------------------------------------------------------------------

    ! !--------------------------------------------------------------------------------------------------
    ! !Wendland C6 kernel function
    ! if (dim==2) then
    !     factor=78/(28*PI*average_smooth_length**2)
    ! elseif (dim==3) then
    !     factor=1365/(512*PI*average_smooth_length**3)
    ! endif
      
         
    ! !w and dwdx
    ! if(distance<=2.0*average_smooth_length) then
        
    !     !Kernel function value 
    !     w=factor*(1+4*q+6.25*q**2+4*q**3)*(1-q/2.0)**8

    !     !Derivates of kernel function
    !     do i=1,dim
    !         temp_dwdx(i)=factor*(-4*(1-q/2.0)**7*(4*q**3+6.25*q**2+4*q+1)+(1-q/2.0)**8*(12*q**2+12.5*q+4))*position_difference(i)/(distance*average_smooth_length)
    !     end do
    
    ! else
    !     w=0.0
    !     temp_dwdx=0.0
    ! end if
    ! !--------------------------------------------------------------------------------------------------
   
end subroutine compute_kernel 
    