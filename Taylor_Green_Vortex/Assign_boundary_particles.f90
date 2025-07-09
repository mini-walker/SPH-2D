!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: Assign_boundary_particles
!
!  PURPOSE: Assign boundary particles information
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

subroutine Assign_boundary_particles(i_time_step)

    use function_module
    use information_module
    use Public_variable_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables

    ! Variables from the mother subroutine
    integer,intent(in)::i_time_step                              !Current time step

    !-----------------------------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------------------------------------------
    ! Variables in subroutine
    integer::i,j,k,l,o,v,r                                                           !内部循环控制变量
    integer::fix_ghost_name                                                          !固定镜像粒子的编号
    integer::mesh_name                                                               !网格名称
    real(kind=8)::temp_w                                                             !核函数值
    real(kind=8),dimension(dim)::temp_dwdx                                           !核函数偏导数 
    integer,dimension(NumberInEachGridPrediction)::effect_particle                   !有效粒子存储数组
    integer::effect_particle_number                                                  !有效粒子数量
    integer::ture_effect_number                                                      !真正有效粒子数量
    integer,dimension(NumberInEachGridPrediction)::ture_effect_name                  !真正有效粒子存储数组
    integer::near_mesh_number                                                        !粒子邻居网格数量
    integer::near_mesh_name                                                          !粒子邻居网格编号
    integer::particle_in_mesh                                                        !网格粒子数量
    real(kind=8)::distance                                                           !i，j粒子点的距离
    real(kind=8)::average_smooth_length                                              !i，j粒子点的平均光滑长度
    real(kind=8),dimension(dim)::position_difference                                 !坐标的差(i-j),1-x,2-y,3-z
    real(kind=8),dimension(dim)::Temp_Velocity                                       !Temp_Velocity
    integer::x_n,y_n,z_n                                                             !粒子所在的链表坐标

    !粒子插值变量
    real(kind=8),dimension(NumberInEachGridPrediction)::W_kernel                     !核函数值
    real(kind=8),dimension(NumberInEachGridPrediction)::Fai                          !权函数值

    integer::file_index
    character(len=40)::file_name
    character(len=4)::char_Current_Processor_ID
    integer::ioerror=0                                                               !open return value
    !----------------------------------------------------------------------------------------------------------
    !==========================================================================================================

    ! Body of subroutine Assign_boundary_particles

    !----------------------------------------------------------------------------------------------------------
    !Initialized the variables
    subdomain_number=0
    subdomain_particle_in_chain=0
    subdomain_particle_chain_number=0

    !Project all the particles (total_particle_number_in_subdomain) to the background grid
    do i=1,total_particle_number_in_subdomain
        
        if(subdomain_particle_type(i)==120 .or. subdomain_particle_type(i)==0) cycle         !No error particles

        !Calculate the coordinates of the particles in the list
        !注意，这里计算的块的坐标，要减去子块的起始坐标,然后要右移一列
        x_n=floor((subdomain_particle_position(i,1)-subdomain_chain_origin_x)/chain_dx+1)+1
        y_n=floor((subdomain_particle_position(i,2)-subdomain_chain_origin_y)/chain_dy+1)

        !Calculate the grid number of the particle in the background grid
        k=(y_n-1)*subdomain_chain_x_number+x_n                                 !Row by row

        if (k<=subdomain_chain_max_number) then
            
            subdomain_particle_chain_number(i)=k                               !Save the particle background grid index

            !write(*,*) k,chain_max_number
            subdomain_number(k)=subdomain_number(k)+1                          !Particle number in current grid index plus 1

            !write(3,*) k
            if(subdomain_number(k)<=NumberInEachGridPrediction) then
                subdomain_particle_in_chain(k,subdomain_number(k))=i           !The particle index in kth grid is i
            else
                write(*,*) "The chain's containablity is not enough! (Find pair)"
            end if

        else
            subdomain_particle_type(i)=120
        endif

        ! if (i==1) then
        !     write(*,*) i,x_n,y_n,k,subdomain_chain_x_number,subdomain_chain_y_number
        ! endif

    end do
    !----------------------------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------------------------
    !Find all the vicinity particles of each particle for total particles(actual+buffer zone)
    
    !Initialized the subdomain effect particle and number
    subdomain_effect_particle=0
    subdomain_effect_particle_number=0

    do i=1,total_particle_number_in_subdomain

        !The background grid index of current particle
        mesh_name=subdomain_particle_chain_number(i)                  !mesh index
        near_mesh_number=subdomain_chain_near_mesh_number(mesh_name)  !vicinity grid index of background grid index of current particle
        
        do k=1,near_mesh_number                         
           
           near_mesh_name=subdomain_chain_near_mesh_name(mesh_name,k) !current background grid index
           particle_in_mesh=subdomain_number(near_mesh_name)          !particle number in current background grid index
           
           do v=1,particle_in_mesh                                    !Search all the particles in current background grid
    
               j=subdomain_particle_in_chain(near_mesh_name,v)

               subdomain_effect_particle_number(i)=subdomain_effect_particle_number(i)+1 
               
               !Save the particle index
               if(subdomain_effect_particle_number(i)<=NumberInEachGridPrediction) then
                   subdomain_effect_particle(i,subdomain_effect_particle_number(i))=j
               else
                   write(*,*) "Memory for effect_particle is not engough, it should be larger than 300! (Assign)"
                   !write(*,*) i,particle_position(i,1),particle_position(i,2)
               end if
               
           end do
        end do

    end do
    !----------------------------------------------------------------------------------------------------------

    !***********************************Assign the boundary particles information******************************* 
    !Assign the boundary particles directly
    if (DistributeFixedBoundaryParticleOrNot==1) then
        
        do i=1,actual_particle_number_in_subdomain

            if(subdomain_particle_type(i)/=-3) cycle

            !Initialized the Variables
            subdomain_particle_velocity(i,:)=0.0d0
            subdomain_particle_press(i)=0.0d0
            subdomain_particle_volume(i)=0.0d0
            
            !--------------------------------------------------------------------------------------------------
            W_kernel=0.0d0
            ture_effect_number=0                                     ! i粒子的真正有效粒子数量
            ture_effect_name=0                                       ! i粒子的真正有效粒子名称
            
            !Search all the vicinity particles
            do L=1,subdomain_effect_particle_number(i)

                !effect particle index
                j=subdomain_effect_particle(i,L)

                if(subdomain_particle_type(j)==-3) cycle            ! No other boundary particles
                
                !position difference and distance between i and j particle
                position_difference(:)=subdomain_particle_position(i,:)-subdomain_particle_position(j,:)
                distance=sqrt(DOT_PRODUCT(position_difference,position_difference))
                
                !计算i，j粒子的平均光滑长度(注意个固定镜像粒子的镜像粒子的涉及到计算的不变的变量进行赋值)
                average_smooth_length=1.2*(subdomain_particle_smooth_lengh(i)+subdomain_particle_smooth_lengh(j))/2.0
                !write(*,*) particle_smooth_lengh(i)
                
                !比较粒子点与(k*平均光滑长度)的大小(并且排除自身点distance==0.0)
                if(distance<=(kernel_scale*average_smooth_length) .and. distance>=1.0E-7) then
                    
                    ture_effect_number=ture_effect_number+1              !ture effect particle number
                    ture_effect_name(ture_effect_number)=j               !ture effect particle index
                    
                    !Get the kernel function
                    call compute_kernel(dim,&                            !dimension
                                        distance,&                       !distance(in)                
                                        position_difference,&            !position difference(in)
                                        average_smooth_length,&          !average smooth length(in)
                                        temp_w,&                         !kernel function value(out)
                                        temp_dwdx&                       !derivates of kernel function(out)
                                        )
                
                    !kernel function value
                    W_kernel(ture_effect_number)=temp_w
      
                end if      
                
            end do
            !--------------------------------------------------------------------------------------------------
            
            !--------------------------------------------------------------------------------------------------
            !Interpolate fixed particle quality
            if(ture_effect_number>=3) then
                     
                !Get the weight value of the effect particle
                call MLS_Interpolation_In_Subdomain_Domain(MLS_order_Boundary,subdomain_particle_position(i,:),ture_effect_number,ture_effect_name,W_kernel,Fai)
           
                !Interpolation
                do L=1,ture_effect_number
            
                    !effect particle index
                    j=ture_effect_name(L)

                    if (trim(adjustl(Gravity_Switch))=='on') then
                    
                        subdomain_particle_press(i)=subdomain_particle_press(i)+subdomain_particle_press(j)*Fai(L)-g*subdomain_particle_rho(j)*(subdomain_particle_position(i,dim)-subdomain_particle_position(j,dim))*Fai(L)
                    
                    else
                        subdomain_particle_press(i)=subdomain_particle_press(i)+subdomain_particle_press(j)*Fai(L)
                    
                    endif

                    subdomain_particle_velocity(i,:)=subdomain_particle_velocity(i,:)+subdomain_particle_velocity(j,:)*Fai(L)

                    !subdomain_particle_volume(i)=subdomain_particle_volume(i)+subdomain_particle_mass(j)/subdomain_particle_rho(j)*Fai(L)
                   
                end do
                
            end if
            !--------------------------------------------------------------------------------------------------

            !--------------------------------------------------------------------------------------------------
            !Assgin the velocity fixed particle quality
            if(subdomain_Boundary_particle_type(i)==-6) then

                subdomain_particle_velocity(i,1)=subdomain_particle_velocity(i,1)
                subdomain_particle_velocity(i,2)=-subdomain_particle_velocity(i,2)

            else if(subdomain_Boundary_particle_type(i)==-7) then
        
                subdomain_particle_velocity(i,1)=-subdomain_particle_velocity(i,1)
                subdomain_particle_velocity(i,2)=-subdomain_particle_velocity(i,2)

            else if(subdomain_Boundary_particle_type(i)==-8) then         !子程序输入的是角度制
                
                ! call cline_free_slip(subdomain_particle_velocity(i,1),subdomain_particle_velocity(i,2),theta_ang,Temp_Velocity(1),Temp_Velocity(2))
                ! subdomain_particle_velocity(i,:)=Temp_Velocity(:)

                subdomain_particle_velocity(i,1)=-subdomain_particle_velocity(i,1)
                subdomain_particle_velocity(i,2)=subdomain_particle_velocity(i,2)
            
            else if(subdomain_Boundary_particle_type(i)==-9) then
                
                !圆心在函数内自动写入，若要封装需要传入圆心坐标
                call radius_free_slip(particle_position(i,1),particle_position(i,2),particle_velocity(i,1),particle_velocity(i,2),particle_velocity(fix_ghost_name,1),particle_velocity(fix_ghost_name,2))

            end if

            !--------------------------------------------------------------------------------------------------
            !Assgin the other qualities of fixed particle 

            if (subdomain_particle_press(i)/=0.0) then

                subdomain_particle_energy(i)=e_0
                subdomain_particle_rho(i)=water_rho_0*(1+7*subdomain_particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                !subdomain_particle_mass(i)=subdomain_particle_rho(i)*subdomain_particle_volume(i)
                subdomain_particle_mass(i)=subdomain_particle_rho(i)*intial_volume_0
                subdomain_particle_c(i)=sqrt(square_c_0)*(subdomain_particle_rho(i)/water_rho_0)**3
            
            else

                subdomain_particle_velocity(i,:)=0.0d0
                subdomain_particle_press(i)=0.0d0
                subdomain_particle_energy(i)=e_0
                subdomain_particle_rho(i)=water_rho_0
                subdomain_particle_mass(i)=0.0d0
                subdomain_particle_c(i)=sqrt(square_c_0)

            endif

            !--------------------------------------------------------------------------------------------------

        end do

        !Exchange the fixed boundary particles in buffer domain
        call exchange_buffer_domain_information(i_time_step)

    endif

    !------------------------------------------------------------------------------------------
    !Don't forget calculate the particle volume
    do k=1,total_particle_number_in_subdomain  

       !if (subdomain_particle_rho(k)<1.0e-8) cycle
        
       !particle volume                          
       subdomain_particle_volume(k)=subdomain_particle_mass(k)/subdomain_particle_rho(k)

       ! if (isnan(subdomain_particle_volume(k))) then
       !      write(*,*) subdomain_particle_mass(k),subdomain_particle_rho(k)
       !     subdomain_particle_volume(k)=intial_volume_0
       ! end if

       ! if (subdomain_particle_volume(k)<1.0e-8) then
       !    subdomain_particle_volume(k)=intial_volume_0
       ! endif

       !write(*,*) subdomain_particle_volume(k)

    end do
    !------------------------------------------------------------------------------------------


!     !***********************************************************************************************************
!     ! Call MPI functions for debugging
!     call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID, ierror_MPI )
!     call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )

!     !***********************************************************************************************************


!     !***********************************************************************************************************
!     !输出检验真正有效的粒子点
!     if(i_time_step==1) then

!       !------------------------------------------------------------------------------------------
!       !open the output tecplot data files
!       file_index=700+Current_Processor_ID

!       !transfer the Current_Processor_ID from integer to character
!       write(char_Current_Processor_ID,'(I4)') Current_Processor_ID
!       !write(*,*) char_Current_Processor_ID

!       file_name="./Subdomain/Boundary_in_subdomain_"//trim(adjustl(char_Current_Processor_ID))//".dat"
!       !write(*,*) file_name

!       open(unit=file_index,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    

!       !tecplot header
!       write(file_index,*) "TITLE='DISTRIBUTION'"
!       write(file_index,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
!       !------------------------------------------------------------------------------------------

!       !------------------------------------------------------------------------------------------
!       write(file_index,*) "ZONE I=",total_particle_number_in_subdomain," F=POINT"
!       do i=1,total_particle_number_in_subdomain
!          write(file_index,170) (subdomain_particle_position(i,j),j=1,dim),(subdomain_particle_velocity(i,j),j=1,dim),subdomain_particle_rho(i),subdomain_particle_press(i)/pressure_0
! 170          format(8F20.10) 
!       end do
!       !------------------------------------------------------------------------------------------

!       close(file_index)
    
!     end if
!     !***********************************************************************************************************
        
    
end subroutine Assign_boundary_particles