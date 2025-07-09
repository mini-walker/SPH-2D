!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: Distribute_Fixed_Boundary_Particle
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

subroutine  Distribute_Fixed_Boundary_Particle()

    use Public_variable_module
    use information_module
    use MPI
    
    implicit none

    !==========================================================================================================
    ! Variables 
    integer::i,j,k,l,O,V                                          !内部循环控制变量
    integer::mirror_particle_name=0                               !正在进行镜像粒子编号
    integer::temp_particle_number=0                               !边界一类虚粒子和内部粒子数总和
    integer::mirror                                               !镜像粒子在总粒子数中的编号
    integer::copy_mirror=0                                        !需要镜像的镜像粒子数
    real(kind=8)::position_x,position_y,position_z                !中间变量
    real(kind=8)::angle_x,angle_y                                 !中间变量
    real(kind=8)::angle                                           !夹角
    real(kind=8)::distance,min_distance                           !点到直线的距离
    integer::x_n,y_n,z_n                                          !粒子所在的链表坐标
    integer::before_fix_ghost_particle_number                     !内部实粒子和固定边界虚粒子数量
    !==========================================================================================================

    ! Body of subroutine Define_fixed_boundary_particle_Position


    !内部实粒子和固定边界虚粒子数量
    before_fix_ghost_particle_number=particle_ture_number+wave_maker_particle_number

    
    !对内部粒子进行检验是否靠近顶部和底部边界，靠近则镜像粒子
    do i=1,fix_ghost_particle_number
        
        mirror=i+before_fix_ghost_particle_number                   !计算正在进行镜像的固定镜像粒子编号
       
        position_x=particle_position(mirror,1)
        position_y=particle_position(mirror,2)
        

        !底部直线段
        if(position_y<=0.0) then

            !左下角中心对称
            if (position_x<0.0) then
               
                Boundary_particle_type(mirror)=-7                                                              !镜像粒子类型

            !右下角中心对称
            elseif(position_x>boundary_size_x) then
               
                Boundary_particle_type(mirror)=-7                                                              !镜像粒子类型

            !中间粒子点(可以镜像到流体内)
            else 
                
                Boundary_particle_type(mirror)=-6                                                      !镜像粒子类型

            end if
                
        !TOP parts
        elseif (position_y>=boundary_size_y) then    

            !左上角中心对称
            if (position_x<0.0) then
               
                Boundary_particle_type(mirror)=-7                                                                !镜像粒子类型

            !右上角中心对称
            elseif(position_x>boundary_size_x) then
               
                Boundary_particle_type(mirror)=-7                                                                !镜像粒子类型

            !中间粒子点(可以镜像到流体内)
            else 
                
                Boundary_particle_type(mirror)=-6                                                                !镜像粒子类型

            end if

        !Vertical Left or right
        else
                  
            Boundary_particle_type(mirror)=-8                                                                                                                                                !镜像粒子类型
                                                                                                         !镜像粒子类型
        end if
             
    end do
    

end subroutine Distribute_Fixed_Boundary_Particle