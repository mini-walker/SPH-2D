!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  PROGRAM: Public_variable_module
!
!  PURPOSE: Public variable module for SPH calculation 
!          (contains: all the public variables used in the SPH program)
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

module Public_variable_module

    use information_module
    use MPI

    implicit none
    
      !*******************************Basic particle information arraies***************************************
      !Note: the dimension of the array should be parameter, can not be a general variable
      real(kind=8),save,dimension(n_total,dim)::particle_position                  !粒子位置坐标(粒子编号，x/y/z)
      real(kind=8),dimension(n_total,dim)::Normal_vector                           !particle Normal vector 
      real(kind=8),save,dimension(n_total)::particle_mass                          !粒子点质量
      integer,save,dimension(n_total)::particle_type                               !粒子点类型
      integer,save,dimension(n_total)::particle_initial_type                       !粒子点初始类型
      real(kind=8),save,dimension(n_total)::particle_rho                           !粒子点密度
      real(kind=8),save,dimension(n_total,dim)::particle_velocity                  !粒子点速度(粒子编号，某方向速度)
      real(kind=8),save,dimension(n_total)::particle_smooth_lengh                  !粒子光滑长度
      real(kind=8),save,dimension(n_total)::particle_c                             !粒子声速
      real(kind=8),save,dimension(n_total)::particle_press                         !粒子压力
      real(kind=8),save,dimension(n_total)::particle_energy                        !粒子能量
      real(kind=8),save,dimension(n_total)::particle_eta                           !粒子动力粘性系数
      integer,save,dimension(n_total)::free_surface_type                           !粒子点是否为自由表面粒子(0-不是自由表面粒子，1-是自由表面粒子)
      integer,save,dimension(n_total)::Boundary_particle_type                      !固定边界粒子的速度赋值类型
      integer,save,dimension(n_total)::wave_maker_particle_layer                   !造波机粒子层数
      real(kind=8),save,dimension(n_total,dim)::Initial_particle_Position          !Initial particle Position at t=0
      real(kind=8),save,dimension(n_total,3)::particle_vorticity                   !particle vorticity

      !Analytical results for particles
      real(kind=8),dimension(n_total,dim)::Analytical_particle_velocity            !Analytical particle velocity
      real(kind=8),dimension(n_total)::Analytical_particle_press                   !Analytical particle pressure
      real(kind=8),dimension(n_total)::Analytical_particle_velocity_magnitude      !Analytical particle pressure magnitude
      real(kind=8),dimension(n_total,dim)::particle_velocity_error                 !Particle velocity error 
      real(kind=8),dimension(n_total)::particle_press_error                        !Particle pressure error
      real(kind=8),dimension(n_total)::particle_velocity_magnitude_error           !particle pressure magnitude error
      
      !----------------------------------------------------------------------------------------------------------
      !Variables for wavemaker damping zone      
      integer,save,dimension(n_total)::In_WaveMakerDampingZone_OrNot               !Particle In WaveMaker Damping Zone Or Not (1 is in; 0 is not in)
      real(kind=8),save,dimension(n_total,dim)::WaveMaker_Analytical_position      !Analytical position in wavemaker damping zone
      real(kind=8),save,dimension(n_total,dim)::WaveMaker_Analytical_velocity      !Analytical velocity in wavemaker damping zone
      real(kind=8),save,dimension(n_total)::WaveMaker_Analytical_rho               !Analytical density in wavemaker damping zone
      real(kind=8),save,dimension(n_total)::WaveMaker_Analytical_press             !Analytical pressure in wavemaker damping zone
      real(kind=8),save,dimension(n_total)::position_soft_coefficient              !position soft coefficient
      !----------------------------------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------------------------
      !Variables for post-proceeding grid 
      !Note: the dimension of the array should be parameter, can not be a general variable
      real(kind=8),dimension(total_element_node_number,dim)::element_node_position      !element node position
      
      real(kind=8),dimension(total_element_node_number,dim)::element_node_velocity      !element node velocity
      real(kind=8),dimension(total_element_node_number,3)::element_node_vorticity       !element node vorticity
      real(kind=8),dimension(total_element_node_number)::element_node_press             !element node pressure
      real(kind=8),dimension(total_element_node_number)::element_node_Fai               !element node level set value

      real(kind=8),dimension(total_element_node_number,dim)::temp_element_node_velocity !temp element node velocity
      real(kind=8),dimension(total_element_node_number)::temp_element_node_press        !temp element node pressure
      real(kind=8),dimension(total_element_node_number)::temp_element_node_Fai          !temp element node Level set value
      real(kind=8),dimension(total_element_node_number,3)::temp_element_node_vorticity  !temp element node vorticity

      real(kind=8),dimension(total_element_node_number)::Reduce_Real_Array_For_Grid     !Reduce Real Array For Grid 

      integer,dimension(total_element_node_number)::element_node_Type                   !element node type
      integer,dimension(total_element_node_number)::element_node_processor_index        !element node processor index
      integer,dimension(total_element_node_number)::element_node_chain_number           !element node chain number
      
      integer,dimension(total_element_number,8)::element_node_index                     !element node index
      integer::NodeNumberInOneElement                                                   !Node Number In One Element
      
      !Analytical results for element nodes
      real(kind=8),dimension(total_element_node_number,dim)::Analytical_element_node_velocity            !Analytical particle velocity
      real(kind=8),dimension(total_element_node_number)::Analytical_element_node_press                   !Analytical particle pressure
      real(kind=8),dimension(total_element_node_number)::Analytical_element_node_velocity_magnitude      !Analytical element node magnitude
      real(kind=8),dimension(total_element_node_number,dim)::element_node_velocity_error                 !element node velocity relative error 
      real(kind=8),dimension(total_element_node_number)::element_node_press_error                        !element node pressure relative error
      real(kind=8),dimension(total_element_node_number)::element_node_velocity_magnitude_error           !element node velocity magnitude

      !----------------------------------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------------------------
      !Temp variables for reduction
      real(kind=8),dimension(n_total)::Temp_Reduce_Real_Array
      integer,dimension(n_total)::Temp_Reduce_Int_Array
      integer::Temp_Reduce_Int_Variable
      !----------------------------------------------------------------------------------------------------------
      
      integer,allocatable,dimension(:)::nearest_particle_index_for_FBP                  !Nearest particle index for FBP

      !*******************************MPI temp particle information arraies*************************************
      real(kind=8),save,dimension(subdomain_ntotal,dim)::subdomain_particle_position                  !粒子位置坐标(粒子编号，x/y/z)
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_particle_mass                          !粒子点质量
      integer,save,dimension(subdomain_ntotal)::subdomain_particle_type                               !粒子点类型
      integer,save,dimension(subdomain_ntotal)::subdomain_particle_initial_type                       !粒子点初始类型
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_particle_rho                           !粒子点密度
      real(kind=8),save,dimension(subdomain_ntotal,dim)::subdomain_particle_velocity                  !粒子点速度(粒子编号，某方向速度)
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_particle_smooth_lengh                  !粒子光滑长度
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_particle_c                             !粒子声速
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_particle_press                         !粒子压力
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_particle_energy                        !粒子能量
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_particle_eta                           !粒子动力粘性系数
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_particle_volume                        !粒子体积
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_sum_w                                  !粒子体积
      integer,save,dimension(subdomain_ntotal)::subdomain_free_surface_type                           !粒子点是否为自由表面粒子(0-不是自由表面粒子，1-是自由表面粒子)
      integer,dimension(subdomain_ntotal)::subdomain_wave_maker_particle_layer                        !造波机粒子层数
      integer,save,dimension(subdomain_ntotal)::subdomain_Boundary_particle_type                      !固定边界粒子的速度赋值类型
      real(kind=8),dimension(subdomain_ntotal,dim)::subdomain_Normal_vector                           !particle Normal vector 
      real(kind=8),dimension(subdomain_ntotal,dim)::subdomain_sum_dwdx                                !subdomain sum dwdx
      integer,save,dimension(subdomain_ntotal)::subdomain_near_free_surface_or_not                    !subdomain particle near free surface or not(0-no，1-yes)

      real(kind=8),save,dimension(subdomain_ntotal,3)::subdomain_particle_vorticity                   !particle vorticity

      !----------------------------------------------------------------------------------------------------------
      !Variables for wavemaker damping zone      
      integer,save,dimension(subdomain_ntotal)::subdomain_In_WaveMakerDampingZone_OrNot               !subdomain Particle In WaveMaker Damping Zone Or Not (1 is in; 0 is not in)
      real(kind=8),save,dimension(subdomain_ntotal,dim)::subdomain_WaveMaker_Analytical_position      !subdomain Analytical position in wavemaker damping zone
      real(kind=8),save,dimension(subdomain_ntotal,dim)::subdomain_WaveMaker_Analytical_velocity      !subdomain Analytical velocity in wavemaker damping zone
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_WaveMaker_Analytical_rho               !subdomain Analytical density in wavemaker damping zone
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_WaveMaker_Analytical_press             !subdomain Analytical pressure in wavemaker damping zone
      real(kind=8),save,dimension(subdomain_ntotal)::subdomain_position_soft_coefficient              !subdomain position soft coefficient
      !----------------------------------------------------------------------------------------------------------

      !粒子控制方程求解量
      real(kind=8),save,dimension(subdomain_ntotal,dim)::internal_acceleration               !内力产生的加速度(粒子编号，某方向加速度)
      real(kind=8),save,dimension(subdomain_ntotal,dim)::artifical_internal_acceleration     !人工内力产生的加速度(粒子编号，某方向加速度)
      real(kind=8),save,dimension(subdomain_ntotal,dim)::external_acceleration               !外力产生的加速度(粒子编号，某方向加速度)
      real(kind=8),save,dimension(subdomain_ntotal,dim)::artifical_viscosity_acceleration    !人工粘性力产生的加速度(粒子编号，某方向加速度)
      real(kind=8),save,dimension(subdomain_ntotal)::artifical_viscosity_dedt                !人工粘性产生的能量变化率
      real(kind=8),save,dimension(subdomain_ntotal)::artifical_heat_dedt                     !人工热量修正的能量变化率
      real(kind=8),save,dimension(subdomain_ntotal,dim)::particle_acceleration               !粒子加速度(粒子编号，某方向加速度)
      real(kind=8),save,dimension(subdomain_ntotal,dim)::average_velocity                    !粒子平均速度(粒子编号，某方向加速度)
      real(kind=8),save,dimension(subdomain_ntotal)::drhodt                                  !密度变化率
      real(kind=8),save,dimension(subdomain_ntotal)::dedt                                    !能量变化率
      real(kind=8),save,dimension(subdomain_ntotal,dim)::SPS_acceleration                    !SPS viscous term
      real(kind=8),save,dimension(subdomain_ntotal,dim)::laplacian_viscous_acceleration      !laplacian viscous term


      !人工密度修正
      real(kind=8),dimension(subdomain_ntotal,dim)::grad_rho                           !粒子点密度修正项中的密度的梯度
      real(kind=8),dimension(subdomain_ntotal)::artifical_rho_correct                  !人工密度修正
      real(kind=8),dimension(subdomain_ntotal,dim,dim)::dyadic_matrix                  !并矢
      real(kind=8),dimension(subdomain_ntotal,dim,dim)::dyadic_matrix_inversion        !并矢的逆

      !Gradience
      real(kind=8),dimension(subdomain_ntotal,dim)::grad_velocity_u                    !gradience of horizontial velocity
      real(kind=8),dimension(subdomain_ntotal,dim)::grad_velocity_v                    !gradience of vertical velocity
      real(kind=8),dimension(subdomain_ntotal,dim)::grad_pressure                      !gradience of pressure

      ! CSPH的张量计算变量
      real(kind=8),save,dimension(subdomain_ntotal,dim,dim)::Tensor_dyadic_matrix                             !并矢
      real(kind=8),save,dimension(subdomain_ntotal,dim,dim)::Tensor_dyadic_inverse_matrix                     !并矢的逆矩阵
      real(kind=8),save,dimension(subdomain_ntotal,dim,dim)::density_right_hand_side                          ! 正则化密度的右边项
      real(kind=8),dimension(subdomain_ntotal,dim,dim)::density_change_matrix                                 ! 正则化密度结果
      real(kind=8),save,dimension(subdomain_ntotal,dim)::Pressure_right_hand_side                             ! 正则化压力梯度的右边项

      ! 采用正则化核函数进行插值采用，而不是用平均值 (最多允许100个采用点，每个采用点邻近粒子最多200个)
      real(kind=8),save,dimension(100,dim)::Sampling_Point_Position                                 ! 采样点位置
      real(kind=8),dimension(100)::Sampling_Point_Normalized_Pressure                               ! 正则化采样点压力
      real(kind=8),dimension(100)::Sampling_Point_Average_Pressure                                  ! 平均采样点压力
      real(kind=8),dimension(100)::Sampling_Point_Kalman_filter_Pressure                            ! Kalman_filter正则化采样点压力
      integer,dimension(100)::Sampling_Position_near_particle_number                                ! 采样点邻近粒子个数
      integer,dimension(100,200)::Sampling_Position_near_particle_name                              ! 采样点邻近粒子编号

      !Variables for wave porbe 
      real(kind=8),dimension(20,dim)::Wave_Probe_Position                                           ! Wave Probe Position
      real(kind=8),dimension(20)::Wave_Probe_elevation                                              ! Wave elevation at Probe Position
      real(kind=8),dimension(20)::distanceBetween_WP_and_P                                          ! Distance between Wave Probe and particles

      !*******************************kalman_filter******************************
      real(kind=8),dimension(Sampling_Point_number)::last_step_P,last_step_X

      integer::FluidParticleNumberInWaveMakerDamping                                           !Fluid Particle Number In WaveMaker Damping Zone

      !*******************************自由表面变量*******************************   
      !子块自由表面变量
      integer,dimension(10000)::subdomain_surface_particle_name                                 !自由表面粒子编号
      integer::subdomain_surface_particle_number                                                !自由表面粒子数量(out)
      integer::New_subdomain_surface_particle_number                                            !新的自由表面粒子数量(out)
      real(kind=8),dimension(subdomain_ntotal)::subdomain_fra_v                                 !粒子体积分数

      !总的自由表面变量
      integer,dimension(30000)::surface_particle_name                                           !自由表面粒子编号
      integer::surface_particle_number                                                          !自由表面粒子数量(out)


      !**********************************迭代变量********************************
      real(kind=8),dimension(subdomain_ntotal)::before_rho                     !粒子点之前密度
      real(kind=8),dimension(subdomain_ntotal)::before_energy                  !粒子点之前能量
      real(kind=8),dimension(subdomain_ntotal,dim)::before_velocity            !粒子点之前速度
      real(kind=8),dimension(subdomain_ntotal,dim)::before_position            !粒子点之前的位置

      real(kind=8),dimension(subdomain_ntotal,dim)::temp_particle_velocity     !欧拉迭代预测步速度
      real(kind=8),dimension(subdomain_ntotal)::temp_drhodt                    !欧拉迭代预测步密度变化率
      real(kind=8),dimension(subdomain_ntotal)::temp_dedt                      !欧拉迭代预测步能量变化率
      real(kind=8),dimension(subdomain_ntotal,dim)::temp_particle_acceleration !欧拉迭代预测步加速度

      real(kind=8),dimension(subdomain_ntotal,dim)::k_1_temp_particle_velocity     !Runge–Kutta预测步速度
      real(kind=8),dimension(subdomain_ntotal)::k_1_temp_drhodt                    !Runge–Kutta预测步密度变化率
      real(kind=8),dimension(subdomain_ntotal)::k_1_temp_dedt                      !Runge–Kutta预测步能量变化率
      real(kind=8),dimension(subdomain_ntotal,dim)::k_1_temp_particle_acceleration !Runge–Kutta预测步加速度

      real(kind=8),dimension(subdomain_ntotal,dim)::k_2_temp_particle_velocity     !Runge–Kutta预测步速度
      real(kind=8),dimension(subdomain_ntotal)::k_2_temp_drhodt                    !Runge–Kutta预测步密度变化率
      real(kind=8),dimension(subdomain_ntotal)::k_2_temp_dedt                      !Runge–Kutta预测步能量变化率
      real(kind=8),dimension(subdomain_ntotal,dim)::k_2_temp_particle_acceleration !Runge–Kutta预测步加速度

      real(kind=8),dimension(subdomain_ntotal,dim)::k_3_temp_particle_velocity     !Runge–Kutta预测步速度
      real(kind=8),dimension(subdomain_ntotal)::k_3_temp_drhodt                    !Runge–Kutta预测步密度变化率
      real(kind=8),dimension(subdomain_ntotal)::k_3_temp_dedt                      !Runge–Kutta预测步能量变化率
      real(kind=8),dimension(subdomain_ntotal,dim)::k_3_temp_particle_acceleration !Runge–Kutta预测步加速度

    
      !总块网格变量
      integer,save,dimension(chain_max_number)::chain_near_mesh_number                   !链表网格的相邻网格数量
      integer,save,dimension(chain_max_number,27)::chain_near_mesh_name                  !链表网格的相邻网格编号(包括自身)
      integer,save,dimension(chain_max_number,dim)::mesh_center_positon                  !网格中心点坐标
      integer,dimension(chain_max_number,NumberInEachGridPrediction)::particle_in_chain  !链表中的粒子编号
      integer,dimension(chain_max_number)::number                                        !链表网格中粒子数
      integer,dimension(chain_x_number)::particle_number_in_each_column_mesh             !每列链表网格中的总粒子数
      integer,dimension(n_total)::particle_chain_number                                  !粒子所在网格的编号


      !********************************起始控制变量*****************************
      integer::actual_start_time_step                                               !实际起始时间步

      !-------------------------------------------------------------------------
      !Variables of the MPI Running
      integer::ierror_MPI                                                           ! MPI Function return value for success runing or not
      integer,dimension(MPI_STATUS_SIZE)::Status_MPI                                ! MPI running status
      integer::Current_Processor_ID                                                 ! Current Processor ID
      integer::Total_Processors_Number                                              ! Total Processors Number (0,1,2,...,Total_Processors_Number-1)
      character(len=MPI_MAX_PROCESSOR_NAME)::Current_Processor_Group_Name           ! The Host Name For Different Computer (When running case on different computers)
      integer::Name_Length                                                          ! The Processor ID length
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !Variables for Subdomain chain                                              
      integer::subdomain_chain_x_number
      integer::subdomain_chain_y_number
      integer::subdomain_chain_z_number
      integer::subdomain_chain_max_number

      real(kind=8)::subdomain_chain_origin_x                                         !子块链表坐标原点x
      real(kind=8)::subdomain_chain_origin_y                                         !子块链表坐标原点y
      real(kind=8)::subdomain_chain_origin_z                                         !子块链表坐标原点z

      integer::actual_subdomain_chain_x_number
      integer::total_subdomain_chain_x_number

      integer,allocatable,dimension(:,:)::subdomain_chain_name                         !子块里的网格编号
      integer,allocatable,dimension(:)::subdomain_column_number                        !子块里的列号
      integer,allocatable,dimension(:)::subdomain_each_column_total_particle_number    !子块每列网格中的所有粒子数

      !子块网格变量
      integer,allocatable,dimension(:)::subdomain_number                               !子块链表网格中粒子数
      integer,allocatable,dimension(:,:)::subdomain_particle_in_chain                  !子块的网格中的粒子编号
      integer,allocatable,dimension(:)::subdomain_chain_near_mesh_number               !子块链表网格的相邻网格数量
      integer,allocatable,dimension(:,:)::subdomain_chain_near_mesh_name               !子块链表网格的相邻网格编号(包括自身)
      integer,allocatable,dimension(:,:)::subdomain_mesh_center_positon                !子块网格中心点坐标

      integer,dimension(subdomain_ntotal)::subdomain_particle_chain_number             !子块中粒子的网格编号

      integer,dimension(subdomain_ntotal,NumberInEachGridPrediction)::subdomain_effect_particle     !子块中粒子的有效粒子存储数组
      integer,dimension(subdomain_ntotal)::subdomain_effect_particle_number                         !子块中粒子的有效粒子数量

      !-------------------------------------------------------------------------

      integer::actual_particle_number_in_subdomain                                     !子块的实际粒子数（2 - k-1列）
      integer::total_particle_number_in_subdomain                                      !子块的总的粒子数（1 - k列）
      integer::particle_number_in_left_buffer_zone                                     !子块左侧缓冲区的总粒子数（1列）
      integer::particle_number_in_right_buffer_zone                                    !子块右侧缓冲区的总粒子数（1列）
      integer::first_particle_index_in_left_buffer                                     !子块左侧缓冲区的粒子存储起始编号
      integer::first_particle_index_in_right_buffer                                    !子块右侧缓冲区的粒子存储起始编号

      integer,dimension(subdomain_ntotal)::actual_particle_name_in_subdomain           !子块的实际粒子数的编号（对应全局粒子编号）
      integer,dimension(subdomain_ntotal)::total_particle_name_in_subdomain            !子块的总的粒子数的编号（对应全局粒子编号）
      integer,dimension(column_ntotal)::particle_name_in_left_buffer_zone              !左侧缓冲区的粒子编号（对应全局粒子编号）
      integer,dimension(column_ntotal)::particle_name_in_right_buffer_zone             !右侧缓冲区的粒子编号（对应全局粒子编号）
      integer,dimension(subdomain_ntotal)::Full_subdomain_particle_name                !整个子块的粒子编号（对应全局粒子编号）

      integer,dimension(Desired_Processor_Number)::start_column,end_column

      integer::total_particle_number_for_MPI                                           ! All the particles 
      integer::average_particle_number_per_process                                     ! average particle number in per process  
      
      !------------------------------------------------------------------------------------
      !Variables for message exchange
      integer::left_processor_ID,right_processor_ID
      integer::left_Message_ID,right_Message_ID

      integer::message_first_particle_index_to_left,message_first_particle_index_to_right
      integer::message_first_particle_index_from_left,message_first_particle_index_from_right

      integer::message_length_send_to_left,message_length_send_to_right
      integer::message_length_from_left,message_length_from_right
      !------------------------------------------------------------------------------------

      !粒子对信息
      integer::pair_number                                                               !粒子对数
      integer,dimension(pair_n)::pair_i                                                  !第k个粒子对中i的粒子编号  
      integer,dimension(pair_n)::pair_j                                                  !第k个粒子对中j的粒子编号 
      real(kind=8),dimension(pair_n)::w                                                  !第k个粒子对中核函数值
      real(kind=8),dimension(pair_n,dim)::dwdx                                           !第k个粒子对中核函数的偏导数值
      real(kind=8),dimension(pair_n,dim)::Normalized_dwdx_i                              !第k个粒子对中核函数的正则化偏导数值
      real(kind=8),dimension(pair_n,dim)::Normalized_dwdx_j                              !第k个粒子对中核函数的正则化偏导数值
      real(kind=8),dimension(pair_n,dim)::position_difference_ij                         !第k个粒子对中位置差(position_difference_i-position_difference_j)
      real(kind=8),dimension(pair_n,dim)::velocity_difference_ij                         !第k个粒子对中速度差(velocity_difference_i-velocity_difference_j)
      real(kind=8),dimension(pair_n)::rho_difference_ij                                  !第k个粒子对中密度差(rho_i-rho_j)
      real(kind=8),dimension(pair_n)::average_c_ij                                       !第k个粒子对声速平均值(particle_c(i)+particle_c(j))/2.0
      real(kind=8),dimension(pair_n)::average_rho_ij                                     !第k个粒子对密度平均值(particle_rho(i)+particle_rho(j))/2.0
      real(kind=8),dimension(pair_n)::average_smooth_length_ij                           !第k个粒子对光滑长度平均值(particle_smooth_lengh(i)+particle_smooth_lengh(j))/2.0

      !Shift particle position Variables
      real(kind=8),dimension(subdomain_ntotal,dim)::subdomain_shift_particle_position    !subdomain shift particle position
      real(kind=8),dimension(subdomain_ntotal)::sum_mass                                 !Sum mass of the near particles
      real(kind=8),dimension(subdomain_ntotal)::sum_distance                             !Sum distance of the near particles
      real(kind=8),dimension(subdomain_ntotal)::near_particle_max_velocity               !Max velocity of the near particles
      integer,dimension(subdomain_ntotal)::shift_effect_particle_number                  !Effect particle number of the near particles

      !wave maker information
      real(kind=8)::wave_maker_position                                       !wave maker position
      real(kind=8)::wave_maker_velocity                                       !wave maker velocity
      real(kind=8)::wave_maker_acceleration                                   !wave maker acceleration
      real(kind=8)::wave_maker_height                                         !wave maker height(It is setted in Subroutine-Distribute_Initial_Particle)
      
      real(kind=8)::left_wave_height                                 !造波机端的波高
      real(kind=8)::right_wave_height                                !右端端的波高

      real(kind=8),save,dimension(subdomain_ntotal)::Reinitialize_subdomain_particle_rho
      real(kind=8),save,dimension(subdomain_ntotal,dim)::Reinitialize_subdomain_particle_velocity

      !--------------------------------------------------------------------------------
      !SPS viscous variables
      !tao
      real(kind=8),dimension(subdomain_ntotal)::tao_xx                !viscous shear term
      real(kind=8),dimension(subdomain_ntotal)::tao_yy
      real(kind=8),dimension(subdomain_ntotal)::tao_zz
      real(kind=8),dimension(subdomain_ntotal)::tao_xy
      real(kind=8),dimension(subdomain_ntotal)::tao_xz
      real(kind=8),dimension(subdomain_ntotal)::tao_yz 

      !SPS strain tensor
      real(kind=8),dimension(subdomain_ntotal)::SPS_Strain_Tensor_xx  !SPS Strain Tensor
      real(kind=8),dimension(subdomain_ntotal)::SPS_Strain_Tensor_yy
      real(kind=8),dimension(subdomain_ntotal)::SPS_Strain_Tensor_zz
      real(kind=8),dimension(subdomain_ntotal)::SPS_Strain_Tensor_xy
      real(kind=8),dimension(subdomain_ntotal)::SPS_Strain_Tensor_xz
      real(kind=8),dimension(subdomain_ntotal)::SPS_Strain_Tensor_yz
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      !Time Variables
      real(kind=8)::Main_Processor_start_time,Main_Processor_end_time            !Variables of time in calculation
      real(kind=8)::Main_Processor_cost_time                                     !Variables of time in calculation
      real(kind=8)::every_step_time                                              !Variables of time in calculation
      CHARACTER(len=10)::DATE,TIME,PRETTY_TIME 
      
      !---------------------------------------------------------------------------------
      !Wave model variables
      real(kind=8)::Wave_celerity_0                                              !linear wave celerity (Wave_celerity_0^2=g/k tanh(kd))
      real(kind=8)::Wave_Length_0                                                !Deep water wave length (Wave_Length_0=g*T^2/(2*PI))
      
      !Coefficients for 5th stokes wave
      real(kind=8),dimension(5,5)::Fifth_stokes_A  
      real(kind=8),dimension(5,5)::Fifth_stokes_B 
      real(kind=8),dimension(5)::Fifth_stokes_C 
      real(kind=8),dimension(5)::Fifth_stokes_D 
      real(kind=8),dimension(5)::Fifth_stokes_E 
      real(kind=8)::Stokes_S,Stokes_C,Stokes_C_0 
      real(kind=8)::Stokes_Lamda,Temp_Stokes_Lamda,Stokes_Epsilon 
      real(kind=8)::temp_wave_number,Last_wave_number     
      real(kind=8)::Stokes_DFai_Dt     
      real(kind=8)::Stokes_pressure_C_0,Stokes_pressure_S,Bernoulli_constant_R
      real(kind=8),dimension(5)::Fifth_stokes_pressure_E    

      !Coefficients for 3rd stokes wave
      real(kind=8),dimension(3)::Third_stokes_D 
      real(kind=8),dimension(3)::Third_stokes_E 
      real(kind=8)::Third_stokes_Alpha                                   !coth(kd)  
      real(kind=8)::water_depth_bar 

      !Coefficients for 3rd solitary wave
      real(kind=8)::Solitary_alpha                                       !alpha
      real(kind=8)::Solitary_Epsilon                                     !h/d
      real(kind=8)::Solitary_s                                           !sech(alphax)  
      real(kind=8)::Solitary_t                                           !tanh(alphax) 
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      !Variables for wave maker damping
      real(kind=8),dimension(dim)::analytical_position                                     !analytical position
      real(kind=8),dimension(dim)::analytical_velocity                                     !analytical velocity
      real(kind=8)::analytical_press                                                       !analytical press              
      real(kind=8)::analytical_rho                                                         !analytical density
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      !Variables for sampling file
      integer,dimension(Sampling_Point_number)::Normalized_Pressure_Sampling_File_index    !Normalized Pressure Sampling File Initial index+sampling order
      integer,dimension(Sampling_Point_number)::Average_Pressure_Sampling_File_index       !Average Pressure Sampling File Initial index+sampling order
      integer,dimension(Sampling_Point_number)::Kalman_Pressure_Sampling_File_index        !Kalman Pressure Sampling File Initial index+sampling order
      integer,dimension(Wave_Probe_number)::Wave_Probe_File_index                          !Wave_Probe_File_Initial_index+Wave Probe Order
      
      real(kind=8),dimension(SamplingLineNumber,dim)::SamplingLineInitialPoint             !Initial Point on the Sampling line
      real(kind=8),dimension(SamplingLineNumber,dim)::SamplingLineEndPoint                 !End Point on the Sampling line
      integer,parameter::SamplingLinePointNumber=41                                        !Sampling Point Number
      
      real(kind=8),dimension(SamplingLineNumber,SamplingLinePointNumber,dim)::SamplingLinePointPosition    !Sampling Line Point Position
      real(kind=8),dimension(SamplingLineNumber,dim)::SamplingLinePointInterval                            !Sampling Line Point Inyerval
      
      real(kind=8),dimension(SamplingLineNumber,SamplingLinePointNumber,dim)::SamplingLinePointVelocity              !Numerical Velocity of Sampling Points on Line
      real(kind=8),dimension(SamplingLineNumber,SamplingLinePointNumber,dim)::SamplingLinePointAnalyticalVelocity    !Analytical Velocity of Sampling Points on Line
      
      real(kind=8),dimension(SamplingLineNumber,SamplingLinePointNumber)::SamplingLinePointMLSPressure        !MLS Pressure of Sampling Points on Line
      real(kind=8),dimension(SamplingLineNumber,SamplingLinePointNumber)::SamplingLinePointAveragePressure    !Average Pressure of Sampling Points on Line
      real(kind=8),dimension(SamplingLineNumber,SamplingLinePointNumber)::SamplingLinePointAnalyticalPressure !Analytical Pressure of Sampling Points on Line
      
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      !CFL Variables
      real(kind=8)::CFL_Morris_dt,CFL_Monaghan_dt,CFL_dt                                   !CFL dt
      real(kind=8),dimension(2)::Min_CFL_dt
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      real(kind=8)::Real_time,Nondimensionalized_real_time
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      !Energy variables
      real(kind=8)::Initial_total_kinetic_energy                                           !Initial total kinetic energy of the system
      real(kind=8)::Total_kinetic_energy                                                   !Total kinetic energy of the system
      real(kind=8)::Norm_1,Norm_2,Norm_Infinity                                            !norm results 
      !---------------------------------------------------------------------------------

end module Public_variable_module