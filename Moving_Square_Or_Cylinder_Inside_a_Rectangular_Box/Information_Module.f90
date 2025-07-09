!**************************************************************************************************************
!  MODULE     : Information_Module
!
!  PURPOSE    : Control information module for SPH calculation
!
!  Programer  : Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location   : MUN
!
!  Time       : 2017.3.18
!
!  Copyright  : Memorial University
!
!  Version    : 1.0
!
!  Note       : MPI version: mpich-3.2
!               Fortran version: Inter Fortran (ifort) 18.0.0
!**************************************************************************************************************

Module Information_Module

    implicit none

    !*******************************************Parameter Variables********************************************
    ! Constant Parameters 
    integer,parameter::Linux_or_Windows=1                                                                  ! 0---Windows; 1---Linux

    ! Case dimension
    integer,parameter::dim=2                                                                               ! Case dimension 
    integer,parameter::Compatibility_Equation_dim=dim+1                                                    ! Dimension of the linear equation of compatibility interpolation
    
    ! Parameter Variables
    real(kind=8),parameter::Sound_Speed_In_Air=340.0d0                                                     ! Speed of sound in the air
    real(kind=8),parameter::Sound_Speed_In_Water=1497.0d0                                                  ! Speed of sound in the water
    real(kind=8),parameter::PI=3.14159265358979324d0                                                       ! PI
    real(kind=8),parameter::air_rho_0=1.29d0                                                               ! Air density          [ unit: kg/(m^3)]
    real(kind=8),parameter::water_rho_0=1000.0d0                                                           ! Water density        [ unit: kg/(m^3)]
    real(kind=8),parameter::g=9.81d0                                                                       ! Gravity acceleration [ unit: m/(s^2)]                                        
    real(kind=8),protected,dimension(dim)::Gravity_grad=0.0d0                                              ! Gravity gardience   
    character(len=5)::Gravity_Switch='off'                                                                 ! Gravity switch 'on' or 'off'

    ! Dynamic viscosity: water :1.01×10-3 Pa·s; air: μ＝17.9×10-6 Pa·s;
    ! Kinematic viscosity: water :v＝1.01×10-6 m2/s ; air:  v＝14.8×10-6 m2/s ;
    real(kind=8),parameter::Water_kinematic_viscosity=0.01d0                             ! 0.0666666667(Re150); 0.01 (Re100); 0.02(Re50); 
    real(kind=8),parameter::Air_kinematic_viscosity=14.8E-6

    !----------------------------------------------------------------------------------------------------------


    ! Parameter for Background grid
    ! 4*4*4 
    integer,parameter::PredictionNumberInSingleGrid       = 200                                            ! Prediction particle number in each grid 
    integer,parameter::PredictionNumberInSupportDomain    = 300                                            ! Prediction particle number in support domain
    integer,parameter::PredictionNumberInAllVicinityGrids = 800                                            ! Prediction particle number in support domain
    
    !----------------------------------------------------------------------------------------------------------

    ! Particle Type
    integer,parameter::Air_Particle_Label         =  1                                                     ! Air particle label
    integer,parameter::Ghost_Air_Particle_Label   = -1                                                     ! Ghost air particle label
    integer,parameter::Water_Particle_Label       =  2                                                     ! Water particle label
    integer,parameter::Ghost_Water_Particle_Label = -2                                                     ! Ghost water particle label
    integer,parameter::Fix_Ghost_Particle_Label   = -3                                                     ! Fix ghost particle label
    integer,parameter::Period_Particle_Label      = -101                                                   ! Period particle label
    integer,parameter::Inlet_Particle_Label       = -201                                                   ! Inlet particle label
    integer,parameter::Outlet_Particle_Label      = -202                                                   ! Outlet particle label
    integer,parameter::Wave_Maker_Particle_Label  =  130                                                   ! Wavemaker particle label  

    ! Attention: This label only assign in the initial particle, and the body particles are also the boundary particles
    ! It is used for the motion and force calculation
    integer,parameter::Body_Particle_Label        = -4                                                     ! Body particle label

    integer,parameter::Ill_Particle_Label         = 120                                                    ! Ill particle label
    !----------------------------------------------------------------------------------------------------------

    !----------------------------------------------------------------------------------------------------------
    ! File Manage System
    integer,parameter::Log_File_Port                                = 10086                                ! Log file port                                       ! Log File Port

    integer,parameter::General_File_Port                            = 1                                    ! General File Port ( All the output data will be collected in one time step )
    integer,parameter::Animation_Result_File_Port                   = 2                                    ! All the animation results in one file
    integer,parameter::Animation_Grid_Result_File_Port              = 3                                    ! All the animation results of grid in one file

    integer,parameter::Tecplot_General_File_Port                    = 4                                    ! Tecplot general File Port ( All the output data will be collected in one time step )
    integer,parameter::Paraview_General_File_Port                   = 5                                    ! Paraview general File Port ( All the output data will be collected in one time step )
    
    integer,parameter::Standby_General_File_Port                    = 6                                    ! Standby general file port, it will be used when general file port is using ( All the output data will be collected in one time step )
    
    integer,parameter::Initial_File_Port_For_Each_Processor         = 100                                  ! Initial file index for checking result of each processor
    integer::Current_Processor_File_Index                                                                  ! Current_Processor_File_Index=Initial_File_Port_For_Each_Processor+Current_Processor_ID

    integer,parameter::Standby_Initial_File_Port_For_Each_Processor = 200                                  ! Standby initial file index for checking result of each processor
    integer::Standby_Current_Processor_File_Index                                                          ! Standby_Current_Processor_File_Index=Standby_Initial_File_Port_For_Each_Processor+Current_Processor_ID


    integer,parameter::Restart_Result_File_Port                     = 20000                                ! Restart result File Port ( All the output data will be collected in one time step )
    integer,parameter::Free_Surface_File_Port                       = 20001                                ! Free surface result File Port 
    integer,parameter::Body_Force_File_Port                         = 20002                                ! Body force result File Port 
    integer,parameter::Smoke_Line_File_Port                         = 20003                                ! Smoke line result File Port 
    
    integer,parameter::Wave_Front_File_Port                         = 20004                                ! Wave front result File Port 
    integer,parameter::Kinetic_Energy_File_Port                     = 20005                                ! Kinetic energy result File Port 
    integer,parameter::Maximum_velocity_magnitude_File_Port         = 20006                                ! Maximum velocity magnitude File Port

    ! File port for Taylor-Green vortex
    integer,parameter::Norm_Horizontal_velocity_File_Port           = 20007                                ! Normalized horizontal velocity file port
    integer,parameter::Norm_Vertical_velocity_File_Port             = 20008                                ! Normalized vertical velocity file port
    integer,parameter::Norm_Velocity_magnitude_File_Port            = 20009                                ! Normalized velocity magnitude file port
    integer,parameter::Norm_Pressure_File_Port                      = 20010                                ! Normalized pressure file port


    ! Maximum pressure probe and wave gauge is 100
    integer,parameter::Normalized_Pressure_Sampling_File_Initial_index = 1000                              ! Normalized Pressure Sampling File Initial index
    integer,parameter::Average_Pressure_Sampling_File_Initial_index    = 1100                              ! Average Pressure Sampling File Initial index
    integer,parameter::Kalman_Pressure_Sampling_File_Initial_index     = 1200                              ! Kalman Pressure Sampling File Initial index
    integer,parameter::Wave_Probe_File_Initial_index                   = 1300                              ! Wave Probe File Initial index
    
    !----------------------------------------------------------------------------------------------------------

    !**********************************************************************************************************







    !*******************************************Computing domain Variables*************************************

    !-------------------------------------Deviding Particle position variables---------------------------------
    integer,parameter::Deviding_from_input_files_or_not=1                                                 ! Deviding from input files or not (0 is yes, others are deviding by computer)
    
    ! Input file type: 0---tecplot data file; 1---gambit data file
    integer,parameter::Initial_file_type=0                                  

    character(len=100)::Tecplot_Initial_Particle_File_name='pretreat_particle.txt'                        ! Input Tecplot file name  
    character(len=100)::Gambit_Initial_Particle_File_name ='particle.neu'                                 ! Input Gambit file name  
    
    integer,parameter::Gambit_element_node_number=4                                                       ! Gambit element node number (3 or 4)


    !-----------------------------------------------Domain Variables-------------------------------------------
    ! Attention: Change the boundary szie when you change the domain size at the same time
    ! U_Module---0.1d0 (100) // 0.4d0 (400) // 1.0d0 (1000) // 3.2d0 (3200) // 10.0d0 (10000) 
    real(kind=8),parameter::U_Module=1.0d0                                                                 ! Velocity Module (|U|)

    real(kind=8),parameter::Domain_size_x=10.0d0                                                           ! Doamin size in x direction (from 0 to domain_size_x)
    real(kind=8),parameter::Domain_size_y=5.0d0                                                            ! Doamin size in y direction (from 0 to domain_size_y)
    real(kind=8),parameter::Boundary_size_x=10.0d0                                                         ! Boundary size in x direction (from 0 to domain_size_x)
    real(kind=8),parameter::Boundary_size_y=5.0d0                                                          ! Boundary size in y direction (from 0 to domain_size_y)
     
    real(kind=8)::Block_Center_X=1.5d0                                                                     ! Block center coordinate: X
    real(kind=8)::Block_Center_Y=Domain_size_y/2.0                                                         ! Block center coordinate: Y
    real(kind=8),parameter::Block_Length=1.0d0                                                             ! Block Length
 
    real(kind=8),parameter::Smoke_Width=Block_Length/4.0                                                   ! Smoke Width for post-proceeding
 
    real(kind=8),parameter::Backward_Step_X=0.5d0                                                          ! Backward Step coordinate: X
    real(kind=8),parameter::Backward_Step_Y=0.5d0                                                          ! Backward Step coordinate: Y

    real(kind=8),parameter::Perfectly_Matched_Layer_Thickness=0.25d0                                       ! Perfectly Matched Layer Thickness

    real(kind=8),parameter::PostProceeding_Grid_size_x=boundary_size_x                                     ! Boundary size in x direction (from 0 to domain_size_x)
    real(kind=8),parameter::PostProceeding_Grid_size_y=boundary_size_y                                     ! Boundary size in y direction (from 0 to domain_size_y)
         
    ! Incline bottom start position
    real(kind=8),parameter::incline_x=boundary_size_x                                                      ! Incline bottom start position x
    real(kind=8),parameter::incline_y=0.0d0                                                                ! Incline bottom start position y
    real(kind=8),parameter::incline_z=0.0d0                                                                ! Incline bottom start position z

    real(kind=8),parameter::theta_ang=90.0d0                                                               ! Incline angle at the end of the tank (unit: degree)
    real(kind=8),parameter::theta_rad=theta_ang*PI/180.0                                                   ! Incline angle at the end of the tank (unit: radian)                                                                  
    real(kind=8)::Tank_square=0.5*domain_size_y*(domain_size_x+domain_size_x+domain_size_y/tan(theta_rad)) ! Tank area
    
    !----------------------------------------------------------------------------------------------------------


    !-------------------------------------------Interior particle Variables------------------------------------
    integer,parameter::interior_particle_number_x=400                                                      ! Interior particle number in x direction 
    integer,parameter::interior_particle_number_y=200                                                      ! Interior particle number in y direction          
    integer,parameter::interior_particle_number=interior_particle_number_x*interior_particle_number_y      ! Total interior particle number
    real(kind=8),parameter::interior_dx=domain_size_x/interior_particle_number_x                           ! Interior particle interval in x direction (dx)
    real(kind=8),parameter::interior_dy=domain_size_y/interior_particle_number_y                           ! Interior particle interval in y direction (dy)
    real(kind=8)::intial_volume_0                                                                          ! Intial particle volume
    real(kind=8)::water_particle_mass_0                                                                    ! Intial particle mass
    !----------------------------------------------------------------------------------------------------------

    !-----------------------------------------Background grid information--------------------------------------
    real(kind=8),parameter::chain_dx=4*interior_dx                                                         ! Grid length dx
    real(kind=8),parameter::chain_dy=4*interior_dx                                                         ! Grid length dy
    real(kind=8),parameter::chain_dz=4*interior_dx                                                         ! Grid length dz
    integer,parameter::chain_x_offset=5                                                                    ! Offset of the grid in x direction
    integer,parameter::chain_y_offset=5                                                                    ! Offset of the grid in y direction
    integer,parameter::chain_z_offset=5                                                                    ! Offset of the grid in z direction
    real(kind=8),parameter::chain_origin_x=-chain_x_offset*chain_dx                                        ! Origin coordinates x of the grid
    real(kind=8),parameter::chain_origin_y=-chain_y_offset*chain_dy                                        ! Origin coordinates y of the grid
    real(kind=8),parameter::chain_origin_z=-chain_z_offset*chain_dz                                        ! Origin coordinates z of the grid
    integer,parameter::chain_x_number=int(boundary_size_x/chain_dx+2*chain_x_offset)                       ! Grid number in x direction
    integer,parameter::chain_y_number=int(boundary_size_y/chain_dy+2*chain_y_offset)                       ! Grid number in y direction
    integer,parameter::chain_z_number=1                                                                    ! Grid number in z direction (It is 1 for dim=2)
    integer,parameter::chain_max_number=int(chain_x_number*chain_y_number*chain_z_number)                  ! Max number of the grid 
    
    real(kind=8),parameter::Background_Grid_x=chain_x_number*chain_dx                                      ! The length of background grid in x 
    real(kind=8),parameter::Background_Grid_y=chain_y_number*chain_dy                                      ! The length of background grid in y
    real(kind=8),parameter::Background_Grid_z=chain_z_number*chain_dz                                      ! The length of background grid in z
    
    !----------------------------------------------------------------------------------------------------------


    !**********************************************************************************************************




















    !*************************************************MPI variables********************************************

    !--------------------------------------------Desired Processor Number--------------------------------------
    integer,parameter::Main_Processor=0                                                                    ! Main Processor ID (default is 0)


    ! integer,parameter::MPI_Dim=1                                                                           ! MPI block devide dimension
    ! integer,parameter::Desired_Processor_Number=32                                                         ! Parallel Processor Number
    

    ! For 2D devide method
    integer,parameter::MPI_Dim=2                                                                           ! MPI block devide dimension
    integer,parameter::Desired_Processor_Number=16                                                         ! Parallel Processor Number 
    
    integer,parameter::Processor_Number_Y=2                                                                ! Parallel Processor Number in y
    integer,parameter::Processor_Number_X=Desired_Processor_Number/Processor_Number_Y                      ! Parallel Processor Number in x
    integer,parameter::Processor_Number_XY=Processor_Number_X*Processor_Number_Y                           ! Parallel Processor Number in xy plane
    integer,parameter::Processor_Number_Z=Desired_Processor_Number/Processor_Number_XY                     ! Parallel Processor Number in z

    !----------------------------------------------------------------------------------------------------------
    

    !--------------------------------------Total Particle Number For Memory------------------------------------
    ! This one is very important, it should be large than the total particles
    ! Some times when the array is out of the boundary, the Program will not output the errors
    ! When the interior particle number is not too large (<100000), the total processor can not be too large (<=16)

    ! When the interior particles number is too much, this factor cannot be too much also, you should keep it in a suitable range

    real(kind=8),parameter::n_total_Factor=1.2d0
    integer,parameter::n_total=INT(interior_particle_number*n_total_Factor)                                ! Total particle number(not ture just make this as large as possible)       


    ! Variables for Subdomain particle number
    integer,parameter::subdomain_ntotal=INT(n_total/Desired_Processor_Number*1.8)                          ! Subdomain particle number for memory
    
    ! IF you need a larger smooth length(>1.33dx), you can define the buffer column larger than 1
    integer,parameter::Buffer_column_number=1                                                              ! Buffer column number
    integer,parameter::column_ntotal=buffer_column_number*chain_y_number*200                               ! Particle number of each column for memory
    integer,parameter::row_ntotal=buffer_column_number*chain_x_number*200                                  ! Particle number of each row for memory
    integer,parameter::corner_ntotal=buffer_column_number*buffer_column_number*600                         ! Particle number of each row for memory
    integer,parameter::pair_n=60*subdomain_ntotal                                                          ! Pair number
    !----------------------------------------------------------------------------------------------------------










    !--------------------------------------Boundary Type Control Variables-------------------------------------
    integer,parameter::Has_FreeSurface_Or_Not=0                                                            ! Has free surface ---1; No free surface---0

    integer,parameter::DistributeFixedBoundaryParticleOrNot=1                                              ! Distribute fixed boundary particle or not(Yes-1; No-0)

    integer::DistributePeriodBoundaryParticleOrNot=0                                                       ! Distribute period boundary particle or not(Yes-1; No-0)
    
    integer::DistributeInletOutletBoundaryParticleOrNot=0                                                  ! Distribute in/outlet boundary particle or not(Yes-1; No-0)

    integer::DistributeTraditionalGhostBoundaryParticleOrNot=0                                             ! Distribute traditional ghost boundary particle or not(Yes-1; No-0)

    integer,parameter::Fix_wave_maker_or_not=0                                                             ! Wave maker with fixed boundary particle or not(0 is no fixed wave maker particles, others are yes)
    
    integer,parameter::Make_wave_or_not=0                                                                  ! Make wave or not (0 is no, 1 is yes)

    integer,parameter::BodyBoundaryCondition=0                                                             ! Boundary Condition of moving body ( 1 is FreeSlip or 0 is No Slip ) 

    integer,parameter::ExternalBoundaryCondition=0                                                         ! External boundary Condition ( 1 is FreeSlip or 0 is No Slip ) 

    integer,parameter::laplacianTerm=1                                                                     ! 1 is Morris or 0 is Monaghan

    integer,parameter::EquationOfState=1                                                                   ! 1 is Tammann or 0 is Tait
    
    integer,parameter::Calculate_Body_Force_or_Not=1                                                       ! Calculate Body Force or Not
    
    integer,parameter::Body_Motion_Solver_On_or_Off=1                                                      ! Update the body motion

    integer,parameter::Free_motion_or_Custom_motion=1                                                      ! 0---Free motion; 1---User-defined motion
    !----------------------------------------------------------------------------------------------------------
 
    !==========================================================================================================










    !**********************************************************************************************************
    
    !--------------------------------------Fixed wave maker particle variables---------------------------------
    integer,parameter::wave_maker_layer=10                                                                 ! The layer of wave maker fixed boundary particle


    !-----------------------------------Fixed ghost boundary particle variables--------------------------------
    real(kind=8),parameter::Fix_ghost_dx=interior_dx                                                       ! Fixed boundary particle interval in x direction (dx)
    real(kind=8),parameter::Fix_ghost_dy=interior_dy                                                       ! Fixed boundary particle interval in y direction (dy)
    integer,parameter::Fix_ghost_layer=4                                                                   ! The layer of fixed ghost boundary particle
    !----------------------------------------------------------------------------------------------------------




    !-----------------------------------InletOutlet boundary particle variables--------------------------------
    real(kind=8),protected::InletOutlet_Boundary_Zone_Width                                                ! InletOutlet Boundary Zone Width
    real(kind=8),protected::InletOutlet_Non_Reflection_Zone_Width                                          ! InletOutlet Non Reflection Boundary Zone Width
    
    integer,parameter::Calos_Non_Reflection_Boundary=0                                                     ! Non-Reflection boundary by Calos; 1---Yes, 0---No;
    
    ! Inflow velocity (It is initialized in the bottom subroutine)
    real(kind=8),protected::Inflow_velocity_magnitude = 0.0d0                                              ! Inflow velocity magnitude
    real(kind=8),dimension(dim)::Inflow_velocity = 0.0d0                                                   ! Inflow velocity

    ! Inflow direction
    integer,parameter::Distribute_Inlet_Boundary_X = 0
    integer,parameter::Distribute_Inlet_Boundary_Y = 0
    integer,parameter::Distribute_Inlet_Boundary_Z = 0

    real(kind=8),parameter::Inlet_Boundary_X=0.0d0
    real(kind=8),parameter::Inlet_Boundary_Y=0.0d0
    real(kind=8),parameter::Inlet_Boundary_Z=0.0d0
    



    ! Outflow direction
    integer,parameter::Distribute_Outlet_Boundary_X = 0
    integer,parameter::Distribute_Outlet_Boundary_Y = 0
    integer,parameter::Distribute_Outlet_Boundary_Z = 0

    real(kind=8),parameter::Outlet_Boundary_X=domain_size_x
    real(kind=8),parameter::Outlet_Boundary_Y=0.0d0
    real(kind=8),parameter::Outlet_Boundary_Z=0.0d0
    !----------------------------------------------------------------------------------------------------------



    !-------------------------------------Period boundary particle variables-----------------------------------
    real(kind=8),protected::Period_Boundary_Zone_Width                                                     ! Period Boundary Zone Width (kernel_scale+2)*interior_dx

    ! The maximun and minimun value in X,Y,Z direction
    integer,parameter::X_circulation=0                                                                     ! Period boundary on x direction or not (0-No, 1-Yes)
    real(kind=8),parameter::Period_Boundary_X_Min=0.0d0
    real(kind=8),parameter::Period_Boundary_X_Max=boundary_size_x
    real(kind=8),parameter::Period_Boundary_Length_X=Period_Boundary_X_Max-Period_Boundary_X_Min

    integer,parameter::Y_circulation=0                                                                     ! Period boundary on y direction or not (0-No, 1-Yes)
    
    real(kind=8),parameter::Period_Boundary_Y_Min=0.0d0
    real(kind=8),parameter::Period_Boundary_Y_Max=boundary_size_y
    real(kind=8),parameter::Period_Boundary_Length_Y=Period_Boundary_Y_Max-Period_Boundary_Y_Min
    
    integer,parameter::Z_circulation=0                                                                     ! Period boundary on z direction or not (0-No, 1-Yes)
    
    real(kind=8),parameter::Period_Boundary_Z_Min=0.0d0
    real(kind=8),parameter::Period_Boundary_Z_Max=0.0d0
    real(kind=8),parameter::Period_Boundary_Length_Z=Period_Boundary_Z_Max-Period_Boundary_Z_Min
    !----------------------------------------------------------------------------------------------------------
   
   

    !--------------------------------Traditional ghost boundary particle variables-----------------------------
    real(kind=8),protected::Traditional_Ghost_Boundary_Zone_Width                                          ! Traditional ghost boundary zone width (kernel_scale+1)*interior_dx

    ! The maximun and minimun value in X,Y,Z direction
    integer,parameter::Traditional_Ghost_X_Direction=0                                                     ! Traditional ghost boundary on x direction or not (0-No, 1-Yes)
    integer,parameter::Traditional_Ghost_X_Left  = 0                                                       ! Traditional ghost boundary on x left or not (0-No, 1-Yes)
    integer,parameter::Traditional_Ghost_X_Right = 0                                                       ! Traditional ghost boundary on x right or not (0-No, 1-Yes)
    
    real(kind=8),parameter::Traditional_Ghost_Boundary_X_Min=0.0d0
    real(kind=8),parameter::Traditional_Ghost_Boundary_X_Max=boundary_size_x
    real(kind=8),parameter::Traditional_Ghost_Boundary_Length_X=Traditional_Ghost_Boundary_X_Max-Traditional_Ghost_Boundary_X_Min



    integer,parameter::Traditional_Ghost_Y_Direction=1                                                     ! Traditional ghost boundary on y direction or not (0-No, 1-Yes)
    integer,parameter::Traditional_Ghost_Y_Front = 1                                                       ! Traditional ghost boundary on y front or not (0-No, 1-Yes)
    integer,parameter::Traditional_Ghost_Y_Back  = 1                                                       ! Traditional ghost boundary on y back or not (0-No, 1-Yes)
    
    real(kind=8),parameter::Traditional_Ghost_Boundary_Y_Min=0.0d0
    real(kind=8),parameter::Traditional_Ghost_Boundary_Y_Max=boundary_size_y
    real(kind=8),parameter::Traditional_Ghost_Boundary_Length_Y=Traditional_Ghost_Boundary_Y_Max-Traditional_Ghost_Boundary_Y_Min
    


    integer,parameter::Traditional_Ghost_Z_Direction=0                                                     ! Traditional ghost boundary on z direction or not (0-No, 1-Yes)
    integer,parameter::Traditional_Ghost_Z_Top    = 1                                                      ! Traditional ghost boundary on y top or not (0-No, 1-Yes)
    integer,parameter::Traditional_Ghost_Z_Bottom = 1                                                      ! Traditional ghost boundary on y bottom or not (0-No, 1-Yes)
    
    real(kind=8),parameter::Traditional_Ghost_Boundary_Z_Min=0.0d0
    real(kind=8),parameter::Traditional_Ghost_Boundary_Z_Max=0.0d0
    real(kind=8),parameter::Traditional_Ghost_Boundary_Length_Z=Traditional_Ghost_Boundary_Z_Max-Traditional_Ghost_Boundary_Z_Min
    !----------------------------------------------------------------------------------------------------------



    !----------------------------------------Particle number Variables-----------------------------------------
    integer::Particle_ture_number=0                                                                        ! Actual interior particle number
    
    integer::Ghost_wave_maker_particle_number=0                                                            ! Wave maker number (by using ghost method)
    integer::Fix_ghost_wave_maker_particle_number=0                                                        ! Wave maker number (by using ghost method)
    integer::Wave_maker_particle_number=0                                                                  ! Wave maker particle number

    integer::Fix_ghost_particle_number=0                                                                   ! Fixed ghost particle number
    integer::New_Fix_ghost_particle_number=0                                                               ! New fixed ghost particle number
    
    integer::Traditional_ghost_particle_number=0                                                           ! Traditional ghost particle number

    integer::Period_Boundary_Particle_Number=0                                                             ! Period Boundary Particle Number

    integer::Inlet_Particle_Number=0                                                                       ! Inlet particle number
    integer::Outlet_Particle_Number=0                                                                      ! Outlet particle number
    integer::InletOutlet_Boundary_Particle_Number=0                                                        ! Inlet particle number+Outlet particle number
    
    integer::Ture_total_particle_number=0                                                                  ! Actual total particle number
    
    !----------------------------------------------------------------------------------------------------------

    !**********************************************************************************************************            




    !************************************Non-dimensionalization variables**************************************
    ! Non-dimensionalization reference variables
    real(kind=8),parameter::NonDimensionalization_Reference_Length   = Block_Length                        ! Reference Length for NonDimension
    real(kind=8),parameter::NonDimensionalization_Reference_Velocity = U_Module                            ! Reference Length for NonDimension

    ! Reynolds Number
    real(kind=8),parameter::Reynolds_Number = NonDimensionalization_Reference_Velocity*NonDimensionalization_Reference_Length/Water_kinematic_viscosity 
    
    ! Mach number              
    real(kind=8),parameter::Mach_number     = NonDimensionalization_Reference_Velocity/Sound_Speed_In_Water                                      


    ! Reference variables  
    real(kind=8),parameter::Pressure_0 = Water_rho_0*NonDimensionalization_Reference_Velocity**2
    real(kind=8),parameter::Velocity_0 = NonDimensionalization_Reference_Velocity
    real(kind=8),parameter::Time_0     = NonDimensionalization_Reference_Length/NonDimensionalization_Reference_Velocity
    real(kind=8),parameter::Force_0    = 0.5*Water_rho_0*NonDimensionalization_Reference_Velocity**2*NonDimensionalization_Reference_Length

    !**********************************************************************************************************


    !******************************************Status equation function****************************************
    real(kind=8),parameter::Max_velocity = U_Module                                                        ! Maximun velocity
    real(kind=8),parameter::Square_c_0   = 200*U_Module                                                    ! c_0^2
    real(kind=8),parameter::Cap_gamma    = 100                                                             ! Second order status equation (100 order)
    real(kind=8)::e_0                    = Square_c_0/(7.0-1.0)                                            ! Initial energy (c_0^2/(gamma-1)) gamma=7
    real(kind=8),parameter::CoefficienceBackgroundPressure=0.00d0                                          ! Coefficience for Background Pressure
    real(kind=8),parameter::Background_Pressure=CoefficienceBackgroundPressure*Square_c_0*water_rho_0      ! Background Pressure                  !
    real(kind=8),parameter::U_Char=sqrt(square_c_0)*Mach_number                                            ! Characterstical velocity
    !**********************************************************************************************************













    !---------------------------------------Distribute grid for postproceeding---------------------------------
    integer,parameter::Postproceeding_Grid_From_Input_Files=1                                              ! Postproceeding Grid from input files or not (1 is yes, others are deviding by computer)
    
    ! Input post-proceeding file type: 0---'Tecplot'; 1---'Gambit'; 2---'Pointwise'.
    integer,parameter::Postproceeding_Grid_Initial_File_Type=2                                             ! 0---'Tecplot'; 1---'Gambit'; 2---'Pointwise'.
    
    character(len=100)::Postproceeding_File_Name='Flow_Over_a_Block_Final_0.01.facet'                      ! File name of post-proceeding grid                        


    integer,parameter::Postproceeding_Grid_element_node_number=4                                           ! Postproceeding_Grid_element_node_number element node number (3 or 4)


    ! For safety, we delete one gird
    real(kind=8),parameter::Factor_Grid_Point=1.0d0                                                        ! Factor between element node and points dx

    integer,parameter::element_number_x=INT(PostProceeding_Grid_size_x/(Factor_Grid_Point*interior_dx))    ! Element number in x                                                       
    integer,parameter::element_number_y=INT(PostProceeding_Grid_size_y/(Factor_Grid_Point*interior_dx))    ! Element number in y

    real(kind=8),parameter::element_dx=boundary_size_x/element_number_x
    real(kind=8),parameter::element_dy=boundary_size_y/element_number_y

    integer,parameter::element_node_number_x=element_number_x+1                                            ! Element node number in x                                                            
    integer,parameter::element_node_number_y=element_number_y+1                                            ! Element node number in y

    ! These variables are used for distributing grid array space
    integer::Total_element_node_number
    integer::Total_element_number
    
    !----------------------------------------------------------------------------------------------------------

    !********************************************************************************************************** 







    !***************************************Wave information variables*****************************************  
    ! Tank information
    real(kind=8),parameter::Water_depth=domain_size_y

    ! Period wave information
    character(len=50)::InputWaveInformation='WaveHeight_and_WavePeriod'                                   ! Input Wave Information Type
    !character(len=50)::InputWaveInformation='WaveHeight_and_WaveLength'                                  ! Input Wave Information Type
    !character(len=50)::InputWaveInformation='WaveHeight_and_WaveFreuency'                                ! Input Wave Information Type
    
    !Solitary wave information ( 0.15d0 )
    real(kind=8)::Wave_height    = 0.15d0                                                                 ! wave height (unit: m)

    real(kind=8)::Wave_Period    = 1.0d0                                                                  ! wave Period (unit: s)
    real(kind=8)::Wave_Frequency = 5.0d0                                                                  ! wave Period (unit: rad/s)
    real(kind=8)::Wave_Length    = 1.0d0                                                                  ! Wave Length (unit: m)
    real(kind=8)::Wave_Number    = 5.0d0                                                                  ! Wave number (unit: rad/m)
    real(kind=8)::Wave_Celerity                                                                           ! Wave velocity moving forward
            
    real(kind=8)::Wave_Amplitue                                                                           ! Wave amplitude (unit: m)
    real(kind=8)::Wave_initial_phase_angle=PI/2                                                           ! Wave initial phase angle (unit: rad)


    !------------------------------------------Wave maker variables--------------------------------------------
    ! Make wave or not
    ! Wave type
    character(len=30)::wave_maker_type='PistonType'                                                       ! For soltary wave or others
    
    !character(len=30)::wave_maker_type='WaveTheory'                                                      ! For regular wave or others
    !character(len=30)::wave_type='Stokes_2nd'                                                            ! For stokes 2nd wave
    !character(len=30)::wave_type='Stokes_3rd'                                                            ! For stokes 3th wave
    !character(len=30)::wave_type='Stokes_5th'                                                            ! For stokes 5th wave
    !character(len=30)::wave_type='Airy'                                                                  ! For Airy wave theory 
    character(len=30)::wave_type='Solitarywave'                                                           ! For Solitary wave combine with Piston Type wavemaker

    ! Damping zone for wave maker zone and end zone
    character(len=5)::WavemakerDamping_Switch='off'                                                       ! Damping zone for wave maker 'on' or 'off'
    
    character(len=5)::EndDampingZone_Switch='off'                                                         ! Damping zone for tank end 'on' or 'off'
    
    !character(len=30)::EndDampingZoneType='Cosin_Type'                                                   ! Damping coefficient from the cosin function (Zheng's Thesis)
    character(len=30)::EndDampingZoneType='OpenFoam_Type'                                                 ! Damping coefficient from the exp function (OpenFoam Type)
    !character(len=30)::EndDampingZoneType='DualSPH_Type'                                                 ! Damping coefficient from the DualSPH(Paper) 

    ! Variables for damping zone
    real(kind=8)::EndZoneDampingLengh=3.0                                                                 ! End Zone Damping Lengh
    real(kind=8)::WavemakerDampingLengh=1.0                                                               ! Here we just mention the initial value and Wave maker Damping Lengh (one wave length if used)
    real(kind=8)::WavemakerDampingLenghScale=1.0                                                          ! The scale of wave length for the Wavemaker Damping Lengh
    real(kind=8)::EndZoneDampingAmplitude=0.3                                                             ! Damping coefficient amplitude              
    !----------------------------------------------------------------------------------------------------------


    !**********************************************************************************************************
 







    !**********************************************Variables for SPH solver************************************ 

    !==========================================================================================================
    ! SPH or ALE
    integer,parameter::SPH_Or_ALE=0                                                                        ! 0---SPH; 1---ALE

    !---------------------------------------------Kernel function variables------------------------------------
    ! Kernel function model: 1---Quintic spline function; 
    !                        2---Cubic spline function; 
    !                        3---Quadratic smooth function;
    !                        4---Improved Gauss Kernel Fcuntion;
    !                        5---C2 wendland kernel function; 
    !                        6---C4 wendland kernel function; 
    !                        7---C6 wendland kernel function;
    integer,parameter::Kernel_function_model=4                                                             ! Kernel function model
    real(kind=8),parameter::Smooth_length_factor=1.15d0                                                    ! Factor: h/dx (1.05 or 1.33)


    ! kernel_scale, smooth_length and the relative variables will be initialized in "Kernel_Function_Initialization" 
    character(len=100),protected::Kernel_function_name                                                     ! Kernel Function name
    real(kind=8),protected::Kernel_scale=3                                                                 ! Kernel Function support size for smooth length
    real(kind=8),protected::Smooth_length=Smooth_length_factor*interior_dx                                 ! Smooth length
    real(kind=8),protected::FreeSurfaceLevelSet_Fai                                                        ! Free Surface Level Set Standard for grid ( -0.25*interior_dx/smooth_length )
    !----------------------------------------------------------------------------------------------------------



    !----------------------------------------------------------------------------------------------------------
    ! Variables for XSPH
    real(kind=8)::epsilon=0.10d0                                                                           ! Don't use for accuracy calculation
      
    ! Variables for Monhan artifical viscity
    real(kind=8)::alpha=0.03d0                                                                             ! (0.0-0.1) it can not be too large, this term is for viscous
    real(kind=8)::beta=2.0d0                                                                               ! (0.8-1.2) this term is for particle collision, no effect on the viscous
    real(kind=8)::fai_coefficient=0.1d0                                                                    ! it should be 0.1
     
    ! Variables for Delta-SPH
    ! Delta SPH will make the surface move up, however the Reinitialization operation with MLS will move the furface down
    ! So it is better Delta=0.03 with MLS 20 steps for wave making
    ! Delta_alpha is also same with the artifical viscous term, it can not be too larger and it will smooth the results
    ! Its function is simular with the artifical viscous, But the viscous is small. Influence is small than  artifical viscous
    ! So Delta_alpha set as 0.02 is Ok!
    real(kind=8)::delta=0.1d0                                                                              ! Delta SPH for inter accleration (Recommand:0.1)                             
    real(kind=8)::Delta_alpha=0.02d0                                                                       ! Delta SPH for inter heat  (Recommand:0.02) 
    real(kind=8)::chi=0.1d0                                                                                ! Delta SPH for viscous heat 
    !----------------------------------------------------------------------------------------------------------



    !----------------------------------------------------------------------------------------------------------
    ! Upwind MUSCL direction ( Position direction or velocity direction as the upwind direction )
    integer,parameter::Position_Or_Velocity_Direction=1                                                    ! Position (Xj>Xi) ---0; Velocity direction ---1
    !----------------------------------------------------------------------------------------------------------



    !----------------------------------------------------------------------------------------------------------
    ! Variables for density refreshing
    ! Refresh the domain by MLS or not
    ! It may casue uncertainty for the regular linear wave surface
    ! So I Recommand turn off this opinion when simulating the regular wave
    character(len=5)::RefreshDensityByMLS_Switch='on'                                                      ! Refresh denisity by MLS 'on' or 'off'
    
    ! 1---Normalized; 3---Liner MLS; 6---Quadric MLS
    ! Recommand: 1, Normalized with every 20-100 steps, results are good and compute comsuming is low;
    ! Recommand: 2, Also you can try Quadric MLS for its high accuracy, the reuslts are meet well with the crest of the pressure;
    ! Recommand: 3, Linear MLS also good in same cases, but it is the last choice in my opinion.       
    ! Attention: The refresh domain density may cause the free surface going up automaticaly
    !            It should be closed for wave making
    integer,parameter::RefreshDensity_MLS_Order=3                                                          ! MLS order for refresh density(m=1, 3, 6) 
    integer,parameter::RefreshDensity_time_step=40                                                         ! Refresh density time step for refresh the domain density


    ! Variables for particle shift 
    character(len=5)::ShiftParticle_Switch='on'                                                            ! Shift Particle 'on' or 'off'
    real(kind=8)::shift_belta=0.5d0                                                                        ! Recommand 0.5     
    real(kind=8)::shift_R=0.2d0                                                                            ! Recommand 0.2    
    integer,parameter::shift_n=4                                                                           ! Recommand 4
    real(kind=8),parameter::max_shift=0.15d0*interior_dx                                                   ! Max shift distance(Recommand 0.2*smooth_length)
    
    !----------------------------------------------------------------------------------------------------------









    !********************************************Time integral Variables***************************************


    !==========================================================================================================                                                    
    ! Time iteration Control
    ! Recommand: 1, Eular is the best choice for debug for its low time consuming
    !            2, The time step, dt, can be 2 or 3 times of which for Eular when using RUGGER KUTA
    !            3, RUGGER KUTA can make the simulation more stable
    integer::Iteration_Method = 2                                                                         ! Time iteration type: 1 is leapforg; 2 is Eular; 3 is Predictor-Corrector; 4 is Beeman; other is RUGGER KUTA

    !----------------------------------------------------------------------------------------------------------
    ! Initial information source
    ! Simulation start from initial time or restart time : (1) 0 --- Initial time;
    !                                                      (2) 1 --- Restart time. 
    integer,parameter::StartFrom = 0                                                                      ! LatestTime or InitialTime
    real(kind=8)::Restart_Time   = 17.0d0                                                                 ! Restart real time 

    real(kind=8)::dt=1.0E-4                                                                               ! Time step length (dt) 
    real(kind=8),parameter::Fixed_dt                = 1.0E-4                                              ! Fixed time step length (Fixed dt) 
    integer,parameter::Max_time_step                = 80000                                               ! Max time step
    
    real(kind=8),parameter::Standing_time           = 0.0d0                                               ! Time for standing simulation
    real(kind=8),parameter::Wavemaker_ramping_time  = 1.0d0                                               ! Wave maker ramping time
    real(kind=8),parameter::Boundary_ramping_time   = 0.0d0                                               ! Boundary ramping time
    real(kind=8),parameter::Shock_wave_release_time = 10.0d0                                               ! Time for release shocking wave


    ! It seems that changing the boundary type in the running period, it will make unstable, shockwave will generated
    integer::BoundaryTransferOpeartion = 0                                                                ! 0---Need transfer; 1---No transfer

    real(kind=8),parameter::ExternalBoundaryTransferTime        = 10.0d0                                  ! Time for changing the external boundary 
    real(kind=8),parameter::InletOutletBoundaryRunningTime      = ExternalBoundaryTransferTime            ! Inlet outlet boundary running time
    real(kind=8),parameter::PeriodBoundaryRunningTime           = ExternalBoundaryTransferTime            ! Period boundary running time
    real(kind=8),parameter::TraditionalGhostBoundaryRunningTime = ExternalBoundaryTransferTime            ! Traditional Ghost boundary running time

    real(kind=8)::InOutlet_Ramping_coefficient                                                            ! In/Outlet boundary ramping coefficient
    real(kind=8)::InOutlet_Ramping_accleration_coefficient                                                ! In/Outlet boundary ramping accleration coefficient
    real(kind=8)::Wavemaker_Ramping_coefficient                                                           ! Wavemaker boundary ramping coefficient
    real(kind=8)::Wavemaker_Ramping_accleration_coefficient                                               ! Wavemaker boundary ramping accleration coefficient


    !----------------------------------------------------------------------------------------------------------
    ! Result saving control

    ! Output domain size (It is initialized in the bottom subroutine 'Kernel_Function_Initialization' );
    ! Only output the desired particle in this domain.
    real(kind=8),protected::OutputDomainXStart,OutputDomainYStart,OutputDomainZStart 
    real(kind=8),protected::OutputDomainXEnd,  OutputDomainYEnd,  OutputDomainZEnd    
   


    character(len=5)::OutputForTecplot_Switch='on'                                                        ! Output Files For Tecplot ('on' or 'off')
    character(len=5)::OutputForParaview_Switch='off'                                                      ! Output Files For Tecplot ('on' or 'off')
    character(len=5)::OutputForPostProceedingGrid_Switch='off'                                            ! Output Files For Post Proceeding Grid ('on' or 'off')
    character(len=5)::OutputForRestartSimulation_Switch='off'                                             ! Output Files For Restart calculation ('on' or 'off')

    integer,parameter::Save_time_step           = 1000                                                    ! Results save time step
    integer,parameter::Sampling_time_step       = 200                                                     ! Sampling frequency
    integer,parameter::Save_time_step_for_debug = 5*save_time_step                                        ! Results save time step for debug
    integer,parameter::Display_time_step        = Sampling_time_step                                      ! Results save time step
    
    !==========================================================================================================
















    !----------------------------------------------------------------------------------------------------------
    ! Interpolation Control
    ! 1---Normalized; 3---Liner MLS; 6---Quadric MLS
    ! Recommand: 1; 
    integer,parameter::MLS_order_Grid=1                                      ! MLS order control for output grid calculation


    ! MLS Order Control for Boundary Interpolation
    ! 1---Normalized; 3---Liner MLS; 6---Quadric MLS
    ! Recommand: 1; 
    integer,parameter::MLS_order_Boundary=1                                  ! MLS order (m=1, 3, 6) for boundary interpolation 
    integer,parameter::MLS_order_InOutlet_Boundary=1                         ! MLS order (m=1, 3, 6) for In/Outlet boundary interpolation 
    integer,parameter::MLS_order_Non_Reflection=1                            ! MLS order (m=1, 3, 6) for non-reflection boundary interpolation 
    real(kind=8)::Slighted_modified_alpha=0.00001d0                          ! Slighted_modified_alpha
    !----------------------------------------------------------------------------------------------------------












    !----------------------------------------------------------------------------------------------------------
    ! Sampling Control
    ! The Pressure sampling and wave probe position is defined in "subroutine Initialize_SPH_Calculation"
    
    integer,parameter::Sampling_Point_number=11                              ! Sampling Point(Maximun is 100)

    integer,parameter::Wave_Probe_number=0                                   ! Wave Probe number(Maximun is 20) 

    integer,parameter::SamplingLineNumber=0                                  ! Sampling Line Number
    integer,parameter::SamplingLinePointNumber=201                           ! Sampling Point Number

    integer,parameter::MLS_order_Sampling=3                                  ! MLS order (m=1, 2, 3, 6) for sampling interpolation 
    
    integer,parameter::Directory_Number=11                                   ! Directory Number and the directory name is initialized in 'Initialize_SPH_Calculation SUBROUTINE'

    !----------------------------------------------------------------------------------------------------------

    
    !********************************************************************************************************** 


































































    contains

    !==========================================================================================================
    ! Subroutine Kernel_Function_Initialization
    ! Initialize the kernel function ( Choose the kernel function and get the smooth length )
    subroutine Kernel_Function_Initialization()

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L

        !------------------------------------------------------------------------------------------------------

        ! Body of Subroutine Kernel_Function_Initialization
        
        !------------------------------------------------------------------------------------------------------
        ! Initialize the variables for weight function
        ! kernel_scale (Kernel Function support size for smooth length)
        if ( Kernel_function_model==1 ) then

            ! Quintic Spline Function
            kernel_scale=3.0d0 
            Kernel_function_name='Quintic Spline Function'

        elseif ( Kernel_function_model==2 ) then

            ! Cubic Spline Function
            kernel_scale=2.0d0 
            Kernel_function_name='Cubic Spline Function'

        elseif ( Kernel_function_model==3 ) then

            ! Quadratic smooth function
            kernel_scale=2.0d0 
            Kernel_function_name='Quadratic smooth function'

        elseif ( Kernel_function_model==4 ) then

            ! Improved Gauss Kernel Fcuntion
            kernel_scale=3.0d0 
            Kernel_function_name='Improved Gauss Kernel Fcuntion'

        elseif ( Kernel_function_model==5 ) then

            ! C2 wendland Kernel Function; 
            kernel_scale=2.0d0 
            Kernel_function_name='Wendland C2 Kernel Function'

        elseif ( Kernel_function_model==6 ) then

            ! C4 wendland Kernel Function; 
            kernel_scale=2.0d0 
            Kernel_function_name='Wendland C4 Kernel Function'

        elseif ( Kernel_function_model==7 ) then

            ! C6 wendland Kernel Function; 
            kernel_scale=2.0d0 
            Kernel_function_name='Wendland C6 Kernel Function'

        endif
        !------------------------------------------------------------------------------------------------------
        
        
        !------------------------------------------------------------------------------------------------------
        ! Initialize kernel_scale, smooth_length and the relative variables
        Smooth_length=smooth_length_factor*interior_dx                             ! Smooth length
        
        FreeSurfaceLevelSet_Fai=-0.25*interior_dx/smooth_length                    ! Free Surface Level Set Standard for grid

        Period_Boundary_Zone_Width=(kernel_scale+2)*interior_dx                    ! Period Boundary Zone Width                     (kernel_scale+2)*interior_dx
        Traditional_Ghost_Boundary_Zone_Width=(kernel_scale+1)*interior_dx         ! Traditional Ghost Boundary Zone Width          (kernel_scale+1)*interior_dx
        InletOutlet_Boundary_Zone_Width=(kernel_scale+1)*interior_dx               ! InletOutlet Boundary Zone Width                (kernel_scale+1)*interior_dx
        InletOutlet_Non_Reflection_Zone_Width=3*InletOutlet_Boundary_Zone_Width    ! InletOutlet Non Reflection Boundary Zone Width (3*InletOutlet_Boundary_Zone_Width)
        
        if ( trim(adjustl(Gravity_Switch))=='on' ) then
            Gravity_grad(dim)=-g                                                   ! Gravity gardience
        endif
        !------------------------------------------------------------------------------------------------------


        !------------------------------------------------------------------------------------------------------
        ! Initialize the output domain size, only output the desired particle in this domain
        OutputDomainXStart = 0.0d0         - kernel_scale*Smooth_length ;
        OutputDomainXEnd   = Domain_size_x + kernel_scale*Smooth_length ; 

        OutputDomainYStart = 0.0d0         - kernel_scale*Smooth_length ;
        OutputDomainYEnd   = Domain_size_y + kernel_scale*Smooth_length ; 

        OutputDomainZStart = 0.0d0 ;
        OutputDomainZEnd   = 0.0d0 ; 
        !------------------------------------------------------------------------------------------------------


        !------------------------------------------------------------------------------------------------------
        ! Initialize the inflow velocity in three direction
        do i=1,dim

            if ( i==1 ) then

                Inflow_velocity(i) = Inflow_velocity_magnitude

            elseif ( i==2 ) then

                Inflow_velocity(i) = 0.0d0

            elseif ( i==3 ) then

                Inflow_velocity(i) = 0.0d0

            endif

        enddo

        !------------------------------------------------------------------------------------------------------


    end subroutine Kernel_Function_Initialization
    !==========================================================================================================

      
end Module Information_Module