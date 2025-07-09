!*****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE: information_module
!
!  PURPOSE: Control information module for SPH calculation
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
!*****************************************************************************

module information_module

    implicit none
    
    !***************************Parameter Variables***************************
    ! Parameter Variables
    real(kind=8),parameter::PI=3.14159265358979324d0                               !PI
    real(kind=8),parameter::air_rho_0=1.29d0                                       !air density          [ unit: kg/(m^3)]
    real(kind=8),parameter::water_rho_0=1000.0d0                                   !water density        [ unit: kg/(m^3)]
    real(kind=8),parameter::g=9.81d0                                               !gravity acceleration [ unit: m/(s^2)]                                        
    character(len=30)::Gravity_Switch='off'                                        !Gravity Switch 'on' or 'off'

    !dynamic viscosity: water :1.01×10-3 Pa·s; air: μ＝17.9×10-6 Pa·s;
    !kinematic viscosity: water :v＝1.01×10-6 m2/s ; air:  v＝14.8×10-6 m2/s ;
    real(kind=8),parameter::water_kinematic_viscosity=1.0E-2
    real(kind=8),parameter::air_kinematic_viscosity=14.8E-6

    integer,parameter::dim=2                                                       !case dimensions 

    !-------------------------------------------------------------------------
    !Parameter for Background grid
    ! 4*4*4*3 
    integer,parameter::NumberInEachGridPrediction=300                              !Prediction Particle Number In Each Grid 

    !*************************************************************************

    !***********************Computing domain Variables************************

    !-----------------------------Domain Variables--------------------------
    !Attentation: Change the boundary szie when you change the domain size at the same time
    real(kind=8),parameter::U_Module=1.0d0                                                                !Velocity Module (|U|)
      
    real(kind=8),parameter::domain_size_x=1.0d0                                                            !Doamin size in x direction (from 0 to domain_size_x)
    real(kind=8),parameter::domain_size_y=1.0d0                                                            !Doamin size in y direction (from 0 to domain_size_y)
    real(kind=8),parameter::boundary_size_x=1.0d0                                                          !Boundary size in x direction (from 0 to domain_size_x)
    real(kind=8),parameter::boundary_size_y=1.0d0                                                          !Boundary size in y direction (from 0 to domain_size_y)
    
    real(kind=8),parameter::Reynolds_Number=U_Module*domain_size_x/water_kinematic_viscosity               !Reynolds Number

    real(kind=8),parameter::PostProceeding_Grid_size_x=boundary_size_x                                     !Boundary size in x direction (from 0 to domain_size_x)
    real(kind=8),parameter::PostProceeding_Grid_size_y=boundary_size_y                                     !Boundary size in y direction (from 0 to domain_size_y)
         
    !incline bottom start position
    real(kind=8),parameter::incline_x=boundary_size_x                                                      !incline bottom start position x
    real(kind=8),parameter::incline_y=0.0d0                                                                !incline bottom start position y
    real(kind=8),parameter::incline_z=0.0d0                                                                !incline bottom start position z

    real(kind=8),parameter::theta_ang=90.0d0                                                               !incline angle at the end of the tank (unit: degree)
    real(kind=8),parameter::theta_rad=theta_ang*PI/180.0                                                   !incline angle at the end of the tank (unit: radian)                                                                  
    real(kind=8)::tank_square=0.5*domain_size_y*(domain_size_x+domain_size_x+domain_size_y/tan(theta_rad)) !tank area
    
    !-----------------------------------------------------------------------
    
    !------------------------Interior particle Variables--------------------
    integer,parameter::interior_particle_number_x=100                                                      !interior particle number in x direction 
    integer,parameter::interior_particle_number_y=100                                                      !interior particle number in y direction          
    integer,parameter::interior_particle_number=interior_particle_number_x*interior_particle_number_y      !total interior particle number
    real(kind=8),parameter::interior_dx=domain_size_x/interior_particle_number_x                           !interior particle interval in x direction (dx)
    real(kind=8),parameter::interior_dy=domain_size_y/interior_particle_number_y                           !interior particle interval in y direction (dy)
    real(kind=8)::intial_volume_0                                                                          !intial particle volume
    real(kind=8)::water_particle_mass_0                                                                    !intial particle mass

    !-----------------Fixed ghost boundary particle Variables---------------
    integer,parameter::DistributeFixedBoundaryParticleOrNot=0                                              !Distribute Fixed Boundary Particle Or Not(Yes-1; No-0)

    real(kind=8),parameter::fix_ghost_dx=interior_dx                                                       !Fixed boundary particle interval in x direction (dx)
    real(kind=8),parameter::fix_ghost_dy=interior_dy                                                       !Fixed boundary particle interval in y direction (dy)
    integer,parameter::fix_ghost_layer=4                                                                   !The layer of fixed ghost boundary particle
    
    !--------------------Period boundary particle Variables-----------------
    integer,parameter::DistributePeriodBoundaryParticleOrNot=1                                             !Distribute Period Boundary Particle Or Not(Yes-1; No-0)
    
    !The maximun and minimun value in X,Y,Z direction
    real(kind=8),parameter::Period_Boundary_X_Min=0.0d0
    real(kind=8),parameter::Period_Boundary_X_Max=boundary_size_x
    real(kind=8),parameter::Period_Boundary_Length_X=Period_Boundary_X_Max-Period_Boundary_X_Min

    real(kind=8),parameter::Period_Boundary_Y_Min=0.0d0
    real(kind=8),parameter::Period_Boundary_Y_Max=boundary_size_y
    real(kind=8),parameter::Period_Boundary_Length_Y=Period_Boundary_Y_Max-Period_Boundary_Y_Min

    real(kind=8),parameter::Period_Boundary_Z_Min=0.0d0
    real(kind=8),parameter::Period_Boundary_Z_Max=0.0d0
    real(kind=8),parameter::Period_Boundary_Length_Z=Period_Boundary_Z_Max-Period_Boundary_Z_Min

    !--------------------------Particle number Variables--------------------
    integer::particle_ture_number=0                                                                        !actual interior particle number
    integer::ghost_wave_maker_particle_number=0                                                            !wave maker number (by using ghost method)
    integer::fix_ghost_wave_maker_particle_number=0                                                        !wave maker number (by using ghost method)
    integer::wave_maker_particle_number=0                                                                  !wave maker particle number

    integer::fix_ghost_particle_number=0                                                                   !Fix ghost particle number
    
    integer::Period_Boundary_Particle_Number=0                                                             !Period Boundary Particle Number
    !**************************************************************************                 
  
    !==========================================================================
    !This one is very important, it should be large than the total particles
    !Some times when the array is out of the boundary, the Program will not output the errors
    !When the interior particle number is not too large, this one should be large enough
    !I cost 2 Days on this bug! Fuck
    real(kind=8),parameter::n_total_Factor=2.5d0
    integer,parameter::n_total=int(interior_particle_number*n_total_Factor)                              !Total particle number(not ture just make this as large as possible)       
    
    !==========================================================================
    integer::ture_total_particle_number=0                                                                  !actual total particle number
    !**************************************************************************               


    !*************************Kernel function variables************************
    integer,parameter::kernel_scale=3                                                   !kernel function length scale 
    real(kind=8),parameter::smooth_length_factor=1.20d0                                 !factor: h/dx (1.05 or 1.33)
    real(kind=8),parameter::smooth_length=smooth_length_factor*interior_dx              !smooth length
    real(kind=8),parameter::FreeSurfaceLevelSet_Fai=-0.25*interior_dx/smooth_length     !Free Surface Level Set Standard for grid

    !**************************************************************************  


    !*************************Status equation function*************************
    real(kind=8),parameter::max_velocity=0.4d0                                           !maximun velocity
    real(kind=8),parameter::square_c_0=100*U_Module                                      !c_0^2
    !real(kind=8),parameter::square_c_0=100                                              !c_0^2
    real(kind=8),parameter::cap_gamma=100                                                !two order status equation (100 order)
    real(kind=8)::e_0=square_c_0/(7.0-1.0)                                               !initial energy (c_0^2/(gamma-1)) gamma=7

    !*************************Non-dimensional variables************************
    !Reference variables  
    real(kind=8),parameter::pressure_0=water_rho_0*U_Module**2
    real(kind=8),parameter::velocity_0=U_Module
    real(kind=8),parameter::time_0=domain_size_x/U_Module
    !************************************************************************** 

    !Background Grid information
    real(kind=8),parameter::chain_dx=(kernel_scale+1)*smooth_length                          !grid length dx
    real(kind=8),parameter::chain_dy=(kernel_scale+1)*smooth_length                          !grid length dy
    real(kind=8),parameter::chain_dz=(kernel_scale+1)*smooth_length                          !grid length dz
    integer,parameter::chain_x_offset=5                                                      !Offset of the grid in x direction
    integer,parameter::chain_y_offset=5                                                      !Offset of the grid in y direction
    integer,parameter::chain_z_offset=5                                                      !Offset of the grid in z direction
    real(kind=8),parameter::chain_origin_x=-chain_x_offset*chain_dx                          !origin coordinates x of the grid
    real(kind=8),parameter::chain_origin_y=-chain_y_offset*chain_dy                          !origin coordinates y of the grid
    real(kind=8),parameter::chain_origin_z=-chain_z_offset*chain_dz                          !origin coordinates z of the grid
    integer,parameter::chain_x_number=int(boundary_size_x/chain_dx+2*chain_x_offset)         !grid number in x direction
    integer,parameter::chain_y_number=int(boundary_size_y/chain_dy+2*chain_y_offset)         !grid number in y direction
    integer,parameter::chain_z_number=1                                                      !grid number in z direction (It is 1 for dim=2)
    integer,parameter::chain_max_number=int(chain_x_number*chain_y_number*chain_z_number)    !max number of the grid 
    
    real(kind=8),parameter::Background_Grid_x=chain_x_number*chain_dx                        !The length of background grid in x 
    real(kind=8),parameter::Background_Grid_y=chain_y_number*chain_dy                        !The length of background grid in y
    real(kind=8),parameter::Background_Grid_z=chain_z_number*chain_dz                        !The length of background grid in z
    
    !---------------------------------------------------------------------------
    !Post-proceeding grid (0-boundary size)
    !For safety, we delete one gird
    real(kind=8),parameter::Factor_Grid_Point=1.0d0                                          !Factor between element node and points dx

    integer,parameter::element_number_x=int(PostProceeding_Grid_size_x/(Factor_Grid_Point*interior_dx))        !Element number in x                                                       
    integer,parameter::element_number_y=int(PostProceeding_Grid_size_y/(Factor_Grid_Point*interior_dx))        !Element number in y

    real(kind=8),parameter::element_dx=boundary_size_x/element_number_x
    real(kind=8),parameter::element_dy=boundary_size_y/element_number_y

    integer,parameter::element_node_number_x=element_number_x+1                              !Element node number in x                                                            
    integer,parameter::element_node_number_y=element_number_y+1                              !Element node number in y

    integer,parameter::total_element_node_number=element_node_number_x*element_node_number_y
    integer,parameter::total_element_number=element_number_x*element_number_y

    !*************************Wave information variables************************  
    
    !Tank information
    real(kind=8),parameter::water_depth=domain_size_y

    !Period wave information
    character(len=50)::InputWaveInformation='WaveHeight_and_WavePeriod'                       ! Input Wave Information Type
    !character(len=50)::InputWaveInformation='WaveHeight_and_WaveLength'                      ! Input Wave Information Type
    !character(len=50)::InputWaveInformation='WaveHeight_and_WaveFreuency'                    ! Input Wave Information Type
    !real(kind=8)::wave_height=0.01d0                                                          ! wave height (unit: m)
    real(kind=8)::wave_Period=1.0d0                                                           ! wave Period (unit: s)
    real(kind=8)::Wave_Frequency=5.0d0                                                        ! wave Period (unit: rad/s)
    real(kind=8)::Wave_Length=1.0d0                                                           ! Wave Length (unit: m)
    real(kind=8)::Wave_number=5.0d0                                                           ! Wave number (unit: rad/m)
    real(kind=8)::Wave_celerity                                                               ! Wave velocity

    real(kind=8)::Wave_amplitue                                                               ! wave amplitude (unit: m)
    real(kind=8)::Wave_initial_phase_angle=PI/2                                               ! Wave initial phase angle (unit: rad)

    ! !Solitary wave information
    real(kind=8)::wave_height=0.15d0

    !--------------------Fixed wave maker particle Variables------------------
    integer,parameter::fix_wave_maker_or_not=0                          !wave maker with fixed boundary particles(0 is no fixed wave maker particles, others are yes)
    integer,parameter::wave_maker_layer=10                              !the layer of wave maker fixed boundary particle

    !************************************************************************* 
    !Variables for SPH Model
    !Variables for XSPH
    real(kind=8)::epsilon=0.0d0                                           !Don't use for accuray calculation
      
    !Variables for Monhan artifical viscity
    real(kind=8)::alpha=0.00d0                                            !(0.0-0.1) it can not be too large, this term is for viscous
    real(kind=8)::beta=2.0d0                                              !(0.8-1.2) this term is for particle collision, no effect on the viscous
    real(kind=8)::fai_coefficient=0.1d0                                   !it should be 0.1
     
    !Variables for Delta-SPH
    !Delta SPH will make the surface move up, however the Reinitialization operation with MLS will move the furface down
    !So it is better Delta=0.03 with MLS 20 steps for wave making
    real(kind=8)::delta=0.1d0                                             !Delta SPH for inter accleration (Recommand:0.1)                             
    real(kind=8)::Delta_alpha=0.02d0                                      !Delta SPH for inter heat  (Recommand:0.1) 
    real(kind=8)::chi=0.1d0                                               !Delta SPH for viscous heat 
    !-----------------------------------------------------------------------

    !****************************Time integral Variables********************

    !-------------------Deviding Particle position variables----------------
    integer,parameter::deviding_from_input_files_or_not=1                    !Deviding from input files or not (0 is yes, others are deviding by computer)
    !character(len=10),parameter::initial_file_type='Tecplot'                !input file type --- Tecplot 
    !character(len=10),parameter::initial_file_type='Gambit'                 !input file type --- Gambit
    integer,parameter::gambit_element_node_number=3                          !Gambit element node number (3 or 4)
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !Initial information source
    !character(len=15)::startFrom='LatestTime'                               ! 'LatestTime' or 'InitialTime' 
    character(len=15)::startFrom='InitialTime'                               ! 'LatestTime' or 'InitialTime' 
    character(len=1)::Change_last_start_time_or_not='y'                      ! 'y'---yes; 'n'---no 
    character(len=15)::Initial_file_type ='Tecplot'                          ! 'Tecplot' or 'Gambit' 
    !character(len=15)::Initial_file_type ='Gambit'                          ! 'Tecplot' or 'Gambit' 
    !-----------------------------------------------------------------------

    !*******************************MPI variables*************************** 
    integer,parameter::Main_Processor=0                                      ! Main Processor ID (default is 0)

    !=======================================================================
    integer,parameter::Desired_Processor_Number=5                            ! Parallel Processor Number
    !=======================================================================

    integer,parameter::buffer_column_number=1                                ! Buffer column number
    !IF you need a larger smooth length(>1.33dx), you can define the buffer column larger than 1

    !------------------------------------------------------------------------
    !Variables for Subdomain particle number
    integer,parameter::subdomain_ntotal=int(n_total/Desired_Processor_Number*1.4)  !Subdomain particle number for memory
    integer,parameter::column_ntotal=buffer_column_number*chain_y_number*100       !particle number of each column for memory
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    !Make wave or not
    integer,parameter::make_wave_or_not=0                                     !make_wave_or_not (0 is no, 1 is yes)

    !Wave type
    character(len=30)::wave_maker_type='PistonType'                           !For soltary wave or others
    !character(len=30)::wave_maker_type='WaveTheory'                          !For regular wave or others
    !character(len=30)::wave_type='Stokes_2nd'                                !For stokes 2nd wave
    !character(len=30)::wave_type='Stokes_3rd'                                !For stokes 3th wave
    !character(len=30)::wave_type='Stokes_5th'                                !For stokes 5th wave
    !character(len=30)::wave_type='Airy'                                        !For Airy wave theory 
    character(len=30)::wave_type='Solitarywave'                               !For Solitary wave combine with Piston Type wavemaker

    !Damping zone for wave maker zone and end zone
    character(len=30)::WavemakerDamping_Switch='off'                           !Damping zone for wave maker 'on' or 'off'
    
    character(len=30)::EndDampingZone_Switch='off'                             !Damping zone for tank end 'on' or 'off'
    !character(len=30)::EndDampingZoneType='Cosin_Type'                        !Damping coefficient from the cosin function (Zheng's Thesis)
    character(len=30)::EndDampingZoneType='OpenFoam_Type'                      !Damping coefficient from the exp function (OpenFoam Type)
    !character(len=30)::EndDampingZoneType='DualSPH_Type'                      !Damping coefficient from the DualSPH(Paper) 

    !Variables for damping zone
    real(kind=8)::EndZoneDampingLengh=3.0                                      !End Zone Damping Lengh
    real(kind=8)::WavemakerDampingLengh=1.0                                    !Here we just mention the initial value and Wave maker Damping Lengh (one wave length if used)
    real(kind=8)::WavemakerDampingLenghScale=1.0                               !The scale of wave length for the Wavemaker Damping Lengh
    real(kind=8)::EndZoneDampingAmplitude=0.3                                  !damping coefficient amplitude              
    !------------------------------------------------------------------------
    
    !************************************************************************                                                       
    !Time iteration Control
    !Recommand: 1, Eular is the best choice for debug for its low time consuming
    !           2, The time step, dt, can be 2 or 3 times of which for Eular when using RUGGER KUTA
    !           3, RUGGER KUTA can make the simulation more stable
    integer,parameter::pair_n=60*subdomain_ntotal                            !pair number
    integer::leapforg_or_eular=3                                             !time iteration type: 1 is leapforg; 2 is Eular; other is RUGGER KUTA
    integer,parameter::max_time_step=50000                                   !max time step
    integer,parameter::standing_time_step=    0                              !time step to release shocking
    integer,parameter::wave_maker_soft_time=10000                            !wave maker soft time
    integer,parameter::save_time_step=500                                    !RESULTS save time step
    integer,parameter::save_time_step_for_debug=10*save_time_step            !RESULTS save_time_step_for_debug
    real(kind=8)::dt=0.0001d0                                                !time step length, dt

    character(len=30)::OutputForTecplot_Switch='on'                          !Output Files For Tecplot ('on' or 'off')
    character(len=30)::OutputForParaview_Switch='off'                        !Output Files For Tecplot ('on' or 'off')
    
    character(len=30)::OutputForPostProcedingGrid_Switch='on'                !Output Files For Post Proceding Grid ('on' or 'off')

    !Interpolation Control
    !1---Normalized; 2---Updated_Spherd; 3---Liner MLS; 6---Quadric MLS
    ! Recommand: 1; 2 maybe have some errors during the calculation
    integer,parameter::MLS_order_Grid=6                                      !MLS order (m=1, 2, 3, 6) 

    !Refresh the domain by MLS or not
    !It may casue uncertainty for the regular linear wave surface
    !I Recommand close this opinion when simulating the regular wave
    character(len=30)::RefreshDensityByMLS_Switch='on'                      !Refresh denisity by MLS 'on' or 'off'
    !1---Normalized; 3---Liner MLS; 6---Quadric MLS
    !Recommand: 1, Normalized with every 100 steps, results are good and compute comsuming is low;
    !Recommand: 2, Also you can try Quadric MLS for its high accuracy, the reuslts are meet well with the crest of the pressure;
    !Recommand: 3, Linear MLS also good in same cases, but it is the last choice in my opinion.       
    !Attention: The refresh domain density may cause the free surface going up automaticaly
    !           It should be closed for wave making
    integer,parameter::RefreshDensity_MLS_Order=6                            !MLS order for refresh density(m=1, 3, 6) 
    integer::MLS_time_step=20                                                !MLS time step for refresh the domain density

    !Sampling Control
    !The Pressure sampling and wave probe position is defined in "subroutine Initialize_SPH_Calculation"
    integer,parameter::sampling_time_step=100                                !Sampling Frequency
    integer,parameter::Sampling_Point_number=0                               !Sampling Point(Maximun is 100)
    integer,parameter::Normalized_Pressure_Sampling_File_Initial_index=1000  !Normalized Pressure Sampling File Initial index
    integer,parameter::Average_Pressure_Sampling_File_Initial_index=2000     !Average Pressure Sampling File Initial index
    integer,parameter::Kalman_Pressure_Sampling_File_Initial_index=3000      !Kalman Pressure Sampling File Initial index
    integer,parameter::Wave_Probe_number=0                                   !Wave Probe number(Maximun is 20) 
    integer,parameter::Wave_Probe_File_Initial_index=5000                    !Wave Probe File Initial index
    integer,parameter::SamplingLineNumber=3                                  !Sampling Line Number
    
    integer,parameter::MLS_order_Sampling=6                                  !MLS order (m=1, 2, 3, 6) for sampling interpolation 
    

    !Variables FOR particle shift SUBROUTINE
    character(len=30)::ShiftParticle_Switch='on'                             !Shift Particle 'on' or 'off'
    real(kind=8)::shift_belta=0.5d0                                          !Recommand 0.5     
    real(kind=8)::shift_R=0.2d0                                              !Recommand 0.2    
    integer,parameter::shift_n=4                                             !Recommand 4
    real(kind=8),parameter::max_shift=0.25d0*interior_dx                     !max shift(Recommand 0.1*smooth_length)
    
    !Interpolation Control
    !1---Normalized; 2---Updated_Spherd; 3---Liner MLS; 6---Quadric MLS
    ! Recommand: 1; 2 maybe have some errors during the calculation
    integer,parameter::MLS_order_Boundary=1                                  !MLS order (m=1, 2, 3, 6) for boundary interpolation 
    real(kind=8)::slighted_modified_alpha=0.0001                             !Slighted_modified_alpha
    !************************************************************************

      
end module information_module