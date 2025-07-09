!**************************************************************************************************************
!  Subroutine : Initialize_SPH_Calculation
!
!  PURPOSE    : Set the initial information from initial or latest domain data for SPH calculation
!               (contains: particles deviding and backgroud mesh)
!
!  Programer  : Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location   : MUN
!
!  Time       : 2017.3.18
!
!  Copyright  : Memorial University
!
!  Version    : 1.1
!
!  Note       : MPI version: mpich-3.2
!               Fortran version: Inter Fortran (ifort) 18.0.0
!**************************************************************************************************************

subroutine Initialize_SPH_Calculation()

    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Wave_Gereration_Module
    USE Function_module
    USE MPI
    
    implicit none

    !==========================================================================================================
    ! Variables in subroutine
    integer::i,j,k,L,m                                              ! Variables for looping

    !==========================================================================================================

    ! Body of Initialize_SPH_Calculation

    ! Call MPI functions
    call MPI_COMM_RANK( MPI_COMM_WORLD, Current_Processor_ID   , ierror_MPI )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, Total_Processors_Number, ierror_MPI )

    !==========================================================================================================
    ! Initialize particle information array
    call Public_Variable_Initialization()

    ! Initialize the directories and files for data saving
    call Directory_and_File_Initialization()

    ! Initialize the kernel function ( Choose the kernel function and get the smooth length )
    call Kernel_Function_Initialization()


    ! Distribute the initial partciles (Both fluid particle and boundary particles are initialized in following subroutine)
    if ( StartFrom == 0 ) then

        call Initial_Particle_From_First_Step()                     ! Satrt from initial step

    elseif ( StartFrom == 1 ) then

        call Initial_Particle_From_Latest_Step()                    ! Satrt from Latest step

    else
        write(*,*) " The 'startFrom' variable is not right. ('InitialTime' or 'LatestTime')"
    end if


    ! Distribute the mesh for post-proceeding
    if( trim(adjustl(OutputForPostProceedingGrid_Switch))=='on' ) then
        call Distribute_Grid_For_Postproceeding()
    endif

    ! Distribute Backgroud Grid
    ! Input: chain_x_number,chain_y_number,chain_z_number
    ! Output: chain_near_mesh_name
    call Distribute_Backgroud_Grid(dim,&                            ! Dimension of grid
                                   chain_x_number,&                 ! Grid number in x direction
                                   chain_y_number,&                 ! Grid number in y direction
                                   chain_z_number,&                 ! Grid number in z direction
                                   chain_max_number,&               ! Max grid number (grid_x_number*grid_y_number*grid_z_number)
                                   mesh_center_positon,&            ! Grid center positon
                                   chain_near_mesh_number,&         ! Near mesh number
                                   chain_near_mesh_name&            ! grid near mesh name
                                   )
    
    ! Initialized the variables for pressure Sampling
    ! Define the sampling points and wave probe position
    ! For 2D simulation, we just need (x,y) position for pressure Sampling and (x) for wave probe 
    ! For 3D simulation, we just need (x,y,z) position for pressure Sampling and (x,y) for wave probe 
    call Sampling_and_Wave_Probe_Initialization()
    

    ! Check the set-up for SPH calculation which contains:(1) CFL checking for SPH time iteration 
    !                                                     (2) Local Reynolds number checking for inlet/outlet boundary condition
    !                                                     (3) Wave model checking for wave making
    call Simulation_Setup_Checking()

    ! Initialized the rest variables such as : (1) Smoke value
    !                                          (2) W(dx) (This part has been move to the initialization subroutine)
    !                                          (3) Identity_matrix
    !                                          (4) Initial total kinetic energy
    !                                          (5) analytical results such as for Green vortex case
    call Additional_Initialization()

    ! Initialization operation results outputting contain the SPH simulation base information and initial particle results
    ! The 'Main_Processor' takes charge of outputting
    if (Current_Processor_ID==Main_Processor) then

        call Initialization_Result_Outputting()

    endif
    !==========================================================================================================


end subroutine Initialize_SPH_Calculation