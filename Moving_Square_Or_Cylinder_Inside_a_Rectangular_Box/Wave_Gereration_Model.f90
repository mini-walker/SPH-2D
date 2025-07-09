!**************************************************************************************************************
!  SUBROUTINE : Wave_Gereration_Module
!
!  PURPOSE    : Module for wave generation
!               (contains: the subroutine for wave generation such as: wave model verification and so on)
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


module Wave_Gereration_Module
       
    USE Information_module
    USE Public_variable_module
    USE SPH_subroutine_module
    USE Function_module
    USE MPI

    contains


    !==========================================================================================================
    ! Subroutine Wave_Model_Verification
    ! Check the wave generation setup is acceptiable or not
    ! Provide a suitable wave theroy for wave generation

    subroutine Wave_Model_Verification()

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables in local subroutine
        integer::i,j,k,L,m                                                         ! Variables for looping

        real(kind=8)::Ursell_Number                                                ! Ursell Number
        real(kind=8)::Wave_Steepnes                                                ! Wave steepnes
        real(kind=8)::WaterDepth_over_WaveLength                                   ! water depth over wave length
        real(kind=8)::Assume_WaveLength                                            ! Assume Wave Length
        real(kind=8)::current_time                                                 ! Current time
        real(kind=8)::current_wave_angle                                           ! Current wave angle
        real(kind=8)::current_wave_elevation                                       ! Wave elevation at current particle z direction
        real(kind=8)::current_wave_probe_position                                  ! Current wave probe position
        real(kind=8)::temp                                                         ! Temp value


        ! Variables for file opeartion
        character(len=40)::File_name                                               ! File name
        character(len=40)::Char_present_wave_probe_location                        ! Char present wave probe location

        !------------------------------------------------------------------------------------------------------


        ! Body of subroutine Wave_model_verification

        !******************************************************************************************************

        Solitary_wave_or_not: if (trim(adjustl(wave_type))=='Solitarywave') then


            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! For Solitary wave
            ! Reference: JOHN FENTON 'A ninth-order solution for the solitary wave'
            Solitary_Epsilon = wave_height/water_depth
            Solitary_alpha   = sqrt(0.75*Solitary_Epsilon)*(1-0.625*Solitary_Epsilon+71/128.0*Solitary_Epsilon**2)

            WavemakerDampingLengh=WavemakerDampingLenghScale*wave_length

            ! The main processor output the data
            if (Current_Processor_ID==Main_Processor) then

                write(*,*) " ====================================================================================="
                write(*,*) "Wave model is 'Third order solitary wave' !"
                write(*,*) "Wave Height   :",wave_Height,"(m)"
                write(*,*) "Wave Length   :",Wave_Length,"(m)---initial value"
                write(*,*) " ====================================================================================="

            endif
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            
        else


            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Linear wave theory 
            ! Calculate the basic infomration for the input wave for checking the wave generation theory
            ! wave length; wave period; wave frequency; wave height; wave number
            if(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WavePeriod') then
                
                Wave_Frequency=2*PI/wave_Period

                !----------------------------------------------------------------------------------------------
                ! Calculate the wave length
                Assume_WaveLength=1.0d0

                do
                    Wave_Length=g*wave_Period**2/(2*PI)*tanh(2*PI*water_depth/Assume_WaveLength)

                    !write(*,*) Assume_WaveLength,Wave_Length
                    ! When the Assume wave length is nearly same with Wave Length result, exit
                    if (abs(Wave_Length-Assume_WaveLength)<1.0E-8) then
                        exit
                    else
                        Assume_WaveLength=Wave_Length
                    endif

                end do
                !----------------------------------------------------------------------------------------------

                Wave_number=2*PI/Wave_Length
                Wave_amplitue=wave_Height/2.0
            
            elseif(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WaveLength') then
                
                !----------------------------------------------------------------------------------------------
                
                Wave_number    = 2*PI/Wave_Length
                Wave_Frequency = sqrt(g*Wave_number*tanh(Wave_number*water_depth))
                wave_Period    = 2*PI/Wave_Frequency
                Wave_amplitue  = wave_Height/2.0

                !----------------------------------------------------------------------------------------------

            elseif(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WaveFreuency') then
                
                !----------------------------------------------------------------------------------------------
                
                wave_Period   = 2*PI/Wave_Frequency
                Wave_amplitue = wave_Height/2.0

                !----------------------------------------------------------------------------------------------
                ! Calculate the wave length
                Assume_WaveLength=1.0d0

                do
                    Wave_Length=g*wave_Period**2/(2*PI)*tanh(2*PI*water_depth/Assume_WaveLength)

                    ! When the Assume wave length is nearly same with Wave Length result, exit
                    if ( abs(Wave_Length-Assume_WaveLength)<1.0E-8 ) then
                        exit
                    else
                        Assume_WaveLength=Wave_Length
                    endif

                end do
                !----------------------------------------------------------------------------------------------

                Wave_number   = 2*PI/Wave_Length
                Wave_amplitue = wave_Height/2.0

                !----------------------------------------------------------------------------------------------

            else

                if (Current_Processor_ID==Main_Processor) then
                     write(*,'(A)') " The Input Wave Information type is not include in the code, please check the type."
                endif

            endif
            !--------------------------------------------------------------------------------------------------
         
            !--------------------------------------------------------------------------------------------------
            ! Calculate the h/L for the wave checking: deep water wave; shallow water wave and finite water wave
            WaterDepth_over_WaveLength=water_depth/Wave_Length
            Ursell_Number=wave_Height/water_depth*(Wave_Length/water_depth)**2
            Wave_Steepnes=wave_Height/Wave_Length

            ! The main processor output the data
            if (Current_Processor_ID==Main_Processor) then

                !----------------------------------------------------------------------------------------------
                write(*,*) "WaterDepth/WaveLength:",WaterDepth_over_WaveLength
                write(*,*) "Ursell Number:",Ursell_Number

                if (WaterDepth_over_WaveLength>=0.5) then

                    write(*,*) "Wave Type: Deep water waves"

                elseif (WaterDepth_over_WaveLength<=0.05) then

                    write(*,*) "Wave Type: Shallow water waves"

                else

                    write(*,*) "Wave Type: Finite water waves"
                    
                endif
                !----------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------
                ! Calculate the Ursell Number for different model
                ! Ursell=H/h*(λ/h)^2
                ! H : the wave height, i.e. the difference between the elevations of the wave crest and trough;
                ! h : the mean water depth
                ! λ : the wavelength, which has to be large compared to the depth

                if (Ursell_Number>=26 .or. WaterDepth_over_WaveLength<0.1) then

                    write(*,'(A)') " Recommand wave theory: Cnoidal Wave Theory. (Based on Ursell Number and Steepness)"

                elseif (Ursell_Number<26) then

                    if (Wave_Steepnes<=0.006) then

                        write(*,'(A)') " Recommand wave theory: Airy Wave Theory. (Based on Ursell Number and Steepness)"

                    elseif (Wave_Steepnes>0.006 .and. Wave_Steepnes<=0.04) then 
                    
                        write(*,'(A)') " Recommand wave theory: 2nd or 3th order Stokes Wave Theory. (Based on Ursell Number and Steepness)"

                    elseif (Wave_Steepnes>0.04 .and. Wave_Steepnes<=0.14) then 
                    
                        write(*,'(A)') " Recommand wave theory: 5th order Stokes Wave Theory. (Based on Ursell Number and Steepness)"
                        
                    endif
                    
                endif
                write(*,*) " ====================================================================================="

                !----------------------------------------------------------------------------------------------


                !----------------------------------------------------------------------------------------------
                ! Recommand wave theory based on Second Theory

                ! For Airy wave
                if (wave_length<2*water_depth .and. wave_Height/wave_length<=0.06) then
                    
                    write(*,'(A)') " Recommand wave theory: Airy Wave Theory. (Based on wave length and period)"

                endif

                ! For stokes wave
                if (wave_length<8*water_depth) then

                    if (wave_Height/(g*wave_Period**2)<0.008 .and. water_depth/(g*wave_Period**2)>=0.002) then

                        write(*,'(A)') " Recommand wave theory: 2nd order Stokes Wave Theory. (Based on wave length and period)"

                    endif

                    if (wave_Height/(g*wave_Period**2)>=0.008 .and. wave_Height/(g*wave_Period**2)<=0.02 .and. water_depth/(g*wave_Period**2)>=0.01) then
                        
                        write(*,'(A)') " Recommand wave theory: 3rd order Stokes Wave Theory. (Based on wave length and period)"

                    endif

                    if (wave_Height/(g*wave_Period**2)>=0.01 .and. wave_Height/(g*wave_Period**2)<=0.05 .and. water_depth/(g*wave_Period**2)>=0.01) then
                        
                        write(*,'(A)') " Recommand wave theory: 5th order Stokes Wave Theory. (Based on wave length and period)"

                    endif

                endif

                ! For Cnoidal wave
                if (wave_length>8*water_depth .and. water_depth/wave_Height>1.8 .and. water_depth/wave_Height<2.2 .and. wave_Height/(g*wave_Period**2)<=0.02 .and. water_depth/(g*wave_Period**2)<=0.03) then
                    
                    write(*,'(A)') " Recommand wave theory: Cnoidal Wave Theory. (Based on wave length and period)"

                endif


                write(*,*) " ====================================================================================="

                !----------------------------------------------------------------------------------------------

            end if
            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^






            !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ! Re-calculate the wave information again based on the wave theory you used
            ! For nonliear wave (3th or 5th Stokes waves), we should calculate them again
            ! For linear wave, we just keep same
            ! For 5th order Stokes, the reference:"Hydrodynamics of Offshore Structures.P59"
            if(trim(adjustl(wave_type))=='Stokes_5th') then

                !----------------------------------------------------------------------------------------------
                ! Calculate the coefficients for 5th order for stokes wave
                Fifth_stokes_A=0.0d0  
                Fifth_stokes_B=0.0d0  
                Fifth_stokes_C=0.0d0 
                Fifth_stokes_D=0.0d0   
                Fifth_stokes_E=0.0d0 
                !----------------------------------------------------------------------------------------------

                if(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WavePeriod') then
                    
                    Wave_Frequency=2*PI/wave_Period

                    !------------------------------------------------------------------------------------------
                    !Determine wave number through Equ.(3.89) and Equ.(3.90)
                    
                    GET_k_1:do i=1,1000000

                        Last_wave_number=wave_number                                      !wave number in last step

                        !---------------------------------------------------------------------------------------------------------
                        !Define the k first
                        !Coefficients (Deep water wave is the initial condition)
                        Stokes_S=sinh(wave_number*water_depth)
                        Stokes_C=cosh(wave_number*water_depth)

                        Fifth_stokes_C(1)=(8*Stokes_C**4-8*Stokes_C**2+9)/(8*Stokes_S**4) 
                        Fifth_stokes_C(2)=(3840*Stokes_C**12-4096*Stokes_C**10-2592*Stokes_C**8-1008*Stokes_C**6+5944*Stokes_C**4-1830*Stokes_C**2+147)/(512*Stokes_S**10*(6*Stokes_C**2-1)) 
                        Fifth_stokes_C(3)=(-1.0)/(4*Stokes_S*Stokes_C) 
                        Fifth_stokes_C(4)=(12*Stokes_C**8+36*Stokes_C**6-162*Stokes_C**4+141*Stokes_C**2-27)/(192*Stokes_C*Stokes_C**9) 
                        

                        !GET the lamda coeffcient in fifth stokes wave theory (Equ.(3.90) )
                        !initial wave number is know, calculate wave number

                        if (i==1) then
                            Stokes_Lamda=0.0d0                                             !initial lamda value
                        endif

                        do  L=1,1000000

                            Temp_wave_number=wave_number

                            wave_number=Wave_Frequency**2/(g*tanh(wave_number*water_depth)*(1+Stokes_Lamda**2*Fifth_stokes_C(1)+Stokes_Lamda**4*Fifth_stokes_C(2)))

                            if (abs(wave_number-Temp_wave_number)<=1.0e-8 ) exit
                            
                        enddo
                        !---------------------------------------------------------------------------------------------------------

                        !---------------------------------------------------------------------------------------------------------
                        !Renew the s and c based on the new k value
                        !Coefficients (Deep water wave is the initial condition)
                        Stokes_S=sinh(wave_number*water_depth)
                        Stokes_C=cosh(wave_number*water_depth)

                        Fifth_stokes_B(2,2)=Stokes_C*(2*Stokes_C**2+1)/(4*Stokes_S**3)
                        Fifth_stokes_B(2,4)=Stokes_C*(272*Stokes_C**8-504*Stokes_C**6+192*Stokes_C**4+322*Stokes_C**2+21)/(384*Stokes_S**9)
                        Fifth_stokes_B(3,3)=3*(8*Stokes_C**6+1)/(64*Stokes_S**6)
                        Fifth_stokes_B(3,5)=(88128*Stokes_C**14-208224*Stokes_C**12+70848*Stokes_C**10+54000*Stokes_C**8-21816*Stokes_C**6+6264*Stokes_C**4-54*Stokes_C**2-81)/(12288*Stokes_S**12*(6*Stokes_C**2-1))
                        Fifth_stokes_B(4,4)=Stokes_C*(768*Stokes_C**10-448*Stokes_C**8-48*Stokes_C**6+48*Stokes_C**4+106*Stokes_C**2-21)/(384*Stokes_S**9*(6*Stokes_C**2-1))
                        Fifth_stokes_B(5,5)=(192000*Stokes_C**16-262720*Stokes_C**14+83680*Stokes_C**12+20160*Stokes_C**10-7280*Stokes_C**8+7160*Stokes_C**6-1800*Stokes_C**4-1050*Stokes_C**2+255)/(12288*Stokes_S**10*(6*Stokes_C**2-1)*(8*Stokes_C**4-11*Stokes_C**2+3))


                        !Re-calculate the Lamda (Equ.(3.89))
                        do m=1,1000000
                           
                            Temp_Stokes_Lamda=Stokes_Lamda
                            Stokes_Lamda=wave_number*wave_Height/2.0-Fifth_stokes_B(3,3)*Stokes_Lamda**3-(Fifth_stokes_B(3,5)+Fifth_stokes_B(5,5))*Stokes_Lamda**5
                            
                            if (abs(Temp_Stokes_Lamda-Stokes_Lamda)<1.0E-8) exit

                        end do


                        if (abs(Last_wave_number-wave_number)<1.0E-8) exit 
                            
                        !---------------------------------------------------------------------------------------------------------

                        !Check the results
                        if (Current_Processor_ID==Main_Processor .and. i==1000000) then

                            write(*,'(A)') ' Errors in 5th order stokes wave information calculation. (Wave_Model_Verification)'

                        end if

                    enddo GET_k_1
                    !------------------------------------------------------------------------------------------

                    Wave_Length=2*PI/Wave_number
                    Wave_celerity=2*PI/(wave_number*wave_Period)                             !c=2PI/(KT)
                    Wave_amplitue=wave_Height/2.0
                    Wave_celerity=2*PI/(wave_number*wave_Period)                             !c=2PI/(KT)
                    
                    Wave_celerity_0=sqrt((g/wave_number)*tanh(wave_number*water_depth))     !c=sqrt(g/K*tanh(kd))
                
                elseif(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WaveLength') then

                    !------------------------------------------------------------------------------------------
                    Wave_number=2*PI/Wave_Length                                             !k=2PI/L
                    Wave_celerity_0=sqrt((g/wave_number)*tanh(wave_number*water_depth))     !c=sqrt(g/K*tanh(kd))

                    !Coefficients (Deep water wave is the initial condition)
                    Stokes_S=sinh(wave_number*water_depth)
                    Stokes_C=cosh(wave_number*water_depth)

                    Fifth_stokes_B(2,2)=Stokes_C*(2*Stokes_C**2+1)/(4*Stokes_S**3)
                    Fifth_stokes_B(2,4)=Stokes_C*(272*Stokes_C**8-504*Stokes_C**6+192*Stokes_C**4+322*Stokes_C**2+21)/(384*Stokes_S**9)
                    Fifth_stokes_B(3,3)=3*(8*Stokes_C**6+1)/(64*Stokes_S**6)
                    Fifth_stokes_B(3,5)=(88128*Stokes_C**14-208224*Stokes_C**12+70848*Stokes_C**10+54000*Stokes_C**8-21816*Stokes_C**6+6264*Stokes_C**4-54*Stokes_C**2-81)/(12288*Stokes_S**12*(6*Stokes_C**2-1))
                    Fifth_stokes_B(4,4)=Stokes_C*(768*Stokes_C**10-448*Stokes_C**8-48*Stokes_C**6+48*Stokes_C**4+106*Stokes_C**2-21)/(384*Stokes_S**9*(6*Stokes_C**2-1))
                    Fifth_stokes_B(5,5)=(192000*Stokes_C**16-262720*Stokes_C**14+83680*Stokes_C**12+20160*Stokes_C**10-7280*Stokes_C**8+7160*Stokes_C**6-1800*Stokes_C**4-1050*Stokes_C**2+255)/(12288*Stokes_S**10*(6*Stokes_C**2-1)*(8*Stokes_C**4-11*Stokes_C**2+3))

                    !GET the lamda coeffcient in fifth stokes wave theory (Equ.(3.89) )
                    
                    Stokes_Lamda=0.0d0                                             !initial lamda value

                    do m=1,1000000
                       
                        Temp_Stokes_Lamda=Stokes_Lamda
                        Stokes_Lamda=wave_number*wave_Height/2.0-Fifth_stokes_B(3,3)*Stokes_Lamda**3-(Fifth_stokes_B(3,5)+Fifth_stokes_B(5,5))*Stokes_Lamda**5
                        
                        if (abs(Temp_Stokes_Lamda-Stokes_Lamda)<1.0E-8) exit

                    end do
                    !------------------------------------------------------------------------------------------


                    Wave_celerity=Wave_celerity_0*sqrt(1+Fifth_stokes_C(1)*Stokes_Lamda**2+Fifth_stokes_C(2)*Stokes_Lamda**4)

                    Wave_Frequency=wave_number*Wave_celerity           ! w=kc

                    wave_Period=2*PI/Wave_Frequency                    ! T=2PI/w

                    !------------------------------------------------------------------------------------------

                elseif(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WaveFreuency') then
                    
                    !------------------------------------------------------------------------------------------

                    wave_Period=2*PI/Wave_Frequency
                    Wave_amplitue=wave_Height/2.0

                    !------------------------------------------------------------------------------------------
                    !Determine wave number through Equ.(3.89) and Equ.(3.90)
                    
                    GET_k_2:do i=1,1000000

                        Last_wave_number=wave_number                                      !wave number in last step
                        
                        !---------------------------------------------------------------------------------------------------------
                        !Define the k first
                        !Coefficients (Deep water wave is the initial condition)
                        Stokes_S=sinh(wave_number*water_depth)
                        Stokes_C=cosh(wave_number*water_depth)

                        Fifth_stokes_C(1)=(8*Stokes_C**4-8*Stokes_C**2+9)/(8*Stokes_S**4) 
                        Fifth_stokes_C(2)=(3840*Stokes_C**12-4096*Stokes_C**10-2592*Stokes_C**8-1008*Stokes_C**6+5944*Stokes_C**4-1830*Stokes_C**2+147)/(512*Stokes_S**10*(6*Stokes_C**2-1)) 
                        Fifth_stokes_C(3)=(-1.0)/(4*Stokes_S*Stokes_C) 
                        Fifth_stokes_C(4)=(12*Stokes_C**8+36*Stokes_C**6-162*Stokes_C**4+141*Stokes_C**2-27)/(192*Stokes_C*Stokes_C**9) 
                        

                        !GET the lamda coeffcient in fifth stokes wave theory (Equ.(3.90) )
                        !initial wave number is know, calculate wave number

                        if (i==1) then
                            Stokes_Lamda=0.0d0                                             !initial lamda value
                        endif

                        do  L=1,1000000

                            Temp_wave_number=wave_number

                            wave_number=Wave_Frequency**2/(g*tanh(wave_number*water_depth)*(1+Stokes_Lamda**2*Fifth_stokes_C(1)+Stokes_Lamda**4*Fifth_stokes_C(2)))

                            if (abs(wave_number-Temp_wave_number)<=1.0e-8 ) exit
                            
                        enddo
                        !---------------------------------------------------------------------------------------------------------

                        !---------------------------------------------------------------------------------------------------------
                        !Renew the s and c based on the new k value
                        !Coefficients (Deep water wave is the initial condition)
                        Stokes_S=sinh(wave_number*water_depth)
                        Stokes_C=cosh(wave_number*water_depth)

                        Fifth_stokes_B(2,2)=Stokes_C*(2*Stokes_C**2+1)/(4*Stokes_S**3)
                        Fifth_stokes_B(2,4)=Stokes_C*(272*Stokes_C**8-504*Stokes_C**6+192*Stokes_C**4+322*Stokes_C**2+21)/(384*Stokes_S**9)
                        Fifth_stokes_B(3,3)=3*(8*Stokes_C**6+1)/(64*Stokes_S**6)
                        Fifth_stokes_B(3,5)=(88128*Stokes_C**14-208224*Stokes_C**12+70848*Stokes_C**10+54000*Stokes_C**8-21816*Stokes_C**6+6264*Stokes_C**4-54*Stokes_C**2-81)/(12288*Stokes_S**12*(6*Stokes_C**2-1))
                        Fifth_stokes_B(4,4)=Stokes_C*(768*Stokes_C**10-448*Stokes_C**8-48*Stokes_C**6+48*Stokes_C**4+106*Stokes_C**2-21)/(384*Stokes_S**9*(6*Stokes_C**2-1))
                        Fifth_stokes_B(5,5)=(192000*Stokes_C**16-262720*Stokes_C**14+83680*Stokes_C**12+20160*Stokes_C**10-7280*Stokes_C**8+7160*Stokes_C**6-1800*Stokes_C**4-1050*Stokes_C**2+255)/(12288*Stokes_S**10*(6*Stokes_C**2-1)*(8*Stokes_C**4-11*Stokes_C**2+3))


                        !Re-calculate the Lamda (Equ.(3.89))
                        do m=1,1000000
                           
                            Temp_Stokes_Lamda=Stokes_Lamda
                            Stokes_Lamda=wave_number*wave_Height/2.0-Fifth_stokes_B(3,3)*Stokes_Lamda**3-(Fifth_stokes_B(3,5)+Fifth_stokes_B(5,5))*Stokes_Lamda**5
                            
                            if (abs(Temp_Stokes_Lamda-Stokes_Lamda)<1.0E-8) exit

                        end do


                        if (abs(Last_wave_number-wave_number)<1.0E-8) exit 
                            
                        !---------------------------------------------------------------------------------------------------------

                        !Check the results
                        if (Current_Processor_ID==Main_Processor .and. i==1000000) then

                            write(*,'(A)') ' Errors in 5th order stokes wave information calculation. (Wave_Model_Verification)'

                        end if

                    enddo GET_k_2
                    !------------------------------------------------------------------------------------------

                    !------------------------------------------------------------------------------------------
                    Wave_Length=2*PI/Wave_number
                    Wave_celerity=2*PI/(wave_number*wave_Period)                             ! c=2PI/(KT)

                    Wave_celerity_0=sqrt((g/wave_number)*tanh(wave_number*water_depth))      ! c=sqrt(g/K*tanh(kd))
                    !------------------------------------------------------------------------------------------

                else

                    if (Current_Processor_ID==Main_Processor) then
                         write(*,*) "The Input Wave Information type is not include in the code, please check the type."
                    endif

                endif
                !----------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------
                ! Coefficients (Deep water wave is the initial condition)
                ! "Hydrodynamics of Offshore Structures.P59" and "Non-linear waves (Brorsen, Michael)"
                Fifth_stokes_A=0.0d0  
                Fifth_stokes_B=0.0d0  
                Fifth_stokes_C=0.0d0 
                Fifth_stokes_D=0.0d0  
                Fifth_stokes_E=0.0d0

                Stokes_S=sinh(wave_number*water_depth)
                Stokes_C=cosh(wave_number*water_depth)

                Fifth_stokes_A(1,1)=1.0/Stokes_S
                Fifth_stokes_A(1,3)=(-Stokes_C**2*(5*Stokes_C**2+1))/(8*Stokes_S**5)
                Fifth_stokes_A(1,5)=-(1184*Stokes_C**10-1440*Stokes_C**8-1992*Stokes_C**6+2641*Stokes_C**4-249*Stokes_C**2+18)/(1536*Stokes_S**11)
                Fifth_stokes_A(2,2)=3.0/(8*Stokes_S**4)
                Fifth_stokes_A(2,4)=(192*Stokes_C**8-424*Stokes_C**6-312*Stokes_C**4+480*Stokes_C**2-17)/(768*Stokes_S**10)
                Fifth_stokes_A(3,3)=(13-4*Stokes_C**2)/(64*Stokes_S**7)
                Fifth_stokes_A(3,5)=(512*Stokes_C**12-4224*Stokes_C**10-6800*Stokes_C**8-12808*Stokes_C**6+16704*Stokes_C**4-3154*Stokes_C**2+107)/(4096*Stokes_S**13*(6*Stokes_C**2-1))
                Fifth_stokes_A(4,4)=(80*Stokes_C**6-816*Stokes_C**4+1338*Stokes_C**2-197)/(1536*Stokes_S**10*(6*Stokes_C**2-1))
                Fifth_stokes_A(5,5)=-(2880*Stokes_C**10-72480*Stokes_C**8+324000*Stokes_C**6-432000*Stokes_C**4+163470*Stokes_C**2-16245)/(61440*Stokes_S**11*(6*Stokes_C**2-1)*(8*Stokes_C**4-11*Stokes_C**2+3))
                
                Fifth_stokes_B(2,2)=Stokes_C*(2*Stokes_C**2+1)/(4*Stokes_S**3)
                Fifth_stokes_B(2,4)=Stokes_C*(272*Stokes_C**8-504*Stokes_C**6+192*Stokes_C**4+322*Stokes_C**2+21)/(384*Stokes_S**9)
                Fifth_stokes_B(3,3)=3*(8*Stokes_C**6+1)/(64*Stokes_S**6)
                Fifth_stokes_B(3,5)=(88128*Stokes_C**14-208224*Stokes_C**12+70848*Stokes_C**10+54000*Stokes_C**8-21816*Stokes_C**6+6264*Stokes_C**4-54*Stokes_C**2-81)/(12288*Stokes_S**12*(6*Stokes_C**2-1))
                Fifth_stokes_B(4,4)=Stokes_C*(768*Stokes_C**10-448*Stokes_C**8-48*Stokes_C**6+48*Stokes_C**4+106*Stokes_C**2-21)/(384*Stokes_S**9*(6*Stokes_C**2-1))
                Fifth_stokes_B(5,5)=(192000*Stokes_C**16-262720*Stokes_C**14+83680*Stokes_C**12+20160*Stokes_C**10-7280*Stokes_C**8+7160*Stokes_C**6-1800*Stokes_C**4-1050*Stokes_C**2+255)/(12288*Stokes_S**10*(6*Stokes_C**2-1)*(8*Stokes_C**4-11*Stokes_C**2+3))

                Fifth_stokes_C(1)=(8*Stokes_C**4-8*Stokes_C**2+9)/(8*Stokes_S**4) 
                Fifth_stokes_C(2)=(3840*Stokes_C**12-4096*Stokes_C**10-2592*Stokes_C**8-1008*Stokes_C**6+5944*Stokes_C**4-1830*Stokes_C**2+147)/(512*Stokes_S**10*(6*Stokes_C**2-1)) 
                Fifth_stokes_C(3)=(-1.0)/(4*Stokes_S*Stokes_C) 
                Fifth_stokes_C(4)=(12*Stokes_C**8+36*Stokes_C**6-162*Stokes_C**4+141*Stokes_C**2-27)/(192*Stokes_C*Stokes_C**9) 

                Fifth_stokes_D(1)=Stokes_Lamda*Fifth_stokes_A(1,1)+Stokes_Lamda**3*Fifth_stokes_A(1,3)+Stokes_Lamda**5*Fifth_stokes_A(1,5)
                Fifth_stokes_D(2)=Stokes_Lamda**2*Fifth_stokes_A(2,2)+Stokes_Lamda**4*Fifth_stokes_A(2,4)
                Fifth_stokes_D(3)=Stokes_Lamda**3*Fifth_stokes_A(3,3)+Stokes_Lamda**5*Fifth_stokes_A(3,5)
                Fifth_stokes_D(4)=Stokes_Lamda**4*Fifth_stokes_A(4,4)
                Fifth_stokes_D(5)=Stokes_Lamda**5*Fifth_stokes_A(5,5)

                Fifth_stokes_E(1)=Stokes_Lamda
                Fifth_stokes_E(2)=Stokes_Lamda**2*Fifth_stokes_B(2,2)+Stokes_Lamda**4*Fifth_stokes_B(2,4)
                Fifth_stokes_E(3)=Stokes_Lamda**3*Fifth_stokes_B(3,3)+Stokes_Lamda**5*Fifth_stokes_B(3,5)
                Fifth_stokes_E(4)=Stokes_Lamda**4*Fifth_stokes_B(4,4)
                Fifth_stokes_E(5)=Stokes_Lamda**5*Fifth_stokes_B(5,5)
                !----------------------------------------------------------------------------------------------

                !----------------------------------------------------------------------------------------------
                ! For pressure calculation
                Stokes_pressure_S=1.0/cosh(2*wave_number*water_depth)

                Stokes_pressure_C_0=sqrt(tanh(wave_number*water_depth))
                Stokes_Epsilon=0.5*wave_number*wave_height
                Fifth_stokes_pressure_E=0

                Fifth_stokes_pressure_E(2)=tanh(wave_number*water_depth)*(2+2*Stokes_pressure_S+5*Stokes_pressure_S**2)/(4*(1-Stokes_pressure_S)**2)
                Fifth_stokes_pressure_E(4)=tanh(wave_number*water_depth)*(8+12*Stokes_pressure_S-152*Stokes_pressure_S**2-308*Stokes_pressure_S**3-42*Stokes_pressure_S**4+77*Stokes_pressure_S**5)/(32*(1-Stokes_pressure_S)**5)
                        
                Bernoulli_constant_R=g/wave_number*(0.5*Stokes_pressure_C_0**2+wave_number*water_depth+Stokes_Epsilon**2*Fifth_stokes_pressure_E(2)+Stokes_Epsilon**4*Fifth_stokes_pressure_E(4))

                !----------------------------------------------------------------------------------------------

            else if(trim(adjustl(wave_type))=='Stokes_3rd') then

                !----------------------------------------------------------------------------------------------
                ! Calculate the wave number based on Equ(3.86) in "Hydrodynamics of Offshore Structures.P58"
                ! c^2=c_0^2*(1+(PI*H/L)^2*((9-8*coshkd^2+8*coshkd^4)/(8*sinhkd^4)))

                !----------------------------------------------------------------------------------------------
                ! Use the wave number based on linear theory
                ! wave celerity from the linear wave theory
                Wave_celerity_0=sqrt(g/wave_number*tanh(wave_number*water_depth))

                Wave_celerity=Wave_celerity_0*sqrt(1+(PI*wave_height/wave_length)**2*(9-8*(cosh(wave_number*water_depth))**2+8*(cosh(wave_number*water_depth))**4)/(8*(sinh(wave_number*water_depth))**4))

                ! Calculate the coefficients for third order stokes wave
                Stokes_Epsilon=(wave_number*wave_height)/2.0

                Stokes_S=sinh(wave_number*water_depth)
                Stokes_C=cosh(wave_number*water_depth)

                !Calculate mean water deepth
                water_depth_bar=water_depth+Stokes_Epsilon**2*Stokes_C/(2*wave_number*Stokes_S)
                
                Stokes_S=sinh(wave_number*water_depth_bar)
                Stokes_C=cosh(wave_number*water_depth_bar)

                Third_stokes_D=0.0d0
                Third_stokes_E=0.0d0

                Third_stokes_E(1)=Stokes_Epsilon*(1-Stokes_Epsilon**2*3*(8*Stokes_S**6+24*Stokes_S**4+24*Stokes_S**2+9)/(64*Stokes_S**6))
                Third_stokes_E(2)=Stokes_Epsilon**2*((Stokes_C*(2*Stokes_S**2+3))/(4*Stokes_S**3))
                Third_stokes_E(3)=Stokes_Epsilon**3*(3*(8*Stokes_S**6+24*Stokes_S**4+24*Stokes_S**2+9)/(64*Stokes_S**6))

                Third_stokes_D(1)=Stokes_Epsilon*(1-Stokes_Epsilon**2*(64*Stokes_S**6+160*Stokes_S**4+120*Stokes_S**2+27)/(64*Stokes_S**6))/Stokes_S
                Third_stokes_D(2)=Stokes_Epsilon**2*(3.0/(8*Stokes_S**4))
                Third_stokes_D(3)=Stokes_Epsilon**3*((-4*Stokes_S**2+9)/(64*Stokes_S**7))

                !----------------------------------------------------------------------------------------------

            else if(trim(adjustl(wave_type))=='Cnoidal_5th') then



            endif

            !--------------------------------------------------------------------------------------------------
            ! Output the basic wave infomration (Main_Processor)
            if (Current_Processor_ID==Main_Processor) then

                write(*,*) "==============================================================="
                
                if (trim(adjustl(wave_type))=='Airy') then
                    write(*,*) "Results based on Linear wave theory. (Airy)"
                else if (trim(adjustl(wave_type))=='Stokes_2nd') then
                    write(*,*) "Results based on Non-linear wave theory. (2nd order stokes)"
                else if (trim(adjustl(wave_type))=='Stokes_3rd') then
                    write(*,*) "Results based on Non-linear wave theory. (3th order stokes)"
                else if (trim(adjustl(wave_type))=='Stokes_5th') then
                    write(*,*) "Results based on Non-linear wave theory. (5th order stokes)"
                endif
                
                write(*,*) "Wave Height   :",wave_Height,"(m)"
                write(*,*) "Wave Amplitude:",Wave_amplitue,"(m)"
                write(*,*) "Wave Length   :",Wave_Length,"(m)"
                write(*,*) "Wave Number   :",Wave_number,"(rad/m)"
                write(*,*) "Wave Period   :",wave_Period,"(s)"
                write(*,*) "Wave Frequency:",Wave_Frequency,"(rad/s)"
                write(*,*) "==============================================================="

            endif
            !--------------------------------------------------------------------------------------------------

        endif Solitary_wave_or_not
        !======================================================================================================












        ! Assign the wave maker damping Lengh based on real wave length calculate from the wavemaker subroutine
        WavemakerDampingLengh=WavemakerDampingLenghScale*wave_length


        !======================================================================================================
        ! Generate the analytical solution for regular waves
        
        if (Current_Processor_ID==Main_Processor) then

            do k=1,Wave_Probe_number

                !transfer the Current Wave Probe Location from real to character
                write(Char_present_wave_probe_location,'(f6.3)') Wave_Probe_Position(k,1)

                File_name=trim(adjustl(wave_type))//'_'//trim(adjustl(Char_present_wave_probe_location))//"_Analytical_Results.dat"

                !Open the saving file
                open(unit=1,file=trim(adjustl(File_name)),status="replace",position="rewind",action="write",iostat=ioerror)    
                
                if(ioerror==0) then
                   
                    write(*,*) "The open action of the "// trim(adjustl(File_name))//" is successful!"
                    
                    !Calculate the wave elevations
                    if(trim(adjustl(wave_type))=='Stokes_2nd') then

                        !---------------------------------------------------------------------------------------------------------
                        !Wave_elevation=H/2*cos(kx-wt+Fai0)+(PI*H)/8*(H/L)*cosh(kh)[cos(2kh)+2]*cos(2(kx-wt+Fai0))/(sinh(kh))**3
                        !Fai0: wave initial phase angle;
                        !H: wave height;
                        !k: wave number;
                        !w: wave frequency;
                        !t: time;
                        !L: wave length;
                        !h: water depth.

                        current_wave_probe_position=Wave_Probe_Position(k,1)               ! Unit: m

                        do i=1,max_time_step,50

                            current_time=dt*i

                            current_wave_angle=wave_number*current_wave_probe_position-wave_frequency*current_time+Wave_initial_phase_angle
                            current_wave_elevation=water_depth+0.5*wave_height*cos(current_wave_angle)+PI*wave_height*wave_height/(8*wave_length)*cosh(wave_number*water_depth)/(sinh(wave_number*water_depth))**3*(2+cosh(2*wave_number*water_depth))*cos(2*current_wave_angle)
                 
                            write(1,*) current_time,current_wave_elevation

                        enddo

                    elseif(trim(adjustl(wave_type))=='Airy') then
                        
                        !---------------------------------------------------------------------------------------------------------
                        !Wave_elevation=H/2*cos(kx-wt+Fai0)

                        current_wave_probe_position=Wave_Probe_Position(k,1)               ! Unit: m

                        do i=1,max_time_step,50

                            current_time=dt*i

                            current_wave_angle=wave_number*current_wave_probe_position-wave_frequency*current_time+Wave_initial_phase_angle
                            current_wave_elevation=water_depth+0.5*wave_height*cos(current_wave_angle)
                 
                            write(1,*) current_time,current_wave_elevation

                        enddo

                    elseif(trim(adjustl(wave_type))=='Stokes_5th') then

                        !---------------------------------------------------------------------------------------------------------
                        ! Wave_elevation=H/2*cos(kx-wt+Fai0)

                        current_wave_probe_position=Wave_Probe_Position(k,1)               ! Unit: m

                        ! !--------------------------------------------------------------------------------------
                        ! !Fenton "Nonlinear Wave Theories"
                        ! do i=1,max_time_step,50

                        !     current_time=dt*i

                        !     current_wave_angle=wave_number*current_wave_probe_position-wave_frequency*current_time+Wave_initial_phase_angle
                            
                        !     current_wave_elevation=0.0d0
                        !     do m=1,5
                        !         temp=0.0d0
                        !         do L=1,m
                        !             temp=temp+Fifth_stokes_B(m,L)*cos(L*current_wave_angle)
                        !         enddo
                        !         current_wave_elevation=current_wave_elevation+Fifth_stokes_Lamda**m*temp
                        !     enddo
                        !     current_wave_elevation=water_depth+current_wave_elevation/wave_number
                 
                        !     write(1,*) current_time,current_wave_elevation

                        ! end do 
                        ! !--------------------------------------------------------------------------------------

                        !---------------------------------------------------------------------------------------------------------
                        ! Brorsen, Michael "Nonlinear Wave Theories"
                        do i=1,max_time_step,50

                            current_time=dt*i

                            current_wave_angle=wave_number*current_wave_probe_position-wave_frequency*current_time+Wave_initial_phase_angle
                            
                            current_wave_elevation=0.0d0
                            do m=1,5
                                current_wave_elevation=current_wave_elevation+Fifth_stokes_E(m)*cos(m*current_wave_angle)
                            enddo
                            current_wave_elevation=water_depth+current_wave_elevation/wave_number
                 
                            write(1,*) current_time,current_wave_elevation

                        enddo 
                        !---------------------------------------------------------------------------------------------------------

                    elseif(trim(adjustl(wave_type))=='Stokes_3rd') then


                        current_wave_probe_position=Wave_Probe_Position(k,1)               ! Unit: m

                        !---------------------------------------------------------------------------------------------------------
                        ! Reference: Song Zhiyao "On the universal third order stokes wave solution"
                        do i=1,max_time_step,50

                            current_time=dt*i

                            current_wave_angle=wave_number*current_wave_probe_position-wave_frequency*current_time+Wave_initial_phase_angle
                            
                            current_wave_elevation=0.0d0
                            do m=1,3
                                current_wave_elevation=current_wave_elevation+Third_stokes_E(m)*cos(m*current_wave_angle)
                            enddo
                            current_wave_elevation=water_depth+current_wave_elevation/wave_number
                 
                            write(1,*) current_time,current_wave_elevation

                        enddo 
                        !---------------------------------------------------------------------------------------------------------

                    elseif(trim(adjustl(wave_type))=='Solitarywave') then

                        current_wave_probe_position=Wave_Probe_Position(k,1)               ! Unit: m

                        !---------------------------------------------------------------------------------------------------------
                        ! Reference: Song Zhiyao "On the universal third order stokes wave solution"
                        do i=1,max_time_step,50

                            current_time=dt*i

                            current_wave_angle=wave_number*current_wave_probe_position-wave_frequency*current_time+Wave_initial_phase_angle
                            
                            current_wave_elevation=0.0d0
                            do m=1,3
                                current_wave_elevation=current_wave_elevation+Third_stokes_E(m)*cos(m*current_wave_angle)
                            enddo
                            current_wave_elevation=water_depth+current_wave_elevation/wave_number
                 
                            write(1,*) current_time,current_wave_elevation

                        enddo 
                        !---------------------------------------------------------------------------------------------------------

                    else

                        write(*,*) "The wave type is not include in the code, please review the type."

                   endif

                else
                    write(*,*) "The open action of the "// trim(adjustl(File_name))//" is fail!"
                endif
               
                close(1)
                !----------------------------------------------------------------------------------------------

            enddo

        endif
        !======================================================================================================



    end subroutine Wave_Model_Verification
    !==========================================================================================================










    !==========================================================================================================
    ! *Plate motion of piston type for solitary wave maker
    ! *FROM THE GORING (1978),
    ! *XP(T)=H/K*[TANH(X(T))+TANH(K/D*LAMDA)],X(T)=K/D*(C*T-XP(T)-LAMDA)
    ! *K=SQRT(3H/4D),C=SQRT(G*(D+H))
    ! *UP(T)=C*H/D*(1/(COSH(X(T))**2.+H/D)
    ! *LAMDA=L/K*D, L=ARCCOSH(1/SQRT(EZ)),EZ=0.002
    subroutine Solitary_wave_maker_plate_motion(WMA,WMV,WMX,TTT,HEIGHT,H0,DT,G)
    
        implicit none
       
        !------------------------------------------------------------------------------------------------------
        ! Variables 
        real(kind=8) :: WMA,WMV,WMX,TTT,HEIGHT,H0,DT,G
        real(kind=8) :: EZ,CL,CK,CK1,C,CLAMDA,CERR,XPI,XP,XT,ERR,UP,VP
        !------------------------------------------------------------------------------------------------------


        ! Body of subroutine Solitary_wave_maker_plate_motion

        !******************************************************************************************************


        !HEIGHT=0.45
        !H0=1.0
        !dt=0.0001

        !write(*,*) TTT,HEIGHT,H0,DT,G

        EZ=0.002
        CL=LOG((1./(EZ**0.5))*(1+SQRT(1.-EZ)))
    !   CL=ARCCOSH(1./EZ**0.5)
        CK=SQRT(3.*HEIGHT/(4.*H0))
          !CK=SQRT(4.*HEIGHT/(3.*H0))
        !CK1=SQRT(3.*HEIGHT/(4.*H0**3))
        C=SQRT(g*(H0+HEIGHT))

        CLAMDA=CL/CK*H0
        CERR=0.0001
          XPI=0.0
        XP=WMX+WMV*DT
          !T0=2*(3.8+HEIGHT/H0)/(C*CK1)
        !SS=SQRT(16*HEIGHT/(3*H0))

    100   XT=CK/H0*(C*TTT-XP-CLAMDA)
          !XPI=CK*H0*(1+TANH(CK1)*(C*(TTT-T0/2)-(XPI-SS/2)))
          XPI=HEIGHT/CK*(TANH(XT)+TANH(CK*CLAMDA/H0))
        ERR=ABS(XPI-XP)
        XP=XPI
        IF(ERR.GT.CERR)GOTO 100


          UP=HEIGHT*C/(H0*(COSH(XT)**2.+HEIGHT/H0))
        VP=2*HEIGHT*C*COSH(XT)*SINH(XT)*CK*C
        !VP=-SQRT(3.)*HEIGHT**1.5*(HEIGHT+H0)*COSH(XT)**3.*SINH(XT)
        VP=VP/(H0*(COSH(XT)**2.+HEIGHT/H0))**2

        WMA=VP  !ACCELERATE
        WMV=UP  !VELOCITY
        WMX=XP  !DISPLACEMENT

        RETURN 
        !******************************************************************************************************


    end subroutine Solitary_wave_maker_plate_motion

    !==========================================================================================================







    !==========================================================================================================

    !**********************************************************************************************************
    !
    !  SUBROUTINE : Wavemaker_Subroutine
    !
    !  PURPOSE    : Assign wavemaker particles information
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
    !**********************************************************************************************************

    subroutine Wavemaker_Subroutine(i_time_step)

        implicit none

        !======================================================================================================
        ! Variables from the superior subroutine
        integer,intent(in)::i_time_step                                            ! Current time step

        !------------------------------------------------------------------------------------------------------
        ! Wave Variables
        integer::i,j,k,l,m,o,v,r                                                   ! Loop variables
        real(kind=8)::along_cline,vertical_cline,angle
        real(kind=8)::damping_coefficient                                          ! Damping coefficient
        real(kind=8)::current_time                                                 ! Current time
        real(kind=8)::body_g                                                       ! Gravity acceleration

        real(kind=8)::current_wave_angle                                           ! Current wave angle
        real(kind=8)::temp,temp_1,temp_2,temp_3,temp_4                             ! Temp value
        real(kind=8),dimension(dim)::temp_analytical_particle_position             ! Current particle position in the loop
        real(kind=8)::current_wave_elevation                                       ! Wave elevation at current particle z direction
        real(kind=8)::m_1                                                          ! m value
        real(kind=8)::piston_stroke                                                ! Piston stroke
        real(kind=8),dimension(dim)::current_particle_position                     ! Current particle position in the loop

        !------------------------------------------------------------------------------------------------------
        
        !======================================================================================================
        

        ! Body of subroutine Wavemaker_Subroutine

        !------------------------------------------------------------------------------------------------------
        ! Define th wavemaker ramping coefficient
        if ( Real_time <= Wavemaker_ramping_time+Standing_time ) then
            Wavemaker_Ramping_coefficient             = 0.5*( sin( ( (Real_time-Standing_time)/Wavemaker_ramping_time-0.5 )*PI )+1 )
            Wavemaker_Ramping_accleration_coefficient = 0.5*  cos( ( (Real_time-Standing_time)/Wavemaker_ramping_time-0.5 )*PI )*PI/Wavemaker_ramping_time
        else
            Wavemaker_Ramping_coefficient             = 1.0d0
            Wavemaker_Ramping_accleration_coefficient = 0.0d0  
        endif
        !------------------------------------------------------------------------------------------------------



        !********************************Assign the wave maker particles information***************************
        ! Piston Type wavemaker for solitary waves or others
        WaveMakerType: if(trim(adjustl(wave_maker_type))=='PistonType') then

            body_g=g
            wave_maker_acceleration = 0.0d0
            wave_maker_velocity     = 0.0d0
            wave_maker_position     = 0.0d0
            
            if ( trim(adjustl(wave_type))=='Solitarywave' ) then

                !----------------------------------------------------------------------------------------------
                ! The solitary wave wavemaker position
                if ( Real_time > Standing_time ) then

                    current_time = Real_time-Standing_time
                    call Solitary_wave_maker_plate_motion(wave_maker_acceleration,wave_maker_velocity,wave_maker_position,current_time,wave_height,water_depth,dt,body_g)
                
                endif
                !----------------------------------------------------------------------------------------------


            elseif ( trim(adjustl(wave_type))=='Strokes_2nd' ) then

                !----------------------------------------------------------------------------------------------
                ! The Strokes 2nd wave wavemaker position
                ! Reference: Long-crested wave generation and absorption for SPH-based DualSPHysics model
                if ( Real_time>Standing_time ) then
                    
                    current_time = Real_time-Standing_time

                    m_1=2*(sinh(wave_number*water_depth))**2/(sinh(wave_number*water_depth)*cosh(wave_number*water_depth)+wave_number*water_depth)

                    current_wave_angle=wave_frequency*current_time+Wave_initial_phase_angle

                    piston_stroke=wave_height/m_1

                    wave_maker_position     = piston_stroke/2.0*sin(current_wave_angle)+((wave_height**2/(32*water_depth))*(3*cosh(water_depth*wave_number)/(sinh(water_depth*wave_number))**3)-2.0/m_1)*sin(2*current_wave_angle)
                    wave_maker_velocity     = piston_stroke/2.0*wave_frequency*cos(current_wave_angle)+((wave_height**2/(32*water_depth))*(3*cosh(water_depth*wave_number)/(sinh(water_depth*wave_number))**3)-2.0/m_1)*2*wave_frequency*cos(2*current_wave_angle)
                    wave_maker_acceleration =-piston_stroke/2.0*wave_frequency**2*sin(current_wave_angle)-((wave_height**2/(32*water_depth))*(3*cosh(water_depth*wave_number)/(sinh(water_depth*wave_number))**3)-2.0/m_1)*4*wave_frequency**2*sin(2*current_wave_angle)

                    wave_maker_position     = Wavemaker_Ramping_coefficient*wave_maker_position
                    wave_maker_velocity     = Wavemaker_Ramping_coefficient*wave_maker_velocity
                    wave_maker_acceleration = Wavemaker_Ramping_coefficient*wave_maker_acceleration

                endif
                !----------------------------------------------------------------------------------------------

            elseif (trim(adjustl(wave_type))=='Airy') then

                !----------------------------------------------------------------------------------------------
                ! The linear wave wavemaker position
                ! Reference: Long-crested wave generation and absorption for SPH-based DualSPHysics model
                if ( Real_time>Standing_time ) then

                    current_time = Real_time-Standing_time

                    m_1=2*(sinh(wave_number*water_depth))**2/(sinh(wave_number*water_depth)*cosh(wave_number*water_depth)+wave_number*water_depth)

                    current_wave_angle=wave_frequency*current_time+Wave_initial_phase_angle

                    piston_stroke=wave_height/m_1

                    wave_maker_position     =  piston_stroke/2.0*sin(current_wave_angle)
                    wave_maker_velocity     =  piston_stroke/2.0*wave_frequency*cos(current_wave_angle)
                    wave_maker_acceleration = -piston_stroke/2.0*wave_frequency**2*sin(current_wave_angle)

                    wave_maker_position     = Wavemaker_Ramping_coefficient*wave_maker_position
                    wave_maker_velocity     = Wavemaker_Ramping_coefficient*wave_maker_velocity
                    wave_maker_acceleration = Wavemaker_Ramping_coefficient*wave_maker_acceleration

                endif
                !----------------------------------------------------------------------------------------------

            endif

            !--------------------------------------------------------------------------------------------------
            ! Check the wave height at the tank ends
            left_wave_height =0.0d0
            right_wave_height=0.0d0

            do k=1,particle_ture_number
                                        
                ! Calculate the distance between the particle and wavemaker plate 
                distance=abs(wave_maker_position-particle_position(k,1))
                
                ! Get the wave height at the left end
                if(distance<=kernel_scale*interior_dx .and. particle_position(k,2)>left_wave_height ) then
                    left_wave_height=particle_position(k,2)
                end if

                ! Calculate the distance between the particle and tank right end
                distance=abs(particle_position(k,2)-tan(theta_rad)*(particle_position(k,1)-incline_x))/sqrt(1+tan(theta_rad)**2)
                
                ! Get the wave height at the right end
                if(distance<=kernel_scale*interior_dx .and. particle_position(k,2)>right_wave_height ) then
                    right_wave_height=particle_position(k,2)
                end if

            end do

            left_wave_height =left_wave_height+0.5*interior_dx
            right_wave_height=right_wave_height+0.5*interior_dx
            !--------------------------------------------------------------------------------------------------
                
            !write(*,*) left_wave_height
            
            !--------------------------------------------------------------------------------------------------
            do j=1,wave_maker_particle_number
               
                i=j+particle_ture_number
                particle_position(i,1)=wave_maker_position-(wave_maker_particle_layer(i)-1)*interior_dx-interior_dx/2.0                      !虚粒子编号，横坐标
                particle_velocity(i,1)=wave_maker_velocity                                                                                                                  !虚粒子编号，x方向速度(初始时刻的位移、速度为0)
                    
                if(particle_position(i,2)>left_wave_height) then
                    
                    particle_rho(i)  =water_rho_0                                                                                              !密度赋值不能出现0
                    particle_press(i)=0.0d0
                    particle_mass(i) =0.0d0                                                       
                
                else

                    particle_rho(i)  =water_rho_0*(-6.0*g*(particle_position(i,2)-left_wave_height)/square_c_0+1.0)**(1.0/6.0)
                    particle_press(i)=square_c_0*water_rho_0*((particle_rho(i)/water_rho_0)**7.0-1.0)/7.0
                    particle_mass(i) =particle_rho(i)*intial_volume_0                             
                    
                end if

            end do
            
            !write(*,*) wave_maker_position,wave_maker_velocity
            !--------------------------------------------------------------------------------------------------
            
            !--------------------------------------------------------------------------------------------------
            ! Drag the particles out the wavemake plate back to the tank
            do i=1,particle_ture_number
                
                if(particle_type(i)==Ill_Particle_Label) cycle        ! No ill particles
                
                if(particle_position(i,1)<wave_maker_position) then
                    particle_position(i,1)=2*wave_maker_position-particle_position(i,1)
                end if
                    
                ! assign the particle velocity
                if((particle_position(i,1)-wave_maker_position)<interior_dx) then
                    particle_velocity(i,1)=wave_maker_velocity
                end if

                ! if(particle_position(i,1)>9.5) then
                   
                !     !计算点到直线的距离
                !     distance=abs(particle_position(i,2)-tan(theta_rad)*(particle_position(i,1)-10))/sqrt(1+tan(theta_rad)**2)    !点到右部倾斜板的距离
               
                !     !对于靠近边界粒子只保留切向速度
                !     if(distance<=0.2*smooth_length ) then
                !         along_cline=particle_velocity(i,1)*cos(theta_rad)+particle_velocity(i,2)*sin(theta_rad)
                !         vertical_cline=-particle_velocity(i,1)*sin(theta_rad)+particle_velocity(i,2)*cos(theta_rad)
            
                !         particle_velocity(i,1)=along_cline*cos(theta_rad)
                !         particle_velocity(i,2)=along_cline*sin(theta_rad)
            
                !     end if
                      
                ! end if
                
                ! !检测底部是否有漏
                ! if(particle_position(i,2)<=0.5*smooth_length ) then         !底部边界
                !    particle_velocity(i,2)=0.5*particle_velocity(i,2)
                ! end if
               
                ! !检测是否越界
                ! if(particle_position(i,2)<0) then         !底部边界拉回底部边界
                !     if(particle_position(i,2)>=-1.5*smooth_length) then
                !         particle_position(i,2)=-particle_position(i,2)
                !         particle_press(i)=particle_press(i)-2*9.81*1000*abs(particle_position(i,2))
                !     else
                !         particle_type(i)=Ill_Particle_Label     
                !     end if   
                    
                ! end if
                
                ! if(particle_position(i,1)>10 .and. particle_position(i,2)>0) then         !底部边界
                !     angle=atan2(particle_position(i,2),(particle_position(i,1)-10))
                !     if(angle<=theta_rad) then
                !         particle_type(i)=Ill_Particle_Label
                !     end if
                ! end if
         
            end do
            !--------------------------------------------------------------------------------------------------

        !======================================================================================================
        ! Assign wavemaker particles directly basedon wave theory
        elseif(trim(adjustl(wave_maker_type))=='WaveTheory') then


            !--------------------------------------------------------------------------------------------------
            ! Check the fluid particles in the domain or not
            do i=1,particle_ture_number
                
                if(particle_type(i)==Ill_Particle_Label) cycle        ! No ill particles
                
                if(particle_position(i,2)<0.0d0) then
                    particle_type(i)=Ill_Particle_Label
                end if
                    
            end do
            !--------------------------------------------------------------------------------------------------


            !--------------------------------------------------------------------------------------------------
            current_time = Real_time-Standing_time


            ! Define the right_wave_height for regular wave
            right_wave_height=water_depth

            !--------------------------------------------------------------------------------------------------
            ! For the time start wave making
            Skip_standing_time: if (current_time>=0.0d0) then

                Wavemaker_zone: if(trim(adjustl(wave_type))=='Stokes_2nd') then

                    !------------------------------------------------------------------------------------------
                    !Wave_elevation=H/2*cos(kx-wt+Fai0)+(PI*H)/8*(H/L)*cosh(kh)[cos(2kh)+2]*cos(2(kx-wt+Fai0))/(sinh(kh))**3
                    !Fai0: wave initial phase angle;
                    !H: wave height;
                    !k: wave number;
                    !w: wave frequency;
                    !t: time;
                    !L: wave length;
                    !h: water depth.

                    ! Define the wavemaker particles velocity 
                    do j=1,wave_maker_particle_number
                   
                        i=j+particle_ture_number

                        current_wave_angle=wave_number*particle_position(i,1)-wave_frequency*current_time+Wave_initial_phase_angle
                        
                        current_wave_elevation=Wavemaker_Ramping_coefficient*(0.5*wave_height*cos(current_wave_angle)+PI*wave_height*wave_height/(8*wave_length)*cosh(wave_number*water_depth)/(sinh(wave_number*water_depth))**3*(2+cosh(2*wave_number*water_depth))*cos(2*current_wave_angle))
                 
                        ! Define the particle base on wave theory
                        particle_velocity(i,1)=PI*wave_height/wave_period*cosh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)+3.0/4.0*(PI**2)*wave_height/wave_period*(wave_height/wave_length)*cosh(2*wave_number*particle_position(i,dim))*cos(2*current_wave_angle)/(sinh(wave_number*water_depth))**4                  
                        particle_velocity(i,dim)=PI*wave_height/wave_period*sinh(wave_number*particle_position(i,dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)+3.0/4.0*(PI**2)*wave_height/wave_period*(wave_height/wave_length)*sinh(2*wave_number*particle_position(i,dim))*sin(2*current_wave_angle)/(sinh(wave_number*water_depth))**4
                        
                        !--------------------------------------------------------------------------------------
                        ! Define the particle base on reamping coefficient
                        particle_velocity(i,1)=Wavemaker_Ramping_coefficient*particle_velocity(i,1)
                        particle_velocity(i,dim)=Wavemaker_Ramping_coefficient*particle_velocity(i,dim)

                        ! Update the wave maker particles position
                        particle_position(i,1)=Initial_particle_Position(i,1)-Wavemaker_Ramping_coefficient*(Wave_amplitue*cosh(wave_number*particle_position(i,dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)+3.0*PI*wave_height**2/(16*wave_length)*cosh(2*wave_number*particle_position(i,dim))*sin(2*current_wave_angle)/(sinh(wave_number*water_depth))**4)
                        particle_position(i,dim)=Initial_particle_Position(i,dim)+Wavemaker_Ramping_coefficient*(Wave_amplitue*sinh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)+3.0*PI*wave_height**2/(16*wave_length)*sinh(2*wave_number*particle_position(i,dim))*cos(2*current_wave_angle)/(sinh(wave_number*water_depth))**4)

                        
                        if (particle_position(i,dim)<=current_wave_elevation+water_depth) then

                            ! Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                            particle_press(i)=water_rho_0*g*(current_wave_elevation+water_depth-particle_position(i,dim))+Wavemaker_Ramping_coefficient*(water_rho_0*g*wave_height/2.0*cosh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/cosh(wave_number*water_depth)+3.0/4.0*water_rho_0*g*wave_height*(PI*wave_height/wave_length)*(1.0/sinh(2*wave_number*water_depth))*(cosh(2*wave_number*particle_position(i,dim))/(sinh(wave_number*water_depth))**2-1.0/3.0)*cos(2*current_wave_angle)-1.0/4.0*water_rho_0*g*wave_height*(PI*wave_height/wave_length)*cosh(2*wave_number*particle_position(i,dim))/sinh(2*wave_number*water_depth))
                            particle_rho(i)=water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                            particle_mass(i)=water_particle_mass_0
                            particle_c(i)=sqrt(square_c_0)

                        else
                           
                            particle_press(i)=0.0d0
                            particle_rho(i)=water_rho_0
                            particle_mass(i)=0.0d0
                            particle_c(i)=sqrt(square_c_0)                        
                            
                        endif
                        !--------------------------------------------------------------------------------------

                    end do
                    !------------------------------------------------------------------------------------------

                elseif(trim(adjustl(wave_type))=='Airy') then
                    
                    !------------------------------------------------------------------------------------------
                    ! Wave_elevation=H/2*cos(kx-wt+Fai0)
                    
                    ! Define the wavemaker particles velocity 
                    do j=1,wave_maker_particle_number
                   
                        i=j+particle_ture_number

                        current_wave_angle=wave_number*particle_position(i,1)-wave_frequency*current_time+Wave_initial_phase_angle
                        current_wave_elevation=Wavemaker_Ramping_coefficient*0.5*wave_height*cos(current_wave_angle)

                        ! Define the particle base on wave theory
                        particle_velocity(i,1)=PI*wave_height/wave_period*cosh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)              
                        particle_velocity(i,dim)=PI*wave_height/wave_period*sinh(wave_number*particle_position(i,dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)
            
                        !--------------------------------------------------------------------------------------
                        ! Define the particle base on soft_coefficient
                        particle_velocity(i,1)=Wavemaker_Ramping_coefficient*particle_velocity(i,1)
                        particle_velocity(i,dim)=Wavemaker_Ramping_coefficient*particle_velocity(i,dim)

                        ! Update the wave maker particles position
                        particle_position(i,1)=Initial_particle_Position(i,1)-Wavemaker_Ramping_coefficient*Wave_amplitue*cosh(wave_number*particle_position(i,dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)
                        particle_position(i,dim)=Initial_particle_Position(i,dim)+Wavemaker_Ramping_coefficient*Wave_amplitue*sinh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)

                        
                        if (particle_position(i,dim)<=current_wave_elevation+water_depth) then

                            ! Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                            particle_press(i)= water_rho_0*g*(water_depth+current_wave_elevation-particle_position(i,dim))+Wavemaker_Ramping_coefficient*(water_rho_0*g*wave_height/2.0*cosh(wave_number*particle_position(i,dim))*cos(current_wave_angle)/cosh(wave_number*water_depth))
                            particle_rho(i)  = water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                            particle_mass(i) = water_particle_mass_0
                            particle_c(i)    = sqrt(square_c_0)

                        else
                           
                            particle_press(i)=0.0d0
                            particle_rho(i)=water_rho_0
                            particle_mass(i)=0.0d0
                            particle_c(i)=sqrt(square_c_0)                        
                            
                        endif
                        !--------------------------------------------------------------------------------------

                    end do
                    !------------------------------------------------------------------------------------------

                elseif(trim(adjustl(wave_type))=='Stokes_5th') then
                    
                    !Define the wavemaker particles velocity 
                    do j=1,wave_maker_particle_number
                   
                        i=j+particle_ture_number

                        current_wave_angle=wave_number*particle_position(i,1)-wave_frequency*current_time+Wave_initial_phase_angle

                        !--------------------------------------------------------------------------------------
                        ! Coefficients (Deep water wave is the initial condition)
                        ! "Hydrodynamics of Offshore Structures.P59"
                        current_wave_elevation=0.0d0
                        do m=1,5
                            current_wave_elevation=current_wave_elevation+Fifth_stokes_E(m)*cos(m*current_wave_angle)
                        enddo
                        current_wave_elevation=Wavemaker_Ramping_coefficient*current_wave_elevation/wave_number

                        !Define the particle base on wave theory
                        particle_velocity(i,:)=0.0d0
                        temp_analytical_particle_position(:)=0.0d0
                        do m=1,5

                            particle_velocity(i,1)=particle_velocity(i,1)+m*Fifth_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)
                            particle_velocity(i,dim)=particle_velocity(i,dim)+m*Fifth_stokes_D(m)*sinh(m*wave_number*particle_position(i,dim))*sin(m*current_wave_angle)

                            temp_analytical_particle_position(1)=temp_analytical_particle_position(1)+Fifth_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*sin(m*current_wave_angle)
                            temp_analytical_particle_position(dim)=temp_analytical_particle_position(dim)+Fifth_stokes_D(m)*sinh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)
                        
                        enddo

                        particle_velocity(i,1)=Wavemaker_Ramping_coefficient*Wave_celerity*particle_velocity(i,1)
                        particle_velocity(i,dim)=Wavemaker_Ramping_coefficient*Wave_celerity*particle_velocity(i,dim)
                        
                        temp_analytical_particle_position(1)=Wavemaker_Ramping_coefficient*(-Wave_celerity/wave_frequency)*temp_analytical_particle_position(1)
                        temp_analytical_particle_position(dim)=Wavemaker_Ramping_coefficient*(Wave_celerity/wave_frequency)*temp_analytical_particle_position(dim)

                        !Update the wave maker particles position
                        particle_position(i,1)=Initial_particle_Position(i,1)+temp_analytical_particle_position(1)
                        particle_position(i,dim)=Initial_particle_Position(i,dim)+temp_analytical_particle_position(dim)

                        ! Stokes_DFai_Dt=0.0
                        ! do m=1,5

                        !     Stokes_DFai_Dt=Stokes_DFai_Dt+m*Fifth_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)

                        ! enddo
                        ! Stokes_DFai_Dt=-Wave_celerity**2*Stokes_DFai_Dt*Wavemaker_Ramping_coefficient

                        if (particle_position(i,dim)<=current_wave_elevation+water_depth) then

                            ! Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                            !particle_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-particle_position(i,dim))-water_rho_0*0.5*(particle_velocity(i,1)**2+particle_velocity(i,dim)**2)+water_rho_0*Stokes_DFai_Dt
                            
                            particle_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-particle_position(i,dim))-water_rho_0*0.5*((particle_velocity(i,1)-Wavemaker_Ramping_coefficient*Wave_celerity)**2+particle_velocity(i,dim)**2)-water_rho_0*Wavemaker_Ramping_coefficient*Bernoulli_constant_R
                            
                            particle_rho(i)=water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                            particle_mass(i)=water_particle_mass_0
                            particle_c(i)=sqrt(square_c_0)

                        else
                           
                            particle_press(i)=0.0d0
                            particle_rho(i)=water_rho_0
                            particle_mass(i)=0.0d0
                            particle_c(i)=sqrt(square_c_0)                        
                            
                        endif
                        !--------------------------------------------------------------------------------------

                    end do
                    !------------------------------------------------------------------------------------------

                elseif(trim(adjustl(wave_type))=='Stokes_3rd') then
                    
                    ! Define the wavemaker particles velocity 
                    do j=1,wave_maker_particle_number
                   
                        i=j+particle_ture_number

                        current_wave_angle=wave_number*particle_position(i,1)-wave_frequency*current_time+Wave_initial_phase_angle

                        !--------------------------------------------------------------------------------------
                        ! Reference: Song Zhiyao "On the universal third order stokes wave solution"
                        current_wave_elevation=0.0d0

                        do m=1,3
                            current_wave_elevation=current_wave_elevation+Third_stokes_E(m)*cos(m*current_wave_angle)
                        enddo
                        current_wave_elevation=Wavemaker_Ramping_coefficient*current_wave_elevation/wave_number

                        ! Define the particle base on wave theory
                        particle_velocity(i,:)=0.0d0
                        temp_analytical_particle_position(:)=0.0d0
                        do m=1,3

                            particle_velocity(i,1)=particle_velocity(i,1)+m*Third_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)
                            particle_velocity(i,dim)=particle_velocity(i,dim)+m*Third_stokes_D(m)*sinh(m*wave_number*particle_position(i,dim))*sin(m*current_wave_angle)

                            temp_analytical_particle_position(1)=temp_analytical_particle_position(1)+Third_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*sin(m*current_wave_angle)
                            temp_analytical_particle_position(dim)=temp_analytical_particle_position(dim)+Third_stokes_D(m)*sinh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)
                        
                        enddo

                        particle_velocity(i,1)=Wavemaker_Ramping_coefficient*Wave_celerity*particle_velocity(i,1)
                        particle_velocity(i,dim)=Wavemaker_Ramping_coefficient*Wave_celerity*particle_velocity(i,dim)
                        
                        temp_analytical_particle_position(1)=Wavemaker_Ramping_coefficient*(-Wave_celerity/wave_frequency)*temp_analytical_particle_position(1)
                        temp_analytical_particle_position(dim)=Wavemaker_Ramping_coefficient*(Wave_celerity/wave_frequency)*temp_analytical_particle_position(dim)

                        ! Update the wave maker particles position
                        particle_position(i,1)=Initial_particle_Position(i,1)+temp_analytical_particle_position(1)
                        particle_position(i,dim)=Initial_particle_Position(i,dim)+temp_analytical_particle_position(dim)

                        Stokes_DFai_Dt=0.0
                        do m=1,3

                            Stokes_DFai_Dt=Stokes_DFai_Dt+m*Third_stokes_D(m)*cosh(m*wave_number*particle_position(i,dim))*cos(m*current_wave_angle)

                        enddo
                        Stokes_DFai_Dt=-Wave_celerity**2*Stokes_DFai_Dt*Wavemaker_Ramping_coefficient

                        if (particle_position(i,dim)<=current_wave_elevation+water_depth) then

                            ! Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                            particle_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-particle_position(i,dim))-water_rho_0*0.5*(particle_velocity(i,1)**2+particle_velocity(i,dim)**2)+water_rho_0*Stokes_DFai_Dt

                            particle_rho(i)=water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                            particle_mass(i)=water_particle_mass_0
                            particle_c(i)=sqrt(square_c_0)

                        else
                           
                            particle_press(i)=0.0d0
                            particle_rho(i)=water_rho_0
                            particle_mass(i)=0.0d0
                            particle_c(i)=sqrt(square_c_0)                        
                            
                        endif
                        !--------------------------------------------------------------------------------------

                    end do
                    !------------------------------------------------------------------------------------------

                elseif(trim(adjustl(wave_type))=='Solitarywave') then
                    
                    ! Define the wavemaker particles velocity 
                    do j=1,wave_maker_particle_number
                   
                        i=j+particle_ture_number

                        Solitary_s=1.0/cosh(Solitary_alpha*particle_position(i,1))
                        Solitary_t=tanh(Solitary_alpha*particle_position(i,1))
                        
                        current_wave_angle=wave_number*particle_position(i,1)-wave_frequency*current_time+Wave_initial_phase_angle

                        !--------------------------------------------------------------------------------------
                        ! Reference: John Fenton "A ninth-order solution for the solitary wave"
                        current_wave_elevation=1+Solitary_Epsilon*Solitary_s**2-3.0/4.0*Solitary_Epsilon**2*Solitary_s**2*Solitary_t**2+Solitary_Epsilon**3*(5.0/8.0*Solitary_s**2*Solitary_t**2-101/80*Solitary_s**4*Solitary_t**2)


                        ! Define the particle base on wave theory
                        particle_velocity(i,1)=1+0.5*Solitary_Epsilon-3.0/20.0*Solitary_Epsilon**2+3/56.0*Solitary_Epsilon**3-Solitary_Epsilon*Solitary_s**2 &
                                                +Solitary_Epsilon**2*(-0.25*Solitary_s**2+Solitary_s**4+particle_position(i,dim)**2*(1.5*Solitary_s**2-2.25*Solitary_s**4))&
                                                +Solitary_Epsilon**3*(19.0/40.0*Solitary_s**2+0.2*Solitary_s**4-1.2*Solitary_s**6+particle_position(i,dim)**2*(-1.5*Solitary_s**2-3.75*Solitary_s**4+7.5*Solitary_s**6)+particle_position(i,dim)**4*(-0.375*Solitary_s**2+45.0/16.0*Solitary_s**4-45.0/16.0*Solitary_s**6))
                        
                        
                        particle_velocity(i,dim)=-Solitary_Epsilon*Solitary_s**2 &
                                                 +Solitary_Epsilon**2*(0.375*Solitary_s**2+2*Solitary_s**4+particle_position(i,dim)**2*(0.5*Solitary_s**2-1.5*Solitary_s**4)) &
                                                 +Solitary_Epsilon**3*(49/640.0*Solitary_s**2-17.0/20.0*Solitary_s**4-18/5.0*Solitary_s**6+particle_position(i,dim)**2*(-13.0/16.0*Solitary_s**2-25.0/16.0*Solitary_s**4+7.5*Solitary_s**6)+particle_position(i,dim)**4*(-3/40.0*Solitary_s**2+9.0/8.0*Solitary_s**4-27/16.0*Solitary_s**6))

                        
                        particle_velocity(i,1)=sqrt(g*water_depth)*particle_velocity(i,1)
                        particle_velocity(i,dim)=sqrt(3*Solitary_Epsilon)*particle_position(i,dim)*Solitary_t*sqrt(g*water_depth)*particle_velocity(i,dim)

                        if (particle_position(i,dim)<=current_wave_elevation+water_depth) then

                            ! Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                            particle_press(i)= water_rho_0*g*water_depth*(1-particle_position(i,dim)+Solitary_Epsilon*Solitary_s**2+Solitary_Epsilon**2*(0.75*Solitary_s**2-1.5*Solitary_s**4+particle_position(i,dim)**2*(-1.5*Solitary_s**2+2.25*Solitary_s**4))&
                                                                          +Solitary_Epsilon**3*(-0.5*Solitary_s**2-19.0/20.0*Solitary_s**4+11.0/5*Solitary_s**6+particle_position(i,dim)**2*(0.75*Solitary_s**2+39.0/8.0*Solitary_s**4-33.0/4.0*Solitary_s**6)+particle_position(i,dim)**4*(0.375*Solitary_s**2-45.0/16.0*Solitary_s**4+45/16.0*Solitary_s**6))&
                                                                         )

                            particle_rho(i)  = water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                            particle_mass(i) = water_particle_mass_0
                            particle_c(i)    = sqrt(square_c_0)

                        else
                           
                            particle_press(i)=0.0d0
                            particle_rho(i)  =water_rho_0
                            particle_mass(i) =0.0d0
                            particle_c(i)    =sqrt(square_c_0)                        
                            
                        endif

                        ! Update the wave maker particles position
                        particle_position(i,1)=particle_position(i,1)+particle_velocity(i,1)*dt
                        particle_position(i,dim)=particle_position(i,dim)+particle_velocity(i,dim)*dt
                        !--------------------------------------------------------------------------------------

                    end do
                    !------------------------------------------------------------------------------------------

                else
                    
                    ! The main processor output the error information
                    if (Current_Processor_ID==Main_Processor) then

                        write(*,'(A)') " The wave type is not include in the code, please review the type."

                    endif

                endif Wavemaker_zone
                !==============================================================================================


                !==============================================================================================
                ! Update the particles in For the wavemaker damping zone
                ! This is an independent data, not included in the SPH calculation
                WavemakerDamping:if (trim(adjustl(WavemakerDamping_Switch))=='on') then

                    do i=1,particle_ture_number

                        if(particle_type(i)==Ill_Particle_Label) cycle                        ! No error particles

                        if (In_WaveMakerDampingZone_OrNot(i)/=1) cycle                        ! Only for the particles in Wavemaker Damping Zone

                        !--------------------------------------------------------------------------------------
                        ! Wave_elevation=H/2*cos(kx-wt+Fai0)+(PI*H)/8*(H/L)*cosh(kh)[cos(2kh)+2]*cos(2(kx-wt+Fai0))/(sinh(kh))**3
                        if(trim(adjustl(wave_type))=='Stokes_2nd') then

                            current_particle_position(:)=WaveMaker_Analytical_position(i,:)

                            current_wave_angle=wave_number*current_particle_position(1)-wave_frequency*current_time+Wave_initial_phase_angle
                            current_wave_elevation=Wavemaker_Ramping_coefficient*(0.5*wave_height*cos(current_wave_angle)+PI*wave_height*wave_height*(1.0/tanh(wave_number*water_depth))*(1+3.0/(2*(sinh(wave_number*water_depth))**3)*cos(2*current_wave_angle)))

                            ! Define the particle base on wave theory
                            WaveMaker_Analytical_velocity(i,1)=PI*wave_height/wave_period*cosh(wave_number*current_particle_position(dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)-3.0/4.0*(PI**2)*wave_height/wave_period*(wave_height/wave_length)*cosh(2*wave_number*current_particle_position(dim))*cos(2*current_wave_angle)/(sinh(wave_number*water_depth))**4                  
                            WaveMaker_Analytical_velocity(i,dim)=PI*wave_height/wave_period*sinh(wave_number*current_particle_position(dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)-3.0/4.0*(PI**2)*wave_height/wave_period*(wave_height/wave_length)*sinh(2*wave_number*current_particle_position(dim))*sin(2*current_wave_angle)/(sinh(wave_number*water_depth))**4
                            
                            !----------------------------------------------------------------------------------
                            ! Define the particle base on soft_coefficient
                            WaveMaker_Analytical_velocity(i,1)=Wavemaker_Ramping_coefficient*WaveMaker_Analytical_velocity(i,1)
                            WaveMaker_Analytical_velocity(i,dim)=Wavemaker_Ramping_coefficient*WaveMaker_Analytical_velocity(i,dim)

                            ! Update the wave maker particles position  
                            WaveMaker_Analytical_position(i,1)=Initial_particle_Position(i,1)-Wavemaker_Ramping_coefficient*(Wave_amplitue*cosh(wave_number*current_particle_position(dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)+3.0*PI*wave_height**2/(16*wave_length)*cosh(2*wave_number*current_particle_position(dim))*sin(2*current_wave_angle)/(sinh(wave_number*water_depth))**4)
                            WaveMaker_Analytical_position(i,dim)=Initial_particle_Position(i,dim)+Wavemaker_Ramping_coefficient*(Wave_amplitue*sinh(wave_number*current_particle_position(dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)+3.0*PI*wave_height**2/(16*wave_length)*sinh(2*wave_number*current_particle_position(dim))*cos(2*current_wave_angle)/(sinh(wave_number*water_depth))**4)

                            if (WaveMaker_Analytical_position(i,dim)<=current_wave_elevation+water_depth) then

                                ! Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                                WaveMaker_Analytical_press(i)=water_rho_0*g*(current_wave_elevation+water_depth-WaveMaker_Analytical_position(i,dim))+Wavemaker_Ramping_coefficient*(water_rho_0*g*wave_height/2.0*cosh(wave_number*WaveMaker_Analytical_position(i,dim))*cos(current_wave_angle)/cosh(wave_number*water_depth)+3.0/4.0*water_rho_0*g*wave_height*(PI*wave_height/wave_length)*(1.0/sinh(2*wave_number*water_depth))*(cosh(2*wave_number*WaveMaker_Analytical_position(i,dim))/(sinh(wave_number*water_depth))**2-1.0/3.0)*cos(2*current_wave_angle)-1.0/4.0*water_rho_0*g*wave_height*(PI*wave_height/wave_length)*cosh(2*wave_number*WaveMaker_Analytical_position(i,dim))/sinh(2*wave_number*water_depth))
                                WaveMaker_Analytical_rho(i)=water_rho_0*(1+7*WaveMaker_Analytical_press(i)/(square_c_0*water_rho_0))**(0.142857142)

                            else
                               
                                WaveMaker_Analytical_press(i)=0.0d0
                                WaveMaker_Analytical_rho(i)=water_rho_0                      
                                
                            endif
                            !----------------------------------------------------------------------------------


                        !--------------------------------------------------------------------------------------
                        ! Wave_elevation=H/2*cos(kx-wt+Fai0)
                        elseif(trim(adjustl(wave_type))=='Airy') then

                            !----------------------------------------------------------------------------------
                            current_particle_position(:)=WaveMaker_Analytical_position(i,:)

                            current_wave_angle=wave_number*current_particle_position(1)-wave_frequency*current_time+Wave_initial_phase_angle
                            current_wave_elevation=Wavemaker_Ramping_coefficient*0.5*wave_height*cos(current_wave_angle)

                            ! Define the particle base on wave theory
                            WaveMaker_Analytical_velocity(i,1)=PI*wave_height/wave_period*cosh(wave_number*current_particle_position(dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)              
                            WaveMaker_Analytical_velocity(i,dim)=PI*wave_height/wave_period*sinh(wave_number*current_particle_position(dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)
                
                            !----------------------------------------------------------------------------------
                            ! Define the particle base on soft_coefficient
                            WaveMaker_Analytical_velocity(i,1)=Wavemaker_Ramping_coefficient*WaveMaker_Analytical_velocity(i,1)
                            WaveMaker_Analytical_velocity(i,dim)=Wavemaker_Ramping_coefficient*WaveMaker_Analytical_velocity(i,dim)

                            ! Update the wave maker particles position
                            !position=initial_position+the ellpside
                            
                            WaveMaker_Analytical_position(i,1)=Initial_particle_Position(i,1)-Wavemaker_Ramping_coefficient*Wave_amplitue*cosh(wave_number*current_particle_position(dim))*sin(current_wave_angle)/sinh(wave_number*water_depth)              
                            WaveMaker_Analytical_position(i,dim)=Initial_particle_Position(i,dim)+Wavemaker_Ramping_coefficient*Wave_amplitue*sinh(wave_number*current_particle_position(dim))*cos(current_wave_angle)/sinh(wave_number*water_depth)              
                           

                            if (WaveMaker_Analytical_position(i,dim)<=current_wave_elevation+water_depth) then

                                ! Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                                WaveMaker_Analytical_press(i)=water_rho_0*g*(current_wave_elevation+water_depth-WaveMaker_Analytical_position(i,dim))+Wavemaker_Ramping_coefficient*(water_rho_0*g*wave_height/2.0*cosh(wave_number*WaveMaker_Analytical_position(i,dim))*cos(current_wave_angle)/cosh(wave_number*water_depth))
                                WaveMaker_Analytical_rho(i)=water_rho_0*(1+7*WaveMaker_Analytical_press(i)/(square_c_0*water_rho_0))**(0.142857142)

                            else
                               
                                WaveMaker_Analytical_press(i)=0.0d0
                                WaveMaker_Analytical_rho(i)=water_rho_0                        
                                
                            endif
                            !----------------------------------------------------------------------------------


                        elseif(trim(adjustl(wave_type))=='Stokes_5th') then

                            !----------------------------------------------------------------------------------
                            current_particle_position(:)=WaveMaker_Analytical_position(i,:)

                            current_wave_angle=wave_number*current_particle_position(1)-wave_frequency*current_time+Wave_initial_phase_angle
                            
                            !----------------------------------------------------------------------------------
                            ! Coefficients (Deep water wave is the initial condition)
                            ! "Hydrodynamics of Offshore Structures.P59"
                            current_wave_elevation=0.0d0
                            do m=1,5
                                current_wave_elevation=current_wave_elevation+Fifth_stokes_E(m)*cos(m*current_wave_angle)
                            enddo
                            current_wave_elevation=Wavemaker_Ramping_coefficient*current_wave_elevation/wave_number


                            ! Define the particle base on wave theory
                            WaveMaker_Analytical_velocity(i,:)=0.0d0
                            temp_analytical_particle_position(:)=0.0d0
                            do m=1,5

                                WaveMaker_Analytical_velocity(i,1)=WaveMaker_Analytical_velocity(i,1)+m*Fifth_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)
                                WaveMaker_Analytical_velocity(i,dim)=WaveMaker_Analytical_velocity(i,dim)+m*Fifth_stokes_D(m)*sinh(m*wave_number*current_particle_position(dim))*sin(m*current_wave_angle)

                                temp_analytical_particle_position(1)=temp_analytical_particle_position(1)+Fifth_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*sin(m*current_wave_angle)
                                temp_analytical_particle_position(dim)=temp_analytical_particle_position(dim)+Fifth_stokes_D(m)*sinh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)
                            
                            enddo

                            WaveMaker_Analytical_velocity(i,1)=Wavemaker_Ramping_coefficient*Wave_celerity*WaveMaker_Analytical_velocity(i,1)
                            WaveMaker_Analytical_velocity(i,dim)=Wavemaker_Ramping_coefficient*Wave_celerity*WaveMaker_Analytical_velocity(i,dim)
                            
                            temp_analytical_particle_position(1)=Wavemaker_Ramping_coefficient*(-Wave_celerity/wave_frequency)*temp_analytical_particle_position(1)
                            temp_analytical_particle_position(dim)=Wavemaker_Ramping_coefficient*(Wave_celerity/wave_frequency)*temp_analytical_particle_position(dim)

                            ! Update the wave maker particles position
                            WaveMaker_Analytical_position(i,1)=Initial_particle_Position(i,1)+temp_analytical_particle_position(1)
                            WaveMaker_Analytical_position(i,dim)=Initial_particle_Position(i,dim)+temp_analytical_particle_position(dim)


                            ! Stokes_DFai_Dt=0.0
                            ! do m=1,5

                            !     Stokes_DFai_Dt=Stokes_DFai_Dt+m*Fifth_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)

                            ! enddo
                            ! Stokes_DFai_Dt=-Wave_celerity**2*Stokes_DFai_Dt*Wavemaker_Ramping_coefficient

                            if (WaveMaker_Analytical_position(i,dim)<=current_wave_elevation+water_depth) then

                                !Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                                !WaveMaker_Analytical_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-current_particle_position(dim))-water_rho_0*0.5*(WaveMaker_Analytical_velocity(i,1)**2+WaveMaker_Analytical_velocity(i,dim)**2)+water_rho_0*Stokes_DFai_Dt
                                
                                WaveMaker_Analytical_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-WaveMaker_Analytical_position(i,dim))-water_rho_0*0.5*((WaveMaker_Analytical_velocity(i,1)-Wavemaker_Ramping_coefficient*Wave_celerity)**2+WaveMaker_Analytical_velocity(i,dim)**2)-Wavemaker_Ramping_coefficient*water_rho_0*Bernoulli_constant_R
                                
                                WaveMaker_Analytical_rho(i)=water_rho_0*(1+7*WaveMaker_Analytical_press(i)/(square_c_0*water_rho_0))**(0.142857142)

                            else
                               
                                WaveMaker_Analytical_press(i)=0.0d0
                                WaveMaker_Analytical_rho(i)=water_rho_0                        
                                
                            endif
                            !----------------------------------------------------------------------------------

                        elseif(trim(adjustl(wave_type))=='Stokes_3rd') then

                            !----------------------------------------------------------------------------------

                            current_particle_position(:)=WaveMaker_Analytical_position(i,:)

                            current_wave_angle=wave_number*current_particle_position(1)-wave_frequency*current_time+Wave_initial_phase_angle
                            
                            !----------------------------------------------------------------------------------
                            ! Coefficients (Deep water wave is the initial condition)
                            ! "Hydrodynamics of Offshore Structures.P59"
                            current_wave_elevation=0.0d0
                            do m=1,3
                                current_wave_elevation=current_wave_elevation+Third_stokes_E(m)*cos(m*current_wave_angle)
                            enddo
                            current_wave_elevation=Wavemaker_Ramping_coefficient*current_wave_elevation/wave_number


                            ! Define the particle base on wave theory
                            WaveMaker_Analytical_velocity(i,:)=0.0d0
                            temp_analytical_particle_position(:)=0.0d0
                            do m=1,3

                                WaveMaker_Analytical_velocity(i,1)=WaveMaker_Analytical_velocity(i,1)+m*Third_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)
                                WaveMaker_Analytical_velocity(i,dim)=WaveMaker_Analytical_velocity(i,dim)+m*Third_stokes_D(m)*sinh(m*wave_number*current_particle_position(dim))*sin(m*current_wave_angle)

                                temp_analytical_particle_position(1)=temp_analytical_particle_position(1)+Third_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*sin(m*current_wave_angle)
                                temp_analytical_particle_position(dim)=temp_analytical_particle_position(dim)+Third_stokes_D(m)*sinh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)
                            
                            enddo

                            WaveMaker_Analytical_velocity(i,1)=Wavemaker_Ramping_coefficient*Wave_celerity*WaveMaker_Analytical_velocity(i,1)
                            WaveMaker_Analytical_velocity(i,dim)=Wavemaker_Ramping_coefficient*Wave_celerity*WaveMaker_Analytical_velocity(i,dim)
                            
                            temp_analytical_particle_position(1)=Wavemaker_Ramping_coefficient*(-Wave_celerity/wave_frequency)*temp_analytical_particle_position(1)
                            temp_analytical_particle_position(dim)=Wavemaker_Ramping_coefficient*(Wave_celerity/wave_frequency)*temp_analytical_particle_position(dim)

                            ! Update the wave maker particles position
                            WaveMaker_Analytical_position(i,1)=Initial_particle_Position(i,1)+temp_analytical_particle_position(1)
                            WaveMaker_Analytical_position(i,dim)=Initial_particle_Position(i,dim)+temp_analytical_particle_position(dim)

                            Stokes_DFai_Dt=0.0
                            do m=1,3

                                Stokes_DFai_Dt=Stokes_DFai_Dt+m*Third_stokes_D(m)*cosh(m*wave_number*current_particle_position(dim))*cos(m*current_wave_angle)

                            enddo
                            Stokes_DFai_Dt=-Wave_celerity**2*Stokes_DFai_Dt*Wavemaker_Ramping_coefficient

                            if (WaveMaker_Analytical_position(i,dim)<=current_wave_elevation+water_depth) then

                                ! Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                                WaveMaker_Analytical_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-current_particle_position(dim))-water_rho_0*0.5*(WaveMaker_Analytical_velocity(i,1)**2+WaveMaker_Analytical_velocity(i,dim)**2)+water_rho_0*Stokes_DFai_Dt
                                
                                !WaveMaker_Analytical_press(i)=water_rho_0*g*(water_depth+current_wave_elevation-WaveMaker_Analytical_position(i,dim))+water_rho_0*0.5*((WaveMaker_Analytical_velocity(i,1)-Wavemaker_Ramping_coefficient*Wave_celerity)**2+WaveMaker_Analytical_velocity(i,dim)**2)-Wavemaker_Ramping_coefficient*water_rho_0*Bernoulli_constant_R
                                
                                WaveMaker_Analytical_rho(i)=water_rho_0*(1+7*WaveMaker_Analytical_press(i)/(square_c_0*water_rho_0))**(0.142857142)

                            else
                               
                                WaveMaker_Analytical_press(i)=0.0d0
                                WaveMaker_Analytical_rho(i)=water_rho_0                        
                                
                            endif
                            !----------------------------------------------------------------------------------

                        !--------------------------------------------------------------------------------------
                        else
                            ! The main processor output the error information
                            if (Current_Processor_ID==Main_Processor) then

                              write(*,'(A)') " The wave type is not include in the code, please review the type.(wavemaker subroutine) "

                            endif

                        endif
                        !--------------------------------------------------------------------------------------

                    enddo

                endif WavemakerDamping
                !==============================================================================================
            
            else

                !==============================================================================================
                ! For the standing time All the results for the hydrostatic 

                !----------------------------------------------------------------------------------------------
                do j=1,wave_maker_particle_number
               
                    i=j+particle_ture_number

                    ! Define the particle base on wave theory
                    particle_velocity(i,:)=0.0d0

                    ! Attention: water_rho_0*g*0.5*interior_dx to strength the pressure ot prevent the point go outside
                    particle_press(i)=water_rho_0*g*(water_depth-particle_position(i,dim))+water_rho_0*g*0.25*interior_dx  
                    particle_rho(i)=water_rho_0*(1+7*particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)
                    particle_mass(i)=water_particle_mass_0
                    particle_c(i)=sqrt(square_c_0)

                end do
                !----------------------------------------------------------------------------------------------


                !----------------------------------------------------------------------------------------------
                do i=1,particle_ture_number

                    if(particle_type(i)==Ill_Particle_Label) cycle                        ! No error particles

                    if (In_WaveMakerDampingZone_OrNot(i)/=1) cycle                        ! Only for the particles in Wavemaker Damping Zone
                            
                    ! Define the particle base on wave theory
                    WaveMaker_Analytical_velocity(i,:)=0.0d0

                    ! Update the wave maker particles position  
                    WaveMaker_Analytical_position(i,:)=Initial_particle_Position(i,:)

                    WaveMaker_Analytical_press(i)=water_rho_0*g*(water_depth-WaveMaker_Analytical_position(i,dim))
                    WaveMaker_Analytical_rho(i)=water_rho_0*(1+7*WaveMaker_Analytical_press(i)/(square_c_0*water_rho_0))**(0.142857142)

                enddo
                !----------------------------------------------------------------------------------------------

            endif Skip_standing_time
            !==================================================================================================


        !======================================================================================================
        !The main processor output the error information
        else

           if (Current_Processor_ID==Main_Processor) then

              write(*,'(A)') " The wave maker type is not include in the code, please review the type."

           endif

        endif WaveMakerType
        !======================================================================================================





    !     !******************************************************************************************************
    !     ! For Debug
    !     ! The main Processor Output the initial particle information
    !     if (Current_Processor_ID==Main_Processor) then

    !         ! Open the Output file
    !         open(unit=General_File_Port,File="Check_Wavemaker.dat")    

    !         ! input tecplot header
    !         write(General_File_Port,*) "TITLE='DISTRIBUTION'"
    !         write(General_File_Port,*) "VARIABLES= 'X' 'Y' 'VX' 'VY' 'RHO' 'P' "
            
    !         ! output fluid particle
    !         write(General_File_Port,*) "ZONE I=",particle_ture_number," F=POINT"
    !         do i=1,particle_ture_number
    !              write(General_File_Port,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    ! 100          format(8F20.10) 
    !         end do
            
    !         ! output wave maker particle
    !         write(General_File_Port,*) "ZONE I=",wave_maker_particle_number," F=POINT"
    !         do k=1,wave_maker_particle_number
    !             i=k+particle_ture_number
    !             write(General_File_Port,100) (particle_position(i,j),j=1,dim),(particle_velocity(i,j),j=1,dim),particle_rho(i),particle_press(i)/pressure_0
    !         end do
                
            
    !         close(General_File_Port)

    !     endif

    !     ! Synchronize all processors calculation
    !     ! MPI_BCAST function has the character "Synchronize", so we don't need Synchronize
    !     call MPI_Barrier(MPI_COMM_WORLD,ierror_MPI)
    !     !******************************************************************************************************


    end subroutine Wavemaker_Subroutine

    !==========================================================================================================








    !==========================================================================================================

    !**********************************************************************************************************
    !
    !  PROGRAM   : wave_maker_damping_zone
    !
    !  Programer : Shanqin Jin (E-mail:sjin@mun.ca)
    !
    !  Location  : MUN
    !
    !  Time      : 2017.3.18
    !
    !  Copyright : Memorial University
    !
    !  Version   : 1.0
    !
    !  Note      : MPI version: mpich-3.2
    !            Fortran version: Inter Fortran (ifort) 18.0.0
    !**********************************************************************************************************

    subroutine Wave_Maker_Damping_Zone(i_time_step)              

        
        implicit none

        !======================================================================================================
        ! Variables form the superior subroutine
        integer,intent(in)::i_time_step                               ! Current time step

        ! Variables in subroutine
        integer::i,j,k,L                                              ! Variables for looping
        real(kind=8)::damping_coefficient                             ! Damping coefficient
        real(kind=8)::current_time                                    ! Current time

        real(kind=8)::current_wave_angle                              ! Current wave angle
        real(kind=8)::current_wave_elevation                          ! Wave elevation at current particle z direction

        ! Variables for file opeartion
        integer::File_index
        character(len=100)::File_name
        character(len=4)::Char_Current_Processor_ID

        !======================================================================================================

        ! Body of subroutine Wave_Maker_Damping_Zone


        !------------------------------------------------------------------------------------------------------
        ! For the wavemaker damping zone
        do i=1,actual_particle_number_in_subdomain

           !if(subdomain_particle_type(i)/=2 .or. subdomain_In_WaveMakerDampingZone_OrNot(i)/=1) cycle      ! Only for the fluid particles

           !! Define the particles in the wavemaker damping zone or not
           if(subdomain_In_WaveMakerDampingZone_OrNot(i)==1) then   

            ! if(subdomain_particle_type(i)/=2 ) cycle
            ! if(subdomain_particle_position(i,1)<=WavemakerDampingLengh) then

                ! write(*,*) subdomain_In_WaveMakerDampingZone_OrNot(i),Current_Processor_ID 

                ! calculate the position damping coefficient
                
                ! if (Current_Processor_ID==1) then
                !    write(*,*) subdomain_particle_position(i,1),WavemakerDampingLengh,subdomain_position_soft_coefficient(i)
                ! endif
                   
                ! !Renew the wavemaker damping zone 
                ! subdomain_particle_press(i)=position_soft_coefficient*subdomain_WaveMaker_Analytical_press(i)+(1-position_soft_coefficient)*subdomain_particle_press(i)   
                ! subdomain_particle_rho(i)=position_soft_coefficient*subdomain_WaveMaker_Analytical_rho(i)+(1-position_soft_coefficient)*subdomain_particle_rho(i)
                
                subdomain_particle_velocity(i,:)=subdomain_position_soft_coefficient(i)*subdomain_WaveMaker_Analytical_velocity(i,:)+(1-subdomain_position_soft_coefficient(i))*subdomain_particle_velocity(i,:)
                subdomain_particle_position(i,:)=subdomain_position_soft_coefficient(i)*subdomain_WaveMaker_Analytical_position(i,:)+(1-subdomain_position_soft_coefficient(i))*subdomain_particle_position(i,:)

                subdomain_particle_press(i)=subdomain_position_soft_coefficient(i)*subdomain_WaveMaker_Analytical_press(i)+(1-subdomain_position_soft_coefficient(i))*subdomain_particle_press(i)   
                subdomain_particle_rho(i)=water_rho_0*(1+7*subdomain_particle_press(i)/(square_c_0*water_rho_0))**(0.142857142)

                !subdomain_particle_rho(i)=subdomain_position_soft_coefficient(i)*subdomain_WaveMaker_Analytical_rho(i)+(1-subdomain_position_soft_coefficient(i))*subdomain_particle_rho(i)
                
                ! subdomain_particle_velocity(i,:)=subdomain_WaveMaker_Analytical_velocity(i,:)
                ! subdomain_particle_position(i,:)=subdomain_WaveMaker_Analytical_position(i,:)

            endif 

        enddo
        !------------------------------------------------------------------------------------------------------

        ! write(*,*) L,FluidParticleNumberInWaveMakerDamping
        
    end subroutine wave_maker_damping_zone


    !==========================================================================================================


end module Wave_Gereration_Module
