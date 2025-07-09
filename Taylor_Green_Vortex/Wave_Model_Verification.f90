!****************************************************************************
!
!  PROGRAM: CSPH_MPI_Parallel_wave_maker_2D
!
!  SUBROUTINE:: Wave_Model_Verification
!
!  PURPOSE: Check the wave model for the calculation based on the Ursell Number
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

subroutine Wave_Model_Verification()

    use Public_variable_module
    use information_module
    use MPI

    implicit none

    !==========================================================================================================
    !Variables in subroutine
    integer::i,j,k,L,m                                           !loop Variables
  
    real(kind=8)::Ursell_Number                                  !Ursell Number
    real(kind=8)::Wave_Steepnes                                  !Ursell Number
    real(kind=8)::WaterDepth_over_WaveLength                     !water depth over wave length
    real(kind=8)::Assume_WaveLength                              !Assume Wave Length
    real(kind=8)::current_time                                   !current time
    real(kind=8)::current_wave_angle                             !current wave angle
    real(kind=8)::current_wave_elevation                         !wave elevation at current particle z direction
    real(kind=8)::current_wave_probe_position                    !current wave probe position
    real(kind=8)::temp                                           !temp value
    real(kind=8)::Function_value,Function_derivate_value         !Function value and its derivates
    integer::Iteration_number                                    !Iteration number

    !Variables for file opeartion
    integer::ioerror=0                                                         !open return value
    integer::stat                                                              !read return value
    integer::status                                                            !allocation return value
    character(len=40)::file_name                                               !file name
    character(len=40)::char_Current_wave_probe_location                        !char Current wave probe location

    !==========================================================================================================

    ! Body of subroutine Wave_Model_Verification

    Solitary_wave_or_not:if (trim(adjustl(wave_type))=='Solitarywave') then
        
        !=======================================================================================================
        !For Solitary wave
        !Reference: JOHN FENTON 'A ninth-order solution for the solitary wave'
        Solitary_Epsilon=wave_height/water_depth
        Solitary_alpha=sqrt(0.75*Solitary_Epsilon)*(1-0.625*Solitary_Epsilon+71/128.0*Solitary_Epsilon**2)

        WavemakerDampingLengh=WavemakerDampingLenghScale*wave_length

        !The main processor output the data
        if (Current_Processor_ID==Main_Processor) then

            write(*,*) "==============================================================="
            !-------------------------------------------------------------------------
            write(*,*) "Wave model is 'Third order solitary wave' !"
            write(*,*) "Wave Height   :",wave_Height,"(m)"
            write(*,*) "Wave Length   :",Wave_Length,"(m)---initial value"

            !-------------------------------------------------------------------------
            write(*,*) "==============================================================="

        endif
        
    else
        !=======================================================================================================
        ! Linear wave theory 
        ! Generate the basic infomration for the input wave
        ! wave length; wave period; wave frequency; wave height; wave number
        if(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WavePeriod') then
            
            Wave_Frequency=2*PI/wave_Period

            !----------------------------------------------------------------------
            !Calculate the wave length
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
            !----------------------------------------------------------------------

            Wave_number=2*PI/Wave_Length
            Wave_amplitue=wave_Height/2.0
        
        elseif(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WaveLength') then
            
            !----------------------------------------------------------------------
            Wave_number=2*PI/Wave_Length
            Wave_Frequency=sqrt(g*Wave_number*tanh(Wave_number*water_depth))
            wave_Period=2*PI/Wave_Frequency
            Wave_amplitue=wave_Height/2.0

            !----------------------------------------------------------------------

        elseif(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WaveFreuency') then
            
            !----------------------------------------------------------------------
            
            wave_Period=2*PI/Wave_Frequency
            Wave_amplitue=wave_Height/2.0

            !----------------------------------------------------------------------
            !Calculate the wave length
            Assume_WaveLength=1.0d0

            do
                Wave_Length=g*wave_Period**2/(2*PI)*tanh(2*PI*water_depth/Assume_WaveLength)

                ! When the Assume wave length is nearly same with Wave Length result, exit
                if (abs(Wave_Length-Assume_WaveLength)<1.0E-8) then
                    exit
                else
                    Assume_WaveLength=Wave_Length
                endif

            end do
            !----------------------------------------------------------------------

            Wave_number=2*PI/Wave_Length
            Wave_amplitue=wave_Height/2.0

            !----------------------------------------------------------------------

        else

            if (Current_Processor_ID==Main_Processor) then
                 write(*,*) "The Input Wave Information type is not include in the code, please check the type."
            endif

        endif
        !-------------------------------------------------------------------------
     
        !-------------------------------------------------------------------------
        !Calculate the h/L for the wave checking: deep water wave; shallow water wave and finite water wave
        WaterDepth_over_WaveLength=water_depth/Wave_Length
        Ursell_Number=wave_Height/water_depth*(Wave_Length/water_depth)**2
        Wave_Steepnes=wave_Height/Wave_Length

        !The main processor output the data
        if (Current_Processor_ID==Main_Processor) then

            !-------------------------------------------------------------------------
            write(*,*) "WaterDepth/WaveLength:",WaterDepth_over_WaveLength
            write(*,*) "Ursell Number:",Ursell_Number

            if (WaterDepth_over_WaveLength>=0.5) then

                write(*,*) "Wave Type: Deep water waves"

            elseif (WaterDepth_over_WaveLength<=0.05) then

                write(*,*) "Wave Type: Shallow water waves"

            else

                write(*,*) "Wave Type: Finite water waves"
                
            endif
            !-------------------------------------------------------------------------

            !-------------------------------------------------------------------------
            ! calculate the Ursell Number for different model
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
            write(*,*) "==============================================================="

            !-------------------------------------------------------------------------

            !-------------------------------------------------------------------------
            !Recommand wave theory based on Second Theory

            !For Airy wave
            if (wave_length<2*water_depth .and. wave_Height/wave_length<=0.06) then
                
                write(*,'(A)') " Recommand wave theory: Airy Wave Theory. (Based on wave length and period)"

            endif

            !For stokes wave
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

            !For Cnoidal wave
            if (wave_length>8*water_depth .and. water_depth/wave_Height>1.8 .and. water_depth/wave_Height<2.2 .and. wave_Height/(g*wave_Period**2)<=0.02 .and. water_depth/(g*wave_Period**2)<=0.03) then
                
                write(*,'(A)') " Recommand wave theory: Cnoidal Wave Theory. (Based on wave length and period)"

            endif


            write(*,*) "==============================================================="
            !-------------------------------------------------------------------------

        end if
        !==========================================================================================================


        !==========================================================================================================
        !Recalculate the wave information again based on the wave theory you used
        !For nonliear wave (3th or 5th Stokes waves), we should calculate them again
        !For linear wave, we just keep same
        !For 5th order Stokes, the reference:"Hydrodynamics of Offshore Structures.P59"
        if(trim(adjustl(wave_type))=='Stokes_5th') then

            !--------------------------------------------------------------------------
            !Calculate the coefficients for 5th order for stokes wave
            Fifth_stokes_A=0.0d0  
            Fifth_stokes_B=0.0d0  
            Fifth_stokes_C=0.0d0 
            Fifth_stokes_D=0.0d0   
            Fifth_stokes_E=0.0d0 
            !--------------------------------------------------------------------------

            if(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WavePeriod') then
                
                Wave_Frequency=2*PI/wave_Period

                !----------------------------------------------------------------------
                !Determine wave number through Equ.(3.89) and Equ.(3.90)
                
                GET_k_1:do i=1,1000000

                    Last_wave_number=wave_number                                      !wave number in last step
                    
                    !------------------------------------------------------------------
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
                    !------------------------------------------------------------------

                    !------------------------------------------------------------------
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
                        
                    !------------------------------------------------------------------

                    !Check the results
                    if (Current_Processor_ID==Main_Processor .and. i==1000000) then

                        write(*,'(A)') ' Errors in 5th order stokes wave information calculation. (Wave_Model_Verification)'

                    end if

                enddo GET_k_1
                !----------------------------------------------------------------------

                Wave_Length=2*PI/Wave_number
                Wave_celerity=2*PI/(wave_number*wave_Period)                             !c=2PI/(KT)
                Wave_amplitue=wave_Height/2.0
                Wave_celerity=2*PI/(wave_number*wave_Period)                             !c=2PI/(KT)
                
                Wave_celerity_0=sqrt((g/wave_number)*tanh(wave_number*water_depth))     !c=sqrt(g/K*tanh(kd))
            
            elseif(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WaveLength') then
                
                !----------------------------------------------------------------------
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
                !------------------------------------------------------------------


                Wave_celerity=Wave_celerity_0*sqrt(1+Fifth_stokes_C(1)*Stokes_Lamda**2+Fifth_stokes_C(2)*Stokes_Lamda**4)

                Wave_Frequency=wave_number*Wave_celerity           !w=kc

                wave_Period=2*PI/Wave_Frequency                    !T=2PI/w

                !----------------------------------------------------------------------

            elseif(trim(adjustl(InputWaveInformation))=='WaveHeight_and_WaveFreuency') then
                
                !----------------------------------------------------------------------
                
                wave_Period=2*PI/Wave_Frequency
                Wave_amplitue=wave_Height/2.0

                !----------------------------------------------------------------------
                !Determine wave number through Equ.(3.89) and Equ.(3.90)
                
                GET_k_2:do i=1,1000000

                    Last_wave_number=wave_number                                      !wave number in last step
                    
                    !------------------------------------------------------------------
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
                    !------------------------------------------------------------------

                    !------------------------------------------------------------------
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
                        
                    !------------------------------------------------------------------

                    !Check the results
                    if (Current_Processor_ID==Main_Processor .and. i==1000000) then

                        write(*,'(A)') ' Errors in 5th order stokes wave information calculation. (Wave_Model_Verification)'

                    end if

                enddo GET_k_2
                !----------------------------------------------------------------------

                !----------------------------------------------------------------------
                Wave_Length=2*PI/Wave_number
                Wave_celerity=2*PI/(wave_number*wave_Period)                             !c=2PI/(KT)

                Wave_celerity_0=sqrt((g/wave_number)*tanh(wave_number*water_depth))     !c=sqrt(g/K*tanh(kd))
                !----------------------------------------------------------------------

            else

                if (Current_Processor_ID==Main_Processor) then
                     write(*,*) "The Input Wave Information type is not include in the code, please check the type."
                endif

            endif
            !-------------------------------------------------------------------------

            !-------------------------------------------------------------------------
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
            !-------------------------------------------------------------------------

            !-------------------------------------------------------------------------
            !For pressure calculation
            Stokes_pressure_S=1.0/cosh(2*wave_number*water_depth)

            Stokes_pressure_C_0=sqrt(tanh(wave_number*water_depth))
            Stokes_Epsilon=0.5*wave_number*wave_height
            Fifth_stokes_pressure_E=0

            Fifth_stokes_pressure_E(2)=tanh(wave_number*water_depth)*(2+2*Stokes_pressure_S+5*Stokes_pressure_S**2)/(4*(1-Stokes_pressure_S)**2)
            Fifth_stokes_pressure_E(4)=tanh(wave_number*water_depth)*(8+12*Stokes_pressure_S-152*Stokes_pressure_S**2-308*Stokes_pressure_S**3-42*Stokes_pressure_S**4+77*Stokes_pressure_S**5)/(32*(1-Stokes_pressure_S)**5)
                    
            Bernoulli_constant_R=g/wave_number*(0.5*Stokes_pressure_C_0**2+wave_number*water_depth+Stokes_Epsilon**2*Fifth_stokes_pressure_E(2)+Stokes_Epsilon**4*Fifth_stokes_pressure_E(4))

            !-------------------------------------------------------------------------

        else if(trim(adjustl(wave_type))=='Stokes_3rd') then

            !-------------------------------------------------------------------------
            !Calculate the wave number based on Equ(3.86) in "Hydrodynamics of Offshore Structures.P58"
            !c^2=c_0^2*(1+(PI*H/L)^2*((9-8*coshkd^2+8*coshkd^4)/(8*sinhkd^4)))

            !-------------------------------------------------------------------------
            !Use the wave number based on linear theory
            !wave celerity from the linear wave theory
            Wave_celerity_0=sqrt(g/wave_number*tanh(wave_number*water_depth))

            Wave_celerity=Wave_celerity_0*sqrt(1+(PI*wave_height/wave_length)**2*(9-8*(cosh(wave_number*water_depth))**2+8*(cosh(wave_number*water_depth))**4)/(8*(sinh(wave_number*water_depth))**4))

            !Calculate the coefficients for third order stokes wave
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

            !-------------------------------------------------------------------------

        else if(trim(adjustl(wave_type))=='Cnoidal_5th') then



        endif

        !-------------------------------------------------------------------------
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
        !-------------------------------------------------------------------------
   
    endif Solitary_wave_or_not
    !==========================================================================================================

    !Assign the Wave maker Damping Lengh based on actual wave length calculate from the wavemaker subroutine
    WavemakerDampingLengh=WavemakerDampingLenghScale*wave_length

    !==========================================================================================================
    !Generate the analytical solution for regular waves
    
    if (Current_Processor_ID==Main_Processor) then

        do k=1,Wave_Probe_number

            !transfer the Current Wave Probe Location from real to character
            write(char_Current_wave_probe_location,'(f6.3)') Wave_Probe_Position(k,1)

            file_name=trim(adjustl(wave_type))//'_'//trim(adjustl(char_Current_wave_probe_location))//"_Analytical_Results.dat"

            !Open the saving file
            open(unit=1,file=trim(adjustl(file_name)),status="replace",position="rewind",action="write",iostat=ioerror)    
            
            if(ioerror==0) then
               
                write(*,*) "The open action of the "// trim(adjustl(file_name))//" is successful!"
                
                !Calculate the wave elevations
                if(trim(adjustl(wave_type))=='Stokes_2nd') then
                    
                    !----------------------------------------------------------------------
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

                    end do

                elseif(trim(adjustl(wave_type))=='Airy') then
                    
                    !----------------------------------------------------------------------
                    !Wave_elevation=H/2*cos(kx-wt+Fai0)

                    current_wave_probe_position=Wave_Probe_Position(k,1)               ! Unit: m

                    do i=1,max_time_step,50

                        current_time=dt*i

                        current_wave_angle=wave_number*current_wave_probe_position-wave_frequency*current_time+Wave_initial_phase_angle
                        current_wave_elevation=water_depth+0.5*wave_height*cos(current_wave_angle)
             
                        write(1,*) current_time,current_wave_elevation

                    end do

                elseif(trim(adjustl(wave_type))=='Stokes_5th') then

                    !----------------------------------------------------------------------
                    !Wave_elevation=H/2*cos(kx-wt+Fai0)

                    current_wave_probe_position=Wave_Probe_Position(k,1)               ! Unit: m

                    ! !----------------------------------------------------------------------
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
                    ! !----------------------------------------------------------------------

                    !----------------------------------------------------------------------
                    !Brorsen, Michael "Nonlinear Wave Theories"
                    do i=1,max_time_step,50

                        current_time=dt*i

                        current_wave_angle=wave_number*current_wave_probe_position-wave_frequency*current_time+Wave_initial_phase_angle
                        
                        current_wave_elevation=0.0d0
                        do m=1,5
                            current_wave_elevation=current_wave_elevation+Fifth_stokes_E(m)*cos(m*current_wave_angle)
                        enddo
                        current_wave_elevation=water_depth+current_wave_elevation/wave_number
             
                        write(1,*) current_time,current_wave_elevation

                    end do 
                    !----------------------------------------------------------------------

                elseif(trim(adjustl(wave_type))=='Stokes_3rd') then


                    current_wave_probe_position=Wave_Probe_Position(k,1)               ! Unit: m

                    !----------------------------------------------------------------------
                    !Reference: Song Zhiyao "On the universal third order stokes wave solution"
                    do i=1,max_time_step,50

                        current_time=dt*i

                        current_wave_angle=wave_number*current_wave_probe_position-wave_frequency*current_time+Wave_initial_phase_angle
                        
                        current_wave_elevation=0.0d0
                        do m=1,3
                            current_wave_elevation=current_wave_elevation+Third_stokes_E(m)*cos(m*current_wave_angle)
                        enddo
                        current_wave_elevation=water_depth+current_wave_elevation/wave_number
             
                        write(1,*) current_time,current_wave_elevation

                    end do 
                    !----------------------------------------------------------------------

                elseif(trim(adjustl(wave_type))=='Solitarywave') then

                    current_wave_probe_position=Wave_Probe_Position(k,1)               ! Unit: m

                    !----------------------------------------------------------------------
                    !Reference: Song Zhiyao "On the universal third order stokes wave solution"
                    do i=1,max_time_step,50

                        current_time=dt*i

                        current_wave_angle=wave_number*current_wave_probe_position-wave_frequency*current_time+Wave_initial_phase_angle
                        
                        current_wave_elevation=0.0d0
                        do m=1,3
                            current_wave_elevation=current_wave_elevation+Third_stokes_E(m)*cos(m*current_wave_angle)
                        enddo
                        current_wave_elevation=water_depth+current_wave_elevation/wave_number
             
                        write(1,*) current_time,current_wave_elevation

                    end do 
                    !----------------------------------------------------------------------

                else

                    write(*,*) "The wave type is not include in the code, please review the type."

               endif

            else
                write(*,*) "The open action of the "// trim(adjustl(file_name))//" is fail!"
            end if
           
            close(1)
            !-------------------------------------------------------------------------------

        end do
    endif
    !==========================================================================================================

end subroutine Wave_Model_Verification

