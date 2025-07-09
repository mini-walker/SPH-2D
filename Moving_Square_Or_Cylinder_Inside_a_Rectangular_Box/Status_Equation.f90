!****************************************************************************
!
!  PROGRAM: status_equation
!
!  Programer: Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location: MUN
!
!  Time£º2017.3.18
!
!  Copyright: Memorial University
!
!  Version: 1.0
!
!  Note: MPI version: mpich-3.2
!        Fortran version: Inter Fortran (ifort) 18.0.0
!****************************************************************************

subroutine   status_equation(particle_rho,&                 !particle density
                             particle_energy,&              !particle energy
                             particle_type,&                !particle type
                             particle_press,&               !particle pressure(out)
                             particle_c&                    !particle sound velocity(out)
                             )

    use information_module
    
    implicit none
    
    ! Variables form mother subroutine

    integer,intent(in)::particle_type                      !particle type
    real(kind=8),intent(in)::particle_rho                  !particle density
    real(kind=8),intent(in)::particle_energy               !particle energy
    real(kind=8),intent(out)::particle_c                   !particle c
    real(kind=8),intent(out)::particle_press               !particle pressure
    
    ! Variables in subroutine
    real(kind=8)::gamma
    real(kind=8)::b                                        !initial reference pressure
    
    ! Body of subroutine status_equation
    
    !Equations of Status
    !Air EOS
    if(abs(particle_type)==1) then

        gamma=1.4
        !particle_press=(gamma-1)*particle_rho*particle_internal_energy
        particle_c=sqrt((gamma-1)*particle_rho)

    !Water EOS
    else 

        ! Tait equations
        if ( EquationOfState==0 ) then

            ! ------------------------------------------------------------------------
            ! 1st order EOS(1)
            ! 1st order EOS: P=B*((rho/rho_0)^7-1.0)
            !                B=C_0^2*rho_0/7.0
            ! ------------------------------------------------------------------------
            gamma=7.0
            b=square_c_0*water_rho_0/7.0
            particle_press=b*((particle_rho/water_rho_0)**7-1.0)+Background_Pressure                                         
            particle_c=sqrt(square_c_0)*(particle_rho/water_rho_0)**3

        ! Tammann equations
        elseif( EquationOfState==1 ) then

            ! ------------------------------------------------------------------------
            ! 1st order EOS(2)
            ! 1st order EOS: P=C_0^2*(rho-rho_0)
            ! ------------------------------------------------------------------------
            particle_press=square_c_0*(particle_rho-water_rho_0)+Background_Pressure                                           
            particle_c=sqrt(square_c_0)*(particle_rho/water_rho_0)**3

            ! ------------------------------------------------------------------------

            ! ------------------------------------------------------------------------
            ! 2nd order EOS: P=C_0^2*(1+cap_gamma*(particle_energy/e_0-1.0))*(rho-rho_0)
            ! ------------------------------------------------------------------------
            ! particle_press=square_c_0*(1+cap_gamma*(particle_energy/e_0-1.0))*(particle_rho-water_rho_0)+Background_Pressure !2nd order EOS
            ! particle_c=sqrt(square_c_0)*sqrt(1+cap_gamma*(particle_energy/e_0-1.0))

        else

            write(*,*) 'The defination of EquationOfState is not right! (EOS)'
            
        endif


        
    end if
  
end subroutine status_equation 
    