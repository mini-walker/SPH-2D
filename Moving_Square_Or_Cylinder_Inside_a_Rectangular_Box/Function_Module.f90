!**************************************************************************************************************
!  SUBROUTINE: Function_Module
!
!  PURPOSE: Function module for SPH calculation 
!          (contains: the linear equation solver, matrix inverse and so on)
!
!  Programer: Shanqin Jin (E-mail:sjin@mun.ca)
!
!  Location: MUN
!
!  Time: 2017.3.18
!
!  Copyright: Memorial University
!
!  Version: 1.0
!
!  Note: MPI version: mpich-3.2
!        Fortran version: Inter Fortran (ifort) 18.0.0
!**************************************************************************************************************

module Function_Module

    use  information_module
 
    implicit none
  
    contains

    !**********************************************************************************************************
    subroutine MUSCL_Limiter(delphi_i,delphi_j,d_phi_ij,del_phi_lim_i,del_phi_lim_j)

      implicit none

      real(kind=8),intent(in)::delphi_i,delphi_j,d_phi_ij
      real(kind=8),intent(out)::del_phi_lim_i,del_phi_lim_j

      real(kind=8)::diff_factor,del_phi_lim,diff_factor_i,diff_factor_j
      real(kind=8)::beta_lim=1.3d0

      ! Body of Subroutine MUSCL_Limiter

      !- MUSCL Limiters - Toro (2001) p. 208
      diff_factor = 0.5d0

      if( delphi_j.gt.0.0d0 )then

        del_phi_lim = max(0.0d0 , min(beta_lim*delphi_i,delphi_j), min(delphi_i,beta_lim*delphi_j) )

        if( del_phi_lim .gt. (d_phi_ij) .and. abs(del_phi_lim).gt. 0.0d0 )then
              
            if( delphi_i .gt. delphi_j )then
              diff_factor_i = diff_factor
              diff_factor_j = 0.0
            else
              diff_factor_i = 0.0
              diff_factor_j = diff_factor
            endif

            del_phi_lim_i = d_phi_ij*diff_factor_i
            del_phi_lim_j = d_phi_ij*diff_factor_j

        else

            del_phi_lim_i = del_phi_lim
            del_phi_lim_j = del_phi_lim

        endif

      else

        del_phi_lim =min(0.0d0 ,max(beta_lim*delphi_i,delphi_j),max(delphi_i,beta_lim*delphi_j))

        if( del_phi_lim .lt. (d_phi_ij) .and. abs(del_phi_lim).gt.0.0d0 )then

            if(delphi_i.lt.delphi_j)then
              diff_factor_i = diff_factor
              diff_factor_j = 0.0
            else
              diff_factor_i = 0.0
              diff_factor_j = diff_factor
            endif

            del_phi_lim_i = d_phi_ij*diff_factor_i
            del_phi_lim_j = d_phi_ij*diff_factor_j

        else

            del_phi_lim_i = del_phi_lim
            del_phi_lim_j = del_phi_lim

        endif

      end if

    end subroutine MUSCL_Limiter
    !**********************************************************************************************************

    !**********************************************************************************************************
    ! Input Variables: Delta_Backward, Delta_Forward
    subroutine Upwind_MUSCL_Limiter(delphi_i,delphi_j,del_phi_lim)

      implicit none

      real(kind=8),intent(in)::delphi_i,delphi_j
      real(kind=8),intent(out)::del_phi_lim

      real(kind=8)::d_phi_ij
      real(kind=8)::diff_factor,diff_factor_i,diff_factor_j

      real(kind=8)::beta_lim=1.3d0                                      !The small belta, more smooth (1.0-1.5 is best)

      ! Body of Subroutine Upwind_MUSCL_Limiter

      if( delphi_j.gt.0.0d0 )then

        del_phi_lim = max(0.0d0 , min(beta_lim*delphi_i,delphi_j), min(delphi_i,beta_lim*delphi_j) )

      else

        del_phi_lim =min(0.0d0 ,max(beta_lim*delphi_i,delphi_j),max(delphi_i,beta_lim*delphi_j))

      end if

    end subroutine Upwind_MUSCL_Limiter
    !**********************************************************************************************************

    !**********************************************************************************************************
    ! Transfer ASCII to character
    subroutine Write_ASCII_to_BINARY(File_Port,ASCII_Character,character_length)
       
      implicit none

        integer,intent(in)::character_length
        integer,intent(in)::File_Port
        character(len=character_length)::ASCII_Character
        
        integer::i,j,k
        
           do i=1,character_length
                j=ICHAR(ASCII_Character(i:i))
            write(File_Port) j
           end do
           write(File_Port) 0

    end subroutine Write_ASCII_to_BINARY
    !**********************************************************************************************************
    
    !**********************************************************************************************************
    ! Matric Cross product
    Function Cross_Product(A,B,n) result(Z)

    ! CROSS_3D computes the cross product of two vectors in 3D.
    !
    !  Definition: X Corss_Product Y
    !
    !    The cross product in 3D can be regarded as the determinant of the
    !    symbolic matrix:
    !
    !          |  i  j  k |
    !      det | x1 y1 z1 |
    !          | x2 y2 z2 |
    !
    !      = ( y1 * z2 - z1 * y2 ) * i
    !      + ( z1 * x2 - x1 * z2 ) * j
    !      + ( x1 * y2 - y1 * x2 ) * k
    !  Parameters:
    !    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
    !    Output, real X3, Y3, Z3, the cross product vector.

      implicit none
      integer::n                                !Matrix dimension
      real(kind=8),dimension(n)::A              !Input A array
      real(kind=8),dimension(n)::B              !Input B array
      
      real(kind=8),dimension(3)::X              !X array
      real(kind=8),dimension(3)::Y              !Y array
      real(kind=8),dimension(3)::Z              !Z array

      ! Body of Function Corss_Product

        if (n==2) then

            X(1) = A(1) ; X(2) = A(2) ; X(3) = 0.0d0 ;
            Y(1) = B(1) ; Y(2) = B(2) ; Y(3) = 0.0d0 ;

        else

          X = A
          Y = B

        endif

        Z(1) = X(2) * Y(3) - X(3) * Y(2)
        Z(2) = X(3) * Y(1) - X(1) * Y(3)
        Z(3) = X(1) * Y(2) - X(2) * Y(1)

    end Function Cross_Product
    !**********************************************************************************************************

    !**********************************************************************************************************
    ! Matric Tensor product
    Function Tensor_Product(A,B,n) result(Z)

    ! CROSS_3D computes the Tensor product of two vectors in 3D.
    !
    !  Definition: A Tensor_Product B
    !
    !    The cross product in 3D can be regarded as the determinant of the
    !    symbolic matrix:
    !
    !          | a1 |
    !          | a2 | Tensor_Product | b1 b2 b3 |
    !          | a3 |
    !          | a1b1  a1b2  a1b3 |
    !      =   | a2b1  a2b2  a2b3 |
    !          | a3b1  a3b2  a3b3 |
    !


      implicit none
      integer::n,i,j                            !Matrix dimension
      real(kind=8),dimension(n)::A              !Input A array
      real(kind=8),dimension(n)::B              !Input B array
      real(kind=8),dimension(n,n)::Z            !Z tensor matrix


      ! Body of Function Tensor_Product

      do i=1,n
        do j=1,n

            Z(i,j)=A(i)*B(j)

        enddo
      enddo

    end Function Tensor_Product
    !**********************************************************************************************************

    !**********************************************************************************************************
    ! Matric Kronecker product
    Function Kronecker_Product(A,B,n) result(Z)

    ! CROSS_3D computes the Kronecker product of two vectors in 3D.
    !
    !  Definition: A Tensor_Product B
    !
    !    The cross product in 3D can be regarded as the determinant of the
    !    symbolic matrix:
    !
    !                                          | b1 |
    !          | a1 a2 a3 | Kronecker_Product  | b2 |  
    !                                          | b3 |
    !          | a1b1  a2b1  a3b1 |
    !      =   | a1b2  a2b2  a3b2 |
    !          | a1b3  a2b3  a3b3 |
    !


      implicit none
      integer::n,i,j                            !Matrix dimension
      real(kind=8),dimension(n)::A              !Input A array
      real(kind=8),dimension(n)::B              !Input B array
      real(kind=8),dimension(n,n)::Z            !Z tensor matrix


      ! Body of Function Kronecker_Product

      do i=1,n
        do j=1,n

            Z(i,j)=A(j)*B(i)

        enddo
      enddo

    end Function Kronecker_Product
    !**********************************************************************************************************

    !**********************************************************************************************************
    subroutine matrix_inversion_less_three(n,input_matrix,result_matrix,check_error)

      implicit none
    
      ! Variables
      integer,intent(in)::n                                             !Input matrix dimension: n*n
      real(kind=8),dimension(n,n),intent(in)::input_matrix              !Input matrix
      real(kind=8),dimension(n,n),intent(out)::result_matrix            !Result matrix
      integer::check_error
      
      !Coefficience for linear equations solution
      integer::k
      real(kind=8),dimension(n,n)::A                                    ! A
      real(kind=8)::martix_value
      real(kind=8),dimension(n,n)::Inverse_A                            ! Inverse martix of A
      real(kind=8)::a_1,a_2,a_3,b_1,b_2,b_3,c_1,c_2,c_3
          
      ! Body of subroutine matrix_inversion_less_three
        
      Inverse_A=0.0d0
      A=input_matrix
      check_error=1
      martix_value=0.0d0
      
      !********************************************************************************************************
      ! Initialized the results martix with identity matrix
      do k=1,n
          result_matrix(k,k)=1.0d0
      enddo
      !********************************************************************************************************
    
      !********************************************************************************************************
      ! 2D martix inversion
      if(n==2) then
           
          a_1=A(1,1); b_1=A(1,2)
          a_2=A(2,1); b_2=A(2,2)

          martix_value=(a_1*b_2-b_1*a_2)
          
          if(martix_value/=0) then

              Inverse_A(1,1)=b_2/martix_value
              Inverse_A(1,2)=-b_1/martix_value
              Inverse_A(2,1)=-a_2/martix_value
              Inverse_A(2,2)=a_1/martix_value
              
          end if
      !********************************************************************************************************

      !********************************************************************************************************
      ! 3D martix inversion
      elseif(n==3) then
       
          a_1=A(1,1); b_1=A(1,2); c_1=A(1,3)
          a_2=A(2,1); b_2=A(2,2); c_2=A(2,3)
          a_3=A(3,1); b_3=A(3,2); c_3=A(3,3)

          martix_value=(a_1*(b_2*c_3-c_2*b_3)-a_2*(b_1*c_3-c_1*b_3)+a_3*(b_1*c_2-c_1*b_2))
          
          if(martix_value/=0) then

              Inverse_A(1,1)=(b_2*c_3-c_2*b_3)/martix_value
              Inverse_A(1,2)=(b_3*c_1-c_3*b_1)/martix_value
              Inverse_A(1,3)=(b_1*c_2-c_1*b_2)/martix_value
              Inverse_A(2,1)=(c_2*a_3-a_2*c_3)/martix_value
              Inverse_A(2,2)=(c_3*a_1-a_3*c_1)/martix_value
              Inverse_A(2,3)=(c_1*a_2-a_1*c_2)/martix_value
              Inverse_A(3,1)=(a_2*b_3-b_2*a_3)/martix_value
              Inverse_A(3,2)=(a_3*b_1-b_3*a_1)/martix_value
              Inverse_A(3,3)=(a_1*b_2-b_1*a_2)/martix_value
          
          end if
      !********************************************************************************************************

      else
       
          write(*,*) "The matrix is larger than three order !"
       
      end if
      !********************************************************************************************************
     
      !Check the inversion opeartion is right or not
      if(martix_value/=0) then
          result_matrix=Inverse_A
      else
          check_error=0
          write(*,*) "The inversion calculation is wrong! "
      end if
          
      !result_matrix=Inverse_A
      
    end subroutine matrix_inversion_less_three
    !**********************************************************************************************************

    !**********************************************************************************************************
    ! Solve the inversion of Coefficient matrix directly, when n less or equal to 3 this is faster
    ! When dimension is larger than 3 this method is not good, however, this is also a choice
    ! Xu's subroutine faster than the Fortran inner code and accuracy
    subroutine matrix_inversion_less_six(n,input_matrix,result_matrix,check_error)

       implicit none
    
        ! Variables
        integer,intent(in)::n                                             !Input matrix dimension: n*n
        real(kind=8),dimension(n,n),intent(in)::input_matrix              !Input matrix
        real(kind=8),dimension(n,n),intent(out)::result_matrix            !Result matrix
        integer::check_error
      
        ! Variables for Coefficient
        integer::k
        real(kind=8),dimension(n,n)::A                                    ! A
        real(kind=8)::martix_value
        real(kind=8),dimension(n,n)::Inverse_A                            ! Inverse martix of A
        real(kind=8)::a_1,a_2,a_3,a_4,a_5,a_6,a_7,a_8,a_9
        real(kind=8)::b_1,b_2,b_3,b_4,b_5,b_6,b_7,b_8,b_9
        real(kind=8)::c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9
        real(kind=8)::d_1,d_2,d_3,d_4,d_5,d_6,d_7,d_8,d_9
        real(kind=8)::e_1,e_2,e_3,e_4,e_5,e_6,e_7,e_8,e_9
        real(kind=8)::f_1,f_2,f_3,f_4,f_5,f_6,f_7,f_8,f_9
        real(kind=8)::g_1,g_2,g_3,g_4,g_5,g_6,g_7,g_8,g_9
        real(kind=8)::h_1,h_2,h_3,h_4,h_5,h_6,h_7,h_8,h_9
        real(kind=8)::i_1,i_2,i_3,i_4,i_5,i_6,i_7,i_8,i_9
          
        Inverse_A=0.0d0
        A=input_matrix
        check_error=1
        martix_value=0.0d0
        
        ! ! Initialized the results martix with identity matrix
        ! do k=1,n
        !     result_matrix(k,k)=1.0d0
        ! end do

        !********************************************************************************************************
        ! 2D martix inversion
        if(n==2) then
             
            a_1=A(1,1); b_1=A(1,2)
            a_2=A(2,1); b_2=A(2,2)
    
            martix_value=(a_1*b_2-b_1*a_2)
            
            if(martix_value/=0) then
    
                Inverse_A(1,1)=b_2/martix_value
                Inverse_A(1,2)=-b_1/martix_value
                Inverse_A(2,1)=-a_2/martix_value
                Inverse_A(2,2)=a_1/martix_value
                
            end if
        !********************************************************************************************************

        !********************************************************************************************************
        ! 3D martix inversion
        else if(n==3) then
         
            a_1=A(1,1); b_1=A(1,2); c_1=A(1,3)
            a_2=A(2,1); b_2=A(2,2); c_2=A(2,3)
            a_3=A(3,1); b_3=A(3,2); c_3=A(3,3)
    
            martix_value=(a_1*(b_2*c_3-c_2*b_3)-a_2*(b_1*c_3-c_1*b_3)+a_3*(b_1*c_2-c_1*b_2))
            
            if(martix_value/=0) then
    
                Inverse_A(1,1)=(b_2*c_3-c_2*b_3)/martix_value
                Inverse_A(1,2)=(b_3*c_1-c_3*b_1)/martix_value
                Inverse_A(1,3)=(b_1*c_2-c_1*b_2)/martix_value
                Inverse_A(2,1)=(c_2*a_3-a_2*c_3)/martix_value
                Inverse_A(2,2)=(c_3*a_1-a_3*c_1)/martix_value
                Inverse_A(2,3)=(c_1*a_2-a_1*c_2)/martix_value
                Inverse_A(3,1)=(a_2*b_3-b_2*a_3)/martix_value
                Inverse_A(3,2)=(a_3*b_1-b_3*a_1)/martix_value
                Inverse_A(3,3)=(a_1*b_2-b_1*a_2)/martix_value
            
            end if
        !********************************************************************************************************
            
        !********************************************************************************************************   
        ! 4D martix inversion
        else if(n==4) then
         
            a_1=A(1,1); b_1=A(1,2); c_1=A(1,3); d_1=A(1,4)
            a_2=A(2,1); b_2=A(2,2); c_2=A(2,3); d_2=A(2,4)
            a_3=A(3,1); b_3=A(3,2); c_3=A(3,3); d_3=A(3,4)
            a_4=A(4,1); b_4=A(4,2); c_4=A(4,3); d_4=A(4,4)

            martix_value=(a_1*b_2*c_3*d_4-a_1*b_2*c_4*d_3-a_1*b_3*c_2*d_4+a_1*b_3*c_4*d_2+a_1*b_4*c_2*d_3-a_1*b_4*c_3*d_2-a_2*b_1*c_3*d_4 &
                          +a_2*b_1*c_4*d_3+a_2*b_3*c_1*d_4-a_2*b_3*c_4*d_1-a_2*b_4*c_1*d_3+a_2*b_4*c_3*d_1+a_3*b_1*c_2*d_4 &
                          -a_3*b_1*c_4*d_2-a_3*b_2*c_1*d_4+a_3*b_2*c_4*d_1+a_3*b_4*c_1*d_2-a_3*b_4*c_2*d_1-a_4*b_1*c_2*d_3 & 
                          +a_4*b_1*c_3*d_2+a_4*b_2*c_1*d_3-a_4*b_2*c_3*d_1-a_4*b_3*c_1*d_2+a_4*b_3*c_2*d_1) 
            
            if(martix_value/=0) then
    
                Inverse_A(1,1)=(b_2*c_3*d_4 - b_2*c_4*d_3 - b_3*c_2*d_4 + b_3*c_4*d_2 + b_4*c_2*d_3 - b_4*c_3*d_2)/martix_value
                Inverse_A(2,1)=-(a_2*c_3*d_4 - a_2*c_4*d_3 - a_3*c_2*d_4 + a_3*c_4*d_2 + a_4*c_2*d_3 - a_4*c_3*d_2)/martix_value
                Inverse_A(3,1)=(a_2*b_3*d_4 - a_2*b_4*d_3 - a_3*b_2*d_4 + a_3*b_4*d_2 + a_4*b_2*d_3 - a_4*b_3*d_2)/martix_value
                Inverse_A(4,1)=-(a_2*b_3*c_4 - a_2*b_4*c_3 - a_3*b_2*c_4 + a_3*b_4*c_2 + a_4*b_2*c_3 - a_4*b_3*c_2)/martix_value
                
                Inverse_A(1,2)=-(b_1*c_3*d_4 - b_1*c_4*d_3 - b_3*c_1*d_4 + b_3*c_4*d_1 + b_4*c_1*d_3 - b_4*c_3*d_1)/martix_value
                Inverse_A(2,2)=(a_1*c_3*d_4 - a_1*c_4*d_3 - a_3*c_1*d_4 + a_3*c_4*d_1 + a_4*c_1*d_3 - a_4*c_3*d_1)/martix_value
                Inverse_A(3,2)=-(a_1*b_3*d_4 - a_1*b_4*d_3 - a_3*b_1*d_4 + a_3*b_4*d_1 + a_4*b_1*d_3 - a_4*b_3*d_1)/martix_value
                Inverse_A(4,2)=(a_1*b_3*c_4 - a_1*b_4*c_3 - a_3*b_1*c_4 + a_3*b_4*c_1 + a_4*b_1*c_3 - a_4*b_3*c_1)/martix_value
                
                Inverse_A(1,3)=(b_1*c_2*d_4 - b_1*c_4*d_2 - b_2*c_1*d_4 + b_2*c_4*d_1 + b_4*c_1*d_2 - b_4*c_2*d_1)/martix_value
                Inverse_A(2,3)=-(a_1*c_2*d_4 - a_1*c_4*d_2 - a_2*c_1*d_4 + a_2*c_4*d_1 + a_4*c_1*d_2 - a_4*c_2*d_1)/martix_value
                Inverse_A(3,3)=(a_1*b_2*d_4 - a_1*b_4*d_2 - a_2*b_1*d_4 + a_2*b_4*d_1 + a_4*b_1*d_2 - a_4*b_2*d_1)/martix_value
                Inverse_A(4,3)=-(a_1*b_2*c_4 - a_1*b_4*c_2 - a_2*b_1*c_4 + a_2*b_4*c_1 + a_4*b_1*c_2 - a_4*b_2*c_1)/martix_value
                
                Inverse_A(1,4)=-(b_1*c_2*d_3 - b_1*c_3*d_2 - b_2*c_1*d_3 + b_2*c_3*d_1 + b_3*c_1*d_2 - b_3*c_2*d_1)/martix_value
                Inverse_A(2,4)=(a_1*c_2*d_3 - a_1*c_3*d_2 - a_2*c_1*d_3 + a_2*c_3*d_1 + a_3*c_1*d_2 - a_3*c_2*d_1)/martix_value
                Inverse_A(3,4)=-(a_1*b_2*d_3 - a_1*b_3*d_2 - a_2*b_1*d_3 + a_2*b_3*d_1 + a_3*b_1*d_2 - a_3*b_2*d_1)/martix_value
                Inverse_A(4,4)=(a_1*b_2*c_3 - a_1*b_3*c_2 - a_2*b_1*c_3 + a_2*b_3*c_1 + a_3*b_1*c_2 - a_3*b_2*c_1)/martix_value
            
            end if
        !********************************************************************************************************

        !********************************************************************************************************
        ! 5D martix inversion
        else if(n==5) then
         
          a_1=A(1,1); b_1=A(1,2); c_1=A(1,3); d_1=A(1,4); e_1=A(1,5)
          a_2=A(2,1); b_2=A(2,2); c_2=A(2,3); d_2=A(2,4); e_2=A(2,5) 
          a_3=A(3,1); b_3=A(3,2); c_3=A(3,3); d_3=A(3,4); e_3=A(3,5)
          a_4=A(4,1); b_4=A(4,2); c_4=A(4,3); d_4=A(4,4); e_4=A(4,5)
          a_5=A(5,1); b_5=A(5,2); c_5=A(5,3); d_5=A(5,4); e_5=A(5,5)

  
          martix_value=(a_1*b_2*c_3*d_4*e_5 - a_1*b_2*c_3*d_5*e_4 - a_1*b_2*c_4*d_3*e_5 + a_1*b_2*c_4*d_5*e_3 + a_1*b_2*c_5*d_3*e_4 - a_1*b_2*c_5*d_4*e_3 - a_1*b_3*c_2*d_4*e_5 + a_1*b_3*c_2*d_5*e_4 + a_1*b_3*c_4*d_2*e_5 - a_1*b_3*c_4*d_5*e_2 - a_1*b_3*c_5*d_2*e_4 + a_1*b_3*c_5*d_4*e_2 + a_1*b_4*c_2*d_3*e_5 - a_1*b_4*c_2*d_5*e_3 - a_1*b_4*c_3*d_2*e_5 + a_1*b_4*c_3*d_5*e_2 + a_1*b_4*c_5*d_2*e_3 - a_1*b_4*c_5*d_3*e_2 - a_1*b_5*c_2*d_3*e_4 + a_1*b_5*c_2*d_4*e_3 &
                        + a_1*b_5*c_3*d_2*e_4 - a_1*b_5*c_3*d_4*e_2 - a_1*b_5*c_4*d_2*e_3 + a_1*b_5*c_4*d_3*e_2 - a_2*b_1*c_3*d_4*e_5 + a_2*b_1*c_3*d_5*e_4 + a_2*b_1*c_4*d_3*e_5 - a_2*b_1*c_4*d_5*e_3 - a_2*b_1*c_5*d_3*e_4 + a_2*b_1*c_5*d_4*e_3 + a_2*b_3*c_1*d_4*e_5 - a_2*b_3*c_1*d_5*e_4 - a_2*b_3*c_4*d_1*e_5 + a_2*b_3*c_4*d_5*e_1 + a_2*b_3*c_5*d_1*e_4 - a_2*b_3*c_5*d_4*e_1 - a_2*b_4*c_1*d_3*e_5 + a_2*b_4*c_1*d_5*e_3 + a_2*b_4*c_3*d_1*e_5 - a_2*b_4*c_3*d_5*e_1 &
                        - a_2*b_4*c_5*d_1*e_3 + a_2*b_4*c_5*d_3*e_1 + a_2*b_5*c_1*d_3*e_4 - a_2*b_5*c_1*d_4*e_3 - a_2*b_5*c_3*d_1*e_4 + a_2*b_5*c_3*d_4*e_1 + a_2*b_5*c_4*d_1*e_3 - a_2*b_5*c_4*d_3*e_1 + a_3*b_1*c_2*d_4*e_5 - a_3*b_1*c_2*d_5*e_4 - a_3*b_1*c_4*d_2*e_5 + a_3*b_1*c_4*d_5*e_2 + a_3*b_1*c_5*d_2*e_4 - a_3*b_1*c_5*d_4*e_2 - a_3*b_2*c_1*d_4*e_5 + a_3*b_2*c_1*d_5*e_4 + a_3*b_2*c_4*d_1*e_5 - a_3*b_2*c_4*d_5*e_1 - a_3*b_2*c_5*d_1*e_4 + a_3*b_2*c_5*d_4*e_1 &
                        + a_3*b_4*c_1*d_2*e_5 - a_3*b_4*c_1*d_5*e_2 - a_3*b_4*c_2*d_1*e_5 + a_3*b_4*c_2*d_5*e_1 + a_3*b_4*c_5*d_1*e_2 - a_3*b_4*c_5*d_2*e_1 - a_3*b_5*c_1*d_2*e_4 + a_3*b_5*c_1*d_4*e_2 + a_3*b_5*c_2*d_1*e_4 - a_3*b_5*c_2*d_4*e_1 - a_3*b_5*c_4*d_1*e_2 + a_3*b_5*c_4*d_2*e_1 - a_4*b_1*c_2*d_3*e_5 + a_4*b_1*c_2*d_5*e_3 + a_4*b_1*c_3*d_2*e_5 - a_4*b_1*c_3*d_5*e_2 - a_4*b_1*c_5*d_2*e_3 + a_4*b_1*c_5*d_3*e_2 + a_4*b_2*c_1*d_3*e_5 - a_4*b_2*c_1*d_5*e_3 &
                        - a_4*b_2*c_3*d_1*e_5 + a_4*b_2*c_3*d_5*e_1 + a_4*b_2*c_5*d_1*e_3 - a_4*b_2*c_5*d_3*e_1 - a_4*b_3*c_1*d_2*e_5 + a_4*b_3*c_1*d_5*e_2 + a_4*b_3*c_2*d_1*e_5 - a_4*b_3*c_2*d_5*e_1 - a_4*b_3*c_5*d_1*e_2 + a_4*b_3*c_5*d_2*e_1 + a_4*b_5*c_1*d_2*e_3 - a_4*b_5*c_1*d_3*e_2 - a_4*b_5*c_2*d_1*e_3 + a_4*b_5*c_2*d_3*e_1 + a_4*b_5*c_3*d_1*e_2 - a_4*b_5*c_3*d_2*e_1 + a_5*b_1*c_2*d_3*e_4 - a_5*b_1*c_2*d_4*e_3 - a_5*b_1*c_3*d_2*e_4 + a_5*b_1*c_3*d_4*e_2 &
                        + a_5*b_1*c_4*d_2*e_3 - a_5*b_1*c_4*d_3*e_2 - a_5*b_2*c_1*d_3*e_4 + a_5*b_2*c_1*d_4*e_3 + a_5*b_2*c_3*d_1*e_4 - a_5*b_2*c_3*d_4*e_1 - a_5*b_2*c_4*d_1*e_3 + a_5*b_2*c_4*d_3*e_1 + a_5*b_3*c_1*d_2*e_4 - a_5*b_3*c_1*d_4*e_2 - a_5*b_3*c_2*d_1*e_4 + a_5*b_3*c_2*d_4*e_1 + a_5*b_3*c_4*d_1*e_2 - a_5*b_3*c_4*d_2*e_1 - a_5*b_4*c_1*d_2*e_3 + a_5*b_4*c_1*d_3*e_2 + a_5*b_4*c_2*d_1*e_3 - a_5*b_4*c_2*d_3*e_1 - a_5*b_4*c_3*d_1*e_2 + a_5*b_4*c_3*d_2*e_1)
          
          if(martix_value/=0) then
  
              Inverse_A(1,1)=(b_2*c_3*d_4*e_5 - b_2*c_3*d_5*e_4 - b_2*c_4*d_3*e_5 + b_2*c_4*d_5*e_3 + b_2*c_5*d_3*e_4 - b_2*c_5*d_4*e_3 - b_3*c_2*d_4*e_5 + b_3*c_2*d_5*e_4 + b_3*c_4*d_2*e_5 - b_3*c_4*d_5*e_2 - b_3*c_5*d_2*e_4 + b_3*c_5*d_4*e_2 + b_4*c_2*d_3*e_5 - b_4*c_2*d_5*e_3 - b_4*c_3*d_2*e_5 + b_4*c_3*d_5*e_2 + b_4*c_5*d_2*e_3 - b_4*c_5*d_3*e_2 - b_5*c_2*d_3*e_4 + b_5*c_2*d_4*e_3 + b_5*c_3*d_2*e_4 - b_5*c_3*d_4*e_2 - b_5*c_4*d_2*e_3 + b_5*c_4*d_3*e_2)/martix_value
              Inverse_A(2,1)=-(a_2*c_3*d_4*e_5 - a_2*c_3*d_5*e_4 - a_2*c_4*d_3*e_5 + a_2*c_4*d_5*e_3 + a_2*c_5*d_3*e_4 - a_2*c_5*d_4*e_3 - a_3*c_2*d_4*e_5 + a_3*c_2*d_5*e_4 + a_3*c_4*d_2*e_5 - a_3*c_4*d_5*e_2 - a_3*c_5*d_2*e_4 + a_3*c_5*d_4*e_2 + a_4*c_2*d_3*e_5 - a_4*c_2*d_5*e_3 - a_4*c_3*d_2*e_5 + a_4*c_3*d_5*e_2 + a_4*c_5*d_2*e_3 - a_4*c_5*d_3*e_2 - a_5*c_2*d_3*e_4 + a_5*c_2*d_4*e_3 + a_5*c_3*d_2*e_4 - a_5*c_3*d_4*e_2 - a_5*c_4*d_2*e_3 + a_5*c_4*d_3*e_2)/martix_value
              Inverse_A(3,1)=(a_2*b_3*d_4*e_5 - a_2*b_3*d_5*e_4 - a_2*b_4*d_3*e_5 + a_2*b_4*d_5*e_3 + a_2*b_5*d_3*e_4 - a_2*b_5*d_4*e_3 - a_3*b_2*d_4*e_5 + a_3*b_2*d_5*e_4 + a_3*b_4*d_2*e_5 - a_3*b_4*d_5*e_2 - a_3*b_5*d_2*e_4 + a_3*b_5*d_4*e_2 + a_4*b_2*d_3*e_5 - a_4*b_2*d_5*e_3 - a_4*b_3*d_2*e_5 + a_4*b_3*d_5*e_2 + a_4*b_5*d_2*e_3 - a_4*b_5*d_3*e_2 - a_5*b_2*d_3*e_4 + a_5*b_2*d_4*e_3 + a_5*b_3*d_2*e_4 - a_5*b_3*d_4*e_2 - a_5*b_4*d_2*e_3 + a_5*b_4*d_3*e_2)/martix_value
              Inverse_A(4,1)=-(a_2*b_3*c_4*e_5 - a_2*b_3*c_5*e_4 - a_2*b_4*c_3*e_5 + a_2*b_4*c_5*e_3 + a_2*b_5*c_3*e_4 - a_2*b_5*c_4*e_3 - a_3*b_2*c_4*e_5 + a_3*b_2*c_5*e_4 + a_3*b_4*c_2*e_5 - a_3*b_4*c_5*e_2 - a_3*b_5*c_2*e_4 + a_3*b_5*c_4*e_2 + a_4*b_2*c_3*e_5 - a_4*b_2*c_5*e_3 - a_4*b_3*c_2*e_5 + a_4*b_3*c_5*e_2 + a_4*b_5*c_2*e_3 - a_4*b_5*c_3*e_2 - a_5*b_2*c_3*e_4 + a_5*b_2*c_4*e_3 + a_5*b_3*c_2*e_4 - a_5*b_3*c_4*e_2 - a_5*b_4*c_2*e_3 + a_5*b_4*c_3*e_2)/martix_value
              Inverse_A(5,1)=(a_2*b_3*c_4*d_5 - a_2*b_3*c_5*d_4 - a_2*b_4*c_3*d_5 + a_2*b_4*c_5*d_3 + a_2*b_5*c_3*d_4 - a_2*b_5*c_4*d_3 - a_3*b_2*c_4*d_5 + a_3*b_2*c_5*d_4 + a_3*b_4*c_2*d_5 - a_3*b_4*c_5*d_2 - a_3*b_5*c_2*d_4 + a_3*b_5*c_4*d_2 + a_4*b_2*c_3*d_5 - a_4*b_2*c_5*d_3 - a_4*b_3*c_2*d_5 + a_4*b_3*c_5*d_2 + a_4*b_5*c_2*d_3 - a_4*b_5*c_3*d_2 - a_5*b_2*c_3*d_4 + a_5*b_2*c_4*d_3 + a_5*b_3*c_2*d_4 - a_5*b_3*c_4*d_2 - a_5*b_4*c_2*d_3 + a_5*b_4*c_3*d_2)/martix_value
              
              Inverse_A(1,2)=-(b_1*c_3*d_4*e_5 - b_1*c_3*d_5*e_4 - b_1*c_4*d_3*e_5 + b_1*c_4*d_5*e_3 + b_1*c_5*d_3*e_4 - b_1*c_5*d_4*e_3 - b_3*c_1*d_4*e_5 + b_3*c_1*d_5*e_4 + b_3*c_4*d_1*e_5 - b_3*c_4*d_5*e_1 - b_3*c_5*d_1*e_4 + b_3*c_5*d_4*e_1 + b_4*c_1*d_3*e_5 - b_4*c_1*d_5*e_3 - b_4*c_3*d_1*e_5 + b_4*c_3*d_5*e_1 + b_4*c_5*d_1*e_3 - b_4*c_5*d_3*e_1 - b_5*c_1*d_3*e_4 + b_5*c_1*d_4*e_3 + b_5*c_3*d_1*e_4 - b_5*c_3*d_4*e_1 - b_5*c_4*d_1*e_3 + b_5*c_4*d_3*e_1)/martix_value
              Inverse_A(2,2)=(a_1*c_3*d_4*e_5 - a_1*c_3*d_5*e_4 - a_1*c_4*d_3*e_5 + a_1*c_4*d_5*e_3 + a_1*c_5*d_3*e_4 - a_1*c_5*d_4*e_3 - a_3*c_1*d_4*e_5 + a_3*c_1*d_5*e_4 + a_3*c_4*d_1*e_5 - a_3*c_4*d_5*e_1 - a_3*c_5*d_1*e_4 + a_3*c_5*d_4*e_1 + a_4*c_1*d_3*e_5 - a_4*c_1*d_5*e_3 - a_4*c_3*d_1*e_5 + a_4*c_3*d_5*e_1 + a_4*c_5*d_1*e_3 - a_4*c_5*d_3*e_1 - a_5*c_1*d_3*e_4 + a_5*c_1*d_4*e_3 + a_5*c_3*d_1*e_4 - a_5*c_3*d_4*e_1 - a_5*c_4*d_1*e_3 + a_5*c_4*d_3*e_1)/martix_value
              Inverse_A(3,2)=-(a_1*b_3*d_4*e_5 - a_1*b_3*d_5*e_4 - a_1*b_4*d_3*e_5 + a_1*b_4*d_5*e_3 + a_1*b_5*d_3*e_4 - a_1*b_5*d_4*e_3 - a_3*b_1*d_4*e_5 + a_3*b_1*d_5*e_4 + a_3*b_4*d_1*e_5 - a_3*b_4*d_5*e_1 - a_3*b_5*d_1*e_4 + a_3*b_5*d_4*e_1 + a_4*b_1*d_3*e_5 - a_4*b_1*d_5*e_3 - a_4*b_3*d_1*e_5 + a_4*b_3*d_5*e_1 + a_4*b_5*d_1*e_3 - a_4*b_5*d_3*e_1 - a_5*b_1*d_3*e_4 + a_5*b_1*d_4*e_3 + a_5*b_3*d_1*e_4 - a_5*b_3*d_4*e_1 - a_5*b_4*d_1*e_3 + a_5*b_4*d_3*e_1)/martix_value
              Inverse_A(4,2)=(a_1*b_3*c_4*e_5 - a_1*b_3*c_5*e_4 - a_1*b_4*c_3*e_5 + a_1*b_4*c_5*e_3 + a_1*b_5*c_3*e_4 - a_1*b_5*c_4*e_3 - a_3*b_1*c_4*e_5 + a_3*b_1*c_5*e_4 + a_3*b_4*c_1*e_5 - a_3*b_4*c_5*e_1 - a_3*b_5*c_1*e_4 + a_3*b_5*c_4*e_1 + a_4*b_1*c_3*e_5 - a_4*b_1*c_5*e_3 - a_4*b_3*c_1*e_5 + a_4*b_3*c_5*e_1 + a_4*b_5*c_1*e_3 - a_4*b_5*c_3*e_1 - a_5*b_1*c_3*e_4 + a_5*b_1*c_4*e_3 + a_5*b_3*c_1*e_4 - a_5*b_3*c_4*e_1 - a_5*b_4*c_1*e_3 + a_5*b_4*c_3*e_1)/martix_value
              Inverse_A(5,2)=-(a_1*b_3*c_4*d_5 - a_1*b_3*c_5*d_4 - a_1*b_4*c_3*d_5 + a_1*b_4*c_5*d_3 + a_1*b_5*c_3*d_4 - a_1*b_5*c_4*d_3 - a_3*b_1*c_4*d_5 + a_3*b_1*c_5*d_4 + a_3*b_4*c_1*d_5 - a_3*b_4*c_5*d_1 - a_3*b_5*c_1*d_4 + a_3*b_5*c_4*d_1 + a_4*b_1*c_3*d_5 - a_4*b_1*c_5*d_3 - a_4*b_3*c_1*d_5 + a_4*b_3*c_5*d_1 + a_4*b_5*c_1*d_3 - a_4*b_5*c_3*d_1 - a_5*b_1*c_3*d_4 + a_5*b_1*c_4*d_3 + a_5*b_3*c_1*d_4 - a_5*b_3*c_4*d_1 - a_5*b_4*c_1*d_3 + a_5*b_4*c_3*d_1)/martix_value
              
              Inverse_A(1,3)=(b_1*c_2*d_4*e_5 - b_1*c_2*d_5*e_4 - b_1*c_4*d_2*e_5 + b_1*c_4*d_5*e_2 + b_1*c_5*d_2*e_4 - b_1*c_5*d_4*e_2 - b_2*c_1*d_4*e_5 + b_2*c_1*d_5*e_4 + b_2*c_4*d_1*e_5 - b_2*c_4*d_5*e_1 - b_2*c_5*d_1*e_4 + b_2*c_5*d_4*e_1 + b_4*c_1*d_2*e_5 - b_4*c_1*d_5*e_2 - b_4*c_2*d_1*e_5 + b_4*c_2*d_5*e_1 + b_4*c_5*d_1*e_2 - b_4*c_5*d_2*e_1 - b_5*c_1*d_2*e_4 + b_5*c_1*d_4*e_2 + b_5*c_2*d_1*e_4 - b_5*c_2*d_4*e_1 - b_5*c_4*d_1*e_2 + b_5*c_4*d_2*e_1)/martix_value
              Inverse_A(2,3)=-(a_1*c_2*d_4*e_5 - a_1*c_2*d_5*e_4 - a_1*c_4*d_2*e_5 + a_1*c_4*d_5*e_2 + a_1*c_5*d_2*e_4 - a_1*c_5*d_4*e_2 - a_2*c_1*d_4*e_5 + a_2*c_1*d_5*e_4 + a_2*c_4*d_1*e_5 - a_2*c_4*d_5*e_1 - a_2*c_5*d_1*e_4 + a_2*c_5*d_4*e_1 + a_4*c_1*d_2*e_5 - a_4*c_1*d_5*e_2 - a_4*c_2*d_1*e_5 + a_4*c_2*d_5*e_1 + a_4*c_5*d_1*e_2 - a_4*c_5*d_2*e_1 - a_5*c_1*d_2*e_4 + a_5*c_1*d_4*e_2 + a_5*c_2*d_1*e_4 - a_5*c_2*d_4*e_1 - a_5*c_4*d_1*e_2 + a_5*c_4*d_2*e_1)/martix_value
              Inverse_A(3,3)=(a_1*b_2*d_4*e_5 - a_1*b_2*d_5*e_4 - a_1*b_4*d_2*e_5 + a_1*b_4*d_5*e_2 + a_1*b_5*d_2*e_4 - a_1*b_5*d_4*e_2 - a_2*b_1*d_4*e_5 + a_2*b_1*d_5*e_4 + a_2*b_4*d_1*e_5 - a_2*b_4*d_5*e_1 - a_2*b_5*d_1*e_4 + a_2*b_5*d_4*e_1 + a_4*b_1*d_2*e_5 - a_4*b_1*d_5*e_2 - a_4*b_2*d_1*e_5 + a_4*b_2*d_5*e_1 + a_4*b_5*d_1*e_2 - a_4*b_5*d_2*e_1 - a_5*b_1*d_2*e_4 + a_5*b_1*d_4*e_2 + a_5*b_2*d_1*e_4 - a_5*b_2*d_4*e_1 - a_5*b_4*d_1*e_2 + a_5*b_4*d_2*e_1)/martix_value
              Inverse_A(4,3)=-(a_1*b_2*c_4*e_5 - a_1*b_2*c_5*e_4 - a_1*b_4*c_2*e_5 + a_1*b_4*c_5*e_2 + a_1*b_5*c_2*e_4 - a_1*b_5*c_4*e_2 - a_2*b_1*c_4*e_5 + a_2*b_1*c_5*e_4 + a_2*b_4*c_1*e_5 - a_2*b_4*c_5*e_1 - a_2*b_5*c_1*e_4 + a_2*b_5*c_4*e_1 + a_4*b_1*c_2*e_5 - a_4*b_1*c_5*e_2 - a_4*b_2*c_1*e_5 + a_4*b_2*c_5*e_1 + a_4*b_5*c_1*e_2 - a_4*b_5*c_2*e_1 - a_5*b_1*c_2*e_4 + a_5*b_1*c_4*e_2 + a_5*b_2*c_1*e_4 - a_5*b_2*c_4*e_1 - a_5*b_4*c_1*e_2 + a_5*b_4*c_2*e_1)/martix_value
              Inverse_A(5,3)=(a_1*b_2*c_4*d_5 - a_1*b_2*c_5*d_4 - a_1*b_4*c_2*d_5 + a_1*b_4*c_5*d_2 + a_1*b_5*c_2*d_4 - a_1*b_5*c_4*d_2 - a_2*b_1*c_4*d_5 + a_2*b_1*c_5*d_4 + a_2*b_4*c_1*d_5 - a_2*b_4*c_5*d_1 - a_2*b_5*c_1*d_4 + a_2*b_5*c_4*d_1 + a_4*b_1*c_2*d_5 - a_4*b_1*c_5*d_2 - a_4*b_2*c_1*d_5 + a_4*b_2*c_5*d_1 + a_4*b_5*c_1*d_2 - a_4*b_5*c_2*d_1 - a_5*b_1*c_2*d_4 + a_5*b_1*c_4*d_2 + a_5*b_2*c_1*d_4 - a_5*b_2*c_4*d_1 - a_5*b_4*c_1*d_2 + a_5*b_4*c_2*d_1)/martix_value

              Inverse_A(1,4)=-(b_1*c_2*d_3*e_5 - b_1*c_2*d_5*e_3 - b_1*c_3*d_2*e_5 + b_1*c_3*d_5*e_2 + b_1*c_5*d_2*e_3 - b_1*c_5*d_3*e_2 - b_2*c_1*d_3*e_5 + b_2*c_1*d_5*e_3 + b_2*c_3*d_1*e_5 - b_2*c_3*d_5*e_1 - b_2*c_5*d_1*e_3 + b_2*c_5*d_3*e_1 + b_3*c_1*d_2*e_5 - b_3*c_1*d_5*e_2 - b_3*c_2*d_1*e_5 + b_3*c_2*d_5*e_1 + b_3*c_5*d_1*e_2 - b_3*c_5*d_2*e_1 - b_5*c_1*d_2*e_3 + b_5*c_1*d_3*e_2 + b_5*c_2*d_1*e_3 - b_5*c_2*d_3*e_1 - b_5*c_3*d_1*e_2 + b_5*c_3*d_2*e_1)/martix_value
              Inverse_A(2,4)=(a_1*c_2*d_3*e_5 - a_1*c_2*d_5*e_3 - a_1*c_3*d_2*e_5 + a_1*c_3*d_5*e_2 + a_1*c_5*d_2*e_3 - a_1*c_5*d_3*e_2 - a_2*c_1*d_3*e_5 + a_2*c_1*d_5*e_3 + a_2*c_3*d_1*e_5 - a_2*c_3*d_5*e_1 - a_2*c_5*d_1*e_3 + a_2*c_5*d_3*e_1 + a_3*c_1*d_2*e_5 - a_3*c_1*d_5*e_2 - a_3*c_2*d_1*e_5 + a_3*c_2*d_5*e_1 + a_3*c_5*d_1*e_2 - a_3*c_5*d_2*e_1 - a_5*c_1*d_2*e_3 + a_5*c_1*d_3*e_2 + a_5*c_2*d_1*e_3 - a_5*c_2*d_3*e_1 - a_5*c_3*d_1*e_2 + a_5*c_3*d_2*e_1)/martix_value
              Inverse_A(3,4)=-(a_1*b_2*d_3*e_5 - a_1*b_2*d_5*e_3 - a_1*b_3*d_2*e_5 + a_1*b_3*d_5*e_2 + a_1*b_5*d_2*e_3 - a_1*b_5*d_3*e_2 - a_2*b_1*d_3*e_5 + a_2*b_1*d_5*e_3 + a_2*b_3*d_1*e_5 - a_2*b_3*d_5*e_1 - a_2*b_5*d_1*e_3 + a_2*b_5*d_3*e_1 + a_3*b_1*d_2*e_5 - a_3*b_1*d_5*e_2 - a_3*b_2*d_1*e_5 + a_3*b_2*d_5*e_1 + a_3*b_5*d_1*e_2 - a_3*b_5*d_2*e_1 - a_5*b_1*d_2*e_3 + a_5*b_1*d_3*e_2 + a_5*b_2*d_1*e_3 - a_5*b_2*d_3*e_1 - a_5*b_3*d_1*e_2 + a_5*b_3*d_2*e_1)/martix_value
              Inverse_A(4,4)=(a_1*b_2*c_3*e_5 - a_1*b_2*c_5*e_3 - a_1*b_3*c_2*e_5 + a_1*b_3*c_5*e_2 + a_1*b_5*c_2*e_3 - a_1*b_5*c_3*e_2 - a_2*b_1*c_3*e_5 + a_2*b_1*c_5*e_3 + a_2*b_3*c_1*e_5 - a_2*b_3*c_5*e_1 - a_2*b_5*c_1*e_3 + a_2*b_5*c_3*e_1 + a_3*b_1*c_2*e_5 - a_3*b_1*c_5*e_2 - a_3*b_2*c_1*e_5 + a_3*b_2*c_5*e_1 + a_3*b_5*c_1*e_2 - a_3*b_5*c_2*e_1 - a_5*b_1*c_2*e_3 + a_5*b_1*c_3*e_2 + a_5*b_2*c_1*e_3 - a_5*b_2*c_3*e_1 - a_5*b_3*c_1*e_2 + a_5*b_3*c_2*e_1)/martix_value
              Inverse_A(5,4)=-(a_1*b_2*c_3*d_5 - a_1*b_2*c_5*d_3 - a_1*b_3*c_2*d_5 + a_1*b_3*c_5*d_2 + a_1*b_5*c_2*d_3 - a_1*b_5*c_3*d_2 - a_2*b_1*c_3*d_5 + a_2*b_1*c_5*d_3 + a_2*b_3*c_1*d_5 - a_2*b_3*c_5*d_1 - a_2*b_5*c_1*d_3 + a_2*b_5*c_3*d_1 + a_3*b_1*c_2*d_5 - a_3*b_1*c_5*d_2 - a_3*b_2*c_1*d_5 + a_3*b_2*c_5*d_1 + a_3*b_5*c_1*d_2 - a_3*b_5*c_2*d_1 - a_5*b_1*c_2*d_3 + a_5*b_1*c_3*d_2 + a_5*b_2*c_1*d_3 - a_5*b_2*c_3*d_1 - a_5*b_3*c_1*d_2 + a_5*b_3*c_2*d_1)/martix_value
              
              Inverse_A(1,5)=(b_1*c_2*d_3*e_4 - b_1*c_2*d_4*e_3 - b_1*c_3*d_2*e_4 + b_1*c_3*d_4*e_2 + b_1*c_4*d_2*e_3 - b_1*c_4*d_3*e_2 - b_2*c_1*d_3*e_4 + b_2*c_1*d_4*e_3 + b_2*c_3*d_1*e_4 - b_2*c_3*d_4*e_1 - b_2*c_4*d_1*e_3 + b_2*c_4*d_3*e_1 + b_3*c_1*d_2*e_4 - b_3*c_1*d_4*e_2 - b_3*c_2*d_1*e_4 + b_3*c_2*d_4*e_1 + b_3*c_4*d_1*e_2 - b_3*c_4*d_2*e_1 - b_4*c_1*d_2*e_3 + b_4*c_1*d_3*e_2 + b_4*c_2*d_1*e_3 - b_4*c_2*d_3*e_1 - b_4*c_3*d_1*e_2 + b_4*c_3*d_2*e_1)/martix_value
              Inverse_A(2,5)=-(a_1*c_2*d_3*e_4 - a_1*c_2*d_4*e_3 - a_1*c_3*d_2*e_4 + a_1*c_3*d_4*e_2 + a_1*c_4*d_2*e_3 - a_1*c_4*d_3*e_2 - a_2*c_1*d_3*e_4 + a_2*c_1*d_4*e_3 + a_2*c_3*d_1*e_4 - a_2*c_3*d_4*e_1 - a_2*c_4*d_1*e_3 + a_2*c_4*d_3*e_1 + a_3*c_1*d_2*e_4 - a_3*c_1*d_4*e_2 - a_3*c_2*d_1*e_4 + a_3*c_2*d_4*e_1 + a_3*c_4*d_1*e_2 - a_3*c_4*d_2*e_1 - a_4*c_1*d_2*e_3 + a_4*c_1*d_3*e_2 + a_4*c_2*d_1*e_3 - a_4*c_2*d_3*e_1 - a_4*c_3*d_1*e_2 + a_4*c_3*d_2*e_1)/martix_value
              Inverse_A(3,5)=(a_1*b_2*d_3*e_4 - a_1*b_2*d_4*e_3 - a_1*b_3*d_2*e_4 + a_1*b_3*d_4*e_2 + a_1*b_4*d_2*e_3 - a_1*b_4*d_3*e_2 - a_2*b_1*d_3*e_4 + a_2*b_1*d_4*e_3 + a_2*b_3*d_1*e_4 - a_2*b_3*d_4*e_1 - a_2*b_4*d_1*e_3 + a_2*b_4*d_3*e_1 + a_3*b_1*d_2*e_4 - a_3*b_1*d_4*e_2 - a_3*b_2*d_1*e_4 + a_3*b_2*d_4*e_1 + a_3*b_4*d_1*e_2 - a_3*b_4*d_2*e_1 - a_4*b_1*d_2*e_3 + a_4*b_1*d_3*e_2 + a_4*b_2*d_1*e_3 - a_4*b_2*d_3*e_1 - a_4*b_3*d_1*e_2 + a_4*b_3*d_2*e_1)/martix_value
              Inverse_A(4,5)=-(a_1*b_2*c_3*e_4 - a_1*b_2*c_4*e_3 - a_1*b_3*c_2*e_4 + a_1*b_3*c_4*e_2 + a_1*b_4*c_2*e_3 - a_1*b_4*c_3*e_2 - a_2*b_1*c_3*e_4 + a_2*b_1*c_4*e_3 + a_2*b_3*c_1*e_4 - a_2*b_3*c_4*e_1 - a_2*b_4*c_1*e_3 + a_2*b_4*c_3*e_1 + a_3*b_1*c_2*e_4 - a_3*b_1*c_4*e_2 - a_3*b_2*c_1*e_4 + a_3*b_2*c_4*e_1 + a_3*b_4*c_1*e_2 - a_3*b_4*c_2*e_1 - a_4*b_1*c_2*e_3 + a_4*b_1*c_3*e_2 + a_4*b_2*c_1*e_3 - a_4*b_2*c_3*e_1 - a_4*b_3*c_1*e_2 + a_4*b_3*c_2*e_1)/martix_value
              Inverse_A(5,5)=(a_1*b_2*c_3*d_4 - a_1*b_2*c_4*d_3 - a_1*b_3*c_2*d_4 + a_1*b_3*c_4*d_2 + a_1*b_4*c_2*d_3 - a_1*b_4*c_3*d_2 - a_2*b_1*c_3*d_4 + a_2*b_1*c_4*d_3 + a_2*b_3*c_1*d_4 - a_2*b_3*c_4*d_1 - a_2*b_4*c_1*d_3 + a_2*b_4*c_3*d_1 + a_3*b_1*c_2*d_4 - a_3*b_1*c_4*d_2 - a_3*b_2*c_1*d_4 + a_3*b_2*c_4*d_1 + a_3*b_4*c_1*d_2 - a_3*b_4*c_2*d_1 - a_4*b_1*c_2*d_3 + a_4*b_1*c_3*d_2 + a_4*b_2*c_1*d_3 - a_4*b_2*c_3*d_1 - a_4*b_3*c_1*d_2 + a_4*b_3*c_2*d_1)/martix_value
          
          end if
        !********************************************************************************************************

        !********************************************************************************************************
        ! 6D martix inversion
        else if(n==6) then
         
          a_1=A(1,1); b_1=A(1,2); c_1=A(1,3); d_1=A(1,4); e_1=A(1,5); f_1=A(1,6)
          a_2=A(2,1); b_2=A(2,2); c_2=A(2,3); d_2=A(2,4); e_2=A(2,5); f_2=A(2,6)
          a_3=A(3,1); b_3=A(3,2); c_3=A(3,3); d_3=A(3,4); e_3=A(3,5); f_3=A(3,6)
          a_4=A(4,1); b_4=A(4,2); c_4=A(4,3); d_4=A(4,4); e_4=A(4,5); f_4=A(4,6)
          a_5=A(5,1); b_5=A(5,2); c_5=A(5,3); d_5=A(5,4); e_5=A(5,5); f_5=A(5,6)
          a_6=A(6,1); b_6=A(6,2); c_6=A(6,3); d_6=A(6,4); e_6=A(6,5); f_6=A(6,6)
  
          martix_value=(a_1*b_2*c_3*d_4*e_5*f_6 - a_1*b_2*c_3*d_4*e_6*f_5 - a_1*b_2*c_3*d_5*e_4*f_6 + a_1*b_2*c_3*d_5*e_6*f_4 + a_1*b_2*c_3*d_6*e_4*f_5 - a_1*b_2*c_3*d_6*e_5*f_4 - a_1*b_2*c_4*d_3*e_5*f_6 + a_1*b_2*c_4*d_3*e_6*f_5 + a_1*b_2*c_4*d_5*e_3*f_6 - a_1*b_2*c_4*d_5*e_6*f_3 - a_1*b_2*c_4*d_6*e_3*f_5 + a_1*b_2*c_4*d_6*e_5*f_3 + a_1*b_2*c_5*d_3*e_4*f_6 - a_1*b_2*c_5*d_3*e_6*f_4 - a_1*b_2*c_5*d_4*e_3*f_6 + a_1*b_2*c_5*d_4*e_6*f_3 + a_1*b_2*c_5*d_6*e_3*f_4 - a_1*b_2*c_5*d_6*e_4*f_3 - a_1*b_2*c_6*d_3*e_4*f_5 + a_1*b_2*c_6*d_3*e_5*f_4 + a_1*b_2*c_6*d_4*e_3*f_5 - a_1*b_2*c_6*d_4*e_5*f_3 - a_1*b_2*c_6*d_5*e_3*f_4 + a_1*b_2*c_6*d_5*e_4*f_3 - a_1*b_3*c_2*d_4*e_5*f_6 + a_1*b_3*c_2*d_4*e_6*f_5 + a_1*b_3*c_2*d_5*e_4*f_6 - a_1*b_3*c_2*d_5*e_6*f_4 - a_1*b_3*c_2*d_6*e_4*f_5 + a_1*b_3*c_2*d_6*e_5*f_4 + a_1*b_3*c_4*d_2*e_5*f_6 - a_1*b_3*c_4*d_2*e_6*f_5 - a_1*b_3*c_4*d_5*e_2*f_6 + a_1*b_3*c_4*d_5*e_6*f_2 + a_1*b_3*c_4*d_6*e_2*f_5 - a_1*b_3*c_4*d_6*e_5*f_2 - a_1*b_3*c_5*d_2*e_4*f_6 + a_1*b_3*c_5*d_2*e_6*f_4 + a_1*b_3*c_5*d_4*e_2*f_6 - a_1*b_3*c_5*d_4*e_6*f_2 - a_1*b_3*c_5*d_6*e_2*f_4 + a_1*b_3*c_5*d_6*e_4*f_2 + a_1*b_3*c_6*d_2*e_4*f_5 - a_1*b_3*c_6*d_2*e_5*f_4 - a_1*b_3*c_6*d_4*e_2*f_5 + a_1*b_3*c_6*d_4*e_5*f_2 + a_1*b_3*c_6*d_5*e_2*f_4 - a_1*b_3*c_6*d_5*e_4*f_2 + a_1*b_4*c_2*d_3*e_5*f_6 - a_1*b_4*c_2*d_3*e_6*f_5 - a_1*b_4*c_2*d_5*e_3*f_6 + a_1*b_4*c_2*d_5*e_6*f_3 + a_1*b_4*c_2*d_6*e_3*f_5 - a_1*b_4*c_2*d_6*e_5*f_3 - a_1*b_4*c_3*d_2*e_5*f_6 + a_1*b_4*c_3*d_2*e_6*f_5 + a_1*b_4*c_3*d_5*e_2*f_6 - a_1*b_4*c_3*d_5*e_6*f_2 - a_1*b_4*c_3*d_6*e_2*f_5 + a_1*b_4*c_3*d_6*e_5*f_2 + a_1*b_4*c_5*d_2*e_3*f_6 - a_1*b_4*c_5*d_2*e_6*f_3 - a_1*b_4*c_5*d_3*e_2*f_6 + a_1*b_4*c_5*d_3*e_6*f_2 + a_1*b_4*c_5*d_6*e_2*f_3 - a_1*b_4*c_5*d_6*e_3*f_2 - a_1*b_4*c_6*d_2*e_3*f_5 + a_1*b_4*c_6*d_2*e_5*f_3 + a_1*b_4*c_6*d_3*e_2*f_5 - a_1*b_4*c_6*d_3*e_5*f_2 - a_1*b_4*c_6*d_5*e_2*f_3 + a_1*b_4*c_6*d_5*e_3*f_2 - a_1*b_5*c_2*d_3*e_4*f_6 + a_1*b_5*c_2*d_3*e_6*f_4 + a_1*b_5*c_2*d_4*e_3*f_6 - a_1*b_5*c_2*d_4*e_6*f_3 - a_1*b_5*c_2*d_6*e_3*f_4 + a_1*b_5*c_2*d_6*e_4*f_3 + a_1*b_5*c_3*d_2*e_4*f_6 - a_1*b_5*c_3*d_2*e_6*f_4 - a_1*b_5*c_3*d_4*e_2*f_6 + a_1*b_5*c_3*d_4*e_6*f_2 + a_1*b_5*c_3*d_6*e_2*f_4 - a_1*b_5*c_3*d_6*e_4*f_2 - a_1*b_5*c_4*d_2*e_3*f_6 + a_1*b_5*c_4*d_2*e_6*f_3 + a_1*b_5*c_4*d_3*e_2*f_6 - a_1*b_5*c_4*d_3*e_6*f_2 - a_1*b_5*c_4*d_6*e_2*f_3 + a_1*b_5*c_4*d_6*e_3*f_2 + a_1*b_5*c_6*d_2*e_3*f_4 - a_1*b_5*c_6*d_2*e_4*f_3 - a_1*b_5*c_6*d_3*e_2*f_4 + a_1*b_5*c_6*d_3*e_4*f_2 + a_1*b_5*c_6*d_4*e_2*f_3 - a_1*b_5*c_6*d_4*e_3*f_2 + a_1*b_6*c_2*d_3*e_4*f_5 - a_1*b_6*c_2*d_3*e_5*f_4 - a_1*b_6*c_2*d_4*e_3*f_5 + a_1*b_6*c_2*d_4*e_5*f_3 + a_1*b_6*c_2*d_5*e_3*f_4 - a_1*b_6*c_2*d_5*e_4*f_3 &
                        - a_1*b_6*c_3*d_2*e_4*f_5 + a_1*b_6*c_3*d_2*e_5*f_4 + a_1*b_6*c_3*d_4*e_2*f_5 - a_1*b_6*c_3*d_4*e_5*f_2 - a_1*b_6*c_3*d_5*e_2*f_4 + a_1*b_6*c_3*d_5*e_4*f_2 + a_1*b_6*c_4*d_2*e_3*f_5 - a_1*b_6*c_4*d_2*e_5*f_3 - a_1*b_6*c_4*d_3*e_2*f_5 + a_1*b_6*c_4*d_3*e_5*f_2 + a_1*b_6*c_4*d_5*e_2*f_3 - a_1*b_6*c_4*d_5*e_3*f_2 - a_1*b_6*c_5*d_2*e_3*f_4 + a_1*b_6*c_5*d_2*e_4*f_3 + a_1*b_6*c_5*d_3*e_2*f_4 - a_1*b_6*c_5*d_3*e_4*f_2 - a_1*b_6*c_5*d_4*e_2*f_3 + a_1*b_6*c_5*d_4*e_3*f_2 - a_2*b_1*c_3*d_4*e_5*f_6 + a_2*b_1*c_3*d_4*e_6*f_5 + a_2*b_1*c_3*d_5*e_4*f_6 - a_2*b_1*c_3*d_5*e_6*f_4 - a_2*b_1*c_3*d_6*e_4*f_5 + a_2*b_1*c_3*d_6*e_5*f_4 + a_2*b_1*c_4*d_3*e_5*f_6 - a_2*b_1*c_4*d_3*e_6*f_5 - a_2*b_1*c_4*d_5*e_3*f_6 + a_2*b_1*c_4*d_5*e_6*f_3 + a_2*b_1*c_4*d_6*e_3*f_5 - a_2*b_1*c_4*d_6*e_5*f_3 - a_2*b_1*c_5*d_3*e_4*f_6 + a_2*b_1*c_5*d_3*e_6*f_4 + a_2*b_1*c_5*d_4*e_3*f_6 - a_2*b_1*c_5*d_4*e_6*f_3 - a_2*b_1*c_5*d_6*e_3*f_4 + a_2*b_1*c_5*d_6*e_4*f_3 + a_2*b_1*c_6*d_3*e_4*f_5 - a_2*b_1*c_6*d_3*e_5*f_4 - a_2*b_1*c_6*d_4*e_3*f_5 + a_2*b_1*c_6*d_4*e_5*f_3 + a_2*b_1*c_6*d_5*e_3*f_4 - a_2*b_1*c_6*d_5*e_4*f_3 + a_2*b_3*c_1*d_4*e_5*f_6 - a_2*b_3*c_1*d_4*e_6*f_5 - a_2*b_3*c_1*d_5*e_4*f_6 + a_2*b_3*c_1*d_5*e_6*f_4 + a_2*b_3*c_1*d_6*e_4*f_5 - a_2*b_3*c_1*d_6*e_5*f_4 - a_2*b_3*c_4*d_1*e_5*f_6 + a_2*b_3*c_4*d_1*e_6*f_5 + a_2*b_3*c_4*d_5*e_1*f_6 - a_2*b_3*c_4*d_5*e_6*f_1 - a_2*b_3*c_4*d_6*e_1*f_5 + a_2*b_3*c_4*d_6*e_5*f_1 + a_2*b_3*c_5*d_1*e_4*f_6 - a_2*b_3*c_5*d_1*e_6*f_4 - a_2*b_3*c_5*d_4*e_1*f_6 + a_2*b_3*c_5*d_4*e_6*f_1 + a_2*b_3*c_5*d_6*e_1*f_4 - a_2*b_3*c_5*d_6*e_4*f_1 - a_2*b_3*c_6*d_1*e_4*f_5 + a_2*b_3*c_6*d_1*e_5*f_4 + a_2*b_3*c_6*d_4*e_1*f_5 - a_2*b_3*c_6*d_4*e_5*f_1 - a_2*b_3*c_6*d_5*e_1*f_4 + a_2*b_3*c_6*d_5*e_4*f_1 - a_2*b_4*c_1*d_3*e_5*f_6 + a_2*b_4*c_1*d_3*e_6*f_5 + a_2*b_4*c_1*d_5*e_3*f_6 - a_2*b_4*c_1*d_5*e_6*f_3 - a_2*b_4*c_1*d_6*e_3*f_5 + a_2*b_4*c_1*d_6*e_5*f_3 + a_2*b_4*c_3*d_1*e_5*f_6 - a_2*b_4*c_3*d_1*e_6*f_5 - a_2*b_4*c_3*d_5*e_1*f_6 + a_2*b_4*c_3*d_5*e_6*f_1 + a_2*b_4*c_3*d_6*e_1*f_5 - a_2*b_4*c_3*d_6*e_5*f_1 - a_2*b_4*c_5*d_1*e_3*f_6 + a_2*b_4*c_5*d_1*e_6*f_3 + a_2*b_4*c_5*d_3*e_1*f_6 - a_2*b_4*c_5*d_3*e_6*f_1 - a_2*b_4*c_5*d_6*e_1*f_3 + a_2*b_4*c_5*d_6*e_3*f_1 + a_2*b_4*c_6*d_1*e_3*f_5 - a_2*b_4*c_6*d_1*e_5*f_3 - a_2*b_4*c_6*d_3*e_1*f_5 + a_2*b_4*c_6*d_3*e_5*f_1 + a_2*b_4*c_6*d_5*e_1*f_3 - a_2*b_4*c_6*d_5*e_3*f_1 + a_2*b_5*c_1*d_3*e_4*f_6 - a_2*b_5*c_1*d_3*e_6*f_4 - a_2*b_5*c_1*d_4*e_3*f_6 + a_2*b_5*c_1*d_4*e_6*f_3 + a_2*b_5*c_1*d_6*e_3*f_4 - a_2*b_5*c_1*d_6*e_4*f_3 - a_2*b_5*c_3*d_1*e_4*f_6 + a_2*b_5*c_3*d_1*e_6*f_4 + a_2*b_5*c_3*d_4*e_1*f_6 - a_2*b_5*c_3*d_4*e_6*f_1 - a_2*b_5*c_3*d_6*e_1*f_4 + a_2*b_5*c_3*d_6*e_4*f_1 &
                        + a_2*b_5*c_4*d_1*e_3*f_6 - a_2*b_5*c_4*d_1*e_6*f_3 - a_2*b_5*c_4*d_3*e_1*f_6 + a_2*b_5*c_4*d_3*e_6*f_1 + a_2*b_5*c_4*d_6*e_1*f_3 - a_2*b_5*c_4*d_6*e_3*f_1 - a_2*b_5*c_6*d_1*e_3*f_4 + a_2*b_5*c_6*d_1*e_4*f_3 + a_2*b_5*c_6*d_3*e_1*f_4 - a_2*b_5*c_6*d_3*e_4*f_1 - a_2*b_5*c_6*d_4*e_1*f_3 + a_2*b_5*c_6*d_4*e_3*f_1 - a_2*b_6*c_1*d_3*e_4*f_5 + a_2*b_6*c_1*d_3*e_5*f_4 + a_2*b_6*c_1*d_4*e_3*f_5 - a_2*b_6*c_1*d_4*e_5*f_3 - a_2*b_6*c_1*d_5*e_3*f_4 + a_2*b_6*c_1*d_5*e_4*f_3 + a_2*b_6*c_3*d_1*e_4*f_5 - a_2*b_6*c_3*d_1*e_5*f_4 - a_2*b_6*c_3*d_4*e_1*f_5 + a_2*b_6*c_3*d_4*e_5*f_1 + a_2*b_6*c_3*d_5*e_1*f_4 - a_2*b_6*c_3*d_5*e_4*f_1 - a_2*b_6*c_4*d_1*e_3*f_5 + a_2*b_6*c_4*d_1*e_5*f_3 + a_2*b_6*c_4*d_3*e_1*f_5 - a_2*b_6*c_4*d_3*e_5*f_1 - a_2*b_6*c_4*d_5*e_1*f_3 + a_2*b_6*c_4*d_5*e_3*f_1 + a_2*b_6*c_5*d_1*e_3*f_4 - a_2*b_6*c_5*d_1*e_4*f_3 - a_2*b_6*c_5*d_3*e_1*f_4 + a_2*b_6*c_5*d_3*e_4*f_1 + a_2*b_6*c_5*d_4*e_1*f_3 - a_2*b_6*c_5*d_4*e_3*f_1 + a_3*b_1*c_2*d_4*e_5*f_6 - a_3*b_1*c_2*d_4*e_6*f_5 - a_3*b_1*c_2*d_5*e_4*f_6 + a_3*b_1*c_2*d_5*e_6*f_4 + a_3*b_1*c_2*d_6*e_4*f_5 - a_3*b_1*c_2*d_6*e_5*f_4 - a_3*b_1*c_4*d_2*e_5*f_6 + a_3*b_1*c_4*d_2*e_6*f_5 + a_3*b_1*c_4*d_5*e_2*f_6 - a_3*b_1*c_4*d_5*e_6*f_2 - a_3*b_1*c_4*d_6*e_2*f_5 + a_3*b_1*c_4*d_6*e_5*f_2 + a_3*b_1*c_5*d_2*e_4*f_6 - a_3*b_1*c_5*d_2*e_6*f_4 - a_3*b_1*c_5*d_4*e_2*f_6 + a_3*b_1*c_5*d_4*e_6*f_2 + a_3*b_1*c_5*d_6*e_2*f_4 - a_3*b_1*c_5*d_6*e_4*f_2 - a_3*b_1*c_6*d_2*e_4*f_5 + a_3*b_1*c_6*d_2*e_5*f_4 + a_3*b_1*c_6*d_4*e_2*f_5 - a_3*b_1*c_6*d_4*e_5*f_2 - a_3*b_1*c_6*d_5*e_2*f_4 + a_3*b_1*c_6*d_5*e_4*f_2 - a_3*b_2*c_1*d_4*e_5*f_6 + a_3*b_2*c_1*d_4*e_6*f_5 + a_3*b_2*c_1*d_5*e_4*f_6 - a_3*b_2*c_1*d_5*e_6*f_4 - a_3*b_2*c_1*d_6*e_4*f_5 + a_3*b_2*c_1*d_6*e_5*f_4 + a_3*b_2*c_4*d_1*e_5*f_6 - a_3*b_2*c_4*d_1*e_6*f_5 - a_3*b_2*c_4*d_5*e_1*f_6 + a_3*b_2*c_4*d_5*e_6*f_1 + a_3*b_2*c_4*d_6*e_1*f_5 - a_3*b_2*c_4*d_6*e_5*f_1 - a_3*b_2*c_5*d_1*e_4*f_6 + a_3*b_2*c_5*d_1*e_6*f_4 + a_3*b_2*c_5*d_4*e_1*f_6 - a_3*b_2*c_5*d_4*e_6*f_1 - a_3*b_2*c_5*d_6*e_1*f_4 + a_3*b_2*c_5*d_6*e_4*f_1 + a_3*b_2*c_6*d_1*e_4*f_5 - a_3*b_2*c_6*d_1*e_5*f_4 - a_3*b_2*c_6*d_4*e_1*f_5 + a_3*b_2*c_6*d_4*e_5*f_1 + a_3*b_2*c_6*d_5*e_1*f_4 - a_3*b_2*c_6*d_5*e_4*f_1 + a_3*b_4*c_1*d_2*e_5*f_6 - a_3*b_4*c_1*d_2*e_6*f_5 - a_3*b_4*c_1*d_5*e_2*f_6 + a_3*b_4*c_1*d_5*e_6*f_2 + a_3*b_4*c_1*d_6*e_2*f_5 - a_3*b_4*c_1*d_6*e_5*f_2 - a_3*b_4*c_2*d_1*e_5*f_6 + a_3*b_4*c_2*d_1*e_6*f_5 + a_3*b_4*c_2*d_5*e_1*f_6 - a_3*b_4*c_2*d_5*e_6*f_1 - a_3*b_4*c_2*d_6*e_1*f_5 + a_3*b_4*c_2*d_6*e_5*f_1 + a_3*b_4*c_5*d_1*e_2*f_6 - a_3*b_4*c_5*d_1*e_6*f_2 - a_3*b_4*c_5*d_2*e_1*f_6 + a_3*b_4*c_5*d_2*e_6*f_1 + a_3*b_4*c_5*d_6*e_1*f_2 - a_3*b_4*c_5*d_6*e_2*f_1 &
                        - a_3*b_4*c_6*d_1*e_2*f_5 + a_3*b_4*c_6*d_1*e_5*f_2 + a_3*b_4*c_6*d_2*e_1*f_5 - a_3*b_4*c_6*d_2*e_5*f_1 - a_3*b_4*c_6*d_5*e_1*f_2 + a_3*b_4*c_6*d_5*e_2*f_1 - a_3*b_5*c_1*d_2*e_4*f_6 + a_3*b_5*c_1*d_2*e_6*f_4 + a_3*b_5*c_1*d_4*e_2*f_6 - a_3*b_5*c_1*d_4*e_6*f_2 - a_3*b_5*c_1*d_6*e_2*f_4 + a_3*b_5*c_1*d_6*e_4*f_2 + a_3*b_5*c_2*d_1*e_4*f_6 - a_3*b_5*c_2*d_1*e_6*f_4 - a_3*b_5*c_2*d_4*e_1*f_6 + a_3*b_5*c_2*d_4*e_6*f_1 + a_3*b_5*c_2*d_6*e_1*f_4 - a_3*b_5*c_2*d_6*e_4*f_1 - a_3*b_5*c_4*d_1*e_2*f_6 + a_3*b_5*c_4*d_1*e_6*f_2 + a_3*b_5*c_4*d_2*e_1*f_6 - a_3*b_5*c_4*d_2*e_6*f_1 - a_3*b_5*c_4*d_6*e_1*f_2 + a_3*b_5*c_4*d_6*e_2*f_1 + a_3*b_5*c_6*d_1*e_2*f_4 - a_3*b_5*c_6*d_1*e_4*f_2 - a_3*b_5*c_6*d_2*e_1*f_4 + a_3*b_5*c_6*d_2*e_4*f_1 + a_3*b_5*c_6*d_4*e_1*f_2 - a_3*b_5*c_6*d_4*e_2*f_1 + a_3*b_6*c_1*d_2*e_4*f_5 - a_3*b_6*c_1*d_2*e_5*f_4 - a_3*b_6*c_1*d_4*e_2*f_5 + a_3*b_6*c_1*d_4*e_5*f_2 + a_3*b_6*c_1*d_5*e_2*f_4 - a_3*b_6*c_1*d_5*e_4*f_2 - a_3*b_6*c_2*d_1*e_4*f_5 + a_3*b_6*c_2*d_1*e_5*f_4 + a_3*b_6*c_2*d_4*e_1*f_5 - a_3*b_6*c_2*d_4*e_5*f_1 - a_3*b_6*c_2*d_5*e_1*f_4 + a_3*b_6*c_2*d_5*e_4*f_1 + a_3*b_6*c_4*d_1*e_2*f_5 - a_3*b_6*c_4*d_1*e_5*f_2 - a_3*b_6*c_4*d_2*e_1*f_5 + a_3*b_6*c_4*d_2*e_5*f_1 + a_3*b_6*c_4*d_5*e_1*f_2 - a_3*b_6*c_4*d_5*e_2*f_1 - a_3*b_6*c_5*d_1*e_2*f_4 + a_3*b_6*c_5*d_1*e_4*f_2 + a_3*b_6*c_5*d_2*e_1*f_4 - a_3*b_6*c_5*d_2*e_4*f_1 - a_3*b_6*c_5*d_4*e_1*f_2 + a_3*b_6*c_5*d_4*e_2*f_1 - a_4*b_1*c_2*d_3*e_5*f_6 + a_4*b_1*c_2*d_3*e_6*f_5 + a_4*b_1*c_2*d_5*e_3*f_6 - a_4*b_1*c_2*d_5*e_6*f_3 - a_4*b_1*c_2*d_6*e_3*f_5 + a_4*b_1*c_2*d_6*e_5*f_3 + a_4*b_1*c_3*d_2*e_5*f_6 - a_4*b_1*c_3*d_2*e_6*f_5 - a_4*b_1*c_3*d_5*e_2*f_6 + a_4*b_1*c_3*d_5*e_6*f_2 + a_4*b_1*c_3*d_6*e_2*f_5 - a_4*b_1*c_3*d_6*e_5*f_2 - a_4*b_1*c_5*d_2*e_3*f_6 + a_4*b_1*c_5*d_2*e_6*f_3 + a_4*b_1*c_5*d_3*e_2*f_6 - a_4*b_1*c_5*d_3*e_6*f_2 - a_4*b_1*c_5*d_6*e_2*f_3 + a_4*b_1*c_5*d_6*e_3*f_2 + a_4*b_1*c_6*d_2*e_3*f_5 - a_4*b_1*c_6*d_2*e_5*f_3 - a_4*b_1*c_6*d_3*e_2*f_5 + a_4*b_1*c_6*d_3*e_5*f_2 + a_4*b_1*c_6*d_5*e_2*f_3 - a_4*b_1*c_6*d_5*e_3*f_2 + a_4*b_2*c_1*d_3*e_5*f_6 - a_4*b_2*c_1*d_3*e_6*f_5 - a_4*b_2*c_1*d_5*e_3*f_6 + a_4*b_2*c_1*d_5*e_6*f_3 + a_4*b_2*c_1*d_6*e_3*f_5 - a_4*b_2*c_1*d_6*e_5*f_3 - a_4*b_2*c_3*d_1*e_5*f_6 + a_4*b_2*c_3*d_1*e_6*f_5 + a_4*b_2*c_3*d_5*e_1*f_6 - a_4*b_2*c_3*d_5*e_6*f_1 - a_4*b_2*c_3*d_6*e_1*f_5 + a_4*b_2*c_3*d_6*e_5*f_1 + a_4*b_2*c_5*d_1*e_3*f_6 - a_4*b_2*c_5*d_1*e_6*f_3 - a_4*b_2*c_5*d_3*e_1*f_6 + a_4*b_2*c_5*d_3*e_6*f_1 + a_4*b_2*c_5*d_6*e_1*f_3 - a_4*b_2*c_5*d_6*e_3*f_1 - a_4*b_2*c_6*d_1*e_3*f_5 + a_4*b_2*c_6*d_1*e_5*f_3 + a_4*b_2*c_6*d_3*e_1*f_5 - a_4*b_2*c_6*d_3*e_5*f_1 - a_4*b_2*c_6*d_5*e_1*f_3 + a_4*b_2*c_6*d_5*e_3*f_1 &
                        - a_4*b_3*c_1*d_2*e_5*f_6 + a_4*b_3*c_1*d_2*e_6*f_5 + a_4*b_3*c_1*d_5*e_2*f_6 - a_4*b_3*c_1*d_5*e_6*f_2 - a_4*b_3*c_1*d_6*e_2*f_5 + a_4*b_3*c_1*d_6*e_5*f_2 + a_4*b_3*c_2*d_1*e_5*f_6 - a_4*b_3*c_2*d_1*e_6*f_5 - a_4*b_3*c_2*d_5*e_1*f_6 + a_4*b_3*c_2*d_5*e_6*f_1 + a_4*b_3*c_2*d_6*e_1*f_5 - a_4*b_3*c_2*d_6*e_5*f_1 - a_4*b_3*c_5*d_1*e_2*f_6 + a_4*b_3*c_5*d_1*e_6*f_2 + a_4*b_3*c_5*d_2*e_1*f_6 - a_4*b_3*c_5*d_2*e_6*f_1 - a_4*b_3*c_5*d_6*e_1*f_2 + a_4*b_3*c_5*d_6*e_2*f_1 + a_4*b_3*c_6*d_1*e_2*f_5 - a_4*b_3*c_6*d_1*e_5*f_2 - a_4*b_3*c_6*d_2*e_1*f_5 + a_4*b_3*c_6*d_2*e_5*f_1 + a_4*b_3*c_6*d_5*e_1*f_2 - a_4*b_3*c_6*d_5*e_2*f_1 + a_4*b_5*c_1*d_2*e_3*f_6 - a_4*b_5*c_1*d_2*e_6*f_3 - a_4*b_5*c_1*d_3*e_2*f_6 + a_4*b_5*c_1*d_3*e_6*f_2 + a_4*b_5*c_1*d_6*e_2*f_3 - a_4*b_5*c_1*d_6*e_3*f_2 - a_4*b_5*c_2*d_1*e_3*f_6 + a_4*b_5*c_2*d_1*e_6*f_3 + a_4*b_5*c_2*d_3*e_1*f_6 - a_4*b_5*c_2*d_3*e_6*f_1 - a_4*b_5*c_2*d_6*e_1*f_3 + a_4*b_5*c_2*d_6*e_3*f_1 + a_4*b_5*c_3*d_1*e_2*f_6 - a_4*b_5*c_3*d_1*e_6*f_2 - a_4*b_5*c_3*d_2*e_1*f_6 + a_4*b_5*c_3*d_2*e_6*f_1 + a_4*b_5*c_3*d_6*e_1*f_2 - a_4*b_5*c_3*d_6*e_2*f_1 - a_4*b_5*c_6*d_1*e_2*f_3 + a_4*b_5*c_6*d_1*e_3*f_2 + a_4*b_5*c_6*d_2*e_1*f_3 - a_4*b_5*c_6*d_2*e_3*f_1 - a_4*b_5*c_6*d_3*e_1*f_2 + a_4*b_5*c_6*d_3*e_2*f_1 - a_4*b_6*c_1*d_2*e_3*f_5 + a_4*b_6*c_1*d_2*e_5*f_3 + a_4*b_6*c_1*d_3*e_2*f_5 - a_4*b_6*c_1*d_3*e_5*f_2 - a_4*b_6*c_1*d_5*e_2*f_3 + a_4*b_6*c_1*d_5*e_3*f_2 + a_4*b_6*c_2*d_1*e_3*f_5 - a_4*b_6*c_2*d_1*e_5*f_3 - a_4*b_6*c_2*d_3*e_1*f_5 + a_4*b_6*c_2*d_3*e_5*f_1 + a_4*b_6*c_2*d_5*e_1*f_3 - a_4*b_6*c_2*d_5*e_3*f_1 - a_4*b_6*c_3*d_1*e_2*f_5 + a_4*b_6*c_3*d_1*e_5*f_2 + a_4*b_6*c_3*d_2*e_1*f_5 - a_4*b_6*c_3*d_2*e_5*f_1 - a_4*b_6*c_3*d_5*e_1*f_2 + a_4*b_6*c_3*d_5*e_2*f_1 + a_4*b_6*c_5*d_1*e_2*f_3 - a_4*b_6*c_5*d_1*e_3*f_2 - a_4*b_6*c_5*d_2*e_1*f_3 + a_4*b_6*c_5*d_2*e_3*f_1 + a_4*b_6*c_5*d_3*e_1*f_2 - a_4*b_6*c_5*d_3*e_2*f_1 + a_5*b_1*c_2*d_3*e_4*f_6 - a_5*b_1*c_2*d_3*e_6*f_4 - a_5*b_1*c_2*d_4*e_3*f_6 + a_5*b_1*c_2*d_4*e_6*f_3 + a_5*b_1*c_2*d_6*e_3*f_4 - a_5*b_1*c_2*d_6*e_4*f_3 - a_5*b_1*c_3*d_2*e_4*f_6 + a_5*b_1*c_3*d_2*e_6*f_4 + a_5*b_1*c_3*d_4*e_2*f_6 - a_5*b_1*c_3*d_4*e_6*f_2 - a_5*b_1*c_3*d_6*e_2*f_4 + a_5*b_1*c_3*d_6*e_4*f_2 + a_5*b_1*c_4*d_2*e_3*f_6 - a_5*b_1*c_4*d_2*e_6*f_3 - a_5*b_1*c_4*d_3*e_2*f_6 + a_5*b_1*c_4*d_3*e_6*f_2 + a_5*b_1*c_4*d_6*e_2*f_3 - a_5*b_1*c_4*d_6*e_3*f_2 - a_5*b_1*c_6*d_2*e_3*f_4 + a_5*b_1*c_6*d_2*e_4*f_3 + a_5*b_1*c_6*d_3*e_2*f_4 - a_5*b_1*c_6*d_3*e_4*f_2 - a_5*b_1*c_6*d_4*e_2*f_3 + a_5*b_1*c_6*d_4*e_3*f_2 - a_5*b_2*c_1*d_3*e_4*f_6 + a_5*b_2*c_1*d_3*e_6*f_4 + a_5*b_2*c_1*d_4*e_3*f_6 - a_5*b_2*c_1*d_4*e_6*f_3 - a_5*b_2*c_1*d_6*e_3*f_4 + a_5*b_2*c_1*d_6*e_4*f_3 &
                        + a_5*b_2*c_3*d_1*e_4*f_6 - a_5*b_2*c_3*d_1*e_6*f_4 - a_5*b_2*c_3*d_4*e_1*f_6 + a_5*b_2*c_3*d_4*e_6*f_1 + a_5*b_2*c_3*d_6*e_1*f_4 - a_5*b_2*c_3*d_6*e_4*f_1 - a_5*b_2*c_4*d_1*e_3*f_6 + a_5*b_2*c_4*d_1*e_6*f_3 + a_5*b_2*c_4*d_3*e_1*f_6 - a_5*b_2*c_4*d_3*e_6*f_1 - a_5*b_2*c_4*d_6*e_1*f_3 + a_5*b_2*c_4*d_6*e_3*f_1 + a_5*b_2*c_6*d_1*e_3*f_4 - a_5*b_2*c_6*d_1*e_4*f_3 - a_5*b_2*c_6*d_3*e_1*f_4 + a_5*b_2*c_6*d_3*e_4*f_1 + a_5*b_2*c_6*d_4*e_1*f_3 - a_5*b_2*c_6*d_4*e_3*f_1 + a_5*b_3*c_1*d_2*e_4*f_6 - a_5*b_3*c_1*d_2*e_6*f_4 - a_5*b_3*c_1*d_4*e_2*f_6 + a_5*b_3*c_1*d_4*e_6*f_2 + a_5*b_3*c_1*d_6*e_2*f_4 - a_5*b_3*c_1*d_6*e_4*f_2 - a_5*b_3*c_2*d_1*e_4*f_6 + a_5*b_3*c_2*d_1*e_6*f_4 + a_5*b_3*c_2*d_4*e_1*f_6 - a_5*b_3*c_2*d_4*e_6*f_1 - a_5*b_3*c_2*d_6*e_1*f_4 + a_5*b_3*c_2*d_6*e_4*f_1 + a_5*b_3*c_4*d_1*e_2*f_6 - a_5*b_3*c_4*d_1*e_6*f_2 - a_5*b_3*c_4*d_2*e_1*f_6 + a_5*b_3*c_4*d_2*e_6*f_1 + a_5*b_3*c_4*d_6*e_1*f_2 - a_5*b_3*c_4*d_6*e_2*f_1 - a_5*b_3*c_6*d_1*e_2*f_4 + a_5*b_3*c_6*d_1*e_4*f_2 + a_5*b_3*c_6*d_2*e_1*f_4 - a_5*b_3*c_6*d_2*e_4*f_1 - a_5*b_3*c_6*d_4*e_1*f_2 + a_5*b_3*c_6*d_4*e_2*f_1 - a_5*b_4*c_1*d_2*e_3*f_6 + a_5*b_4*c_1*d_2*e_6*f_3 + a_5*b_4*c_1*d_3*e_2*f_6 - a_5*b_4*c_1*d_3*e_6*f_2 - a_5*b_4*c_1*d_6*e_2*f_3 + a_5*b_4*c_1*d_6*e_3*f_2 + a_5*b_4*c_2*d_1*e_3*f_6 - a_5*b_4*c_2*d_1*e_6*f_3 - a_5*b_4*c_2*d_3*e_1*f_6 + a_5*b_4*c_2*d_3*e_6*f_1 + a_5*b_4*c_2*d_6*e_1*f_3 - a_5*b_4*c_2*d_6*e_3*f_1 - a_5*b_4*c_3*d_1*e_2*f_6 + a_5*b_4*c_3*d_1*e_6*f_2 + a_5*b_4*c_3*d_2*e_1*f_6 - a_5*b_4*c_3*d_2*e_6*f_1 - a_5*b_4*c_3*d_6*e_1*f_2 + a_5*b_4*c_3*d_6*e_2*f_1 + a_5*b_4*c_6*d_1*e_2*f_3 - a_5*b_4*c_6*d_1*e_3*f_2 - a_5*b_4*c_6*d_2*e_1*f_3 + a_5*b_4*c_6*d_2*e_3*f_1 + a_5*b_4*c_6*d_3*e_1*f_2 - a_5*b_4*c_6*d_3*e_2*f_1 + a_5*b_6*c_1*d_2*e_3*f_4 - a_5*b_6*c_1*d_2*e_4*f_3 - a_5*b_6*c_1*d_3*e_2*f_4 + a_5*b_6*c_1*d_3*e_4*f_2 + a_5*b_6*c_1*d_4*e_2*f_3 - a_5*b_6*c_1*d_4*e_3*f_2 - a_5*b_6*c_2*d_1*e_3*f_4 + a_5*b_6*c_2*d_1*e_4*f_3 + a_5*b_6*c_2*d_3*e_1*f_4 - a_5*b_6*c_2*d_3*e_4*f_1 - a_5*b_6*c_2*d_4*e_1*f_3 + a_5*b_6*c_2*d_4*e_3*f_1 + a_5*b_6*c_3*d_1*e_2*f_4 - a_5*b_6*c_3*d_1*e_4*f_2 - a_5*b_6*c_3*d_2*e_1*f_4 + a_5*b_6*c_3*d_2*e_4*f_1 + a_5*b_6*c_3*d_4*e_1*f_2 - a_5*b_6*c_3*d_4*e_2*f_1 - a_5*b_6*c_4*d_1*e_2*f_3 + a_5*b_6*c_4*d_1*e_3*f_2 + a_5*b_6*c_4*d_2*e_1*f_3 - a_5*b_6*c_4*d_2*e_3*f_1 - a_5*b_6*c_4*d_3*e_1*f_2 + a_5*b_6*c_4*d_3*e_2*f_1 - a_6*b_1*c_2*d_3*e_4*f_5 + a_6*b_1*c_2*d_3*e_5*f_4 + a_6*b_1*c_2*d_4*e_3*f_5 - a_6*b_1*c_2*d_4*e_5*f_3 - a_6*b_1*c_2*d_5*e_3*f_4 + a_6*b_1*c_2*d_5*e_4*f_3 + a_6*b_1*c_3*d_2*e_4*f_5 - a_6*b_1*c_3*d_2*e_5*f_4 - a_6*b_1*c_3*d_4*e_2*f_5 + a_6*b_1*c_3*d_4*e_5*f_2 + a_6*b_1*c_3*d_5*e_2*f_4 - a_6*b_1*c_3*d_5*e_4*f_2 &
                        - a_6*b_1*c_4*d_2*e_3*f_5 + a_6*b_1*c_4*d_2*e_5*f_3 + a_6*b_1*c_4*d_3*e_2*f_5 - a_6*b_1*c_4*d_3*e_5*f_2 - a_6*b_1*c_4*d_5*e_2*f_3 + a_6*b_1*c_4*d_5*e_3*f_2 + a_6*b_1*c_5*d_2*e_3*f_4 - a_6*b_1*c_5*d_2*e_4*f_3 - a_6*b_1*c_5*d_3*e_2*f_4 + a_6*b_1*c_5*d_3*e_4*f_2 + a_6*b_1*c_5*d_4*e_2*f_3 - a_6*b_1*c_5*d_4*e_3*f_2 + a_6*b_2*c_1*d_3*e_4*f_5 - a_6*b_2*c_1*d_3*e_5*f_4 - a_6*b_2*c_1*d_4*e_3*f_5 + a_6*b_2*c_1*d_4*e_5*f_3 + a_6*b_2*c_1*d_5*e_3*f_4 - a_6*b_2*c_1*d_5*e_4*f_3 - a_6*b_2*c_3*d_1*e_4*f_5 + a_6*b_2*c_3*d_1*e_5*f_4 + a_6*b_2*c_3*d_4*e_1*f_5 - a_6*b_2*c_3*d_4*e_5*f_1 - a_6*b_2*c_3*d_5*e_1*f_4 + a_6*b_2*c_3*d_5*e_4*f_1 + a_6*b_2*c_4*d_1*e_3*f_5 - a_6*b_2*c_4*d_1*e_5*f_3 - a_6*b_2*c_4*d_3*e_1*f_5 + a_6*b_2*c_4*d_3*e_5*f_1 + a_6*b_2*c_4*d_5*e_1*f_3 - a_6*b_2*c_4*d_5*e_3*f_1 - a_6*b_2*c_5*d_1*e_3*f_4 + a_6*b_2*c_5*d_1*e_4*f_3 + a_6*b_2*c_5*d_3*e_1*f_4 - a_6*b_2*c_5*d_3*e_4*f_1 - a_6*b_2*c_5*d_4*e_1*f_3 + a_6*b_2*c_5*d_4*e_3*f_1 - a_6*b_3*c_1*d_2*e_4*f_5 + a_6*b_3*c_1*d_2*e_5*f_4 + a_6*b_3*c_1*d_4*e_2*f_5 - a_6*b_3*c_1*d_4*e_5*f_2 - a_6*b_3*c_1*d_5*e_2*f_4 + a_6*b_3*c_1*d_5*e_4*f_2 + a_6*b_3*c_2*d_1*e_4*f_5 - a_6*b_3*c_2*d_1*e_5*f_4 - a_6*b_3*c_2*d_4*e_1*f_5 + a_6*b_3*c_2*d_4*e_5*f_1 + a_6*b_3*c_2*d_5*e_1*f_4 - a_6*b_3*c_2*d_5*e_4*f_1 - a_6*b_3*c_4*d_1*e_2*f_5 + a_6*b_3*c_4*d_1*e_5*f_2 + a_6*b_3*c_4*d_2*e_1*f_5 - a_6*b_3*c_4*d_2*e_5*f_1 - a_6*b_3*c_4*d_5*e_1*f_2 + a_6*b_3*c_4*d_5*e_2*f_1 + a_6*b_3*c_5*d_1*e_2*f_4 - a_6*b_3*c_5*d_1*e_4*f_2 - a_6*b_3*c_5*d_2*e_1*f_4 + a_6*b_3*c_5*d_2*e_4*f_1 + a_6*b_3*c_5*d_4*e_1*f_2 - a_6*b_3*c_5*d_4*e_2*f_1 + a_6*b_4*c_1*d_2*e_3*f_5 - a_6*b_4*c_1*d_2*e_5*f_3 - a_6*b_4*c_1*d_3*e_2*f_5 + a_6*b_4*c_1*d_3*e_5*f_2 + a_6*b_4*c_1*d_5*e_2*f_3 - a_6*b_4*c_1*d_5*e_3*f_2 - a_6*b_4*c_2*d_1*e_3*f_5 + a_6*b_4*c_2*d_1*e_5*f_3 + a_6*b_4*c_2*d_3*e_1*f_5 - a_6*b_4*c_2*d_3*e_5*f_1 - a_6*b_4*c_2*d_5*e_1*f_3 + a_6*b_4*c_2*d_5*e_3*f_1 + a_6*b_4*c_3*d_1*e_2*f_5 - a_6*b_4*c_3*d_1*e_5*f_2 - a_6*b_4*c_3*d_2*e_1*f_5 + a_6*b_4*c_3*d_2*e_5*f_1 + a_6*b_4*c_3*d_5*e_1*f_2 - a_6*b_4*c_3*d_5*e_2*f_1 - a_6*b_4*c_5*d_1*e_2*f_3 + a_6*b_4*c_5*d_1*e_3*f_2 + a_6*b_4*c_5*d_2*e_1*f_3 - a_6*b_4*c_5*d_2*e_3*f_1 - a_6*b_4*c_5*d_3*e_1*f_2 + a_6*b_4*c_5*d_3*e_2*f_1 - a_6*b_5*c_1*d_2*e_3*f_4 + a_6*b_5*c_1*d_2*e_4*f_3 + a_6*b_5*c_1*d_3*e_2*f_4 - a_6*b_5*c_1*d_3*e_4*f_2 - a_6*b_5*c_1*d_4*e_2*f_3 + a_6*b_5*c_1*d_4*e_3*f_2 + a_6*b_5*c_2*d_1*e_3*f_4 - a_6*b_5*c_2*d_1*e_4*f_3 - a_6*b_5*c_2*d_3*e_1*f_4 + a_6*b_5*c_2*d_3*e_4*f_1 + a_6*b_5*c_2*d_4*e_1*f_3 - a_6*b_5*c_2*d_4*e_3*f_1 - a_6*b_5*c_3*d_1*e_2*f_4 + a_6*b_5*c_3*d_1*e_4*f_2 + a_6*b_5*c_3*d_2*e_1*f_4 - a_6*b_5*c_3*d_2*e_4*f_1 - a_6*b_5*c_3*d_4*e_1*f_2 + a_6*b_5*c_3*d_4*e_2*f_1 &
                        + a_6*b_5*c_4*d_1*e_2*f_3 - a_6*b_5*c_4*d_1*e_3*f_2 - a_6*b_5*c_4*d_2*e_1*f_3 + a_6*b_5*c_4*d_2*e_3*f_1 + a_6*b_5*c_4*d_3*e_1*f_2 - a_6*b_5*c_4*d_3*e_2*f_1)
          
          if(martix_value/=0) then
  
              Inverse_A(1,1)=(b_2*c_3*d_4*e_5*f_6 - b_2*c_3*d_4*e_6*f_5 - b_2*c_3*d_5*e_4*f_6 + b_2*c_3*d_5*e_6*f_4 + b_2*c_3*d_6*e_4*f_5 - b_2*c_3*d_6*e_5*f_4 - b_2*c_4*d_3*e_5*f_6 + b_2*c_4*d_3*e_6*f_5 + b_2*c_4*d_5*e_3*f_6 - b_2*c_4*d_5*e_6*f_3 - b_2*c_4*d_6*e_3*f_5 + b_2*c_4*d_6*e_5*f_3 + b_2*c_5*d_3*e_4*f_6 - b_2*c_5*d_3*e_6*f_4 - b_2*c_5*d_4*e_3*f_6 + b_2*c_5*d_4*e_6*f_3 + b_2*c_5*d_6*e_3*f_4 - b_2*c_5*d_6*e_4*f_3 - b_2*c_6*d_3*e_4*f_5 + b_2*c_6*d_3*e_5*f_4 + b_2*c_6*d_4*e_3*f_5 - b_2*c_6*d_4*e_5*f_3 - b_2*c_6*d_5*e_3*f_4 + b_2*c_6*d_5*e_4*f_3 - b_3*c_2*d_4*e_5*f_6 + b_3*c_2*d_4*e_6*f_5 + b_3*c_2*d_5*e_4*f_6 - b_3*c_2*d_5*e_6*f_4 - b_3*c_2*d_6*e_4*f_5 + b_3*c_2*d_6*e_5*f_4 + b_3*c_4*d_2*e_5*f_6 - b_3*c_4*d_2*e_6*f_5 - b_3*c_4*d_5*e_2*f_6 + b_3*c_4*d_5*e_6*f_2 + b_3*c_4*d_6*e_2*f_5 - b_3*c_4*d_6*e_5*f_2 - b_3*c_5*d_2*e_4*f_6 + b_3*c_5*d_2*e_6*f_4 + b_3*c_5*d_4*e_2*f_6 - b_3*c_5*d_4*e_6*f_2 - b_3*c_5*d_6*e_2*f_4 + b_3*c_5*d_6*e_4*f_2 + b_3*c_6*d_2*e_4*f_5 - b_3*c_6*d_2*e_5*f_4 - b_3*c_6*d_4*e_2*f_5 + b_3*c_6*d_4*e_5*f_2 + b_3*c_6*d_5*e_2*f_4 - b_3*c_6*d_5*e_4*f_2 + b_4*c_2*d_3*e_5*f_6 - b_4*c_2*d_3*e_6*f_5 - b_4*c_2*d_5*e_3*f_6 + b_4*c_2*d_5*e_6*f_3 + b_4*c_2*d_6*e_3*f_5 - b_4*c_2*d_6*e_5*f_3 - b_4*c_3*d_2*e_5*f_6 + b_4*c_3*d_2*e_6*f_5 + b_4*c_3*d_5*e_2*f_6 - b_4*c_3*d_5*e_6*f_2 - b_4*c_3*d_6*e_2*f_5 + b_4*c_3*d_6*e_5*f_2 + b_4*c_5*d_2*e_3*f_6 - b_4*c_5*d_2*e_6*f_3 - b_4*c_5*d_3*e_2*f_6 + b_4*c_5*d_3*e_6*f_2 + b_4*c_5*d_6*e_2*f_3 - b_4*c_5*d_6*e_3*f_2 - b_4*c_6*d_2*e_3*f_5 + b_4*c_6*d_2*e_5*f_3 + b_4*c_6*d_3*e_2*f_5 - b_4*c_6*d_3*e_5*f_2 - b_4*c_6*d_5*e_2*f_3 + b_4*c_6*d_5*e_3*f_2 - b_5*c_2*d_3*e_4*f_6 + b_5*c_2*d_3*e_6*f_4 + b_5*c_2*d_4*e_3*f_6 - b_5*c_2*d_4*e_6*f_3 - b_5*c_2*d_6*e_3*f_4 + b_5*c_2*d_6*e_4*f_3 + b_5*c_3*d_2*e_4*f_6 - b_5*c_3*d_2*e_6*f_4 - b_5*c_3*d_4*e_2*f_6 + b_5*c_3*d_4*e_6*f_2 + b_5*c_3*d_6*e_2*f_4 - b_5*c_3*d_6*e_4*f_2 - b_5*c_4*d_2*e_3*f_6 + b_5*c_4*d_2*e_6*f_3 + b_5*c_4*d_3*e_2*f_6 - b_5*c_4*d_3*e_6*f_2 - b_5*c_4*d_6*e_2*f_3 + b_5*c_4*d_6*e_3*f_2 + b_5*c_6*d_2*e_3*f_4 - b_5*c_6*d_2*e_4*f_3 - b_5*c_6*d_3*e_2*f_4 + b_5*c_6*d_3*e_4*f_2 + b_5*c_6*d_4*e_2*f_3 - b_5*c_6*d_4*e_3*f_2 + b_6*c_2*d_3*e_4*f_5 - b_6*c_2*d_3*e_5*f_4 - b_6*c_2*d_4*e_3*f_5 + b_6*c_2*d_4*e_5*f_3 + b_6*c_2*d_5*e_3*f_4 - b_6*c_2*d_5*e_4*f_3 - b_6*c_3*d_2*e_4*f_5 + b_6*c_3*d_2*e_5*f_4 + b_6*c_3*d_4*e_2*f_5 - b_6*c_3*d_4*e_5*f_2 - b_6*c_3*d_5*e_2*f_4 + b_6*c_3*d_5*e_4*f_2 + b_6*c_4*d_2*e_3*f_5 - b_6*c_4*d_2*e_5*f_3 - b_6*c_4*d_3*e_2*f_5 + b_6*c_4*d_3*e_5*f_2 + b_6*c_4*d_5*e_2*f_3 - b_6*c_4*d_5*e_3*f_2 - b_6*c_5*d_2*e_3*f_4 + b_6*c_5*d_2*e_4*f_3 + b_6*c_5*d_3*e_2*f_4 - b_6*c_5*d_3*e_4*f_2 - b_6*c_5*d_4*e_2*f_3 + b_6*c_5*d_4*e_3*f_2)/martix_value
              Inverse_A(2,1)=-(a_2*c_3*d_4*e_5*f_6 - a_2*c_3*d_4*e_6*f_5 - a_2*c_3*d_5*e_4*f_6 + a_2*c_3*d_5*e_6*f_4 + a_2*c_3*d_6*e_4*f_5 - a_2*c_3*d_6*e_5*f_4 - a_2*c_4*d_3*e_5*f_6 + a_2*c_4*d_3*e_6*f_5 + a_2*c_4*d_5*e_3*f_6 - a_2*c_4*d_5*e_6*f_3 - a_2*c_4*d_6*e_3*f_5 + a_2*c_4*d_6*e_5*f_3 + a_2*c_5*d_3*e_4*f_6 - a_2*c_5*d_3*e_6*f_4 - a_2*c_5*d_4*e_3*f_6 + a_2*c_5*d_4*e_6*f_3 + a_2*c_5*d_6*e_3*f_4 - a_2*c_5*d_6*e_4*f_3 - a_2*c_6*d_3*e_4*f_5 + a_2*c_6*d_3*e_5*f_4 + a_2*c_6*d_4*e_3*f_5 - a_2*c_6*d_4*e_5*f_3 - a_2*c_6*d_5*e_3*f_4 + a_2*c_6*d_5*e_4*f_3 - a_3*c_2*d_4*e_5*f_6 + a_3*c_2*d_4*e_6*f_5 + a_3*c_2*d_5*e_4*f_6 - a_3*c_2*d_5*e_6*f_4 - a_3*c_2*d_6*e_4*f_5 + a_3*c_2*d_6*e_5*f_4 + a_3*c_4*d_2*e_5*f_6 - a_3*c_4*d_2*e_6*f_5 - a_3*c_4*d_5*e_2*f_6 + a_3*c_4*d_5*e_6*f_2 + a_3*c_4*d_6*e_2*f_5 - a_3*c_4*d_6*e_5*f_2 - a_3*c_5*d_2*e_4*f_6 + a_3*c_5*d_2*e_6*f_4 + a_3*c_5*d_4*e_2*f_6 - a_3*c_5*d_4*e_6*f_2 - a_3*c_5*d_6*e_2*f_4 + a_3*c_5*d_6*e_4*f_2 + a_3*c_6*d_2*e_4*f_5 - a_3*c_6*d_2*e_5*f_4 - a_3*c_6*d_4*e_2*f_5 + a_3*c_6*d_4*e_5*f_2 + a_3*c_6*d_5*e_2*f_4 - a_3*c_6*d_5*e_4*f_2 + a_4*c_2*d_3*e_5*f_6 - a_4*c_2*d_3*e_6*f_5 - a_4*c_2*d_5*e_3*f_6 + a_4*c_2*d_5*e_6*f_3 + a_4*c_2*d_6*e_3*f_5 - a_4*c_2*d_6*e_5*f_3 - a_4*c_3*d_2*e_5*f_6 + a_4*c_3*d_2*e_6*f_5 + a_4*c_3*d_5*e_2*f_6 - a_4*c_3*d_5*e_6*f_2 - a_4*c_3*d_6*e_2*f_5 + a_4*c_3*d_6*e_5*f_2 + a_4*c_5*d_2*e_3*f_6 - a_4*c_5*d_2*e_6*f_3 - a_4*c_5*d_3*e_2*f_6 + a_4*c_5*d_3*e_6*f_2 + a_4*c_5*d_6*e_2*f_3 - a_4*c_5*d_6*e_3*f_2 - a_4*c_6*d_2*e_3*f_5 + a_4*c_6*d_2*e_5*f_3 + a_4*c_6*d_3*e_2*f_5 - a_4*c_6*d_3*e_5*f_2 - a_4*c_6*d_5*e_2*f_3 + a_4*c_6*d_5*e_3*f_2 - a_5*c_2*d_3*e_4*f_6 + a_5*c_2*d_3*e_6*f_4 + a_5*c_2*d_4*e_3*f_6 - a_5*c_2*d_4*e_6*f_3 - a_5*c_2*d_6*e_3*f_4 + a_5*c_2*d_6*e_4*f_3 + a_5*c_3*d_2*e_4*f_6 - a_5*c_3*d_2*e_6*f_4 - a_5*c_3*d_4*e_2*f_6 + a_5*c_3*d_4*e_6*f_2 + a_5*c_3*d_6*e_2*f_4 - a_5*c_3*d_6*e_4*f_2 - a_5*c_4*d_2*e_3*f_6 + a_5*c_4*d_2*e_6*f_3 + a_5*c_4*d_3*e_2*f_6 - a_5*c_4*d_3*e_6*f_2 - a_5*c_4*d_6*e_2*f_3 + a_5*c_4*d_6*e_3*f_2 + a_5*c_6*d_2*e_3*f_4 - a_5*c_6*d_2*e_4*f_3 - a_5*c_6*d_3*e_2*f_4 + a_5*c_6*d_3*e_4*f_2 + a_5*c_6*d_4*e_2*f_3 - a_5*c_6*d_4*e_3*f_2 + a_6*c_2*d_3*e_4*f_5 - a_6*c_2*d_3*e_5*f_4 - a_6*c_2*d_4*e_3*f_5 + a_6*c_2*d_4*e_5*f_3 + a_6*c_2*d_5*e_3*f_4 - a_6*c_2*d_5*e_4*f_3 - a_6*c_3*d_2*e_4*f_5 + a_6*c_3*d_2*e_5*f_4 + a_6*c_3*d_4*e_2*f_5 - a_6*c_3*d_4*e_5*f_2 - a_6*c_3*d_5*e_2*f_4 + a_6*c_3*d_5*e_4*f_2 + a_6*c_4*d_2*e_3*f_5 - a_6*c_4*d_2*e_5*f_3 - a_6*c_4*d_3*e_2*f_5 + a_6*c_4*d_3*e_5*f_2 + a_6*c_4*d_5*e_2*f_3 - a_6*c_4*d_5*e_3*f_2 - a_6*c_5*d_2*e_3*f_4 + a_6*c_5*d_2*e_4*f_3 + a_6*c_5*d_3*e_2*f_4 - a_6*c_5*d_3*e_4*f_2 - a_6*c_5*d_4*e_2*f_3 + a_6*c_5*d_4*e_3*f_2)/martix_value
              Inverse_A(3,1)=(a_2*b_3*d_4*e_5*f_6 - a_2*b_3*d_4*e_6*f_5 - a_2*b_3*d_5*e_4*f_6 + a_2*b_3*d_5*e_6*f_4 + a_2*b_3*d_6*e_4*f_5 - a_2*b_3*d_6*e_5*f_4 - a_2*b_4*d_3*e_5*f_6 + a_2*b_4*d_3*e_6*f_5 + a_2*b_4*d_5*e_3*f_6 - a_2*b_4*d_5*e_6*f_3 - a_2*b_4*d_6*e_3*f_5 + a_2*b_4*d_6*e_5*f_3 + a_2*b_5*d_3*e_4*f_6 - a_2*b_5*d_3*e_6*f_4 - a_2*b_5*d_4*e_3*f_6 + a_2*b_5*d_4*e_6*f_3 + a_2*b_5*d_6*e_3*f_4 - a_2*b_5*d_6*e_4*f_3 - a_2*b_6*d_3*e_4*f_5 + a_2*b_6*d_3*e_5*f_4 + a_2*b_6*d_4*e_3*f_5 - a_2*b_6*d_4*e_5*f_3 - a_2*b_6*d_5*e_3*f_4 + a_2*b_6*d_5*e_4*f_3 - a_3*b_2*d_4*e_5*f_6 + a_3*b_2*d_4*e_6*f_5 + a_3*b_2*d_5*e_4*f_6 - a_3*b_2*d_5*e_6*f_4 - a_3*b_2*d_6*e_4*f_5 + a_3*b_2*d_6*e_5*f_4 + a_3*b_4*d_2*e_5*f_6 - a_3*b_4*d_2*e_6*f_5 - a_3*b_4*d_5*e_2*f_6 + a_3*b_4*d_5*e_6*f_2 + a_3*b_4*d_6*e_2*f_5 - a_3*b_4*d_6*e_5*f_2 - a_3*b_5*d_2*e_4*f_6 + a_3*b_5*d_2*e_6*f_4 + a_3*b_5*d_4*e_2*f_6 - a_3*b_5*d_4*e_6*f_2 - a_3*b_5*d_6*e_2*f_4 + a_3*b_5*d_6*e_4*f_2 + a_3*b_6*d_2*e_4*f_5 - a_3*b_6*d_2*e_5*f_4 - a_3*b_6*d_4*e_2*f_5 + a_3*b_6*d_4*e_5*f_2 + a_3*b_6*d_5*e_2*f_4 - a_3*b_6*d_5*e_4*f_2 + a_4*b_2*d_3*e_5*f_6 - a_4*b_2*d_3*e_6*f_5 - a_4*b_2*d_5*e_3*f_6 + a_4*b_2*d_5*e_6*f_3 + a_4*b_2*d_6*e_3*f_5 - a_4*b_2*d_6*e_5*f_3 - a_4*b_3*d_2*e_5*f_6 + a_4*b_3*d_2*e_6*f_5 + a_4*b_3*d_5*e_2*f_6 - a_4*b_3*d_5*e_6*f_2 - a_4*b_3*d_6*e_2*f_5 + a_4*b_3*d_6*e_5*f_2 + a_4*b_5*d_2*e_3*f_6 - a_4*b_5*d_2*e_6*f_3 - a_4*b_5*d_3*e_2*f_6 + a_4*b_5*d_3*e_6*f_2 + a_4*b_5*d_6*e_2*f_3 - a_4*b_5*d_6*e_3*f_2 - a_4*b_6*d_2*e_3*f_5 + a_4*b_6*d_2*e_5*f_3 + a_4*b_6*d_3*e_2*f_5 - a_4*b_6*d_3*e_5*f_2 - a_4*b_6*d_5*e_2*f_3 + a_4*b_6*d_5*e_3*f_2 - a_5*b_2*d_3*e_4*f_6 + a_5*b_2*d_3*e_6*f_4 + a_5*b_2*d_4*e_3*f_6 - a_5*b_2*d_4*e_6*f_3 - a_5*b_2*d_6*e_3*f_4 + a_5*b_2*d_6*e_4*f_3 + a_5*b_3*d_2*e_4*f_6 - a_5*b_3*d_2*e_6*f_4 - a_5*b_3*d_4*e_2*f_6 + a_5*b_3*d_4*e_6*f_2 + a_5*b_3*d_6*e_2*f_4 - a_5*b_3*d_6*e_4*f_2 - a_5*b_4*d_2*e_3*f_6 + a_5*b_4*d_2*e_6*f_3 + a_5*b_4*d_3*e_2*f_6 - a_5*b_4*d_3*e_6*f_2 - a_5*b_4*d_6*e_2*f_3 + a_5*b_4*d_6*e_3*f_2 + a_5*b_6*d_2*e_3*f_4 - a_5*b_6*d_2*e_4*f_3 - a_5*b_6*d_3*e_2*f_4 + a_5*b_6*d_3*e_4*f_2 + a_5*b_6*d_4*e_2*f_3 - a_5*b_6*d_4*e_3*f_2 + a_6*b_2*d_3*e_4*f_5 - a_6*b_2*d_3*e_5*f_4 - a_6*b_2*d_4*e_3*f_5 + a_6*b_2*d_4*e_5*f_3 + a_6*b_2*d_5*e_3*f_4 - a_6*b_2*d_5*e_4*f_3 - a_6*b_3*d_2*e_4*f_5 + a_6*b_3*d_2*e_5*f_4 + a_6*b_3*d_4*e_2*f_5 - a_6*b_3*d_4*e_5*f_2 - a_6*b_3*d_5*e_2*f_4 + a_6*b_3*d_5*e_4*f_2 + a_6*b_4*d_2*e_3*f_5 - a_6*b_4*d_2*e_5*f_3 - a_6*b_4*d_3*e_2*f_5 + a_6*b_4*d_3*e_5*f_2 + a_6*b_4*d_5*e_2*f_3 - a_6*b_4*d_5*e_3*f_2 - a_6*b_5*d_2*e_3*f_4 + a_6*b_5*d_2*e_4*f_3 + a_6*b_5*d_3*e_2*f_4 - a_6*b_5*d_3*e_4*f_2 - a_6*b_5*d_4*e_2*f_3 + a_6*b_5*d_4*e_3*f_2)/martix_value
              Inverse_A(4,1)=-(a_2*b_3*c_4*e_5*f_6 - a_2*b_3*c_4*e_6*f_5 - a_2*b_3*c_5*e_4*f_6 + a_2*b_3*c_5*e_6*f_4 + a_2*b_3*c_6*e_4*f_5 - a_2*b_3*c_6*e_5*f_4 - a_2*b_4*c_3*e_5*f_6 + a_2*b_4*c_3*e_6*f_5 + a_2*b_4*c_5*e_3*f_6 - a_2*b_4*c_5*e_6*f_3 - a_2*b_4*c_6*e_3*f_5 + a_2*b_4*c_6*e_5*f_3 + a_2*b_5*c_3*e_4*f_6 - a_2*b_5*c_3*e_6*f_4 - a_2*b_5*c_4*e_3*f_6 + a_2*b_5*c_4*e_6*f_3 + a_2*b_5*c_6*e_3*f_4 - a_2*b_5*c_6*e_4*f_3 - a_2*b_6*c_3*e_4*f_5 + a_2*b_6*c_3*e_5*f_4 + a_2*b_6*c_4*e_3*f_5 - a_2*b_6*c_4*e_5*f_3 - a_2*b_6*c_5*e_3*f_4 + a_2*b_6*c_5*e_4*f_3 - a_3*b_2*c_4*e_5*f_6 + a_3*b_2*c_4*e_6*f_5 + a_3*b_2*c_5*e_4*f_6 - a_3*b_2*c_5*e_6*f_4 - a_3*b_2*c_6*e_4*f_5 + a_3*b_2*c_6*e_5*f_4 + a_3*b_4*c_2*e_5*f_6 - a_3*b_4*c_2*e_6*f_5 - a_3*b_4*c_5*e_2*f_6 + a_3*b_4*c_5*e_6*f_2 + a_3*b_4*c_6*e_2*f_5 - a_3*b_4*c_6*e_5*f_2 - a_3*b_5*c_2*e_4*f_6 + a_3*b_5*c_2*e_6*f_4 + a_3*b_5*c_4*e_2*f_6 - a_3*b_5*c_4*e_6*f_2 - a_3*b_5*c_6*e_2*f_4 + a_3*b_5*c_6*e_4*f_2 + a_3*b_6*c_2*e_4*f_5 - a_3*b_6*c_2*e_5*f_4 - a_3*b_6*c_4*e_2*f_5 + a_3*b_6*c_4*e_5*f_2 + a_3*b_6*c_5*e_2*f_4 - a_3*b_6*c_5*e_4*f_2 + a_4*b_2*c_3*e_5*f_6 - a_4*b_2*c_3*e_6*f_5 - a_4*b_2*c_5*e_3*f_6 + a_4*b_2*c_5*e_6*f_3 + a_4*b_2*c_6*e_3*f_5 - a_4*b_2*c_6*e_5*f_3 - a_4*b_3*c_2*e_5*f_6 + a_4*b_3*c_2*e_6*f_5 + a_4*b_3*c_5*e_2*f_6 - a_4*b_3*c_5*e_6*f_2 - a_4*b_3*c_6*e_2*f_5 + a_4*b_3*c_6*e_5*f_2 + a_4*b_5*c_2*e_3*f_6 - a_4*b_5*c_2*e_6*f_3 - a_4*b_5*c_3*e_2*f_6 + a_4*b_5*c_3*e_6*f_2 + a_4*b_5*c_6*e_2*f_3 - a_4*b_5*c_6*e_3*f_2 - a_4*b_6*c_2*e_3*f_5 + a_4*b_6*c_2*e_5*f_3 + a_4*b_6*c_3*e_2*f_5 - a_4*b_6*c_3*e_5*f_2 - a_4*b_6*c_5*e_2*f_3 + a_4*b_6*c_5*e_3*f_2 - a_5*b_2*c_3*e_4*f_6 + a_5*b_2*c_3*e_6*f_4 + a_5*b_2*c_4*e_3*f_6 - a_5*b_2*c_4*e_6*f_3 - a_5*b_2*c_6*e_3*f_4 + a_5*b_2*c_6*e_4*f_3 + a_5*b_3*c_2*e_4*f_6 - a_5*b_3*c_2*e_6*f_4 - a_5*b_3*c_4*e_2*f_6 + a_5*b_3*c_4*e_6*f_2 + a_5*b_3*c_6*e_2*f_4 - a_5*b_3*c_6*e_4*f_2 - a_5*b_4*c_2*e_3*f_6 + a_5*b_4*c_2*e_6*f_3 + a_5*b_4*c_3*e_2*f_6 - a_5*b_4*c_3*e_6*f_2 - a_5*b_4*c_6*e_2*f_3 + a_5*b_4*c_6*e_3*f_2 + a_5*b_6*c_2*e_3*f_4 - a_5*b_6*c_2*e_4*f_3 - a_5*b_6*c_3*e_2*f_4 + a_5*b_6*c_3*e_4*f_2 + a_5*b_6*c_4*e_2*f_3 - a_5*b_6*c_4*e_3*f_2 + a_6*b_2*c_3*e_4*f_5 - a_6*b_2*c_3*e_5*f_4 - a_6*b_2*c_4*e_3*f_5 + a_6*b_2*c_4*e_5*f_3 + a_6*b_2*c_5*e_3*f_4 - a_6*b_2*c_5*e_4*f_3 - a_6*b_3*c_2*e_4*f_5 + a_6*b_3*c_2*e_5*f_4 + a_6*b_3*c_4*e_2*f_5 - a_6*b_3*c_4*e_5*f_2 - a_6*b_3*c_5*e_2*f_4 + a_6*b_3*c_5*e_4*f_2 + a_6*b_4*c_2*e_3*f_5 - a_6*b_4*c_2*e_5*f_3 - a_6*b_4*c_3*e_2*f_5 + a_6*b_4*c_3*e_5*f_2 + a_6*b_4*c_5*e_2*f_3 - a_6*b_4*c_5*e_3*f_2 - a_6*b_5*c_2*e_3*f_4 + a_6*b_5*c_2*e_4*f_3 + a_6*b_5*c_3*e_2*f_4 - a_6*b_5*c_3*e_4*f_2 - a_6*b_5*c_4*e_2*f_3 + a_6*b_5*c_4*e_3*f_2)/martix_value
              Inverse_A(5,1)=(a_2*b_3*c_4*d_5*f_6 - a_2*b_3*c_4*d_6*f_5 - a_2*b_3*c_5*d_4*f_6 + a_2*b_3*c_5*d_6*f_4 + a_2*b_3*c_6*d_4*f_5 - a_2*b_3*c_6*d_5*f_4 - a_2*b_4*c_3*d_5*f_6 + a_2*b_4*c_3*d_6*f_5 + a_2*b_4*c_5*d_3*f_6 - a_2*b_4*c_5*d_6*f_3 - a_2*b_4*c_6*d_3*f_5 + a_2*b_4*c_6*d_5*f_3 + a_2*b_5*c_3*d_4*f_6 - a_2*b_5*c_3*d_6*f_4 - a_2*b_5*c_4*d_3*f_6 + a_2*b_5*c_4*d_6*f_3 + a_2*b_5*c_6*d_3*f_4 - a_2*b_5*c_6*d_4*f_3 - a_2*b_6*c_3*d_4*f_5 + a_2*b_6*c_3*d_5*f_4 + a_2*b_6*c_4*d_3*f_5 - a_2*b_6*c_4*d_5*f_3 - a_2*b_6*c_5*d_3*f_4 + a_2*b_6*c_5*d_4*f_3 - a_3*b_2*c_4*d_5*f_6 + a_3*b_2*c_4*d_6*f_5 + a_3*b_2*c_5*d_4*f_6 - a_3*b_2*c_5*d_6*f_4 - a_3*b_2*c_6*d_4*f_5 + a_3*b_2*c_6*d_5*f_4 + a_3*b_4*c_2*d_5*f_6 - a_3*b_4*c_2*d_6*f_5 - a_3*b_4*c_5*d_2*f_6 + a_3*b_4*c_5*d_6*f_2 + a_3*b_4*c_6*d_2*f_5 - a_3*b_4*c_6*d_5*f_2 - a_3*b_5*c_2*d_4*f_6 + a_3*b_5*c_2*d_6*f_4 + a_3*b_5*c_4*d_2*f_6 - a_3*b_5*c_4*d_6*f_2 - a_3*b_5*c_6*d_2*f_4 + a_3*b_5*c_6*d_4*f_2 + a_3*b_6*c_2*d_4*f_5 - a_3*b_6*c_2*d_5*f_4 - a_3*b_6*c_4*d_2*f_5 + a_3*b_6*c_4*d_5*f_2 + a_3*b_6*c_5*d_2*f_4 - a_3*b_6*c_5*d_4*f_2 + a_4*b_2*c_3*d_5*f_6 - a_4*b_2*c_3*d_6*f_5 - a_4*b_2*c_5*d_3*f_6 + a_4*b_2*c_5*d_6*f_3 + a_4*b_2*c_6*d_3*f_5 - a_4*b_2*c_6*d_5*f_3 - a_4*b_3*c_2*d_5*f_6 + a_4*b_3*c_2*d_6*f_5 + a_4*b_3*c_5*d_2*f_6 - a_4*b_3*c_5*d_6*f_2 - a_4*b_3*c_6*d_2*f_5 + a_4*b_3*c_6*d_5*f_2 + a_4*b_5*c_2*d_3*f_6 - a_4*b_5*c_2*d_6*f_3 - a_4*b_5*c_3*d_2*f_6 + a_4*b_5*c_3*d_6*f_2 + a_4*b_5*c_6*d_2*f_3 - a_4*b_5*c_6*d_3*f_2 - a_4*b_6*c_2*d_3*f_5 + a_4*b_6*c_2*d_5*f_3 + a_4*b_6*c_3*d_2*f_5 - a_4*b_6*c_3*d_5*f_2 - a_4*b_6*c_5*d_2*f_3 + a_4*b_6*c_5*d_3*f_2 - a_5*b_2*c_3*d_4*f_6 + a_5*b_2*c_3*d_6*f_4 + a_5*b_2*c_4*d_3*f_6 - a_5*b_2*c_4*d_6*f_3 - a_5*b_2*c_6*d_3*f_4 + a_5*b_2*c_6*d_4*f_3 + a_5*b_3*c_2*d_4*f_6 - a_5*b_3*c_2*d_6*f_4 - a_5*b_3*c_4*d_2*f_6 + a_5*b_3*c_4*d_6*f_2 + a_5*b_3*c_6*d_2*f_4 - a_5*b_3*c_6*d_4*f_2 - a_5*b_4*c_2*d_3*f_6 + a_5*b_4*c_2*d_6*f_3 + a_5*b_4*c_3*d_2*f_6 - a_5*b_4*c_3*d_6*f_2 - a_5*b_4*c_6*d_2*f_3 + a_5*b_4*c_6*d_3*f_2 + a_5*b_6*c_2*d_3*f_4 - a_5*b_6*c_2*d_4*f_3 - a_5*b_6*c_3*d_2*f_4 + a_5*b_6*c_3*d_4*f_2 + a_5*b_6*c_4*d_2*f_3 - a_5*b_6*c_4*d_3*f_2 + a_6*b_2*c_3*d_4*f_5 - a_6*b_2*c_3*d_5*f_4 - a_6*b_2*c_4*d_3*f_5 + a_6*b_2*c_4*d_5*f_3 + a_6*b_2*c_5*d_3*f_4 - a_6*b_2*c_5*d_4*f_3 - a_6*b_3*c_2*d_4*f_5 + a_6*b_3*c_2*d_5*f_4 + a_6*b_3*c_4*d_2*f_5 - a_6*b_3*c_4*d_5*f_2 - a_6*b_3*c_5*d_2*f_4 + a_6*b_3*c_5*d_4*f_2 + a_6*b_4*c_2*d_3*f_5 - a_6*b_4*c_2*d_5*f_3 - a_6*b_4*c_3*d_2*f_5 + a_6*b_4*c_3*d_5*f_2 + a_6*b_4*c_5*d_2*f_3 - a_6*b_4*c_5*d_3*f_2 - a_6*b_5*c_2*d_3*f_4 + a_6*b_5*c_2*d_4*f_3 + a_6*b_5*c_3*d_2*f_4 - a_6*b_5*c_3*d_4*f_2 - a_6*b_5*c_4*d_2*f_3 + a_6*b_5*c_4*d_3*f_2)/martix_value
              Inverse_A(6,1)=-(a_2*b_3*c_4*d_5*e_6 - a_2*b_3*c_4*d_6*e_5 - a_2*b_3*c_5*d_4*e_6 + a_2*b_3*c_5*d_6*e_4 + a_2*b_3*c_6*d_4*e_5 - a_2*b_3*c_6*d_5*e_4 - a_2*b_4*c_3*d_5*e_6 + a_2*b_4*c_3*d_6*e_5 + a_2*b_4*c_5*d_3*e_6 - a_2*b_4*c_5*d_6*e_3 - a_2*b_4*c_6*d_3*e_5 + a_2*b_4*c_6*d_5*e_3 + a_2*b_5*c_3*d_4*e_6 - a_2*b_5*c_3*d_6*e_4 - a_2*b_5*c_4*d_3*e_6 + a_2*b_5*c_4*d_6*e_3 + a_2*b_5*c_6*d_3*e_4 - a_2*b_5*c_6*d_4*e_3 - a_2*b_6*c_3*d_4*e_5 + a_2*b_6*c_3*d_5*e_4 + a_2*b_6*c_4*d_3*e_5 - a_2*b_6*c_4*d_5*e_3 - a_2*b_6*c_5*d_3*e_4 + a_2*b_6*c_5*d_4*e_3 - a_3*b_2*c_4*d_5*e_6 + a_3*b_2*c_4*d_6*e_5 + a_3*b_2*c_5*d_4*e_6 - a_3*b_2*c_5*d_6*e_4 - a_3*b_2*c_6*d_4*e_5 + a_3*b_2*c_6*d_5*e_4 + a_3*b_4*c_2*d_5*e_6 - a_3*b_4*c_2*d_6*e_5 - a_3*b_4*c_5*d_2*e_6 + a_3*b_4*c_5*d_6*e_2 + a_3*b_4*c_6*d_2*e_5 - a_3*b_4*c_6*d_5*e_2 - a_3*b_5*c_2*d_4*e_6 + a_3*b_5*c_2*d_6*e_4 + a_3*b_5*c_4*d_2*e_6 - a_3*b_5*c_4*d_6*e_2 - a_3*b_5*c_6*d_2*e_4 + a_3*b_5*c_6*d_4*e_2 + a_3*b_6*c_2*d_4*e_5 - a_3*b_6*c_2*d_5*e_4 - a_3*b_6*c_4*d_2*e_5 + a_3*b_6*c_4*d_5*e_2 + a_3*b_6*c_5*d_2*e_4 - a_3*b_6*c_5*d_4*e_2 + a_4*b_2*c_3*d_5*e_6 - a_4*b_2*c_3*d_6*e_5 - a_4*b_2*c_5*d_3*e_6 + a_4*b_2*c_5*d_6*e_3 + a_4*b_2*c_6*d_3*e_5 - a_4*b_2*c_6*d_5*e_3 - a_4*b_3*c_2*d_5*e_6 + a_4*b_3*c_2*d_6*e_5 + a_4*b_3*c_5*d_2*e_6 - a_4*b_3*c_5*d_6*e_2 - a_4*b_3*c_6*d_2*e_5 + a_4*b_3*c_6*d_5*e_2 + a_4*b_5*c_2*d_3*e_6 - a_4*b_5*c_2*d_6*e_3 - a_4*b_5*c_3*d_2*e_6 + a_4*b_5*c_3*d_6*e_2 + a_4*b_5*c_6*d_2*e_3 - a_4*b_5*c_6*d_3*e_2 - a_4*b_6*c_2*d_3*e_5 + a_4*b_6*c_2*d_5*e_3 + a_4*b_6*c_3*d_2*e_5 - a_4*b_6*c_3*d_5*e_2 - a_4*b_6*c_5*d_2*e_3 + a_4*b_6*c_5*d_3*e_2 - a_5*b_2*c_3*d_4*e_6 + a_5*b_2*c_3*d_6*e_4 + a_5*b_2*c_4*d_3*e_6 - a_5*b_2*c_4*d_6*e_3 - a_5*b_2*c_6*d_3*e_4 + a_5*b_2*c_6*d_4*e_3 + a_5*b_3*c_2*d_4*e_6 - a_5*b_3*c_2*d_6*e_4 - a_5*b_3*c_4*d_2*e_6 + a_5*b_3*c_4*d_6*e_2 + a_5*b_3*c_6*d_2*e_4 - a_5*b_3*c_6*d_4*e_2 - a_5*b_4*c_2*d_3*e_6 + a_5*b_4*c_2*d_6*e_3 + a_5*b_4*c_3*d_2*e_6 - a_5*b_4*c_3*d_6*e_2 - a_5*b_4*c_6*d_2*e_3 + a_5*b_4*c_6*d_3*e_2 + a_5*b_6*c_2*d_3*e_4 - a_5*b_6*c_2*d_4*e_3 - a_5*b_6*c_3*d_2*e_4 + a_5*b_6*c_3*d_4*e_2 + a_5*b_6*c_4*d_2*e_3 - a_5*b_6*c_4*d_3*e_2 + a_6*b_2*c_3*d_4*e_5 - a_6*b_2*c_3*d_5*e_4 - a_6*b_2*c_4*d_3*e_5 + a_6*b_2*c_4*d_5*e_3 + a_6*b_2*c_5*d_3*e_4 - a_6*b_2*c_5*d_4*e_3 - a_6*b_3*c_2*d_4*e_5 + a_6*b_3*c_2*d_5*e_4 + a_6*b_3*c_4*d_2*e_5 - a_6*b_3*c_4*d_5*e_2 - a_6*b_3*c_5*d_2*e_4 + a_6*b_3*c_5*d_4*e_2 + a_6*b_4*c_2*d_3*e_5 - a_6*b_4*c_2*d_5*e_3 - a_6*b_4*c_3*d_2*e_5 + a_6*b_4*c_3*d_5*e_2 + a_6*b_4*c_5*d_2*e_3 - a_6*b_4*c_5*d_3*e_2 - a_6*b_5*c_2*d_3*e_4 + a_6*b_5*c_2*d_4*e_3 + a_6*b_5*c_3*d_2*e_4 - a_6*b_5*c_3*d_4*e_2 - a_6*b_5*c_4*d_2*e_3 + a_6*b_5*c_4*d_3*e_2)/martix_value
              
              Inverse_A(1,2)=-(b_1*c_3*d_4*e_5*f_6 - b_1*c_3*d_4*e_6*f_5 - b_1*c_3*d_5*e_4*f_6 + b_1*c_3*d_5*e_6*f_4 + b_1*c_3*d_6*e_4*f_5 - b_1*c_3*d_6*e_5*f_4 - b_1*c_4*d_3*e_5*f_6 + b_1*c_4*d_3*e_6*f_5 + b_1*c_4*d_5*e_3*f_6 - b_1*c_4*d_5*e_6*f_3 - b_1*c_4*d_6*e_3*f_5 + b_1*c_4*d_6*e_5*f_3 + b_1*c_5*d_3*e_4*f_6 - b_1*c_5*d_3*e_6*f_4 - b_1*c_5*d_4*e_3*f_6 + b_1*c_5*d_4*e_6*f_3 + b_1*c_5*d_6*e_3*f_4 - b_1*c_5*d_6*e_4*f_3 - b_1*c_6*d_3*e_4*f_5 + b_1*c_6*d_3*e_5*f_4 + b_1*c_6*d_4*e_3*f_5 - b_1*c_6*d_4*e_5*f_3 - b_1*c_6*d_5*e_3*f_4 + b_1*c_6*d_5*e_4*f_3 - b_3*c_1*d_4*e_5*f_6 + b_3*c_1*d_4*e_6*f_5 + b_3*c_1*d_5*e_4*f_6 - b_3*c_1*d_5*e_6*f_4 - b_3*c_1*d_6*e_4*f_5 + b_3*c_1*d_6*e_5*f_4 + b_3*c_4*d_1*e_5*f_6 - b_3*c_4*d_1*e_6*f_5 - b_3*c_4*d_5*e_1*f_6 + b_3*c_4*d_5*e_6*f_1 + b_3*c_4*d_6*e_1*f_5 - b_3*c_4*d_6*e_5*f_1 - b_3*c_5*d_1*e_4*f_6 + b_3*c_5*d_1*e_6*f_4 + b_3*c_5*d_4*e_1*f_6 - b_3*c_5*d_4*e_6*f_1 - b_3*c_5*d_6*e_1*f_4 + b_3*c_5*d_6*e_4*f_1 + b_3*c_6*d_1*e_4*f_5 - b_3*c_6*d_1*e_5*f_4 - b_3*c_6*d_4*e_1*f_5 + b_3*c_6*d_4*e_5*f_1 + b_3*c_6*d_5*e_1*f_4 - b_3*c_6*d_5*e_4*f_1 + b_4*c_1*d_3*e_5*f_6 - b_4*c_1*d_3*e_6*f_5 - b_4*c_1*d_5*e_3*f_6 + b_4*c_1*d_5*e_6*f_3 + b_4*c_1*d_6*e_3*f_5 - b_4*c_1*d_6*e_5*f_3 - b_4*c_3*d_1*e_5*f_6 + b_4*c_3*d_1*e_6*f_5 + b_4*c_3*d_5*e_1*f_6 - b_4*c_3*d_5*e_6*f_1 - b_4*c_3*d_6*e_1*f_5 + b_4*c_3*d_6*e_5*f_1 + b_4*c_5*d_1*e_3*f_6 - b_4*c_5*d_1*e_6*f_3 - b_4*c_5*d_3*e_1*f_6 + b_4*c_5*d_3*e_6*f_1 + b_4*c_5*d_6*e_1*f_3 - b_4*c_5*d_6*e_3*f_1 - b_4*c_6*d_1*e_3*f_5 + b_4*c_6*d_1*e_5*f_3 + b_4*c_6*d_3*e_1*f_5 - b_4*c_6*d_3*e_5*f_1 - b_4*c_6*d_5*e_1*f_3 + b_4*c_6*d_5*e_3*f_1 - b_5*c_1*d_3*e_4*f_6 + b_5*c_1*d_3*e_6*f_4 + b_5*c_1*d_4*e_3*f_6 - b_5*c_1*d_4*e_6*f_3 - b_5*c_1*d_6*e_3*f_4 + b_5*c_1*d_6*e_4*f_3 + b_5*c_3*d_1*e_4*f_6 - b_5*c_3*d_1*e_6*f_4 - b_5*c_3*d_4*e_1*f_6 + b_5*c_3*d_4*e_6*f_1 + b_5*c_3*d_6*e_1*f_4 - b_5*c_3*d_6*e_4*f_1 - b_5*c_4*d_1*e_3*f_6 + b_5*c_4*d_1*e_6*f_3 + b_5*c_4*d_3*e_1*f_6 - b_5*c_4*d_3*e_6*f_1 - b_5*c_4*d_6*e_1*f_3 + b_5*c_4*d_6*e_3*f_1 + b_5*c_6*d_1*e_3*f_4 - b_5*c_6*d_1*e_4*f_3 - b_5*c_6*d_3*e_1*f_4 + b_5*c_6*d_3*e_4*f_1 + b_5*c_6*d_4*e_1*f_3 - b_5*c_6*d_4*e_3*f_1 + b_6*c_1*d_3*e_4*f_5 - b_6*c_1*d_3*e_5*f_4 - b_6*c_1*d_4*e_3*f_5 + b_6*c_1*d_4*e_5*f_3 + b_6*c_1*d_5*e_3*f_4 - b_6*c_1*d_5*e_4*f_3 - b_6*c_3*d_1*e_4*f_5 + b_6*c_3*d_1*e_5*f_4 + b_6*c_3*d_4*e_1*f_5 - b_6*c_3*d_4*e_5*f_1 - b_6*c_3*d_5*e_1*f_4 + b_6*c_3*d_5*e_4*f_1 + b_6*c_4*d_1*e_3*f_5 - b_6*c_4*d_1*e_5*f_3 - b_6*c_4*d_3*e_1*f_5 + b_6*c_4*d_3*e_5*f_1 + b_6*c_4*d_5*e_1*f_3 - b_6*c_4*d_5*e_3*f_1 - b_6*c_5*d_1*e_3*f_4 + b_6*c_5*d_1*e_4*f_3 + b_6*c_5*d_3*e_1*f_4 - b_6*c_5*d_3*e_4*f_1 - b_6*c_5*d_4*e_1*f_3 + b_6*c_5*d_4*e_3*f_1)/martix_value
              Inverse_A(2,2)=(a_1*c_3*d_4*e_5*f_6 - a_1*c_3*d_4*e_6*f_5 - a_1*c_3*d_5*e_4*f_6 + a_1*c_3*d_5*e_6*f_4 + a_1*c_3*d_6*e_4*f_5 - a_1*c_3*d_6*e_5*f_4 - a_1*c_4*d_3*e_5*f_6 + a_1*c_4*d_3*e_6*f_5 + a_1*c_4*d_5*e_3*f_6 - a_1*c_4*d_5*e_6*f_3 - a_1*c_4*d_6*e_3*f_5 + a_1*c_4*d_6*e_5*f_3 + a_1*c_5*d_3*e_4*f_6 - a_1*c_5*d_3*e_6*f_4 - a_1*c_5*d_4*e_3*f_6 + a_1*c_5*d_4*e_6*f_3 + a_1*c_5*d_6*e_3*f_4 - a_1*c_5*d_6*e_4*f_3 - a_1*c_6*d_3*e_4*f_5 + a_1*c_6*d_3*e_5*f_4 + a_1*c_6*d_4*e_3*f_5 - a_1*c_6*d_4*e_5*f_3 - a_1*c_6*d_5*e_3*f_4 + a_1*c_6*d_5*e_4*f_3 - a_3*c_1*d_4*e_5*f_6 + a_3*c_1*d_4*e_6*f_5 + a_3*c_1*d_5*e_4*f_6 - a_3*c_1*d_5*e_6*f_4 - a_3*c_1*d_6*e_4*f_5 + a_3*c_1*d_6*e_5*f_4 + a_3*c_4*d_1*e_5*f_6 - a_3*c_4*d_1*e_6*f_5 - a_3*c_4*d_5*e_1*f_6 + a_3*c_4*d_5*e_6*f_1 + a_3*c_4*d_6*e_1*f_5 - a_3*c_4*d_6*e_5*f_1 - a_3*c_5*d_1*e_4*f_6 + a_3*c_5*d_1*e_6*f_4 + a_3*c_5*d_4*e_1*f_6 - a_3*c_5*d_4*e_6*f_1 - a_3*c_5*d_6*e_1*f_4 + a_3*c_5*d_6*e_4*f_1 + a_3*c_6*d_1*e_4*f_5 - a_3*c_6*d_1*e_5*f_4 - a_3*c_6*d_4*e_1*f_5 + a_3*c_6*d_4*e_5*f_1 + a_3*c_6*d_5*e_1*f_4 - a_3*c_6*d_5*e_4*f_1 + a_4*c_1*d_3*e_5*f_6 - a_4*c_1*d_3*e_6*f_5 - a_4*c_1*d_5*e_3*f_6 + a_4*c_1*d_5*e_6*f_3 + a_4*c_1*d_6*e_3*f_5 - a_4*c_1*d_6*e_5*f_3 - a_4*c_3*d_1*e_5*f_6 + a_4*c_3*d_1*e_6*f_5 + a_4*c_3*d_5*e_1*f_6 - a_4*c_3*d_5*e_6*f_1 - a_4*c_3*d_6*e_1*f_5 + a_4*c_3*d_6*e_5*f_1 + a_4*c_5*d_1*e_3*f_6 - a_4*c_5*d_1*e_6*f_3 - a_4*c_5*d_3*e_1*f_6 + a_4*c_5*d_3*e_6*f_1 + a_4*c_5*d_6*e_1*f_3 - a_4*c_5*d_6*e_3*f_1 - a_4*c_6*d_1*e_3*f_5 + a_4*c_6*d_1*e_5*f_3 + a_4*c_6*d_3*e_1*f_5 - a_4*c_6*d_3*e_5*f_1 - a_4*c_6*d_5*e_1*f_3 + a_4*c_6*d_5*e_3*f_1 - a_5*c_1*d_3*e_4*f_6 + a_5*c_1*d_3*e_6*f_4 + a_5*c_1*d_4*e_3*f_6 - a_5*c_1*d_4*e_6*f_3 - a_5*c_1*d_6*e_3*f_4 + a_5*c_1*d_6*e_4*f_3 + a_5*c_3*d_1*e_4*f_6 - a_5*c_3*d_1*e_6*f_4 - a_5*c_3*d_4*e_1*f_6 + a_5*c_3*d_4*e_6*f_1 + a_5*c_3*d_6*e_1*f_4 - a_5*c_3*d_6*e_4*f_1 - a_5*c_4*d_1*e_3*f_6 + a_5*c_4*d_1*e_6*f_3 + a_5*c_4*d_3*e_1*f_6 - a_5*c_4*d_3*e_6*f_1 - a_5*c_4*d_6*e_1*f_3 + a_5*c_4*d_6*e_3*f_1 + a_5*c_6*d_1*e_3*f_4 - a_5*c_6*d_1*e_4*f_3 - a_5*c_6*d_3*e_1*f_4 + a_5*c_6*d_3*e_4*f_1 + a_5*c_6*d_4*e_1*f_3 - a_5*c_6*d_4*e_3*f_1 + a_6*c_1*d_3*e_4*f_5 - a_6*c_1*d_3*e_5*f_4 - a_6*c_1*d_4*e_3*f_5 + a_6*c_1*d_4*e_5*f_3 + a_6*c_1*d_5*e_3*f_4 - a_6*c_1*d_5*e_4*f_3 - a_6*c_3*d_1*e_4*f_5 + a_6*c_3*d_1*e_5*f_4 + a_6*c_3*d_4*e_1*f_5 - a_6*c_3*d_4*e_5*f_1 - a_6*c_3*d_5*e_1*f_4 + a_6*c_3*d_5*e_4*f_1 + a_6*c_4*d_1*e_3*f_5 - a_6*c_4*d_1*e_5*f_3 - a_6*c_4*d_3*e_1*f_5 + a_6*c_4*d_3*e_5*f_1 + a_6*c_4*d_5*e_1*f_3 - a_6*c_4*d_5*e_3*f_1 - a_6*c_5*d_1*e_3*f_4 + a_6*c_5*d_1*e_4*f_3 + a_6*c_5*d_3*e_1*f_4 - a_6*c_5*d_3*e_4*f_1 - a_6*c_5*d_4*e_1*f_3 + a_6*c_5*d_4*e_3*f_1)/martix_value
              Inverse_A(3,2)=-(a_1*b_3*d_4*e_5*f_6 - a_1*b_3*d_4*e_6*f_5 - a_1*b_3*d_5*e_4*f_6 + a_1*b_3*d_5*e_6*f_4 + a_1*b_3*d_6*e_4*f_5 - a_1*b_3*d_6*e_5*f_4 - a_1*b_4*d_3*e_5*f_6 + a_1*b_4*d_3*e_6*f_5 + a_1*b_4*d_5*e_3*f_6 - a_1*b_4*d_5*e_6*f_3 - a_1*b_4*d_6*e_3*f_5 + a_1*b_4*d_6*e_5*f_3 + a_1*b_5*d_3*e_4*f_6 - a_1*b_5*d_3*e_6*f_4 - a_1*b_5*d_4*e_3*f_6 + a_1*b_5*d_4*e_6*f_3 + a_1*b_5*d_6*e_3*f_4 - a_1*b_5*d_6*e_4*f_3 - a_1*b_6*d_3*e_4*f_5 + a_1*b_6*d_3*e_5*f_4 + a_1*b_6*d_4*e_3*f_5 - a_1*b_6*d_4*e_5*f_3 - a_1*b_6*d_5*e_3*f_4 + a_1*b_6*d_5*e_4*f_3 - a_3*b_1*d_4*e_5*f_6 + a_3*b_1*d_4*e_6*f_5 + a_3*b_1*d_5*e_4*f_6 - a_3*b_1*d_5*e_6*f_4 - a_3*b_1*d_6*e_4*f_5 + a_3*b_1*d_6*e_5*f_4 + a_3*b_4*d_1*e_5*f_6 - a_3*b_4*d_1*e_6*f_5 - a_3*b_4*d_5*e_1*f_6 + a_3*b_4*d_5*e_6*f_1 + a_3*b_4*d_6*e_1*f_5 - a_3*b_4*d_6*e_5*f_1 - a_3*b_5*d_1*e_4*f_6 + a_3*b_5*d_1*e_6*f_4 + a_3*b_5*d_4*e_1*f_6 - a_3*b_5*d_4*e_6*f_1 - a_3*b_5*d_6*e_1*f_4 + a_3*b_5*d_6*e_4*f_1 + a_3*b_6*d_1*e_4*f_5 - a_3*b_6*d_1*e_5*f_4 - a_3*b_6*d_4*e_1*f_5 + a_3*b_6*d_4*e_5*f_1 + a_3*b_6*d_5*e_1*f_4 - a_3*b_6*d_5*e_4*f_1 + a_4*b_1*d_3*e_5*f_6 - a_4*b_1*d_3*e_6*f_5 - a_4*b_1*d_5*e_3*f_6 + a_4*b_1*d_5*e_6*f_3 + a_4*b_1*d_6*e_3*f_5 - a_4*b_1*d_6*e_5*f_3 - a_4*b_3*d_1*e_5*f_6 + a_4*b_3*d_1*e_6*f_5 + a_4*b_3*d_5*e_1*f_6 - a_4*b_3*d_5*e_6*f_1 - a_4*b_3*d_6*e_1*f_5 + a_4*b_3*d_6*e_5*f_1 + a_4*b_5*d_1*e_3*f_6 - a_4*b_5*d_1*e_6*f_3 - a_4*b_5*d_3*e_1*f_6 + a_4*b_5*d_3*e_6*f_1 + a_4*b_5*d_6*e_1*f_3 - a_4*b_5*d_6*e_3*f_1 - a_4*b_6*d_1*e_3*f_5 + a_4*b_6*d_1*e_5*f_3 + a_4*b_6*d_3*e_1*f_5 - a_4*b_6*d_3*e_5*f_1 - a_4*b_6*d_5*e_1*f_3 + a_4*b_6*d_5*e_3*f_1 - a_5*b_1*d_3*e_4*f_6 + a_5*b_1*d_3*e_6*f_4 + a_5*b_1*d_4*e_3*f_6 - a_5*b_1*d_4*e_6*f_3 - a_5*b_1*d_6*e_3*f_4 + a_5*b_1*d_6*e_4*f_3 + a_5*b_3*d_1*e_4*f_6 - a_5*b_3*d_1*e_6*f_4 - a_5*b_3*d_4*e_1*f_6 + a_5*b_3*d_4*e_6*f_1 + a_5*b_3*d_6*e_1*f_4 - a_5*b_3*d_6*e_4*f_1 - a_5*b_4*d_1*e_3*f_6 + a_5*b_4*d_1*e_6*f_3 + a_5*b_4*d_3*e_1*f_6 - a_5*b_4*d_3*e_6*f_1 - a_5*b_4*d_6*e_1*f_3 + a_5*b_4*d_6*e_3*f_1 + a_5*b_6*d_1*e_3*f_4 - a_5*b_6*d_1*e_4*f_3 - a_5*b_6*d_3*e_1*f_4 + a_5*b_6*d_3*e_4*f_1 + a_5*b_6*d_4*e_1*f_3 - a_5*b_6*d_4*e_3*f_1 + a_6*b_1*d_3*e_4*f_5 - a_6*b_1*d_3*e_5*f_4 - a_6*b_1*d_4*e_3*f_5 + a_6*b_1*d_4*e_5*f_3 + a_6*b_1*d_5*e_3*f_4 - a_6*b_1*d_5*e_4*f_3 - a_6*b_3*d_1*e_4*f_5 + a_6*b_3*d_1*e_5*f_4 + a_6*b_3*d_4*e_1*f_5 - a_6*b_3*d_4*e_5*f_1 - a_6*b_3*d_5*e_1*f_4 + a_6*b_3*d_5*e_4*f_1 + a_6*b_4*d_1*e_3*f_5 - a_6*b_4*d_1*e_5*f_3 - a_6*b_4*d_3*e_1*f_5 + a_6*b_4*d_3*e_5*f_1 + a_6*b_4*d_5*e_1*f_3 - a_6*b_4*d_5*e_3*f_1 - a_6*b_5*d_1*e_3*f_4 + a_6*b_5*d_1*e_4*f_3 + a_6*b_5*d_3*e_1*f_4 - a_6*b_5*d_3*e_4*f_1 - a_6*b_5*d_4*e_1*f_3 + a_6*b_5*d_4*e_3*f_1)/martix_value
              Inverse_A(4,2)=(a_1*b_3*c_4*e_5*f_6 - a_1*b_3*c_4*e_6*f_5 - a_1*b_3*c_5*e_4*f_6 + a_1*b_3*c_5*e_6*f_4 + a_1*b_3*c_6*e_4*f_5 - a_1*b_3*c_6*e_5*f_4 - a_1*b_4*c_3*e_5*f_6 + a_1*b_4*c_3*e_6*f_5 + a_1*b_4*c_5*e_3*f_6 - a_1*b_4*c_5*e_6*f_3 - a_1*b_4*c_6*e_3*f_5 + a_1*b_4*c_6*e_5*f_3 + a_1*b_5*c_3*e_4*f_6 - a_1*b_5*c_3*e_6*f_4 - a_1*b_5*c_4*e_3*f_6 + a_1*b_5*c_4*e_6*f_3 + a_1*b_5*c_6*e_3*f_4 - a_1*b_5*c_6*e_4*f_3 - a_1*b_6*c_3*e_4*f_5 + a_1*b_6*c_3*e_5*f_4 + a_1*b_6*c_4*e_3*f_5 - a_1*b_6*c_4*e_5*f_3 - a_1*b_6*c_5*e_3*f_4 + a_1*b_6*c_5*e_4*f_3 - a_3*b_1*c_4*e_5*f_6 + a_3*b_1*c_4*e_6*f_5 + a_3*b_1*c_5*e_4*f_6 - a_3*b_1*c_5*e_6*f_4 - a_3*b_1*c_6*e_4*f_5 + a_3*b_1*c_6*e_5*f_4 + a_3*b_4*c_1*e_5*f_6 - a_3*b_4*c_1*e_6*f_5 - a_3*b_4*c_5*e_1*f_6 + a_3*b_4*c_5*e_6*f_1 + a_3*b_4*c_6*e_1*f_5 - a_3*b_4*c_6*e_5*f_1 - a_3*b_5*c_1*e_4*f_6 + a_3*b_5*c_1*e_6*f_4 + a_3*b_5*c_4*e_1*f_6 - a_3*b_5*c_4*e_6*f_1 - a_3*b_5*c_6*e_1*f_4 + a_3*b_5*c_6*e_4*f_1 + a_3*b_6*c_1*e_4*f_5 - a_3*b_6*c_1*e_5*f_4 - a_3*b_6*c_4*e_1*f_5 + a_3*b_6*c_4*e_5*f_1 + a_3*b_6*c_5*e_1*f_4 - a_3*b_6*c_5*e_4*f_1 + a_4*b_1*c_3*e_5*f_6 - a_4*b_1*c_3*e_6*f_5 - a_4*b_1*c_5*e_3*f_6 + a_4*b_1*c_5*e_6*f_3 + a_4*b_1*c_6*e_3*f_5 - a_4*b_1*c_6*e_5*f_3 - a_4*b_3*c_1*e_5*f_6 + a_4*b_3*c_1*e_6*f_5 + a_4*b_3*c_5*e_1*f_6 - a_4*b_3*c_5*e_6*f_1 - a_4*b_3*c_6*e_1*f_5 + a_4*b_3*c_6*e_5*f_1 + a_4*b_5*c_1*e_3*f_6 - a_4*b_5*c_1*e_6*f_3 - a_4*b_5*c_3*e_1*f_6 + a_4*b_5*c_3*e_6*f_1 + a_4*b_5*c_6*e_1*f_3 - a_4*b_5*c_6*e_3*f_1 - a_4*b_6*c_1*e_3*f_5 + a_4*b_6*c_1*e_5*f_3 + a_4*b_6*c_3*e_1*f_5 - a_4*b_6*c_3*e_5*f_1 - a_4*b_6*c_5*e_1*f_3 + a_4*b_6*c_5*e_3*f_1 - a_5*b_1*c_3*e_4*f_6 + a_5*b_1*c_3*e_6*f_4 + a_5*b_1*c_4*e_3*f_6 - a_5*b_1*c_4*e_6*f_3 - a_5*b_1*c_6*e_3*f_4 + a_5*b_1*c_6*e_4*f_3 + a_5*b_3*c_1*e_4*f_6 - a_5*b_3*c_1*e_6*f_4 - a_5*b_3*c_4*e_1*f_6 + a_5*b_3*c_4*e_6*f_1 + a_5*b_3*c_6*e_1*f_4 - a_5*b_3*c_6*e_4*f_1 - a_5*b_4*c_1*e_3*f_6 + a_5*b_4*c_1*e_6*f_3 + a_5*b_4*c_3*e_1*f_6 - a_5*b_4*c_3*e_6*f_1 - a_5*b_4*c_6*e_1*f_3 + a_5*b_4*c_6*e_3*f_1 + a_5*b_6*c_1*e_3*f_4 - a_5*b_6*c_1*e_4*f_3 - a_5*b_6*c_3*e_1*f_4 + a_5*b_6*c_3*e_4*f_1 + a_5*b_6*c_4*e_1*f_3 - a_5*b_6*c_4*e_3*f_1 + a_6*b_1*c_3*e_4*f_5 - a_6*b_1*c_3*e_5*f_4 - a_6*b_1*c_4*e_3*f_5 + a_6*b_1*c_4*e_5*f_3 + a_6*b_1*c_5*e_3*f_4 - a_6*b_1*c_5*e_4*f_3 - a_6*b_3*c_1*e_4*f_5 + a_6*b_3*c_1*e_5*f_4 + a_6*b_3*c_4*e_1*f_5 - a_6*b_3*c_4*e_5*f_1 - a_6*b_3*c_5*e_1*f_4 + a_6*b_3*c_5*e_4*f_1 + a_6*b_4*c_1*e_3*f_5 - a_6*b_4*c_1*e_5*f_3 - a_6*b_4*c_3*e_1*f_5 + a_6*b_4*c_3*e_5*f_1 + a_6*b_4*c_5*e_1*f_3 - a_6*b_4*c_5*e_3*f_1 - a_6*b_5*c_1*e_3*f_4 + a_6*b_5*c_1*e_4*f_3 + a_6*b_5*c_3*e_1*f_4 - a_6*b_5*c_3*e_4*f_1 - a_6*b_5*c_4*e_1*f_3 + a_6*b_5*c_4*e_3*f_1)/martix_value
              Inverse_A(5,2)=-(a_1*b_3*c_4*d_5*f_6 - a_1*b_3*c_4*d_6*f_5 - a_1*b_3*c_5*d_4*f_6 + a_1*b_3*c_5*d_6*f_4 + a_1*b_3*c_6*d_4*f_5 - a_1*b_3*c_6*d_5*f_4 - a_1*b_4*c_3*d_5*f_6 + a_1*b_4*c_3*d_6*f_5 + a_1*b_4*c_5*d_3*f_6 - a_1*b_4*c_5*d_6*f_3 - a_1*b_4*c_6*d_3*f_5 + a_1*b_4*c_6*d_5*f_3 + a_1*b_5*c_3*d_4*f_6 - a_1*b_5*c_3*d_6*f_4 - a_1*b_5*c_4*d_3*f_6 + a_1*b_5*c_4*d_6*f_3 + a_1*b_5*c_6*d_3*f_4 - a_1*b_5*c_6*d_4*f_3 - a_1*b_6*c_3*d_4*f_5 + a_1*b_6*c_3*d_5*f_4 + a_1*b_6*c_4*d_3*f_5 - a_1*b_6*c_4*d_5*f_3 - a_1*b_6*c_5*d_3*f_4 + a_1*b_6*c_5*d_4*f_3 - a_3*b_1*c_4*d_5*f_6 + a_3*b_1*c_4*d_6*f_5 + a_3*b_1*c_5*d_4*f_6 - a_3*b_1*c_5*d_6*f_4 - a_3*b_1*c_6*d_4*f_5 + a_3*b_1*c_6*d_5*f_4 + a_3*b_4*c_1*d_5*f_6 - a_3*b_4*c_1*d_6*f_5 - a_3*b_4*c_5*d_1*f_6 + a_3*b_4*c_5*d_6*f_1 + a_3*b_4*c_6*d_1*f_5 - a_3*b_4*c_6*d_5*f_1 - a_3*b_5*c_1*d_4*f_6 + a_3*b_5*c_1*d_6*f_4 + a_3*b_5*c_4*d_1*f_6 - a_3*b_5*c_4*d_6*f_1 - a_3*b_5*c_6*d_1*f_4 + a_3*b_5*c_6*d_4*f_1 + a_3*b_6*c_1*d_4*f_5 - a_3*b_6*c_1*d_5*f_4 - a_3*b_6*c_4*d_1*f_5 + a_3*b_6*c_4*d_5*f_1 + a_3*b_6*c_5*d_1*f_4 - a_3*b_6*c_5*d_4*f_1 + a_4*b_1*c_3*d_5*f_6 - a_4*b_1*c_3*d_6*f_5 - a_4*b_1*c_5*d_3*f_6 + a_4*b_1*c_5*d_6*f_3 + a_4*b_1*c_6*d_3*f_5 - a_4*b_1*c_6*d_5*f_3 - a_4*b_3*c_1*d_5*f_6 + a_4*b_3*c_1*d_6*f_5 + a_4*b_3*c_5*d_1*f_6 - a_4*b_3*c_5*d_6*f_1 - a_4*b_3*c_6*d_1*f_5 + a_4*b_3*c_6*d_5*f_1 + a_4*b_5*c_1*d_3*f_6 - a_4*b_5*c_1*d_6*f_3 - a_4*b_5*c_3*d_1*f_6 + a_4*b_5*c_3*d_6*f_1 + a_4*b_5*c_6*d_1*f_3 - a_4*b_5*c_6*d_3*f_1 - a_4*b_6*c_1*d_3*f_5 + a_4*b_6*c_1*d_5*f_3 + a_4*b_6*c_3*d_1*f_5 - a_4*b_6*c_3*d_5*f_1 - a_4*b_6*c_5*d_1*f_3 + a_4*b_6*c_5*d_3*f_1 - a_5*b_1*c_3*d_4*f_6 + a_5*b_1*c_3*d_6*f_4 + a_5*b_1*c_4*d_3*f_6 - a_5*b_1*c_4*d_6*f_3 - a_5*b_1*c_6*d_3*f_4 + a_5*b_1*c_6*d_4*f_3 + a_5*b_3*c_1*d_4*f_6 - a_5*b_3*c_1*d_6*f_4 - a_5*b_3*c_4*d_1*f_6 + a_5*b_3*c_4*d_6*f_1 + a_5*b_3*c_6*d_1*f_4 - a_5*b_3*c_6*d_4*f_1 - a_5*b_4*c_1*d_3*f_6 + a_5*b_4*c_1*d_6*f_3 + a_5*b_4*c_3*d_1*f_6 - a_5*b_4*c_3*d_6*f_1 - a_5*b_4*c_6*d_1*f_3 + a_5*b_4*c_6*d_3*f_1 + a_5*b_6*c_1*d_3*f_4 - a_5*b_6*c_1*d_4*f_3 - a_5*b_6*c_3*d_1*f_4 + a_5*b_6*c_3*d_4*f_1 + a_5*b_6*c_4*d_1*f_3 - a_5*b_6*c_4*d_3*f_1 + a_6*b_1*c_3*d_4*f_5 - a_6*b_1*c_3*d_5*f_4 - a_6*b_1*c_4*d_3*f_5 + a_6*b_1*c_4*d_5*f_3 + a_6*b_1*c_5*d_3*f_4 - a_6*b_1*c_5*d_4*f_3 - a_6*b_3*c_1*d_4*f_5 + a_6*b_3*c_1*d_5*f_4 + a_6*b_3*c_4*d_1*f_5 - a_6*b_3*c_4*d_5*f_1 - a_6*b_3*c_5*d_1*f_4 + a_6*b_3*c_5*d_4*f_1 + a_6*b_4*c_1*d_3*f_5 - a_6*b_4*c_1*d_5*f_3 - a_6*b_4*c_3*d_1*f_5 + a_6*b_4*c_3*d_5*f_1 + a_6*b_4*c_5*d_1*f_3 - a_6*b_4*c_5*d_3*f_1 - a_6*b_5*c_1*d_3*f_4 + a_6*b_5*c_1*d_4*f_3 + a_6*b_5*c_3*d_1*f_4 - a_6*b_5*c_3*d_4*f_1 - a_6*b_5*c_4*d_1*f_3 + a_6*b_5*c_4*d_3*f_1)/martix_value
              Inverse_A(6,2)=(a_1*b_3*c_4*d_5*e_6 - a_1*b_3*c_4*d_6*e_5 - a_1*b_3*c_5*d_4*e_6 + a_1*b_3*c_5*d_6*e_4 + a_1*b_3*c_6*d_4*e_5 - a_1*b_3*c_6*d_5*e_4 - a_1*b_4*c_3*d_5*e_6 + a_1*b_4*c_3*d_6*e_5 + a_1*b_4*c_5*d_3*e_6 - a_1*b_4*c_5*d_6*e_3 - a_1*b_4*c_6*d_3*e_5 + a_1*b_4*c_6*d_5*e_3 + a_1*b_5*c_3*d_4*e_6 - a_1*b_5*c_3*d_6*e_4 - a_1*b_5*c_4*d_3*e_6 + a_1*b_5*c_4*d_6*e_3 + a_1*b_5*c_6*d_3*e_4 - a_1*b_5*c_6*d_4*e_3 - a_1*b_6*c_3*d_4*e_5 + a_1*b_6*c_3*d_5*e_4 + a_1*b_6*c_4*d_3*e_5 - a_1*b_6*c_4*d_5*e_3 - a_1*b_6*c_5*d_3*e_4 + a_1*b_6*c_5*d_4*e_3 - a_3*b_1*c_4*d_5*e_6 + a_3*b_1*c_4*d_6*e_5 + a_3*b_1*c_5*d_4*e_6 - a_3*b_1*c_5*d_6*e_4 - a_3*b_1*c_6*d_4*e_5 + a_3*b_1*c_6*d_5*e_4 + a_3*b_4*c_1*d_5*e_6 - a_3*b_4*c_1*d_6*e_5 - a_3*b_4*c_5*d_1*e_6 + a_3*b_4*c_5*d_6*e_1 + a_3*b_4*c_6*d_1*e_5 - a_3*b_4*c_6*d_5*e_1 - a_3*b_5*c_1*d_4*e_6 + a_3*b_5*c_1*d_6*e_4 + a_3*b_5*c_4*d_1*e_6 - a_3*b_5*c_4*d_6*e_1 - a_3*b_5*c_6*d_1*e_4 + a_3*b_5*c_6*d_4*e_1 + a_3*b_6*c_1*d_4*e_5 - a_3*b_6*c_1*d_5*e_4 - a_3*b_6*c_4*d_1*e_5 + a_3*b_6*c_4*d_5*e_1 + a_3*b_6*c_5*d_1*e_4 - a_3*b_6*c_5*d_4*e_1 + a_4*b_1*c_3*d_5*e_6 - a_4*b_1*c_3*d_6*e_5 - a_4*b_1*c_5*d_3*e_6 + a_4*b_1*c_5*d_6*e_3 + a_4*b_1*c_6*d_3*e_5 - a_4*b_1*c_6*d_5*e_3 - a_4*b_3*c_1*d_5*e_6 + a_4*b_3*c_1*d_6*e_5 + a_4*b_3*c_5*d_1*e_6 - a_4*b_3*c_5*d_6*e_1 - a_4*b_3*c_6*d_1*e_5 + a_4*b_3*c_6*d_5*e_1 + a_4*b_5*c_1*d_3*e_6 - a_4*b_5*c_1*d_6*e_3 - a_4*b_5*c_3*d_1*e_6 + a_4*b_5*c_3*d_6*e_1 + a_4*b_5*c_6*d_1*e_3 - a_4*b_5*c_6*d_3*e_1 - a_4*b_6*c_1*d_3*e_5 + a_4*b_6*c_1*d_5*e_3 + a_4*b_6*c_3*d_1*e_5 - a_4*b_6*c_3*d_5*e_1 - a_4*b_6*c_5*d_1*e_3 + a_4*b_6*c_5*d_3*e_1 - a_5*b_1*c_3*d_4*e_6 + a_5*b_1*c_3*d_6*e_4 + a_5*b_1*c_4*d_3*e_6 - a_5*b_1*c_4*d_6*e_3 - a_5*b_1*c_6*d_3*e_4 + a_5*b_1*c_6*d_4*e_3 + a_5*b_3*c_1*d_4*e_6 - a_5*b_3*c_1*d_6*e_4 - a_5*b_3*c_4*d_1*e_6 + a_5*b_3*c_4*d_6*e_1 + a_5*b_3*c_6*d_1*e_4 - a_5*b_3*c_6*d_4*e_1 - a_5*b_4*c_1*d_3*e_6 + a_5*b_4*c_1*d_6*e_3 + a_5*b_4*c_3*d_1*e_6 - a_5*b_4*c_3*d_6*e_1 - a_5*b_4*c_6*d_1*e_3 + a_5*b_4*c_6*d_3*e_1 + a_5*b_6*c_1*d_3*e_4 - a_5*b_6*c_1*d_4*e_3 - a_5*b_6*c_3*d_1*e_4 + a_5*b_6*c_3*d_4*e_1 + a_5*b_6*c_4*d_1*e_3 - a_5*b_6*c_4*d_3*e_1 + a_6*b_1*c_3*d_4*e_5 - a_6*b_1*c_3*d_5*e_4 - a_6*b_1*c_4*d_3*e_5 + a_6*b_1*c_4*d_5*e_3 + a_6*b_1*c_5*d_3*e_4 - a_6*b_1*c_5*d_4*e_3 - a_6*b_3*c_1*d_4*e_5 + a_6*b_3*c_1*d_5*e_4 + a_6*b_3*c_4*d_1*e_5 - a_6*b_3*c_4*d_5*e_1 - a_6*b_3*c_5*d_1*e_4 + a_6*b_3*c_5*d_4*e_1 + a_6*b_4*c_1*d_3*e_5 - a_6*b_4*c_1*d_5*e_3 - a_6*b_4*c_3*d_1*e_5 + a_6*b_4*c_3*d_5*e_1 + a_6*b_4*c_5*d_1*e_3 - a_6*b_4*c_5*d_3*e_1 - a_6*b_5*c_1*d_3*e_4 + a_6*b_5*c_1*d_4*e_3 + a_6*b_5*c_3*d_1*e_4 - a_6*b_5*c_3*d_4*e_1 - a_6*b_5*c_4*d_1*e_3 + a_6*b_5*c_4*d_3*e_1)/martix_value
      
              Inverse_A(1,3)=(b_1*c_2*d_4*e_5*f_6 - b_1*c_2*d_4*e_6*f_5 - b_1*c_2*d_5*e_4*f_6 + b_1*c_2*d_5*e_6*f_4 + b_1*c_2*d_6*e_4*f_5 - b_1*c_2*d_6*e_5*f_4 - b_1*c_4*d_2*e_5*f_6 + b_1*c_4*d_2*e_6*f_5 + b_1*c_4*d_5*e_2*f_6 - b_1*c_4*d_5*e_6*f_2 - b_1*c_4*d_6*e_2*f_5 + b_1*c_4*d_6*e_5*f_2 + b_1*c_5*d_2*e_4*f_6 - b_1*c_5*d_2*e_6*f_4 - b_1*c_5*d_4*e_2*f_6 + b_1*c_5*d_4*e_6*f_2 + b_1*c_5*d_6*e_2*f_4 - b_1*c_5*d_6*e_4*f_2 - b_1*c_6*d_2*e_4*f_5 + b_1*c_6*d_2*e_5*f_4 + b_1*c_6*d_4*e_2*f_5 - b_1*c_6*d_4*e_5*f_2 - b_1*c_6*d_5*e_2*f_4 + b_1*c_6*d_5*e_4*f_2 - b_2*c_1*d_4*e_5*f_6 + b_2*c_1*d_4*e_6*f_5 + b_2*c_1*d_5*e_4*f_6 - b_2*c_1*d_5*e_6*f_4 - b_2*c_1*d_6*e_4*f_5 + b_2*c_1*d_6*e_5*f_4 + b_2*c_4*d_1*e_5*f_6 - b_2*c_4*d_1*e_6*f_5 - b_2*c_4*d_5*e_1*f_6 + b_2*c_4*d_5*e_6*f_1 + b_2*c_4*d_6*e_1*f_5 - b_2*c_4*d_6*e_5*f_1 - b_2*c_5*d_1*e_4*f_6 + b_2*c_5*d_1*e_6*f_4 + b_2*c_5*d_4*e_1*f_6 - b_2*c_5*d_4*e_6*f_1 - b_2*c_5*d_6*e_1*f_4 + b_2*c_5*d_6*e_4*f_1 + b_2*c_6*d_1*e_4*f_5 - b_2*c_6*d_1*e_5*f_4 - b_2*c_6*d_4*e_1*f_5 + b_2*c_6*d_4*e_5*f_1 + b_2*c_6*d_5*e_1*f_4 - b_2*c_6*d_5*e_4*f_1 + b_4*c_1*d_2*e_5*f_6 - b_4*c_1*d_2*e_6*f_5 - b_4*c_1*d_5*e_2*f_6 + b_4*c_1*d_5*e_6*f_2 + b_4*c_1*d_6*e_2*f_5 - b_4*c_1*d_6*e_5*f_2 - b_4*c_2*d_1*e_5*f_6 + b_4*c_2*d_1*e_6*f_5 + b_4*c_2*d_5*e_1*f_6 - b_4*c_2*d_5*e_6*f_1 - b_4*c_2*d_6*e_1*f_5 + b_4*c_2*d_6*e_5*f_1 + b_4*c_5*d_1*e_2*f_6 - b_4*c_5*d_1*e_6*f_2 - b_4*c_5*d_2*e_1*f_6 + b_4*c_5*d_2*e_6*f_1 + b_4*c_5*d_6*e_1*f_2 - b_4*c_5*d_6*e_2*f_1 - b_4*c_6*d_1*e_2*f_5 + b_4*c_6*d_1*e_5*f_2 + b_4*c_6*d_2*e_1*f_5 - b_4*c_6*d_2*e_5*f_1 - b_4*c_6*d_5*e_1*f_2 + b_4*c_6*d_5*e_2*f_1 - b_5*c_1*d_2*e_4*f_6 + b_5*c_1*d_2*e_6*f_4 + b_5*c_1*d_4*e_2*f_6 - b_5*c_1*d_4*e_6*f_2 - b_5*c_1*d_6*e_2*f_4 + b_5*c_1*d_6*e_4*f_2 + b_5*c_2*d_1*e_4*f_6 - b_5*c_2*d_1*e_6*f_4 - b_5*c_2*d_4*e_1*f_6 + b_5*c_2*d_4*e_6*f_1 + b_5*c_2*d_6*e_1*f_4 - b_5*c_2*d_6*e_4*f_1 - b_5*c_4*d_1*e_2*f_6 + b_5*c_4*d_1*e_6*f_2 + b_5*c_4*d_2*e_1*f_6 - b_5*c_4*d_2*e_6*f_1 - b_5*c_4*d_6*e_1*f_2 + b_5*c_4*d_6*e_2*f_1 + b_5*c_6*d_1*e_2*f_4 - b_5*c_6*d_1*e_4*f_2 - b_5*c_6*d_2*e_1*f_4 + b_5*c_6*d_2*e_4*f_1 + b_5*c_6*d_4*e_1*f_2 - b_5*c_6*d_4*e_2*f_1 + b_6*c_1*d_2*e_4*f_5 - b_6*c_1*d_2*e_5*f_4 - b_6*c_1*d_4*e_2*f_5 + b_6*c_1*d_4*e_5*f_2 + b_6*c_1*d_5*e_2*f_4 - b_6*c_1*d_5*e_4*f_2 - b_6*c_2*d_1*e_4*f_5 + b_6*c_2*d_1*e_5*f_4 + b_6*c_2*d_4*e_1*f_5 - b_6*c_2*d_4*e_5*f_1 - b_6*c_2*d_5*e_1*f_4 + b_6*c_2*d_5*e_4*f_1 + b_6*c_4*d_1*e_2*f_5 - b_6*c_4*d_1*e_5*f_2 - b_6*c_4*d_2*e_1*f_5 + b_6*c_4*d_2*e_5*f_1 + b_6*c_4*d_5*e_1*f_2 - b_6*c_4*d_5*e_2*f_1 - b_6*c_5*d_1*e_2*f_4 + b_6*c_5*d_1*e_4*f_2 + b_6*c_5*d_2*e_1*f_4 - b_6*c_5*d_2*e_4*f_1 - b_6*c_5*d_4*e_1*f_2 + b_6*c_5*d_4*e_2*f_1)/martix_value
              Inverse_A(2,3)=-(a_1*c_2*d_4*e_5*f_6 - a_1*c_2*d_4*e_6*f_5 - a_1*c_2*d_5*e_4*f_6 + a_1*c_2*d_5*e_6*f_4 + a_1*c_2*d_6*e_4*f_5 - a_1*c_2*d_6*e_5*f_4 - a_1*c_4*d_2*e_5*f_6 + a_1*c_4*d_2*e_6*f_5 + a_1*c_4*d_5*e_2*f_6 - a_1*c_4*d_5*e_6*f_2 - a_1*c_4*d_6*e_2*f_5 + a_1*c_4*d_6*e_5*f_2 + a_1*c_5*d_2*e_4*f_6 - a_1*c_5*d_2*e_6*f_4 - a_1*c_5*d_4*e_2*f_6 + a_1*c_5*d_4*e_6*f_2 + a_1*c_5*d_6*e_2*f_4 - a_1*c_5*d_6*e_4*f_2 - a_1*c_6*d_2*e_4*f_5 + a_1*c_6*d_2*e_5*f_4 + a_1*c_6*d_4*e_2*f_5 - a_1*c_6*d_4*e_5*f_2 - a_1*c_6*d_5*e_2*f_4 + a_1*c_6*d_5*e_4*f_2 - a_2*c_1*d_4*e_5*f_6 + a_2*c_1*d_4*e_6*f_5 + a_2*c_1*d_5*e_4*f_6 - a_2*c_1*d_5*e_6*f_4 - a_2*c_1*d_6*e_4*f_5 + a_2*c_1*d_6*e_5*f_4 + a_2*c_4*d_1*e_5*f_6 - a_2*c_4*d_1*e_6*f_5 - a_2*c_4*d_5*e_1*f_6 + a_2*c_4*d_5*e_6*f_1 + a_2*c_4*d_6*e_1*f_5 - a_2*c_4*d_6*e_5*f_1 - a_2*c_5*d_1*e_4*f_6 + a_2*c_5*d_1*e_6*f_4 + a_2*c_5*d_4*e_1*f_6 - a_2*c_5*d_4*e_6*f_1 - a_2*c_5*d_6*e_1*f_4 + a_2*c_5*d_6*e_4*f_1 + a_2*c_6*d_1*e_4*f_5 - a_2*c_6*d_1*e_5*f_4 - a_2*c_6*d_4*e_1*f_5 + a_2*c_6*d_4*e_5*f_1 + a_2*c_6*d_5*e_1*f_4 - a_2*c_6*d_5*e_4*f_1 + a_4*c_1*d_2*e_5*f_6 - a_4*c_1*d_2*e_6*f_5 - a_4*c_1*d_5*e_2*f_6 + a_4*c_1*d_5*e_6*f_2 + a_4*c_1*d_6*e_2*f_5 - a_4*c_1*d_6*e_5*f_2 - a_4*c_2*d_1*e_5*f_6 + a_4*c_2*d_1*e_6*f_5 + a_4*c_2*d_5*e_1*f_6 - a_4*c_2*d_5*e_6*f_1 - a_4*c_2*d_6*e_1*f_5 + a_4*c_2*d_6*e_5*f_1 + a_4*c_5*d_1*e_2*f_6 - a_4*c_5*d_1*e_6*f_2 - a_4*c_5*d_2*e_1*f_6 + a_4*c_5*d_2*e_6*f_1 + a_4*c_5*d_6*e_1*f_2 - a_4*c_5*d_6*e_2*f_1 - a_4*c_6*d_1*e_2*f_5 + a_4*c_6*d_1*e_5*f_2 + a_4*c_6*d_2*e_1*f_5 - a_4*c_6*d_2*e_5*f_1 - a_4*c_6*d_5*e_1*f_2 + a_4*c_6*d_5*e_2*f_1 - a_5*c_1*d_2*e_4*f_6 + a_5*c_1*d_2*e_6*f_4 + a_5*c_1*d_4*e_2*f_6 - a_5*c_1*d_4*e_6*f_2 - a_5*c_1*d_6*e_2*f_4 + a_5*c_1*d_6*e_4*f_2 + a_5*c_2*d_1*e_4*f_6 - a_5*c_2*d_1*e_6*f_4 - a_5*c_2*d_4*e_1*f_6 + a_5*c_2*d_4*e_6*f_1 + a_5*c_2*d_6*e_1*f_4 - a_5*c_2*d_6*e_4*f_1 - a_5*c_4*d_1*e_2*f_6 + a_5*c_4*d_1*e_6*f_2 + a_5*c_4*d_2*e_1*f_6 - a_5*c_4*d_2*e_6*f_1 - a_5*c_4*d_6*e_1*f_2 + a_5*c_4*d_6*e_2*f_1 + a_5*c_6*d_1*e_2*f_4 - a_5*c_6*d_1*e_4*f_2 - a_5*c_6*d_2*e_1*f_4 + a_5*c_6*d_2*e_4*f_1 + a_5*c_6*d_4*e_1*f_2 - a_5*c_6*d_4*e_2*f_1 + a_6*c_1*d_2*e_4*f_5 - a_6*c_1*d_2*e_5*f_4 - a_6*c_1*d_4*e_2*f_5 + a_6*c_1*d_4*e_5*f_2 + a_6*c_1*d_5*e_2*f_4 - a_6*c_1*d_5*e_4*f_2 - a_6*c_2*d_1*e_4*f_5 + a_6*c_2*d_1*e_5*f_4 + a_6*c_2*d_4*e_1*f_5 - a_6*c_2*d_4*e_5*f_1 - a_6*c_2*d_5*e_1*f_4 + a_6*c_2*d_5*e_4*f_1 + a_6*c_4*d_1*e_2*f_5 - a_6*c_4*d_1*e_5*f_2 - a_6*c_4*d_2*e_1*f_5 + a_6*c_4*d_2*e_5*f_1 + a_6*c_4*d_5*e_1*f_2 - a_6*c_4*d_5*e_2*f_1 - a_6*c_5*d_1*e_2*f_4 + a_6*c_5*d_1*e_4*f_2 + a_6*c_5*d_2*e_1*f_4 - a_6*c_5*d_2*e_4*f_1 - a_6*c_5*d_4*e_1*f_2 + a_6*c_5*d_4*e_2*f_1)/martix_value
              Inverse_A(3,3)=(a_1*b_2*d_4*e_5*f_6 - a_1*b_2*d_4*e_6*f_5 - a_1*b_2*d_5*e_4*f_6 + a_1*b_2*d_5*e_6*f_4 + a_1*b_2*d_6*e_4*f_5 - a_1*b_2*d_6*e_5*f_4 - a_1*b_4*d_2*e_5*f_6 + a_1*b_4*d_2*e_6*f_5 + a_1*b_4*d_5*e_2*f_6 - a_1*b_4*d_5*e_6*f_2 - a_1*b_4*d_6*e_2*f_5 + a_1*b_4*d_6*e_5*f_2 + a_1*b_5*d_2*e_4*f_6 - a_1*b_5*d_2*e_6*f_4 - a_1*b_5*d_4*e_2*f_6 + a_1*b_5*d_4*e_6*f_2 + a_1*b_5*d_6*e_2*f_4 - a_1*b_5*d_6*e_4*f_2 - a_1*b_6*d_2*e_4*f_5 + a_1*b_6*d_2*e_5*f_4 + a_1*b_6*d_4*e_2*f_5 - a_1*b_6*d_4*e_5*f_2 - a_1*b_6*d_5*e_2*f_4 + a_1*b_6*d_5*e_4*f_2 - a_2*b_1*d_4*e_5*f_6 + a_2*b_1*d_4*e_6*f_5 + a_2*b_1*d_5*e_4*f_6 - a_2*b_1*d_5*e_6*f_4 - a_2*b_1*d_6*e_4*f_5 + a_2*b_1*d_6*e_5*f_4 + a_2*b_4*d_1*e_5*f_6 - a_2*b_4*d_1*e_6*f_5 - a_2*b_4*d_5*e_1*f_6 + a_2*b_4*d_5*e_6*f_1 + a_2*b_4*d_6*e_1*f_5 - a_2*b_4*d_6*e_5*f_1 - a_2*b_5*d_1*e_4*f_6 + a_2*b_5*d_1*e_6*f_4 + a_2*b_5*d_4*e_1*f_6 - a_2*b_5*d_4*e_6*f_1 - a_2*b_5*d_6*e_1*f_4 + a_2*b_5*d_6*e_4*f_1 + a_2*b_6*d_1*e_4*f_5 - a_2*b_6*d_1*e_5*f_4 - a_2*b_6*d_4*e_1*f_5 + a_2*b_6*d_4*e_5*f_1 + a_2*b_6*d_5*e_1*f_4 - a_2*b_6*d_5*e_4*f_1 + a_4*b_1*d_2*e_5*f_6 - a_4*b_1*d_2*e_6*f_5 - a_4*b_1*d_5*e_2*f_6 + a_4*b_1*d_5*e_6*f_2 + a_4*b_1*d_6*e_2*f_5 - a_4*b_1*d_6*e_5*f_2 - a_4*b_2*d_1*e_5*f_6 + a_4*b_2*d_1*e_6*f_5 + a_4*b_2*d_5*e_1*f_6 - a_4*b_2*d_5*e_6*f_1 - a_4*b_2*d_6*e_1*f_5 + a_4*b_2*d_6*e_5*f_1 + a_4*b_5*d_1*e_2*f_6 - a_4*b_5*d_1*e_6*f_2 - a_4*b_5*d_2*e_1*f_6 + a_4*b_5*d_2*e_6*f_1 + a_4*b_5*d_6*e_1*f_2 - a_4*b_5*d_6*e_2*f_1 - a_4*b_6*d_1*e_2*f_5 + a_4*b_6*d_1*e_5*f_2 + a_4*b_6*d_2*e_1*f_5 - a_4*b_6*d_2*e_5*f_1 - a_4*b_6*d_5*e_1*f_2 + a_4*b_6*d_5*e_2*f_1 - a_5*b_1*d_2*e_4*f_6 + a_5*b_1*d_2*e_6*f_4 + a_5*b_1*d_4*e_2*f_6 - a_5*b_1*d_4*e_6*f_2 - a_5*b_1*d_6*e_2*f_4 + a_5*b_1*d_6*e_4*f_2 + a_5*b_2*d_1*e_4*f_6 - a_5*b_2*d_1*e_6*f_4 - a_5*b_2*d_4*e_1*f_6 + a_5*b_2*d_4*e_6*f_1 + a_5*b_2*d_6*e_1*f_4 - a_5*b_2*d_6*e_4*f_1 - a_5*b_4*d_1*e_2*f_6 + a_5*b_4*d_1*e_6*f_2 + a_5*b_4*d_2*e_1*f_6 - a_5*b_4*d_2*e_6*f_1 - a_5*b_4*d_6*e_1*f_2 + a_5*b_4*d_6*e_2*f_1 + a_5*b_6*d_1*e_2*f_4 - a_5*b_6*d_1*e_4*f_2 - a_5*b_6*d_2*e_1*f_4 + a_5*b_6*d_2*e_4*f_1 + a_5*b_6*d_4*e_1*f_2 - a_5*b_6*d_4*e_2*f_1 + a_6*b_1*d_2*e_4*f_5 - a_6*b_1*d_2*e_5*f_4 - a_6*b_1*d_4*e_2*f_5 + a_6*b_1*d_4*e_5*f_2 + a_6*b_1*d_5*e_2*f_4 - a_6*b_1*d_5*e_4*f_2 - a_6*b_2*d_1*e_4*f_5 + a_6*b_2*d_1*e_5*f_4 + a_6*b_2*d_4*e_1*f_5 - a_6*b_2*d_4*e_5*f_1 - a_6*b_2*d_5*e_1*f_4 + a_6*b_2*d_5*e_4*f_1 + a_6*b_4*d_1*e_2*f_5 - a_6*b_4*d_1*e_5*f_2 - a_6*b_4*d_2*e_1*f_5 + a_6*b_4*d_2*e_5*f_1 + a_6*b_4*d_5*e_1*f_2 - a_6*b_4*d_5*e_2*f_1 - a_6*b_5*d_1*e_2*f_4 + a_6*b_5*d_1*e_4*f_2 + a_6*b_5*d_2*e_1*f_4 - a_6*b_5*d_2*e_4*f_1 - a_6*b_5*d_4*e_1*f_2 + a_6*b_5*d_4*e_2*f_1)/martix_value
              Inverse_A(4,3)=-(a_1*b_2*c_4*e_5*f_6 - a_1*b_2*c_4*e_6*f_5 - a_1*b_2*c_5*e_4*f_6 + a_1*b_2*c_5*e_6*f_4 + a_1*b_2*c_6*e_4*f_5 - a_1*b_2*c_6*e_5*f_4 - a_1*b_4*c_2*e_5*f_6 + a_1*b_4*c_2*e_6*f_5 + a_1*b_4*c_5*e_2*f_6 - a_1*b_4*c_5*e_6*f_2 - a_1*b_4*c_6*e_2*f_5 + a_1*b_4*c_6*e_5*f_2 + a_1*b_5*c_2*e_4*f_6 - a_1*b_5*c_2*e_6*f_4 - a_1*b_5*c_4*e_2*f_6 + a_1*b_5*c_4*e_6*f_2 + a_1*b_5*c_6*e_2*f_4 - a_1*b_5*c_6*e_4*f_2 - a_1*b_6*c_2*e_4*f_5 + a_1*b_6*c_2*e_5*f_4 + a_1*b_6*c_4*e_2*f_5 - a_1*b_6*c_4*e_5*f_2 - a_1*b_6*c_5*e_2*f_4 + a_1*b_6*c_5*e_4*f_2 - a_2*b_1*c_4*e_5*f_6 + a_2*b_1*c_4*e_6*f_5 + a_2*b_1*c_5*e_4*f_6 - a_2*b_1*c_5*e_6*f_4 - a_2*b_1*c_6*e_4*f_5 + a_2*b_1*c_6*e_5*f_4 + a_2*b_4*c_1*e_5*f_6 - a_2*b_4*c_1*e_6*f_5 - a_2*b_4*c_5*e_1*f_6 + a_2*b_4*c_5*e_6*f_1 + a_2*b_4*c_6*e_1*f_5 - a_2*b_4*c_6*e_5*f_1 - a_2*b_5*c_1*e_4*f_6 + a_2*b_5*c_1*e_6*f_4 + a_2*b_5*c_4*e_1*f_6 - a_2*b_5*c_4*e_6*f_1 - a_2*b_5*c_6*e_1*f_4 + a_2*b_5*c_6*e_4*f_1 + a_2*b_6*c_1*e_4*f_5 - a_2*b_6*c_1*e_5*f_4 - a_2*b_6*c_4*e_1*f_5 + a_2*b_6*c_4*e_5*f_1 + a_2*b_6*c_5*e_1*f_4 - a_2*b_6*c_5*e_4*f_1 + a_4*b_1*c_2*e_5*f_6 - a_4*b_1*c_2*e_6*f_5 - a_4*b_1*c_5*e_2*f_6 + a_4*b_1*c_5*e_6*f_2 + a_4*b_1*c_6*e_2*f_5 - a_4*b_1*c_6*e_5*f_2 - a_4*b_2*c_1*e_5*f_6 + a_4*b_2*c_1*e_6*f_5 + a_4*b_2*c_5*e_1*f_6 - a_4*b_2*c_5*e_6*f_1 - a_4*b_2*c_6*e_1*f_5 + a_4*b_2*c_6*e_5*f_1 + a_4*b_5*c_1*e_2*f_6 - a_4*b_5*c_1*e_6*f_2 - a_4*b_5*c_2*e_1*f_6 + a_4*b_5*c_2*e_6*f_1 + a_4*b_5*c_6*e_1*f_2 - a_4*b_5*c_6*e_2*f_1 - a_4*b_6*c_1*e_2*f_5 + a_4*b_6*c_1*e_5*f_2 + a_4*b_6*c_2*e_1*f_5 - a_4*b_6*c_2*e_5*f_1 - a_4*b_6*c_5*e_1*f_2 + a_4*b_6*c_5*e_2*f_1 - a_5*b_1*c_2*e_4*f_6 + a_5*b_1*c_2*e_6*f_4 + a_5*b_1*c_4*e_2*f_6 - a_5*b_1*c_4*e_6*f_2 - a_5*b_1*c_6*e_2*f_4 + a_5*b_1*c_6*e_4*f_2 + a_5*b_2*c_1*e_4*f_6 - a_5*b_2*c_1*e_6*f_4 - a_5*b_2*c_4*e_1*f_6 + a_5*b_2*c_4*e_6*f_1 + a_5*b_2*c_6*e_1*f_4 - a_5*b_2*c_6*e_4*f_1 - a_5*b_4*c_1*e_2*f_6 + a_5*b_4*c_1*e_6*f_2 + a_5*b_4*c_2*e_1*f_6 - a_5*b_4*c_2*e_6*f_1 - a_5*b_4*c_6*e_1*f_2 + a_5*b_4*c_6*e_2*f_1 + a_5*b_6*c_1*e_2*f_4 - a_5*b_6*c_1*e_4*f_2 - a_5*b_6*c_2*e_1*f_4 + a_5*b_6*c_2*e_4*f_1 + a_5*b_6*c_4*e_1*f_2 - a_5*b_6*c_4*e_2*f_1 + a_6*b_1*c_2*e_4*f_5 - a_6*b_1*c_2*e_5*f_4 - a_6*b_1*c_4*e_2*f_5 + a_6*b_1*c_4*e_5*f_2 + a_6*b_1*c_5*e_2*f_4 - a_6*b_1*c_5*e_4*f_2 - a_6*b_2*c_1*e_4*f_5 + a_6*b_2*c_1*e_5*f_4 + a_6*b_2*c_4*e_1*f_5 - a_6*b_2*c_4*e_5*f_1 - a_6*b_2*c_5*e_1*f_4 + a_6*b_2*c_5*e_4*f_1 + a_6*b_4*c_1*e_2*f_5 - a_6*b_4*c_1*e_5*f_2 - a_6*b_4*c_2*e_1*f_5 + a_6*b_4*c_2*e_5*f_1 + a_6*b_4*c_5*e_1*f_2 - a_6*b_4*c_5*e_2*f_1 - a_6*b_5*c_1*e_2*f_4 + a_6*b_5*c_1*e_4*f_2 + a_6*b_5*c_2*e_1*f_4 - a_6*b_5*c_2*e_4*f_1 - a_6*b_5*c_4*e_1*f_2 + a_6*b_5*c_4*e_2*f_1)/martix_value
              Inverse_A(5,3)=(a_1*b_2*c_4*d_5*f_6 - a_1*b_2*c_4*d_6*f_5 - a_1*b_2*c_5*d_4*f_6 + a_1*b_2*c_5*d_6*f_4 + a_1*b_2*c_6*d_4*f_5 - a_1*b_2*c_6*d_5*f_4 - a_1*b_4*c_2*d_5*f_6 + a_1*b_4*c_2*d_6*f_5 + a_1*b_4*c_5*d_2*f_6 - a_1*b_4*c_5*d_6*f_2 - a_1*b_4*c_6*d_2*f_5 + a_1*b_4*c_6*d_5*f_2 + a_1*b_5*c_2*d_4*f_6 - a_1*b_5*c_2*d_6*f_4 - a_1*b_5*c_4*d_2*f_6 + a_1*b_5*c_4*d_6*f_2 + a_1*b_5*c_6*d_2*f_4 - a_1*b_5*c_6*d_4*f_2 - a_1*b_6*c_2*d_4*f_5 + a_1*b_6*c_2*d_5*f_4 + a_1*b_6*c_4*d_2*f_5 - a_1*b_6*c_4*d_5*f_2 - a_1*b_6*c_5*d_2*f_4 + a_1*b_6*c_5*d_4*f_2 - a_2*b_1*c_4*d_5*f_6 + a_2*b_1*c_4*d_6*f_5 + a_2*b_1*c_5*d_4*f_6 - a_2*b_1*c_5*d_6*f_4 - a_2*b_1*c_6*d_4*f_5 + a_2*b_1*c_6*d_5*f_4 + a_2*b_4*c_1*d_5*f_6 - a_2*b_4*c_1*d_6*f_5 - a_2*b_4*c_5*d_1*f_6 + a_2*b_4*c_5*d_6*f_1 + a_2*b_4*c_6*d_1*f_5 - a_2*b_4*c_6*d_5*f_1 - a_2*b_5*c_1*d_4*f_6 + a_2*b_5*c_1*d_6*f_4 + a_2*b_5*c_4*d_1*f_6 - a_2*b_5*c_4*d_6*f_1 - a_2*b_5*c_6*d_1*f_4 + a_2*b_5*c_6*d_4*f_1 + a_2*b_6*c_1*d_4*f_5 - a_2*b_6*c_1*d_5*f_4 - a_2*b_6*c_4*d_1*f_5 + a_2*b_6*c_4*d_5*f_1 + a_2*b_6*c_5*d_1*f_4 - a_2*b_6*c_5*d_4*f_1 + a_4*b_1*c_2*d_5*f_6 - a_4*b_1*c_2*d_6*f_5 - a_4*b_1*c_5*d_2*f_6 + a_4*b_1*c_5*d_6*f_2 + a_4*b_1*c_6*d_2*f_5 - a_4*b_1*c_6*d_5*f_2 - a_4*b_2*c_1*d_5*f_6 + a_4*b_2*c_1*d_6*f_5 + a_4*b_2*c_5*d_1*f_6 - a_4*b_2*c_5*d_6*f_1 - a_4*b_2*c_6*d_1*f_5 + a_4*b_2*c_6*d_5*f_1 + a_4*b_5*c_1*d_2*f_6 - a_4*b_5*c_1*d_6*f_2 - a_4*b_5*c_2*d_1*f_6 + a_4*b_5*c_2*d_6*f_1 + a_4*b_5*c_6*d_1*f_2 - a_4*b_5*c_6*d_2*f_1 - a_4*b_6*c_1*d_2*f_5 + a_4*b_6*c_1*d_5*f_2 + a_4*b_6*c_2*d_1*f_5 - a_4*b_6*c_2*d_5*f_1 - a_4*b_6*c_5*d_1*f_2 + a_4*b_6*c_5*d_2*f_1 - a_5*b_1*c_2*d_4*f_6 + a_5*b_1*c_2*d_6*f_4 + a_5*b_1*c_4*d_2*f_6 - a_5*b_1*c_4*d_6*f_2 - a_5*b_1*c_6*d_2*f_4 + a_5*b_1*c_6*d_4*f_2 + a_5*b_2*c_1*d_4*f_6 - a_5*b_2*c_1*d_6*f_4 - a_5*b_2*c_4*d_1*f_6 + a_5*b_2*c_4*d_6*f_1 + a_5*b_2*c_6*d_1*f_4 - a_5*b_2*c_6*d_4*f_1 - a_5*b_4*c_1*d_2*f_6 + a_5*b_4*c_1*d_6*f_2 + a_5*b_4*c_2*d_1*f_6 - a_5*b_4*c_2*d_6*f_1 - a_5*b_4*c_6*d_1*f_2 + a_5*b_4*c_6*d_2*f_1 + a_5*b_6*c_1*d_2*f_4 - a_5*b_6*c_1*d_4*f_2 - a_5*b_6*c_2*d_1*f_4 + a_5*b_6*c_2*d_4*f_1 + a_5*b_6*c_4*d_1*f_2 - a_5*b_6*c_4*d_2*f_1 + a_6*b_1*c_2*d_4*f_5 - a_6*b_1*c_2*d_5*f_4 - a_6*b_1*c_4*d_2*f_5 + a_6*b_1*c_4*d_5*f_2 + a_6*b_1*c_5*d_2*f_4 - a_6*b_1*c_5*d_4*f_2 - a_6*b_2*c_1*d_4*f_5 + a_6*b_2*c_1*d_5*f_4 + a_6*b_2*c_4*d_1*f_5 - a_6*b_2*c_4*d_5*f_1 - a_6*b_2*c_5*d_1*f_4 + a_6*b_2*c_5*d_4*f_1 + a_6*b_4*c_1*d_2*f_5 - a_6*b_4*c_1*d_5*f_2 - a_6*b_4*c_2*d_1*f_5 + a_6*b_4*c_2*d_5*f_1 + a_6*b_4*c_5*d_1*f_2 - a_6*b_4*c_5*d_2*f_1 - a_6*b_5*c_1*d_2*f_4 + a_6*b_5*c_1*d_4*f_2 + a_6*b_5*c_2*d_1*f_4 - a_6*b_5*c_2*d_4*f_1 - a_6*b_5*c_4*d_1*f_2 + a_6*b_5*c_4*d_2*f_1)/martix_value
              Inverse_A(6,3)=-(a_1*b_2*c_4*d_5*e_6 - a_1*b_2*c_4*d_6*e_5 - a_1*b_2*c_5*d_4*e_6 + a_1*b_2*c_5*d_6*e_4 + a_1*b_2*c_6*d_4*e_5 - a_1*b_2*c_6*d_5*e_4 - a_1*b_4*c_2*d_5*e_6 + a_1*b_4*c_2*d_6*e_5 + a_1*b_4*c_5*d_2*e_6 - a_1*b_4*c_5*d_6*e_2 - a_1*b_4*c_6*d_2*e_5 + a_1*b_4*c_6*d_5*e_2 + a_1*b_5*c_2*d_4*e_6 - a_1*b_5*c_2*d_6*e_4 - a_1*b_5*c_4*d_2*e_6 + a_1*b_5*c_4*d_6*e_2 + a_1*b_5*c_6*d_2*e_4 - a_1*b_5*c_6*d_4*e_2 - a_1*b_6*c_2*d_4*e_5 + a_1*b_6*c_2*d_5*e_4 + a_1*b_6*c_4*d_2*e_5 - a_1*b_6*c_4*d_5*e_2 - a_1*b_6*c_5*d_2*e_4 + a_1*b_6*c_5*d_4*e_2 - a_2*b_1*c_4*d_5*e_6 + a_2*b_1*c_4*d_6*e_5 + a_2*b_1*c_5*d_4*e_6 - a_2*b_1*c_5*d_6*e_4 - a_2*b_1*c_6*d_4*e_5 + a_2*b_1*c_6*d_5*e_4 + a_2*b_4*c_1*d_5*e_6 - a_2*b_4*c_1*d_6*e_5 - a_2*b_4*c_5*d_1*e_6 + a_2*b_4*c_5*d_6*e_1 + a_2*b_4*c_6*d_1*e_5 - a_2*b_4*c_6*d_5*e_1 - a_2*b_5*c_1*d_4*e_6 + a_2*b_5*c_1*d_6*e_4 + a_2*b_5*c_4*d_1*e_6 - a_2*b_5*c_4*d_6*e_1 - a_2*b_5*c_6*d_1*e_4 + a_2*b_5*c_6*d_4*e_1 + a_2*b_6*c_1*d_4*e_5 - a_2*b_6*c_1*d_5*e_4 - a_2*b_6*c_4*d_1*e_5 + a_2*b_6*c_4*d_5*e_1 + a_2*b_6*c_5*d_1*e_4 - a_2*b_6*c_5*d_4*e_1 + a_4*b_1*c_2*d_5*e_6 - a_4*b_1*c_2*d_6*e_5 - a_4*b_1*c_5*d_2*e_6 + a_4*b_1*c_5*d_6*e_2 + a_4*b_1*c_6*d_2*e_5 - a_4*b_1*c_6*d_5*e_2 - a_4*b_2*c_1*d_5*e_6 + a_4*b_2*c_1*d_6*e_5 + a_4*b_2*c_5*d_1*e_6 - a_4*b_2*c_5*d_6*e_1 - a_4*b_2*c_6*d_1*e_5 + a_4*b_2*c_6*d_5*e_1 + a_4*b_5*c_1*d_2*e_6 - a_4*b_5*c_1*d_6*e_2 - a_4*b_5*c_2*d_1*e_6 + a_4*b_5*c_2*d_6*e_1 + a_4*b_5*c_6*d_1*e_2 - a_4*b_5*c_6*d_2*e_1 - a_4*b_6*c_1*d_2*e_5 + a_4*b_6*c_1*d_5*e_2 + a_4*b_6*c_2*d_1*e_5 - a_4*b_6*c_2*d_5*e_1 - a_4*b_6*c_5*d_1*e_2 + a_4*b_6*c_5*d_2*e_1 - a_5*b_1*c_2*d_4*e_6 + a_5*b_1*c_2*d_6*e_4 + a_5*b_1*c_4*d_2*e_6 - a_5*b_1*c_4*d_6*e_2 - a_5*b_1*c_6*d_2*e_4 + a_5*b_1*c_6*d_4*e_2 + a_5*b_2*c_1*d_4*e_6 - a_5*b_2*c_1*d_6*e_4 - a_5*b_2*c_4*d_1*e_6 + a_5*b_2*c_4*d_6*e_1 + a_5*b_2*c_6*d_1*e_4 - a_5*b_2*c_6*d_4*e_1 - a_5*b_4*c_1*d_2*e_6 + a_5*b_4*c_1*d_6*e_2 + a_5*b_4*c_2*d_1*e_6 - a_5*b_4*c_2*d_6*e_1 - a_5*b_4*c_6*d_1*e_2 + a_5*b_4*c_6*d_2*e_1 + a_5*b_6*c_1*d_2*e_4 - a_5*b_6*c_1*d_4*e_2 - a_5*b_6*c_2*d_1*e_4 + a_5*b_6*c_2*d_4*e_1 + a_5*b_6*c_4*d_1*e_2 - a_5*b_6*c_4*d_2*e_1 + a_6*b_1*c_2*d_4*e_5 - a_6*b_1*c_2*d_5*e_4 - a_6*b_1*c_4*d_2*e_5 + a_6*b_1*c_4*d_5*e_2 + a_6*b_1*c_5*d_2*e_4 - a_6*b_1*c_5*d_4*e_2 - a_6*b_2*c_1*d_4*e_5 + a_6*b_2*c_1*d_5*e_4 + a_6*b_2*c_4*d_1*e_5 - a_6*b_2*c_4*d_5*e_1 - a_6*b_2*c_5*d_1*e_4 + a_6*b_2*c_5*d_4*e_1 + a_6*b_4*c_1*d_2*e_5 - a_6*b_4*c_1*d_5*e_2 - a_6*b_4*c_2*d_1*e_5 + a_6*b_4*c_2*d_5*e_1 + a_6*b_4*c_5*d_1*e_2 - a_6*b_4*c_5*d_2*e_1 - a_6*b_5*c_1*d_2*e_4 + a_6*b_5*c_1*d_4*e_2 + a_6*b_5*c_2*d_1*e_4 - a_6*b_5*c_2*d_4*e_1 - a_6*b_5*c_4*d_1*e_2 + a_6*b_5*c_4*d_2*e_1)/martix_value
              
              Inverse_A(1,4)=-(b_1*c_2*d_3*e_5*f_6 - b_1*c_2*d_3*e_6*f_5 - b_1*c_2*d_5*e_3*f_6 + b_1*c_2*d_5*e_6*f_3 + b_1*c_2*d_6*e_3*f_5 - b_1*c_2*d_6*e_5*f_3 - b_1*c_3*d_2*e_5*f_6 + b_1*c_3*d_2*e_6*f_5 + b_1*c_3*d_5*e_2*f_6 - b_1*c_3*d_5*e_6*f_2 - b_1*c_3*d_6*e_2*f_5 + b_1*c_3*d_6*e_5*f_2 + b_1*c_5*d_2*e_3*f_6 - b_1*c_5*d_2*e_6*f_3 - b_1*c_5*d_3*e_2*f_6 + b_1*c_5*d_3*e_6*f_2 + b_1*c_5*d_6*e_2*f_3 - b_1*c_5*d_6*e_3*f_2 - b_1*c_6*d_2*e_3*f_5 + b_1*c_6*d_2*e_5*f_3 + b_1*c_6*d_3*e_2*f_5 - b_1*c_6*d_3*e_5*f_2 - b_1*c_6*d_5*e_2*f_3 + b_1*c_6*d_5*e_3*f_2 - b_2*c_1*d_3*e_5*f_6 + b_2*c_1*d_3*e_6*f_5 + b_2*c_1*d_5*e_3*f_6 - b_2*c_1*d_5*e_6*f_3 - b_2*c_1*d_6*e_3*f_5 + b_2*c_1*d_6*e_5*f_3 + b_2*c_3*d_1*e_5*f_6 - b_2*c_3*d_1*e_6*f_5 - b_2*c_3*d_5*e_1*f_6 + b_2*c_3*d_5*e_6*f_1 + b_2*c_3*d_6*e_1*f_5 - b_2*c_3*d_6*e_5*f_1 - b_2*c_5*d_1*e_3*f_6 + b_2*c_5*d_1*e_6*f_3 + b_2*c_5*d_3*e_1*f_6 - b_2*c_5*d_3*e_6*f_1 - b_2*c_5*d_6*e_1*f_3 + b_2*c_5*d_6*e_3*f_1 + b_2*c_6*d_1*e_3*f_5 - b_2*c_6*d_1*e_5*f_3 - b_2*c_6*d_3*e_1*f_5 + b_2*c_6*d_3*e_5*f_1 + b_2*c_6*d_5*e_1*f_3 - b_2*c_6*d_5*e_3*f_1 + b_3*c_1*d_2*e_5*f_6 - b_3*c_1*d_2*e_6*f_5 - b_3*c_1*d_5*e_2*f_6 + b_3*c_1*d_5*e_6*f_2 + b_3*c_1*d_6*e_2*f_5 - b_3*c_1*d_6*e_5*f_2 - b_3*c_2*d_1*e_5*f_6 + b_3*c_2*d_1*e_6*f_5 + b_3*c_2*d_5*e_1*f_6 - b_3*c_2*d_5*e_6*f_1 - b_3*c_2*d_6*e_1*f_5 + b_3*c_2*d_6*e_5*f_1 + b_3*c_5*d_1*e_2*f_6 - b_3*c_5*d_1*e_6*f_2 - b_3*c_5*d_2*e_1*f_6 + b_3*c_5*d_2*e_6*f_1 + b_3*c_5*d_6*e_1*f_2 - b_3*c_5*d_6*e_2*f_1 - b_3*c_6*d_1*e_2*f_5 + b_3*c_6*d_1*e_5*f_2 + b_3*c_6*d_2*e_1*f_5 - b_3*c_6*d_2*e_5*f_1 - b_3*c_6*d_5*e_1*f_2 + b_3*c_6*d_5*e_2*f_1 - b_5*c_1*d_2*e_3*f_6 + b_5*c_1*d_2*e_6*f_3 + b_5*c_1*d_3*e_2*f_6 - b_5*c_1*d_3*e_6*f_2 - b_5*c_1*d_6*e_2*f_3 + b_5*c_1*d_6*e_3*f_2 + b_5*c_2*d_1*e_3*f_6 - b_5*c_2*d_1*e_6*f_3 - b_5*c_2*d_3*e_1*f_6 + b_5*c_2*d_3*e_6*f_1 + b_5*c_2*d_6*e_1*f_3 - b_5*c_2*d_6*e_3*f_1 - b_5*c_3*d_1*e_2*f_6 + b_5*c_3*d_1*e_6*f_2 + b_5*c_3*d_2*e_1*f_6 - b_5*c_3*d_2*e_6*f_1 - b_5*c_3*d_6*e_1*f_2 + b_5*c_3*d_6*e_2*f_1 + b_5*c_6*d_1*e_2*f_3 - b_5*c_6*d_1*e_3*f_2 - b_5*c_6*d_2*e_1*f_3 + b_5*c_6*d_2*e_3*f_1 + b_5*c_6*d_3*e_1*f_2 - b_5*c_6*d_3*e_2*f_1 + b_6*c_1*d_2*e_3*f_5 - b_6*c_1*d_2*e_5*f_3 - b_6*c_1*d_3*e_2*f_5 + b_6*c_1*d_3*e_5*f_2 + b_6*c_1*d_5*e_2*f_3 - b_6*c_1*d_5*e_3*f_2 - b_6*c_2*d_1*e_3*f_5 + b_6*c_2*d_1*e_5*f_3 + b_6*c_2*d_3*e_1*f_5 - b_6*c_2*d_3*e_5*f_1 - b_6*c_2*d_5*e_1*f_3 + b_6*c_2*d_5*e_3*f_1 + b_6*c_3*d_1*e_2*f_5 - b_6*c_3*d_1*e_5*f_2 - b_6*c_3*d_2*e_1*f_5 + b_6*c_3*d_2*e_5*f_1 + b_6*c_3*d_5*e_1*f_2 - b_6*c_3*d_5*e_2*f_1 - b_6*c_5*d_1*e_2*f_3 + b_6*c_5*d_1*e_3*f_2 + b_6*c_5*d_2*e_1*f_3 - b_6*c_5*d_2*e_3*f_1 - b_6*c_5*d_3*e_1*f_2 + b_6*c_5*d_3*e_2*f_1)/martix_value
              Inverse_A(2,4)=(a_1*c_2*d_3*e_5*f_6 - a_1*c_2*d_3*e_6*f_5 - a_1*c_2*d_5*e_3*f_6 + a_1*c_2*d_5*e_6*f_3 + a_1*c_2*d_6*e_3*f_5 - a_1*c_2*d_6*e_5*f_3 - a_1*c_3*d_2*e_5*f_6 + a_1*c_3*d_2*e_6*f_5 + a_1*c_3*d_5*e_2*f_6 - a_1*c_3*d_5*e_6*f_2 - a_1*c_3*d_6*e_2*f_5 + a_1*c_3*d_6*e_5*f_2 + a_1*c_5*d_2*e_3*f_6 - a_1*c_5*d_2*e_6*f_3 - a_1*c_5*d_3*e_2*f_6 + a_1*c_5*d_3*e_6*f_2 + a_1*c_5*d_6*e_2*f_3 - a_1*c_5*d_6*e_3*f_2 - a_1*c_6*d_2*e_3*f_5 + a_1*c_6*d_2*e_5*f_3 + a_1*c_6*d_3*e_2*f_5 - a_1*c_6*d_3*e_5*f_2 - a_1*c_6*d_5*e_2*f_3 + a_1*c_6*d_5*e_3*f_2 - a_2*c_1*d_3*e_5*f_6 + a_2*c_1*d_3*e_6*f_5 + a_2*c_1*d_5*e_3*f_6 - a_2*c_1*d_5*e_6*f_3 - a_2*c_1*d_6*e_3*f_5 + a_2*c_1*d_6*e_5*f_3 + a_2*c_3*d_1*e_5*f_6 - a_2*c_3*d_1*e_6*f_5 - a_2*c_3*d_5*e_1*f_6 + a_2*c_3*d_5*e_6*f_1 + a_2*c_3*d_6*e_1*f_5 - a_2*c_3*d_6*e_5*f_1 - a_2*c_5*d_1*e_3*f_6 + a_2*c_5*d_1*e_6*f_3 + a_2*c_5*d_3*e_1*f_6 - a_2*c_5*d_3*e_6*f_1 - a_2*c_5*d_6*e_1*f_3 + a_2*c_5*d_6*e_3*f_1 + a_2*c_6*d_1*e_3*f_5 - a_2*c_6*d_1*e_5*f_3 - a_2*c_6*d_3*e_1*f_5 + a_2*c_6*d_3*e_5*f_1 + a_2*c_6*d_5*e_1*f_3 - a_2*c_6*d_5*e_3*f_1 + a_3*c_1*d_2*e_5*f_6 - a_3*c_1*d_2*e_6*f_5 - a_3*c_1*d_5*e_2*f_6 + a_3*c_1*d_5*e_6*f_2 + a_3*c_1*d_6*e_2*f_5 - a_3*c_1*d_6*e_5*f_2 - a_3*c_2*d_1*e_5*f_6 + a_3*c_2*d_1*e_6*f_5 + a_3*c_2*d_5*e_1*f_6 - a_3*c_2*d_5*e_6*f_1 - a_3*c_2*d_6*e_1*f_5 + a_3*c_2*d_6*e_5*f_1 + a_3*c_5*d_1*e_2*f_6 - a_3*c_5*d_1*e_6*f_2 - a_3*c_5*d_2*e_1*f_6 + a_3*c_5*d_2*e_6*f_1 + a_3*c_5*d_6*e_1*f_2 - a_3*c_5*d_6*e_2*f_1 - a_3*c_6*d_1*e_2*f_5 + a_3*c_6*d_1*e_5*f_2 + a_3*c_6*d_2*e_1*f_5 - a_3*c_6*d_2*e_5*f_1 - a_3*c_6*d_5*e_1*f_2 + a_3*c_6*d_5*e_2*f_1 - a_5*c_1*d_2*e_3*f_6 + a_5*c_1*d_2*e_6*f_3 + a_5*c_1*d_3*e_2*f_6 - a_5*c_1*d_3*e_6*f_2 - a_5*c_1*d_6*e_2*f_3 + a_5*c_1*d_6*e_3*f_2 + a_5*c_2*d_1*e_3*f_6 - a_5*c_2*d_1*e_6*f_3 - a_5*c_2*d_3*e_1*f_6 + a_5*c_2*d_3*e_6*f_1 + a_5*c_2*d_6*e_1*f_3 - a_5*c_2*d_6*e_3*f_1 - a_5*c_3*d_1*e_2*f_6 + a_5*c_3*d_1*e_6*f_2 + a_5*c_3*d_2*e_1*f_6 - a_5*c_3*d_2*e_6*f_1 - a_5*c_3*d_6*e_1*f_2 + a_5*c_3*d_6*e_2*f_1 + a_5*c_6*d_1*e_2*f_3 - a_5*c_6*d_1*e_3*f_2 - a_5*c_6*d_2*e_1*f_3 + a_5*c_6*d_2*e_3*f_1 + a_5*c_6*d_3*e_1*f_2 - a_5*c_6*d_3*e_2*f_1 + a_6*c_1*d_2*e_3*f_5 - a_6*c_1*d_2*e_5*f_3 - a_6*c_1*d_3*e_2*f_5 + a_6*c_1*d_3*e_5*f_2 + a_6*c_1*d_5*e_2*f_3 - a_6*c_1*d_5*e_3*f_2 - a_6*c_2*d_1*e_3*f_5 + a_6*c_2*d_1*e_5*f_3 + a_6*c_2*d_3*e_1*f_5 - a_6*c_2*d_3*e_5*f_1 - a_6*c_2*d_5*e_1*f_3 + a_6*c_2*d_5*e_3*f_1 + a_6*c_3*d_1*e_2*f_5 - a_6*c_3*d_1*e_5*f_2 - a_6*c_3*d_2*e_1*f_5 + a_6*c_3*d_2*e_5*f_1 + a_6*c_3*d_5*e_1*f_2 - a_6*c_3*d_5*e_2*f_1 - a_6*c_5*d_1*e_2*f_3 + a_6*c_5*d_1*e_3*f_2 + a_6*c_5*d_2*e_1*f_3 - a_6*c_5*d_2*e_3*f_1 - a_6*c_5*d_3*e_1*f_2 + a_6*c_5*d_3*e_2*f_1)/martix_value
              Inverse_A(3,4)=-(a_1*b_2*d_3*e_5*f_6 - a_1*b_2*d_3*e_6*f_5 - a_1*b_2*d_5*e_3*f_6 + a_1*b_2*d_5*e_6*f_3 + a_1*b_2*d_6*e_3*f_5 - a_1*b_2*d_6*e_5*f_3 - a_1*b_3*d_2*e_5*f_6 + a_1*b_3*d_2*e_6*f_5 + a_1*b_3*d_5*e_2*f_6 - a_1*b_3*d_5*e_6*f_2 - a_1*b_3*d_6*e_2*f_5 + a_1*b_3*d_6*e_5*f_2 + a_1*b_5*d_2*e_3*f_6 - a_1*b_5*d_2*e_6*f_3 - a_1*b_5*d_3*e_2*f_6 + a_1*b_5*d_3*e_6*f_2 + a_1*b_5*d_6*e_2*f_3 - a_1*b_5*d_6*e_3*f_2 - a_1*b_6*d_2*e_3*f_5 + a_1*b_6*d_2*e_5*f_3 + a_1*b_6*d_3*e_2*f_5 - a_1*b_6*d_3*e_5*f_2 - a_1*b_6*d_5*e_2*f_3 + a_1*b_6*d_5*e_3*f_2 - a_2*b_1*d_3*e_5*f_6 + a_2*b_1*d_3*e_6*f_5 + a_2*b_1*d_5*e_3*f_6 - a_2*b_1*d_5*e_6*f_3 - a_2*b_1*d_6*e_3*f_5 + a_2*b_1*d_6*e_5*f_3 + a_2*b_3*d_1*e_5*f_6 - a_2*b_3*d_1*e_6*f_5 - a_2*b_3*d_5*e_1*f_6 + a_2*b_3*d_5*e_6*f_1 + a_2*b_3*d_6*e_1*f_5 - a_2*b_3*d_6*e_5*f_1 - a_2*b_5*d_1*e_3*f_6 + a_2*b_5*d_1*e_6*f_3 + a_2*b_5*d_3*e_1*f_6 - a_2*b_5*d_3*e_6*f_1 - a_2*b_5*d_6*e_1*f_3 + a_2*b_5*d_6*e_3*f_1 + a_2*b_6*d_1*e_3*f_5 - a_2*b_6*d_1*e_5*f_3 - a_2*b_6*d_3*e_1*f_5 + a_2*b_6*d_3*e_5*f_1 + a_2*b_6*d_5*e_1*f_3 - a_2*b_6*d_5*e_3*f_1 + a_3*b_1*d_2*e_5*f_6 - a_3*b_1*d_2*e_6*f_5 - a_3*b_1*d_5*e_2*f_6 + a_3*b_1*d_5*e_6*f_2 + a_3*b_1*d_6*e_2*f_5 - a_3*b_1*d_6*e_5*f_2 - a_3*b_2*d_1*e_5*f_6 + a_3*b_2*d_1*e_6*f_5 + a_3*b_2*d_5*e_1*f_6 - a_3*b_2*d_5*e_6*f_1 - a_3*b_2*d_6*e_1*f_5 + a_3*b_2*d_6*e_5*f_1 + a_3*b_5*d_1*e_2*f_6 - a_3*b_5*d_1*e_6*f_2 - a_3*b_5*d_2*e_1*f_6 + a_3*b_5*d_2*e_6*f_1 + a_3*b_5*d_6*e_1*f_2 - a_3*b_5*d_6*e_2*f_1 - a_3*b_6*d_1*e_2*f_5 + a_3*b_6*d_1*e_5*f_2 + a_3*b_6*d_2*e_1*f_5 - a_3*b_6*d_2*e_5*f_1 - a_3*b_6*d_5*e_1*f_2 + a_3*b_6*d_5*e_2*f_1 - a_5*b_1*d_2*e_3*f_6 + a_5*b_1*d_2*e_6*f_3 + a_5*b_1*d_3*e_2*f_6 - a_5*b_1*d_3*e_6*f_2 - a_5*b_1*d_6*e_2*f_3 + a_5*b_1*d_6*e_3*f_2 + a_5*b_2*d_1*e_3*f_6 - a_5*b_2*d_1*e_6*f_3 - a_5*b_2*d_3*e_1*f_6 + a_5*b_2*d_3*e_6*f_1 + a_5*b_2*d_6*e_1*f_3 - a_5*b_2*d_6*e_3*f_1 - a_5*b_3*d_1*e_2*f_6 + a_5*b_3*d_1*e_6*f_2 + a_5*b_3*d_2*e_1*f_6 - a_5*b_3*d_2*e_6*f_1 - a_5*b_3*d_6*e_1*f_2 + a_5*b_3*d_6*e_2*f_1 + a_5*b_6*d_1*e_2*f_3 - a_5*b_6*d_1*e_3*f_2 - a_5*b_6*d_2*e_1*f_3 + a_5*b_6*d_2*e_3*f_1 + a_5*b_6*d_3*e_1*f_2 - a_5*b_6*d_3*e_2*f_1 + a_6*b_1*d_2*e_3*f_5 - a_6*b_1*d_2*e_5*f_3 - a_6*b_1*d_3*e_2*f_5 + a_6*b_1*d_3*e_5*f_2 + a_6*b_1*d_5*e_2*f_3 - a_6*b_1*d_5*e_3*f_2 - a_6*b_2*d_1*e_3*f_5 + a_6*b_2*d_1*e_5*f_3 + a_6*b_2*d_3*e_1*f_5 - a_6*b_2*d_3*e_5*f_1 - a_6*b_2*d_5*e_1*f_3 + a_6*b_2*d_5*e_3*f_1 + a_6*b_3*d_1*e_2*f_5 - a_6*b_3*d_1*e_5*f_2 - a_6*b_3*d_2*e_1*f_5 + a_6*b_3*d_2*e_5*f_1 + a_6*b_3*d_5*e_1*f_2 - a_6*b_3*d_5*e_2*f_1 - a_6*b_5*d_1*e_2*f_3 + a_6*b_5*d_1*e_3*f_2 + a_6*b_5*d_2*e_1*f_3 - a_6*b_5*d_2*e_3*f_1 - a_6*b_5*d_3*e_1*f_2 + a_6*b_5*d_3*e_2*f_1)/martix_value
              Inverse_A(4,4)=(a_1*b_2*c_3*e_5*f_6 - a_1*b_2*c_3*e_6*f_5 - a_1*b_2*c_5*e_3*f_6 + a_1*b_2*c_5*e_6*f_3 + a_1*b_2*c_6*e_3*f_5 - a_1*b_2*c_6*e_5*f_3 - a_1*b_3*c_2*e_5*f_6 + a_1*b_3*c_2*e_6*f_5 + a_1*b_3*c_5*e_2*f_6 - a_1*b_3*c_5*e_6*f_2 - a_1*b_3*c_6*e_2*f_5 + a_1*b_3*c_6*e_5*f_2 + a_1*b_5*c_2*e_3*f_6 - a_1*b_5*c_2*e_6*f_3 - a_1*b_5*c_3*e_2*f_6 + a_1*b_5*c_3*e_6*f_2 + a_1*b_5*c_6*e_2*f_3 - a_1*b_5*c_6*e_3*f_2 - a_1*b_6*c_2*e_3*f_5 + a_1*b_6*c_2*e_5*f_3 + a_1*b_6*c_3*e_2*f_5 - a_1*b_6*c_3*e_5*f_2 - a_1*b_6*c_5*e_2*f_3 + a_1*b_6*c_5*e_3*f_2 - a_2*b_1*c_3*e_5*f_6 + a_2*b_1*c_3*e_6*f_5 + a_2*b_1*c_5*e_3*f_6 - a_2*b_1*c_5*e_6*f_3 - a_2*b_1*c_6*e_3*f_5 + a_2*b_1*c_6*e_5*f_3 + a_2*b_3*c_1*e_5*f_6 - a_2*b_3*c_1*e_6*f_5 - a_2*b_3*c_5*e_1*f_6 + a_2*b_3*c_5*e_6*f_1 + a_2*b_3*c_6*e_1*f_5 - a_2*b_3*c_6*e_5*f_1 - a_2*b_5*c_1*e_3*f_6 + a_2*b_5*c_1*e_6*f_3 + a_2*b_5*c_3*e_1*f_6 - a_2*b_5*c_3*e_6*f_1 - a_2*b_5*c_6*e_1*f_3 + a_2*b_5*c_6*e_3*f_1 + a_2*b_6*c_1*e_3*f_5 - a_2*b_6*c_1*e_5*f_3 - a_2*b_6*c_3*e_1*f_5 + a_2*b_6*c_3*e_5*f_1 + a_2*b_6*c_5*e_1*f_3 - a_2*b_6*c_5*e_3*f_1 + a_3*b_1*c_2*e_5*f_6 - a_3*b_1*c_2*e_6*f_5 - a_3*b_1*c_5*e_2*f_6 + a_3*b_1*c_5*e_6*f_2 + a_3*b_1*c_6*e_2*f_5 - a_3*b_1*c_6*e_5*f_2 - a_3*b_2*c_1*e_5*f_6 + a_3*b_2*c_1*e_6*f_5 + a_3*b_2*c_5*e_1*f_6 - a_3*b_2*c_5*e_6*f_1 - a_3*b_2*c_6*e_1*f_5 + a_3*b_2*c_6*e_5*f_1 + a_3*b_5*c_1*e_2*f_6 - a_3*b_5*c_1*e_6*f_2 - a_3*b_5*c_2*e_1*f_6 + a_3*b_5*c_2*e_6*f_1 + a_3*b_5*c_6*e_1*f_2 - a_3*b_5*c_6*e_2*f_1 - a_3*b_6*c_1*e_2*f_5 + a_3*b_6*c_1*e_5*f_2 + a_3*b_6*c_2*e_1*f_5 - a_3*b_6*c_2*e_5*f_1 - a_3*b_6*c_5*e_1*f_2 + a_3*b_6*c_5*e_2*f_1 - a_5*b_1*c_2*e_3*f_6 + a_5*b_1*c_2*e_6*f_3 + a_5*b_1*c_3*e_2*f_6 - a_5*b_1*c_3*e_6*f_2 - a_5*b_1*c_6*e_2*f_3 + a_5*b_1*c_6*e_3*f_2 + a_5*b_2*c_1*e_3*f_6 - a_5*b_2*c_1*e_6*f_3 - a_5*b_2*c_3*e_1*f_6 + a_5*b_2*c_3*e_6*f_1 + a_5*b_2*c_6*e_1*f_3 - a_5*b_2*c_6*e_3*f_1 - a_5*b_3*c_1*e_2*f_6 + a_5*b_3*c_1*e_6*f_2 + a_5*b_3*c_2*e_1*f_6 - a_5*b_3*c_2*e_6*f_1 - a_5*b_3*c_6*e_1*f_2 + a_5*b_3*c_6*e_2*f_1 + a_5*b_6*c_1*e_2*f_3 - a_5*b_6*c_1*e_3*f_2 - a_5*b_6*c_2*e_1*f_3 + a_5*b_6*c_2*e_3*f_1 + a_5*b_6*c_3*e_1*f_2 - a_5*b_6*c_3*e_2*f_1 + a_6*b_1*c_2*e_3*f_5 - a_6*b_1*c_2*e_5*f_3 - a_6*b_1*c_3*e_2*f_5 + a_6*b_1*c_3*e_5*f_2 + a_6*b_1*c_5*e_2*f_3 - a_6*b_1*c_5*e_3*f_2 - a_6*b_2*c_1*e_3*f_5 + a_6*b_2*c_1*e_5*f_3 + a_6*b_2*c_3*e_1*f_5 - a_6*b_2*c_3*e_5*f_1 - a_6*b_2*c_5*e_1*f_3 + a_6*b_2*c_5*e_3*f_1 + a_6*b_3*c_1*e_2*f_5 - a_6*b_3*c_1*e_5*f_2 - a_6*b_3*c_2*e_1*f_5 + a_6*b_3*c_2*e_5*f_1 + a_6*b_3*c_5*e_1*f_2 - a_6*b_3*c_5*e_2*f_1 - a_6*b_5*c_1*e_2*f_3 + a_6*b_5*c_1*e_3*f_2 + a_6*b_5*c_2*e_1*f_3 - a_6*b_5*c_2*e_3*f_1 - a_6*b_5*c_3*e_1*f_2 + a_6*b_5*c_3*e_2*f_1)/martix_value
              Inverse_A(5,4)=-(a_1*b_2*c_3*d_5*f_6 - a_1*b_2*c_3*d_6*f_5 - a_1*b_2*c_5*d_3*f_6 + a_1*b_2*c_5*d_6*f_3 + a_1*b_2*c_6*d_3*f_5 - a_1*b_2*c_6*d_5*f_3 - a_1*b_3*c_2*d_5*f_6 + a_1*b_3*c_2*d_6*f_5 + a_1*b_3*c_5*d_2*f_6 - a_1*b_3*c_5*d_6*f_2 - a_1*b_3*c_6*d_2*f_5 + a_1*b_3*c_6*d_5*f_2 + a_1*b_5*c_2*d_3*f_6 - a_1*b_5*c_2*d_6*f_3 - a_1*b_5*c_3*d_2*f_6 + a_1*b_5*c_3*d_6*f_2 + a_1*b_5*c_6*d_2*f_3 - a_1*b_5*c_6*d_3*f_2 - a_1*b_6*c_2*d_3*f_5 + a_1*b_6*c_2*d_5*f_3 + a_1*b_6*c_3*d_2*f_5 - a_1*b_6*c_3*d_5*f_2 - a_1*b_6*c_5*d_2*f_3 + a_1*b_6*c_5*d_3*f_2 - a_2*b_1*c_3*d_5*f_6 + a_2*b_1*c_3*d_6*f_5 + a_2*b_1*c_5*d_3*f_6 - a_2*b_1*c_5*d_6*f_3 - a_2*b_1*c_6*d_3*f_5 + a_2*b_1*c_6*d_5*f_3 + a_2*b_3*c_1*d_5*f_6 - a_2*b_3*c_1*d_6*f_5 - a_2*b_3*c_5*d_1*f_6 + a_2*b_3*c_5*d_6*f_1 + a_2*b_3*c_6*d_1*f_5 - a_2*b_3*c_6*d_5*f_1 - a_2*b_5*c_1*d_3*f_6 + a_2*b_5*c_1*d_6*f_3 + a_2*b_5*c_3*d_1*f_6 - a_2*b_5*c_3*d_6*f_1 - a_2*b_5*c_6*d_1*f_3 + a_2*b_5*c_6*d_3*f_1 + a_2*b_6*c_1*d_3*f_5 - a_2*b_6*c_1*d_5*f_3 - a_2*b_6*c_3*d_1*f_5 + a_2*b_6*c_3*d_5*f_1 + a_2*b_6*c_5*d_1*f_3 - a_2*b_6*c_5*d_3*f_1 + a_3*b_1*c_2*d_5*f_6 - a_3*b_1*c_2*d_6*f_5 - a_3*b_1*c_5*d_2*f_6 + a_3*b_1*c_5*d_6*f_2 + a_3*b_1*c_6*d_2*f_5 - a_3*b_1*c_6*d_5*f_2 - a_3*b_2*c_1*d_5*f_6 + a_3*b_2*c_1*d_6*f_5 + a_3*b_2*c_5*d_1*f_6 - a_3*b_2*c_5*d_6*f_1 - a_3*b_2*c_6*d_1*f_5 + a_3*b_2*c_6*d_5*f_1 + a_3*b_5*c_1*d_2*f_6 - a_3*b_5*c_1*d_6*f_2 - a_3*b_5*c_2*d_1*f_6 + a_3*b_5*c_2*d_6*f_1 + a_3*b_5*c_6*d_1*f_2 - a_3*b_5*c_6*d_2*f_1 - a_3*b_6*c_1*d_2*f_5 + a_3*b_6*c_1*d_5*f_2 + a_3*b_6*c_2*d_1*f_5 - a_3*b_6*c_2*d_5*f_1 - a_3*b_6*c_5*d_1*f_2 + a_3*b_6*c_5*d_2*f_1 - a_5*b_1*c_2*d_3*f_6 + a_5*b_1*c_2*d_6*f_3 + a_5*b_1*c_3*d_2*f_6 - a_5*b_1*c_3*d_6*f_2 - a_5*b_1*c_6*d_2*f_3 + a_5*b_1*c_6*d_3*f_2 + a_5*b_2*c_1*d_3*f_6 - a_5*b_2*c_1*d_6*f_3 - a_5*b_2*c_3*d_1*f_6 + a_5*b_2*c_3*d_6*f_1 + a_5*b_2*c_6*d_1*f_3 - a_5*b_2*c_6*d_3*f_1 - a_5*b_3*c_1*d_2*f_6 + a_5*b_3*c_1*d_6*f_2 + a_5*b_3*c_2*d_1*f_6 - a_5*b_3*c_2*d_6*f_1 - a_5*b_3*c_6*d_1*f_2 + a_5*b_3*c_6*d_2*f_1 + a_5*b_6*c_1*d_2*f_3 - a_5*b_6*c_1*d_3*f_2 - a_5*b_6*c_2*d_1*f_3 + a_5*b_6*c_2*d_3*f_1 + a_5*b_6*c_3*d_1*f_2 - a_5*b_6*c_3*d_2*f_1 + a_6*b_1*c_2*d_3*f_5 - a_6*b_1*c_2*d_5*f_3 - a_6*b_1*c_3*d_2*f_5 + a_6*b_1*c_3*d_5*f_2 + a_6*b_1*c_5*d_2*f_3 - a_6*b_1*c_5*d_3*f_2 - a_6*b_2*c_1*d_3*f_5 + a_6*b_2*c_1*d_5*f_3 + a_6*b_2*c_3*d_1*f_5 - a_6*b_2*c_3*d_5*f_1 - a_6*b_2*c_5*d_1*f_3 + a_6*b_2*c_5*d_3*f_1 + a_6*b_3*c_1*d_2*f_5 - a_6*b_3*c_1*d_5*f_2 - a_6*b_3*c_2*d_1*f_5 + a_6*b_3*c_2*d_5*f_1 + a_6*b_3*c_5*d_1*f_2 - a_6*b_3*c_5*d_2*f_1 - a_6*b_5*c_1*d_2*f_3 + a_6*b_5*c_1*d_3*f_2 + a_6*b_5*c_2*d_1*f_3 - a_6*b_5*c_2*d_3*f_1 - a_6*b_5*c_3*d_1*f_2 + a_6*b_5*c_3*d_2*f_1)/martix_value
              Inverse_A(6,4)=(a_1*b_2*c_3*d_5*e_6 - a_1*b_2*c_3*d_6*e_5 - a_1*b_2*c_5*d_3*e_6 + a_1*b_2*c_5*d_6*e_3 + a_1*b_2*c_6*d_3*e_5 - a_1*b_2*c_6*d_5*e_3 - a_1*b_3*c_2*d_5*e_6 + a_1*b_3*c_2*d_6*e_5 + a_1*b_3*c_5*d_2*e_6 - a_1*b_3*c_5*d_6*e_2 - a_1*b_3*c_6*d_2*e_5 + a_1*b_3*c_6*d_5*e_2 + a_1*b_5*c_2*d_3*e_6 - a_1*b_5*c_2*d_6*e_3 - a_1*b_5*c_3*d_2*e_6 + a_1*b_5*c_3*d_6*e_2 + a_1*b_5*c_6*d_2*e_3 - a_1*b_5*c_6*d_3*e_2 - a_1*b_6*c_2*d_3*e_5 + a_1*b_6*c_2*d_5*e_3 + a_1*b_6*c_3*d_2*e_5 - a_1*b_6*c_3*d_5*e_2 - a_1*b_6*c_5*d_2*e_3 + a_1*b_6*c_5*d_3*e_2 - a_2*b_1*c_3*d_5*e_6 + a_2*b_1*c_3*d_6*e_5 + a_2*b_1*c_5*d_3*e_6 - a_2*b_1*c_5*d_6*e_3 - a_2*b_1*c_6*d_3*e_5 + a_2*b_1*c_6*d_5*e_3 + a_2*b_3*c_1*d_5*e_6 - a_2*b_3*c_1*d_6*e_5 - a_2*b_3*c_5*d_1*e_6 + a_2*b_3*c_5*d_6*e_1 + a_2*b_3*c_6*d_1*e_5 - a_2*b_3*c_6*d_5*e_1 - a_2*b_5*c_1*d_3*e_6 + a_2*b_5*c_1*d_6*e_3 + a_2*b_5*c_3*d_1*e_6 - a_2*b_5*c_3*d_6*e_1 - a_2*b_5*c_6*d_1*e_3 + a_2*b_5*c_6*d_3*e_1 + a_2*b_6*c_1*d_3*e_5 - a_2*b_6*c_1*d_5*e_3 - a_2*b_6*c_3*d_1*e_5 + a_2*b_6*c_3*d_5*e_1 + a_2*b_6*c_5*d_1*e_3 - a_2*b_6*c_5*d_3*e_1 + a_3*b_1*c_2*d_5*e_6 - a_3*b_1*c_2*d_6*e_5 - a_3*b_1*c_5*d_2*e_6 + a_3*b_1*c_5*d_6*e_2 + a_3*b_1*c_6*d_2*e_5 - a_3*b_1*c_6*d_5*e_2 - a_3*b_2*c_1*d_5*e_6 + a_3*b_2*c_1*d_6*e_5 + a_3*b_2*c_5*d_1*e_6 - a_3*b_2*c_5*d_6*e_1 - a_3*b_2*c_6*d_1*e_5 + a_3*b_2*c_6*d_5*e_1 + a_3*b_5*c_1*d_2*e_6 - a_3*b_5*c_1*d_6*e_2 - a_3*b_5*c_2*d_1*e_6 + a_3*b_5*c_2*d_6*e_1 + a_3*b_5*c_6*d_1*e_2 - a_3*b_5*c_6*d_2*e_1 - a_3*b_6*c_1*d_2*e_5 + a_3*b_6*c_1*d_5*e_2 + a_3*b_6*c_2*d_1*e_5 - a_3*b_6*c_2*d_5*e_1 - a_3*b_6*c_5*d_1*e_2 + a_3*b_6*c_5*d_2*e_1 - a_5*b_1*c_2*d_3*e_6 + a_5*b_1*c_2*d_6*e_3 + a_5*b_1*c_3*d_2*e_6 - a_5*b_1*c_3*d_6*e_2 - a_5*b_1*c_6*d_2*e_3 + a_5*b_1*c_6*d_3*e_2 + a_5*b_2*c_1*d_3*e_6 - a_5*b_2*c_1*d_6*e_3 - a_5*b_2*c_3*d_1*e_6 + a_5*b_2*c_3*d_6*e_1 + a_5*b_2*c_6*d_1*e_3 - a_5*b_2*c_6*d_3*e_1 - a_5*b_3*c_1*d_2*e_6 + a_5*b_3*c_1*d_6*e_2 + a_5*b_3*c_2*d_1*e_6 - a_5*b_3*c_2*d_6*e_1 - a_5*b_3*c_6*d_1*e_2 + a_5*b_3*c_6*d_2*e_1 + a_5*b_6*c_1*d_2*e_3 - a_5*b_6*c_1*d_3*e_2 - a_5*b_6*c_2*d_1*e_3 + a_5*b_6*c_2*d_3*e_1 + a_5*b_6*c_3*d_1*e_2 - a_5*b_6*c_3*d_2*e_1 + a_6*b_1*c_2*d_3*e_5 - a_6*b_1*c_2*d_5*e_3 - a_6*b_1*c_3*d_2*e_5 + a_6*b_1*c_3*d_5*e_2 + a_6*b_1*c_5*d_2*e_3 - a_6*b_1*c_5*d_3*e_2 - a_6*b_2*c_1*d_3*e_5 + a_6*b_2*c_1*d_5*e_3 + a_6*b_2*c_3*d_1*e_5 - a_6*b_2*c_3*d_5*e_1 - a_6*b_2*c_5*d_1*e_3 + a_6*b_2*c_5*d_3*e_1 + a_6*b_3*c_1*d_2*e_5 - a_6*b_3*c_1*d_5*e_2 - a_6*b_3*c_2*d_1*e_5 + a_6*b_3*c_2*d_5*e_1 + a_6*b_3*c_5*d_1*e_2 - a_6*b_3*c_5*d_2*e_1 - a_6*b_5*c_1*d_2*e_3 + a_6*b_5*c_1*d_3*e_2 + a_6*b_5*c_2*d_1*e_3 - a_6*b_5*c_2*d_3*e_1 - a_6*b_5*c_3*d_1*e_2 + a_6*b_5*c_3*d_2*e_1)/martix_value
              
              Inverse_A(1,5)=(b_1*c_2*d_3*e_4*f_6 - b_1*c_2*d_3*e_6*f_4 - b_1*c_2*d_4*e_3*f_6 + b_1*c_2*d_4*e_6*f_3 + b_1*c_2*d_6*e_3*f_4 - b_1*c_2*d_6*e_4*f_3 - b_1*c_3*d_2*e_4*f_6 + b_1*c_3*d_2*e_6*f_4 + b_1*c_3*d_4*e_2*f_6 - b_1*c_3*d_4*e_6*f_2 - b_1*c_3*d_6*e_2*f_4 + b_1*c_3*d_6*e_4*f_2 + b_1*c_4*d_2*e_3*f_6 - b_1*c_4*d_2*e_6*f_3 - b_1*c_4*d_3*e_2*f_6 + b_1*c_4*d_3*e_6*f_2 + b_1*c_4*d_6*e_2*f_3 - b_1*c_4*d_6*e_3*f_2 - b_1*c_6*d_2*e_3*f_4 + b_1*c_6*d_2*e_4*f_3 + b_1*c_6*d_3*e_2*f_4 - b_1*c_6*d_3*e_4*f_2 - b_1*c_6*d_4*e_2*f_3 + b_1*c_6*d_4*e_3*f_2 - b_2*c_1*d_3*e_4*f_6 + b_2*c_1*d_3*e_6*f_4 + b_2*c_1*d_4*e_3*f_6 - b_2*c_1*d_4*e_6*f_3 - b_2*c_1*d_6*e_3*f_4 + b_2*c_1*d_6*e_4*f_3 + b_2*c_3*d_1*e_4*f_6 - b_2*c_3*d_1*e_6*f_4 - b_2*c_3*d_4*e_1*f_6 + b_2*c_3*d_4*e_6*f_1 + b_2*c_3*d_6*e_1*f_4 - b_2*c_3*d_6*e_4*f_1 - b_2*c_4*d_1*e_3*f_6 + b_2*c_4*d_1*e_6*f_3 + b_2*c_4*d_3*e_1*f_6 - b_2*c_4*d_3*e_6*f_1 - b_2*c_4*d_6*e_1*f_3 + b_2*c_4*d_6*e_3*f_1 + b_2*c_6*d_1*e_3*f_4 - b_2*c_6*d_1*e_4*f_3 - b_2*c_6*d_3*e_1*f_4 + b_2*c_6*d_3*e_4*f_1 + b_2*c_6*d_4*e_1*f_3 - b_2*c_6*d_4*e_3*f_1 + b_3*c_1*d_2*e_4*f_6 - b_3*c_1*d_2*e_6*f_4 - b_3*c_1*d_4*e_2*f_6 + b_3*c_1*d_4*e_6*f_2 + b_3*c_1*d_6*e_2*f_4 - b_3*c_1*d_6*e_4*f_2 - b_3*c_2*d_1*e_4*f_6 + b_3*c_2*d_1*e_6*f_4 + b_3*c_2*d_4*e_1*f_6 - b_3*c_2*d_4*e_6*f_1 - b_3*c_2*d_6*e_1*f_4 + b_3*c_2*d_6*e_4*f_1 + b_3*c_4*d_1*e_2*f_6 - b_3*c_4*d_1*e_6*f_2 - b_3*c_4*d_2*e_1*f_6 + b_3*c_4*d_2*e_6*f_1 + b_3*c_4*d_6*e_1*f_2 - b_3*c_4*d_6*e_2*f_1 - b_3*c_6*d_1*e_2*f_4 + b_3*c_6*d_1*e_4*f_2 + b_3*c_6*d_2*e_1*f_4 - b_3*c_6*d_2*e_4*f_1 - b_3*c_6*d_4*e_1*f_2 + b_3*c_6*d_4*e_2*f_1 - b_4*c_1*d_2*e_3*f_6 + b_4*c_1*d_2*e_6*f_3 + b_4*c_1*d_3*e_2*f_6 - b_4*c_1*d_3*e_6*f_2 - b_4*c_1*d_6*e_2*f_3 + b_4*c_1*d_6*e_3*f_2 + b_4*c_2*d_1*e_3*f_6 - b_4*c_2*d_1*e_6*f_3 - b_4*c_2*d_3*e_1*f_6 + b_4*c_2*d_3*e_6*f_1 + b_4*c_2*d_6*e_1*f_3 - b_4*c_2*d_6*e_3*f_1 - b_4*c_3*d_1*e_2*f_6 + b_4*c_3*d_1*e_6*f_2 + b_4*c_3*d_2*e_1*f_6 - b_4*c_3*d_2*e_6*f_1 - b_4*c_3*d_6*e_1*f_2 + b_4*c_3*d_6*e_2*f_1 + b_4*c_6*d_1*e_2*f_3 - b_4*c_6*d_1*e_3*f_2 - b_4*c_6*d_2*e_1*f_3 + b_4*c_6*d_2*e_3*f_1 + b_4*c_6*d_3*e_1*f_2 - b_4*c_6*d_3*e_2*f_1 + b_6*c_1*d_2*e_3*f_4 - b_6*c_1*d_2*e_4*f_3 - b_6*c_1*d_3*e_2*f_4 + b_6*c_1*d_3*e_4*f_2 + b_6*c_1*d_4*e_2*f_3 - b_6*c_1*d_4*e_3*f_2 - b_6*c_2*d_1*e_3*f_4 + b_6*c_2*d_1*e_4*f_3 + b_6*c_2*d_3*e_1*f_4 - b_6*c_2*d_3*e_4*f_1 - b_6*c_2*d_4*e_1*f_3 + b_6*c_2*d_4*e_3*f_1 + b_6*c_3*d_1*e_2*f_4 - b_6*c_3*d_1*e_4*f_2 - b_6*c_3*d_2*e_1*f_4 + b_6*c_3*d_2*e_4*f_1 + b_6*c_3*d_4*e_1*f_2 - b_6*c_3*d_4*e_2*f_1 - b_6*c_4*d_1*e_2*f_3 + b_6*c_4*d_1*e_3*f_2 + b_6*c_4*d_2*e_1*f_3 - b_6*c_4*d_2*e_3*f_1 - b_6*c_4*d_3*e_1*f_2 + b_6*c_4*d_3*e_2*f_1)/martix_value
              Inverse_A(2,5)=-(a_1*c_2*d_3*e_4*f_6 - a_1*c_2*d_3*e_6*f_4 - a_1*c_2*d_4*e_3*f_6 + a_1*c_2*d_4*e_6*f_3 + a_1*c_2*d_6*e_3*f_4 - a_1*c_2*d_6*e_4*f_3 - a_1*c_3*d_2*e_4*f_6 + a_1*c_3*d_2*e_6*f_4 + a_1*c_3*d_4*e_2*f_6 - a_1*c_3*d_4*e_6*f_2 - a_1*c_3*d_6*e_2*f_4 + a_1*c_3*d_6*e_4*f_2 + a_1*c_4*d_2*e_3*f_6 - a_1*c_4*d_2*e_6*f_3 - a_1*c_4*d_3*e_2*f_6 + a_1*c_4*d_3*e_6*f_2 + a_1*c_4*d_6*e_2*f_3 - a_1*c_4*d_6*e_3*f_2 - a_1*c_6*d_2*e_3*f_4 + a_1*c_6*d_2*e_4*f_3 + a_1*c_6*d_3*e_2*f_4 - a_1*c_6*d_3*e_4*f_2 - a_1*c_6*d_4*e_2*f_3 + a_1*c_6*d_4*e_3*f_2 - a_2*c_1*d_3*e_4*f_6 + a_2*c_1*d_3*e_6*f_4 + a_2*c_1*d_4*e_3*f_6 - a_2*c_1*d_4*e_6*f_3 - a_2*c_1*d_6*e_3*f_4 + a_2*c_1*d_6*e_4*f_3 + a_2*c_3*d_1*e_4*f_6 - a_2*c_3*d_1*e_6*f_4 - a_2*c_3*d_4*e_1*f_6 + a_2*c_3*d_4*e_6*f_1 + a_2*c_3*d_6*e_1*f_4 - a_2*c_3*d_6*e_4*f_1 - a_2*c_4*d_1*e_3*f_6 + a_2*c_4*d_1*e_6*f_3 + a_2*c_4*d_3*e_1*f_6 - a_2*c_4*d_3*e_6*f_1 - a_2*c_4*d_6*e_1*f_3 + a_2*c_4*d_6*e_3*f_1 + a_2*c_6*d_1*e_3*f_4 - a_2*c_6*d_1*e_4*f_3 - a_2*c_6*d_3*e_1*f_4 + a_2*c_6*d_3*e_4*f_1 + a_2*c_6*d_4*e_1*f_3 - a_2*c_6*d_4*e_3*f_1 + a_3*c_1*d_2*e_4*f_6 - a_3*c_1*d_2*e_6*f_4 - a_3*c_1*d_4*e_2*f_6 + a_3*c_1*d_4*e_6*f_2 + a_3*c_1*d_6*e_2*f_4 - a_3*c_1*d_6*e_4*f_2 - a_3*c_2*d_1*e_4*f_6 + a_3*c_2*d_1*e_6*f_4 + a_3*c_2*d_4*e_1*f_6 - a_3*c_2*d_4*e_6*f_1 - a_3*c_2*d_6*e_1*f_4 + a_3*c_2*d_6*e_4*f_1 + a_3*c_4*d_1*e_2*f_6 - a_3*c_4*d_1*e_6*f_2 - a_3*c_4*d_2*e_1*f_6 + a_3*c_4*d_2*e_6*f_1 + a_3*c_4*d_6*e_1*f_2 - a_3*c_4*d_6*e_2*f_1 - a_3*c_6*d_1*e_2*f_4 + a_3*c_6*d_1*e_4*f_2 + a_3*c_6*d_2*e_1*f_4 - a_3*c_6*d_2*e_4*f_1 - a_3*c_6*d_4*e_1*f_2 + a_3*c_6*d_4*e_2*f_1 - a_4*c_1*d_2*e_3*f_6 + a_4*c_1*d_2*e_6*f_3 + a_4*c_1*d_3*e_2*f_6 - a_4*c_1*d_3*e_6*f_2 - a_4*c_1*d_6*e_2*f_3 + a_4*c_1*d_6*e_3*f_2 + a_4*c_2*d_1*e_3*f_6 - a_4*c_2*d_1*e_6*f_3 - a_4*c_2*d_3*e_1*f_6 + a_4*c_2*d_3*e_6*f_1 + a_4*c_2*d_6*e_1*f_3 - a_4*c_2*d_6*e_3*f_1 - a_4*c_3*d_1*e_2*f_6 + a_4*c_3*d_1*e_6*f_2 + a_4*c_3*d_2*e_1*f_6 - a_4*c_3*d_2*e_6*f_1 - a_4*c_3*d_6*e_1*f_2 + a_4*c_3*d_6*e_2*f_1 + a_4*c_6*d_1*e_2*f_3 - a_4*c_6*d_1*e_3*f_2 - a_4*c_6*d_2*e_1*f_3 + a_4*c_6*d_2*e_3*f_1 + a_4*c_6*d_3*e_1*f_2 - a_4*c_6*d_3*e_2*f_1 + a_6*c_1*d_2*e_3*f_4 - a_6*c_1*d_2*e_4*f_3 - a_6*c_1*d_3*e_2*f_4 + a_6*c_1*d_3*e_4*f_2 + a_6*c_1*d_4*e_2*f_3 - a_6*c_1*d_4*e_3*f_2 - a_6*c_2*d_1*e_3*f_4 + a_6*c_2*d_1*e_4*f_3 + a_6*c_2*d_3*e_1*f_4 - a_6*c_2*d_3*e_4*f_1 - a_6*c_2*d_4*e_1*f_3 + a_6*c_2*d_4*e_3*f_1 + a_6*c_3*d_1*e_2*f_4 - a_6*c_3*d_1*e_4*f_2 - a_6*c_3*d_2*e_1*f_4 + a_6*c_3*d_2*e_4*f_1 + a_6*c_3*d_4*e_1*f_2 - a_6*c_3*d_4*e_2*f_1 - a_6*c_4*d_1*e_2*f_3 + a_6*c_4*d_1*e_3*f_2 + a_6*c_4*d_2*e_1*f_3 - a_6*c_4*d_2*e_3*f_1 - a_6*c_4*d_3*e_1*f_2 + a_6*c_4*d_3*e_2*f_1)/martix_value
              Inverse_A(3,5)=(a_1*b_2*d_3*e_4*f_6 - a_1*b_2*d_3*e_6*f_4 - a_1*b_2*d_4*e_3*f_6 + a_1*b_2*d_4*e_6*f_3 + a_1*b_2*d_6*e_3*f_4 - a_1*b_2*d_6*e_4*f_3 - a_1*b_3*d_2*e_4*f_6 + a_1*b_3*d_2*e_6*f_4 + a_1*b_3*d_4*e_2*f_6 - a_1*b_3*d_4*e_6*f_2 - a_1*b_3*d_6*e_2*f_4 + a_1*b_3*d_6*e_4*f_2 + a_1*b_4*d_2*e_3*f_6 - a_1*b_4*d_2*e_6*f_3 - a_1*b_4*d_3*e_2*f_6 + a_1*b_4*d_3*e_6*f_2 + a_1*b_4*d_6*e_2*f_3 - a_1*b_4*d_6*e_3*f_2 - a_1*b_6*d_2*e_3*f_4 + a_1*b_6*d_2*e_4*f_3 + a_1*b_6*d_3*e_2*f_4 - a_1*b_6*d_3*e_4*f_2 - a_1*b_6*d_4*e_2*f_3 + a_1*b_6*d_4*e_3*f_2 - a_2*b_1*d_3*e_4*f_6 + a_2*b_1*d_3*e_6*f_4 + a_2*b_1*d_4*e_3*f_6 - a_2*b_1*d_4*e_6*f_3 - a_2*b_1*d_6*e_3*f_4 + a_2*b_1*d_6*e_4*f_3 + a_2*b_3*d_1*e_4*f_6 - a_2*b_3*d_1*e_6*f_4 - a_2*b_3*d_4*e_1*f_6 + a_2*b_3*d_4*e_6*f_1 + a_2*b_3*d_6*e_1*f_4 - a_2*b_3*d_6*e_4*f_1 - a_2*b_4*d_1*e_3*f_6 + a_2*b_4*d_1*e_6*f_3 + a_2*b_4*d_3*e_1*f_6 - a_2*b_4*d_3*e_6*f_1 - a_2*b_4*d_6*e_1*f_3 + a_2*b_4*d_6*e_3*f_1 + a_2*b_6*d_1*e_3*f_4 - a_2*b_6*d_1*e_4*f_3 - a_2*b_6*d_3*e_1*f_4 + a_2*b_6*d_3*e_4*f_1 + a_2*b_6*d_4*e_1*f_3 - a_2*b_6*d_4*e_3*f_1 + a_3*b_1*d_2*e_4*f_6 - a_3*b_1*d_2*e_6*f_4 - a_3*b_1*d_4*e_2*f_6 + a_3*b_1*d_4*e_6*f_2 + a_3*b_1*d_6*e_2*f_4 - a_3*b_1*d_6*e_4*f_2 - a_3*b_2*d_1*e_4*f_6 + a_3*b_2*d_1*e_6*f_4 + a_3*b_2*d_4*e_1*f_6 - a_3*b_2*d_4*e_6*f_1 - a_3*b_2*d_6*e_1*f_4 + a_3*b_2*d_6*e_4*f_1 + a_3*b_4*d_1*e_2*f_6 - a_3*b_4*d_1*e_6*f_2 - a_3*b_4*d_2*e_1*f_6 + a_3*b_4*d_2*e_6*f_1 + a_3*b_4*d_6*e_1*f_2 - a_3*b_4*d_6*e_2*f_1 - a_3*b_6*d_1*e_2*f_4 + a_3*b_6*d_1*e_4*f_2 + a_3*b_6*d_2*e_1*f_4 - a_3*b_6*d_2*e_4*f_1 - a_3*b_6*d_4*e_1*f_2 + a_3*b_6*d_4*e_2*f_1 - a_4*b_1*d_2*e_3*f_6 + a_4*b_1*d_2*e_6*f_3 + a_4*b_1*d_3*e_2*f_6 - a_4*b_1*d_3*e_6*f_2 - a_4*b_1*d_6*e_2*f_3 + a_4*b_1*d_6*e_3*f_2 + a_4*b_2*d_1*e_3*f_6 - a_4*b_2*d_1*e_6*f_3 - a_4*b_2*d_3*e_1*f_6 + a_4*b_2*d_3*e_6*f_1 + a_4*b_2*d_6*e_1*f_3 - a_4*b_2*d_6*e_3*f_1 - a_4*b_3*d_1*e_2*f_6 + a_4*b_3*d_1*e_6*f_2 + a_4*b_3*d_2*e_1*f_6 - a_4*b_3*d_2*e_6*f_1 - a_4*b_3*d_6*e_1*f_2 + a_4*b_3*d_6*e_2*f_1 + a_4*b_6*d_1*e_2*f_3 - a_4*b_6*d_1*e_3*f_2 - a_4*b_6*d_2*e_1*f_3 + a_4*b_6*d_2*e_3*f_1 + a_4*b_6*d_3*e_1*f_2 - a_4*b_6*d_3*e_2*f_1 + a_6*b_1*d_2*e_3*f_4 - a_6*b_1*d_2*e_4*f_3 - a_6*b_1*d_3*e_2*f_4 + a_6*b_1*d_3*e_4*f_2 + a_6*b_1*d_4*e_2*f_3 - a_6*b_1*d_4*e_3*f_2 - a_6*b_2*d_1*e_3*f_4 + a_6*b_2*d_1*e_4*f_3 + a_6*b_2*d_3*e_1*f_4 - a_6*b_2*d_3*e_4*f_1 - a_6*b_2*d_4*e_1*f_3 + a_6*b_2*d_4*e_3*f_1 + a_6*b_3*d_1*e_2*f_4 - a_6*b_3*d_1*e_4*f_2 - a_6*b_3*d_2*e_1*f_4 + a_6*b_3*d_2*e_4*f_1 + a_6*b_3*d_4*e_1*f_2 - a_6*b_3*d_4*e_2*f_1 - a_6*b_4*d_1*e_2*f_3 + a_6*b_4*d_1*e_3*f_2 + a_6*b_4*d_2*e_1*f_3 - a_6*b_4*d_2*e_3*f_1 - a_6*b_4*d_3*e_1*f_2 + a_6*b_4*d_3*e_2*f_1)/martix_value
              Inverse_A(4,5)=-(a_1*b_2*c_3*e_4*f_6 - a_1*b_2*c_3*e_6*f_4 - a_1*b_2*c_4*e_3*f_6 + a_1*b_2*c_4*e_6*f_3 + a_1*b_2*c_6*e_3*f_4 - a_1*b_2*c_6*e_4*f_3 - a_1*b_3*c_2*e_4*f_6 + a_1*b_3*c_2*e_6*f_4 + a_1*b_3*c_4*e_2*f_6 - a_1*b_3*c_4*e_6*f_2 - a_1*b_3*c_6*e_2*f_4 + a_1*b_3*c_6*e_4*f_2 + a_1*b_4*c_2*e_3*f_6 - a_1*b_4*c_2*e_6*f_3 - a_1*b_4*c_3*e_2*f_6 + a_1*b_4*c_3*e_6*f_2 + a_1*b_4*c_6*e_2*f_3 - a_1*b_4*c_6*e_3*f_2 - a_1*b_6*c_2*e_3*f_4 + a_1*b_6*c_2*e_4*f_3 + a_1*b_6*c_3*e_2*f_4 - a_1*b_6*c_3*e_4*f_2 - a_1*b_6*c_4*e_2*f_3 + a_1*b_6*c_4*e_3*f_2 - a_2*b_1*c_3*e_4*f_6 + a_2*b_1*c_3*e_6*f_4 + a_2*b_1*c_4*e_3*f_6 - a_2*b_1*c_4*e_6*f_3 - a_2*b_1*c_6*e_3*f_4 + a_2*b_1*c_6*e_4*f_3 + a_2*b_3*c_1*e_4*f_6 - a_2*b_3*c_1*e_6*f_4 - a_2*b_3*c_4*e_1*f_6 + a_2*b_3*c_4*e_6*f_1 + a_2*b_3*c_6*e_1*f_4 - a_2*b_3*c_6*e_4*f_1 - a_2*b_4*c_1*e_3*f_6 + a_2*b_4*c_1*e_6*f_3 + a_2*b_4*c_3*e_1*f_6 - a_2*b_4*c_3*e_6*f_1 - a_2*b_4*c_6*e_1*f_3 + a_2*b_4*c_6*e_3*f_1 + a_2*b_6*c_1*e_3*f_4 - a_2*b_6*c_1*e_4*f_3 - a_2*b_6*c_3*e_1*f_4 + a_2*b_6*c_3*e_4*f_1 + a_2*b_6*c_4*e_1*f_3 - a_2*b_6*c_4*e_3*f_1 + a_3*b_1*c_2*e_4*f_6 - a_3*b_1*c_2*e_6*f_4 - a_3*b_1*c_4*e_2*f_6 + a_3*b_1*c_4*e_6*f_2 + a_3*b_1*c_6*e_2*f_4 - a_3*b_1*c_6*e_4*f_2 - a_3*b_2*c_1*e_4*f_6 + a_3*b_2*c_1*e_6*f_4 + a_3*b_2*c_4*e_1*f_6 - a_3*b_2*c_4*e_6*f_1 - a_3*b_2*c_6*e_1*f_4 + a_3*b_2*c_6*e_4*f_1 + a_3*b_4*c_1*e_2*f_6 - a_3*b_4*c_1*e_6*f_2 - a_3*b_4*c_2*e_1*f_6 + a_3*b_4*c_2*e_6*f_1 + a_3*b_4*c_6*e_1*f_2 - a_3*b_4*c_6*e_2*f_1 - a_3*b_6*c_1*e_2*f_4 + a_3*b_6*c_1*e_4*f_2 + a_3*b_6*c_2*e_1*f_4 - a_3*b_6*c_2*e_4*f_1 - a_3*b_6*c_4*e_1*f_2 + a_3*b_6*c_4*e_2*f_1 - a_4*b_1*c_2*e_3*f_6 + a_4*b_1*c_2*e_6*f_3 + a_4*b_1*c_3*e_2*f_6 - a_4*b_1*c_3*e_6*f_2 - a_4*b_1*c_6*e_2*f_3 + a_4*b_1*c_6*e_3*f_2 + a_4*b_2*c_1*e_3*f_6 - a_4*b_2*c_1*e_6*f_3 - a_4*b_2*c_3*e_1*f_6 + a_4*b_2*c_3*e_6*f_1 + a_4*b_2*c_6*e_1*f_3 - a_4*b_2*c_6*e_3*f_1 - a_4*b_3*c_1*e_2*f_6 + a_4*b_3*c_1*e_6*f_2 + a_4*b_3*c_2*e_1*f_6 - a_4*b_3*c_2*e_6*f_1 - a_4*b_3*c_6*e_1*f_2 + a_4*b_3*c_6*e_2*f_1 + a_4*b_6*c_1*e_2*f_3 - a_4*b_6*c_1*e_3*f_2 - a_4*b_6*c_2*e_1*f_3 + a_4*b_6*c_2*e_3*f_1 + a_4*b_6*c_3*e_1*f_2 - a_4*b_6*c_3*e_2*f_1 + a_6*b_1*c_2*e_3*f_4 - a_6*b_1*c_2*e_4*f_3 - a_6*b_1*c_3*e_2*f_4 + a_6*b_1*c_3*e_4*f_2 + a_6*b_1*c_4*e_2*f_3 - a_6*b_1*c_4*e_3*f_2 - a_6*b_2*c_1*e_3*f_4 + a_6*b_2*c_1*e_4*f_3 + a_6*b_2*c_3*e_1*f_4 - a_6*b_2*c_3*e_4*f_1 - a_6*b_2*c_4*e_1*f_3 + a_6*b_2*c_4*e_3*f_1 + a_6*b_3*c_1*e_2*f_4 - a_6*b_3*c_1*e_4*f_2 - a_6*b_3*c_2*e_1*f_4 + a_6*b_3*c_2*e_4*f_1 + a_6*b_3*c_4*e_1*f_2 - a_6*b_3*c_4*e_2*f_1 - a_6*b_4*c_1*e_2*f_3 + a_6*b_4*c_1*e_3*f_2 + a_6*b_4*c_2*e_1*f_3 - a_6*b_4*c_2*e_3*f_1 - a_6*b_4*c_3*e_1*f_2 + a_6*b_4*c_3*e_2*f_1)/martix_value
              Inverse_A(5,5)=(a_1*b_2*c_3*d_4*f_6 - a_1*b_2*c_3*d_6*f_4 - a_1*b_2*c_4*d_3*f_6 + a_1*b_2*c_4*d_6*f_3 + a_1*b_2*c_6*d_3*f_4 - a_1*b_2*c_6*d_4*f_3 - a_1*b_3*c_2*d_4*f_6 + a_1*b_3*c_2*d_6*f_4 + a_1*b_3*c_4*d_2*f_6 - a_1*b_3*c_4*d_6*f_2 - a_1*b_3*c_6*d_2*f_4 + a_1*b_3*c_6*d_4*f_2 + a_1*b_4*c_2*d_3*f_6 - a_1*b_4*c_2*d_6*f_3 - a_1*b_4*c_3*d_2*f_6 + a_1*b_4*c_3*d_6*f_2 + a_1*b_4*c_6*d_2*f_3 - a_1*b_4*c_6*d_3*f_2 - a_1*b_6*c_2*d_3*f_4 + a_1*b_6*c_2*d_4*f_3 + a_1*b_6*c_3*d_2*f_4 - a_1*b_6*c_3*d_4*f_2 - a_1*b_6*c_4*d_2*f_3 + a_1*b_6*c_4*d_3*f_2 - a_2*b_1*c_3*d_4*f_6 + a_2*b_1*c_3*d_6*f_4 + a_2*b_1*c_4*d_3*f_6 - a_2*b_1*c_4*d_6*f_3 - a_2*b_1*c_6*d_3*f_4 + a_2*b_1*c_6*d_4*f_3 + a_2*b_3*c_1*d_4*f_6 - a_2*b_3*c_1*d_6*f_4 - a_2*b_3*c_4*d_1*f_6 + a_2*b_3*c_4*d_6*f_1 + a_2*b_3*c_6*d_1*f_4 - a_2*b_3*c_6*d_4*f_1 - a_2*b_4*c_1*d_3*f_6 + a_2*b_4*c_1*d_6*f_3 + a_2*b_4*c_3*d_1*f_6 - a_2*b_4*c_3*d_6*f_1 - a_2*b_4*c_6*d_1*f_3 + a_2*b_4*c_6*d_3*f_1 + a_2*b_6*c_1*d_3*f_4 - a_2*b_6*c_1*d_4*f_3 - a_2*b_6*c_3*d_1*f_4 + a_2*b_6*c_3*d_4*f_1 + a_2*b_6*c_4*d_1*f_3 - a_2*b_6*c_4*d_3*f_1 + a_3*b_1*c_2*d_4*f_6 - a_3*b_1*c_2*d_6*f_4 - a_3*b_1*c_4*d_2*f_6 + a_3*b_1*c_4*d_6*f_2 + a_3*b_1*c_6*d_2*f_4 - a_3*b_1*c_6*d_4*f_2 - a_3*b_2*c_1*d_4*f_6 + a_3*b_2*c_1*d_6*f_4 + a_3*b_2*c_4*d_1*f_6 - a_3*b_2*c_4*d_6*f_1 - a_3*b_2*c_6*d_1*f_4 + a_3*b_2*c_6*d_4*f_1 + a_3*b_4*c_1*d_2*f_6 - a_3*b_4*c_1*d_6*f_2 - a_3*b_4*c_2*d_1*f_6 + a_3*b_4*c_2*d_6*f_1 + a_3*b_4*c_6*d_1*f_2 - a_3*b_4*c_6*d_2*f_1 - a_3*b_6*c_1*d_2*f_4 + a_3*b_6*c_1*d_4*f_2 + a_3*b_6*c_2*d_1*f_4 - a_3*b_6*c_2*d_4*f_1 - a_3*b_6*c_4*d_1*f_2 + a_3*b_6*c_4*d_2*f_1 - a_4*b_1*c_2*d_3*f_6 + a_4*b_1*c_2*d_6*f_3 + a_4*b_1*c_3*d_2*f_6 - a_4*b_1*c_3*d_6*f_2 - a_4*b_1*c_6*d_2*f_3 + a_4*b_1*c_6*d_3*f_2 + a_4*b_2*c_1*d_3*f_6 - a_4*b_2*c_1*d_6*f_3 - a_4*b_2*c_3*d_1*f_6 + a_4*b_2*c_3*d_6*f_1 + a_4*b_2*c_6*d_1*f_3 - a_4*b_2*c_6*d_3*f_1 - a_4*b_3*c_1*d_2*f_6 + a_4*b_3*c_1*d_6*f_2 + a_4*b_3*c_2*d_1*f_6 - a_4*b_3*c_2*d_6*f_1 - a_4*b_3*c_6*d_1*f_2 + a_4*b_3*c_6*d_2*f_1 + a_4*b_6*c_1*d_2*f_3 - a_4*b_6*c_1*d_3*f_2 - a_4*b_6*c_2*d_1*f_3 + a_4*b_6*c_2*d_3*f_1 + a_4*b_6*c_3*d_1*f_2 - a_4*b_6*c_3*d_2*f_1 + a_6*b_1*c_2*d_3*f_4 - a_6*b_1*c_2*d_4*f_3 - a_6*b_1*c_3*d_2*f_4 + a_6*b_1*c_3*d_4*f_2 + a_6*b_1*c_4*d_2*f_3 - a_6*b_1*c_4*d_3*f_2 - a_6*b_2*c_1*d_3*f_4 + a_6*b_2*c_1*d_4*f_3 + a_6*b_2*c_3*d_1*f_4 - a_6*b_2*c_3*d_4*f_1 - a_6*b_2*c_4*d_1*f_3 + a_6*b_2*c_4*d_3*f_1 + a_6*b_3*c_1*d_2*f_4 - a_6*b_3*c_1*d_4*f_2 - a_6*b_3*c_2*d_1*f_4 + a_6*b_3*c_2*d_4*f_1 + a_6*b_3*c_4*d_1*f_2 - a_6*b_3*c_4*d_2*f_1 - a_6*b_4*c_1*d_2*f_3 + a_6*b_4*c_1*d_3*f_2 + a_6*b_4*c_2*d_1*f_3 - a_6*b_4*c_2*d_3*f_1 - a_6*b_4*c_3*d_1*f_2 + a_6*b_4*c_3*d_2*f_1)/martix_value
              Inverse_A(6,5)=-(a_1*b_2*c_3*d_4*e_6 - a_1*b_2*c_3*d_6*e_4 - a_1*b_2*c_4*d_3*e_6 + a_1*b_2*c_4*d_6*e_3 + a_1*b_2*c_6*d_3*e_4 - a_1*b_2*c_6*d_4*e_3 - a_1*b_3*c_2*d_4*e_6 + a_1*b_3*c_2*d_6*e_4 + a_1*b_3*c_4*d_2*e_6 - a_1*b_3*c_4*d_6*e_2 - a_1*b_3*c_6*d_2*e_4 + a_1*b_3*c_6*d_4*e_2 + a_1*b_4*c_2*d_3*e_6 - a_1*b_4*c_2*d_6*e_3 - a_1*b_4*c_3*d_2*e_6 + a_1*b_4*c_3*d_6*e_2 + a_1*b_4*c_6*d_2*e_3 - a_1*b_4*c_6*d_3*e_2 - a_1*b_6*c_2*d_3*e_4 + a_1*b_6*c_2*d_4*e_3 + a_1*b_6*c_3*d_2*e_4 - a_1*b_6*c_3*d_4*e_2 - a_1*b_6*c_4*d_2*e_3 + a_1*b_6*c_4*d_3*e_2 - a_2*b_1*c_3*d_4*e_6 + a_2*b_1*c_3*d_6*e_4 + a_2*b_1*c_4*d_3*e_6 - a_2*b_1*c_4*d_6*e_3 - a_2*b_1*c_6*d_3*e_4 + a_2*b_1*c_6*d_4*e_3 + a_2*b_3*c_1*d_4*e_6 - a_2*b_3*c_1*d_6*e_4 - a_2*b_3*c_4*d_1*e_6 + a_2*b_3*c_4*d_6*e_1 + a_2*b_3*c_6*d_1*e_4 - a_2*b_3*c_6*d_4*e_1 - a_2*b_4*c_1*d_3*e_6 + a_2*b_4*c_1*d_6*e_3 + a_2*b_4*c_3*d_1*e_6 - a_2*b_4*c_3*d_6*e_1 - a_2*b_4*c_6*d_1*e_3 + a_2*b_4*c_6*d_3*e_1 + a_2*b_6*c_1*d_3*e_4 - a_2*b_6*c_1*d_4*e_3 - a_2*b_6*c_3*d_1*e_4 + a_2*b_6*c_3*d_4*e_1 + a_2*b_6*c_4*d_1*e_3 - a_2*b_6*c_4*d_3*e_1 + a_3*b_1*c_2*d_4*e_6 - a_3*b_1*c_2*d_6*e_4 - a_3*b_1*c_4*d_2*e_6 + a_3*b_1*c_4*d_6*e_2 + a_3*b_1*c_6*d_2*e_4 - a_3*b_1*c_6*d_4*e_2 - a_3*b_2*c_1*d_4*e_6 + a_3*b_2*c_1*d_6*e_4 + a_3*b_2*c_4*d_1*e_6 - a_3*b_2*c_4*d_6*e_1 - a_3*b_2*c_6*d_1*e_4 + a_3*b_2*c_6*d_4*e_1 + a_3*b_4*c_1*d_2*e_6 - a_3*b_4*c_1*d_6*e_2 - a_3*b_4*c_2*d_1*e_6 + a_3*b_4*c_2*d_6*e_1 + a_3*b_4*c_6*d_1*e_2 - a_3*b_4*c_6*d_2*e_1 - a_3*b_6*c_1*d_2*e_4 + a_3*b_6*c_1*d_4*e_2 + a_3*b_6*c_2*d_1*e_4 - a_3*b_6*c_2*d_4*e_1 - a_3*b_6*c_4*d_1*e_2 + a_3*b_6*c_4*d_2*e_1 - a_4*b_1*c_2*d_3*e_6 + a_4*b_1*c_2*d_6*e_3 + a_4*b_1*c_3*d_2*e_6 - a_4*b_1*c_3*d_6*e_2 - a_4*b_1*c_6*d_2*e_3 + a_4*b_1*c_6*d_3*e_2 + a_4*b_2*c_1*d_3*e_6 - a_4*b_2*c_1*d_6*e_3 - a_4*b_2*c_3*d_1*e_6 + a_4*b_2*c_3*d_6*e_1 + a_4*b_2*c_6*d_1*e_3 - a_4*b_2*c_6*d_3*e_1 - a_4*b_3*c_1*d_2*e_6 + a_4*b_3*c_1*d_6*e_2 + a_4*b_3*c_2*d_1*e_6 - a_4*b_3*c_2*d_6*e_1 - a_4*b_3*c_6*d_1*e_2 + a_4*b_3*c_6*d_2*e_1 + a_4*b_6*c_1*d_2*e_3 - a_4*b_6*c_1*d_3*e_2 - a_4*b_6*c_2*d_1*e_3 + a_4*b_6*c_2*d_3*e_1 + a_4*b_6*c_3*d_1*e_2 - a_4*b_6*c_3*d_2*e_1 + a_6*b_1*c_2*d_3*e_4 - a_6*b_1*c_2*d_4*e_3 - a_6*b_1*c_3*d_2*e_4 + a_6*b_1*c_3*d_4*e_2 + a_6*b_1*c_4*d_2*e_3 - a_6*b_1*c_4*d_3*e_2 - a_6*b_2*c_1*d_3*e_4 + a_6*b_2*c_1*d_4*e_3 + a_6*b_2*c_3*d_1*e_4 - a_6*b_2*c_3*d_4*e_1 - a_6*b_2*c_4*d_1*e_3 + a_6*b_2*c_4*d_3*e_1 + a_6*b_3*c_1*d_2*e_4 - a_6*b_3*c_1*d_4*e_2 - a_6*b_3*c_2*d_1*e_4 + a_6*b_3*c_2*d_4*e_1 + a_6*b_3*c_4*d_1*e_2 - a_6*b_3*c_4*d_2*e_1 - a_6*b_4*c_1*d_2*e_3 + a_6*b_4*c_1*d_3*e_2 + a_6*b_4*c_2*d_1*e_3 - a_6*b_4*c_2*d_3*e_1 - a_6*b_4*c_3*d_1*e_2 + a_6*b_4*c_3*d_2*e_1)/martix_value
      
              Inverse_A(1,6)=-(b_1*c_2*d_3*e_4*f_5 - b_1*c_2*d_3*e_5*f_4 - b_1*c_2*d_4*e_3*f_5 + b_1*c_2*d_4*e_5*f_3 + b_1*c_2*d_5*e_3*f_4 - b_1*c_2*d_5*e_4*f_3 - b_1*c_3*d_2*e_4*f_5 + b_1*c_3*d_2*e_5*f_4 + b_1*c_3*d_4*e_2*f_5 - b_1*c_3*d_4*e_5*f_2 - b_1*c_3*d_5*e_2*f_4 + b_1*c_3*d_5*e_4*f_2 + b_1*c_4*d_2*e_3*f_5 - b_1*c_4*d_2*e_5*f_3 - b_1*c_4*d_3*e_2*f_5 + b_1*c_4*d_3*e_5*f_2 + b_1*c_4*d_5*e_2*f_3 - b_1*c_4*d_5*e_3*f_2 - b_1*c_5*d_2*e_3*f_4 + b_1*c_5*d_2*e_4*f_3 + b_1*c_5*d_3*e_2*f_4 - b_1*c_5*d_3*e_4*f_2 - b_1*c_5*d_4*e_2*f_3 + b_1*c_5*d_4*e_3*f_2 - b_2*c_1*d_3*e_4*f_5 + b_2*c_1*d_3*e_5*f_4 + b_2*c_1*d_4*e_3*f_5 - b_2*c_1*d_4*e_5*f_3 - b_2*c_1*d_5*e_3*f_4 + b_2*c_1*d_5*e_4*f_3 + b_2*c_3*d_1*e_4*f_5 - b_2*c_3*d_1*e_5*f_4 - b_2*c_3*d_4*e_1*f_5 + b_2*c_3*d_4*e_5*f_1 + b_2*c_3*d_5*e_1*f_4 - b_2*c_3*d_5*e_4*f_1 - b_2*c_4*d_1*e_3*f_5 + b_2*c_4*d_1*e_5*f_3 + b_2*c_4*d_3*e_1*f_5 - b_2*c_4*d_3*e_5*f_1 - b_2*c_4*d_5*e_1*f_3 + b_2*c_4*d_5*e_3*f_1 + b_2*c_5*d_1*e_3*f_4 - b_2*c_5*d_1*e_4*f_3 - b_2*c_5*d_3*e_1*f_4 + b_2*c_5*d_3*e_4*f_1 + b_2*c_5*d_4*e_1*f_3 - b_2*c_5*d_4*e_3*f_1 + b_3*c_1*d_2*e_4*f_5 - b_3*c_1*d_2*e_5*f_4 - b_3*c_1*d_4*e_2*f_5 + b_3*c_1*d_4*e_5*f_2 + b_3*c_1*d_5*e_2*f_4 - b_3*c_1*d_5*e_4*f_2 - b_3*c_2*d_1*e_4*f_5 + b_3*c_2*d_1*e_5*f_4 + b_3*c_2*d_4*e_1*f_5 - b_3*c_2*d_4*e_5*f_1 - b_3*c_2*d_5*e_1*f_4 + b_3*c_2*d_5*e_4*f_1 + b_3*c_4*d_1*e_2*f_5 - b_3*c_4*d_1*e_5*f_2 - b_3*c_4*d_2*e_1*f_5 + b_3*c_4*d_2*e_5*f_1 + b_3*c_4*d_5*e_1*f_2 - b_3*c_4*d_5*e_2*f_1 - b_3*c_5*d_1*e_2*f_4 + b_3*c_5*d_1*e_4*f_2 + b_3*c_5*d_2*e_1*f_4 - b_3*c_5*d_2*e_4*f_1 - b_3*c_5*d_4*e_1*f_2 + b_3*c_5*d_4*e_2*f_1 - b_4*c_1*d_2*e_3*f_5 + b_4*c_1*d_2*e_5*f_3 + b_4*c_1*d_3*e_2*f_5 - b_4*c_1*d_3*e_5*f_2 - b_4*c_1*d_5*e_2*f_3 + b_4*c_1*d_5*e_3*f_2 + b_4*c_2*d_1*e_3*f_5 - b_4*c_2*d_1*e_5*f_3 - b_4*c_2*d_3*e_1*f_5 + b_4*c_2*d_3*e_5*f_1 + b_4*c_2*d_5*e_1*f_3 - b_4*c_2*d_5*e_3*f_1 - b_4*c_3*d_1*e_2*f_5 + b_4*c_3*d_1*e_5*f_2 + b_4*c_3*d_2*e_1*f_5 - b_4*c_3*d_2*e_5*f_1 - b_4*c_3*d_5*e_1*f_2 + b_4*c_3*d_5*e_2*f_1 + b_4*c_5*d_1*e_2*f_3 - b_4*c_5*d_1*e_3*f_2 - b_4*c_5*d_2*e_1*f_3 + b_4*c_5*d_2*e_3*f_1 + b_4*c_5*d_3*e_1*f_2 - b_4*c_5*d_3*e_2*f_1 + b_5*c_1*d_2*e_3*f_4 - b_5*c_1*d_2*e_4*f_3 - b_5*c_1*d_3*e_2*f_4 + b_5*c_1*d_3*e_4*f_2 + b_5*c_1*d_4*e_2*f_3 - b_5*c_1*d_4*e_3*f_2 - b_5*c_2*d_1*e_3*f_4 + b_5*c_2*d_1*e_4*f_3 + b_5*c_2*d_3*e_1*f_4 - b_5*c_2*d_3*e_4*f_1 - b_5*c_2*d_4*e_1*f_3 + b_5*c_2*d_4*e_3*f_1 + b_5*c_3*d_1*e_2*f_4 - b_5*c_3*d_1*e_4*f_2 - b_5*c_3*d_2*e_1*f_4 + b_5*c_3*d_2*e_4*f_1 + b_5*c_3*d_4*e_1*f_2 - b_5*c_3*d_4*e_2*f_1 - b_5*c_4*d_1*e_2*f_3 + b_5*c_4*d_1*e_3*f_2 + b_5*c_4*d_2*e_1*f_3 - b_5*c_4*d_2*e_3*f_1 - b_5*c_4*d_3*e_1*f_2 + b_5*c_4*d_3*e_2*f_1)/martix_value
              Inverse_A(2,6)=(a_1*c_2*d_3*e_4*f_5 - a_1*c_2*d_3*e_5*f_4 - a_1*c_2*d_4*e_3*f_5 + a_1*c_2*d_4*e_5*f_3 + a_1*c_2*d_5*e_3*f_4 - a_1*c_2*d_5*e_4*f_3 - a_1*c_3*d_2*e_4*f_5 + a_1*c_3*d_2*e_5*f_4 + a_1*c_3*d_4*e_2*f_5 - a_1*c_3*d_4*e_5*f_2 - a_1*c_3*d_5*e_2*f_4 + a_1*c_3*d_5*e_4*f_2 + a_1*c_4*d_2*e_3*f_5 - a_1*c_4*d_2*e_5*f_3 - a_1*c_4*d_3*e_2*f_5 + a_1*c_4*d_3*e_5*f_2 + a_1*c_4*d_5*e_2*f_3 - a_1*c_4*d_5*e_3*f_2 - a_1*c_5*d_2*e_3*f_4 + a_1*c_5*d_2*e_4*f_3 + a_1*c_5*d_3*e_2*f_4 - a_1*c_5*d_3*e_4*f_2 - a_1*c_5*d_4*e_2*f_3 + a_1*c_5*d_4*e_3*f_2 - a_2*c_1*d_3*e_4*f_5 + a_2*c_1*d_3*e_5*f_4 + a_2*c_1*d_4*e_3*f_5 - a_2*c_1*d_4*e_5*f_3 - a_2*c_1*d_5*e_3*f_4 + a_2*c_1*d_5*e_4*f_3 + a_2*c_3*d_1*e_4*f_5 - a_2*c_3*d_1*e_5*f_4 - a_2*c_3*d_4*e_1*f_5 + a_2*c_3*d_4*e_5*f_1 + a_2*c_3*d_5*e_1*f_4 - a_2*c_3*d_5*e_4*f_1 - a_2*c_4*d_1*e_3*f_5 + a_2*c_4*d_1*e_5*f_3 + a_2*c_4*d_3*e_1*f_5 - a_2*c_4*d_3*e_5*f_1 - a_2*c_4*d_5*e_1*f_3 + a_2*c_4*d_5*e_3*f_1 + a_2*c_5*d_1*e_3*f_4 - a_2*c_5*d_1*e_4*f_3 - a_2*c_5*d_3*e_1*f_4 + a_2*c_5*d_3*e_4*f_1 + a_2*c_5*d_4*e_1*f_3 - a_2*c_5*d_4*e_3*f_1 + a_3*c_1*d_2*e_4*f_5 - a_3*c_1*d_2*e_5*f_4 - a_3*c_1*d_4*e_2*f_5 + a_3*c_1*d_4*e_5*f_2 + a_3*c_1*d_5*e_2*f_4 - a_3*c_1*d_5*e_4*f_2 - a_3*c_2*d_1*e_4*f_5 + a_3*c_2*d_1*e_5*f_4 + a_3*c_2*d_4*e_1*f_5 - a_3*c_2*d_4*e_5*f_1 - a_3*c_2*d_5*e_1*f_4 + a_3*c_2*d_5*e_4*f_1 + a_3*c_4*d_1*e_2*f_5 - a_3*c_4*d_1*e_5*f_2 - a_3*c_4*d_2*e_1*f_5 + a_3*c_4*d_2*e_5*f_1 + a_3*c_4*d_5*e_1*f_2 - a_3*c_4*d_5*e_2*f_1 - a_3*c_5*d_1*e_2*f_4 + a_3*c_5*d_1*e_4*f_2 + a_3*c_5*d_2*e_1*f_4 - a_3*c_5*d_2*e_4*f_1 - a_3*c_5*d_4*e_1*f_2 + a_3*c_5*d_4*e_2*f_1 - a_4*c_1*d_2*e_3*f_5 + a_4*c_1*d_2*e_5*f_3 + a_4*c_1*d_3*e_2*f_5 - a_4*c_1*d_3*e_5*f_2 - a_4*c_1*d_5*e_2*f_3 + a_4*c_1*d_5*e_3*f_2 + a_4*c_2*d_1*e_3*f_5 - a_4*c_2*d_1*e_5*f_3 - a_4*c_2*d_3*e_1*f_5 + a_4*c_2*d_3*e_5*f_1 + a_4*c_2*d_5*e_1*f_3 - a_4*c_2*d_5*e_3*f_1 - a_4*c_3*d_1*e_2*f_5 + a_4*c_3*d_1*e_5*f_2 + a_4*c_3*d_2*e_1*f_5 - a_4*c_3*d_2*e_5*f_1 - a_4*c_3*d_5*e_1*f_2 + a_4*c_3*d_5*e_2*f_1 + a_4*c_5*d_1*e_2*f_3 - a_4*c_5*d_1*e_3*f_2 - a_4*c_5*d_2*e_1*f_3 + a_4*c_5*d_2*e_3*f_1 + a_4*c_5*d_3*e_1*f_2 - a_4*c_5*d_3*e_2*f_1 + a_5*c_1*d_2*e_3*f_4 - a_5*c_1*d_2*e_4*f_3 - a_5*c_1*d_3*e_2*f_4 + a_5*c_1*d_3*e_4*f_2 + a_5*c_1*d_4*e_2*f_3 - a_5*c_1*d_4*e_3*f_2 - a_5*c_2*d_1*e_3*f_4 + a_5*c_2*d_1*e_4*f_3 + a_5*c_2*d_3*e_1*f_4 - a_5*c_2*d_3*e_4*f_1 - a_5*c_2*d_4*e_1*f_3 + a_5*c_2*d_4*e_3*f_1 + a_5*c_3*d_1*e_2*f_4 - a_5*c_3*d_1*e_4*f_2 - a_5*c_3*d_2*e_1*f_4 + a_5*c_3*d_2*e_4*f_1 + a_5*c_3*d_4*e_1*f_2 - a_5*c_3*d_4*e_2*f_1 - a_5*c_4*d_1*e_2*f_3 + a_5*c_4*d_1*e_3*f_2 + a_5*c_4*d_2*e_1*f_3 - a_5*c_4*d_2*e_3*f_1 - a_5*c_4*d_3*e_1*f_2 + a_5*c_4*d_3*e_2*f_1)/martix_value
              Inverse_A(3,6)=-(a_1*b_2*d_3*e_4*f_5 - a_1*b_2*d_3*e_5*f_4 - a_1*b_2*d_4*e_3*f_5 + a_1*b_2*d_4*e_5*f_3 + a_1*b_2*d_5*e_3*f_4 - a_1*b_2*d_5*e_4*f_3 - a_1*b_3*d_2*e_4*f_5 + a_1*b_3*d_2*e_5*f_4 + a_1*b_3*d_4*e_2*f_5 - a_1*b_3*d_4*e_5*f_2 - a_1*b_3*d_5*e_2*f_4 + a_1*b_3*d_5*e_4*f_2 + a_1*b_4*d_2*e_3*f_5 - a_1*b_4*d_2*e_5*f_3 - a_1*b_4*d_3*e_2*f_5 + a_1*b_4*d_3*e_5*f_2 + a_1*b_4*d_5*e_2*f_3 - a_1*b_4*d_5*e_3*f_2 - a_1*b_5*d_2*e_3*f_4 + a_1*b_5*d_2*e_4*f_3 + a_1*b_5*d_3*e_2*f_4 - a_1*b_5*d_3*e_4*f_2 - a_1*b_5*d_4*e_2*f_3 + a_1*b_5*d_4*e_3*f_2 - a_2*b_1*d_3*e_4*f_5 + a_2*b_1*d_3*e_5*f_4 + a_2*b_1*d_4*e_3*f_5 - a_2*b_1*d_4*e_5*f_3 - a_2*b_1*d_5*e_3*f_4 + a_2*b_1*d_5*e_4*f_3 + a_2*b_3*d_1*e_4*f_5 - a_2*b_3*d_1*e_5*f_4 - a_2*b_3*d_4*e_1*f_5 + a_2*b_3*d_4*e_5*f_1 + a_2*b_3*d_5*e_1*f_4 - a_2*b_3*d_5*e_4*f_1 - a_2*b_4*d_1*e_3*f_5 + a_2*b_4*d_1*e_5*f_3 + a_2*b_4*d_3*e_1*f_5 - a_2*b_4*d_3*e_5*f_1 - a_2*b_4*d_5*e_1*f_3 + a_2*b_4*d_5*e_3*f_1 + a_2*b_5*d_1*e_3*f_4 - a_2*b_5*d_1*e_4*f_3 - a_2*b_5*d_3*e_1*f_4 + a_2*b_5*d_3*e_4*f_1 + a_2*b_5*d_4*e_1*f_3 - a_2*b_5*d_4*e_3*f_1 + a_3*b_1*d_2*e_4*f_5 - a_3*b_1*d_2*e_5*f_4 - a_3*b_1*d_4*e_2*f_5 + a_3*b_1*d_4*e_5*f_2 + a_3*b_1*d_5*e_2*f_4 - a_3*b_1*d_5*e_4*f_2 - a_3*b_2*d_1*e_4*f_5 + a_3*b_2*d_1*e_5*f_4 + a_3*b_2*d_4*e_1*f_5 - a_3*b_2*d_4*e_5*f_1 - a_3*b_2*d_5*e_1*f_4 + a_3*b_2*d_5*e_4*f_1 + a_3*b_4*d_1*e_2*f_5 - a_3*b_4*d_1*e_5*f_2 - a_3*b_4*d_2*e_1*f_5 + a_3*b_4*d_2*e_5*f_1 + a_3*b_4*d_5*e_1*f_2 - a_3*b_4*d_5*e_2*f_1 - a_3*b_5*d_1*e_2*f_4 + a_3*b_5*d_1*e_4*f_2 + a_3*b_5*d_2*e_1*f_4 - a_3*b_5*d_2*e_4*f_1 - a_3*b_5*d_4*e_1*f_2 + a_3*b_5*d_4*e_2*f_1 - a_4*b_1*d_2*e_3*f_5 + a_4*b_1*d_2*e_5*f_3 + a_4*b_1*d_3*e_2*f_5 - a_4*b_1*d_3*e_5*f_2 - a_4*b_1*d_5*e_2*f_3 + a_4*b_1*d_5*e_3*f_2 + a_4*b_2*d_1*e_3*f_5 - a_4*b_2*d_1*e_5*f_3 - a_4*b_2*d_3*e_1*f_5 + a_4*b_2*d_3*e_5*f_1 + a_4*b_2*d_5*e_1*f_3 - a_4*b_2*d_5*e_3*f_1 - a_4*b_3*d_1*e_2*f_5 + a_4*b_3*d_1*e_5*f_2 + a_4*b_3*d_2*e_1*f_5 - a_4*b_3*d_2*e_5*f_1 - a_4*b_3*d_5*e_1*f_2 + a_4*b_3*d_5*e_2*f_1 + a_4*b_5*d_1*e_2*f_3 - a_4*b_5*d_1*e_3*f_2 - a_4*b_5*d_2*e_1*f_3 + a_4*b_5*d_2*e_3*f_1 + a_4*b_5*d_3*e_1*f_2 - a_4*b_5*d_3*e_2*f_1 + a_5*b_1*d_2*e_3*f_4 - a_5*b_1*d_2*e_4*f_3 - a_5*b_1*d_3*e_2*f_4 + a_5*b_1*d_3*e_4*f_2 + a_5*b_1*d_4*e_2*f_3 - a_5*b_1*d_4*e_3*f_2 - a_5*b_2*d_1*e_3*f_4 + a_5*b_2*d_1*e_4*f_3 + a_5*b_2*d_3*e_1*f_4 - a_5*b_2*d_3*e_4*f_1 - a_5*b_2*d_4*e_1*f_3 + a_5*b_2*d_4*e_3*f_1 + a_5*b_3*d_1*e_2*f_4 - a_5*b_3*d_1*e_4*f_2 - a_5*b_3*d_2*e_1*f_4 + a_5*b_3*d_2*e_4*f_1 + a_5*b_3*d_4*e_1*f_2 - a_5*b_3*d_4*e_2*f_1 - a_5*b_4*d_1*e_2*f_3 + a_5*b_4*d_1*e_3*f_2 + a_5*b_4*d_2*e_1*f_3 - a_5*b_4*d_2*e_3*f_1 - a_5*b_4*d_3*e_1*f_2 + a_5*b_4*d_3*e_2*f_1)/martix_value
              Inverse_A(4,6)=(a_1*b_2*c_3*e_4*f_5 - a_1*b_2*c_3*e_5*f_4 - a_1*b_2*c_4*e_3*f_5 + a_1*b_2*c_4*e_5*f_3 + a_1*b_2*c_5*e_3*f_4 - a_1*b_2*c_5*e_4*f_3 - a_1*b_3*c_2*e_4*f_5 + a_1*b_3*c_2*e_5*f_4 + a_1*b_3*c_4*e_2*f_5 - a_1*b_3*c_4*e_5*f_2 - a_1*b_3*c_5*e_2*f_4 + a_1*b_3*c_5*e_4*f_2 + a_1*b_4*c_2*e_3*f_5 - a_1*b_4*c_2*e_5*f_3 - a_1*b_4*c_3*e_2*f_5 + a_1*b_4*c_3*e_5*f_2 + a_1*b_4*c_5*e_2*f_3 - a_1*b_4*c_5*e_3*f_2 - a_1*b_5*c_2*e_3*f_4 + a_1*b_5*c_2*e_4*f_3 + a_1*b_5*c_3*e_2*f_4 - a_1*b_5*c_3*e_4*f_2 - a_1*b_5*c_4*e_2*f_3 + a_1*b_5*c_4*e_3*f_2 - a_2*b_1*c_3*e_4*f_5 + a_2*b_1*c_3*e_5*f_4 + a_2*b_1*c_4*e_3*f_5 - a_2*b_1*c_4*e_5*f_3 - a_2*b_1*c_5*e_3*f_4 + a_2*b_1*c_5*e_4*f_3 + a_2*b_3*c_1*e_4*f_5 - a_2*b_3*c_1*e_5*f_4 - a_2*b_3*c_4*e_1*f_5 + a_2*b_3*c_4*e_5*f_1 + a_2*b_3*c_5*e_1*f_4 - a_2*b_3*c_5*e_4*f_1 - a_2*b_4*c_1*e_3*f_5 + a_2*b_4*c_1*e_5*f_3 + a_2*b_4*c_3*e_1*f_5 - a_2*b_4*c_3*e_5*f_1 - a_2*b_4*c_5*e_1*f_3 + a_2*b_4*c_5*e_3*f_1 + a_2*b_5*c_1*e_3*f_4 - a_2*b_5*c_1*e_4*f_3 - a_2*b_5*c_3*e_1*f_4 + a_2*b_5*c_3*e_4*f_1 + a_2*b_5*c_4*e_1*f_3 - a_2*b_5*c_4*e_3*f_1 + a_3*b_1*c_2*e_4*f_5 - a_3*b_1*c_2*e_5*f_4 - a_3*b_1*c_4*e_2*f_5 + a_3*b_1*c_4*e_5*f_2 + a_3*b_1*c_5*e_2*f_4 - a_3*b_1*c_5*e_4*f_2 - a_3*b_2*c_1*e_4*f_5 + a_3*b_2*c_1*e_5*f_4 + a_3*b_2*c_4*e_1*f_5 - a_3*b_2*c_4*e_5*f_1 - a_3*b_2*c_5*e_1*f_4 + a_3*b_2*c_5*e_4*f_1 + a_3*b_4*c_1*e_2*f_5 - a_3*b_4*c_1*e_5*f_2 - a_3*b_4*c_2*e_1*f_5 + a_3*b_4*c_2*e_5*f_1 + a_3*b_4*c_5*e_1*f_2 - a_3*b_4*c_5*e_2*f_1 - a_3*b_5*c_1*e_2*f_4 + a_3*b_5*c_1*e_4*f_2 + a_3*b_5*c_2*e_1*f_4 - a_3*b_5*c_2*e_4*f_1 - a_3*b_5*c_4*e_1*f_2 + a_3*b_5*c_4*e_2*f_1 - a_4*b_1*c_2*e_3*f_5 + a_4*b_1*c_2*e_5*f_3 + a_4*b_1*c_3*e_2*f_5 - a_4*b_1*c_3*e_5*f_2 - a_4*b_1*c_5*e_2*f_3 + a_4*b_1*c_5*e_3*f_2 + a_4*b_2*c_1*e_3*f_5 - a_4*b_2*c_1*e_5*f_3 - a_4*b_2*c_3*e_1*f_5 + a_4*b_2*c_3*e_5*f_1 + a_4*b_2*c_5*e_1*f_3 - a_4*b_2*c_5*e_3*f_1 - a_4*b_3*c_1*e_2*f_5 + a_4*b_3*c_1*e_5*f_2 + a_4*b_3*c_2*e_1*f_5 - a_4*b_3*c_2*e_5*f_1 - a_4*b_3*c_5*e_1*f_2 + a_4*b_3*c_5*e_2*f_1 + a_4*b_5*c_1*e_2*f_3 - a_4*b_5*c_1*e_3*f_2 - a_4*b_5*c_2*e_1*f_3 + a_4*b_5*c_2*e_3*f_1 + a_4*b_5*c_3*e_1*f_2 - a_4*b_5*c_3*e_2*f_1 + a_5*b_1*c_2*e_3*f_4 - a_5*b_1*c_2*e_4*f_3 - a_5*b_1*c_3*e_2*f_4 + a_5*b_1*c_3*e_4*f_2 + a_5*b_1*c_4*e_2*f_3 - a_5*b_1*c_4*e_3*f_2 - a_5*b_2*c_1*e_3*f_4 + a_5*b_2*c_1*e_4*f_3 + a_5*b_2*c_3*e_1*f_4 - a_5*b_2*c_3*e_4*f_1 - a_5*b_2*c_4*e_1*f_3 + a_5*b_2*c_4*e_3*f_1 + a_5*b_3*c_1*e_2*f_4 - a_5*b_3*c_1*e_4*f_2 - a_5*b_3*c_2*e_1*f_4 + a_5*b_3*c_2*e_4*f_1 + a_5*b_3*c_4*e_1*f_2 - a_5*b_3*c_4*e_2*f_1 - a_5*b_4*c_1*e_2*f_3 + a_5*b_4*c_1*e_3*f_2 + a_5*b_4*c_2*e_1*f_3 - a_5*b_4*c_2*e_3*f_1 - a_5*b_4*c_3*e_1*f_2 + a_5*b_4*c_3*e_2*f_1)/martix_value
              Inverse_A(5,6)=-(a_1*b_2*c_3*d_4*f_5 - a_1*b_2*c_3*d_5*f_4 - a_1*b_2*c_4*d_3*f_5 + a_1*b_2*c_4*d_5*f_3 + a_1*b_2*c_5*d_3*f_4 - a_1*b_2*c_5*d_4*f_3 - a_1*b_3*c_2*d_4*f_5 + a_1*b_3*c_2*d_5*f_4 + a_1*b_3*c_4*d_2*f_5 - a_1*b_3*c_4*d_5*f_2 - a_1*b_3*c_5*d_2*f_4 + a_1*b_3*c_5*d_4*f_2 + a_1*b_4*c_2*d_3*f_5 - a_1*b_4*c_2*d_5*f_3 - a_1*b_4*c_3*d_2*f_5 + a_1*b_4*c_3*d_5*f_2 + a_1*b_4*c_5*d_2*f_3 - a_1*b_4*c_5*d_3*f_2 - a_1*b_5*c_2*d_3*f_4 + a_1*b_5*c_2*d_4*f_3 + a_1*b_5*c_3*d_2*f_4 - a_1*b_5*c_3*d_4*f_2 - a_1*b_5*c_4*d_2*f_3 + a_1*b_5*c_4*d_3*f_2 - a_2*b_1*c_3*d_4*f_5 + a_2*b_1*c_3*d_5*f_4 + a_2*b_1*c_4*d_3*f_5 - a_2*b_1*c_4*d_5*f_3 - a_2*b_1*c_5*d_3*f_4 + a_2*b_1*c_5*d_4*f_3 + a_2*b_3*c_1*d_4*f_5 - a_2*b_3*c_1*d_5*f_4 - a_2*b_3*c_4*d_1*f_5 + a_2*b_3*c_4*d_5*f_1 + a_2*b_3*c_5*d_1*f_4 - a_2*b_3*c_5*d_4*f_1 - a_2*b_4*c_1*d_3*f_5 + a_2*b_4*c_1*d_5*f_3 + a_2*b_4*c_3*d_1*f_5 - a_2*b_4*c_3*d_5*f_1 - a_2*b_4*c_5*d_1*f_3 + a_2*b_4*c_5*d_3*f_1 + a_2*b_5*c_1*d_3*f_4 - a_2*b_5*c_1*d_4*f_3 - a_2*b_5*c_3*d_1*f_4 + a_2*b_5*c_3*d_4*f_1 + a_2*b_5*c_4*d_1*f_3 - a_2*b_5*c_4*d_3*f_1 + a_3*b_1*c_2*d_4*f_5 - a_3*b_1*c_2*d_5*f_4 - a_3*b_1*c_4*d_2*f_5 + a_3*b_1*c_4*d_5*f_2 + a_3*b_1*c_5*d_2*f_4 - a_3*b_1*c_5*d_4*f_2 - a_3*b_2*c_1*d_4*f_5 + a_3*b_2*c_1*d_5*f_4 + a_3*b_2*c_4*d_1*f_5 - a_3*b_2*c_4*d_5*f_1 - a_3*b_2*c_5*d_1*f_4 + a_3*b_2*c_5*d_4*f_1 + a_3*b_4*c_1*d_2*f_5 - a_3*b_4*c_1*d_5*f_2 - a_3*b_4*c_2*d_1*f_5 + a_3*b_4*c_2*d_5*f_1 + a_3*b_4*c_5*d_1*f_2 - a_3*b_4*c_5*d_2*f_1 - a_3*b_5*c_1*d_2*f_4 + a_3*b_5*c_1*d_4*f_2 + a_3*b_5*c_2*d_1*f_4 - a_3*b_5*c_2*d_4*f_1 - a_3*b_5*c_4*d_1*f_2 + a_3*b_5*c_4*d_2*f_1 - a_4*b_1*c_2*d_3*f_5 + a_4*b_1*c_2*d_5*f_3 + a_4*b_1*c_3*d_2*f_5 - a_4*b_1*c_3*d_5*f_2 - a_4*b_1*c_5*d_2*f_3 + a_4*b_1*c_5*d_3*f_2 + a_4*b_2*c_1*d_3*f_5 - a_4*b_2*c_1*d_5*f_3 - a_4*b_2*c_3*d_1*f_5 + a_4*b_2*c_3*d_5*f_1 + a_4*b_2*c_5*d_1*f_3 - a_4*b_2*c_5*d_3*f_1 - a_4*b_3*c_1*d_2*f_5 + a_4*b_3*c_1*d_5*f_2 + a_4*b_3*c_2*d_1*f_5 - a_4*b_3*c_2*d_5*f_1 - a_4*b_3*c_5*d_1*f_2 + a_4*b_3*c_5*d_2*f_1 + a_4*b_5*c_1*d_2*f_3 - a_4*b_5*c_1*d_3*f_2 - a_4*b_5*c_2*d_1*f_3 + a_4*b_5*c_2*d_3*f_1 + a_4*b_5*c_3*d_1*f_2 - a_4*b_5*c_3*d_2*f_1 + a_5*b_1*c_2*d_3*f_4 - a_5*b_1*c_2*d_4*f_3 - a_5*b_1*c_3*d_2*f_4 + a_5*b_1*c_3*d_4*f_2 + a_5*b_1*c_4*d_2*f_3 - a_5*b_1*c_4*d_3*f_2 - a_5*b_2*c_1*d_3*f_4 + a_5*b_2*c_1*d_4*f_3 + a_5*b_2*c_3*d_1*f_4 - a_5*b_2*c_3*d_4*f_1 - a_5*b_2*c_4*d_1*f_3 + a_5*b_2*c_4*d_3*f_1 + a_5*b_3*c_1*d_2*f_4 - a_5*b_3*c_1*d_4*f_2 - a_5*b_3*c_2*d_1*f_4 + a_5*b_3*c_2*d_4*f_1 + a_5*b_3*c_4*d_1*f_2 - a_5*b_3*c_4*d_2*f_1 - a_5*b_4*c_1*d_2*f_3 + a_5*b_4*c_1*d_3*f_2 + a_5*b_4*c_2*d_1*f_3 - a_5*b_4*c_2*d_3*f_1 - a_5*b_4*c_3*d_1*f_2 + a_5*b_4*c_3*d_2*f_1)/martix_value
              Inverse_A(6,6)=(a_1*b_2*c_3*d_4*e_5 - a_1*b_2*c_3*d_5*e_4 - a_1*b_2*c_4*d_3*e_5 + a_1*b_2*c_4*d_5*e_3 + a_1*b_2*c_5*d_3*e_4 - a_1*b_2*c_5*d_4*e_3 - a_1*b_3*c_2*d_4*e_5 + a_1*b_3*c_2*d_5*e_4 + a_1*b_3*c_4*d_2*e_5 - a_1*b_3*c_4*d_5*e_2 - a_1*b_3*c_5*d_2*e_4 + a_1*b_3*c_5*d_4*e_2 + a_1*b_4*c_2*d_3*e_5 - a_1*b_4*c_2*d_5*e_3 - a_1*b_4*c_3*d_2*e_5 + a_1*b_4*c_3*d_5*e_2 + a_1*b_4*c_5*d_2*e_3 - a_1*b_4*c_5*d_3*e_2 - a_1*b_5*c_2*d_3*e_4 + a_1*b_5*c_2*d_4*e_3 + a_1*b_5*c_3*d_2*e_4 - a_1*b_5*c_3*d_4*e_2 - a_1*b_5*c_4*d_2*e_3 + a_1*b_5*c_4*d_3*e_2 - a_2*b_1*c_3*d_4*e_5 + a_2*b_1*c_3*d_5*e_4 + a_2*b_1*c_4*d_3*e_5 - a_2*b_1*c_4*d_5*e_3 - a_2*b_1*c_5*d_3*e_4 + a_2*b_1*c_5*d_4*e_3 + a_2*b_3*c_1*d_4*e_5 - a_2*b_3*c_1*d_5*e_4 - a_2*b_3*c_4*d_1*e_5 + a_2*b_3*c_4*d_5*e_1 + a_2*b_3*c_5*d_1*e_4 - a_2*b_3*c_5*d_4*e_1 - a_2*b_4*c_1*d_3*e_5 + a_2*b_4*c_1*d_5*e_3 + a_2*b_4*c_3*d_1*e_5 - a_2*b_4*c_3*d_5*e_1 - a_2*b_4*c_5*d_1*e_3 + a_2*b_4*c_5*d_3*e_1 + a_2*b_5*c_1*d_3*e_4 - a_2*b_5*c_1*d_4*e_3 - a_2*b_5*c_3*d_1*e_4 + a_2*b_5*c_3*d_4*e_1 + a_2*b_5*c_4*d_1*e_3 - a_2*b_5*c_4*d_3*e_1 + a_3*b_1*c_2*d_4*e_5 - a_3*b_1*c_2*d_5*e_4 - a_3*b_1*c_4*d_2*e_5 + a_3*b_1*c_4*d_5*e_2 + a_3*b_1*c_5*d_2*e_4 - a_3*b_1*c_5*d_4*e_2 - a_3*b_2*c_1*d_4*e_5 + a_3*b_2*c_1*d_5*e_4 + a_3*b_2*c_4*d_1*e_5 - a_3*b_2*c_4*d_5*e_1 - a_3*b_2*c_5*d_1*e_4 + a_3*b_2*c_5*d_4*e_1 + a_3*b_4*c_1*d_2*e_5 - a_3*b_4*c_1*d_5*e_2 - a_3*b_4*c_2*d_1*e_5 + a_3*b_4*c_2*d_5*e_1 + a_3*b_4*c_5*d_1*e_2 - a_3*b_4*c_5*d_2*e_1 - a_3*b_5*c_1*d_2*e_4 + a_3*b_5*c_1*d_4*e_2 + a_3*b_5*c_2*d_1*e_4 - a_3*b_5*c_2*d_4*e_1 - a_3*b_5*c_4*d_1*e_2 + a_3*b_5*c_4*d_2*e_1 - a_4*b_1*c_2*d_3*e_5 + a_4*b_1*c_2*d_5*e_3 + a_4*b_1*c_3*d_2*e_5 - a_4*b_1*c_3*d_5*e_2 - a_4*b_1*c_5*d_2*e_3 + a_4*b_1*c_5*d_3*e_2 + a_4*b_2*c_1*d_3*e_5 - a_4*b_2*c_1*d_5*e_3 - a_4*b_2*c_3*d_1*e_5 + a_4*b_2*c_3*d_5*e_1 + a_4*b_2*c_5*d_1*e_3 - a_4*b_2*c_5*d_3*e_1 - a_4*b_3*c_1*d_2*e_5 + a_4*b_3*c_1*d_5*e_2 + a_4*b_3*c_2*d_1*e_5 - a_4*b_3*c_2*d_5*e_1 - a_4*b_3*c_5*d_1*e_2 + a_4*b_3*c_5*d_2*e_1 + a_4*b_5*c_1*d_2*e_3 - a_4*b_5*c_1*d_3*e_2 - a_4*b_5*c_2*d_1*e_3 + a_4*b_5*c_2*d_3*e_1 + a_4*b_5*c_3*d_1*e_2 - a_4*b_5*c_3*d_2*e_1 + a_5*b_1*c_2*d_3*e_4 - a_5*b_1*c_2*d_4*e_3 - a_5*b_1*c_3*d_2*e_4 + a_5*b_1*c_3*d_4*e_2 + a_5*b_1*c_4*d_2*e_3 - a_5*b_1*c_4*d_3*e_2 - a_5*b_2*c_1*d_3*e_4 + a_5*b_2*c_1*d_4*e_3 + a_5*b_2*c_3*d_1*e_4 - a_5*b_2*c_3*d_4*e_1 - a_5*b_2*c_4*d_1*e_3 + a_5*b_2*c_4*d_3*e_1 + a_5*b_3*c_1*d_2*e_4 - a_5*b_3*c_1*d_4*e_2 - a_5*b_3*c_2*d_1*e_4 + a_5*b_3*c_2*d_4*e_1 + a_5*b_3*c_4*d_1*e_2 - a_5*b_3*c_4*d_2*e_1 - a_5*b_4*c_1*d_2*e_3 + a_5*b_4*c_1*d_3*e_2 + a_5*b_4*c_2*d_1*e_3 - a_5*b_4*c_2*d_3*e_1 - a_5*b_4*c_3*d_1*e_2 + a_5*b_4*c_3*d_2*e_1)/martix_value
          
          end if
        !********************************************************************************************************
                        
        else
         
            write(*,*) "the matrix is not less than six order !"
         
        end if
     

        !Check the inversion opeartion is right or not
        if(martix_value/=0) then
            result_matrix=Inverse_A
        else
            check_error=0
            write(*,*) "the compute_liner_equations is wrong! "
        end if
      
    end subroutine matrix_inversion_less_six
    !**********************************************************************************
    
    !**********************************************************************************
    !Kalman filter
    Function Kalman_filter(Total_Sampling,current_time_number,input_data,Q,R,x_0,P_0)
  
      ! Variables
      integer::current_time_number,Total_Sampling
      real(kind=8)::Kalman_filter,Q,R,x_0,P_0
      real(kind=8),dimension(Total_Sampling)::input_data
  
      integer::i,j
      real(kind=8),dimension(current_time_number)::K
      real(kind=8),dimension(current_time_number)::X
      real(kind=8),dimension(current_time_number)::P
  
      X(1) =x_0; 
      P(1) =P_0; 
  
      do i=2,current_time_number
              
          K(i)=P(i-1)/(P(i-1)+R); 
          X(i)=X(i-1)+K(i)*(input_data(i)-X(i-1)); 
          P(i)=P(i-1)-K(i)*P(i-1)+Q; 
      
      end do
      
      Kalman_filter=X(current_time_number)

    end Function Kalman_filter
    !**********************************************************************************
    
    
    
    !**********************************************************************************
    SUBROUTINE BRINV(A,N,L) ! Calculate the matrix inversion

    implicit none
    
    integer::N,L                                                 ! Input matrix dimension: n*n
    real(kind=8),dimension(N,N)::A                               ! Input matrix
    integer,dimension(N)::IS,JS                                  ! Output result matrix
    real(kind=8)::T,D
    integer::I,J,K

    L=1
    DO 100 K=1,N
      D=0.0
      DO 10 I=K,N
      DO 10 J=K,N
        IF (ABS(A(I,J)).GT.D) THEN
          D=ABS(A(I,J))
          IS(K)=I
          JS(K)=J
        END IF
10    CONTINUE
      IF (D+1.0.EQ.1.0) THEN
        L=0
        WRITE(*,20)
        RETURN
      END IF
20    FORMAT(1X,'ERR**NOT INV')
      DO 30 J=1,N
        T=A(K,J)
        A(K,J)=A(IS(K),J)
        A(IS(K),J)=T
30    CONTINUE
      DO 40 I=1,N
        T=A(I,K)
        A(I,K)=A(I,JS(K))
        A(I,JS(K))=T
40    CONTINUE
      A(K,K)=1/A(K,K)
      DO 50 J=1,N
        IF (J.NE.K) THEN
          A(K,J)=A(K,J)*A(K,K)
        END IF
50    CONTINUE
      DO 70 I=1,N
        IF (I.NE.K) THEN
          DO 60 J=1,N
            IF (J.NE.K) THEN
              A(I,J)=A(I,J)-A(I,K)*A(K,J)
            END IF
60        CONTINUE
        END IF
70    CONTINUE
      DO 80 I=1,N
        IF (I.NE.K) THEN
          A(I,K)=-A(I,K)*A(K,K)
        END IF
80    CONTINUE
100 CONTINUE
    DO 130 K=N,1,-1
      DO 110 J=1,N
        T=A(K,J)
        A(K,J)=A(JS(K),J)
        A(JS(K),J)=T
110   CONTINUE
      DO 120 I=1,N
        T=A(I,K)
        A(I,K)=A(I,IS(K))
        A(I,IS(K))=T
120   CONTINUE
130 CONTINUE
    RETURN
    END SUBROUTINE BRINV
    !**********************************************************************************



    !**********************************************************************************************************
    !Standard gaussian elimination slover
    !Input data: (a)equation_dim---the dimension of the linear system--------n
    !            (b)extern_dim---the dimension of the linear system----------n+1
    !            (c)coefficient_matrix---the coefficient matrix (n*n)
    !            (d)right_hand_vector---the right hand matrix (n*1)
    !            (e)result_x---the return result array (n*1)
    !**********************************************************************************************************

    SUBROUTINE Gaussian_Elimination_Solver(equation_dim,extern_dim,coefficient_matrix,right_hand_vector,result_x)
    
        implicit none
    
        !------------------------------------------------------------------------------------------------------
        ! Variables input 
        integer,intent(in)::equation_dim                                                        ! equation dim
        integer,intent(in)::extern_dim                                                          ! equation dim+1
        real(kind=8),intent(in),dimension(equation_dim,equation_dim)::coefficient_matrix        ! coefficient matrix
        real(kind=8),intent(in),dimension(equation_dim)::right_hand_vector                      ! RHS matrix
        real(kind=8),intent(out),dimension(equation_dim)::result_x                              ! return results matrix        

        ! Variables in Gaussian_Elimination_Solver
        integer::i,j,k,m,n                                                                      ! loop control variable
        real(kind=8),dimension(equation_dim,extern_dim)::extern_matrix                          ! the matrix contain coefficient and rhs
        real(kind=8),dimension(equation_dim)::temp_list                                         ! tempture list array 
        real(kind=8),dimension(extern_dim)::temp_line                                           ! tempture line array 
        real(kind=8)::max_value                                                                 ! the value of the max in list
        integer,dimension(1)::max_location                                                      ! the location of the max in list
        real(kind=8)::temp                                                                      ! temp value
        real(kind=8)::elimination
        !------------------------------------------------------------------------------------------------------


        !------------------------------------------------------------------------------------------------------
        ! Body of Guassian_Elimination
        
        result_x=0                                                                              ! initial the return results array
        
        ! Assemble the extern matrix which contains the coefficient and right hand side
        do i=1,equation_dim
            extern_matrix(:,i)=coefficient_matrix(:,i)
        end do
        extern_matrix(:,extern_dim)=right_hand_vector(:)


        do j=1,equation_dim
            
            ! Save the abs of list of the extern array
            temp_list=-1
            do i=j,equation_dim                                  ! only need to start from j
                temp_list(i)=abs(coefficient_matrix(i,j))
            end do
            
            ! Search the max abs value and its loaction in the array
            max_location=MaxLoc(temp_list)                       ! as the return is an array in MacLoc
            m=max_location(1)
            max_value=coefficient_matrix(m,j)                    ! the ture data is in the original matrix as the temp list is just used for find the largest value
            
            ! Check the coefficient array
            if(max_value==0) then
                write(*,*) "There are something wrong in the linear system---coefficient matrix is singular"
            end if
                
            if(m/=j) then                                        ! need the array moving
                   
                ! change the line       
                temp_line=extern_matrix(j,:)                     ! save j line
                extern_matrix(j,:)=extern_matrix(m,:)            ! change j and m line
                extern_matrix(m,:)=temp_line                     ! change j and m line
                
            end if
            
            ! Doing the element elimination
            do i=j+1,equation_dim                                ! Doing elimination in every line
               elimination=-extern_matrix(i,j)/extern_matrix(j,j)
                    
               ! Doing adding action in every list number
               do k=1,extern_dim
                  extern_matrix(i,k)=extern_matrix(i,k)+extern_matrix(j,k)*elimination
               end do
             end do   
            
        enddo
        
        ! check the coefficient array
        if(extern_matrix(equation_dim,equation_dim)==0) then
            write(*,*) "There are something wrong in the linear system"
        endif

        ! compute the results
        do i=equation_dim,1,-1
            
            temp=0
            do j=i,equation_dim
               temp=temp+extern_matrix(i,j)*result_x(j)
            end do
            
            result_x(i)=(extern_matrix(i,extern_dim)-temp)/extern_matrix(i,i)
            
        enddo
        !------------------------------------------------------------------------------------------------------

        
    END SUBROUTINE Gaussian_Elimination_Solver
    !**********************************************************************************************************



    !**********************************************************************************************************
    !This is right
    !cline angle (anti-clockwise) free slip boundary condition 
    !coordinates transfer from global to local, apply the boundary condition , then transfer from local to gobal
    subroutine cline_free_slip(particle_velocity_x,particle_velocity_y,beta,ghost_velocity_x,ghost_velocity_y)

        implicit none

        ! Variables
              
        integer::i,j,k,l                                                       !Variables for looping
        real(kind=8),intent(in)::particle_velocity_x,particle_velocity_y       !Particle velocity in gobal coordinate
        real(kind=8),intent(in)::beta                                          !Incline angle in degree
        real(kind=8),intent(out)::ghost_velocity_x,ghost_velocity_y            !Ghost particle velocity (output) 

        !Variables for subroutine
        real(kind=8)::particle_u,particle_v                                    !Particle velocity in local coordinate
        real(kind=8)::ghost_u,ghost_v                                          !velocity of ghost particle
        real(kind=8)::theta                                                    !angle
        
        !body of subroutine cline_free_slip
        !transfer angle from degree to radian
        theta=beta*PI/180.0
        
        !transfer velocity from gobal to local
        particle_u=particle_velocity_x*cos(theta)+particle_velocity_y*sin(theta)
        particle_v=-particle_velocity_x*sin(theta)+particle_velocity_y*cos(theta)
        
        !Free slip 
        ghost_u=particle_u
        ghost_v=-particle_v
        
        !transfer velocity from local to gobal
        ghost_velocity_x=ghost_u*cos(theta)-ghost_v*sin(theta)
        ghost_velocity_y=ghost_u*sin(theta)+ghost_v*cos(theta)
    
    end subroutine cline_free_slip
    !**********************************************************************************************************
    
    !**********************************************************************************************************
    !cline angle (anti-clockwise) no slip boundary condition 
    subroutine cline_no_slip(particle_velocity_x,particle_velocity_y,beta,ghost_velocity_x,ghost_velocity_y)

        implicit none

        ! Variables
              
        integer::i,j,k,l                                                       !Variables for looping
        real(kind=8),intent(in)::particle_velocity_x,particle_velocity_y       !Particle velocity in gobal coordinate
        real(kind=8),intent(in)::beta                                          !Incline angle in degree
        real(kind=8),intent(out)::ghost_velocity_x,ghost_velocity_y            !Ghost particle velocity (output) 

        !Variables for subroutine
        real(kind=8)::particle_u,particle_v                                    !Particle velocity in local coordinate
        real(kind=8)::ghost_u,ghost_v                                          !velocity of ghost particle
        real(kind=8)::theta                                                    !angle
     
        !body of subroutine cline_no_slip
        !transfer velocity from gobal to local
        particle_u=particle_velocity_x*cos(theta)+particle_velocity_y*sin(theta)
        particle_v=-particle_velocity_x*sin(theta)+particle_velocity_y*cos(theta)
        
        !No slip 
        ghost_u=-particle_u
        ghost_v=-particle_v
        
        !transfer velocity from gobal to local
        ghost_velocity_x=ghost_u*cos(theta)-ghost_v*sin(theta)
        ghost_velocity_y=ghost_u*sin(theta)+ghost_v*cos(theta)
    
    end subroutine cline_no_slip
    !**********************************************************************************************************


    !**********************************************************************************************************
    !Free slip boundary condition for circle
    subroutine radius_free_slip(particle_position_x,particle_position_y,particle_velocity_x,particle_velocity_y,ghost_velocity_x,ghost_velocity_y)

        implicit none

        ! Variables      
        integer::i,j,k,l                                                       !Variables for looping
        real(kind=8),intent(in)::particle_position_x,particle_position_y       !Particle position in gobal coordinate
        real(kind=8),intent(in)::particle_velocity_x,particle_velocity_y       !Particle velocity in gobal coordinate
        real(kind=8),intent(out)::ghost_velocity_x,ghost_velocity_y            !Ghost particle velocity (output) 

        !Variables for subroutine
        real(kind=8)::particle_u,particle_v                                    !Particle velocity in local coordinate
        real(kind=8)::ghost_u,ghost_v                                          !velocity of ghost particle
        real(kind=8)::theta                                                    !angle
        
        !body of subroutine cline_free_slip
        !calculate the radius angle
        theta=atan2((0.5-particle_position_y),(0.5-particle_position_x))
        
        !transfer velocity from gobal to local
        particle_u=-particle_velocity_x*sin(theta)+particle_velocity_y*cos(theta)
        particle_v=particle_velocity_x*cos(theta)+particle_velocity_y*sin(theta)
        
        !Free slip
        ghost_u=particle_u
        ghost_v=-particle_v
        
        !transfer velocity from gobal to local
        ghost_velocity_x=-ghost_u*sin(theta)+ghost_v*cos(theta)
        ghost_velocity_y=ghost_u*cos(theta)+ghost_v*sin(theta)
    
    end subroutine radius_free_slip
    !**********************************************************************************************************


    !**********************************************************************************************************
    !No slip boundary condition for circle
    subroutine radius_no_slip(particle_position_x,particle_position_y,particle_velocity_x,particle_velocity_y,ghost_velocity_x,ghost_velocity_y)

        implicit none

        ! Variables      
        integer::i,j,k,l                                                       !Variables for looping
        real(kind=8),intent(in)::particle_position_x,particle_position_y       !Particle position in gobal coordinate
        real(kind=8),intent(in)::particle_velocity_x,particle_velocity_y       !Particle velocity in gobal coordinate
        real(kind=8),intent(out)::ghost_velocity_x,ghost_velocity_y            !Ghost particle velocity (output) 

        !Variables for subroutine
        real(kind=8)::particle_u,particle_v                                    !Particle velocity in local coordinate
        real(kind=8)::ghost_u,ghost_v                                          !velocity of ghost particle
        real(kind=8)::theta                                                    !angle

        !body of subroutine cline_free_slip
        !calculate the radius angle
        theta=atan2((0.5-particle_position_y),(0.5-particle_position_x))
        
        !transfer velocity from gobal to local
        particle_u=-particle_velocity_x*sin(theta)+particle_velocity_y*cos(theta)
        particle_v=particle_velocity_x*cos(theta)+particle_velocity_y*sin(theta)
        
        !No slip
        ghost_u=-particle_u
        ghost_v=-particle_v
        
        !transfer velocity from gobal to local
        ghost_velocity_x=-ghost_u*sin(theta)+ghost_v*cos(theta)
        ghost_velocity_y=ghost_u*cos(theta)+ghost_v*sin(theta)
    
    end subroutine radius_no_slip
    !**********************************************************************************************************























    !==========================================================================================================
    !       ┌───┐   ┌───┬───┬───┬───┐ ┌───┬───┬───┬───┐ ┌───┬───┬───┬───┐ ┌───┬───┬───┐
    !       │Esc│   │ F1│ F2│ F3│ F4│ │ F5│ F6│ F7│ F8│ │ F9│F10│F11│F12│ │P/S│S L│P/B│  ┌┐    ┌┐    ┌┐
    !       └───┘   └───┴───┴───┴───┘ └───┴───┴───┴───┘ └───┴───┴───┴───┘ └───┴───┴───┘  └┘    └┘    └┘
    !       ┌───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───────┐ ┌───┬───┬───┐ ┌───┬───┬───┬───┐
    !       │~ `│! 1│@ 2│# 3│$ 4│% 5│^ 6│& 7│* 8│( 9│) 0│_ -│+ =│ BacSp │ │Ins│Hom│PUp│ │N L│ / │ * │ - │
    !       ├───┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─────┤ ├───┼───┼───┤ ├───┼───┼───┼───┤
    !       │ Tab │ Q │ W │ E │ R │ T │ Y │ U │ I │ O │ P │{ [│} ]│ | \ │ │Del│End│PDn│ │ 7 │ 8 │ 9 │   │
    !       ├─────┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴─────┤ └───┴───┴───┘ ├───┼───┼───┤ + │
    !       │ Caps │ A │ S │ D │ F │ G │ H │ J │ K │ L │: ;│" '│ Enter  │               │ 4 │ 5 │ 6 │   │
    !       ├──────┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴────────┤     ┌───┐     ├───┼───┼───┼───┤
    !       │ Shift  │ Z │ X │ C │ V │ B │ N │ M │< ,│> .│? /│  Shift   │     │ ↑ │     │ 1 │ 2 │ 3 │   │
    !       ├─────┬──┴─┬─┴──┬┴───┴───┴───┴───┴───┴──┬┴───┼───┴┬────┬────┤ ┌───┼───┼───┐ ├───┴───┼───┤ E││
    !       │ Ctrl│    │Alt │         Space         │ Alt│    │    │Ctrl│ │ ← │ ↓ │ → │ │   0   │ . │←─┘│
    !       └─────┴────┴────┴───────────────────────┴────┴────┴────┴────┘ └───┴───┴───┘ └───────┴───┴───┘
    !==========================================================================================================




    !==========================================================================================================
    subroutine Initialziting_Writing_File(File_index,File_name,ioerror)

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables from superior-subroutine
        integer,intent(in)::File_index                ! File index
        character(len=100),intent(in)::File_name      ! File name
        integer,intent(out)::ioerror                  ! File operation output value

        ! Variables in local subroutine
        !------------------------------------------------------------------------------------------------------

        ! Body of subroutine Initialziting_Writing_File
        open(unit=File_index,file=trim(adjustl(File_name)),status="replace",position="rewind",action="write",iostat=ioerror)    
        
        if(ioerror/=0) then
            write(*,*) "The open action of the "//trim(adjustl(File_name))//" is not successful!"
        end if

    end subroutine Initialziting_Writing_File
    !==========================================================================================================


    !==========================================================================================================
    subroutine Initialziting_Existed_Writing_File(File_index,File_name,ioerror)

        implicit none

        !------------------------------------------------------------------------------------------------------
        ! Variables from superior-subroutine
        integer,intent(in)::File_index                ! Input file index
        character(len=100),intent(in)::File_name      ! Input file name
        integer,intent(out)::ioerror                  ! File operation return value

        ! Variables in local subroutine
        Logical::file_exist                           ! File exist or not 
        !------------------------------------------------------------------------------------------------------

        ! Body of subroutine Initialziting_Reading_File

        ! Inquire the reading file exist or not
        Inquire(File=adjustl(trim(File_name)), Exist=file_exist)

        if(file_exist) then

            open(unit=File_index,file=trim(adjustl(File_name)),status="old",position="append",action="write",iostat=ioerror)       
            
            if(ioerror/=0) then
                write(*,*) "The open action of the "//trim(adjustl(File_name))//" is not successful!"
            end if

        else

            write(*,*) "The file "//trim(adjustl(File_name))//" is not exist!"

        endif

    end subroutine Initialziting_Existed_Writing_File
    !==========================================================================================================


    !==========================================================================================================
    subroutine Initialziting_Reading_File(File_index,File_name,ioerror)

        implicit none
        
        !------------------------------------------------------------------------------------------------------
        ! Variables from superior-subroutine
        integer,intent(in)::File_index                ! Input file index
        character(len=100),intent(in)::File_name      ! Input file name
        integer,intent(out)::ioerror                  ! File operation return value

        ! Variables in local subroutine
        Logical::File_exist                           ! File exist or not 
        !------------------------------------------------------------------------------------------------------

        ! Body of subroutine Initialziting_Reading_File

        ! Inquire the reading file exist or not
        Inquire(File=adjustl(trim(File_name)), Exist=file_exist)

        if(file_exist) then

            open(unit=File_index,file=trim(adjustl(File_name)),status="old",position="rewind",action="read",iostat=ioerror)       
            
            if(ioerror/=0) then
                write(*,*) "The open action of the "//trim(adjustl(File_name))//" is not successful!"
            end if

        else

            write(*,*) "The file "//trim(adjustl(File_name))//" is not exist!"

        endif    

    end subroutine Initialziting_Reading_File
    !==========================================================================================================













    !==========================================================================================================
end module Function_Module