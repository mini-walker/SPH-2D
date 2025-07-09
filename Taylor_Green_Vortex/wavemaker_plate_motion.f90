!*FROM THE GORING (1978),
!*XP(T)=H/K*[TANH(X(T))+TANH(K/D*LAMDA)],X(T)=K/D*(C*T-XP(T)-LAMDA)
!*K=SQRT(3H/4D),C=SQRT(G*(D+H))
!*UP(T)=C*H/D*(1/(COSH(X(T))**2.+H/D)
!*LAMDA=L/K*D, L=ARCCOSH(1/SQRT(EZ)),EZ=0.002
subroutine wavemaker_plate_motion(WMA,WMV,WMX,TTT,HEIGHT,H0,DT,G)
    
    implicit none
   
    !精度类型常量
    integer,parameter::s=selected_real_kind(p=6,r=37)                      !单精度常量
    integer,parameter::d=selected_real_kind(p=13)                          !双精度常量
    real(kind=8) :: WMA,WMV,WMX,TTT,HEIGHT,H0,DT,G
    real(kind=8) :: EZ,CL,CK,CK1,C,CLAMDA,CERR,XPI,XP,XT,ERR,UP,VP
      
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
	end subroutine wavemaker_plate_motion





