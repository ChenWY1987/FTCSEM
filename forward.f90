!================================================================
!================================================================
!	   1D CSAMT(Ex/t)、SOTEM(Vy\Vz) OCCMA直接反演程序
!      基于频率向时间域的转化时间 SOTEM正演
!================================================================
!      多模型 同时正演计算 实现灵敏度并行
!================================================================           
!================================================================  
    MODULE D1_CSAMTINV
	USE FREQUENCY_FORWARD
    CONTAINS
!=================================================================================
    SUBROUTINE D1_CSOT_MF(TEDM,DATE_XZ,POTNS,CDBLF,LSI,NF,NL,FJ,MODS,MODSR, &
	                      EH_F,EH_ERSY,MAX_FDN,LLL,MOD_RLZ,BTAS)
!=================================================================================
	IMPLICIT NONE
!=======================================================================
	INTEGER :: I,J,K,L,S,IL,LLL,NRPD,IHDS,IHDF,LH
	INTEGER :: NF,NL,NK,TEDM,POTNS
	INTEGER :: NN,NIVSN,MN,MAX_FDN,PSPA_LH,MOD_RLZ
	INTEGER :: NKS,NHS,YNL,DATE_XZ(6)
	REAL(8) :: MOD_ZYS,MOD_YYS
	REAL(8) :: BHS,PLNHS
   !----------------------------------------------------------------------------------------
	REAL(8) :: FJ(NF)
	TYPE(TYPE_CDINF)  :: CDBLF
	TYPE(TYPE_ISINF)  :: LSI(POTNS)
	TYPE(TYPE_MODINF) :: MODS(NL),MODF(NL),MODSR(NL)
	TYPE(TYPE_EHINF)  :: EH_F(NF),EH_ERSY(NF)
   !----------------------------------------------------------------------------------------
	TYPE(TYPE_EHINF),ALLOCATABLE  :: EH_PL(:),EH_PJ(:),EH_PK(:)
	REAL(8),ALLOCATABLE :: WD(:),WS(:),DELT(:)
	REAL(8),ALLOCATABLE :: JAC(:,:),ATA(:,:),SR(:)
	REAL(8),ALLOCATABLE :: DATSX0(:),DATSX1(:),DATSX2(:),DATSX3(:),DATSX4(:),DATSX5(:)
   !----------------------------------------------------------------------------------------
	REAL(8) :: WC,WCD,WCL,WC1,WC2,WC3,WC4,WC5,WC6,DRP,SMH 
	REAL(8) :: BTAS,SQSM,FCXS0,FCXS1,FCXS2
	REAL(8) :: TIME0,TIME1
	CHARACTER(100) :: OUT_SS,OUT_DS
!=====================================================================
    BHS=CDBLF%N
    I=INT( BHS/1000 )
	J=INT( ( BHS - I*1000 )/100 )
	K=INT( ( BHS - I*1000 - J*100 )/10 )
	L=INT(   BHS - I*1000 - J*100 - K*10 )
	OUT_SS(1:1)=CHAR(I+48)
	OUT_SS(2:2)=CHAR(J+48)
	OUT_SS(3:3)=CHAR(K+48)
	OUT_SS(4:4)=CHAR(L+48)
	OUT_SS(5:12)="拟合.TXT"
!=====================================================================
    NHS=0
	NKS=NF
	IF (DATE_XZ(1)==1) NHS=NHS + NF
	IF (DATE_XZ(2)==1) NHS=NHS + NF
	IF (DATE_XZ(3)==1) NHS=NHS + NF
	IF (DATE_XZ(4)==1) NHS=NHS + NF
	IF (DATE_XZ(5)==1) NHS=NHS + NF
	IF (DATE_XZ(6)==1) NHS=NHS + NF
	MOD_ZYS=0.5D0
	MOD_YYS=1.0D0
!====================================================================================
    ALLOCATE( DATSX0(NF),DATSX1(NF),DATSX2(NF),DATSX3(NF),DATSX4(NF),DATSX5(NF) )
    ALLOCATE( EH_PJ(NKS),DELT(NHS),WD(NHS),WS(NHS) )
!====================================================================================
	DO I=1,NKS
	IF (DATE_XZ(1)==0) DATSX0(I)=1.0D0
	IF (DATE_XZ(2)==0) DATSX1(I)=1.0D0
	IF (DATE_XZ(3)==0) DATSX2(I)=1.0D0
	IF (DATE_XZ(4)==0) DATSX3(I)=1.0D0
	IF (DATE_XZ(5)==0) DATSX4(I)=1.0D0
	IF (DATE_XZ(6)==0) DATSX5(I)=1.0D0
	IF (DATE_XZ(1)==1) DATSX0(I)=10.0D0/EH_F(I)%EX
	IF (DATE_XZ(2)==1) DATSX1(I)=10.0D0/EH_F(I)%EY
	IF (DATE_XZ(3)==1) DATSX2(I)=10.0D0/EH_F(I)%EZ
	IF (DATE_XZ(4)==1) DATSX3(I)=10.0D0/EH_F(I)%HX
	IF (DATE_XZ(5)==1) DATSX4(I)=10.0D0/EH_F(I)%HY
	IF (DATE_XZ(6)==1) DATSX5(I)=10.0D0/EH_F(I)%HZ
	IF (DATE_XZ(1)==0) EH_PJ(I)%EX=EH_F(I)%EX
	IF (DATE_XZ(2)==0) EH_PJ(I)%EY=EH_F(I)%EY
	IF (DATE_XZ(3)==0) EH_PJ(I)%EZ=EH_F(I)%EZ
	IF (DATE_XZ(4)==0) EH_PJ(I)%HX=EH_F(I)%HX
	IF (DATE_XZ(5)==0) EH_PJ(I)%HY=EH_F(I)%HY
	IF (DATE_XZ(6)==0) EH_PJ(I)%HZ=EH_F(I)%HZ
	IF (DATE_XZ(1)==1) EH_PJ(I)%EX=10.0D0
	IF (DATE_XZ(2)==1) EH_PJ(I)%EY=10.0D0
	IF (DATE_XZ(3)==1) EH_PJ(I)%EZ=10.0D0
	IF (DATE_XZ(4)==1) EH_PJ(I)%HX=10.0D0
	IF (DATE_XZ(5)==1) EH_PJ(I)%HY=10.0D0
	IF (DATE_XZ(6)==1) EH_PJ(I)%HZ=10.0D0
	END DO
!=====================================================================
    L=0
	DO J=1,6
	IF (DATE_XZ(J)==1) THEN
	L=L+1
	K=(L-1)*NKS
    IF (J==1) WD(K+1:K+NF)=EH_ERSY(1:NF)%EX
    IF (J==2) WD(K+1:K+NF)=EH_ERSY(1:NF)%EY
    IF (J==3) WD(K+1:K+NF)=EH_ERSY(1:NF)%EZ
    IF (J==4) WD(K+1:K+NF)=EH_ERSY(1:NF)%HX
    IF (J==5) WD(K+1:K+NF)=EH_ERSY(1:NF)%HY
    IF (J==6) WD(K+1:K+NF)=EH_ERSY(1:NF)%HZ
	END IF
	END DO
!====================================================
    DO I=1,NHS
	 IF ( WD(I)>99.9D0 ) WD(I)=99.9D0   !防止误差都为100.0导致WD都为0.0 重而出现奇异点
	 IF ( WD(I)< 0.0D0 ) WD(I)=0.0D0
     WD(I)=1.0D0- WD(I)/100.0D0
	END DO
!=====================================================================
    CALL CPU_TIME(TIME0)
!==========================================================================================
!******************************************************************************************
!==========================================================================================
    LH=0
	NN=1
    WCL=100.0D0
	DO 2000 NIVSN=1,MAX_FDN	               !开始反演
!==========================================================================================
!       对修改好的模型进行正演计算	!第一个正演
!=====================================================================
    ALLOCATE( EH_PL(NKS) )
!---------------------------------------------------------------------
	CALL D1_SOTEMZ(TEDM,CDBLF,LSI,POTNS,NKS,NL,FJ,MODS,EH_PL)
!---------------------------------------------------------------------
    DO I=1,NKS
	EH_PL(I)%EX=EH_PL(I)%EX*DATSX0(I)
	EH_PL(I)%EY=EH_PL(I)%EY*DATSX1(I)
	EH_PL(I)%EZ=EH_PL(I)%EZ*DATSX2(I)
	EH_PL(I)%HX=EH_PL(I)%HX*DATSX3(I)
	EH_PL(I)%HY=EH_PL(I)%HY*DATSX4(I)
	EH_PL(I)%HZ=EH_PL(I)%HZ*DATSX5(I)
	END DO
!==================================================
	OPEN(Unit=51,File="数据拟合\"//OUT_SS(1:12))
    IF (TEDM==1) THEN
	  WRITE(51,"(5X,A,1X,12(11X,A,1X),1X,6(8X,A))")"频率","SEx_RD","SEy_RD","SEz_RD","SHx_RD","SHy_RD","SHz_RD",  &
	                                                      "NEx_RD","NEy_RD","NEz_RD","NHx_RD","NHy_RD","NHz_RD",  &
														  "Ex误差","Ey误差","Ez误差","Hx误差","Hy误差","Hz误差"
	ELSEIF (TEDM==2) THEN
	  WRITE(51,"(5X,A,1X,12(11X,A,1X),1X,6(8X,A))")"时间","SEx_T","SEy_T","SEz_T","SVx_T","SVy_T","SVz_T",  &
	                                                      "NEx_T","NEy_T","NEz_T","NVx_T","NVy_T","NVz_T",  &
														  "Ex误差","Ey误差","Ez误差","Vx误差","Vy误差","Vz误差"
	END IF
!==================================================
    K=0
	DO I=1,NF
!============================================================================
       WC1=1.0d0 - DLOG10( EH_PL(I)%EX )/DLOG10( EH_PJ(I)%EX )
       WC2=1.0d0 - DLOG10( EH_PL(I)%EY )/DLOG10( EH_PJ(I)%EY )
       WC3=1.0d0 - DLOG10( EH_PL(I)%EZ )/DLOG10( EH_PJ(I)%EZ )
       WC4=1.0d0 - DLOG10( EH_PL(I)%HX )/DLOG10( EH_PJ(I)%HX )
       WC5=1.0d0 - DLOG10( EH_PL(I)%HY )/DLOG10( EH_PJ(I)%HY )
       WC6=1.0d0 - DLOG10( EH_PL(I)%HZ )/DLOG10( EH_PJ(I)%HZ )
	   WRITE(51,'(E12.4,12E18.10,6F14.7)')FJ(I), &
	    EH_PJ(I)%EX/DATSX0(I) ,EH_PJ(I)%EY/DATSX1(I) ,EH_PJ(I)%EZ/DATSX2(I) ,   &
	    EH_PJ(I)%HX/DATSX3(I) ,EH_PJ(I)%HY/DATSX4(I) ,EH_PJ(I)%HZ/DATSX5(I) ,   &
	    EH_PL(I)%EX/DATSX0(I) ,EH_PL(I)%EY/DATSX1(I) ,EH_PL(I)%EZ/DATSX2(I) ,   &
	    EH_PL(I)%HX/DATSX3(I) ,EH_PL(I)%HY/DATSX4(I) ,EH_PL(I)%HZ/DATSX5(I) ,   &
	    WC1,WC2,WC3,WC4,WC5,WC6
!============================================================================
       IF (DATE_XZ(1)==1) THEN
	    K=K+1
        WC1=DLOG10( EH_PJ(I)%EX ) - DLOG10( EH_PL(I)%EX )
        DELT(K)=WC1
	   END IF
       IF (DATE_XZ(2)==1) THEN
	    K=K+1
        WC1=DLOG10( EH_PJ(I)%EY ) - DLOG10( EH_PL(I)%EY )
        DELT(K)=WC1
	   END IF
       IF (DATE_XZ(3)==1) THEN
	    K=K+1
        WC1=DLOG10( EH_PJ(I)%EZ ) - DLOG10( EH_PL(I)%EZ )
        DELT(K)=WC1
	   END IF
       IF (DATE_XZ(4)==1) THEN
	    K=K+1
        WC1=DLOG10( EH_PJ(I)%HX ) - DLOG10( EH_PL(I)%HX )
        DELT(K)=WC1
	   END IF
       IF (DATE_XZ(5)==1) THEN
	    K=K+1
        WC1=DLOG10( EH_PJ(I)%HY ) - DLOG10( EH_PL(I)%HY )
        DELT(K)=WC1
	   END IF
       IF (DATE_XZ(6)==1) THEN
	    K=K+1
        WC1=DLOG10( EH_PJ(I)%HZ ) - DLOG10( EH_PL(I)%HZ )
        DELT(K)=WC1
	   END IF
!============================================================================
    END DO
!============================================================================
    CLOSE(51)
!============================================================================
    FCXS0=0.0D0
    FCXS1=0.0D0
	DO I=1,NHS
    FCXS0=FCXS0 + DELT(I)*DELT(I)
    FCXS1=FCXS1 + DELT(I)*WD(I)*DELT(I)
    END DO
    FCXS2=FCXS0/FCXS1
	DO I=1,NHS
      WS(I)=WD(I)*FCXS2
    END DO
!==================================================      DBLE
	WC=0.0D0
    PLNHS=0.0D0
	K=0
    DO I=1,NKS
      IF (DATE_XZ(1)==1) THEN
	  K=K+1
      WC1=DLOG10( EH_PL(I)%EX )/DLOG10( EH_PJ(I)%EX )
	  WC1=1.0D0 - WC1
	  WC=WC + WC1*WC1*WS(K)*WS(K)
	  PLNHS=PLNHS + 1.0D0
	  END IF
      IF (DATE_XZ(2)==1) THEN
	  K=K+1
      WC1=DLOG10( EH_PL(I)%EY )/DLOG10( EH_PJ(I)%EY )
	  WC1=1.0D0 - WC1
	  WC=WC + WC1*WC1*WS(K)*WS(K)
	  PLNHS=PLNHS + 1.0D0
	  END IF
      IF (DATE_XZ(3)==1) THEN
	  K=K+1
      WC1=DLOG10( EH_PL(I)%EZ )/DLOG10( EH_PJ(I)%EZ )
	  WC1=1.0D0 - WC1
	  WC=WC + WC1*WC1*WS(K)*WS(K)
	  PLNHS=PLNHS + 1.0D0
	  END IF
      IF (DATE_XZ(4)==1) THEN
	  K=K+1
      WC1=DLOG10( EH_PL(I)%HX )/DLOG10( EH_PJ(I)%HX )
	  WC1=1.0D0 - WC1
	  WC=WC + WC1*WC1*WS(K)*WS(K)
	  PLNHS=PLNHS + 1.0D0
	  END IF
      IF (DATE_XZ(5)==1) THEN
	  K=K+1
      WC1=DLOG10( EH_PL(I)%HY )/DLOG10( EH_PJ(I)%HY )
	  WC1=1.0D0 - WC1
	  WC=WC + WC1*WC1*WS(K)*WS(K)
	  PLNHS=PLNHS + 1.0D0
	  END IF
      IF (DATE_XZ(6)==1) THEN
	  K=K+1
      WC1=DLOG10( EH_PL(I)%HZ )/DLOG10( EH_PJ(I)%HZ )
	  WC1=1.0D0 - WC1
	  WC=WC + WC1*WC1*WS(K)*WS(K)
	  PLNHS=PLNHS + 1.0D0
	  END IF
	END DO
!==================================================
	WC=DSQRT(WC/PLNHS)	      !计算残差
!=====================================================================
    CALL CPU_TIME(TIME1)
!=====================================================================
    WRITE(*,"(1X,A,A,I2,A,F8.3,A,F6.2,A)")OUT_SS(1:4),'号测点,第',NN,'次反演迭代的残差为:',WC*100.0D0,"%  用时:",TIME1-TIME0,"秒"
	IF (WC<=0.001) GOTO 1001  !达到反演拟合精度，退出反演
    TIME0=TIME1
!==================================================
!   拟合残差及残差下降情况
!==================================================
    WCD=WC-WCL
    IF (WCD<0.0D0) BTAS=BTAS*0.682D0
    IF (WCD>0.0D0) BTAS=BTAS*1.318D0
    IF (BTAS>=1.0E-1) BTAS=1.0E-1
    IF (BTAS<=5.0E-4) BTAS=5.0E-4
!==================================================
!   使用拟合残差及残差下降情况来判断 层厚度 是否参与反演
!==================================================
    WCD=DABS( WCD )
	  YNL=NL
	 NRPD=0
    IF ( LH==0 .AND. (WC<0.05D0 .OR. (WCD<0.01D0 .AND. WC<WCL)) ) THEN
	  YNL=NL-1
	  NRPD=1
	  LH=1
	END IF
	WCL=WC
!************************************************************************************************
!   视电阻率对模型电阻率的偏导数
!   改变多个参数，建立分别独立的模型同时正演
!************************************************************************************************
	ALLOCATE( JAC(NHS,YNL) )
!==================================================
    IF (WC<=0.1D0) THEN
    DRP=DLOG(0.8D0) / DLOG(0.1D0)/10.0D0   !扰动量
    ELSE IF (WC>=0.8D0) THEN
    DRP=DLOG(0.8D0) / DLOG(0.8D0)/10.0D0   !扰动量
	ELSE
    DRP=DLOG(0.8D0) / DLOG(WC)/10.0D0      !扰动量
	END IF
!************************************************************************************************
    IF (NRPD==1) GOTO 22
!************************************************************************************************    
    ALLOCATE( EH_PK(NKS) )
!============================================================================
    DO 10 I=1,YNL                !改变层电阻率求偏导数
!============================================================================
     MODF(1:NL)=MODS(1:NL)
!---------------------------------------------------------------------------
	 IF (MODS(I)%R .LE. 1.D0) THEN     !电阻率参数   !以免DLOG10(1)=0的出现。
	 IF ( WC<=0.1D0 ) THEN        
	   MODF(I)%R=MODF(I)%R*( 1.0D0+DLOG10(MODF(I)%R*0.1D0)*DRP ) !修改一个参数
     ELSEIF ( WC>=0.8D0 ) THEN
	   MODF(I)%R=MODF(I)%R*( 1.0D0+DLOG10(MODF(I)%R*0.8D0)*DRP ) !修改一个参数
	 ELSE
	   MODF(I)%R=MODF(I)%R*( 1.0D0+DLOG10(MODF(I)%R*WC)*DRP )    !修改一个参数
	 END IF
	 ELSE
	   MODF(I)%R=10.0D0**( DLOG10(MODF(I)%R)*( 1.0D0+DRP ) )  !修改一个参数
	 END IF
!---------------------------------------------------------------------------
	 CALL D1_SOTEMZ(TEDM,CDBLF,LSI,POTNS,NKS,NL,FJ,MODF,EH_PK)	!第二个正演
!============================================================================
!    计算偏导数矩阵
!============================================================================
     DO J=1,NKS
	 EH_PK(J)%EX=EH_PK(J)%EX*DATSX0(J)
	 EH_PK(J)%EY=EH_PK(J)%EY*DATSX1(J)
	 EH_PK(J)%EZ=EH_PK(J)%EZ*DATSX2(J)
	 EH_PK(J)%HX=EH_PK(J)%HX*DATSX3(J)
	 EH_PK(J)%HY=EH_PK(J)%HY*DATSX4(J)
	 EH_PK(J)%HZ=EH_PK(J)%HZ*DATSX5(J)
	 END DO
!============================================================================
	 L=0
	 DO J=1,NKS
	   IF (DATE_XZ(1)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Ex)-DLOG10(EH_PL(J)%EX))/( DLOG10(MODF(I)%R)-DLOG10(MODS(I)%R) )   !幅值偏导数
	   END IF
	   IF (DATE_XZ(2)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Ey)-DLOG10(EH_PL(J)%Ey))/( DLOG10(MODF(I)%R)-DLOG10(MODS(I)%R) )   !幅值偏导数
	   END IF
	   IF (DATE_XZ(3)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Ez)-DLOG10(EH_PL(J)%Ez))/( DLOG10(MODF(I)%R)-DLOG10(MODS(I)%R) )   !幅值偏导数
	   END IF
	   IF (DATE_XZ(4)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Hx)-DLOG10(EH_PL(J)%Hx))/( DLOG10(MODF(I)%R)-DLOG10(MODS(I)%R) )   !幅值偏导数
	   END IF
	   IF (DATE_XZ(5)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Hy)-DLOG10(EH_PL(J)%Hy))/( DLOG10(MODF(I)%R)-DLOG10(MODS(I)%R) )   !幅值偏导数
	   END IF
	   IF (DATE_XZ(6)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Hz)-DLOG10(EH_PL(J)%Hz))/( DLOG10(MODF(I)%R)-DLOG10(MODS(I)%R) )   !幅值偏导数
	   END IF
	END DO
!============================================================================
10  CONTINUE
    DEALLOCATE( EH_PK )
!============================================================================
22  CONTINUE
!============================================================================
    IF (NRPD==0) GOTO 44
!============================================================================
!************************************************************************************************
    ALLOCATE( EH_PK(NKS) )
!============================================================================
	DO 20 I=1,YNL                    !改变层厚度求偏导数
!============================================================================
     MODF(1:NL)=MODS(1:NL)
!============================================================================
	 IF (MODF(I)%H .LE. 1.D0) THEN     !电阻率参数   !以免DLOG10(1)=0的出现。
	 IF ( WC<=0.1D0 ) THEN        
	   MODF(I)%H=MODF(I)%H*( 1.0D0+DLOG10(MODF(I)%H*0.1D0)*DRP ) !修改一个参数
     ELSEIF ( WC>=0.8D0 ) THEN
	   MODF(I)%H=MODF(I)%H*( 1.0D0+DLOG10(MODF(I)%H*0.8D0)*DRP ) !修改一个参数
	 ELSE
	   MODF(I)%H=MODF(I)%H*( 1.0D0+DLOG10(MODF(I)%H*WC)*DRP )    !修改一个参数
	 END IF
	 ELSE
	   MODF(I)%H=10.0D0**( DLOG10(MODF(I)%H)*( 1.0D0+DRP ) )  !修改一个参数
	 END IF
!---------------------------------------------------------------------------
	 CALL D1_SOTEMZ(TEDM,CDBLF,LSI,POTNS,NKS,NL,FJ,MODF,EH_PK)	!第二个正演
!============================================================================
!    计算偏导数矩阵
!============================================================================
     DO J=1,NKS
	 EH_PK(J)%EX=EH_PK(J)%EX*DATSX0(J)
	 EH_PK(J)%EY=EH_PK(J)%EY*DATSX1(J)
	 EH_PK(J)%EZ=EH_PK(J)%EZ*DATSX2(J)
	 EH_PK(J)%HX=EH_PK(J)%HX*DATSX3(J)
	 EH_PK(J)%HY=EH_PK(J)%HY*DATSX4(J)
	 EH_PK(J)%HZ=EH_PK(J)%HZ*DATSX5(J)
	 END DO
!============================================================================
     L=0
	 DO J=1,NKS
	   IF (DATE_XZ(1)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Ex)-DLOG10(EH_PL(J)%EX))/( DLOG10(MODF(I)%H)-DLOG10(MODS(I)%H) )   !幅值偏导数
	   END IF
	   IF (DATE_XZ(2)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Ey)-DLOG10(EH_PL(J)%Ey))/( DLOG10(MODF(I)%H)-DLOG10(MODS(I)%H) )   !幅值偏导数
	   END IF
	   IF (DATE_XZ(3)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Ez)-DLOG10(EH_PL(J)%Ez))/( DLOG10(MODF(I)%H)-DLOG10(MODS(I)%H) )   !幅值偏导数
	   END IF
	   IF (DATE_XZ(4)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Hx)-DLOG10(EH_PL(J)%Hx))/( DLOG10(MODF(I)%H)-DLOG10(MODS(I)%H) )   !幅值偏导数
	   END IF
	   IF (DATE_XZ(5)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Hy)-DLOG10(EH_PL(J)%Hy))/( DLOG10(MODF(I)%H)-DLOG10(MODS(I)%H) )   !幅值偏导数
	   END IF
	   IF (DATE_XZ(6)==1) THEN
	   L=L+1
	   JAC(L,I)=(DLOG10(EH_PK(J)%Hz)-DLOG10(EH_PL(J)%Hz))/( DLOG10(MODF(I)%H)-DLOG10(MODS(I)%H) )   !幅值偏导数
	   END IF
	 END DO
!============================================================================
20   CONTINUE
     DEALLOCATE( EH_PK )
!************************************************************************************************
44   CONTINUE
     DEALLOCATE( EH_PL )
!==================================================================
!=================电阻率的圆滑约束正则化最小二乘反演===============
!=============================================Wd*A WdPs============
      ALLOCATE( ATA( NHS+YNL ,YNL+1) )
!==================================================================
      DO 657 I=1,NHS
!==================================================================
      DO 653 IL=1,YNL
         ATA(I,IL )=JAC(I,IL )*WS(I)  !误差限制
!--------------------------------------------------------
653   CONTINUE
       ATA(I,YNL+1)= DelT(I)*WS(I)
657   CONTINUE
!=========================================================================
      DEALLOCATE( JAC )
!=========================================================================
      MN=NHS
      IHDS=1
      IHDF=YNL
      SQSM=0.0D0
!****************************************************************************************
!****************************************************************************************
      DO 790 I=IHDS,IHDF
!================================================
      DO 780 J=1,YNL+1
         ATA(MN+I,J)=0.0D0    !光滑约束矩阵
780   CONTINUE
!================================================
	  WC1=0.0D0
!================================================
	  IF (MOD_RLZ==1) THEN
        IF (I/=IHDS) ATA(MN+I,I-1)= -MOD_YYS   !光滑约束矩阵
        IF (I/=IHDS) ATA(MN+I,I  )=  MOD_YYS   !光滑约束矩阵
		!---------------------------------------------------------------------
	    IF (NRPD==0) THEN
        IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I-1)*( DLOG10( MODS(I-1)%R ) )
        IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I  )*( DLOG10( MODS(I  )%R ) )
		ELSE
        IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I-1)*( DLOG10( MODS(I-1)%H ) )
        IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I  )*( DLOG10( MODS(I  )%H ) )
		END IF
		!---------------------------------------------------------------------
	  ELSE IF (MOD_RLZ==2) THEN
        IF (I/=IHDS) ATA(MN+I,I-1)=  MOD_YYS   !光滑约束矩阵
        IF (I/=IHDF) ATA(MN+I,I+1)=  MOD_YYS   !光滑约束矩阵
		IF (I/=IHDS .AND. I/=IHDF)  &
                     ATA(MN+I,I  )= -(ATA(MN+I,I-1)+ATA(MN+I,I+1)) !光滑约束矩阵
		!---------------------------------------------------------------------
	    IF (NRPD==0) THEN
        IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I-1)*( DLOG10( MODS(I-1)%R ) )
        IF (I/=IHDF) WC1=WC1 + ATA(MN+I,I+1)*( DLOG10( MODS(I+1)%R ) )
		             WC1=WC1 + ATA(MN+I,I  )*( DLOG10( MODS(I  )%R ) )
		ELSE
        IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I-1)*( DLOG10( MODS(I-1)%H ) )
        IF (I/=IHDF) WC1=WC1 + ATA(MN+I,I+1)*( DLOG10( MODS(I+1)%H ) )
		             WC1=WC1 + ATA(MN+I,I  )*( DLOG10( MODS(I  )%H ) )
		END IF
		!---------------------------------------------------------------------
	  ELSE IF (MOD_RLZ==3) THEN
		!---------------------------------------------------------------------
        IF (I/=IHDS) ATA(MN+I,I-1)= -MOD_YYS/DSQRT( ( DLOG10( MODS(I-1)%R ) )**2 + BTAS**2 ) !光滑约束矩阵
        IF (I/=IHDS) ATA(MN+I,I  )=  MOD_YYS/DSQRT( ( DLOG10( MODS(I  )%R ) )**2 + BTAS**2 ) !光滑约束矩阵
		!---------------------------------------------------------------------
	    IF (NRPD==0) THEN
		IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I-1)*( DLOG10( MODS(I-1)%R ) )
        IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I  )*( DLOG10( MODS(I  )%R ) )
		ELSE
		IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I-1)*( DLOG10( MODS(I-1)%H ) )
        IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I  )*( DLOG10( MODS(I  )%H ) )
		END IF
		!---------------------------------------------------------------------
	  ELSE IF (MOD_RLZ==4) THEN
        IF (I/=IHDS) ATA(MN+I,I-1)= MOD_YYS  !光滑约束矩阵
        IF (I/=IHDF) ATA(MN+I,I+1)= MOD_YYS  !光滑约束矩阵
		!---------------------------------------------------------------------
		IF (I/=IHDS .AND. I/=IHDF)  &
                     ATA(MN+I,I  )=-(ATA(MN+I,I-1)+ATA(MN+I,I+1))/DSQRT( ( DLOG10( MODS(I)%R ))**2 + BTAS**2 )
		IF (I/=IHDS) ATA(MN+I,I-1)= MOD_YYS/DSQRT( (DLOG10( MODS(I-1)%R ))**2 + BTAS**2 ) !光滑约束矩阵
        IF (I/=IHDF) ATA(MN+I,I+1)= MOD_YYS/DSQRT( (DLOG10( MODS(I+1)%R ))**2 + BTAS**2 ) !光滑约束矩阵
		!---------------------------------------------------------------------
	    IF (NRPD==0) THEN
        IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I-1)*( DLOG10( MODS(I-1)%R ) )
        IF (I/=IHDF) WC1=WC1 + ATA(MN+I,I+1)*( DLOG10( MODS(I+1)%R ) )
		             WC1=WC1 + ATA(MN+I,I  )*( DLOG10( MODS(I  )%R ) )
		ELSE
        IF (I/=IHDS) WC1=WC1 + ATA(MN+I,I-1)*( DLOG10( MODS(I-1)%H ) )
        IF (I/=IHDF) WC1=WC1 + ATA(MN+I,I+1)*( DLOG10( MODS(I+1)%H ) )
		             WC1=WC1 + ATA(MN+I,I  )*( DLOG10( MODS(I  )%H ) )
		END IF
		!---------------------------------------------------------------------
	  END IF
!=====================================
      SQSM=SQSM + WC1*WC1
      ATA(MN+I,YNL+1)=0.0D0
790   CONTINUE
!****************************************************************************************
      FCXS0=DSQRT( FCXS0/DBLE(NHS) )
       SQSM=DSQRT( SQSM/DBLE(YNL) )
       SQSM=DSQRT( FCXS0/(FCXS0+SQSM) )
!================================================================
      DO 990 I=1,YNL
	  DO 980 J=1,YNL+1
	  ATA(MN+I,J)=ATA(MN+I,J)*SQSM*MOD_ZYS
980   CONTINUE
990   CONTINUE
      MN=MN + YNL
!=======================================================
!=================正则化解奇异方程算法 =================
!=======================================================
      ALLOCATE( SR(YNL) )
	  CALL GRAMS(ATA,MN,YNL,SR)
      DEALLOCATE( ATA )
!=======================================================
	  IF (NRPD==0) THEN
	  DO I=1,YNL
	    SR(I)=SR(I) + DLOG10(MODS(I)%R)
	    IF ( DABS(10.0D0**( SR(i) )- MODS(I)%R ) .LT. MODS(I)%R/2.0D0) THEN
	     MODS(I)%R=10.0D0**SR(i)
	    ELSE
	     MODS(I)%R=MODS(I)%R + SIGN( MODS(I)%R/2.0D0 , 10.0D0**( SR(i) )-MODS(I)%R )
	    END IF
	  END DO
	  ELSE
	  DO I=1,YNL
	    SR(I)=SR(I) + DLOG10(MODS(I)%H)
	    IF ( DABS(10.0D0**( SR(i) )- MODS(I)%H ) .LT. MODS(I)%H/2.0D0) THEN
	      MODS(I)%H=10.0D0**SR(i)
	    ELSE
	      MODS(I)%H=MODS(I)%H + SIGN( MODS(I)%H/2.0D0 , 10.0D0**( SR(i) )-MODS(I)%H )
	    END IF
	  END DO
      END IF
!=========================================================================
    DEALLOCATE( SR )
!=========================================================================
    I=INT(NN/10)
	J=NN-I*10
	OUT_DS(1:2)="第"
	OUT_DS(3:3)=CHAR(48+I)
	OUT_DS(4:4)=CHAR(48+J)
	OUT_DS(5:18)="次反演模型.TXT"
	open(Unit=51,File="反演过程/"//OUT_DS(1:18))
	WRITE(51,"(3X,A,3X,A)")"深度","电阻率"
     SmH=MODS(1)%H*0.1D0
	 WRITE(51,'(F10.2,F12.3)')SmH,MODS(1)%R
     SmH=MODS(1)%H
	 WRITE(51,'(F10.2,F12.3)')SmH,MODS(1)%R
	 WRITE(51,'(F10.2,F12.3)')SmH,MODS(2)%R
	DO I=2,NL-1
     SmH=SmH + MODS(I)%H
	 WRITE(51,'(F10.2,F12.3)')SmH,MODS(I  )%R
	 WRITE(51,'(F10.2,F12.3)')SmH,MODS(I+1)%R
	END DO
     SmH=SmH + MODS(NL-1)%H*0.5D0
	 WRITE(51,'(F10.2,F12.3)')SmH,MODS(NL)%R
	CLOSE(51)
!======================================================
	 NN=NN+1	   !反演次数加1
!======================================================
2000 CONTINUE
     GOTO 1002
!==========================================================================================
!******************************************************************************************
!==========================================================================================
!    反演达到最大迭代次数时---推出反演---释放动态数组
!=====================================================================
1001  DEALLOCATE( EH_PL )
1002  DEALLOCATE( DELT,EH_PJ )
      DEALLOCATE( WD,WS )
      DEALLOCATE( DATSX0,DATSX1,DATSX2,DATSX3,DATSX4,DATSX5 )
!=====================================================================
!	  反演如果预先达到要求计算的精度时，则会直接转到 1001 位置
!=====================================================================
      RETURN
!=====================================================================
      END SUBROUTINE


!====================================================================
!=======================正则化解奇异方程算法=========================
!====================================================================
      SUBROUTINE GRAMS(A,N,M,R)
	  IMPLICIT NONE
	  INTEGER :: I,J,K,N,M
      REAL(8) :: A(N,M+1),RW(M,M),R(M)
	  REAL(8) :: D,RKI,X
!==========================================
      DO 5 K=1,M
      D=0.0D0
!==========================================
      DO 1 J=1,N
1       D=D+A(J,K)**2  !正则化后的主对角元素
!==========================================
!==========================================
        RW(K,K)=D      !正则化后的主对角元素
      DO 4 I=K+1,M+1
        RKI=0.0D-14
      DO 2 J=1,N
2       RKI=RKI+A(J,K)*A(J,I)
        RKI=RKI/D
        IF (I/=M+1) RW(K,I)=RKI !第K行，第I列
        IF (I/=M+1) RW(I,K)=RKI !第K行，第I列
        IF (I==M+1) R(K)=RKI    !右边向量
      DO 3 J=1,N
3       A(J,I)=A(J,I)-RKI*A(J,K)
4     CONTINUE
!==========================================
5     CONTINUE
!==================================================
      DO 7 K=2,M
         J=M-K+1
         X=R(J)
      DO 6 I=J+1,M
6       X=X-RW(J,I)*R(I)
        R(J)=X
7     CONTINUE
!==================================================
      RETURN                                                         
      END SUBROUTINE
!==================================================

!=========================================================================
!    小波分析
!=========================================================================
       SUBROUTINE WAVE_M( NM,WN,NK,RS_YS)
       IMPLICIT NONE
       INTEGER :: I,J,K,LK,LP,N,NM,WN,NK,IX,IY,IJ,IK
       REAL(8) :: RS_YS(NM),RS_UN(NM)
       REAL(8) :: SX,SY,W,U,X,PI,WH,WV,WD
!===============================================
       PI=3.141592653589793E+00
!===============================================
       DO I=1,WN
       SX=2.0D0**(I-1)  !尺度
!===============================================
       N=2
20     W=0.0D0
       DO K=-N,N
       X=DBLE(K)
       W=W + DEXP(-X*X/2.0D0/SX)/DSQRT(2.0*PI*SX)
       END DO
       IF ((1.0-W)>1.0E-6) THEN 
       N=N+1
       GOTO 20
       ELSEIF ((1.0-W)<-1.0E-6) THEN 
       N=N-1
       IF (N<=1) N=1
       GOTO 20
       END IF
!===============================================
       DO IX=1,NK
       U=0.0
       WV=0.0
       DO J=-N,N
        LP=IX+J
      	IF (LP<(1 )) LP=1
      	IF (LP>(NK)) LP=NK
      	X=DBLE(J)
	    W=DEXP(-X*X/2.0D0/SX)
	    WH=W/DSQRT(2.0D0*PI*SX)
	    WV=WV+WH
	    IF (I==1) THEN
	    U=U + WH*RS_YS(LP)
      	ELSE
      	U=U + WH*RS_UN(LP)
	    END IF
       END DO
       U=U/WV
       RS_UN(IX)=U  !I阶逼近
       END DO
!===============================================
       END DO
       RS_YS(1:NK)=RS_UN(1:NK)
!===============================================
       RETURN
!===============================================
       END SUBROUTINE
!===============================================

	END MODULE D1_CSAMTINV

