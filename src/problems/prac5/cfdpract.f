	PROGRAM CFD_5

C Version November 2001
C
C October 2003: convergence criterion Poisson equation changed
C
C***************************************************************************
C
c  --- Most important variables ---
c  NF(I,J)    : cell label - 0=Ful 7=Out 8=In 9=Bnd/Obs
c  P/PN(I,J)  : pressure at new/old time level
c  U/UN(I,J)  : radial velocity at new/old time level
c  V/VN(I,J)  : axial velocity at new/old time level

c  --- Subroutines ---
c  BCBND : sets boundary conditions for velocity components.
c  CFLCHK: monitors CFL number and sets flag for time-step adjustment
c  COEFF : defines coefficient matrix for Poisson equation
c          including boundary conditions at wall and obstacles.
c  DTADJ : halves/doubles time step (old time step will be repeated).
c  GRID  : makes (non-uniform) grid.
c  INIT  : starts new time step.
c  MILU  : solves Poisson equation.
c  PRT   : prints and writes results.
c  SETFLD: initializes fluid configuration.
c  SETPAR: reads input file.
c  SOLVEP: organizes pressure calculation and updates velocity.
c  SOR   : alternative Poisson solver by SOR. 
c  TILDE : integrates momentum equations.
c  PRTCNF, etc: generate plots and plotfiles.
c
c  --- Main calling sequence, excluding postprocessing routines ---
c  initialisation  SETPAR  GRID
c                  SETFLD 
c  time step       INIT    BCBND
c                  TILDE   
c                  SOLVEP  COEFF
c                          MILU
c                          SOR
c                          BCBND
c                  CFLCHK  DTADJ
c***********************************************************************

      IMPLICIT NONE

      INTEGER NSMAX,SDTMAX
      PARAMETER (NSMAX=20, SDTMAX=10)

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

      common /rel/ relatief
      real relatief

      COMMON /PRINTS/ TPrt,PrtDT,TAdPrt,TCnf,NP
     >   ,PGnu,PAvs,Puvpf,PVelop,PForce
      REAL TPrt, PrtDT, TAdPrt, TCnf(2001)
      INTEGER NP,PGnu,PAvs,Puvpf,PVelop,PForce

      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
     >               StrtP
      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
      REAL  Epsi,OmStrt,Alpha,NrmRhs

c  local variables
      REAL HalfDT
      INTEGER FlgAdj

c**   most units opened and closed in PRT(1) and PRT(3), resp.

      OPEN (UNIT=5, FILE='cfd.in')

C**** PROGRAM SETUP ***************************************************

      CALL SETPAR
      CALL SETFLD
      CALL DIVERR
C
C**** PLOT/RESTART CLOCK
C
      TStart = T
      TPrt = T
      TAdPrt = T

      CALL PRT(1)

C**** START TIME CYCLE ************************************************
      FlgAdj=0
   10 CONTINUE

c--- initialization of new time step
         T = T + DelT
         Cycle = Cycle + 1
         HalfDT = 0.5*DelT
         CALL INIT
         IF (FlgAdj.NE.0) relatief = 1.0         

c--- momentum equation
         CALL TILDE

c--- Poisson equation
         CALL SOLVEP

c--- check CFL criterion and adjust time step
         CALL CFLCHK(FlgAdj)


c--- finish time step when OK
         IF (FlgAdj.EQ.0) THEN         
            CALL DIVERR
c...    print and save results
            IF (T+HalfDT .GT. TPrt) CALL PRT(2)
         ENDIF

c--- check for end of computation
         IF (CYCLE.LT.4) relatief = 1.0
         IF (T+HalfDT .LT. TFin .AND.
     >       relatief.GE. -1e-6*delT) GOTO 10         


C**** END OF TIME CYCLE ***********************************************

C--- stop program
 20   CALL PRT(3)
      CLOSE (5)
      STOP

      END
c  End of main program

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE BCBND
      IMPLICIT NONE

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX, RDY,
     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

      COMMON /ORGA/ NF
      INTEGER NF(IMAX,JMAX)

      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
     1     , POutN, POutE, POutS, POutW
      INTEGER WB, WL, WR, WT

c  local variables
      INTEGER i,j

c... velocity conditions: SgnL=1 slip wall; SgnL=-1 no-slip wall

      DO 100 J=1,JMaxUs
c-- western boundary
         V(1,J)=SgnL*V(2,J)
         IF (nf(1,j).EQ.9) THEN
            U(1,J)=0.0
            P(1,J)=P(2,J)
         ELSEIF (nf(1,j).EQ.8) THEN
            u(1,j)=UIn*COS(FreqIn*t)
            u(0,j)=u(1,j)
            v(1,j)=2.*VIn*COS(FreqIn*t) - v(2,j)
            P(1,J)=P(2,J)
         ELSEIF (nf(1,j).EQ.7) THEN
            IF(nf(2,j).NE.0) u(1,j)=u(2,j)
            u(0,j)=2.0*u(1,j)-u(2,j)
c... negative outflow
            IF (u(1,j).GT.0.0) u(0,j)=u(1,j)
            P(1,J)=POutW 
         ENDIF
C
c-- eastern boundary
         V(IMaxUs,J)=SgnR*V(IM1Us,J)
         IF (nf(IMaxUs,j).EQ.9) THEN
            U(IM1Us,J)=0.0
            P(IMaxUs,J)=P(IM1Us,J)
         ELSEIF (nf(IMaxUs,j).EQ.8) THEN
            u(IM1Us,j)=UIn*COS(FreqIn*t)
            u(IMaxUs,j)=u(IM1Us,j)
            v(IMaxUs,j)=2.0*VIn*COS(FreqIn*t)-v(IM1Us,j)
            P(IMaxUs,J)=P(IM1Us,J)
         ELSEIF (nf(IMaxUs,j).EQ.7) THEN
            IF (nf(IM1Us,j).NE.0) u(IM1Us,j)=u(IM2Us,j)
            u(IMaxUs,j)=2.0*u(IM1Us,j)-u(IM2Us,j)
c... negative outflow
            IF (u(IM1Us,j).LT.0.0) u(IMaxUs,j)=u(IM1Us,j)
            P(IMaxUS,J)=POutE 
C            P(IMaxUS,J)=P(IM1Us,J)
         ENDIF
C
 100  CONTINUE
C
      DO 200 I=1,IMaxUs
c-- southern boundary
         U(I,1)=SgnB*U(I,2)
         IF (nf(i,1).EQ.9) THEN
            V(I,1)=0.0
            P(I,1)=P(I,2)
         ELSEIF (nf(i,1).EQ.8) THEN
            U(I,1)=2.0*UIn*COS(FreqIn*t)-U(I,2)
            v(i,1)=VIn*COS(FreqIn*t)
            v(i,0)=v(i,1)
            P(I,1)=P(I,2)
         ELSEIF (nf(i,1).EQ.7) THEN
            IF (nf(i,2).NE.0) v(i,1)=v(i,2)
            v(i,0)=2.0*v(i,1)-v(i,2)
c... negative outflow
            IF (v(i,1).GT.0.0) v(i,0)=v(i,1)
            P(I,1)=POutS
         ENDIF
C
c-- northern boundary
         U(I,JMaxUs)=SgnT*U(I,JM1Us)
         IF (nf(i,JMaxUs).EQ.9) THEN
            V(I,JM1Us)=0.0
            P(I,JMaxUs)=P(I,JM1Us)
         ELSEIF (nf(i,JMaxUs).EQ.8) THEN
            U(I,JMaxUs)=2.0*UIn*COS(FreqIn*t)-U(I,JM1Us)
            v(i,JM1Us)=VIn*COS(FreqIn*t)
            v(i,JMaxUs)=v(i,JM1Us)
            P(I,JMaxUs)=P(I,JM1Us)
         ELSEIF (nf(i,JMaxUs).EQ.7) THEN
            IF (nf(i,JM1Us).NE.0) v(i,JM1Us)=v(i,JM2Us)
            v(i,JMaxUs)=2.0*v(i,JM1Us)-v(i,JM2Us)
c... negative outflow
            IF (v(i,JM1Us).LT.0.0) v(i,JMaxUs)=v(i,JM1Us)
            P(I,JMaxUs)=POutW
         ENDIF
C
  200 CONTINUE

      RETURN
      END
c End of subroutine BCBND
c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      SUBROUTINE CFLCHK(FlgAdj)
      IMPLICIT NONE

      INTEGER FlgAdj

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX, RDY,
     1   CX,CY, IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

      COMMON /ADJUST/ CFLMin, CFLMax, NDTAd, MxDTAd, SlwFlw
      REAL CFLMin, CFLMax
      INTEGER NDTAd, MxDTAd, SlwFlw

c-- local variables
      INTEGER i,j
      REAL kladx,klady,CFL,CFLx,CFLy,CFLtmp

      FlgAdj = 0
      CFL=0.0
      DO 20 J=2,JM1Us
         klady=delt/dely(j)
      DO 20 I=2,IM1Us
         kladx=delt/delx(i)
         CFLx=ABS(u(i,j))*kladx
         CFLy=ABS(v(i,j))*klady
         CFLtmp=AMAX1(CFLx,CFLy)
         CFL=AMAX1(CFLtmp,CFL)
 20   CONTINUE

c--- check for large time step

      IF (CFL.GT.5.0) THEN
         WRITE(6,*) 'CFL is larger than 5'
         STOP 'CFL is larger than 5'
      ENDIF

      IF (CFL.GT.CFLMax) THEN
c--- reduce time step
         FlgAdj=1
         IF (NDTAd.GE.ABS(MxDTAd)) THEN
            FlgAdj=0
         ENDIF
      ENDIF

c--- check for small time step;
c    if 25 times in a row  CFL < CFLMin  then enlarge time step

      IF (CFL.LT.CFLMin) THEN
         SlwFlw=SlwFlw+1
      ELSE
         SlwFlw=0
      ENDIF

      IF (SlwFlw.EQ.25 .AND. DelT.LT.0.5*DelTMx) THEN
         SlwFlw=0
         FlgAdj=-1
      ENDIF


c--- adjust time step if required

      IF (FlgAdj .NE. 0) CALL DTADJ(FlgAdj)

      RETURN
      END
c  End of subroutine CFLCHK

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DTADJ(FlgAdj)
      IMPLICIT NONE

      INTEGER FlgAdj

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX, RDY,
     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

      COMMON /ADJUST/ CFLMin, CFLMax, NDTAd, MxDTAd, SlwFlw
      REAL CFLMin, CFLMax
      INTEGER NDTAd, MxDTAd, SlwFlw

c local variables
      INTEGER i,j

c    when MxDTAd > 0 then DelT can be reduced and enlarged
c    when MxDTAd < 0 then DelT may only be reduced

c--- NDTAd is effective number of time-step reductions
c    FlgAdj = 1 reduction of DelT; FlgAdj = -1 enlargement

      NDTAd = NDTAd + FlgAdj

c-- reset time level

      t = t - delt

c-- time step is adjusted

      IF (FlgAdj.EQ.1) THEN
         delt = 0.5*delt
         WRITE(6,4000) t, delt
      ELSEIF (FlgAdj.EQ.-1 .AND. MxDTAd.GT.0) THEN
         delt = 2.0*delt
         WRITE(6,4001) t, delt
      ENDIF

c-- reset all variables

      DO 10 J=1,JMaxUs
         U(0,J)=UN(0,J)
      DO 10 I=1,IMaxUs
         P(I,J)=PN(I,J)
         U(I,J)=UN(I,J)
         V(I,J)=VN(I,J)
 10   CONTINUE
      DO 11 I=1,IMaxUs
         V(I,0)=VN(I,0)
 11   CONTINUE

 4000 FORMAT (/,2X,' T =',E12.5,'  TIME STEP REDUCED: DelT =',E11.4,/)
 4001 FORMAT (/,2X,' T =',E12.5,'  TIME STEP DOUBLED: DelT =',E11.4,/)

      RETURN
      END
c  End of subroutine DTADJ

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GRID
      IMPLICIT NONE

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX,RDY,
     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /DXRL/ DXRI,DXLI,DYTJ,DYBJ
      REAL DXRI(IMAX),DXLI(IMAX),DYTJ(JMAX),DYBJ(JMAX)

      COMMON /coef2/ CWI,CEI,CNJ,CSJ
      REAL CWI(IMAX),CEI(IMAX),CNJ(JMAX),CSJ(JMAX)

      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI
      REAL EMF,EMF1,EM6,EM10,EP10,PI

      COMMON /TANK/ XMax, XMin, YMax, YMin
      REAL  XMax, XMin, YMax, YMin

      COMMON /REK/ xpos,ypos
      REAL xpos,ypos

c  local variables
      REAL stepx,stepy,scalex,scaley,ymid,xmid,ccx,ccy
     >      ,reli,kladl,kladr,AbsCX,AbsCY
      INTEGER i,j,jmid,imid
C
C**********     GRID IN X-DIRECTION     ***************************
C
C... definition of stretched grid
C    when cx = 0 uniform grid
C         cx > 0 refinement near position x=Xpos
C         cx < 0 refinement near both endpoints

         xpos=AMAX1(xpos,xmin)
         xpos=AMIN1(xpos,xmax)
         AbsCX=ABS(CX)

      IF (AbsCX.GE.EM6) THEN
c.. determine number of gridpoints per interval (logarithmic split)
         reli=ALOG(AbsCX*(xpos-xmin)+1.)/
     >     (ALOG(AbsCX*(xpos-xmin)+1.)+ALOG(AbsCX*(xmax-xpos)+1.))
         IMID=NINT(reli*float(IM2Us))+1
         XMID=xpos
c.. stretching factors for both intervals
         kladl=1.0+2.0*(reli-0.5)/(0.5*AbsCX**2+1.0)
         kladr=1.0+2.0*(0.5-reli)/(0.5*AbsCX**2+1.0)
      ENDIF

      IF (AbsCX.LT.EM6) THEN
         STEPX=(XMax-XMin)/FLOAT(IM2Us)
         DO 2 I=1,IM1Us
    2    X(I)=XMin+(I-1)*STEPX

      ELSEIF (CX.GT.0.0) THEN

c.. refinement near X=Xpos
c.. left interval
         CCX=AbsCX*kladl
         SCALEX=1.0
         IF (CCX.GT.em6) SCALEX=(XMID-XMin)/TANH(CCX)
         DO 4 I=1,IMID-1
 4          X(I)=XMin+SCALEX*TANH(CCX*FLOAT(I-1)/FLOAT(IMID-1))
c.. `mid'point
         X(IMID)=XMID
c.. right interval
         CCX=AbsCX*kladr
         SCALEX=1.0
         IF (CCX.GT.em6) SCALEX=(XMax-XMID)/TANH(CCX)
         DO 5 I=IMID+1,IM1Us
 5          X(I)=XMax-SCALEX*TANH(CCX*FLOAT(IM1Us-I)/FLOAT(IM1Us-IMID))

      ELSEIF (CX.LT.0.0) THEN

c.. refinement near both endpoints, away from X=Xpos
c         IMID=INT(0.5*float(IM2Us))+1
c         XMID=(xmax+xmin)*0.5
c.. left interval
         CCX=AbsCX
         SCALEX=(XMID-XMin)/TANH(CCX)
         DO 7 I=1,IMID-1
 7          X(I)=XMID-SCALEX*TANH(CCX*FLOAT(IMID-I)/FLOAT(IMID-1))
c.. midpoint
         X(IMID)=XMID
c.. right interval
         CCX=AbsCX
         SCALEX=(XMax-XMID)/TANH(CCX)
         DO 8 I=IMID+1,IM1Us
 8          X(I)=XMID+SCALEX*TANH(CCX*FLOAT(I-IMID)/FLOAT(IM1Us-IMID))

      ENDIF
C
C... auxiliary quantities
C
      DO 10 I=2,IM1Us
         DelX(I)=X(I)-X(I-1)
         RDX(I)=1./DelX(I)
         XI(I)=0.5*(X(I-1)+X(I))
 10   CONTINUE
C
      DelX(1)=DelX(2)
      RDX(1)=1.0/DelX(1)
      XI(1)=XI(2)-DelX(2)
      DelX(IMaxUs)=DelX(IM1Us)
      RDX(IMaxUs)=1.0/DelX(IMaxUs)
      XI(IMaxUs)=XI(IM1Us)+DelX(IM1Us)
      X(IMaxUs)=X(IM1Us)+DelX(IM1Us)
C

      DO 11 I=2,IM1Us
         DXLI(I)=0.5*(DelX(I)+DelX(I-1))
         DXRI(I)=0.5*(DelX(I)+DelX(I+1))
         CWI(I)=RDX(I)/DXLI(I)
         CEI(I)=RDX(I)/DXRI(I)
 11   CONTINUE
      DXLI(1)=DXLI(2)
      DXRI(1)=DXRI(2)
C
C***********     GRID IN Y-DIRECTION      ***************************
C
C... definition of stretched grid
C    when cy = 0 uniform grid
C         cy > 0 refinement near Ypos
C         cy < 0 refinement away from Ypos
C
         ypos=AMAX1(ypos,ymin)
         ypos=AMIN1(ypos,ymax)
         AbsCY=ABS(CY)

      IF (AbsCY .GE. EM6) THEN
c.. determine number of gridpoints per interval (logarithmic split)
         reli=ALOG(AbsCY*(ypos-ymin)+1.)/
     >     (ALOG(AbsCY*(ypos-ymin)+1.)+ALOG(AbsCY*(ymax-ypos)+1.))
         JMID=NINT(reli*float(JM2Us))+1
         YMID=ypos
c.. stretching factors for both intervals
         kladl=1.0+2.0*(reli-0.5)/(0.5*AbsCY**2+1.0)
         kladr=1.0+2.0*(0.5-reli)/(0.5*AbsCY**2+1.0)
      ENDIF
C
      IF (AbsCY.LT.EM6) THEN

         STEPY=(YMax-YMin)/FLOAT(JM2Us)
         DO 12 J=1,JM1Us
   12    Y(J)=YMin + (J-1)*STEPY

      ELSEIF (CY.GT.0.0) THEN

c.. refinement near Y=Ypos
c.. lower interval
         CCY=AbsCY*kladl
         SCALEY=1.0
         IF (CCY.GT.em6) SCALEY=(YMID-YMin)/TANH(CCY)
         DO 14 J=1,JMID-1
 14         Y(J)=YMin+SCALEY*TANH(CCY*FLOAT(J-1)/FLOAT(JMID-1))
c.. `mid'point
         Y(JMID)=YMID
c.. upper interval
         CCY=AbsCY*kladr
         SCALEY=1.0
         IF (CCY.GT.em6) SCALEY=(YMax-YMID)/TANH(CCY)
         DO 15 J=JMID+1,JM1Us
 15         Y(J)=YMax-SCALEY*TANH(CCY*FLOAT(JM1Us-J)/FLOAT(JM1Us-JMID))

      ELSEIF (CY.LT.0.0) THEN

c.. refinement near both endpoints, away from Y=Ypos
c         JMID=INT(0.5*float(JM2Us))+1
c         YMID=(ymax+ymin)*0.5
c.. lower interval
         SCALEY=(YMID-YMin)/TANH(CY)
         DO 17 J=1,JMID-1
 17         Y(J)=YMID-SCALEY*TANH(CY*FLOAT(JMID-J)/FLOAT(JMID-1))
c.. midpoint
         Y(JMID)=YMID
c.. upper interval
         SCALEY=(YMax-YMID)/TANH(CY)
         DO 18 J=JMID+1,JM1Us
 18         Y(J)=YMID+SCALEY*TANH(CY*FLOAT(J-JMID)/FLOAT(JM1Us-JMID))

      ENDIF

C
C... auxiliary quantities
C
      DO 20 J=2,JM1Us
         DelY(J)=Y(J)-Y(J-1)
         RDY(J)=1.0/DelY(J)
         YJ(J)=0.5*(Y(J-1)+Y(J))
 20   CONTINUE
C
      DelY(1)=DelY(2)
      RDY(1)=RDY(2)
      YJ(1)=YJ(2)-DelY(2)
      DelY(JMaxUs)=DelY(JM1Us)
      RDY(JMaxUs)=1.0/DelY(JMaxUs)
      YJ(JMaxUs)=YJ(JM1Us)+DelY(JM1Us)
      Y(JMaxUs)=Y(JM1Us)+DelY(JM1Us)
C
      DO 25 J=2,JM1Us
         DYBJ(J)=0.5*(DelY(J)+DelY(J-1))
         DYTJ(J)=0.5*(DelY(J)+DelY(J+1))
         CSJ(J)=RDY(J)/DYBJ(J)
         CNJ(J)=RDY(J)/DYTJ(J)
   25 CONTINUE
      DYBJ(1)=DYBJ(2)
      DYTJ(1)=DYTJ(2)
C
      RETURN
      END
c  End of subroutine GRID

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE INIT
      IMPLICIT NONE

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX, RDY,
     1   CX,CY, IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      common /rel/ relatief
      real relatief      

      COMMON /CONDIF/ RU,RV,divuv,divmax
      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
      INTEGER divuv

c local variables
      REAL POld,verschil,absoluut
      INTEGER i,j

      verschil=0.0
      absoluut=0.0
      DO 5 J=2,JM1Us
      DO 5 I=2,IM1Us
         verschil = verschil +
     >            (UN(I,J)-U(I,J))**2 + (VN(I,J)-V(I,J))**2
         absoluut = absoluut + U(I,J)**2 + V(I,J)**2
 5    CONTINUE   
      verschil=sqrt(verschil/(IM2Us*JM2Us))
      absoluut=sqrt(absoluut/(IM2Us*JM2Us))
      relatief =  verschil / (absoluut + 1e-10)

      DO 10 J=1,JMaxUs
         UN(0,J)=U(0,J)
         U(0,J)=0.0
      DO 10 I=1,IMaxUs
         UN(I,J) = U(I,J)
         VN(I,J) = V(I,J)
         U(I,J) = 0.0
         V(I,J) = 0.0
         POLD = P(I,J)
c         P(I,J)=2.0*P(I,J)-PN(I,J)
         PN(I,J)=POLD
 10   CONTINUE
      DO 11 I=1,IMaxUs
        VN(I,0)=V(I,0)
        V(I,0)=0.0
 11   CONTINUE

c-- restore non-zero boundary conditions
      CALL BCBND

      RETURN
      END
c End of subroutine INIT

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE DIVERR

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX, RDY,
     1   CX,CY, IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /CONDIF/ RU,RV,divuv,divmax
      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
      INTEGER divuv

      INTEGER i,j
      REAL divu

      divmax=0.0
      DO 15 J=2,JM1Us
      DO 15 I=2,IM1Us
	 divu=RDX(I)*(U(I,J)-U(I-1,J)) + 
     >        RDY(J)*(V(I,J)-V(I,J-1))
	 divmax=AMAX1(ABS(divu),divmax)
 15    CONTINUE

       RETURN
       END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SETFLD
      IMPLICIT NONE

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX, RDY,
     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /ORGA/ NF
      INTEGER NF(IMAX,JMAX)

      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
     1     , POutN, POutE, POutS, POutW
      INTEGER WB, WL, WR, WT

      COMMON /TANK/ XMax, XMin, YMax, YMin
      REAL  XMax, XMin, YMax, YMin

c local variables
      INTEGER i,j
C
C**** PRESSURE and VELOCITIES
C
      DO 10 J=1,JMaxUs
      DO 10 I=1,IMaxUs
         PN(I,J)=0.0
         P(I,J)=0.0
         U(I,J)=0.0
         V(I,J)=0.0
         UN(I,J)=U(I,J)
         VN(I,J)=V(I,J)
         F(I,J)=1.0
         NF(I,J)=0
 10   CONTINUE

      DO J=1,JMaxUs
         U(0,J)=U(1,J)
         UN(0,J)=U(0,J)
      ENDDO
      DO I=1,IMaxUs
         V(I,0)=V(I,1)
         VN(I,0)=V(I,0)
      ENDDO
C
c... define boundary conditions for F
c... solid: NF=9, F=0;  inflow: NF=8, F=1;  outflow: NF=7, F=1
c
      DO 70 I=2,IM1Us
         NF(I,1)=9
         NF(I,JMaxUs)=9
 70   CONTINUE
C
       DO 71 J=1,JMaxUs
         NF(1,J)=9
         NF(IMaxUs,J)=9
 71   CONTINUE

c...  INFLOW / OUTFLOW boundaries
c       
      IF (wb.EQ.7) THEN
         DO i=1,IMaxUs
            nf(i,1)=wb
         ENDDO
      ELSEIF (wb.EQ.8) THEN
         DO i=2,IM1Us
            nf(i,1)=wb
         ENDDO
      ENDIF

      IF (wt.EQ.7) THEN
         DO i=1,IMaxUs
            nf(i,JMaxUs)=wt
         ENDDO
      ELSEIF (wt.EQ.8) THEN
         DO i=2,IM1Us
            nf(i,JMaxUs)=wt
         ENDDO
      ENDIF

      IF (wl.EQ.7) THEN
         DO j=1,JMaxUs
            nf(1,j)=wl
         ENDDO
      ELSEIF (wl.EQ.8) THEN
         DO j=2,JM1Us
            nf(1,j)=wl
         ENDDO
      ENDIF

      IF (wr.EQ.7) THEN
         DO j=1,JMaxUs
            nf(IMaxUs,j)=wr
         ENDDO
      ELSEIF (wr.EQ.8) THEN
         DO j=2,JM1Us
            nf(IMaxUs,j)=wr
         ENDDO
      ENDIF
c
c... adaption of vof-function to cell labels
c
      DO j=1,JMaxUs
      DO i=1,IMaxUs
         IF (nf(i,j).EQ.9) THEN
            f(i,j)=0.0
         ELSEIF (nf(i,j).EQ.7 .OR. nf(i,j).EQ.8) THEN
            f(i,j)=1.0
         ENDIF
      ENDDO
      ENDDO

c      call bcbnd
C
      RETURN
      END
c  End of subroutine SETFLD

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SETPAR
      IMPLICIT NONE

      INTEGER IMAX,JMAX,NSMAX,SDTMAX,NRFMAX,NRHMAX,mmpnts
      PARAMETER (IMAX=130, JMAX=130, NSMAX=20, SDTMAX=10, 
     >           NRFMAX=10, NRHMAX=10, mmpnts=25)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX, RDY,
     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
     1     , POutN, POutE, POutS, POutW
      INTEGER WB, WL, WR, WT

      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
     >               StrtP
      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
      REAL  Epsi,OmStrt,Alpha,NrmRhs

      COMMON /CONDIF/ RU,RV,divuv,divmax
      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
      INTEGER divuv

      COMMON /CONST/ EMF, EMF1, EM6, EM10, EP10, PI
      REAL EMF,EMF1,EM6,EM10,EP10,PI

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

      COMMON /ADJUST/ CFLMin, CFLMax, NDTAd, MxDTAd, SlwFlw
      REAL CFLMin, CFLMax
      INTEGER NDTAd, MxDTAd, SlwFlw

      COMMON /PRINTS/ TPrt,PrtDT,TAdPrt,TCnf,NP
     >   ,PGnu,PAvs,Puvpf,PVelop,PForce
      REAL TPrt, PrtDT, TAdPrt, TCnf(2001)
      INTEGER NP,PGnu,PAvs,Puvpf,PVelop,PForce

      COMMON /TANK/ XMax, XMin, YMax, YMin
      REAL  XMax, XMin, YMax, YMin

      COMMON /FLUID/ Nu
      REAL  Nu

      COMMON /REK/ posx,posy
      REAL posx,posy

c local variables
      REAL dx,dy,h1,PIn
      INTEGER i,j,IPIn

C**** CONSTANT PARAMETERS
      PI = 3.14159265359
      EP10 = 1.0e10
      EM6  = 1.0e-6
      EM10 = 1.0e-10
      EMF  = 1.0e-5
      EMF1 = 1.0-EMF

C**** TANK GEOMETRY
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) XMIN,XMAX,YMIN,YMAX
      READ(5,*)

C**** GRID DEFINITION
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) IMaxUs,JMaxUs,CX,CY,posx,posy
      READ(5,*)
      IF (IMaxUs.GT.IMAX .OR. JMaxUs.GT.JMAX)
     >     STOP  'gridsize too large'
      IM1Us = IMaxUs -1
      JM1Us = JMaxUs -1
      IM2Us = IMaxUs -2
      JM2Us = JMaxUs -2

      CALL GRID

C**** MATERIAL QUANTITIES
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) Nu
      READ(5,*)

C**** BOUNDARY CONDITIONS AND IN-/OUTFLOW
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) WL, WR, WT, WB, UIn, VIn
      READ(5,*)
      IF ((WL.NE.1).and.(WL.NE.2).and.(WL.NE.7).and.(WL.NE.8)) THEN
         write(6,*) 'error in cfd.in at boundary condition left'
         stop
      ENDIF
      IF ((WR.NE.1).and.(WR.NE.2).and.(WR.NE.7).and.(WR.NE.8)) THEN
         write(6,*) 'error in cfd.in at boundary condition right'
         stop
      ENDIF
      IF ((WT.NE.1).and.(WT.NE.2).and.(WT.NE.7).and.(WT.NE.8)) THEN
         write(6,*) 'error in cfd.in at boundary condition top'
         stop
      ENDIF
      IF ((WB.NE.1).and.(WB.NE.2).and.(WB.NE.7).and.(WB.NE.8)) THEN
         write(6,*) 'error in cfd.in at boundary condition bottom'
         stop
      ENDIF
      FreqIn=0
      IPIn=0
      PIn=-1.0
      POutN=0.0
      POutE=0.0
      POutS=0.0
      POutW=0.0
      IF (IPIn.EQ.1) POutW=PIn
      IF (IPIn.EQ.2) POutE=PIn
      IF (IPIn.EQ.3) POutN=PIn
      IF (IPIn.EQ.4) POutS=PIn
c
      SgnL=1.0
      IF (WL.EQ.2) SgnL=-1.0
      SgnR=1.0
      IF (WR.EQ.2) SgnR=-1.0
      SgnT=1.0
      IF (wt.EQ.2) SgnT=-1.0
      SgnB=1.0
      IF (wb.EQ.2) SgnB=-1.0

C**** NUMERICAL MODEL PARAMETERS
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) Alpha, Divuv
      READ(5,*)

      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) Epsi, ItMax, StrtP
      ItSum = 0
      OmStrt=1.0
      IMilu=1
      READ(5,*)

C**** TIME PARAMETERS
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) TFin,DelT,CFLMax
      READ(5,*)
      T=0.0
C--- timestep DelT adjustment to stability limits
c    look for smallest mesh size
      DX = XMax-XMin
      DO I=1,IM1Us
         IF (DelX(I).LT.DX) DX=DelX(I)
      ENDDO
      DY = YMax-YMin
      DO J=1,JM1Us
         IF (DelY(J).LT.DY) DY=DelY(J)
      ENDDO
      H1 = 1./(DX*DX) + 1./(DY*DY)
      DelTMx = .5/(H1*(Nu+em10))
      DelT = AMIN1(0.9*DelTMx, DelT)

C**** TIME-STEP ADJUSTMENT CONTROL
      CFLMin=0.4*CFLMax
      MxDTAd=5
      NDTAd=0
      SlwFlw=0
C
C**** PRINT/PLOT CONTROL
      READ(5,*)
      READ(5,*)
      READ(5,*,ERR=9999) PrtDT,Puvpf
      READ(5,*)
C
C**** ORGANISATIONAL PARAMETERS
      NP = 0
      Cycle = 0
      VMAX=0.0
      pmin=0.0
      pmax=0.0

      RETURN

 9999 STOP 'error in reading CFD.IN'

      END
c  End of subroutine SETPAR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TILDE
      IMPLICIT NONE

C NEW version -- conservative symmetry preserving (November 1998)
C Only plane geometry -- axisymmetry has still to be implemented!!!!!!!!

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX,RDY,
     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /ORGA/ NF
      INTEGER NF(IMAX,JMAX)

      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
     >               StrtP
      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
      REAL  Epsi,OmStrt,Alpha,NrmRhs

      COMMON /FLUID/ Nu
      REAL  Nu

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

      COMMON /CONDIF/ RU,RV,divuv,divmax
      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
      INTEGER divuv

c local variables
      REAL dxc,dyc,rdxc,rdxe,rdxw,rdyc,rdyn,rdys
      REAL hue,huw,hvn,hvs,csom,ce,cw,cn,cs,de,dw,dn,ds,klad
      INTEGER i,j

c      double precision dsecnd,t1,t2
C
C**** CALCULATION OF TEMPORARY VELOCITIES
C

C....       TILDE IN X-DIRECTION
C
      DO 10 J=2,JM1Us
         RDYN=1./(YJ(J+1)-YJ(J))
         RDYS=1./(YJ(J)-YJ(J-1))
      DO 10 I=1,IM1Us

c--- zeroise all velocities but do not change IN and BND cells
         IF (nf(i,j).LE.7 .AND. nf(i+1,j).LE.7) u(i,j)=0.0

         IF ( (NF(I,J).EQ.0 .OR. nf(i,j).EQ.7) .AND.
     >           (NF(I+1,J).EQ.0 .OR. nf(i+1,j).EQ.7) ) THEN
c.. both i,j and i+1,j FUL or OUT cell

            DXC=XI(I+1)-XI(I)
            RDXC=1./DXC
            RDXE=RDX(I+1)
            RDXW=RDX(I)

C
c.. nieuwe formulering
            HUE=0.25*(UN(I,J)+UN(I+1,J))
            HUW=0.25*(UN(I,J)+UN(I-1,J))
            HVN=0.25*(DelX(I)*VN(I,J)+DelX(I+1)*VN(I+1,J))*RDXC
            HVS=0.25*(DelX(I)*VN(I,J-1)+DelX(I+1)*VN(I+1,J-1))*RDXC

c coefficients to be scaled by DXC*DELY(J) for (skew)symmetry

            DE=Nu*RDXE+ALPHA*ABS(HUE)
            CE=( HUE - DE) * RDXC
            DW=Nu*RDXW+ALPHA*ABS(HUW)
            CW=(-HUW - DW) * RDXC
            DN=Nu*RDYN+ALPHA*ABS(HVN)
            CN=( HVN - DN) * RDY(J)
            DS=Nu*RDYS+ALPHA*ABS(HVS)
            CS=(-HVS - DS) * RDY(J)
            CSom=-(DE+DW)*RDXC-(DN+DS)*RDY(J)

c expliciet
           klad=-CE*UN(I+1,J)-CW*UN(I-1,J)-
     >            CN*UN(I,J+1)-CS*UN(I,J-1)+Csom*UN(I,J)
           RU(I,J)=KLAD
           U(I,J)=UN(I,J)+DelT*KLAD

         ENDIF
 10   CONTINUE

C
C....       TILDE IN Y-DIRECTION
C
      DO 20 J=1,JM1Us
         DYC=YJ(J+1)-YJ(J)
         RDYC=1./DYC
         RDYN=RDY(J+1)
         RDYS=RDY(J)
      DO 20 I=2,IM1Us

c--- zeroise all velocities but do not change IN and BND cells
         IF (nf(i,j).LE.7 .AND. nf(i,j+1).LE.7) v(i,j)=0.0

         IF ( (NF(I,J).EQ.0 .OR. nf(i,j).EQ.7) .AND.
     >           (NF(I,J+1).EQ.0 .OR. nf(i,j+1).EQ.7) ) THEN
c.. both i,j and i,j+1 FUL or OUT cell

            RDXE=1./(XI(I+1)-XI(I))
            RDXW=1./(XI(I)-XI(I-1))

            HUE=0.25*(DelY(J)*UN(I,J)+DelY(J+1)*UN(I,J+1))*RDYC
            HUW=0.25*(DelY(J)*UN(I-1,J)+DelY(J+1)*UN(I-1,J+1))*RDYC
            HVN=0.25*(VN(I,J)+VN(I,J+1))
            HVS=0.25*(VN(I,J)+VN(I,J-1))

c coefficients to be scaled by DELX(I)*DYC for (skew)symmetry
            DE=Nu*RDXE+ALPHA*ABS(HUE)
            CE=( HUE - DE) * RDX(I)
            DW=Nu*RDXW+ALPHA*ABS(HUW)
            CW=(-HUW - DW) * RDX(I)
            DN=Nu*RDYN+ALPHA*ABS(HVN)
            CN=( HVN - DN) * RDYC
            DS=Nu*RDYS+ALPHA*ABS(HVS)
            CS=(-HVS - DS) * RDYC
            CSom=-(DE+DW)*RDX(I)-(DN+DS)*RDYC

c expliciet
            KLAD=-CE*VN(I+1,J)-CW*VN(I-1,J)-
     >            CN*VN(I,J+1)-CS*VN(I,J-1)+Csom*VN(I,J)
            RV(I,J)=KLAD
            V(I,J)=VN(I,J)+DelT*KLAD


         ENDIF
 20   CONTINUE
C
      RETURN
      END
c  End of subroutine TILDE
C

c*****************************************************************
c**********************poisson solver *****************************
c*****************************************************************

      SUBROUTINE SOLVEP
      IMPLICIT NONE

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX,Y,YJ,DelY,RDX,RDY,
     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
 
      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /ORGA/ NF
      INTEGER NF(IMAX,JMAX)

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

c local variables
      INTEGER I,J
C
      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
     >               StrtP
      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
      REAL  Epsi,OmStrt,Alpha,NrmRhs

      CALL COEFF

      IF (IMilu.EQ.0) THEN
c--- solve Poisson equation with SOR
          CALL SOR
      ELSE
c--- solve with MILU
         CALL MILU
      ENDIF

C
C**** VELOCITY-UPDATE FOR PRESSURE GRADIENT
C
      DO 180 J=1,JM1Us
      DO 180 I=1,IM1Us
        IF (NF(I,J).LE.7 .AND. NF(I+1,J).LE.7)
     1    U(I,J)=U(I,J)-2.0*DelT*(P(I+1,J)-P(I,J))/(DelX(I+1)+DelX(I))
        IF (NF(I,J).LE.7 .AND. NF(I,J+1).LE.7)
     1    V(I,J)=V(I,J)-2.0*DelT*(P(I,J+1)-P(I,J))/(DelY(J+1)+DelY(J))
 180  CONTINUE

      CALL BCBND

      RETURN
      END
c  End of subroutine SOLVEP

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE COEFF
      IMPLICIT NONE

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY,RDX,RDY,
     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
c
      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /ORGA/ NF
      INTEGER NF(IMAX,JMAX)

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

      COMMON /COEFP/ CN, CS, CE, CW, CC, DIV
      REAL CN(IMAX,JMAX), CS(IMAX,JMAX), CE(IMAX,JMAX),
     1     CW(IMAX,JMAX), CC(IMAX,JMAX), DIV(IMAX,JMAX)

      COMMON /coef2/ CWI,CEI,CNJ,CSJ
      REAL CWI(IMAX),CEI(IMAX),CNJ(JMAX),CSJ(JMAX)

      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
     1     , POutN, POutE, POutS, POutW
      INTEGER WB, WL, WR, WT

      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
     >               StrtP
      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
      REAL  Epsi,OmStrt,Alpha,NrmRhs

      COMMON /CONDIF/ RU, RV,divuv,divmax
      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
      INTEGER divuv

c local variables

      INTEGER I,J
      REAL schaal
C
C***  only first time step: set coefficients for Poisson equation
C     scaling by SCHAAL makes coefficient matrix symmetric 
C
      IF (Cycle.EQ.1) THEN

         DO 100 J=2,JM1Us
         DO 100 I=2,IM1Us
            schaal=-DelX(I)*DelY(J)
            IF (NF(I,J).EQ.0) THEN
c--- full cell
               CS(I,J)=CSJ(J)*schaal
               CN(I,J)=CNJ(J)*schaal
               CE(I,J)=CEI(I)*schaal
               CW(I,J)=CWI(I)*schaal
               CC(I,J)=-CW(I,J)-CE(I,J)-CN(I,J)-CS(I,J)
c
c--- eliminate solid walls, inflow (u=prescribed) and outflow
c
               IF (NF(I+1,J).GE.8) CC(I,J)=CC(I,J)+CE(I,J)
               IF (NF(I+1,J).GE.7) CE(I,J)=0.0
 
               IF (NF(I-1,J).GE.8) CC(I,J)=CC(I,J)+CW(I,J)
               IF (NF(I-1,J).GE.7) CW(I,J)=0.0

               IF (NF(I,J+1).GE.8) CC(I,J)=CC(I,J)+CN(I,J)
               IF (NF(I,J+1).GE.7) CN(I,J)=0.0
  
               IF (NF(I,J-1).GE.8) CC(I,J)=CC(I,J)+CS(I,J)
               IF (NF(I,J-1).GE.7) CS(I,J)=0.0

            ELSE
c--- obstacle cell
               CC(I,J) = -schaal
               CN(I,J) = 0.0
               CS(I,J) = 0.0
               CE(I,J) = 0.0
               CW(I,J) = 0.0
               DIV(I,J) = 0.0
            ENDIF

100      CONTINUE

         IF (IMilu.EQ.0) THEN
c--- scale diagonal at 1.0 for SOR
            DO 500 J=2,JM1Us
            DO 500 I=2,IM1Us
               CN(I,J)=CN(I,J)/CC(I,J)
               CS(I,J)=CS(I,J)/CC(I,J)
               CE(I,J)=CE(I,J)/CC(I,J)
               CW(I,J)=CW(I,J)/CC(I,J)
 500        CONTINUE
         ENDIF

      ENDIF

C
C*** each time step compute right-hand side of Poisson equation
C
      DO 200 J=2,JM1Us
      DO 200 I=2,IM1Us
         IF (NF(I,J).EQ.0) THEN

C.... only full cells are considered
C

          IF (divuv.eq.1) THEN
             DIV(I,J)=( RDX(I)*(U(I,J)-U(I-1,J)) + 
     >                             RDY(J)*(V(I,J)-V(I,J-1)) )/DelT  
          ELSE
             DIV(I,J)=( RDX(I)*(RU(I,J)-RU(I-1,J)) + 
     >                             RDY(J)*(RV(I,J)-RV(I,J-1)) )   
          ENDIF

c--- outflow boundary condition p=POut (set in SETPAR)
c
            IF (nf(i-1,j).EQ.7) div(i,j)=div(i,j)-cwi(i)*POutW
            IF (nf(i+1,j).EQ.7) div(i,j)=div(i,j)-cei(i)*POutE
            IF (nf(i,j-1).EQ.7) div(i,j)=div(i,j)-csj(j)*POutS
            IF (nf(i,j+1).EQ.7) div(i,j)=div(i,j)-cnj(j)*POutN
c 
            schaal=-SIGN(1.0, cc(i,j)) *DelX(I)*DelY(J) 
            DIV(I,J)=DIV(I,J)*schaal

            IF (IMilu.EQ.0) DIV(I,J)=DIV(I,J)/CC(I,J)

         ENDIF
  200 CONTINUE 

C*** Compute norm of right-hand side;
C    to be used in SOR convergence criterion for Poisson equation.
C
      NrmRhs=0.0
      DO J=2,JM1Us
      DO I=2,IM1Us
         NrmRhs=NrmRhs + DIV(I,J)*DIV(I,J)
      ENDDO
      ENDDO
      NrmRhs = SQRT(NrmRhs)

      RETURN
      END
c  End of subroutine COEFF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SOR
      IMPLICIT NONE

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)      

      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX, RDY,
     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
     >               StrtP
      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
      REAL  Epsi,OmStrt,Alpha,NrmRhs      

      COMMON /COEFP/ CN, CS, CE, CW, CC, DIV
      REAL CN(IMAX,JMAX), CS(IMAX,JMAX), CE(IMAX,JMAX),
     1     CW(IMAX,JMAX), CC(IMAX,JMAX), DIV(IMAX,JMAX)

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

c local variables
      REAL OMEG,DELTA,DIFF,Epsi2,DEL0
      INTEGER i,j

C     EUCLIDEAN NORM IS USED and SCALING WITH RHS
C october 2003: only relative criterion 
C      Epsi2=Epsi*SQRT(FLOAT(IM2Us*JM2Us))
C      Epsi2=AMAX1(Epsi2, Epsi*NrmRhs)
       Epsi2=Epsi*NrmRhs

      OMEG=OmStrt
      Iter=0

      IF (StrtP.eq.0) then
        DO 10 j=1,JMaxUs
        DO 10 i=1,IMaxUs
	   p(i,j)=0.0
 10	CONTINUE
      ENDIF

 140  Iter=Iter+1
      DELTA=0.0

      DO 150 j=2,Jm1Us
      DO 150 i=2,Im1Us        
         DIFF = - P(i,j) + DIV(i,j)
     &           - CW(i,j)*P(i-1,j) - CE(i,j)*P(i+1,j )
     &           - CS(i,j)*P(i,j-1) - CN(i,j)*P(i,j+1 )
         P(i,j) = P(i,j) + OMEG*DIFF
         DELTA = DELTA + DIFF**2
         IF (Iter.EQ.1) DEL0=DELTA
 150  CONTINUE

      DELTA=SQRT(DELTA)
      
      IF (DELTA.GT.Epsi2.and.Iter.LT.ItMax) goto 140

      ItSum = ItSum + Iter

c      IF (DELTA.GT.0.1*NrmRhs .AND. DELTA.GT.0.1*DEL0) THEN  
c           WRITE(6,'(A17,2E13.4,1F8.4,3E13.4)')
c     >       'slow convergence ', NrmRhs, DEL0, OMEG,
c     >         DELTA, DELTA/NrmRhs, DELTA/DEL0
c      ENDIF

      RETURN
      END
c  End of subroutine SOR


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE MILU
      IMPLICIT NONE

      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX, RDY,
     1   CX,CY,IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /COEFP/ CN, CS, CE, CW, CC, DIV
      REAL CN(IMAX,JMAX), CS(IMAX,JMAX), CE(IMAX,JMAX),
     1     CW(IMAX,JMAX), CC(IMAX,JMAX), DIV(IMAX,JMAX)

      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
     >               StrtP
      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
      REAL  Epsi,OmStrt,Alpha,NrmRhs

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

      REAL rkr,rkrold,zaz,maxkmr,maxp,somw,somp,
     >     alf,beta,delta,cee,cnn,css,cww

      REAL r(IMAX,JMAX),z(IMAX,JMAX),az(IMAX,JMAX),
     >     w(IMAX,JMAX),d1(IMAX,JMAX),dsi(IMAX,JMAX)
      INTEGER i,j,M,N
      SAVE d1,dsi

C ##################        INIT         ###################

      M=IMaxUs
      N=JMaxUs
      DO 5 I=1,IMaxUs
        r(i,1)=0.0
        r(i,N)=0.0
        z(i,1)=0.0
        z(i,N)=0.0
        az(i,1)=0.0
        az(i,N)=0.0
        w(i,1)=0.0
        w(i,N)=0.0
        dsi(i,1)=0.0
        dsi(i,N)=0.0
    5 CONTINUE
      DO 10 J=1,JMaxUs
        r(1,j)=0.0
        r(M,j)=0.0
        z(1,j)=0.0
        z(M,j)=0.0
        az(1,j)=0.0
        az(M,j)=0.0
        w(1,j)=0.0
        w(M,j)=0.0
        dsi(1,j)=0.0
        dsi(M,j)=0.0
   10 CONTINUE

      IF (StrtP.eq.0) then
        DO 15 j=1,JMaxUs
        DO 15 i=1,IMaxUs
	   p(i,j)=0.0
   15   CONTINUE
      ENDIF

C ##################   Preconditioning   ###################
C Onderstaande loop kan veel efficienter
C Ook verderop is winst te behalen. Vectorisatie is mogelijk door
C langs diagonalen (zuidoost-noordwest) te nummeren

      IF (CYCLE.EQ.1) THEN

         beta=(1.0+5./(IM1Us*JM1Us))
         DO 20 j=2,JM1Us
         DO 20 i=2,IM1Us 
           cee=0.0
           cww=0.0
           cnn=0.0
           css=0.0
           IF (i.GT.2) THEN
              cee=ce(i-1,j)
              IF (j.LT.JM1Us) css=cs(i-1,j+1)
           ENDIF
           IF (j.GT.2) THEN
             cnn=cn(i,j-1)
             IF (i.LT.IM1Us) cww=cw(i+1,j-1)
           ENDIF
           Dsi(i,j)=1./(beta*Cc(i,j)
     >                      -Cw(i,j)*(cee+css)*Dsi(i-1,j)
     >                      -Cs(i,j)*(cnn+cww)*Dsi(i,j-1))
           d1(i,j)=Cc(i,j)-2./Dsi(i,j)
   20    CONTINUE

      ENDIF

C ##################  Conjugate gradient  ###################

      Iter = 0
      rkr=0.0
      DO 40 j=2,JM1Us
      DO 40 i=2,IM1Us
        r(i,j)=div(i,j)-(Cc(i,j)*p(i,j)+Cs(i,j)*p(i,j-1)+
     >         Cw(i,j)*p(i-1,j)+Cn(i,j)*p(i,j+1)+Ce(i,j)*p(i+1,j))
        r(i,j)=(r(i,j)-Cs(i,j)*r(i,j-1)-Cw(i,j)*r(i-1,j))*Dsi(i,j)
        z(i,j)=r(i,j)/Dsi(i,j)
        rkr=rkr+r(i,j)*z(i,j)
   40 CONTINUE

      somw=0.0
      somp=0.0
      DO j=2,JM1Us
      DO i=2,IM1Us
         somw=somw+z(i,j)*z(i,j)
         somp=somp+p(i,j)*p(i,j)
      ENDDO
      ENDDO
      maxkmr=SQRT(somw)
      maxp=SQRT(somp)
      delta=maxkmr/(maxp+1.0e-10)
      IF (delta.LT.1e-6*Epsi) GOTO 2000

      DO 45 j=2,JM1Us
      DO 45 i=2,IM1Us
        p(i,j)=p(i,j)/Dsi(i,j)+Cn(i,j)*p(i,j+1)+Ce(i,j)*p(i+1,j)
 45   CONTINUE

 1000 Iter=Iter+1
      rkrold=rkr

      DO 50 j=JM1Us,2,-1
      DO 50 i=IM1Us,2,-1
        az(i,j)=(z(i,j)-Cn(i,j)*az(i,j+1)-Ce(i,j)*az(i+1,j))*Dsi(i,j)
        w(i,j)=z(i,j)+d1(i,j)*az(i,j)
 50   CONTINUE

      zaz=0.0

      DO 60 j=2,JM1Us
      DO 60 i=2,IM1Us
        w(i,j)=(w(i,j)-Cs(i,j)*w(i,j-1)-Cw(i,j)*w(i-1,j))*Dsi(i,j)
        az(i,j)=az(i,j)+w(i,j)
        zaz=zaz+az(i,j)*z(i,j)
 60   CONTINUE

      alf=rkr/zaz
      rkr=0.0
      DO 70 j=2,JM1Us
      DO 70 i=2,IM1Us
        p(i,j)=p(i,j)+alf*z(i,j)
        r(i,j)=r(i,j)-alf*az(i,j)
        w(i,j)=r(i,j)/Dsi(i,j)
        rkr=rkr+r(i,j)*w(i,j)
   70 CONTINUE

      somw=0.0
      somp=0.0
      DO j=2,JM1Us
      DO i=2,IM1Us
         somw=somw+w(i,j)*w(i,j)
         somp=somp+p(i,j)*p(i,j)
      ENDDO
      ENDDO
      maxkmr=SQRT(somw)
      maxp=SQRT(somp)
      delta=maxkmr/(maxp+1.0e-10)

      beta=rkr/rkrold
      DO 80 j=2,JM1Us
      DO 80 i=2,IM1Us
        z(i,j)=w(i,j)+beta*z(i,j)
   80 CONTINUE

      IF (delta.GT.Epsi .AND. Iter.LT.ItMax) GOTO 1000

      DO 90 j=JM1Us,2,-1
      DO 90 i=IM1Us,2,-1
        p(i,j)=(p(i,j)-Cn(i,j)*p(i,j+1)-Ce(i,j)*p(i+1,j))*Dsi(i,j)
   90 CONTINUE

 2000 ItSum = ItSum + Iter

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE PRT(N)
      IMPLICIT NONE

      INTEGER IMAX,JMAX,NSMAX,SDTMAX,mmpnts,NRFMAX,NRHMAX
      PARAMETER (IMAX=130, JMAX=130, NSMAX=20, SDTMAX=10, 
     >           mmpnts=25, NRFMAX=10, NRHMAX=10)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX,RDY,
     1   CX,CY, IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
     1     , POutN, POutE, POutS, POutW
      INTEGER WB, WL, WR, WT

      COMMON /ORGA/ NF
      INTEGER NF(IMAX,JMAX)

      COMMON /ORGA2/ IMilu,Iter,ItMax,ItSum,Epsi,OmStrt,Alpha,NrmRhs,
     >               StrtP
      INTEGER IMilu,Iter,ItMax,ItSum,StrtP
      REAL  Epsi,OmStrt,Alpha,NrmRhs

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

      COMMON /PRINTS/ TPrt,PrtDT,TAdPrt,TCnf,NP
     >   ,PGnu,PAvs,Puvpf,PVelop,PForce
      REAL TPrt, PrtDT, TAdPrt, TCnf(2001)
      INTEGER NP,PGnu,PAvs,Puvpf,PVelop,PForce

      COMMON /TANK/ XMax, XMin, YMax, YMin
      REAL  XMax, XMin, YMax, YMin

      COMMON /FLUID/ Nu
      REAL  Nu

      COMMON /CONDIF/ RU,RV,divuv,divmax
      REAL RU(IMAX,JMAX),RV(IMAX,JMAX),divmax
      iNTEGER divuv

c local variables
      INTEGER i,j,n,imid,jstap
      REAL dx,dy

c----------------------------------------------------------------------
    1 FORMAT (/)
    2 FORMAT ('XMin=',F7.3,3X,'XMax=',F7.3,3X,'YMin=',F7.3,3X,
     C    'YMax=',F7.3)    
    3 FORMAT (I3,' + 2 X-INTERVALS , MEAN MESHSIZE =',F8.4)
    4 FORMAT (I3,' + 2 Y-INTERVALS , MEAN MESHSIZE =',F8.4)
    5 FORMAT ('TIME STEP DelT           =',E12.4)
    7 FORMAT ('BOUNDARY CONDITIONS (1=SLIP, 2=NO-SLIP, 7=OUTFLOW, ',
     C '8=INFLOW)',/,
     C '  RIGHT',i3,3X,'LEFT',i3,3X,'TOP',i3,3X,'BOTTOM',i3)
   71 FORMAT ('INFLOW AMPLITUDE UIn =',E11.4, ' VIn =',E11.4 )
    8 FORMAT ('KINEMATIC VISCOSITY Nu   =',E10.3)
   12 FORMAT ('PRESSURE CONVERGENCE Epsi=',E9.2)
   13 FORMAT ('UPWIND PARAMETER Alpha   =',F4.1)
   15 FORMAT (8F10.4)
   22 FORMAT ('T=',E12.5,2X,'div=',E12.5,2X,I4,3X,14(F9.5))
   41 FORMAT (1H ,48I2)
   42 FORMAT (1H ,17F9.4)
   43 FORMAT (1H ,20E12.4)

C--------------------------------------------------------------------
      IF (N .EQ. 1) THEN
C
C--- initializing print actions
C
         WRITE(6,1)
         WRITE(6,2) XMin,XMax,YMin,YMax
         WRITE(6,*)
         DX=(XMax-XMin)/IM2Us
         WRITE(6,3) IM2Us,DX
         WRITE(6,15) (X(I), I=1,IM1Us)
         WRITE(6,*)
         DY=(YMax-YMin)/JM2Us
         WRITE(6,4) JM2Us, DY
         WRITE(6,15) (Y(J), J=1,JM1Us)
         WRITE(6,*)
         WRITE(6,5) DelT
         WRITE(6,7) WR, WL, WT, WB
         WRITE(6,71) UIn,VIn
         WRITE(6,8) Nu
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,12) Epsi
         WRITE(6,13) Alpha
         WRITE(6,*)

C.... print initial configuration
         CALL PRTCNF(1)
         WRITE(6,1)

C--------------------------------------------------------------------
C        regular printout every PrtDt time steps
C--------------------------------------------------------------------
      ELSEIF (N .EQ. 2) THEN
C 
c** info for screen and 'cfd.out'
c--- horizontal velocity along x-mid plane is chosen here
c
         jstap=MAX0(JM2Us/6,1)
         IMID=IMaxUs/2+1
c         WRITE(6,22) T,divmax,Iter,(U(IMID,J),J=2,JM1Us,jstap)
         WRITE(6,22) T,divmax,Iter,U(IMID,2),U(IMID,JM2Us/2),
     +               U(IMID,JM2Us)
        
         IF (Puvpf.EQ.1) THEN
c... increase counter for number of calls of postpro-routines
        	 NP = NP + 1
C--- velocities, forces, fluxes and streamlines
c--  output of solution in all grid points
	         CALL PRTCNF(2)
         ENDIF   

c-- compute next small print time
         TPrt = TPrt + PrtDT
C
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C
      ELSEIF (N .EQ. 3) THEN
C
C*****  afsluitende uitvoer
C
c-- averages for time step and iterations
         WRITE(6,*)
         WRITE(6,'(a21,i6)') 'No. of time steps  = ', Cycle
         WRITE(6,'(a21,e9.2)') 'average time step  = ',
     >                                      (TFin-TStart)/Cycle
         WRITE(6,'(a21,f6.1)') 'average iterations = ',
     >                            FLOAT(ItSum)/FLOAT(Cycle)
         WRITE(6,*)         

c... increase counter for number of calls of postpro-routines
         NP = NP + 1 
c--- closing of files
C--- velocities, forces, fluxes and streamfunctions
         CALL PRTCNF(2)

      ENDIF

      RETURN
      END
c  End of subroutine PRT
c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      SUBROUTINE PRTCNF(N)
      IMPLICIT NONE
      INTEGER N
c
      INTEGER IMAX,JMAX
      PARAMETER (IMAX=130, JMAX=130)
c
      COMMON /GRIDAR/  X,XI,DelX, Y,YJ,DelY, RDX,RDY,
     1   CX,CY, IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us
      REAL X(IMAX),XI(IMAX),DelX(IMAX), Y(JMAX),YJ(JMAX),DelY(JMAX),
     1   RDX(IMAX),RDY(JMAX),CX,CY
      INTEGER IMaxUs,JMaxUs,IM1Us,JM1Us,IM2Us,JM2Us

      COMMON /TANK/ XMax, XMin, YMax, YMin
      REAL  XMax, XMin, YMax, YMin

      COMMON /TIMES/ Cycle, T, DelT, DelTMx, TFin, TStart
      REAL  T, DelT, DelTMx, TFin, TStart
      INTEGER Cycle

      COMMON /PRINTS/ TPrt,PrtDT,TAdPrt,TCnf,NP
     >   ,PGnu,PAvs,Puvpf,PVelop,PForce
      REAL TPrt, PrtDT, TAdPrt, TCnf(2001)
      INTEGER NP,PGnu,PAvs,Puvpf,PVelop,PForce

      COMMON /PHYS/  F,U,UN,V,VN,P,PN,VMAX,PMIN,PMAX
      REAL F(IMAX,JMAX),U(0:IMAX,JMAX),UN(0:IMAX,JMAX),
     1     V(IMAX,0:JMAX),VN(IMAX,0:JMAX),
     2     P(IMAX,JMAX),PN(IMAX,JMAX),VMAX,PMIN,PMAX

      COMMON /DXRL/ DXRI,DXLI,DYTJ,DYBJ
      REAL DXRI(IMAX),DXLI(IMAX),DYTJ(JMAX),DYBJ(JMAX)

      COMMON /BNDCON/ WB, WL, WR, WT, SgnB, SgnL, SgnR, SgnT,
     1        UIn, VIn, FreqIn, POutN, POutE, POutS, POutW
      REAL  SgnB, SgnL, SgnR, SgnT, UIn, VIn, FreqIn
     1     , POutN, POutE, POutS, POutW
      INTEGER WB, WL, WR, WT

      COMMON /ORGA/ NF
      INTEGER NF(IMAX,JMAX)

      COMMON /WORKAR/ PSI
      REAL PSI(IMAX,JMAX)

c  local variables
      REAL UC(IMAX,JMAX),VC(IMAX,JMAX)
      INTEGER I,J
      CHARACTER*4 nummer
c
      IF (N.EQ.1) THEN

         OPEN (UNIT=28, FILE='config.dat')
         WRITE(28,'(I3,4F12.5)') 0, XMin, XMax, YMin, YMax
         WRITE(28,'(3I5)') IMaxUs, JMaxUs
         WRITE(28,'(10F12.5)') (XI(I),I=1,IMaxUs)
         WRITE(28,'(10F12.5)') (YJ(J),J=1,JMaxUs)
         DO J=1,JMaxUs
            WRITE(28,'(32I3)') (NF(I,J),I=1,IMaxUs)
         ENDDO         
         WRITE(28,'(6F15.5)') 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
         CLOSE(28)

      ELSEIF (N.EQ.2) THEN

         IF (NP .GT. 2001) THEN
            WRITE(6,*) 'too many uvpf-files created'
            RETURN
         ENDIF

         WRITE (nummer,'(I4)') NP
         DO i=1,4
            IF (nummer(i:i).EQ.' ') nummer(i:i)='0'
         ENDDO

c--- average velocities to cell centers

         DO J=1,JMaxUs
         DO I=1,IMaxUs
            UC(I,J)=0.0
            VC(I,J)=0.0
            IF (NF(I,J).NE.9) THEN
               UC(I,J)=0.5*(U(I-1,J)+U(I,J))
               VC(I,J)=0.5*(V(I,J-1)+V(I,J))
               vmax=AMAX1(SQRT(UC(I,J)**2+VC(I,J)**2),vmax)
               pmin=AMIN1(p(i,j),pmin)
               pmax=AMAX1(p(i,j),pmax)
            ENDIF
            PSI(I,J)=0.0
         ENDDO
         ENDDO

c-- slip boundary velocities
         DO I=2,IM1Us
            J=1
            IF (NF(I,J).EQ.9 .AND. WB.EQ.1) UC(I,J)=UC(I,J+1)
            J=JMaxUs
            IF (NF(I,J).EQ.9 .AND. WT.EQ.1) UC(I,J)=UC(I,J-1)
         ENDDO

         DO J=2,JM1Us
            I=1
            IF (NF(I,J).EQ.9 .AND. WL.EQ.1) VC(I,J)=VC(I+1,J)
            I=IMaxUs
            IF (NF(I,J).EQ.9 .AND. WR.EQ.1) VC(I,J)=VC(I-1,J)
         ENDDO

	 PSI(1,1)=0.0
	 DO I=2,IM1Us
	   PSI(I,1)=PSI(I-1,1)-DelX(I)*V(I,1)
	 ENDDO

	 DO I=1,IM1Us
	 DO J=2,JM1Us
	   PSI(I,J)=PSI(I,J-1)+DelY(J)*U(I,J)
	 ENDDO
	 ENDDO

         OPEN (UNIT=28, FILE='uvpf'//nummer//'.dat')

         WRITE(28,'(F12.5)') T
         DO J=1,JMaxUs
         DO I=1,IMaxUs
c.. small numbers set to 0 to avoid format problems on Cray
            IF (ABS(UC(I,J)) .LE. 1.E-20) UC(I,J)=0.0
            IF (ABS(VC(I,J)) .LE. 1.E-20) VC(I,J)=0.0
            IF (ABS(P(I,J)) .LE. 1.E-20) P(I,J)=0.0
            IF (ABS(PSI(I,J)) .LE. 1.E-20) PSI(I,J)=0.0

            WRITE(28,'(3E12.4, F7.3, E12.4)') 
     >         UC(I,J),VC(I,J),P(I,J),F(I,J),PSI(I,J)
         ENDDO
         ENDDO

         CLOSE(28)

         WRITE(6,'(A13,A4,A21)') 'plotfile UVPF', nummer, 
     >                                 '.DAT has been created '

      ENDIF

      RETURN
      END





















