C $Header: /home/ross/cvsroot/parnum/xc.f,v 1.1 2005/03/16 14:25:22 ross Exp $
      SUBROUTINE GASCOR(RHO1,RHO2,ECOR,VC1,VC2,UCOR)
      USE CONSTANTS
      IMPLICIT NONE
      REAL(KIND(1D0)),INTENT(IN) :: RHO1,RHO2
      REAL(KIND(1D0)),INTENT(OUT):: ECOR,VC1,VC2,UCOR
      REAL(KIND(1D0)),PARAMETER :: PGAM = 0.5198421D0
      REAL(KIND(1D0)),PARAMETER :: PFZZ = 1.709921D0
      REAL(KIND(1D0)) PRS,PZET,PF
      REAL(KIND(1D0)) PA,PA1,PB1,PB2,PB3,PB4,PP
      REAL(KIND(1D0)) PP1,PQ0,PRS12,PRS32,PRSP,PQ1,PQ2,PQ3
      REAL(KIND(1D0)) PEU,PEURS,PEP,PEPRS
      REAL(KIND(1D0)) PALFM,ALFRSM,PALFC,PZ4,PEC,PECRS,PFZ,PECZET,PCOMM
      REAL(KIND(1D0)) PTC,PTCOR
C**** PW91 CORRELATION LSDA,
      PRS = (0.75D0/PI/(RHO1+RHO2+EPS))**THRD
      PZET = (RHO1-RHO2)/(RHO1+RHO2+EPS)/(1.D0+EPS)
      PF = ((1.D0+PZET)**THRD4+(1.D0-PZET)**THRD4-2.D0)/PGAM
      PA = 0.0310907D0
      PA1 = 0.21370D0
      PB1 = 7.5957D0
      PB2 = 3.5876D0
      PB3 = 1.6382D0
      PB4 = 0.49294D0
      PP = 1.00D0
      PP1 = PP + 1.D0
      PQ0 = -2.D0*PA*(1.D0+PA1*PRS)
      PRS12 = DSQRT(PRS)
      PRS32 = PRS12**3
      PRSP = PRS**PP
      PQ1 = 2.D0*PA*(PB1*PRS12+PB2*PRS+PB3*PRS32+PB4*PRS*PRSP)
      PQ2 = DLOG(1.D0+1.D0/PQ1)
      PEU = PQ0*PQ2
      PQ3 = PA*(PB1/PRS12+2.D0*PB2+3.D0*PB3*PRS12+2.D0*PB4*PP1*PRSP)
      PEURS = -2.D0*PA*PA1*PQ2-PQ0*PQ3/(PQ1**2+PQ1)
      PA = 0.01554535D0
      PA1 = 0.20548D0
      PB1 = 14.1189D0
      PB2 = 6.1977D0
      PB3 = 3.3662D0
      PB4 = 0.62517D0
      PP = 1.00D0
      PP1 = PP + 1.D0
      PQ0 = -2.D0*PA*(1.D0+PA1*PRS)
      PRSP = PRS**PP
      PQ1 = 2.D0*PA*(PB1*PRS12+PB2*PRS+PB3*PRS32+PB4*PRS*PRSP)
      PQ2 = DLOG(1.D0+1.D0/PQ1)
      PEP = PQ0*PQ2
      PQ3 = PA*(PB1/PRS12+2.D0*PB2+3.D0*PB3*PRS12+2.D0*PB4*PP1*PRSP)
      PEPRS = -2.D0*PA*PA1*PQ2-PQ0*PQ3/(PQ1**2+PQ1)
      PA = 0.0168869D0
      PA1 = 0.11125D0
      PB1 = 10.357D0
      PB2 = 3.6231D0
      PB3 = 0.88026D0
      PB4 = 0.49671D0
      PP = 1.00D0
      PP1 = PP + 1.D0
      PQ0 = -2.D0*PA*(1.D0+PA1*PRS)
      PRSP = PRS**PP
      PQ1 = 2.D0*PA*(PB1*PRS12+PB2*PRS+PB3*PRS32+PB4*PRS*PRSP)
      PQ2 = DLOG(1.D0+1.D0/PQ1)
      PALFM = PQ0*PQ2
      PQ3 = PA*(PB1/PRS12+2.D0*PB2+3.D0*PB3*PRS12+2.D0*PB4*PP1*PRSP)
      ALFRSM = -2.D0*PA*PA1*PQ2-PQ0*PQ3/(PQ1**2+PQ1)
      PALFC = -PALFM
      PZ4 = PZET**4
      PEC = PEU*(1.D0-PF*PZ4)+PEP*PF*PZ4-PALFM*PF*(1.D0-PZ4)/PFZZ
       ECOR = (RHO1+RHO2) * PEC
      PECRS = PEURS*(1.D0-PF*PZ4)+PEPRS*PF*PZ4-ALFRSM*PF*(1.D0-PZ4)/PFZZ
      PFZ = THRD4*((1.D0+PZET)**THRD-(1.D0-PZET)**THRD)/PGAM
      PECZET = 4.D0*(PZET**3)*PF*(PEP-PEU+PALFM/PFZZ)+PFZ*(PZ4*PEP-PZ4
     1         *PEU-(1.D0-PZ4)*PALFM/PFZZ)
      PCOMM = PEC -PRS*PECRS/3.D0-PZET*PECZET
       VC1 = PCOMM + PECZET
       VC2 = PCOMM - PECZET
C**** END PW91 CORRELATION LSDA.
      PTC = -4.D0*PEC + 1.5D0*((1.D0+PZET)* VC1+(1.D0-PZET)* VC2)
      PTCOR = (RHO1+RHO2) * PTC
       UCOR =   ECOR-PTCOR
C**** END PW91 CORRELATION KINETIC AND POTENTIAL ENERGIES.
      RETURN
      END
C
C
C
      SUBROUTINE XCGGA(RHOS,EXL,ECL,UCL,EXGGA,ECGGA,VXC)
C
C     GENERALIZED GRADIENT APPROXIMATION EXCHANGE-CORRELATION POTENTIAL
C
      USE CONSTANTS
      USE GRID
      USE KEYWORDS
      USE SIZES
      IMPLICIT NONE
      REAL(KIND(1D0)),INTENT(IN) :: RHOS(NNN,2)
      REAL(KIND(1D0)),INTENT(OUT):: EXL,ECL,UCL,EXGGA,ECGGA,VXC(NNN,2)
      REAL(KIND(1D0)),PARAMETER :: BET=0.0036D0, GAM=0.004D0
      REAL(KIND(1D0)),PARAMETER :: C12=0.0032D0, C11=0.20D0
      REAL(KIND(1D0)),PARAMETER :: THRD8=8.D0/3.D0, THRD16=16.D0/3.D0
      REAL(KIND(1D0)) RO1XYZ(0:9,NNN),RO2XYZ(0:9,NNN)
      REAL(KIND(1D0)) RO1MIN,RO2MIN,RAISE1,RAISE2
      REAL(KIND(1D0)) RHO1,DRHO12,W1,D2RHO1,D41
      REAL(KIND(1D0)) RHO2,DRHO22,W2,D2RHO2,D42
      REAL(KIND(1D0)) F0,F1,DENOM,G0,G1,G2
      REAL(KIND(1D0)) RHO43,EC01,VC01,DUM,UC,EC02,VC02,ECOR,VC1,VC2
      REAL(KIND(1D0)) ECOR0,ECOR1,ECOR2,DUM1,DUM2,UCOR0,EC12
      INTEGER KKK,I
C
      CALL FDRHO(RHOS(1,1),3,RO1XYZ)
      IF(NSPINS.EQ.2)THEN
      CALL FDRHO(RHOS(1,2),3,RO2XYZ)
      ELSE
      DO 10 KKK=1,NNN
      DO 10 I=0,9
10    RO2XYZ(I,KKK)=RO1XYZ(I,KKK)
      END IF
      RO1MIN=1.D11
      RO2MIN=1.D11
      DO 20 KKK=1,NNN
      VXC(KKK,1)=0.D0
      VXC(KKK,2)=0.D0
      IF(RO1XYZ(0,KKK).LT.RO1MIN)RO1MIN=RO1XYZ(0,KKK)
      IF(RO2XYZ(0,KKK).LT.RO2MIN)RO2MIN=RO2XYZ(0,KKK)
20    CONTINUE
      IF(RO1MIN.GE.EPS)RAISE1=0.D0
      IF(RO1MIN.LT.EPS)RAISE1=EPS-RO1MIN
      IF(RO2MIN.GE.EPS)RAISE2=0.D0
      IF(RO2MIN.LT.EPS)RAISE2=EPS-RO2MIN
C
      EXL=0.D0
      ECL=0.D0
      UCL=0.D0
      EXGGA=0.D0
      ECGGA=0.D0
      DO 100 KKK=1,NNN
      RHO1=RO1XYZ(0,KKK)+RAISE1
      DRHO12=RO1XYZ(1,KKK)**2+RO1XYZ(2,KKK)**2+RO1XYZ(3,KKK)**2
      W1=DRHO12/RHO1**THRD8
      D2RHO1=RO1XYZ(4,KKK)+RO1XYZ(5,KKK)+RO1XYZ(6,KKK)
      D41 = 2.D0*( RO1XYZ(1,KKK)**2*RO1XYZ(4,KKK)
     +           + RO1XYZ(2,KKK)**2*RO1XYZ(5,KKK)
     +           + RO1XYZ(3,KKK)**2*RO1XYZ(6,KKK) )
     +    + 4.D0*( RO1XYZ(2,KKK)*RO1XYZ(3,KKK)*RO1XYZ(7,KKK)
     +           + RO1XYZ(1,KKK)*RO1XYZ(3,KKK)*RO1XYZ(8,KKK)
     +           + RO1XYZ(1,KKK)*RO1XYZ(2,KKK)*RO1XYZ(9,KKK) )
      RHO2=RO2XYZ(0,KKK)+RAISE2
      DRHO22=RO2XYZ(1,KKK)**2+RO2XYZ(2,KKK)**2+RO2XYZ(3,KKK)**2
      W2=DRHO22/RHO2**THRD8
      D2RHO2=RO2XYZ(4,KKK)+RO2XYZ(5,KKK)+RO2XYZ(6,KKK)
      D42 = 2.D0*( RO2XYZ(1,KKK)**2*RO2XYZ(4,KKK)
     +           + RO2XYZ(2,KKK)**2*RO2XYZ(5,KKK)
     +           + RO2XYZ(3,KKK)**2*RO2XYZ(6,KKK) )
     +    + 4.D0*( RO2XYZ(2,KKK)*RO2XYZ(3,KKK)*RO2XYZ(7,KKK)
     +           + RO2XYZ(1,KKK)*RO2XYZ(3,KKK)*RO2XYZ(8,KKK)
     +           + RO2XYZ(1,KKK)*RO2XYZ(2,KKK)*RO2XYZ(9,KKK) )
C     EXCHANGE:
      F0=RHO1**THRD4
      F1=THRD4*RHO1**THRD
      DENOM=1.D0+GAM*W1
      G0=-CEX-BET*W1/DENOM
      G1=-BET/DENOM**2
      G2= 2.D0*BET*GAM/DENOM**3
      VXC(KKK,1) = VXC(KKK,1)
     +           + F1*G0 + THRD8*F0*G1*W1/RHO1 - 2.D0*F1*G1*W1
     +           - 2.D0*F0*G1*D2RHO1/RHO1**THRD8
     +           + THRD16*F0*G2*W1*W1/RHO1 - 2.D0*F0*G2*D41/RHO1**THRD16
      RHO43=(RHOS(KKK,1)+EPS)**THRD4
      EXL=EXL-CEX*WINTS(KKK)*RHO43
      EXGGA=EXGGA+WINTS(KKK)*RHO43*G0
      F0=RHO2**THRD4
      F1=THRD4*RHO2**THRD
      DENOM=1.D0+GAM*W2
      G0=-CEX-BET*W2/DENOM
      G1=-BET/DENOM**2
      G2= 2.D0*BET*GAM/DENOM**3
      VXC(KKK,2) = VXC(KKK,2)
     +           + F1*G0 + THRD8*F0*G1*W2/RHO2 - 2.D0*F1*G1*W2
     +           - 2.D0*F0*G1*D2RHO2/RHO2**THRD8
     +           + THRD16*F0*G2*W2*W2/RHO2 - 2.D0*F0*G2*D42/RHO2**THRD16
      RHO43=(RHOS(KKK,2)+EPS)**THRD4
      EXL=EXL-CEX*WINTS(KKK)*RHO43
      EXGGA=EXGGA+WINTS(KKK)*RHO43*G0
C     SAME-SPIN CORRELATION:
      CALL GASCOR(RHO1,0.D0,EC01,VC01,DUM,UC)
      F0=EC01
      F1=VC01
      DENOM=1.D0+C11*W1
      G0= 1.D0/DENOM
      G1=-C11/DENOM**2
      G2= 2.D0*C11**2/DENOM**3
      VXC(KKK,1) = VXC(KKK,1)
     +           + F1*G0 + THRD8*F0*G1*W1/RHO1 - 2.D0*F1*G1*W1
     +           - 2.D0*F0*G1*D2RHO1/RHO1**THRD8
     +           + THRD16*F0*G2*W1*W1/RHO1 - 2.D0*F0*G2*D41/RHO1**THRD16
      CALL GASCOR(RHOS(KKK,1),0.D0,ECOR1,DUM1,DUM2,UC)
      ECGGA=ECGGA+WINTS(KKK)*ECOR1*G0
      CALL GASCOR(RHO2,0.D0,EC02,VC02,DUM,UC)
      F0=EC02
      F1=VC02
      DENOM=1.D0+C11*W2
      G0= 1.D0/DENOM
      G1=-C11/DENOM**2
      G2= 2.D0*C11**2/DENOM**3
      VXC(KKK,2) = VXC(KKK,2)
     +           + F1*G0 + THRD8*F0*G1*W2/RHO2 - 2.D0*F1*G1*W2
     +           - 2.D0*F0*G1*D2RHO2/RHO2**THRD8
     +           + THRD16*F0*G2*W2*W2/RHO2 - 2.D0*F0*G2*D42/RHO2**THRD16
      CALL GASCOR(RHOS(KKK,2),0.D0,ECOR2,DUM1,DUM2,UC)
      ECGGA=ECGGA+WINTS(KKK)*ECOR2*G0
C     OPPOSITE-SPIN CORRELATION:
      CALL GASCOR(RHO1,RHO2,ECOR,VC1,VC2,UC)
      EC12=ECOR-EC01-EC02
      F0=EC12
      F1=VC1-VC01
      DENOM=1.D0+C11*(W1+W2)
      G0= 1.D0/DENOM
      G1=-C11/DENOM**2
      G2= 2.D0*C11**2/DENOM**3
      VXC(KKK,1) = VXC(KKK,1)
     +           + F1*G0 + THRD8*F0*G1*W1/RHO1 - 2.D0*F1*G1*W1
     +           - 2.D0*F0*G1*D2RHO1/RHO1**THRD8
     +           + THRD16*F0*G2*W1*W1/RHO1 - 2.D0*F0*G2*D41/RHO1**THRD16
      F0=EC12
      F1=VC2-VC02
      DENOM=1.D0+C11*(W1+W2)
      G0= 1.D0/DENOM
      G1=-C11/DENOM**2
      G2= 2.D0*C11**2/DENOM**3
      VXC(KKK,2) = VXC(KKK,2)
     +           + F1*G0 + THRD8*F0*G1*W2/RHO2 - 2.D0*F1*G1*W2
     +           - 2.D0*F0*G1*D2RHO2/RHO2**THRD8
     +           + THRD16*F0*G2*W2*W2/RHO2 - 2.D0*F0*G2*D42/RHO2**THRD16
      CALL GASCOR(RHOS(KKK,1),RHOS(KKK,2),ECOR0,DUM1,DUM2,UCOR0)
      ECL=ECL+WINTS(KKK)*ECOR0
      UCL=UCL+WINTS(KKK)*UCOR0
      ECGGA=ECGGA+WINTS(KKK)*(ECOR0-ECOR1-ECOR2)*G0
100   CONTINUE
C
      RETURN
      END
C
C
C$PRAGMA SUN OPT=1
      SUBROUTINE XCMODS(RHO,RHOS,TAUS,EXDENS,
     +                  EX88,EC8812,EC8811,
     +                  EXBR,ECBR12,ECBR11,
     +                  FXBR,FCBR12,FCBR11,
     +                  EC95OPP,EC95PAR,
     +                  CND12,CND21,
     +                  CD,CDC,
     +                  CNDP,CDP,CDPC,
     +                  BRDINT,BRDIP2)
C
C     NEW HARTREE-FOCK BASED XC MODELS.
C     2005 Mar 12 --- added Bc95 in so we can compute BB1K --- rmd
C     2005 Mar 16 --- Bc95 might not be right yet; literature comparisons off
C
C     TODO:  Deallocate the members of the Potentials module
C     before we enter here?
C
      USE CONSTANTS
      USE ERRORHANDLING
      USE GRID
      USE KEYWORDS
      USE MPI
      USE SIZES
      IMPLICIT NONE
      ! Arguments
      REAL(KIND(1D0)),INTENT(IN) :: RHO(NNN),RHOS(NNN,2),TAUS(NNN,2)
      REAL(KIND(1D0)),INTENT(IN) :: EXDENS(NNN,2)
      REAL(KIND(1D0)),INTENT(OUT):: EX88,EC8812,EC8811
      REAL(KIND(1D0)),INTENT(OUT):: EC95OPP,EC95PAR
      REAL(KIND(1D0)),INTENT(OUT):: EXBR,ECBR12,ECBR11
      REAL(KIND(1D0)),INTENT(OUT):: FXBR,FCBR12,FCBR11
      REAL(KIND(1D0)),INTENT(OUT):: CND12,CND21,CD,CDC,CNDP,CDP,CDPC
      REAL(KIND(1D0)),INTENT(OUT):: BRDINT,BRDIP2(NNN)
      ! Locals
      REAL(KIND(1D0)) RHOXYZ(0:9,NNN),TAUPD(NNN),WEIZS(NNN,2)
      REAL(KIND(1D0)) D2RHO(NNN),DSIGS(NNN,2),QUADS(NNN,2)
      REAL(KIND(1D0)) UXPS(NNN,2),XLNS(NNN,2)
      REAL(KIND(1D0)) RHOSIG,XLN,DRHOSQ,DRHO,RHO43,XL,DLSS,ASNH
      REAL(KIND(1D0)) X88,RFG1,Z11,C8811,DUM,XBR1,RBR1,C1BR11,XBRG,SBR1
      REAL(KIND(1D0)) C2BR11,RLH1,C0BR11,DELN,X,B,F1,F2,F3,F4,F5
      REAL(KIND(1D0)) COEF1,COEF2,RFG2,Z22,C8822,RBR2,C1BR22,SBR2,C2BR22
      REAL(KIND(1D0)) RLH2,C0BR22,FLIP12,FLIP21,FLIP,Z12,C8812,CBR12
      REAL(KIND(1D0)) CD12
      INTEGER KKK
      ! Locals (Bc95 components)
      REAL(KIND(1D0)),PARAMETER :: CC95OPP=0.0031, CC95PAR=0.038
      REAL(KIND(1D0)),PARAMETER :: THRD5=5.D0/3.D0
      REAL(KIND(1D0)) EC12UEG,EC11UEG,EC22UEG,DSIGUEG,EC95D
      REAL(KIND(1D0)) CDSIGUEG,DUM1,DUM2,DUM3,CHI1,CHI2
      CDSIGUEG = 0.6D0*(6.D0*PI**2)**(2.D0/3.D0) ! could be made constant
C
      EX88=0.D0
      EC8812=0.D0
      EC8811=0.D0
      EXBR=0.D0
      ECBR12=0.D0
      ECBR11=0.D0
      FXBR=0.D0
      FCBR12=0.D0
      FCBR11=0.D0
      EC95OPP=0.D0
      EC95PAR=0.D0
      CND12=0.D0
      CND21=0.D0
      CD=0.D0
      CDC=0.D0
      CNDP=0.D0
      CDP=0.D0
      CDPC=0.D0
      BRDINT=0.D0
      BRDIP2(:)=0.D0
C
      CALL TAUWZL(1,TAUS(:,1),TAUPD,WEIZS(:,1),D2RHO)
      if (myRank.eq.0) then
      DO 100 KKK=1,NNN
      RHOSIG=RHOS(KKK,1)
      IF(RHOSIG.LT.SMALL)GO TO 100
      DSIGS(KKK,1)=TAUPD(KKK)-WEIZS(KKK,1)
      UXPS(KKK,1)=-2.D0*EXDENS(KKK,1)/RHOSIG
      QUADS(KKK,1)=(D2RHO(KKK)-2.D0*DSIGS(KKK,1))/6.D0
      CALL XLNORM(RHOSIG,QUADS(KKK,1),UXPS(KKK,1),XLN)
      if (ERRFLG) return
      XLNS(KKK,1)=XLN
100   CONTINUE
      end if ! myRank==0
C
      IF(NSPINS.EQ.2)THEN
      CALL TAUWZL(2,TAUS(:,2),TAUPD,WEIZS(:,2),D2RHO)
      ELSE
      DO 199 KKK=1,NNN
199   WEIZS(KKK,2)=WEIZS(KKK,1)
      END IF
C
      if (myRank.eq.0) then
      DO 200 KKK=1,NNN
      RHOSIG=RHOS(KKK,2)
      IF(RHOSIG.LT.SMALL)GO TO 200
      DSIGS(KKK,2)=TAUPD(KKK)-WEIZS(KKK,2)
      UXPS(KKK,2)=-2.D0*EXDENS(KKK,2)/RHOSIG
      QUADS(KKK,2)=(D2RHO(KKK)-2.D0*DSIGS(KKK,2))/6.D0
      CALL XLNORM(RHOSIG,QUADS(KKK,2),UXPS(KKK,2),XLN)
      if (ERRFLG) return
      XLNS(KKK,2)=XLN
200   CONTINUE
C
      DO 300 KKK=1,NNN
C
      IF (RHOS(KKK,1).LT.SMALL.AND.RHOS(KKK,2).LT.SMALL) THEN
      GO TO 300
C
      ELSE IF (RHOS(KKK,2).LT.SMALL) THEN
      DRHOSQ=4.D0*RHOS(KKK,1)*WEIZS(KKK,1)
      DRHO=DSQRT(DRHOSQ)
      RHO43=RHOS(KKK,1)**THRD4
      XL=-CEX*RHO43
      DLSS=DRHO/RHO43
      ASNH=DLOG(DLSS+DSQRT(DLSS**2+1.D0))
      X88=XL-0.0042D0*DRHO*DLSS/(1.D0+6.D0*0.0042D0*DLSS*ASNH)
      EX88=EX88+WINTS(KKK)*X88
      RFG1=-0.5D0*RHOS(KKK,1)/X88
      Z11=0.88D0*2.D0*RFG1
      C8811=-0.01D0*RHOS(KKK,1)*DSIGS(KKK,1)*Z11**4*(1.D0-2.D0/Z11*
     C                                          DLOG(1.D0+0.5D0*Z11))
      EC8811=EC8811+WINTS(KKK)*C8811
      CALL XMODEL(RHOS(KKK,1),QUADS(KKK,1),1.D0,DUM,XBR1)
      if (ERRFLG) return
      EXBR=EXBR+WINTS(KKK)*XBR1
      RBR1=-0.5D0*RHOS(KKK,1)/XBR1
      Z11=0.88D0*2.D0*RBR1
      C1BR11=-0.01D0*RHOS(KKK,1)*DSIGS(KKK,1)*Z11**4
     +              *(1.D0-2.D0/Z11*DLOG(1.D0+0.5D0*Z11))
      ECBR11=ECBR11+WINTS(KKK)*C1BR11
      CALL XMODEL(RHOS(KKK,1),QUADS(KKK,1),XLNS(KKK,1),
     +                                     UXPS(KKK,1),XBRG)
      if (ERRFLG) return
      FXBR=FXBR+WINTS(KKK)*XBRG
      SBR1=-0.5D0*RHOS(KKK,1)/XBRG
      Z11=0.88D0*2.D0*SBR1
      C2BR11=-0.01D0*RHOS(KKK,1)*DSIGS(KKK,1)*Z11**4
     +              *(1.D0-2.D0/Z11*DLOG(1.D0+0.5D0*Z11))
      FCBR11=FCBR11+WINTS(KKK)*C2BR11
      RLH1=XLNS(KKK,1)*SBR1
      Z11=0.88D0*2.D0*RLH1
      C0BR11=-0.01D0*RHOS(KKK,1)*DSIGS(KKK,1)*Z11**4
     +              *(1.D0-2.D0/Z11*DLOG(1.D0+0.5D0*Z11))
      CDP=CDP+WINTS(KKK)*C0BR11
      DELN=1.D0-XLNS(KKK,1)
      CALL BHOLE(RHOS(KKK,1),QUADS(KKK,1),XLNS(KKK,1),X,B)
      if (ERRFLG) return
      CALL MOMS(XLNS(KKK,1),X,B,F1,F2,F3,F4,F5)
      BRDIP2(KKK)=BRDIP2(KKK)+RHOS(KKK,1)*B**2
      BRDINT=BRDINT+WINTS(KKK)*RHOS(KKK,1)*B**2
      F3=F3/RHOS(KKK,1)
      F4=F4/RHOS(KKK,1)
      COEF1=DELN/(4.D0*PI*F4)
      COEF2=DMIN1(DSIGS(KKK,1)/3.D0,COEF1)
      CNDP=CNDP-WINTS(KKK)*2.D0*PI*COEF2*RHOS(KKK,1)*F3
      ! Bc95:
      call gascor(rhos(kkk,1),0.d0,ec11ueg,dum1,dum2,dum3)
      dsigueg = cdsigueg*rhos(kkk,1)**thrd5
      ec95d = ec11ueg*dsigs(kkk,1)/dsigueg/(1.d0+cc95par*dlss**2)**2
      ec95par = ec95par + wints(kkk)*ec95d
      ! end of Bc95 for RHOS(KKK,2).LT.SMALL
C
      ELSE IF (RHOS(KKK,1).LT.SMALL) THEN
      DRHOSQ=4.D0*RHOS(KKK,2)*WEIZS(KKK,2)
      DRHO=DSQRT(DRHOSQ)
      RHO43=RHOS(KKK,2)**THRD4
      XL=-CEX*RHO43
      DLSS=DRHO/RHO43
      ASNH=DLOG(DLSS+DSQRT(DLSS**2+1.D0))
      X88=XL-0.0042D0*DRHO*DLSS/(1.D0+6.D0*0.0042D0*DLSS*ASNH)
      EX88=EX88+WINTS(KKK)*X88
      RFG2=-0.5D0*RHOS(KKK,2)/X88
      Z22=0.88D0*2.D0*RFG2
      C8822=-0.01D0*RHOS(KKK,2)*DSIGS(KKK,2)*Z22**4*(1.D0-2.D0/Z22*
     C                                          DLOG(1.D0+0.5D0*Z22))
      EC8811=EC8811+WINTS(KKK)*C8822
      CALL XMODEL(RHOS(KKK,2),QUADS(KKK,2),1.D0,DUM,XBR1)
      if (ERRFLG) return
      EXBR=EXBR+WINTS(KKK)*XBR1
      RBR2=-0.5D0*RHOS(KKK,2)/XBR1
      Z22=0.88D0*2.D0*RBR2
      C1BR22=-0.01D0*RHOS(KKK,2)*DSIGS(KKK,2)*Z22**4
     +              *(1.D0-2.D0/Z22*DLOG(1.D0+0.5D0*Z22))
      ECBR11=ECBR11+WINTS(KKK)*C1BR22
      CALL XMODEL(RHOS(KKK,2),QUADS(KKK,2),XLNS(KKK,2),
     +                                     UXPS(KKK,2),XBRG)
      FXBR=FXBR+WINTS(KKK)*XBRG

      Z22=0.88D0*2.D0*SBR2
      C2BR22=-0.01D0*RHOS(KKK,2)*DSIGS(KKK,2)*Z22**4
     +              *(1.D0-2.D0/Z22*DLOG(1.D0+0.5D0*Z22))
      FCBR11=FCBR11+WINTS(KKK)*C2BR22
      RLH2=XLNS(KKK,2)*SBR2
      Z22=0.88D0*2.D0*RLH2
      C0BR22=-0.01D0*RHOS(KKK,2)*DSIGS(KKK,2)*Z22**4
     +              *(1.D0-2.D0/Z22*DLOG(1.D0+0.5D0*Z22))
      CDP=CDP+WINTS(KKK)*C0BR22
      DELN=1.D0-XLNS(KKK,2)
      CALL BHOLE(RHOS(KKK,2),QUADS(KKK,2),XLNS(KKK,2),X,B)
      CALL MOMS(XLNS(KKK,2),X,B,F1,F2,F3,F4,F5)
      BRDIP2(KKK)=BRDIP2(KKK)+RHOS(KKK,2)*B**2
      BRDINT=BRDINT+WINTS(KKK)*RHOS(KKK,2)*B**2
      F3=F3/RHOS(KKK,2)
      F4=F4/RHOS(KKK,2)
      COEF1=DELN/(4.D0*PI*F4)
      COEF2=DMIN1(DSIGS(KKK,2)/3.D0,COEF1)
      CNDP=CNDP-WINTS(KKK)*2.D0*PI*COEF2*RHOS(KKK,2)*F3
      ! Bc95:
      call gascor(rhos(kkk,2),0.d0,ec22ueg,dum1,dum2,dum3)
      dsigueg = cdsigueg*rhos(kkk,2)**thrd5
      ec95d = ec22ueg*dsigs(kkk,2)/dsigueg/(1.d0+cc95par*dlss**2)**2
      ec95par = ec95par + wints(kkk)*ec95d
      ! end of Bc95 for RHOS(KKK,1).LT.SMALL
C
      ELSE ! RHOS(KKK,1).GE.SMALL .AND. RHOS(KKK,2).GE.SMALL
      FLIP12=(1.D0-XLNS(KKK,1))/XLNS(KKK,2)
      FLIP21=(1.D0-XLNS(KKK,2))/XLNS(KKK,1)
      FLIP=DMIN1( FLIP12 , FLIP21 , 1.D0 )
C     SPIN 1:
      DRHOSQ=4.D0*RHOS(KKK,1)*WEIZS(KKK,1) ! a backwards route to del rho
      DRHO=DSQRT(DRHOSQ)                   ! del rho, gradient
      RHO43=RHOS(KKK,1)**THRD4
      XL=-CEX*RHO43                        ! LSDA exchange
      DLSS=DRHO/RHO43                      ! chi, "dimensionless gradient"
      ASNH=DLOG(DLSS+DSQRT(DLSS**2+1.D0))  ! arcsinh
      X88=XL-0.0042D0*DRHO*DLSS/(1.D0+6.D0*0.0042D0*DLSS*ASNH)
      EX88=EX88+WINTS(KKK)*X88
      RFG1=-0.5D0*RHOS(KKK,1)/X88
      Z11=0.88D0*2.D0*RFG1
      C8811=-0.01D0*RHOS(KKK,1)*DSIGS(KKK,1)*Z11**4*(1.D0-2.D0/Z11*
     C                                          DLOG(1.D0+0.5D0*Z11))
      EC8811=EC8811+WINTS(KKK)*C8811
      CALL XMODEL(RHOS(KKK,1),QUADS(KKK,1),1.D0,DUM,XBR1)
      if (ERRFLG) return
      EXBR=EXBR+WINTS(KKK)*XBR1
      RBR1=-0.5D0*RHOS(KKK,1)/XBR1
      Z11=0.88D0*2.D0*RBR1
      C1BR11=-0.01D0*RHOS(KKK,1)*DSIGS(KKK,1)*Z11**4
     +              *(1.D0-2.D0/Z11*DLOG(1.D0+0.5D0*Z11))
      ECBR11=ECBR11+WINTS(KKK)*C1BR11
      CALL XMODEL(RHOS(KKK,1),QUADS(KKK,1),XLNS(KKK,1),
     +                                     UXPS(KKK,1),XBRG)
      if (ERRFLG) return
      FXBR=FXBR+WINTS(KKK)*XBRG
      SBR1=-0.5D0*RHOS(KKK,1)/XBRG
      Z11=0.88D0*2.D0*SBR1
      C2BR11=-0.01D0*RHOS(KKK,1)*DSIGS(KKK,1)*Z11**4
     +              *(1.D0-2.D0/Z11*DLOG(1.D0+0.5D0*Z11))
      FCBR11=FCBR11+WINTS(KKK)*C2BR11
      RLH1=XLNS(KKK,1)*SBR1
      Z11=0.88D0*2.D0*RLH1
      C0BR11=-0.01D0*RHOS(KKK,1)*DSIGS(KKK,1)*Z11**4
     +              *(1.D0-2.D0/Z11*DLOG(1.D0+0.5D0*Z11))
      CDP=CDP+WINTS(KKK)*C0BR11
      DELN=1.D0-XLNS(KKK,1)-FLIP*XLNS(KKK,2)
      CALL BHOLE(RHOS(KKK,1),QUADS(KKK,1),XLNS(KKK,1),X,B)
      if (ERRFLG) return
      CALL MOMS(XLNS(KKK,1),X,B,F1,F2,F3,F4,F5)
      BRDIP2(KKK)=BRDIP2(KKK)+RHOS(KKK,1)*B**2
      BRDINT=BRDINT+WINTS(KKK)*RHOS(KKK,1)*B**2
      F3=F3/RHOS(KKK,1)
      F4=F4/RHOS(KKK,1)
      COEF1=DELN/(4.D0*PI*F4)
      COEF2=DMIN1(DSIGS(KKK,1)/3.D0,COEF1)
      CNDP=CNDP-WINTS(KKK)*2.D0*PI*COEF2*RHOS(KKK,1)*F3
      ! Bc95:
      call gascor(rhos(kkk,1),0.d0,ec11ueg,dum1,dum2,dum3)
      dsigueg = cdsigueg * rhos(kkk,1)**thrd5
      ec95d   = ec11ueg*dsigs(kkk,1)/dsigueg/(1.d0+cc95par*dlss**2)**2
      ec95par = ec95par + wints(kkk)*ec95d
      chi1 = dlss
      ! end Bc95 (spin 1)
C     SPIN 2:
      DRHOSQ=4.D0*RHOS(KKK,2)*WEIZS(KKK,2)
      DRHO=DSQRT(DRHOSQ)
      RHO43=RHOS(KKK,2)**THRD4
      XL=-CEX*RHO43
      DLSS=DRHO/RHO43
      ASNH=DLOG(DLSS+DSQRT(DLSS**2+1.D0))
      X88=XL-0.0042D0*DRHO*DLSS/(1.D0+6.D0*0.0042D0*DLSS*ASNH)
      EX88=EX88+WINTS(KKK)*X88
      RFG2=-0.5D0*RHOS(KKK,2)/X88
      Z22=0.88D0*2.D0*RFG2
      C8822=-0.01D0*RHOS(KKK,2)*DSIGS(KKK,2)*Z22**4*(1.D0-2.D0/Z22*
     C                                          DLOG(1.D0+0.5D0*Z22))
      EC8811=EC8811+WINTS(KKK)*C8822
      CALL XMODEL(RHOS(KKK,2),QUADS(KKK,2),1.D0,DUM,XBR1)
      if (ERRFLG) return
      EXBR=EXBR+WINTS(KKK)*XBR1
      RBR2=-0.5D0*RHOS(KKK,2)/XBR1
      Z22=0.88D0*2.D0*RBR2
      C1BR22=-0.01D0*RHOS(KKK,2)*DSIGS(KKK,2)*Z22**4
     +              *(1.D0-2.D0/Z22*DLOG(1.D0+0.5D0*Z22))
      ECBR11=ECBR11+WINTS(KKK)*C1BR22
      CALL XMODEL(RHOS(KKK,2),QUADS(KKK,2),XLNS(KKK,2),
     +                                     UXPS(KKK,2),XBRG)
      if (ERRFLG) return
      FXBR=FXBR+WINTS(KKK)*XBRG
      SBR2=-0.5D0*RHOS(KKK,2)/XBRG
      Z22=0.88D0*2.D0*SBR2
      C2BR22=-0.01D0*RHOS(KKK,2)*DSIGS(KKK,2)*Z22**4
     +              *(1.D0-2.D0/Z22*DLOG(1.D0+0.5D0*Z22))
      FCBR11=FCBR11+WINTS(KKK)*C2BR22
      RLH2=XLNS(KKK,2)*SBR2
      Z22=0.88D0*2.D0*RLH2
      C0BR22=-0.01D0*RHOS(KKK,2)*DSIGS(KKK,2)*Z22**4
     +              *(1.D0-2.D0/Z22*DLOG(1.D0+0.5D0*Z22))
      CDP=CDP+WINTS(KKK)*C0BR22
      DELN=1.D0-XLNS(KKK,2)-FLIP*XLNS(KKK,1)
      CALL BHOLE(RHOS(KKK,2),QUADS(KKK,2),XLNS(KKK,2),X,B)
      if (ERRFLG) return
      CALL MOMS(XLNS(KKK,2),X,B,F1,F2,F3,F4,F5)
      BRDIP2(KKK)=BRDIP2(KKK)+RHOS(KKK,2)*B**2
      BRDINT=BRDINT+WINTS(KKK)*RHOS(KKK,2)*B**2
      F3=F3/RHOS(KKK,2)
      F4=F4/RHOS(KKK,2)
      COEF1=DELN/(4.D0*PI*F4)
      COEF2=DMIN1(DSIGS(KKK,2)/3.D0,COEF1)
      CNDP=CNDP-WINTS(KKK)*2.D0*PI*COEF2*RHOS(KKK,2)*F3
      ! Bc95:
      call gascor(rhos(kkk,2),0.d0,ec22ueg,dum1,dum2,dum3)
      dsigueg = cdsigueg * rhos(kkk,2)**thrd5
      ec95d   = ec22ueg*dsigs(kkk,2)/dsigueg/(1.d0+cc95par*dlss**2)**2
      ec95par = ec95par + wints(kkk)*ec95d
      chi2 = dlss
      ! end Bc95 (spin 2)
C     OPP SPINS:
      Z12=0.63D0*(RFG1+RFG2)
      C8812=-0.8D0*RHOS(KKK,1)*RHOS(KKK,2)*Z12**2*(1.D0-DLOG(1.D0+Z12)
     C        /Z12)
      EC8812=EC8812+WINTS(KKK)*C8812
      Z12=0.63D0*RBR1+0.63D0*RBR2
      CBR12=-0.8D0*RHOS(KKK,1)*RHOS(KKK,2)*Z12**2
     +            *(1.D0-DLOG(1.D0+Z12)/Z12)
      ECBR12=ECBR12+WINTS(KKK)*CBR12
      Z12=0.63D0*SBR1+0.63D0*SBR2
      CBR12=-0.8D0*RHOS(KKK,1)*RHOS(KKK,2)*Z12**2
     +            *(1.D0-DLOG(1.D0+Z12)/Z12)
      FCBR12=FCBR12+WINTS(KKK)*CBR12
      Z12=0.63D0*RLH1+0.63D0*RLH2
      CD12=-0.8D0*RHOS(KKK,1)*RHOS(KKK,2)*Z12**2
     +            *(1.D0-DLOG(1.D0+Z12)/Z12)
      CD=CD+WINTS(KKK)*CD12
      CDC=CDC+WINTS(KKK)*CD12*FLIP
      CND12=CND12-0.5D0*WINTS(KKK)*RHOS(KKK,1)*FLIP*UXPS(KKK,2)
      CND21=CND21-0.5D0*WINTS(KKK)*RHOS(KKK,2)*FLIP*UXPS(KKK,1)
      ! Bc95:
      call gascor(rhos(kkk,1),rhos(kkk,2),ec12ueg,dum1,dum2,dum3)
      ec12ueg = ec12ueg - ec11ueg - ec22ueg
      ec95d = ec12ueg/(1.d0+cc95opp*(chi1**2+chi2**2))
      ec95opp = ec95opp + wints(kkk)*ec95d
      ! end Bc95 (opposite spins)
      END IF ! RHOS<SMALL
300   CONTINUE
      end if ! myRank==0
      RETURN
      END
C
C
C$PRAGMA SUN OPT=1
      SUBROUTINE TAUWZL(ISPIN,TAU,TAUPD,WEIZ,D2RHO)
C
C     POSITIVE-DEFINITE TAU (NO 1/2 FACTOR),
C     WEIZ = 0.25*(GRAD-RHO)**2/RHO,
C     AND LAPLACIAN OF SPIN-DENSITY (SPIN=ISPIN).
C
      USE CONSTANTS
      USE IOUNITS
      USE MPI
      USE OccupationNumbers
      USE SIZES
      IMPLICIT NONE
      INTEGER,        INTENT(IN) :: ISPIN
      REAL(KIND(1D0)),INTENT(IN) :: TAU(NNN)
      REAL(KIND(1D0)),INTENT(OUT):: TAUPD(NNN),WEIZ(NNN),D2RHO(NNN)
      REAL(KIND(1D0)) PSIXYZ(0:3,NNN),PSIDC(NNN)
      REAL(KIND(1D0)) RHO(NNN),DRHOX(NNN),DRHOY(NNN),DRHOZ(NNN)
      REAL(KIND(1D0)) PSI,OCC
      INTEGER I,ICORB,IODC,IMORB,KKK,IREC
C
      DO 10 I=1,NNN
      RHO(I)=0.D0
      TAUPD(I)=0.D0
      DRHOX(I)=0.D0
      DRHOY(I)=0.D0
10    DRHOZ(I)=0.D0
C
      IF(NCORBS.NE.0.and.myRank.eq.0)THEN
        REWIND(ICRDC)
        DO 100 ICORB=1,NCORBS
        READ(ICRDC)(PSIDC(I),I=1,NNN)
        CALL FDPSI(PSIDC,1,PSIXYZ,ICORB)
        DO 110 KKK=1,NNN
        PSI=PSIXYZ(0,KKK)
        RHO(KKK)=RHO(KKK)+PSI*PSI
        DRHOX(KKK)=DRHOX(KKK)+PSI*PSIXYZ(1,KKK)
        DRHOY(KKK)=DRHOY(KKK)+PSI*PSIXYZ(2,KKK)
        DRHOZ(KKK)=DRHOZ(KKK)+PSI*PSIXYZ(3,KKK)
        TAUPD(KKK)=TAUPD(KKK)
     +          +(PSIXYZ(1,KKK)**2+PSIXYZ(2,KKK)**2+PSIXYZ(3,KKK)**2)
110     CONTINUE
100     CONTINUE
      END IF ! nCOrbs>0
C
      IF(ISPIN.EQ.1)THEN
        IODC=MODCA
      ELSE
        IODC=MODCB
      END IF
C
      DO 200 IMORB=1,NMORBS
      OCC=OCCS(IMORB,ISPIN)
      if (myRank.eq.IOwner(iMOrb)) then
        iRec=MOSeq(iMOrb)
        READ(IODC,REC=iRec)(PSIDC(I),I=1,NNN)
        CALL FDPSI(PSIDC,1,PSIXYZ,0)
        DO 210 KKK=1,NNN
        PSI=PSIXYZ(0,KKK)
        RHO(KKK)=RHO(KKK)+OCC*PSI*PSI
        DRHOX(KKK)=DRHOX(KKK)+OCC*PSI*PSIXYZ(1,KKK)
        DRHOY(KKK)=DRHOY(KKK)+OCC*PSI*PSIXYZ(2,KKK)
        DRHOZ(KKK)=DRHOZ(KKK)+OCC*PSI*PSIXYZ(3,KKK)
        TAUPD(KKK)=TAUPD(KKK)+OCC*
     +           (PSIXYZ(1,KKK)**2+PSIXYZ(2,KKK)**2+PSIXYZ(3,KKK)**2)
210     CONTINUE
      end if ! myRank==IOwner
200   CONTINUE
C
C     Reduce RHO,TAUPD,DRHOX/Y/Z re-using PSIDC for the accumulator
      CALL MPI_REDUCE(RHO,PSIDC,NNN,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                0,MPI_COMM_WORLD,IERROR)
      DO 291 KKK=1,NNN
291     RHO(KKK)=PSIDC(KKK)
C
      CALL MPI_REDUCE(TAUPD,PSIDC,NNN,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                0,MPI_COMM_WORLD,IERROR)
      DO 292 KKK=1,NNN
292     TAUPD(KKK)=PSIDC(KKK)
C
      CALL MPI_REDUCE(DRHOX,PSIDC,NNN,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                0,MPI_COMM_WORLD,IERROR)
      DO 293 KKK=1,NNN
293     DRHOX(KKK)=PSIDC(KKK)
C
      CALL MPI_REDUCE(DRHOY,PSIDC,NNN,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                0,MPI_COMM_WORLD,IERROR)
      DO 294 KKK=1,NNN
294     DRHOY(KKK)=PSIDC(KKK)
C
      CALL MPI_REDUCE(DRHOZ,PSIDC,NNN,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                0,MPI_COMM_WORLD,IERROR)
      DO 295 KKK=1,NNN
295     DRHOZ(KKK)=PSIDC(KKK)
C
      DO 300 KKK=1,NNN
      D2RHO(KKK)=2.D0*(TAU(KKK)+TAUPD(KKK))
300   WEIZ(KKK)=(DRHOX(KKK)**2+DRHOY(KKK)**2+DRHOZ(KKK)**2)
     +                                               /(RHO(KKK)+EPS)
C
      RETURN
      END
C
C
C
      SUBROUTINE XLNORM(RHO,QUAD,UXPOS,XLNRM)
      USE Constants
      USE IOUNITS
      USE ERRORHANDLING
      IMPLICIT NONE
      REAL(KIND(1D0)),INTENT(IN) :: RHO,QUAD,UXPOS
      REAL(KIND(1D0)),INTENT(OUT):: XLNRM
      REAL(KIND(1D0)) RHS,X0,SHIFT,X,F,DF,X1,ALF,A
      INTEGER I
      IF(RHO.LT.SMALL)THEN
        XLNRM=1.D0
        RETURN
      END IF
      RHS=4.D0*PI/3.D0*RHO*RHO/QUAD/UXPOS
      X0=2.D0
      SHIFT=1.D0
      IF(RHS.LT.0.D0)GO TO 10
      IF(RHS.GT.0.D0)GO TO 20
10    DO 11 I=1,16
      X=X0-SHIFT
      CALL XLFUNS(X,RHS,F,DF)
      IF(F.LT.0.D0)GO TO 88
11    SHIFT=0.1D0*SHIFT
      WRITE(IOUT,1002)
      ERRFLG = .TRUE.
      RETURN
20    DO 21 I=1,16
      X=X0+SHIFT
      CALL XLFUNS(X,RHS,F,DF)
      IF(F.GT.0.D0)GO TO 88
21    SHIFT=0.1D0*SHIFT
      WRITE(IOUT,1002)
      ERRFLG = .TRUE.
      RETURN
88    CONTINUE
      DO 100 I=1,100
      CALL XLFUNS(X,RHS,F,DF)
      X1=X-F/DF
      IF(DABS(X1-X).LT.SMALL)GO TO 111
      X=X1
100   CONTINUE
      WRITE(IOUT,1001)
      ERRFLG = .TRUE.
      RETURN
111   X=X1
      ALF=DSQRT(6.D0*QUAD*X/RHO/(X-2.D0))
      A=RHO*DEXP(X)
      XLNRM=DMIN1(8.D0*PI*A/ALF**3,2.D0)
      RETURN
1001  FORMAT(' XLNORM: NEWTON ALGORITHM FAILS TO CONVERGE!')
1002  FORMAT(' XLNORM: NEWTON ALGORITHM FAILS TO INITIALIZE!')
1003  FORMAT(' RHO:',F16.8,' QUAD:',F16.8,' UXPOS:',F16.8)
      END
C
      SUBROUTINE XLFUNS(X,RHS,F,DF)
      IMPLICIT NONE
      REAL(KIND(1D0)),INTENT(IN) :: X,RHS
      REAL(KIND(1D0)),INTENT(OUT):: F,DF
      REAL(KIND(1D0)) EXPO,BOT
      EXPO=DEXP(X)
      BOT=(X-2.D0)*(EXPO-1.D0-0.5D0*X)
      F = X*X/BOT - RHS
      DF=4.D0*X-(4.D0*X-3.D0*X*X+X**3)*EXPO
      DF=DF/BOT**2
      RETURN
      END
C
C
C
      SUBROUTINE XMODEL(RHO,QUAD,HNORM,UXPOS,EX)
      USE Constants
      USE IOUNITS
      USE ERRORHANDLING
      IMPLICIT NONE
      REAL(KIND(1D0)),INTENT(IN) :: RHO,QUAD,HNORM
      REAL(KIND(1D0)),INTENT(OUT):: UXPOS,EX
      REAL(KIND(1D0)) RHS,X0,SHIFT,X,F,DF,X1,EXPO,PREFAC,ALF
      INTEGER I
      IF(RHO.LT.SMALL)THEN
        EX=0.D0
        UXPOS=SMALL
        RETURN
      END IF
      RHS=THRD2*(PI*RHO/HNORM)**THRD2*RHO/QUAD
      X0=2.D0
      SHIFT=1.D0
      IF(RHS.LT.0.D0)GO TO 10
      IF(RHS.GT.0.D0)GO TO 20
10    DO 11 I=1,16
      X=X0-SHIFT
      CALL XFUNCS(X,RHS,F,DF)
      IF(F.LT.0.D0)GO TO 88
11    SHIFT=0.1D0*SHIFT
      WRITE(IOUT,1002)
      WRITE(IOUT,1003)RHO,QUAD,HNORM
      ERRFLG = .TRUE.
      RETURN
20    DO 21 I=1,16
      X=X0+SHIFT
      CALL XFUNCS(X,RHS,F,DF)
      IF(F.GT.0.D0)GO TO 88
21    SHIFT=0.1D0*SHIFT
      WRITE(IOUT,1002)
      WRITE(IOUT,1003)RHO,QUAD,HNORM
      ERRFLG = .TRUE.
      RETURN
88    CONTINUE
      DO 100 I=1,100
      CALL XFUNCS(X,RHS,F,DF)
      X1=X-F/DF
      IF(DABS(X1-X).LT.SMALL)GO TO 111
      X=X1
100   CONTINUE
      WRITE(IOUT,1001)
      WRITE(IOUT,1003)RHO,QUAD,HNORM
      ERRFLG = .TRUE.
      RETURN
111   X=X1
      EXPO=DEXP(-X)
      PREFAC=RHO/EXPO
      ALF=(8.D0*PI*PREFAC/HNORM)**THRD
      UXPOS=4.D0*PI*PREFAC/ALF/ALF/X*(2.D0-2.D0*EXPO-X*EXPO)
      EX=-0.5D0*RHO*UXPOS
      RETURN
1001  FORMAT(' XMODEL: NEWTON ALGORITHM FAILS TO CONVERGE!')
1002  FORMAT(' XMODEL: NEWTON ALGORITHM FAILS TO INITIALIZE!')
1003  FORMAT(' RHO:',F16.8,' QUAD:',F16.8,' HNORM:',F16.8)
      END
C
C
C
      SUBROUTINE XFUNCS(X,RHS,F,DF)
      IMPLICIT NONE
      REAL(KIND(1D0)) X,RHS,F,DF,EXPO23
      EXPO23=DEXP(-2.D0/3.D0*X)
      F = X*EXPO23/(X-2.D0) - RHS
      DF=2.D0/3.D0*(2.D0*X-X**2-3.D0)/(X-2.D0)**2*EXPO23
      RETURN
      END
C
C
C
      SUBROUTINE BHOLE(RHO,QUAD,HNORM,X,B)
      USE Constants
      USE IOUNITS
      USE ERRORHANDLING
      IMPLICIT NONE
      ! Arguments
      REAL(KIND(1D0)),INTENT(IN) :: RHO,QUAD,HNORM
      REAL(KIND(1D0)),INTENT(OUT):: X,B
      ! Locals
      REAL(KIND(1D0)) RHS,X0,SHIFT,F,DF,X1,EXPO,PREFAC,ALF
      INTEGER I
      RHS=THRD2*(PI*RHO/HNORM)**THRD2*RHO/QUAD
      X0=2.D0
      SHIFT=1.D0
      IF(RHS.LT.0.D0)GO TO 10
      IF(RHS.GT.0.D0)GO TO 20
10    DO 11 I=1,16
      X=X0-SHIFT
      CALL XFUNCS(X,RHS,F,DF)
      IF(F.LT.0.D0)GO TO 88
11    SHIFT=0.1D0*SHIFT
      goto 1002
20    DO 21 I=1,16
      X=X0+SHIFT
      CALL XFUNCS(X,RHS,F,DF)
      IF(F.GT.0.D0)GO TO 88
21    SHIFT=0.1D0*SHIFT
      goto 1002
88    CONTINUE
      DO 100 I=1,100
      CALL XFUNCS(X,RHS,F,DF)
      X1=X-F/DF
      IF(DABS(X1-X).LT.SMALL)GO TO 111
      X=X1
100   CONTINUE
      goto 1001
111   X=X1
      EXPO=DEXP(-X)
      PREFAC=RHO/EXPO
      ALF=(8.D0*PI*PREFAC/HNORM)**THRD
      B=X/ALF
      RETURN
1001  WRITE(IOUT,*)' BHOLE: NEWTON ALGORITHM FAILS TO CONVERGE!'
      GOTO 9999
1002  WRITE(IOUT,*)' BHOLE: NEWTON ALGORITHM FAILS TO INITIALIZE!'
 9999 CONTINUE
      ERRFLG = .TRUE.
      RETURN
      END
C
C
C
      SUBROUTINE MOMS(HNORM,X,B,F1,F2,F3,F4,F5)
      USE Constants, only : PI
      IMPLICIT NONE
      ! Arguments
      REAL(KIND(1D0)),INTENT(IN) :: HNORM,X,B
      REAL(KIND(1D0)),INTENT(OUT):: F1,F2,F3,F4,F5
      ! Locals
      REAL(KIND(1D0)) PIREC,BREC,BSQR,BCUB,EXPO,X1,X2,X3,X4
      PIREC=1.0D0/PI
      BREC=1.0D0/B
      BSQR=B**2
      BCUB=B*BSQR
      EXPO=-DEXP(-X)
      X1=1.0D0/X
      X2=1.0D0/X**2
      X3=1.0D0/X**3
      X4=1.0D0/X**4
      F1=EXPO*(0.125D0*X+0.25D0)
      F1=BREC*PIREC*(0.25D0+F1)
      F2=0.25D0*PIREC
      F3=0.25D0+X2
      F3=F3+EXPO*(0.25D0*X1+X2)
      F3=B*PIREC*F3
      F4=0.25D0+3.0D0*X2
      F4=BSQR*PIREC*F4
      F5=0.25D0+6.0D0*X2+18.0D0*X4
      F5=F5+3.0D0*EXPO*(X3+6.0D0*X4)
      F5=BCUB*PIREC*F5
      F1=HNORM*F1
      F2=HNORM*F2
      F3=HNORM*F3
      F4=HNORM*F4
      F5=HNORM*F5
      RETURN
      END
