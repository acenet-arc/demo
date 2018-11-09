C $Header: /home/ross/cvsroot/parnum/parnum.f,v 1.32 2005/04/26 17:03:00 ross Exp $
C$PRAGMA SUN OPT=5
      PROGRAM NUMOL
C**********************************************************************C
C                                                                      C
C     A NUMERICAL DENSITY-FUNCTIONAL MOLECULAR STRUCTURE PROGRAM.      C
C                                                                      C
C  WRITTEN BY A.D. BECKE AND R.M. DICKSON, QUEEN'S UNIVERSITY, 1988.   C
C     THIS VERSION EXTENSIVELY REVISED BY A.D. BECKE, 1988-1996.       C
C                                                                      C
C          MPI parallelized version by Ross Dickson, 2004.             C
C          This is a custom version for Sun Fortran 90 2.0.            C
C                                                                      C
C**********************************************************************C
C                                                                      C
C                 COPYRIGHT (C) 1996, BY A.D. BECKE.                   C
C                                                                      C
C    DISTRIBUTION OF THIS PROGRAM, OR ANY SEGMENT OF THIS PROGRAM,     C
C                            IS PROHIBITED.                            C
C        THE COPYRIGHT OWNER RESERVES ALL DISTRIBUTION RIGHTS.         C
C                                                                      C
C              NO WARRANTIES ARE EXPRESSED OR IMPLIED                  C
C                 THAT THE PROGRAM IS FREE OF ERROR.                   C
C             THE USER SHALL EMPLOY AT HIS/HER OWN RISK,               C
C         AND WE DISCLAIM ALL LIABILITY FOR THE USE OR MISUSE          C
C               OF THIS PROGRAM, OR ANY RESULTS THEREOF.               C
C                                                                      C
C**********************************************************************C
C
C   EXTERNAL SUBROUTINES FROM LIBRARIES  LINPACK AND EISPACK.
C
      USE CONSTANTS
      USE CORES
      USE ERRORHANDLING
      USE GeometryOptimization
      USE GRID
      USE IOUNITS
      USE KEYWORDS
      USE MPI
      USE NuclearFramework
      USE NuclearMeshes
      USE NondefaultOccupancies
      USE OccupationNumbers
      USE Potentials
      USE PROMOLECULE
      USE SCFConvergence
      USE SIZES
      USE SphericalHarmonics ! not used this unit, but needs global scope
      IMPLICIT NONE
C***********************************************************************
C     Historical note:  Dynamic arrays have largely replaced COMMON/WORK/
C
C     Dynamic arrays might actually increase memory usage
C     since Absoft Fortran needs static allocation to co-operate with
C     MPI.  If this turns out to be a problem we might have to re-instate
C     a WORKSPACE module which shares space between various routines via
C     renaming.
C
C     TODO:  Does explicit allocation/deallocation (as opposed to automatic
C     arrays) dodge the static usage when using the -s flag in Absoft
C     Fortran?  If so that provides an alternative approach to space saving.
C
C***********************************************************************
C extra comment
      LOGICAL GCONVG
      REAL(KIND(1D0)) DIPOLE(3)
      REAL(KIND(1D0)), ALLOCATABLE :: RHO(:),RHOS(:,:),TAUS(:,:)
      REAL(KIND(1D0)), ALLOCATABLE :: HIRSH(:)
      REAL(KIND(1D0)), ALLOCATABLE :: EXDENS(:,:)
      REAL(KIND(1D0)), ALLOCATABLE :: BRDIP2(:),XXDIP2(:)
      REAL(KIND(1D0)) TCORE,ENUC,R,CHGTST,VX1,VX2,ECOR,VC1,VC2,UCOR
      REAL(KIND(1D0)) EKIN,EXL,ECL,UCL,EPOT,RHO1,RHO2,VNEW1,VNEW2,ELDA
      REAL(KIND(1D0)) VXNEW1,VXNEW2
      REAL(KIND(1D0)) EXGGA,ECGGA,EGGA,FXSUM,FYSUM,FZSUM,VTNUC,ECLASS
      REAL(KIND(1D0)) CHRG,XCHRG,YCHRG,ZCHRG,SPINN,FUN,XMEAN,YMEAN,ZMEAN
      REAL(KIND(1D0)) EXX,EX88,EC8812,EC8811,EXBR,ECBR12,ECBR11,FXBR
      REAL(KIND(1D0)) EC95,EC95OPP,EC95PAR
      REAL(KIND(1D0)) FCBR12,FCBR11,CND12,CND21,CD,CDC,CNDP,CDP,CDPC
      REAL(KIND(1D0)) BRDINT,XXDINT,ATMXXDIP2,ATMBRDIP2,FUNXX,FUNBR
      REAL(KIND(1D0)) ETOTAL

      INTEGER NELS,NV0S,ITGEOM,INUC,JNUC,KKK,IT,I
C
      call MPI_Init(iError)
      call MPI_Comm_Rank(MPI_COMM_WORLD,myRank,iError)
      call MPI_Comm_Size(MPI_COMM_WORLD,nProcs,iError)
      call MPI_Get_processor_name(myHost,lenHost,iError)
      ERRFLG = .FALSE.
      write(IOUT,*) 'This is MPI process ',myRank,' of ',nProcs,
     &              ' on ',trim(myHost)
      call MPI_Barrier(MPI_COMM_WORLD,iError) ! just to keep the output tidy
      OPEN(UNIT=IATMS,FILE='atomic.data',STATUS='OLD',ERR=9998)
C
      CALL InitializeConstants
C
      IF (myRank.EQ.0) THEN
         CALL INPUT(NELS,NV0S) ! calls ALLOCATE_FRAMEWORK on root node
      END IF
      call flush(IOUT)  ! get Sun status out to the user
      CALL BCAST_INPUT_SCALARS(NELS,NV0S)
      IF (myRank.NE.0) CALL ALLOCATE_FRAMEWORK()
      CALL BCAST_INPUT_ARRAYS()
C
      OPEN(UNIT=IHRSH,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(UNIT=ICOR0,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(UNIT=ICOR2,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(IBAS0,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      OPEN(IBAS2,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      OPEN(IBSDC,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      OPEN(MOSA0,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      OPEN(MOSA2,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      OPEN(MOSB0,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      OPEN(MOSB2,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      OPEN(UNIT=ICRDC,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(MODCA,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      OPEN(MODCB,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      IF (IPOT.EQ.XX_POT) THEN
      OPEN(UNIT=30,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(UNIT=40,STATUS='SCRATCH',FORM='UNFORMATTED')
      open(31   ,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      open(41   ,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      open(42   ,STATUS='SCRATCH',ACCESS='DIRECT',RECL=NNN*8)
      END IF
C
C     Allocate major arrays
C     TODO:  Is there a way to trap out-of-memory failures here?
C     There doesn't appear to be.  Absoft manual makes no mention of one,
C     and a test program using 'allocated' inquiry function doesn't do it.
      allocate(XMESH(NNN))
      allocate(YMESH(NNN))
      allocate(ZMESH(NNN))
      allocate(WINTS(NNN))
      allocate(WNUC(NNN))
      allocate(V(NNN,2))
      allocate(VNUC(NNN))
      allocate(VEL(NNN))
      allocate(VXC(NNN,2))
      allocate(RHOCOR(NNN))
      allocate(TAUCOR(NNN))
      allocate(RHO0(NNN))
      allocate(VEL0(NNN))
      allocate(RHO(NNN))
      allocate(RHOS(NNN,2))
      allocate(TAUS(NNN,2))
C
      ITGEOM=1
      GCONVG=.FALSE.
888   CALL MESH0
      if (ERRFLG) goto 9999
      IF (I_PRINT_MESH.GT.0) THEN
         ! Just print the integration points and quit.
         WRITE(IOUT,'(/I7,'' molecular mesh points:'')') NNN
         DO KKK=1,NNN
            WRITE(IOUT,*) XMESH(KKK),YMESH(KKK),ZMESH(KKK)
         END DO
         GOTO 9999
      END IF
C
      CALL ATOMIN(TCORE,RHO,RHOS,VEL,VNUC,FORCE0)
      if (ERRFLG) goto 9999
C     ON RETURN:
C     TCORE   =  CORE KINETIC ENERGY (ZERO IF CORES UNFROZEN).
C     RHO/VEL =  PROMOLECULE DENSITY AND COULOMB POTENTIAL.
C     RHOS    =  USER-MODIFIABLE STARTING SPIN-DENSITIES.
C     VNUC    =  NUCLEAR POTENTIAL.
C     FORCE0   =  REFERENCE FORCES ARISING FROM PROMOLECULE DENSITY.
C
      IF(IOPTG.EQ.-1)THEN
         CALL NEWGEO(ITGEOM,GCONVG)
         if (ERRFLG.OR.GCONVG) goto 9999
      END IF
C
      IF(ITGEOM.EQ.1)THEN
         CALL PopulateOrbitals(NELS,NV0S)
         if (ERRFLG) GOTO 9999
      END IF ! ITGEOM==1
C
      ENUC=0.D0
      IF(NUCS.NE.1)THEN
      DO 160 INUC=1,NUCS-1
      DO 161 JNUC=INUC+1,NUCS
      R=DSQRT((XNUC(JNUC)-XNUC(INUC))**2
     +       +(YNUC(JNUC)-YNUC(INUC))**2
     +       +(ZNUC(JNUC)-ZNUC(INUC))**2)
161   ENUC=ENUC+NUCZ(INUC)*NUCZ(JNUC)/R
160   CONTINUE
      END IF
C
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'NUMBER OF MESH POINTS =',NNN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'INTEGRATION TEST ON PROMOLECULE DENSITY:'
      END IF
C
      CHGTST=0.D0
      DO 200 KKK=1,NNN
      CHGTST=CHGTST+WINTS(KKK)*RHO(KKK)
C     REPLACE PROMOLECULE DENSITY RHO
C     WITH USER-MODIFIABLE STARTING DENSITY,
200   RHO(KKK) = RHOS(KKK,1)+RHOS(KKK,2)
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)'TOTAL CHARGE  =    ',CHGTST
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'VALENCE MOLECULAR ORBITALS  =',NMORBS
      WRITE(IOUT,*)'LCAO VALENCE BASIS FUNCTIONS  =',NBASIS
      END IF
      call flush(IOUT)   ! for Sun again
C
C     INITIALIZE SPIN-POTENTIALS:
C
      CALL POISS(RHO,VEL,0)
C
      IF (IPOT.EQ.LDA_POT.OR.IPOT.EQ.GGA_POT) THEN
        DO 210 KKK=1,NNN
        VX1=-CVX*(RHOS(KKK,1)+EPS)**THRD
        VX2=-CVX*(RHOS(KKK,2)+EPS)**THRD
        CALL GASCOR(RHOS(KKK,1),RHOS(KKK,2),ECOR,VC1,VC2,UCOR)
        V(KKK,1)=VNUC(KKK)+VEL(KKK)+VX1+VC1
210     V(KKK,2)=VNUC(KKK)+VEL(KKK)+VX2+VC2
      ELSE                ! exact-exchange
        DO 211 KKK=1,NNN
        V(KKK,1)=VNUC(KKK)+VEL(KKK)
        V(KKK,2)=VNUC(KKK)+VEL(KKK)
        VXC(KKK,1)=-CVX*(RHOS(KKK,1)+EPS)**THRD
211     VXC(KKK,2)=-CVX*(RHOS(KKK,2)+EPS)**THRD
      END IF              ! ipot
C
      call MPI_Barrier( MPI_COMM_WORLD, iError )
      CALL LCAO0
      call flush(IOUT)    ! for Sun again
C
C
C     SCF ITERATIONS:
C
      DO 300 IT=1,NITS
C
      IF (IPOT.EQ.LDA_POT.OR.IPOT.EQ.GGA_POT) THEN

      IF(IT.EQ.1)THEN
        IF(NV0S.EQ.1)CALL LCAO1(IT,V,EKIN,RHO,RHOS,TAUS)
        IF(NV0S.EQ.2)CALL LCAO2(IT,V,EKIN,RHO,RHOS,TAUS)
      ELSE IF(IT.GT.1.AND.IT.LE.NITS0)THEN
        IF(NSPINS.EQ.1)CALL LCAO1(IT,V,EKIN,RHO,RHOS,TAUS)
        IF(NSPINS.EQ.2)CALL LCAO2(IT,V,EKIN,RHO,RHOS,TAUS)
      ELSE
        CALL SCFNUM(IT,V,EKIN,RHO,RHOS,TAUS)
      END IF

      ELSE  ! exact-exchange potential

      IF(IT.EQ.1)THEN
        IF(NV0S.EQ.1)CALL LCAO1X(IT,EKIN,RHO,RHOS,TAUS)
        IF(NV0S.EQ.2)CALL LCAO2X(IT,EKIN,RHO,RHOS,TAUS)
      ELSE IF(IT.GT.1.AND.IT.LE.NITS0)THEN
        IF(NSPINS.EQ.1)CALL LCAO1X(IT,EKIN,RHO,RHOS,TAUS)
        IF(NSPINS.EQ.2)CALL LCAO2X(IT,EKIN,RHO,RHOS,TAUS)
      ELSE
        CALL SCFNUMX(IT,EKIN,RHO,RHOS,TAUS)
      END IF

      END IF ! ipot
      if (ERRFLG) goto 9999

C     At this point, EIGS and RHO have been globalized.
C     We need to compute the new V and globalize it.
C     Currently, everybody computes the new V.
C     TODO:  Parallelize or localize the computation.
C     TODO:  Parallelize POISS
      CALL POISS(RHO,VEL,0)
C
      IF (IPOT.EQ.LDA_POT) THEN
      EXL=0.D0
      ECL=0.D0
      UCL=0.D0
      EPOT=0.D0
      DO 310 KKK=1,NNN
      RHO1=RHOS(KKK,1)
      RHO2=RHOS(KKK,2)
      VX1=-CVX*(RHO1+EPS)**THRD
      VX2=-CVX*(RHO2+EPS)**THRD
      CALL GASCOR(RHO1,RHO2,ECOR,VC1,VC2,UCOR)
      VNEW1=VNUC(KKK)+VEL(KKK)+VX1+VC1
      VNEW2=VNUC(KKK)+VEL(KKK)+VX2+VC2
      V(KKK,1)  =  FMIX*VNEW1 + (1.D0-FMIX)*V(KKK,1)
      V(KKK,2)  =  FMIX*VNEW2 + (1.D0-FMIX)*V(KKK,2)
      ECL=ECL+WINTS(KKK)*ECOR
      UCL=UCL+WINTS(KKK)*UCOR
      EPOT=EPOT+WINTS(KKK)*(RHO(KKK)*(VNUC(KKK)+0.5D0*VEL(KKK)))
310   EXL=EXL-CEX*WINTS(KKK)*((RHO1+EPS)**THRD4+(RHO2+EPS)**THRD4)
      ELDA=EKIN+EPOT+EXL+ECL+ENUC+TCORE
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,811)ELDA
 811  FORMAT(/' TOTAL LDA ENERGY  =   ',F23.6)
      END IF ! myrank==0

      ELSE IF (IPOT.EQ.GGA_POT) THEN
      CALL XCGGA(RHOS,EXL,ECL,UCL,EXGGA,ECGGA,VXC)
      EPOT=0.D0
      DO 320 KKK=1,NNN
      VNEW1=VNUC(KKK)+VEL(KKK)+VXC(KKK,1)
      VNEW2=VNUC(KKK)+VEL(KKK)+VXC(KKK,2)
      V(KKK,1)  =  FMIX*VNEW1 + (1.D0-FMIX)*V(KKK,1)
      V(KKK,2)  =  FMIX*VNEW2 + (1.D0-FMIX)*V(KKK,2)
320   EPOT=EPOT+WINTS(KKK)*(RHO(KKK)*(VNUC(KKK)+0.5D0*VEL(KKK)))
      EGGA=EKIN+EPOT+EXGGA+ECGGA+ENUC+TCORE
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,812)EGGA
      END IF ! myrank==0
 812  FORMAT(/' TOTAL GGA ENERGY  =   ',F23.6)

      ELSE IF (IPOT.EQ.XX_POT) THEN
      EXL=0.D0
      EPOT=0.D0
      DO 340 KKK=1,NNN
      RHO1=RHOS(KKK,1)
      RHO2=RHOS(KKK,2)
      VNEW1=VNUC(KKK)+VEL(KKK)
      VNEW2=VNUC(KKK)+VEL(KKK)
      V(KKK,1)=FMIX*VNEW1+(1.D0-FMIX)*V(KKK,1)
      V(KKK,2)=FMIX*VNEW2+(1.D0-FMIX)*V(KKK,2)
      VXNEW1=-CVX*(RHO1+EPS)**THRD
      VXNEW2=-CVX*(RHO2+EPS)**THRD
      VXC(KKK,1)=FMIX*VXNEW1+(1.D0-FMIX)*VXC(KKK,1)
      VXC(KKK,2)=FMIX*VXNEW2+(1.D0-FMIX)*VXC(KKK,2)
      EPOT=EPOT+WINTS(KKK)*(RHO(KKK)*(VNUC(KKK)+0.5D0*VEL(KKK)))
340   EXL=EXL-CEX*WINTS(KKK)*((RHO1+EPS)**THRD4+(RHO2+EPS)**THRD4)
      ELDA=EKIN+EPOT+EXL+ENUC
C TODO:  Replace this LDA energy with a proper SCF energy
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,813)ELDA
 813  FORMAT(/' LDA ENERGY OF SCF DENSITY =  ',F23.6)
      END IF ! myrank==0

      END IF ! ipot
C
      IF(NUCS.NE.1)THEN
      CALL FORCES(RHO,V,FORCE0,FORCE,DIPOLE)
      IF (myRank.EQ.0) THEN
      IF((MUTE.EQ.0.AND.IT.GT.NITS0).OR.IT.EQ.NITS)THEN
      IF(NETCHG.EQ.0)THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'DIPOLE MOMENT (A.U.):'
      WRITE(IOUT,'(5X,3E15.4)')(DIPOLE(I),I=1,3)
      END IF ! netchg==0
      IF(IT.GE.4)THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'HELLMANN-FEYNMAN FORCES'
      WRITE(IOUT,*)'AVERAGED OVER 4 PREVIOUS ITERATIONS:'
      END IF ! it>=4
      END IF ! it>nits0
      END IF ! myrank==0
      DO 330 INUC=1,NUCS
      DO 331 I=1,3
      FSTOR4(INUC,1,5-I)=FSTOR4(INUC,1,4-I)
      FSTOR4(INUC,2,5-I)=FSTOR4(INUC,2,4-I)
331   FSTOR4(INUC,3,5-I)=FSTOR4(INUC,3,4-I)
      FSTOR4(INUC,1,1)=FORCE(INUC,1)
      FSTOR4(INUC,2,1)=FORCE(INUC,2)
      FSTOR4(INUC,3,1)=FORCE(INUC,3)
      IF(IT.LT.4)GO TO 330
      FXSUM=0.D0
      FYSUM=0.D0
      FZSUM=0.D0
      DO 332 I=1,4
      FXSUM=FXSUM+FSTOR4(INUC,1,I)
      FYSUM=FYSUM+FSTOR4(INUC,2,I)
332   FZSUM=FZSUM+FSTOR4(INUC,3,I)
      FORCE(INUC,1)=FXSUM/4.D0
      FORCE(INUC,2)=FYSUM/4.D0
      FORCE(INUC,3)=FZSUM/4.D0
      IF (myRank.EQ.0) THEN
      IF((MUTE.EQ.0.AND.IT.GT.NITS0).OR.IT.EQ.NITS) THEN
      WRITE(IOUT,'(1X,I4,3E15.4)')INUC,(FORCE(INUC,I),I=1,3)
      END IF ! mute
      END IF ! myrank
330   CONTINUE
      END IF ! nucs>1
C
      call flush(IOUT) ! get up-to-date output from Sun
300   CONTINUE  ! next SCF iteration
C
      VTNUC=0.D0
      IF(NUCS.NE.1)THEN
      CALL FCLEAN(FORCE)
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'HELLMANN-FEYNMAN FORCES'
      WRITE(IOUT,*)'CORRECTED FOR TRANS.ROT. INVARIANCE:'
      END IF
      DO 400 INUC=1,NUCS
C*****VIRIAL THEOREM STUFF
      VTNUC=VTNUC+0.5D0*XNUC(INUC)*FORCE(INUC,1)
      VTNUC=VTNUC+0.5D0*YNUC(INUC)*FORCE(INUC,2)
      VTNUC=VTNUC+0.5D0*ZNUC(INUC)*FORCE(INUC,3)
C*****END VIRIAL THEOREM STUFF
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,'(1X,I4,3E15.4)')INUC,(FORCE(INUC,I),I=1,3)
      END IF
400   CONTINUE
      END IF
C
      IF(IOPTG.EQ.0)THEN
C
      IF(NUCS.GT.1)THEN
C      CALL VIBRAN(TFDZPE)
C      if (ERRFLG) goto 9999
C
C     HIRSHFELD ATOMIC PROPERTIES:
C     (F.L. HIRSHFELD, THEORET.CHIM.ACTA 44, 129 (1977))
C
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'ATOMIC CHARGES: MEAN X,Y,Z COORDINATES: NET SPINS:'
      REWIND(IHRSH)
      allocate(HIRSH(NNN))
      DO 410 INUC=1,NUCS
      READ(IHRSH,ERR=9997)(HIRSH(I),I=1,NNN)
      CHRG=0.D0
      XCHRG=0.D0
      YCHRG=0.D0
      ZCHRG=0.D0
      SPINN=0.D0
      DO 411 KKK=1,NNN
      FUN=WINTS(KKK)*HIRSH(KKK)*RHO(KKK)
      CHRG=CHRG+FUN
      XCHRG=XCHRG+XMESH(KKK)*FUN
      YCHRG=YCHRG+YMESH(KKK)*FUN
      ZCHRG=ZCHRG+ZMESH(KKK)*FUN
      SPINN=SPINN+WINTS(KKK)*HIRSH(KKK)*(RHOS(KKK,1)-RHOS(KKK,2))
411   CONTINUE
      XMEAN=XCHRG/CHRG/DistanceScale
      YMEAN=YCHRG/CHRG/DistanceScale
      ZMEAN=ZCHRG/CHRG/DistanceScale
      WRITE(IOUT,82)INUC,ZSYMB(INUC),CHRG,XMEAN,YMEAN,ZMEAN,SPINN
  82  FORMAT(1X,I3,3X,A2,F11.4,'  : ',3F11.4,'  : ',F11.4)
410   CONTINUE
      END IF ! myrank==0
      END IF ! nucs>1
C
      allocate(EXDENS(NNN,2))
      CALL EXACTX(EXX,EXDENS)
      allocate(BRDIP2(NNN))
      if (ERRFLG) goto 9999
      CALL XCMODS(RHO,RHOS,TAUS,EXDENS,
     +            EX88,EC8812,EC8811,
     +            EXBR,ECBR12,ECBR11,
     +            FXBR,FCBR12,FCBR11,
     +            EC95OPP,EC95PAR,
     +            CND12,CND21,
     +            CD,CDC,
     +            CNDP,CDP,CDPC,
     +            BRDINT,BRDIP2)
      if (ERRFLG) goto 9999

      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'POST-SCF ENERGIES:'
      WRITE(IOUT,*)'---------------------------------------------------'
      ECLASS=EKIN+EPOT+ENUC+TCORE
      WRITE(IOUT,83)'CLASSICAL',                      ECLASS
      WRITE(IOUT,83)'EXACT EXCHANGE',                 EXX
      WRITE(IOUT,83)'BECKE (1988) EXCHANGE',          EX88
      WRITE(IOUT,83)'TOTAL HF (ECLASS+EXX)',          ECLASS+EXX
      WRITE(IOUT,83)'NON-DYN CORRELATION (OPP,2004)', CND12+CND21
      WRITE(IOUT,83)'NON-DYN CORRELATION (PAR,2004)', CNDP
      WRITE(IOUT,83)'DYNAMICAL CORR (OPP,2004)',      ECBR12
      WRITE(IOUT,83)'DYNAMICAL CORR (PAR,2004)',      ECBR11
      WRITE(IOUT,83)'CORRELATION (OPP,1995) ???',      EC95OPP
      WRITE(IOUT,83)'CORRELATION (PAR,1995) ???',      EC95PAR
      ETOTAL = ECLASS+0.42*EXX+0.58*EX88+EC95OPP+EC95PAR
      WRITE(IOUT,83)'BB1K TOTAL ENERGY ???', ETOTAL
      WRITE(IOUT,83)'HF+ECBR ENERGY',    ECLASS+EXX+ECBR12+ECBR11
      IF (IPOT.NE.XX_POT) THEN ! post-LDA parameters
         ETOTAL=ECLASS+EXX+0.514D0*(CND12+CND21)+0.651D0*CNDP
     +                        +1.075D0*ECBR12+1.113D0*ECBR11
      ELSE ! post-HF parameters
         ETOTAL=ECLASS+EXX+0.487D0*(CND12+CND21)+0.518D0*CNDP
     +                        +1.117D0*ECBR12+0.991D0*ECBR11
      END IF
      WRITE(IOUT,83)'TOTAL (XC MODEL OF SEPT 2004)',ETOTAL
  83  FORMAT(1X,A32,F18.6)
      WRITE(IOUT,*)' '
      END IF ! myrank==0
C
C     Model C6.
C
C     Break <d**2> into atomic contributions using Hirshfeld partitioning.
C     Added by ERJ, Dec 17, 2004.
C
      allocate(XXDIP2(NNN))
      CALL DISPER(RHOS,XXDINT,XXDIP2)

      if (myRank.eq.0) then
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'INTEGRATED SQUARED XX-HOLE DIPOLE MOMENT = ',XXDINT
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'INTEGRATED SQUARED BR-HOLE DIPOLE MOMENT = ',BRDINT
      WRITE(IOUT,*)' '
      IF (NUCS.GT.1) THEN
      WRITE(IOUT,*)'Atomic contributions to integrated ',
     +             'squared dipole moment: '
      WRITE(IOUT,*)'               XX <d**2>         BR <d**2> '
      REWIND(IHRSH)
      DO 415 INUC=1,NUCS
        READ(IHRSH,ERR=9997)(HIRSH(I),I=1,NNN)
        ATMXXDIP2=0.0
        ATMBRDIP2=0.0
        DO 416 KKK=1,NNN
          FUNXX=WINTS(KKK)*HIRSH(KKK)*XXDIP2(KKK)
          FUNBR=WINTS(KKK)*HIRSH(KKK)*BRDIP2(KKK)
          ATMXXDIP2=ATMXXDIP2+FUNXX
          ATMBRDIP2=ATMBRDIP2+FUNBR
416     CONTINUE
        WRITE(IOUT,4161)INUC,ZSYMB(INUC),ATMXXDIP2,ATMBRDIP2
4161    FORMAT(I4,A5,2F18.9)
415   CONTINUE
      WRITE(IOUT,*)' '
      deallocate(HIRSH)
      END IF ! nucs>1
      end if ! myRank==0
      deallocate(XXDIP2)
      deallocate(BRDIP2)
C
      GOTO 9999
      END IF
C
      if (myrank.eq.0) then
         CALL NEWGEO(ITGEOM,GCONVG)
         if (ERRFLG) goto 9999
      end if
      CALL MPI_BCAST( GCONVG, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, iError )
C     TODO Do post-SCF calculations on converged geometry.
      if (GCONVG) goto 9999
      CALL MPI_BCAST( XNUC, NUCS, MPI_DOUBLE_PRECISION, 0,
     &                MPI_COMM_WORLD, iError )
      CALL MPI_BCAST( YNUC, NUCS, MPI_DOUBLE_PRECISION, 0,
     &                MPI_COMM_WORLD, iError )
      CALL MPI_BCAST( ZNUC, NUCS, MPI_DOUBLE_PRECISION, 0,
     &                MPI_COMM_WORLD, iError )
      ITGEOM=ITGEOM+1
      GO TO 888
C
 9997 CONTINUE
      write(IOUT,*) '*** I/O error on IHRSH. ***'
      ERRFLG = .TRUE.
      goto 9999
C
 9998 CONTINUE
      write(IOUT,*) '*** Cannot find atomic.data file. ***'
      goto 9999
C
 9999 CONTINUE
      if (ERRFLG) then
         write(IOUT,*) '*** Error termination for process',myRank,' ***'
      end if
      call MPI_Finalize(iError)
      END
C
C
C
      SUBROUTINE MESH0
C
C     MESH POINTS, NUCLEAR WEIGHTS, AND INTEGRATION WEIGHTS.
C     TODO:  Incorporate Stratmann/Scuseria alternatives
C            per routine MESH (see archived code)
C
      USE ERRORHANDLING
      USE GRID
      USE IOUNITS
      USE KEYWORDS
      USE MPI, ONLY : myRank
      USE NuclearFramework
      USE NuclearMeshes
      USE SecondDerivativeMatrices
      USE SIZES
      USE SphericalHarmonics, ONLY : MYLM
      IMPLICIT NONE
      REAL(KIND(1D0)) RR(NUCS,NUCS),CUTOFF(NUCS,NUCS)
      REAL(KIND(1D0)) R,R1,R2,X,Y,Z,HYPR,VP0,VPSUM,VPI
      INTEGER INUC,JNUC,NUC,KKK,IR,IL,NL_SUPPLIED,LMAX_SUPPLIED
C
      if (.not.allocated(RMID)) then
      allocate(RMID(NUCS))
      allocate(RADS (MR,NUCS))
      allocate(RQ   (MR,NUCS))
      allocate(RQQ  (MR,NUCS))
      allocate(WRADS(MR,NUCS))
      allocate(XANGS(ML,NUCS))
      allocate(YANGS(ML,NUCS))
      allocate(ZANGS(ML,NUCS))
      allocate(WANGS(ML,NUCS))
      allocate(KKKTAB(NUCS))
      allocate(D2R   (MR,7,NUCS))
      allocate(TRID2R(MR,3,NUCS))
      end if
C
      IF(NUCS.EQ.1)THEN
      CALL MESH1
      RETURN
      END IF
C
      DO 10 INUC=1,NUCS
      DO 10 JNUC=1,NUCS
10    RR(INUC,JNUC)=DSQRT((XNUC(JNUC)-XNUC(INUC))**2
     +                   +(YNUC(JNUC)-YNUC(INUC))**2
     +                   +(ZNUC(JNUC)-ZNUC(INUC))**2)
C
      KKK=0
      DO 100 NUC=1,NUCS
      KKKTAB(NUC)=KKK
C
C     INITIALIZE RADIAL MESH:
C
      RMID(NUC)=1.D0/DFLOAT(NUCZ(NUC))**(1.D0/3.D0)
      CALL RMESH(NUC,NR(NUC),RMID(NUC),
     +      RADS(1,NUC),RQ(1,NUC),RQQ(1,NUC),WRADS(1,NUC))
C
C     INITIALIZE ANGULAR MESH:
C
      CALL LEBQUA(NL(NUC),XANGS(:,NUC),YANGS(:,NUC),ZANGS(:,NUC),
     +            WANGS(:,NUC),NL_SUPPLIED,LMAX_SUPPLIED)
      IF (NL(NUC).NE.NL_SUPPLIED) THEN
        WRITE(IOUT,*)'*** Trouble with angular mesh code!!!'
        ERRFLG = .TRUE.
        RETURN
      END IF
C
      DO 101 IR=1,NR(NUC)
      R=RADS(IR,NUC)
      DO 102 IL=1,NL(NUC)
      KKK=KKK+1
      X=XNUC(NUC)+R*XANGS(IL,NUC)
      Y=YNUC(NUC)+R*YANGS(IL,NUC)
      Z=ZNUC(NUC)+R*ZANGS(IL,NUC)
      DO 110 INUC=2,NUCS
      DO 111 JNUC=1,INUC-1
      R1=DSQRT((X-XNUC(INUC))**2+(Y-YNUC(INUC))**2+(Z-ZNUC(INUC))**2)
      R2=DSQRT((X-XNUC(JNUC))**2+(Y-YNUC(JNUC))**2+(Z-ZNUC(JNUC))**2)
      HYPR=(R1-R2)/RR(INUC,JNUC)
      HYPR=1.5D0*HYPR-0.5D0*HYPR**3
      HYPR=1.5D0*HYPR-0.5D0*HYPR**3
      HYPR=1.5D0*HYPR-0.5D0*HYPR**3
      HYPR=1.5D0*HYPR-0.5D0*HYPR**3
      CUTOFF(INUC,JNUC)=(1.D0-HYPR)/2.D0
111   CUTOFF(JNUC,INUC)=(1.D0+HYPR)/2.D0
110   CUTOFF(INUC,INUC) = 1.D0
      CUTOFF(   1,   1) = 1.D0
      VP0=1.D0
      VPSUM=0.D0
      DO 120 INUC=1,NUCS
      VP0=VP0*CUTOFF(NUC,INUC)
      VPI=1.D0
      DO 121 JNUC=1,NUCS
121   VPI=VPI*CUTOFF(INUC,JNUC)
120   VPSUM=VPSUM+VPI
      WNUC(KKK)=VP0/VPSUM
      XMESH(KKK)=X
      YMESH(KKK)=Y
      ZMESH(KKK)=Z
      WINTS(KKK)=WNUC(KKK)*WRADS(IR,NUC)*WANGS(IL,NUC)
102   CONTINUE
101   CONTINUE
100   CONTINUE
C
      MYLM = MAXVAL(LMAX)
C
      IF(KKK.NE.NNN)THEN
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** MESH CHECK FAILURE IN ROUTINE MESH0 ***'
      WRITE(IOUT,*)'*** KKK:',KKK,' NNN:',NNN
      ERRFLG = .TRUE.
      END IF
C
      IF (MINVAL(WINTS(:)).LT.0.D0.and.myRank.EQ.0) THEN
      WRITE(IOUT,'(/''*** WARNING! Negative integration weights!'')')
      END IF
C
      RETURN
      END
C
C
C
      SUBROUTINE MESH1
C
C     MESH POINTS AND INTEGRATION WEIGHTS FOR SINGLE ATOM.
C
      USE ERRORHANDLING
      USE GRID
      USE IOUNITS
      USE NuclearFramework
      USE NuclearMeshes
      USE SIZES
      USE SphericalHarmonics, ONLY : MYLM
      IMPLICIT NONE
      INTEGER KKK,IR,IL,NL_SUPPLIED,LMAX_SUPPLIED
      REAL(KIND(1D0)) R
C
C     INITIALIZE RADIAL MESH:
C
      RMID(1)=1.D0/DFLOAT(NUCZ(1))**(1.D0/3.D0)
      CALL RMESH(1,NR(1),RMID(1),RADS(1,1),RQ(1,1),RQQ(1,1),WRADS(1,1))
C
C     INITIALIZE ANGULAR MESH:
C
      CALL LEBQUA(NL(1),XANGS(:,1),YANGS(:,1),ZANGS(:,1),
     +            WANGS(:,1),NL_SUPPLIED,LMAX_SUPPLIED)
      IF (NL(1).NE.NL_SUPPLIED) THEN
        WRITE(IOUT,*)'*** Trouble with angular mesh code!!!'
        ERRFLG = .TRUE.
        RETURN
      END IF
      MYLM = LMAX(1)
C
      KKK=0
      DO 101 IR=1,NR(1)
        R=RADS(IR,1)
        DO 102 IL=1,NL(1)
          KKK=KKK+1
          WNUC(KKK)=1.D0
          XMESH(KKK)=XNUC(1)+R*XANGS(IL,1)
          YMESH(KKK)=YNUC(1)+R*YANGS(IL,1)
          ZMESH(KKK)=ZNUC(1)+R*ZANGS(IL,1)
          WINTS(KKK)=WRADS(IR,1)*WANGS(IL,1)
102     CONTINUE
101   CONTINUE
      KKKTAB(1)=0
C
      IF(KKK.NE.NNN)THEN
        WRITE(IOUT,'(/''*** MESH CHECK FAILURE IN ROUTINE MESH1 ***'')')
        ERRFLG = .TRUE.
      END IF
C
      RETURN
      END
C
C
C
      SUBROUTINE RMESH(INUC,N,RMID,R,RQ,RQQ,WINTR)
C
C     RADIAL MESH, INTEGRATION WEIGHTS,
C     DERIVATIVES OF VARIABLE R WITH RESPECT TO VARIABLE Q.
C
C     THE Q-MESH IS UNIFORM ON THE INTERVAL (0,+1). TRANSFORMATION IS
C
C                        Q
C         R  =  RMID ---------
C                    ( 1 - Q )
C
C     ALSO,
C     3 AND 7-POINT FINITE DIFFERENCE MATRICES FOR D2(WRT)R ON Q-MESH.
C
      USE CONSTANTS
      USE SecondDerivativeMatrices
      ! TODO refactor so that slices of D2R and TRID2R are passed in as
      ! TODO arguments instead of by use-association.
      IMPLICIT NONE
      ! Arguments
      INTEGER INUC,N
      REAL(KIND(1D0)):: RMID,R(N),RQ(N),RQQ(N),WINTR(N)
      ! Locals
      INTEGER I,II,K
      REAL(KIND(1D0)) A1,A2,C1,C2,CR1,CR2,H,Q
C
C     FINITE DIFFERENCE COEFFICIENTS:
C
      REAL(KIND(1D0)),PARAMETER :: FD1(-3:3,-3:3)=reshape((/
     +    0.D0,   0.D0,  -12.D0,  -12.D0,  -12.D0,   -4.D0,  -2.D0,
     +    0.D0,   6.D0,  108.D0,  108.D0,  108.D0,   30.D0,  12.D0,
     +   -6.D0, -60.D0, -540.D0, -540.D0, -540.D0, -120.D0, -36.D0,
     +  -20.D0, -40.D0,    0.D0,    0.D0,    0.D0,   40.D0,  20.D0,
     +   36.D0, 120.D0,  540.D0,  540.D0,  540.D0,   60.D0,   6.D0,
     +  -12.D0, -30.D0, -108.D0, -108.D0, -108.D0,   -6.D0,   0.D0,
     +    2.D0,   4.D0,   12.D0,   12.D0,   12.D0,    0.D0,   0.D0/),
     +  (/7,7/))
      REAL(KIND(1D0)),PARAMETER :: FD2(-3:3,-3:3)=reshape((/
     +    0.D0,    0.D0,    4.D0,    4.D0,    4.D0,    0.D0,  -1.D0,
     +    0.D0,   -5.D0,  -54.D0,  -54.D0,  -54.D0,   -5.D0,   4.D0,
     +   11.D0,   80.D0,  540.D0,  540.D0,  540.D0,   80.D0,   6.D0,
     +  -20.D0, -150.D0, -980.D0, -980.D0, -980.D0, -150.D0, -20.D0,
     +    6.D0,   80.D0,  540.D0,  540.D0,  540.D0,   80.D0,  11.D0,
     +    4.D0,   -5.D0,  -54.D0,  -54.D0,  -54.D0,   -5.D0,   0.D0,
     +   -1.D0,    0.D0,    4.D0,    4.D0,    4.D0,    0.D0,   0.D0/),
     +  (/7,7/))
C
      H=1.D0/(N+1)
      DO 100 I=1,N
      Q=H*I
      R(I)=RMID*Q/(1.D0-Q)
      RQ(I)=RMID/(1.D0-Q)**2
      RQQ(I)=2.D0*RMID/(1.D0-Q)**3
100   WINTR(I)=FOURPI*H*R(I)**2*RQ(I)
C
C     CONSTRUCT SECOND DERIVATIVE MATRIX (D2R):
C
      DO 200 I=1,N
      IF(I.EQ.1.OR.I.EQ.N)THEN
      C1=1.D0/24.D0/H
      C2=1.D0/12.D0/H**2
      ELSE IF(I.EQ.2.OR.I.EQ.N-1)THEN
      C1=1.D0/120.D0/H
      C2=1.D0/ 60.D0/H**2
      ELSE
      C1=1.D0/720.D0/H
      C2=1.D0/360.D0/H**2
      END IF
      CR2=C2/RQ(I)**2
      CR1=-C1*RQQ(I)/RQ(I)**3
      II=0
      IF(I.LE.3)II=I-4
      IF(I.GE.N-2)II=I-N+3
      DO 210 K=1,7
210   D2R(I,K,INUC)=CR1*FD1(II,K-4)+CR2*FD2(II,K-4)
200   CONTINUE
C
      DO 300 I=1,N
      A2=1.D0/RQ(I)**2/H**2
      A1=-0.5D0*RQQ(I)/RQ(I)**3/H
      TRID2R(I,1,INUC)=A2-A1
      TRID2R(I,2,INUC)=-2.D0*A2
      TRID2R(I,3,INUC)=A2+A1
300   CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE SEARCH(ZSYMB,IATM,ZCORE)
C
      USE AtomicData
      USE CORES
      USE ERRORHANDLING
      USE IOUNITS
      USE KEYWORDS
      IMPLICIT NONE
      ! Arguments
      REAL(KIND(1D0)) ZCORE
      INTEGER IATM
      CHARACTER*2 ZSYMB
      ! Locals
      INTEGER L,IS,NS
      CHARACTER*2 COMPAR
C
      REWIND(IATMS)
100   CONTINUE
      READ(IATMS,'(A2)',END=300)COMPAR
      IF(COMPAR.NE.ZSYMB)GO TO 100
      ZCORE=0.D0
      READ(IATMS,*)LMAX0(IATM)
      DO 200 L=0,LMAX0(IATM)
      READ(IATMS,*)NORBS0(IATM,L)
      NS=NORBS0(IATM,L)
      READ(IATMS,*)(ICORE(IATM,L,IS),IS=1,NS)
      READ(IATMS,*) (OCC1(IATM,L,IS),IS=1,NS)
      READ(IATMS,*) (OCC2(IATM,L,IS),IS=1,NS)
      DO 210 IS=1,NS
      ICORE(IATM,L,IS)=IFRZC*ICORE(IATM,L,IS)
      OCC12(IATM,L,IS)=OCC1(IATM,L,IS)+OCC2(IATM,L,IS)
210   ZCORE = ZCORE + ICORE(IATM,L,IS)*OCC12(IATM,L,IS)
200   CONTINUE
      RETURN
300   CONTINUE
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** UNRECOGNIZED ATOMIC SYMBOL: ',ZSYMB,' ***'
      WRITE(IOUT,*)'*** NUMATM ATOMIC DATA NOT PRESENT ***'
      ERRFLG = .TRUE.
      RETURN
      END
C
C
C
      SUBROUTINE ATOMIN(TCORE,RHO,RHOS,VEL,VNUC,FORCE0)
C
      USE AtomicData
      USE CONSTANTS
      USE CORES
      USE ERRORHANDLING
      USE GRID
      USE IOUNITS
      USE KEYWORDS
      USE NuclearFramework
      USE NuclearMeshes
      USE PROMOLECULE
      USE SIZES
      USE SpinPolarization
      USE SPLINES
      IMPLICIT NONE
C     Arguments
      REAL(KIND(1D0)),INTENT(OUT):: TCORE
      REAL(KIND(1D0)),INTENT(OUT):: RHO(NNN),RHOS(NNN,2)
      REAL(KIND(1D0)),INTENT(OUT):: VEL(NNN),VNUC(NNN),FORCE0(NUCS,3)
C     Locals
      REAL(KIND(1D0)) WORK0(NNN),WORK2(NNN),WORKDC(NNN)
      REAL(KIND(1D0)) RDATA(MDATA)
      REAL(KIND(1D0)) F1(0:MDATA),A1(0:MDATA),B1(0:MDATA),C1(0:MDATA)
      REAL(KIND(1D0)) F2(0:MDATA),A2(0:MDATA),B2(0:MDATA),C2(0:MDATA)
      REAL(KIND(1D0)) F3(0:MDATA),A3(0:MDATA),B3(0:MDATA),C3(0:MDATA)
      REAL(KIND(1D0)) F4(0:MDATA),A4(0:MDATA),B4(0:MDATA),C4(0:MDATA)
      REAL(KIND(1D0)) F5(0:MDATA),A5(0:MDATA),B5(0:MDATA),C5(0:MDATA)
      REAL(KIND(1D0)) CHRG,RMIDZ,ZRMD,ZCORE,ZNORM,SPOL1,SPOL2
      REAL(KIND(1D0)) H,Q,X,Y,Z,R,DQ,ARHO1,ARHO2,AVEL,ARHOC,ATAUC
      REAL(KIND(1D0)) COREH,VALEH,RHOP1,RHOP2,AVNUC,ARHO,SEPMIN
      REAL(KIND(1D0)) DRQ,DVEL,FABS,DARHOC,OCCKIN,CON,FANG,PSI0,PSI2
      LOGICAL CORE
      INTEGER I,IATM,NDATA,NDATA1,KKK,INTQ,INUC,IR,L,LDEG,M,IORB
      INTEGER KKKBOT,KKKTOP,IBASIS
C     Functions
      INTEGER MaxCoreOrbs
      REAL(KIND(1D0)) SPDF
C
      F1(0)=0.D0; F2(0)=0.D0; F3(0)=0.D0; F4(0)=0.D0; F5(0)=0.D0
      REWIND(ICOR0)
      REWIND(ICOR2)
      REWIND(ICRDC)
      if (.not.allocated(RARHO)) then
         allocate(NSPL(NUCS))
         allocate(RARHO(MDATA,NUCS))
         allocate(RAVEL(MDATA,NUCS))
         allocate(DCOR(MR,NUCS))
         allocate(ICORE(NUCS,0:3,8))
         allocate(LMAX0(NUCS))
         allocate(NORBS0(NUCS,0:3))
         allocate(OCC1(NUCS,0:3,8))
         allocate(OCC2(NUCS,0:3,8))
         allocate(OCC12(NUCS,0:3,8))
         MCORBS = MaxCoreOrbs(NUCZ,NUCS)
         allocate(ICNTR(MCORBS))
         allocate(LCOR(MCORBS))
      end if
C
      NCORBS=0
      NBASIS=0
      TCORE=0.D0
      DO 10 I=1,NNN
      RHO(I)=0.D0
      VEL(I)=0.D0
      VNUC(I)=0.D0
      RHOS(I,1)=0.D0
      RHOS(I,2)=0.D0
      RHOCOR(I)=0.D0
10    TAUCOR(I)=0.D0
      DO 11 I=1,NUCS
      FORCE0(I,1)=0.D0
      FORCE0(I,2)=0.D0
11    FORCE0(I,3)=0.D0
C
      DO 100 IATM=1,NUCS
C
      CHRG=DFLOAT(NUCZ(IATM))
      RMIDZ=1.D0/CHRG**THRD
      ZRMD=CHRG*RMIDZ
      CALL SEARCH(ZSYMB(IATM),IATM,ZCORE)
      if (ERRFLG) return
      ZNORM=0.5D0*(CHRG-ZCORE)
      SPOL1=1.D0+SPOL(IATM,1)/ZNORM
      SPOL2=1.D0+SPOL(IATM,2)/ZNORM
C
      READ(IATMS,*)NDATA
      IF(NDATA.GT.MDATA)THEN
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** NDATA GREATER THAN MDATA ***'
      ERRFLG = .TRUE.
      RETURN
      END IF
      NDATA1=NDATA+1
      H=1.D0/NDATA1
      DO 101 I=1,NDATA
      Q=H*I
101   RDATA(I)=RMIDZ*Q/(1.D0-Q)
      NSPL(IATM)=NDATA
C
      READ(IATMS,*)(F1(I),I=1,NDATA)
      READ(IATMS,*)(F2(I),I=1,NDATA)
      READ(IATMS,*)(F3(I),I=1,NDATA)
      READ(IATMS,*)(F4(I),I=1,NDATA)
      READ(IATMS,*)(F5(I),I=1,NDATA)
      DO 200 I=1,NDATA
      F1(I)=RDATA(I)**2*F1(I)
      F2(I)=RDATA(I)**2*F2(I)
      F3(I)=RDATA(I)*F3(I)
      F4(I)=RDATA(I)**2*F4(I)
      F5(I)=RDATA(I)**2*F5(I)
      RAVEL(I,IATM)=F3(I)
200   RARHO(I,IATM)=(F1(I)+F2(I))/RDATA(I)
      IF(IFRZC.EQ.0)THEN
      DO 201 I=1,NDATA
      F4(I)=0.D0
201   F5(I)=0.D0
      END IF
      IF(IOPTG.EQ.-1)GO TO 100
      CALL SPLINE(H,F1,A1,B1,C1,NDATA,1,0.D0,0.D0)
      CALL SPLINE(H,F2,A2,B2,C2,NDATA,1,0.D0,0.D0)
      CALL SPLINE(H,F3,A3,B3,C3,NDATA,2,0.D0,CHRG)
      CALL SPLINE(H,F4,A4,B4,C4,NDATA,1,0.D0,0.D0)
      CALL SPLINE(H,F5,A5,B5,C5,NDATA,2,0.D0,0.D0)
C
      DO 210 KKK=1,NNN
      X=XMESH(KKK)-XNUC(IATM)
      Y=YMESH(KKK)-YNUC(IATM)
      Z=ZMESH(KKK)-ZNUC(IATM)
      R=DSQRT(X**2+Y**2+Z**2)
      Q=R/(R+RMIDZ)
      INTQ=IDINT(NDATA1*Q)
      DQ=Q-INTQ*H
      ARHO1=(F1(INTQ)+DQ*(A1(INTQ)+DQ*(B1(INTQ)+DQ*C1(INTQ))))/R**2
      ARHO2=(F2(INTQ)+DQ*(A2(INTQ)+DQ*(B2(INTQ)+DQ*C2(INTQ))))/R**2
      AVEL =(F3(INTQ)+DQ*(A3(INTQ)+DQ*(B3(INTQ)+DQ*C3(INTQ))))/R
      ARHOC=(F4(INTQ)+DQ*(A4(INTQ)+DQ*(B4(INTQ)+DQ*C4(INTQ))))/R**2
      ATAUC=(F5(INTQ)+DQ*(A5(INTQ)+DQ*(B5(INTQ)+DQ*C5(INTQ))))/R**2
C
      COREH=0.5D0*ARHOC
      VALEH=0.5D0*(ARHO1+ARHO2-ARHOC)
      RHOP1=SPOL1*VALEH+COREH
      RHOP2=SPOL2*VALEH+COREH
C
      AVNUC=-CHRG/R
      ARHO1=DABS(ARHO1)
      ARHO2=DABS(ARHO2)
      ARHOC=DABS(ARHOC)
      RHOP1=DABS(RHOP1)
      RHOP2=DABS(RHOP2)
      AVEL =DABS(AVEL )
      ARHO =ARHO1+ARHO2
C
      WORK0(KKK)=ARHO
      RHO(KKK)=RHO(KKK)+ARHO
      VEL(KKK)=VEL(KKK)+AVEL
      VNUC(KKK)=VNUC(KKK)+AVNUC
      RHOS(KKK,1)=RHOS(KKK,1)+RHOP1
      RHOS(KKK,2)=RHOS(KKK,2)+RHOP2
      RHOCOR(KKK)=RHOCOR(KKK)+ARHOC
      TAUCOR(KKK)=TAUCOR(KKK)+ATAUC
      TCORE=TCORE-0.5D0*WINTS(KKK)*ATAUC
210   CONTINUE
C
      WRITE(IBAS0,REC=IATM,ERR=9001)(WORK0(I),I=1,NNN)
C
      SEPMIN=1.D11
      DO 300 INUC=1,NUCS
      IF(INUC.EQ.IATM)GO TO 300
      X=XNUC(INUC)-XNUC(IATM)
      Y=YNUC(INUC)-YNUC(IATM)
      Z=ZNUC(INUC)-ZNUC(IATM)
      R=DSQRT(X**2+Y**2+Z**2)
      IF(R.LT.SEPMIN)SEPMIN=R
      Q=R/(R+RMIDZ)
      DRQ=RMIDZ/(1.D0-Q)**2
      INTQ=IDINT(NDATA1*Q)
      DQ=Q-INTQ*H
      AVEL=(F3(INTQ)+DQ*(A3(INTQ)+DQ*(B3(INTQ)+DQ*C3(INTQ))))/R
      DVEL=(A3(INTQ)+2.D0*B3(INTQ)*DQ+3.D0*C3(INTQ)*DQ*DQ)/DRQ
      DVEL=(DVEL-AVEL)/R
      FABS=NUCZ(INUC)*(DVEL+CHRG/R**2)
      FORCE0(INUC,1)=FORCE0(INUC,1)+FABS*X/R
      FORCE0(INUC,2)=FORCE0(INUC,2)+FABS*Y/R
      FORCE0(INUC,3)=FORCE0(INUC,3)+FABS*Z/R
300   CONTINUE
      DO 310 IR=1,NR(IATM)
      R=RADS(IR,IATM)
      Q=R/(R+RMIDZ)
      DRQ=RMIDZ/(1.D0-Q)**2
      INTQ=IDINT(NDATA1*Q)
      DQ=Q-INTQ*H
      ARHOC=(F4(INTQ)+DQ*(A4(INTQ)+DQ*(B4(INTQ)+DQ*C4(INTQ))))/R**2
      DARHOC=(A4(INTQ)+2.D0*B4(INTQ)*DQ+3.D0*C4(INTQ)*DQ*DQ)/DRQ
      DCOR(IR,IATM)=(DARHOC-2.D0*R*ARHOC)/R**2
      IF(R.GT.(SEPMIN-0.1D0))DCOR(IR,IATM)=0.D0
310   CONTINUE
C
      DO 400 L=0,LMAX0(IATM)
      LDEG=2*L+1
C
      DO 400 IORB=1,NORBS0(IATM,L)
      READ(IATMS,*)(F1(I),I=1,NDATA)
      READ(IATMS,*)(F2(I),I=1,NDATA)
      CORE=.FALSE.
      IF(ICORE(IATM,L,IORB).EQ.1)CORE=.TRUE.
      IF(CORE)THEN
      OCCKIN=0.D0
      ELSE
      OCCKIN=OCC12(IATM,L,IORB)/LDEG
      END IF
      CON=DSQRT(FOURPI)
      DO 401 I=1,NDATA
      F1(I)=CON*RDATA(I)*F1(I)
401   F2(I)=CON*RDATA(I)**2*F2(I)
      IF(L.EQ.0)THEN
      CALL SPLINE(H,F1,A1,B1,C1,NDATA,0,ZRMD,0.D0)
      CALL SPLINE(H,F2,A2,B2,C2,NDATA,2,0.D0,0.D0)
      ELSE
      CALL SPLINE(H,F1,A1,B1,C1,NDATA,1,0.D0,0.D0)
      CALL SPLINE(H,F2,A2,B2,C2,NDATA,1,0.D0,0.D0)
      END IF
C
      DO 410 M=1,LDEG
      IF(CORE)THEN
        NCORBS=NCORBS+1
        IF (NCORBS>MCORBS) THEN
          WRITE(IOUT,*) '*** Too many core orbitals?!? ***'
          GOTO 9999
        END IF
        LCOR(NCORBS)=L
        ICNTR(NCORBS)=IATM
      ELSE
        NBASIS=NBASIS+1
      END IF
C
      DO 420 KKK=1,NNN
      X=XMESH(KKK)-XNUC(IATM)
      Y=YMESH(KKK)-YNUC(IATM)
      Z=ZMESH(KKK)-ZNUC(IATM)
      R=DSQRT(X**2+Y**2+Z**2)
      X=X/R
      Y=Y/R
      Z=Z/R
      FANG=SPDF(L,M,X,Y,Z)
      Q=R/(R+RMIDZ)
      INTQ=IDINT(NDATA1*Q)
      DQ=Q-INTQ*H
      PSI0=(F1(INTQ)+DQ*(A1(INTQ)+DQ*(B1(INTQ)+DQ*C1(INTQ))))/R
      PSI2=(F2(INTQ)+DQ*(A2(INTQ)+DQ*(B2(INTQ)+DQ*C2(INTQ))))/R**2
      WORK0(KKK)=PSI0*FANG
      WORK2(KKK)=PSI2*FANG
      WORKDC(KKK)=0.D0
420   CONTINUE
C
      KKKBOT=KKKTAB(IATM)+1
      IF(IATM.LT.NUCS)THEN
      KKKTOP=KKKTAB(IATM+1)
      ELSE
      KKKTOP=NNN
      END IF
      DO 421 KKK=KKKBOT,KKKTOP
421   WORKDC(KKK)=WORK0(KKK)
C
      IF(CORE)THEN
      WRITE(ICOR0,ERR=9002)(WORK0(I),I=1,NNN)
      WRITE(ICOR2,ERR=9003)(WORK2(I),I=1,NNN)
      WRITE(ICRDC,ERR=9004)(WORKDC(I),I=1,NNN)
      ELSE
      WRITE(MOSA0,REC=NBASIS,ERR=9005)(WORK0(I),I=1,NNN)
      WRITE(MOSA2,REC=NBASIS,ERR=9006)(WORK2(I),I=1,NNN)
      WRITE(MODCA,REC=NBASIS,ERR=9007)(WORKDC(I),I=1,NNN)
      END IF
410   CONTINUE
400   CONTINUE
C
100   CONTINUE
      IF(IOPTG.EQ.-1)RETURN
C
      DO 500 KKK=1,NNN
      RHO0(KKK)=RHO(KKK)
500   VEL0(KKK)=VEL(KKK)
C
      REWIND(IHRSH)
      DO 510 IATM=1,NUCS
      READ(IBAS0,REC=IATM)(WORK0(I),I=1,NNN)
      DO 511 KKK=1,NNN
511   WORK0(KKK)=WORK0(KKK)/RHO(KKK)
      WRITE(IHRSH,ERR=9008)(WORK0(I),I=1,NNN)
510   CONTINUE
C
C     TODO: Can we simplify out MOS* since we're now using direct access?
      DO 520 IBASIS=1,NBASIS
      READ(MOSA0,REC=IBASIS)(WORK0(I),I=1,NNN)
      READ(MOSA2,REC=IBASIS)(WORK2(I),I=1,NNN)
      READ(MODCA,REC=IBASIS)(WORKDC(I),I=1,NNN)
      IF(NCORBS.NE.0)CALL ORTHC(WORK0,WORK2,WORKDC)
      WRITE(IBAS0,REC=IBASIS,ERR=9009)(WORK0(I),I=1,NNN)
      WRITE(IBAS2,REC=IBASIS,ERR=9010)(WORK2(I),I=1,NNN)
      WRITE(IBSDC,REC=IBASIS,ERR=9011)(WORKDC(I),I=1,NNN)
520   CONTINUE
C
      IF (IFRZC.EQ.0) TCORE=0.D0
C
      RETURN
 9001 WRITE(IOUT,*) '*** Error writing IBAS0 in ATOMIN (1)'
      GOTO 9999
 9002 WRITE(IOUT,*) '*** Error writing ICOR0 in ATOMIN'
      GOTO 9999
 9003 WRITE(IOUT,*) '*** Error writing ICOR2 in ATOMIN'
      GOTO 9999
 9004 WRITE(IOUT,*) '*** Error writing ICRDC in ATOMIN'
      GOTO 9999
 9005 WRITE(IOUT,*) '*** Error writing MOSA0 in ATOMIN'
      GOTO 9999
 9006 WRITE(IOUT,*) '*** Error writing MOSA2 in ATOMIN'
      GOTO 9999
 9007 WRITE(IOUT,*) '*** Error writing MODCA in ATOMIN'
      GOTO 9999
 9008 WRITE(IOUT,*) '*** Error writing IHRSH in ATOMIN'
      GOTO 9999
 9009 WRITE(IOUT,*) '*** Error writing IBAS0 in ATOMIN (2)'
      GOTO 9999
 9010 WRITE(IOUT,*) '*** Error writing IBAS2 in ATOMIN'
      GOTO 9999
 9011 WRITE(IOUT,*) '*** Error writing IBSDC in ATOMIN'
      GOTO 9999
 9999 ERRFLG = .TRUE.
      RETURN
      END
C
C
C
      FUNCTION MaxCoreOrbs(NUCZ,NUCS)
      ! Compute the typical maximum number of core orbitals.
      ! The actual number of core orbitals is determined later in
      ! ATOMIN while reading the atomic data file IATMS.
      IMPLICIT NONE
      INTEGER MaxCoreOrbs
      INTEGER NUCS, NUCZ(NUCS)
      INTEGER INUC
      MaxCoreOrbs = 0
      DO INUC=1,NUCS
         SELECT CASE (NUCZ(INUC))
         CASE (1:2)
            ! no core orbitals for Z <= 2
         CASE (3:10)
            MaxCoreOrbs = MaxCoreOrbs + 1   ! 1s
         CASE (11:18)
            MaxCoreOrbs = MaxCoreOrbs + 5   ! +2s+2p
         CASE (19:30)
            MaxCoreOrbs = MaxCoreOrbs + 9   ! +3s+3p
         CASE (31:36)
            MaxCoreOrbs = MaxCoreOrbs + 14  ! +3d
         CASE (37:48)
            MaxCoreOrbs = MaxCoreOrbs + 18  ! +4s+4p
         CASE (49:54)
            MaxCoreOrbs = MaxCoreOrbs + 23  ! +4d
         CASE (55:80)
            MaxCoreOrbs = MaxCoreOrbs + 27  ! +5s+5p
         CASE (81:86)
            MaxCoreOrbs = MaxCoreOrbs + 39  ! +5d+4f
         CASE (87:112)
            MaxCoreOrbs = MaxCoreOrbs + 43  ! +6s+6p
         CASE DEFAULT
            ! ought to throw an exception here
         END SELECT
      END DO
      RETURN
      END FUNCTION
C
C
C
      FUNCTION SPDF(L,M,X,Y,Z)
C
C     CALCULATION OF ATOMIC ANGULAR FUNCTIONS.
C     FROM COTTON: CHEMICAL APPLICATIONS OF GROUP THEORY.
C
      USE CONSTANTS
      IMPLICIT NONE
      INTEGER L,M
      REAL(KIND(1D0)) SPDF,X,Y,Z
C
      GO TO (10,20,30,40),L+1
C
10    SPDF=0.5D0/DSQRT(PI)
      RETURN
C
20    GO TO (201,202,203),M
201   SPDF=0.5D0*DSQRT(3.D0/PI)*Z
      RETURN
202   SPDF=0.5D0*DSQRT(3.D0/PI)*X
      RETURN
203   SPDF=0.5D0*DSQRT(3.D0/PI)*Y
      RETURN
C
30    GO TO (301,302,303,304,305),M
301   SPDF=0.25D0*DSQRT(5.D0/PI)*(3.D0*Z**2-1.D0)
      RETURN
302   SPDF=0.5D0*DSQRT(15.D0/PI)*X*Z
      RETURN
303   SPDF=0.5D0*DSQRT(15.D0/PI)*Y*Z
      RETURN
304   SPDF=0.25D0*DSQRT(15.D0/PI)*(X**2-Y**2)
      RETURN
305   SPDF=0.5D0*DSQRT(15.D0/PI)*X*Y
      RETURN
C
40    GO TO (401,402,403,404,405,406,407),M
401   SPDF=0.25D0*DSQRT(7.D0/PI)*Z*(5.D0*Z**2-3.D0)
      RETURN
402   SPDF=0.125D0*DSQRT(42.D0/PI)*X*(5.D0*Z**2-1.D0)
      RETURN
403   SPDF=0.125D0*DSQRT(42.D0/PI)*Y*(5.D0*Z**2-1.D0)
      RETURN
404   SPDF=0.25D0*DSQRT(105.D0/PI)*Z*(X**2-Y**2)
      RETURN
405   SPDF=0.5D0*DSQRT(105.D0/PI)*Z*X*Y
      RETURN
406   SPDF=0.125D0*DSQRT(70.D0/PI)*X*(X**2-3.D0*Y**2)
      RETURN
407   SPDF=0.125D0*DSQRT(70.D0/PI)*Y*(3.D0*X**2-Y**2)
      RETURN
C
      END
C
C
C
      SUBROUTINE ORTHC(F0,F2,FDC)
C
C     ORTHOGONALIZATION WITH RESPECT TO CORE ORBITALS.
C
      USE GRID
      USE IOUNITS
      USE SIZES
      IMPLICIT NONE
      REAL(KIND(1D0)) F0(NNN),F2(NNN),FDC(NNN)
      REAL(KIND(1D0)) SUM,DUM(NNN,3),STOR(NNN)
      REAL(KIND(1D0)) CORE0(NNN),CORE2(NNN),COREDC(NNN)
      INTEGER KKK,ICORB,I
C
      DO 10 KKK=1,NNN
10    STOR(KKK)=F0(KKK)
      REWIND(ICOR0)
      REWIND(ICOR2)
      REWIND(ICRDC)
      DO 100 ICORB=1,NCORBS
      READ(ICOR0)(CORE0(I),I=1,NNN)
      READ(ICOR2)(CORE2(I),I=1,NNN)
      READ(ICRDC)(COREDC(I),I=1,NNN)
      SUM=0.D0
      DO 110 KKK=1,NNN
110   SUM=SUM+WINTS(KKK)*CORE0(KKK)*STOR(KKK)
      DO 120 KKK=1,NNN
      F0(KKK)=F0(KKK)-SUM*CORE0(KKK)
      F2(KKK)=F2(KKK)-SUM*CORE2(KKK)
120   FDC(KKK)=FDC(KKK)-SUM*COREDC(KKK)
100   CONTINUE
C
      RETURN
      END
C
C
C
      subroutine PopulateOrbitals(NELS,NV0S)
      use ErrorHandling
      use IOUnits, only : IOUT
      use KeyWords
      use MPI, only : myrank
      use NondefaultOccupancies
      use OccupationNumbers
      use Sizes
      implicit none
      integer NELS, NV0S
      integer I, IMORB, NCORE
      real(kind(1d0)) CHKSUM, DIFFQ
C
      if (.not.allocated(OCCS)) THEN
         allocate(OCCS(NBASIS,2))
         OCCS(:,:) = 0.D0
      end if
C
      NCORE=2*NCORBS
      NMORBS=(NELS-NCORE)/2
      DO 120 I=1,NMORBS
         OCCS(I,1)=1.D0
         OCCS(I,2)=1.D0
  120 CONTINUE
      IF(MOD(NELS,2).GT.0)THEN
         NMORBS=NMORBS+1
         OCCS(I,1)=1.D0
         OCCS(I,2)=0.D0
      END IF
C
C     TODO: Change input standard so that nondefault occupancies always
C     TODO  include the core orbitals, and we can subtract them out here
      IF(IOCCS.EQ.1)THEN
         if(myrank.eq.0) WRITE(IOUT,128)
  128    format(//'NONDEFAULT OCCUPANCIES:')
         DO 130 I=1,NCHNGS
            IMORB=ICHANGE(I)
            OCCS(IMORB,1)=CHNGOCC(I,1)
            OCCS(IMORB,2)=CHNGOCC(I,2)
            IF(IMORB.GT.NMORBS)NMORBS=IMORB
            if(myrank.eq.0) 
     &               WRITE(IOUT,129)IMORB,OCCS(IMORB,1),OCCS(IMORB,2)
  129       format(1X,I4,4X,2F10.5)
  130    CONTINUE
      END IF
      CHKSUM=DFLOAT(NCORE)
      DO 140 I=1,NMORBS
  140 CHKSUM=CHKSUM+OCCS(I,1)+OCCS(I,2)
      DIFFQ=CHKSUM-DFLOAT(NELS)
      IF(DABS(DIFFQ).GT.1.D-03)THEN
         WRITE(IOUT,'(/)')
         WRITE(IOUT,*)'*** OCCUPANCY CHECKSUM NOT EQUAL TO NELS ***'
         ERRFLG = .TRUE.
         RETURN
      END IF
C
      NSPINS=1
      DO 150 I=1,NMORBS
      IF (OCCS(I,1).NE.OCCS(I,2)) NSPINS=2
150   CONTINUE
      IF(NV0S.EQ.2)NSPINS=2
C
      RETURN
      END
C
C
C
      SUBROUTINE LCAO0
C
C     INITIALIZATION OF LCAO S AND T MATRICES.
C
      USE GRID
      USE IOUNITS
      USE LCAOMatrices
      USE SIZES
      IMPLICIT NONE
      REAL(KIND(1D0)) BAS0A(NNN),BAS2A(NNN),BAS0B(NNN),BAS2B(NNN)
      INTEGER I,II,JJ,KKK
      REAL(KIND(1D0)) SSUM,TSUM
C
      if (.not.allocated(EIGS)) then
         allocate(EIGS(NBASIS,2)) ! TODO could we make this 1 if NSPINS==1?
         allocate(RATSUM(NBASIS,2))  ! TODO make this NMORBS instead?
         allocate(SMAT(NBASIS,NBASIS))
         allocate(TMAT(NBASIS,NBASIS))
      end if
C
      DO 100 II=1,NBASIS
      READ(IBAS0,REC=II)(BAS0A(I),I=1,NNN)
      READ(IBAS2,REC=II)(BAS2A(I),I=1,NNN)
      DO 100 JJ=1,II
      READ(IBAS0,REC=JJ)(BAS0B(I),I=1,NNN)
      READ(IBAS2,REC=JJ)(BAS2B(I),I=1,NNN)
      SSUM=0.D0
      TSUM=0.D0
      DO 110 KKK=1,NNN
      SSUM=SSUM+WINTS(KKK)*BAS0A(KKK)*BAS0B(KKK)
      TSUM=TSUM-0.25D0*WINTS(KKK)*(BAS2A(KKK)*BAS0B(KKK)
     +                            +BAS0A(KKK)*BAS2B(KKK))
110   CONTINUE
      SMAT(JJ,II)=SSUM
      TMAT(JJ,II)=TSUM
100   TMAT(II,JJ)=TSUM
C
      RETURN
      END
C
C
C
C$PRAGMA SUN OPT=5
      SUBROUTINE POISS(RHO,VEL,IREF)
C
C     SOLUTION OF POISSON'S EQUATION FOR THE COULOMB POTENTIAL
C     OF AN ARBITRARY SOURCE DENSITY DEFINED ON THE MOLECULAR MESH.
C
C     IF IREF.EQ.1, SUBTRACT REFERENCE DENSITY.
C     We determined in Sept 2003 that using the difference density lead to
C     significant orientation dependence in the total energy.  Therefore
C     we now use only the total density and this option is not currently
C     exercised anywhere in the program. --- rmd
C
      USE CONSTANTS
      USE ERRORHANDLING
      USE GRID
      USE IOUNITS
      USE NuclearFramework
      USE NuclearMeshes
      USE PROMOLECULE
      USE SecondDerivativeMatrices
      USE SIZES
      USE SphericalHarmonics
      IMPLICIT NONE
C     Arguments
      REAL(KIND(1D0)),INTENT(OUT):: RHO(NNN)
      REAL(KIND(1D0)),INTENT(OUT):: VEL(NNN)
      INTEGER,        INTENT(IN) :: IREF
C     Locals
      REAL(KIND(1D0)) ULM(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) A(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) B(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) C(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) BND(10,-2:MR+3)
      REAL(KIND(1D0)) RHOLM(MR),XX(MR)
      REAL(KIND(1D0)) RS(ML),XUNIT(ML),YUNIT(ML),ZUNIT(ML),DQS(ML)
      REAL(KIND(1D0)) H,CHARGE,SUM,X,Y,Z,R,Q,DQ,RFAC
      INTEGER IPVT(MR),INTQS(ML),KKK,I,ICEL,NRADS,NANGS,NWAVES,IY,JY
      INTEGER IR,IL,L,L2,K,IERR,INUC,INTQ
C
      IF(IREF.EQ.1)THEN
      DO 10 KKK=1,NNN
10    RHO(KKK)=RHO(KKK)-RHO0(KKK)
      END IF
C
      DO 11 I=1,NNN
11    VEL(I) = 0.D0
C
      DO 100 ICEL=1,NUCS
C
      NRADS=NR(ICEL)
      NANGS=NL(ICEL)
      NWAVES=LMAX(ICEL)
      H=1.D0/(NRADS+1)
C
      KKK=KKKTAB(ICEL)
      CHARGE=0.D0
      DO 101 IR=1,NRADS
      DO 101 IL=1,NANGS
      KKK=KKK+1
101   CHARGE=CHARGE+WINTS(KKK)*RHO(KKK)
C
C     SPHERICAL EXPANSION OF THE DENSITY:
C
      CALL YCALC(XANGS(:,ICEL),YANGS(:,ICEL),ZANGS(:,ICEL),NANGS,NWAVES)
C
      DO 200 IY=0,NWAVES
      DO 200 JY=0,NWAVES
C
      KKK=KKKTAB(ICEL)
      DO 210 IR=1,NRADS
      SUM=0.D0
      DO 211 IL=1,NANGS
211   SUM=SUM+WANGS(IL,ICEL)*YLM(IL,IY,JY)*WNUC(KKK+IL)*RHO(KKK+IL)
      RHOLM(IR)=SUM
210   KKK=KKK+NANGS
C
C     SOLVE RADIAL LM EQUATION:
C
      L=MAX0(IY,JY)
      L2=L*(L+1)
C
C     L.H.S. (LINPACK BAND MATRIX):
C
      DO 220 K=1,7
      DO 220 IR=1,NRADS
220   BND(11-K,IR+K-4)=D2R(IR,K,ICEL)
      DO 221 IR=1,NRADS
      XX(IR)=-FOURPI*RADS(IR,ICEL)*RHOLM(IR)
221   BND(7,IR)=BND(7,IR)-L2/RADS(IR,ICEL)**2
C
C     INFINITE-R BOUNDARY CONDITION:
C
      IF(L.EQ.0)THEN
      XX(NRADS)=XX(NRADS)-CHARGE*D2R(NRADS,5,ICEL)
      XX(NRADS-1)=XX(NRADS-1)-CHARGE*D2R(NRADS-1,6,ICEL)
      XX(NRADS-2)=XX(NRADS-2)-CHARGE*D2R(NRADS-2,7,ICEL)
      END IF
C
      CALL DGBSV(NRADS,3,3,1,BND(1:10,1:),10,IPVT,XX,MR,IERR)
      IF(IERR.NE.0)THEN
      WRITE(IOUT,'(/''*** POISS: ERROR IN ROUTINE DGBSV ***'')')
      ERRFLG = .TRUE.
      RETURN
      END IF
C
      ULM(0,IY,JY)=0.D0
      DO 230 IR=1,NRADS
230   ULM(IR,IY,JY)=XX(IR)
C
      KKK=KKKTAB(ICEL)
      DO 240 IR=1,NRADS
      R=RADS(IR,ICEL)
      DO 241 IL=1,NANGS
241   VEL(KKK+IL)=VEL(KKK+IL)+ULM(IR,IY,JY)*YLM(IL,IY,JY)/R
240   KKK=KKK+NANGS
C
      IF(L.EQ.0)THEN
      CALL SPLINE(H,ULM(0,IY,JY),A(0,IY,JY),B(0,IY,JY),C(0,IY,JY),
     +                                        NRADS,2,0.D0,CHARGE)
      ELSE
      CALL SPLINE(H,ULM(0,IY,JY),A(0,IY,JY),B(0,IY,JY),C(0,IY,JY),
     +                                          NRADS,1,0.D0,0.D0)
      END IF
200   CONTINUE
C
      DO 300 INUC=1,NUCS
      IF(INUC.EQ.ICEL)GO TO 300
C
      KKK=KKKTAB(INUC)
      DO 310 IR=1,NR(INUC)
      DO 320 IL=1,NL(INUC)
C     CALC COORDS ON UNIT SPHERE ABOUT ICEL:
      X=XNUC(INUC)+XANGS(IL,INUC)*RADS(IR,INUC)-XNUC(ICEL)
      Y=YNUC(INUC)+YANGS(IL,INUC)*RADS(IR,INUC)-YNUC(ICEL)
      Z=ZNUC(INUC)+ZANGS(IL,INUC)*RADS(IR,INUC)-ZNUC(ICEL)
      RS(IL)=DSQRT(X**2+Y**2+Z**2)
      XUNIT(IL)=X/RS(IL)
      YUNIT(IL)=Y/RS(IL)
      ZUNIT(IL)=Z/RS(IL)
      Q=RS(IL)/(RS(IL)+RMID(ICEL))
      INTQS(IL)=IDINT((NRADS+1)*Q)
      DQS(IL)=Q-INTQS(IL)*H
320   CONTINUE
      CALL YCALC(XUNIT,YUNIT,ZUNIT,NL(INUC),NWAVES)
      DO 330 IY=0,NWAVES
      DO 330 JY=0,NWAVES
      DO 330 IL=1,NL(INUC)
      INTQ=INTQS(IL)
      DQ=DQS(IL)
      Y=ULM(INTQ,IY,JY)
      RFAC=Y+DQ*(A(INTQ,IY,JY)+DQ*(B(INTQ,IY,JY)+DQ*C(INTQ,IY,JY)))
330   VEL(KKK+IL)=VEL(KKK+IL)+RFAC*YLM(IL,IY,JY)/RS(IL)
310   KKK=KKK+NL(INUC)
C
300   CONTINUE
C
100   CONTINUE
C
C     IF IREF.EQ.1, ADD REFERENCE SYSTEM:
      IF(IREF.EQ.1)THEN
      DO 20 KKK=1,NNN
      RHO(KKK)=RHO(KKK)+RHO0(KKK)
20    VEL(KKK)=VEL(KKK)+VEL0(KKK)
      END IF
C
      RETURN
      END
C
C
C
C$PRAGMA SUN OPT=5
      SUBROUTINE YCALC(X,Y,Z,NANGS,NY)
C
C     GENERATES REAL SPHERICAL HARMONICS UP TO ORDER NY
C     FOR A SET OF POINTS X(I),Y(I),Z(I) ON THE UNIT SPHERE.
C
C     NORMALIZED SUCH THAT  Y  = 1,   I.E. UNIT NORMALIZATION
C                            00
C     WITH RESPECT TO ANGULAR AVERAGE ON SURFACE OF UNIT SPHERE.
C                     ***************
C
C     STORAGE:
C
C             Y(0,0)  Y(1,1)  Y(2,2)
C             Y(1,-1) Y(1,0)  Y(2,1)   .  .  .
C             Y(2,-2) Y(2,-1) Y(2,0)
C
C     ** NOTE ** DIMENSIONS OF ARRAY YLM START FROM ZERO, NOT ONE,
C                 SO THAT:
C
C                 =  YLM(K,L-M,L)     WHEN M > 0 (UPPER TRIANGLE)
C         Y
C          LM     =  YLM(K,L,L+M)     WHEN M < 0 (LOWER TRIANGLE)
C
C     GENERATED USING RECURRENCE RELATIONS (SEE "NUMERICAL RECIPES").
C
      USE CONSTANTS
      USE NuclearMeshes, ONLY : ML
      USE SphericalHarmonics
C TODO Eliminate NuclearMeshes here by removing YLM allocation.
C TODO Consider making YLM an argument instead of a module member.
C TODO Or make YLM a function return value!
      IMPLICIT NONE
      ! Arguments
      REAL(KIND(1D0)) X(NANGS),Y(NANGS),Z(NANGS)
      INTEGER NANGS, NY
      ! Locals
      REAL(KIND(1D0)) SINTH(NANGS),SINPH(NANGS),COSPH(NANGS)
      REAL(KIND(1D0)) SINMPH(NANGS,NY),COSMPH(NANGS,NY)
      REAL(KIND(1D0)), PARAMETER :: REALI(20) =
     + (/1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,
     + 11.D0,12.D0,13.D0,14.D0,15.D0,16.D0,17.D0,18.D0,19.D0,20.D0/)
      REAL(KIND(1D0)) C,CN
      INTEGER I,J,K,L,N
C
      if (.not.allocated(YLM)) then
         allocate(YLM(ML,0:MYLM,0:MYLM))
      end if
C
      DO 10 K=1,NANGS
      YLM(K,0,0)=1.D0
      SINTH(K)=DSQRT(X(K)**2+Y(K)**2)+EPS
      SINPH(K)=Y(K)/SINTH(K)
      COSPH(K)=X(K)/SINTH(K)
      SINMPH(K,1)=SINPH(K)
10    COSMPH(K,1)=COSPH(K)
C
      DO 20 J=2,NY
      DO 20 K=1,NANGS
      SINMPH(K,J)=COSPH(K)*SINMPH(K,J-1)+SINPH(K)*COSMPH(K,J-1)
20    COSMPH(K,J)=COSPH(K)*COSMPH(K,J-1)-SINPH(K)*SINMPH(K,J-1)
C
C     FIRST ROW:
C
      C=DSQRT(2.D0)
      DO 100 J=1,NY
      DO 100 K=1,NANGS
      SINMPH(K,J)=C*SINMPH(K,J)
      COSMPH(K,J)=C*COSMPH(K,J)
100   YLM(K,0,J)=DFAC(J)*SINTH(K)**J
C
C     SECOND ROW:
C
      DO 200 I=0,NY-1
      DO 200 K=1,NANGS
200   YLM(K,1,I+1)=(2*I+1)*Z(K)*YLM(K,0,I)
C
C     REMAINING UPPER TRIANGLE:
C
      DO 300 I=2,NY
      DO 300 J=I,NY
      DO 300 K=1,NANGS
300   YLM(K,I,J)=((2*J-1)*Z(K)*YLM(K,I-1,J-1)
     +              -(2*J-I-1)*YLM(K,I-2,J-2))/REALI(I)
C
C     NORMALIZE UPPER TRIANGLE WITH RESPECT TO THETA:
C
      DO 310 L=1,NY
      CN=DSQRT(2.D0*L+1.D0)
      DO 310 K=1,NANGS
310   YLM(K,L,L)=CN*YLM(K,L,L)
      DO 320 J=1,NY
      DO 320 I=0,J-1
      CN=DSQRT((2*J+1)*FAC(I)/FAC(2*J-I))
      DO 320 K=1,NANGS
320   YLM(K,I,J)=CN*YLM(K,I,J)
C
C     COPY UPPER TRIANGLE INTO LOWER TRIANGLE (IGNORE PHASE):
C
      DO 400 I=1,NY
      DO 400 J=0,I-1
      DO 400 K=1,NANGS
400   YLM(K,I,J)=YLM(K,J,I)
C
C     MULTIPLY BY NORMALIZED PHI FACTORS:
C
      DO 500 J=1,NY
      DO 500 I=0,J-1
      N=J-I
      DO 500 K=1,NANGS
      YLM(K,I,J)=YLM(K,I,J)*SINMPH(K,N)
500   YLM(K,J,I)=YLM(K,J,I)*COSMPH(K,N)
C
      RETURN
      END
C
C
C
C$PRAGMA SUN OPT=5
      SUBROUTINE SPLINE (H,Y,A,B,C,N,IBC,BCL,BCR)
C
C     FITS A CUBIC SPLINE TO DATA POINTS ON THE UNIFORM R-MESH.
C     RETURNS COEFFICIENTS A,B,C DEFINING THE CUBIC POLYNOMIAL TO THE
C     RIGHT OF EACH KNOT (INCLUDING A(0),B(0),C(0)).
C
C     INPUT:
C         H    - INTERVAL BETWEEN KNOTS
C         Y    - ARRAY OF ORDINATE VALUES
C         N    - THE NUMBER OF DATA POINTS
C         IBC  - LEFT BOUNDARY CONDITION (R=0);
C           0  - ZERO VALUE, S-WAVE CUSP
C           1  - ZERO VALUE, ZERO 1ST DERIVATIVE
C           2  - ZERO VALUE, ZERO 2ND DERIVATIVE
C         BCL  - CUSP INFORMATION IF IBC=0
C         BCR  - VALUE OF FUNCTION AT R=INFINITY (SLOPE=0)
C
C     OUTPUT:
C         A    - ARRAY OF COEFFICIENTS OF X
C         B    - ARRAY OF COEFFICIENTS OF X**2
C         C    - ARRAY OF COEFFICIENTS OF X**3
C
      USE IOUNITS
      USE ERRORHANDLING
      USE SPLINES, ONLY : MDATA
      IMPLICIT NONE
      ! Arguments
      REAL(KIND(1D0)), INTENT(IN) :: H, Y(0:MDATA)
      REAL(KIND(1D0)), INTENT(OUT):: A(0:MDATA),B(0:MDATA),C(0:MDATA)
      INTEGER,         INTENT(IN) :: N,IBC
      REAL(KIND(1D0)), INTENT(IN) :: BCL,BCR
      ! Locals
      REAL(KIND(1D0)) :: BS(-1:MDATA+1), CUSP
      REAL(KIND(1D0)) :: Y1(0:MDATA+1),T1(MDATA),T2(MDATA),T3(MDATA)
      INTEGER I,IERR
C
      Y1(0)=Y(0)
      DO 10 I=1,N
      BS(I)=Y(I)
      Y1(I)=Y(I)
      T1(I)=1.D0
      T2(I)=4.D0
10    T3(I)=1.D0
C
      IF(IBC.EQ.0)THEN
      CUSP=(1.D0-BCL)*H
      T2(1)=4.D0-CUSP/(2.D0*CUSP+3.D0)
      ELSE IF(IBC.EQ.1)THEN
      T2(1)=3.5D0
      ELSE IF(IBC.EQ.2)THEN
      T2(1)=4.0D0
      ELSE
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** SPLINE: UNRECOGNIZED VALUE OF IBC ***'
      ERRFLG = .TRUE.
      RETURN
      END IF
      T2(N)=3.5D0
      BS(N)=BS(N)-0.25D0*BCR
C
      CALL DGTSV(N,1,T1,T2,T3,BS(1:),MDATA,IERR)
      IF(IERR.NE.0)THEN
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** SPLINE: ERROR IN ROUTINE DGTSL ***'
      ERRFLG = .TRUE.
      RETURN
      END IF
C
      IF(IBC.EQ.0)THEN
      BS(0)=-CUSP/(2.D0*CUSP+3.D0)*BS(1)
      BS(-1)=(2.D0*CUSP-3.D0)/(2.D0*CUSP+3.D0)*BS(1)
      ELSE IF(IBC.EQ.1)THEN
      BS(0)=-0.5D0*BS(1)
      BS(-1)=BS(1)
      ELSE IF(IBC.EQ.2)THEN
      BS(0)=0.D0
      BS(-1)=-BS(1)
      ELSE
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** SPLINE: UNRECOGNIZED VALUE OF IBC ***'
      ERRFLG = .TRUE.
      RETURN
      END IF
      BS(N+1)=-0.5D0*BS(N)+0.25D0*BCR
      Y1(N+1)= BCR
C
      DO 20 I=0,N
      A(I)=3.D0*(BS(I+1)-BS(I-1))/H
      B(I)=3.D0*(BS(I-1)-2.D0*BS(I)+BS(I+1))/H**2
20    C(I)= (Y1(I+1)-Y1(I)-A(I)*H-B(I)*H**2)/H**3
C
      RETURN
      END
C
C
C
      SUBROUTINE FDRHO(RHO,IORDER,RHOXYZ)
C
C     FIRST AND SECOND PARTIAL DERIVATIVES (CARTESIAN)
C     OF DENSITY RHO BY FINITE DIFFERENCES.
C
C     ON OUTPUT, RHOXYZ(I,***) CONTAINS
C     FOR I=0      :  THE MODEL DENSITY
C     FOR I=1,2,3  :  1ST DERIVS  X,Y,Z     (IF IORDER.GE.1)
C     FOR I=4,5,6  :  2ND DERIVS  XX,YY,ZZ  (IF IORDER.GE.2)
C     FOR I=7,8,9  :  2ND DERIVS  YZ,XZ,XY  (IF IORDER.EQ.3)
C
      USE ERRORHANDLING
      USE GRID
      USE NuclearFramework
      USE NuclearMeshes
      USE SIZES
      USE SphericalHarmonics
      IMPLICIT NONE
C     Arguments
      REAL(KIND(1D0)),INTENT(IN) :: RHO(NNN)
      INTEGER,        INTENT(IN) :: IORDER
      REAL(KIND(1D0)),INTENT(OUT):: RHOXYZ(0:9,NNN)
C     Locals
      REAL(KIND(1D0)) RHOLM(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) A(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) B(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) C(0:MR,0:MYLM,0:MYLM)
      INTEGER INTQS(ML),NXYZ,NDEL,I,K,ICEL,NRADS,NANGS,NWAVES,IY,JY,KKK
      INTEGER IR,IL,IDEL,INUC,INTQ,IXYZ
      REAL(KIND(1D0)) RS(ML),XUNIT(ML),YUNIT(ML),ZUNIT(ML),DQS(ML)
      REAL(KIND(1D0)) DELX(0:12),DELY(0:12),DELZ(0:12),COEF(0:9,0:12)
      REAL(KIND(1D0)) H,Q,DQ,RFAC,PRODRL,SUM,X,Y,Z
      INTEGER IDELX(0:12),IDELY(0:12),IDELZ(0:12),ICOEF(0:9,0:12)
      REAL(KIND(1D0)),PARAMETER :: FDH=1.D-04
      DATA (IDELX(I),I=0,12)/ 0,+1,-1, 0, 0, 0, 0, 0, 0,+1,-1,+1,-1/
      DATA (IDELY(I),I=0,12)/ 0, 0, 0,+1,-1, 0, 0,+1,-1, 0, 0,+1,-1/
      DATA (IDELZ(I),I=0,12)/ 0, 0, 0, 0, 0,+1,-1,+1,-1,+1,-1, 0, 0/
      DATA (COEF(0,I),I=0,12)/1.D0,12*0.D0/
      DATA (ICOEF(1,I),I=0,12)/ 0,+1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(2,I),I=0,12)/ 0, 0, 0,+1,-1, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(3,I),I=0,12)/ 0, 0, 0, 0, 0,+1,-1, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(4,I),I=0,12)/-2,+1,+1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(5,I),I=0,12)/-2, 0, 0,+1,+1, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(6,I),I=0,12)/-2, 0, 0, 0, 0,+1,+1, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(7,I),I=0,12)/-2, 0, 0,+1,+1,+1,+1,-1,-1, 0, 0, 0, 0/
      DATA (ICOEF(8,I),I=0,12)/-2,+1,+1, 0, 0,+1,+1, 0, 0,-1,-1, 0, 0/
      DATA (ICOEF(9,I),I=0,12)/-2,+1,+1,+1,+1, 0, 0, 0, 0, 0, 0,-1,-1/
C
      IF(IORDER.EQ.1)THEN
      NXYZ=3
      NDEL=6
      ELSE IF(IORDER.EQ.2)THEN
      NXYZ=6
      NDEL=6
      ELSE IF(IORDER.EQ.3)THEN
      NXYZ=9
      NDEL=12
      END IF
      DO 10 I=0,NDEL
      DELX(I)=FDH*IDELX(I)
      DELY(I)=FDH*IDELY(I)
      DELZ(I)=FDH*IDELZ(I)
      COEF(1,I)=0.5D0*ICOEF(1,I)/FDH
      COEF(2,I)=0.5D0*ICOEF(2,I)/FDH
      COEF(3,I)=0.5D0*ICOEF(3,I)/FDH
      COEF(4,I)=ICOEF(4,I)/FDH**2
      COEF(5,I)=ICOEF(5,I)/FDH**2
      COEF(6,I)=ICOEF(6,I)/FDH**2
      COEF(7,I)=-0.5D0*ICOEF(7,I)/FDH**2
      COEF(8,I)=-0.5D0*ICOEF(8,I)/FDH**2
10    COEF(9,I)=-0.5D0*ICOEF(9,I)/FDH**2
      DO 20 K=1,NNN
      DO 20 I=0,NXYZ
20    RHOXYZ(I,K)=0.D0
C
      DO 100 ICEL=1,NUCS
C
      NRADS=NR(ICEL)
      NANGS=NL(ICEL)
      NWAVES=LMAX(ICEL)
      H=1.D0/(NRADS+1)
C
C     SPHERICAL EXPANSION OF RHO:
C
      CALL YCALC(XANGS(:,ICEL),YANGS(:,ICEL),ZANGS(:,ICEL),NANGS,NWAVES)
C
      DO 200 IY=0,NWAVES
      DO 200 JY=0,NWAVES
C
      RHOLM(0,IY,JY)=0.D0
      KKK=KKKTAB(ICEL)
      DO 210 IR=1,NRADS
      SUM=0.D0
      DO 211 IL=1,NANGS
211   SUM=SUM+WANGS(IL,ICEL)*YLM(IL,IY,JY)*WNUC(KKK+IL)*RHO(KKK+IL)
      RHOLM(IR,IY,JY)=RADS(IR,ICEL)**2*SUM
210   KKK=KKK+NANGS
C
      CALL SPLINE(H,RHOLM(0,IY,JY),A(0,IY,JY),B(0,IY,JY),C(0,IY,JY),
     +                                            NRADS,1,0.D0,0.D0)
      if (ERRFLG) return
200   CONTINUE
C
      DO 12 IDEL=0,NDEL
C
      KKK=0
      DO 300 INUC=1,NUCS
C
      DO 310 IR=1,NR(INUC)
      DO 320 IL=1,NL(INUC)
C     CALC COORDS ON UNIT SPHERE ABOUT ICEL:
      X=XNUC(INUC)+XANGS(IL,INUC)*RADS(IR,INUC)-XNUC(ICEL)+DELX(IDEL)
      Y=YNUC(INUC)+YANGS(IL,INUC)*RADS(IR,INUC)-YNUC(ICEL)+DELY(IDEL)
      Z=ZNUC(INUC)+ZANGS(IL,INUC)*RADS(IR,INUC)-ZNUC(ICEL)+DELZ(IDEL)
      RS(IL)=DSQRT(X**2+Y**2+Z**2)
      XUNIT(IL)=X/RS(IL)
      YUNIT(IL)=Y/RS(IL)
      ZUNIT(IL)=Z/RS(IL)
      Q=RS(IL)/(RS(IL)+RMID(ICEL))
      INTQS(IL)=IDINT((NRADS+1)*Q)
      DQS(IL)=Q-INTQS(IL)*H
320   CONTINUE
      CALL YCALC(XUNIT,YUNIT,ZUNIT,NL(INUC),NWAVES)
      DO 330 IY=0,NWAVES
      DO 330 JY=0,NWAVES
      DO 330 IL=1,NL(INUC)
      INTQ=INTQS(IL)
      DQ=DQS(IL)
      Y=RHOLM(INTQ,IY,JY)
      RFAC=Y+DQ*(A(INTQ,IY,JY)+DQ*(B(INTQ,IY,JY)+DQ*C(INTQ,IY,JY)))
      PRODRL=RFAC*YLM(IL,IY,JY)/RS(IL)**2
      DO 330 IXYZ=0,NXYZ
      RHOXYZ(IXYZ,KKK+IL)=RHOXYZ(IXYZ,KKK+IL)+COEF(IXYZ,IDEL)*PRODRL
330   CONTINUE
      KKK=KKK+NL(INUC)
310   CONTINUE
C
300   CONTINUE
C
12    CONTINUE
C
100   CONTINUE
C
      RETURN
      END
C
C
C
C$PRAGMA SUN OPT=2
      SUBROUTINE SCHROD(V,VX,EN,RHS,RES0,RES2,RESDC,RNORM)
C
C     APPROXIMATE CALCULATION OF DELTA-PSI.
C     ASSUMES INDEPENDENT SINGLE-CENTRE EQUATIONS, EACH WITH
C     RHS(N)=WNUC(N)*RHS AND A SPHERICALIZED MODEL POTENTIAL VSPH(R).
C
C     RADIAL EQTNS SOLVED BY MINIMUM "SQUARED-RESIDUAL + NORM" SCHEME.
C
      USE CONSTANTS
      USE ERRORHANDLING
      USE GRID
      USE IOUNITS
      USE NuclearFramework
      USE NuclearMeshes
      USE SCFConvergence
      USE SecondDerivativeMatrices
      USE SIZES
      USE SphericalHarmonics
      IMPLICIT NONE
C TODO Refactor the interpolant storage here -- can all this be
C TODO locally allocated and deallocated?
C TODO Beware also of name conflicts with module SphericalInterpolants
C     Arguments
      REAL(KIND(1D0)),INTENT(IN) :: V(NNN),VX(NNN),EN,RHS(NNN)
      REAL(KIND(1D0)),INTENT(OUT):: RES0(NNN),RES2(NNN),RESDC(NNN),RNORM
C     Locals
      REAL(KIND(1D0)) RES0LM(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) RES2LM(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) A0(0:MR,0:MYLM,0:MYLM),A2(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) B0(0:MR,0:MYLM,0:MYLM),B2(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) C0(0:MR,0:MYLM,0:MYLM),C2(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) VSPH(MR),W(0:MR+1),C(0:MR+1),D(0:MR+1),E(0:MR+1)
      REAL(KIND(1D0)) RS(ML),XUNIT(ML),YUNIT(ML),ZUNIT(ML),DQS(ML)
      REAL(KIND(1D0)) XX(MR),AA(3,MR)
      REAL(KIND(1D0)) H,SUM,R,ZRMD,X,Y,Z,Q,DQ,Y0,Y2,RFAC0,RFAC2
      INTEGER INTQS(ML),I,ICEL,NRADS,NANGS,NWAVES,KKK,IR,IL,IY,JY,L,L2
      INTEGER IERR,K,II,INUC,INTQ
C
      C(0) = 0.D0
      D(0) = 0.D0
      E(0) = 0.D0
      W(0) = 0.D0
      C(1) = 0.D0
      DO 10 I=1,NNN
      RES0(I) = 0.D0
      RES2(I) = 0.D0
      RESDC(I) = 0.D0
10    CONTINUE
C
      DO 100 ICEL=1,NUCS
C
      NRADS=NR(ICEL)
      NANGS=NL(ICEL)
      NWAVES=NY(ICEL)
      E(NRADS)=0.D0
      C(NRADS+1)=0.D0
      D(NRADS+1)=0.D0
      E(NRADS+1)=0.D0
      W(NRADS+1)=0.D0
      H=1.D0/(NRADS+1)
C
      KKK=KKKTAB(ICEL)
      DO 101 IR=1,NRADS
      VSPH(IR)=0.D0
      W(IR)=WRADS(IR,ICEL)/RADS(IR,ICEL)**2
      DO 102 IL=1,NANGS
102   VSPH(IR)=VSPH(IR)+WANGS(IL,ICEL)*WNUC(KKK+IL)*
     +                        (V(KKK+IL)+VX(KKK+IL))
101   KKK=KKK+NANGS
C
C     SPHERICAL EXPANSION OF WNUC*RHS:
C
      CALL YCALC(XANGS(:,ICEL),YANGS(:,ICEL),ZANGS(:,ICEL),NANGS,NWAVES)
C
      DO 200 IY=0,NWAVES
      DO 200 JY=0,NWAVES
C
      KKK=KKKTAB(ICEL)
      DO 210 IR=1,NRADS
      SUM=0.D0
      DO 211 IL=1,NANGS
211   SUM=SUM+WANGS(IL,ICEL)*YLM(IL,IY,JY)*WNUC(KKK+IL)*RHS(KKK+IL)
      RES0LM(IR,IY,JY)=RADS(IR,ICEL)*SUM
210   KKK=KKK+NANGS
C
C     SOLVE RADIAL LM EQUATION:
C
      L=MAX0(IY,JY)
      L2=L*(L+1)
C
      D(1)=-0.5D0*TRID2R(1,2,ICEL)+0.5D0*L2/RADS(1,ICEL)**2+VSPH(1)-EN
      E(1)=-0.5D0*TRID2R(1,3,ICEL)
      DO 220 I=2,NRADS-1
      C(I)=-0.5D0*TRID2R(I,1,ICEL)
      D(I)=-0.5D0*TRID2R(I,2,ICEL)+0.5D0*L2/RADS(I,ICEL)**2+VSPH(I)-EN
220   E(I)=-0.5D0*TRID2R(I,3,ICEL)
      C(NRADS)=-0.5D0*TRID2R(NRADS,1,ICEL)
      D(NRADS)=-0.5D0*TRID2R(NRADS,2,ICEL)
     +                    +0.5D0*L2/RADS(NRADS,ICEL)**2+VSPH(NRADS)-EN
C
      DO 230 I=1,NRADS
      AA(1,I)=C(I-1)*W(I-1)*E(I-1)
      AA(2,I)=D(I-1)*W(I-1)*E(I-1)+C(I)*W(I)*D(I)
230   AA(3,I)=W(I-1)*E(I-1)**2+W(I)*D(I)**2+W(I+1)*C(I+1)**2 + SING*W(I)
C
      DO 240 I=1,NRADS
240   XX(I) =  E(I-1)*W(I-1)*RES0LM(I-1,IY,JY)
     +        +D(I)*W(I)*RES0LM(I,IY,JY)
     +        +C(I+1)*W(I+1)*RES0LM(I+1,IY,JY)
C
      CALL DPBSV('U',NRADS,2,1,AA,3,XX,MR,IERR)
      IF(IERR.NE.0)THEN
      WRITE(IOUT,900) IERR,'DPBSV'
  900 FORMAT(/'*** SCHROD: ERROR ',I3,' IN ROUTINE ',A6,' ***')
      ERRFLG = .TRUE.
      RETURN
      END IF
C
      RES0LM(0,IY,JY)=0.D0
      DO 250 IR=1,NRADS
250   RES0LM(IR,IY,JY)=XX(IR)
C
      RES2LM(0,IY,JY)=0.D0
      DO 260 IR=1,NRADS
      SUM=0.D0
      DO 261 K=1,7
      II=IR+K-4
      IF(II.LT.1.OR.II.GT.NRADS)GO TO 261
      SUM=SUM+D2R(IR,K,ICEL)*RES0LM(II,IY,JY)
261   CONTINUE
260   RES2LM(IR,IY,JY)=RADS(IR,ICEL)*SUM
     +                 -L2/RADS(IR,ICEL)*RES0LM(IR,IY,JY)
C
      KKK=KKKTAB(ICEL)
      DO 280 IR=1,NRADS
      R=RADS(IR,ICEL)
      DO 281 IL=1,NANGS
      RES0(KKK+IL)=RES0(KKK+IL)
     +             +RES0LM(IR,IY,JY)*YLM(IL,IY,JY)/R
      RES2(KKK+IL)=RES2(KKK+IL)
     +             +RES2LM(IR,IY,JY)*YLM(IL,IY,JY)/R**2
281   RESDC(KKK+IL)=RESDC(KKK+IL)
     +             +RES0LM(IR,IY,JY)*YLM(IL,IY,JY)/R
280   KKK=KKK+NANGS
C
      IF(L.EQ.0)THEN
      ZRMD=DFLOAT(NUCZ(ICEL))*RMID(ICEL)
      CALL SPLINE(H,RES0LM(0,IY,JY),A0(0,IY,JY),B0(0,IY,JY),C0(0,IY,JY),
     +                                                NRADS,0,ZRMD,0.D0)
      CALL SPLINE(H,RES2LM(0,IY,JY),A2(0,IY,JY),B2(0,IY,JY),C2(0,IY,JY),
     +                                                NRADS,2,0.D0,0.D0)
      ELSE
      CALL SPLINE(H,RES0LM(0,IY,JY),A0(0,IY,JY),B0(0,IY,JY),C0(0,IY,JY),
     +                                                NRADS,1,0.D0,0.D0)
      CALL SPLINE(H,RES2LM(0,IY,JY),A2(0,IY,JY),B2(0,IY,JY),C2(0,IY,JY),
     +                                                NRADS,1,0.D0,0.D0)
      END IF
      if (ERRFLG) return
200   CONTINUE
C
      DO 300 INUC=1,NUCS
      IF(INUC.EQ.ICEL)GO TO 300
C
      KKK=KKKTAB(INUC)
      DO 310 IR=1,NR(INUC)
      DO 320 IL=1,NL(INUC)
C     CALC COORDS ON UNIT SPHERE ABOUT ICEL:
      X=XNUC(INUC)+XANGS(IL,INUC)*RADS(IR,INUC)-XNUC(ICEL)
      Y=YNUC(INUC)+YANGS(IL,INUC)*RADS(IR,INUC)-YNUC(ICEL)
      Z=ZNUC(INUC)+ZANGS(IL,INUC)*RADS(IR,INUC)-ZNUC(ICEL)
      RS(IL)=DSQRT(X**2+Y**2+Z**2)
      XUNIT(IL)=X/RS(IL)
      YUNIT(IL)=Y/RS(IL)
      ZUNIT(IL)=Z/RS(IL)
      Q=RS(IL)/(RS(IL)+RMID(ICEL))
      INTQS(IL)=IDINT((NRADS+1)*Q)
      DQS(IL)=Q-INTQS(IL)*H
320   CONTINUE
      CALL YCALC(XUNIT,YUNIT,ZUNIT,NL(INUC),NWAVES)
      DO 330 IY=0,NWAVES
      DO 330 JY=0,NWAVES
      DO 330 IL=1,NL(INUC)
      INTQ=INTQS(IL)
      DQ=DQS(IL)
      Y0=RES0LM(INTQ,IY,JY)
      Y2=RES2LM(INTQ,IY,JY)
      RFAC0=Y0+DQ*(A0(INTQ,IY,JY)+DQ*(B0(INTQ,IY,JY)+DQ*C0(INTQ,IY,JY)))
      RFAC2=Y2+DQ*(A2(INTQ,IY,JY)+DQ*(B2(INTQ,IY,JY)+DQ*C2(INTQ,IY,JY)))
      RES0(KKK+IL)=RES0(KKK+IL)+RFAC0*YLM(IL,IY,JY)/RS(IL)
330   RES2(KKK+IL)=RES2(KKK+IL)+RFAC2*YLM(IL,IY,JY)/RS(IL)**2
310   KKK=KKK+NL(INUC)
C
300   CONTINUE
C
100   CONTINUE
C
      RNORM=0.D0
      DO 400 KKK=1,NNN
400   RNORM=RNORM+WINTS(KKK)*RES0(KKK)**2
      RNORM=DSQRT(RNORM)
C
      RETURN
      END
C
C
C
      SUBROUTINE DISPER(RHOS,XXDINT,XXDIP2)
C
      USE Constants
      USE Grid
      USE IOUnits
      USE KeyWords
      USE MPI
      USE OccupationNumbers
      USE Sizes
      IMPLICIT NONE
      ! Arguments
      REAL(KIND(1D0)) :: RHOS(NNN,2),XXDINT,XXDIP2(NNN)
      ! Locals
      REAL(KIND(1D0)) :: PSII(NNN),PSIJ(NNN)
      REAL(KIND(1D0)) :: DIPX(NNN),DIPY(NNN),DIPZ(NNN)
      REAL(KIND(1D0)) :: DIPX_sum(NNN),DIPY_sum(NNN),DIPZ_sum(NNN)
      INTEGER ISPIN,KKK,IO,IMORB,JMORB,I,IREC,JREC
      REAL(KIND(1D0)) FACTOR,WPROD,XINT,YINT,ZINT,FPROD,SDPOL2
C
      XXDINT=0.D0
      if(myRank.eq.0)then
      DO 1 KKK=1,NNN
1     XXDIP2(KKK)=0.D0
      end if ! myRank==0
C
      DO 1000 ISPIN=1,NSPINS
C
      DO 3 KKK=1,NNN
      DIPX(KKK)=0.D0
      DIPY(KKK)=0.D0
3     DIPZ(KKK)=0.D0
C
      IF(ISPIN.EQ.1)IO=MOSA0
      IF(ISPIN.EQ.2)IO=MOSB0
      DO 10 IMORB=1,NMORBS
        IF (myRank.EQ.IOWNER(IMORB)) THEN
          IREC=MOSEQ(IMORB)
          READ(IO,REC=IREC)(PSII(I),I=1,NNN)
        END IF ! myRank==iOwner
        CALL MPI_BCAST(PSII,NNN,MPI_DOUBLE_PRECISION,
     +                 IOWNER(IMORB),MPI_COMM_WORLD,iError)
        DO 20 JMORB=1,IMORB
          IF (myRank.EQ.IOWNER(JMORB)) THEN
            JREC=MOSEQ(JMORB)
            READ(IO,REC=JREC)(PSIJ(I),I=1,NNN)
            FACTOR=2.D0*OCCS(IMORB,ISPIN)*OCCS(JMORB,ISPIN)
            IF(JMORB.EQ.IMORB)FACTOR=0.5D0*FACTOR
            XINT=0.D0
            YINT=0.D0
            ZINT=0.D0
            DO 200 KKK=1,NNN
              WPROD=WINTS(KKK)*PSII(KKK)*PSIJ(KKK)
              XINT=XINT+XMESH(KKK)*WPROD
              YINT=YINT+YMESH(KKK)*WPROD
              ZINT=ZINT+ZMESH(KKK)*WPROD
200         CONTINUE
            DO 100 KKK=1,NNN
              FPROD=FACTOR*PSII(KKK)*PSIJ(KKK)/(RHOS(KKK,ISPIN)+SMALL)
              DIPX(KKK)=DIPX(KKK)+FPROD*XINT
              DIPY(KKK)=DIPY(KKK)+FPROD*YINT
              DIPZ(KKK)=DIPZ(KKK)+FPROD*ZINT
100         CONTINUE
          END IF ! myRank==iOwner
20      CONTINUE ! jmorb
10    CONTINUE   ! imorb
C
C     Sum up DIPX/Y/Z on root
      CALL MPI_REDUCE(DIPX, DIPX_sum, NNN, MPI_DOUBLE_PRECISION,
     +                MPI_SUM, 0, MPI_COMM_WORLD, iError)
      CALL MPI_REDUCE(DIPY, DIPY_sum, NNN, MPI_DOUBLE_PRECISION,
     +                MPI_SUM, 0, MPI_COMM_WORLD, iError)
      CALL MPI_REDUCE(DIPZ, DIPZ_sum, NNN, MPI_DOUBLE_PRECISION,
     +                MPI_SUM, 0, MPI_COMM_WORLD, iError)
C
      if (myRank.eq.0) then
      DO 300 KKK=1,NNN
      SDPOL2=RHOS(KKK,ISPIN)*
     +           ((DIPX_sum(KKK)-XMESH(KKK))**2
     +           +(DIPY_sum(KKK)-YMESH(KKK))**2
     +           +(DIPZ_sum(KKK)-ZMESH(KKK))**2)
      XXDINT=XXDINT+WINTS(KKK)*SDPOL2
300   XXDIP2(KKK)=XXDIP2(KKK)+SDPOL2
      end if ! myRank==0
C
1000  CONTINUE
C
      if (myRank.eq.0) then
      IF(NSPINS.EQ.1)THEN
      XXDINT=2.D0*XXDINT
      DO 1001 KKK=1,NNN
1001  XXDIP2(KKK)=2.D0*XXDIP2(KKK)
      END IF ! nSpins==1
      end if ! myRank==0
C
      RETURN
      END
C
C
C
      SUBROUTINE FORCES(RHO,V,FORCE0,FORCE,DIPOLE)
C
C     HELLMANN-FEYNMAN FORCES.
C
      USE CORES
      USE GRID
      USE IOUNITS
      USE NuclearFramework
      USE NuclearMeshes
      USE PROMOLECULE
      USE SIZES
      IMPLICIT NONE
C     Arguments
      REAL(KIND(1D0)) RHO(NNN),V(NNN,2),DIPOLE(3)
      REAL(KIND(1D0)) FORCE0(NUCS,3),FORCE(NUCS,3)
C     Locals
      REAL(KIND(1D0)) FRZCOR(3)
      INTEGER I,INUC,NRADS,NANGS,KKK,IR,IL
      REAL(KIND(1D0)) CHARGE,DCORE,FUNC,FX,FY,FZ,R,R3

C     SUBTRACT REFERENCE DENSITY
C     AND CALCULATE DIPOLE MOMENT:
      DIPOLE(1)=0.D0
      DIPOLE(2)=0.D0
      DIPOLE(3)=0.D0
      DO 10 I=1,NNN
      RHO(I)=RHO(I)-RHO0(I)
      DIPOLE(1)=DIPOLE(1)-WINTS(I)*XMESH(I)*RHO(I)
      DIPOLE(2)=DIPOLE(2)-WINTS(I)*YMESH(I)*RHO(I)
10    DIPOLE(3)=DIPOLE(3)-WINTS(I)*ZMESH(I)*RHO(I)
C
      DO 100 INUC=1,NUCS
      NRADS=NR(INUC)
      NANGS=NL(INUC)
      CHARGE=DFLOAT(NUCZ(INUC))
C
C     FROZEN-CORE CORRECTION:
      FRZCOR(1)=0.D0
      FRZCOR(2)=0.D0
      FRZCOR(3)=0.D0
      KKK=KKKTAB(INUC)
      DO 110 IR=1,NRADS
      DCORE=DCOR(IR,INUC)
      DO 111 IL=1,NANGS
      FUNC=WRADS(IR,INUC)*WANGS(IL,INUC)
     +     *DCORE*0.5D0*(V(KKK+IL,1)+V(KKK+IL,2))
      FRZCOR(1)=FRZCOR(1)+XANGS(IL,INUC)*FUNC
      FRZCOR(2)=FRZCOR(2)+YANGS(IL,INUC)*FUNC
      FRZCOR(3)=FRZCOR(3)+ZANGS(IL,INUC)*FUNC
111   CONTINUE
      KKK=KKK+NANGS
110   CONTINUE
C
      FX=0.D0
      FY=0.D0
      FZ=0.D0
      DO 200 I=1,NNN
      R=DSQRT((XMESH(I)-XNUC(INUC))**2
     +       +(YMESH(I)-YNUC(INUC))**2
     +       +(ZMESH(I)-ZNUC(INUC))**2)
      R3=R**3
      FX = FX + WINTS(I)*RHO(I)*(XMESH(I)-XNUC(INUC))/R3
      FY = FY + WINTS(I)*RHO(I)*(YMESH(I)-YNUC(INUC))/R3
      FZ = FZ + WINTS(I)*RHO(I)*(ZMESH(I)-ZNUC(INUC))/R3
200   CONTINUE
      FORCE(INUC,1) = FORCE0(INUC,1) + CHARGE*FX + FRZCOR(1)
      FORCE(INUC,2) = FORCE0(INUC,2) + CHARGE*FY + FRZCOR(2)
      FORCE(INUC,3) = FORCE0(INUC,3) + CHARGE*FZ + FRZCOR(3)
100   CONTINUE
C
C     RESTORE INPUT DENSITY:
      DO 300 I=1,NNN
300   RHO(I)=RHO(I)+RHO0(I)
C
      RETURN
      END
C
C
C
C$PRAGMA SUN OPT=2
      SUBROUTINE MIXRES(V,E,PSI0,PSI2,PSIDC,RES0,RES2,RESDC,RATIO,
     &                  XPSI,XRES)
C
C     CALCULATES LINEAR COMBINATION OF PSI AND RES
C     THAT MINIMIZES RESIDUAL FOR "FULL" V AND ENERGY E.
C
      USE ERRORHANDLING
      USE GRID
      USE IOUNITS
      USE SIZES
      IMPLICIT NONE
      REAL(KIND(1D0)),INTENT(IN)   :: V(NNN),E,RES0(NNN)
      REAL(KIND(1D0)),INTENT(IN)   :: RES2(NNN),RESDC(NNN),XRES(NNN)
      REAL(KIND(1D0)),INTENT(INOUT):: PSI0(NNN),PSI2(NNN),PSIDC(NNN)
      REAL(KIND(1D0)),INTENT(INOUT):: XPSI(NNN)
      REAL(KIND(1D0)),INTENT(OUT)  :: RATIO
      REAL(KIND(1D0)) R(2,2),S(2,2),W(2),WORK(5)
      REAL(KIND(1D0)) S1,S2,R1,R2,WINT,CPSI,CRES
      INTEGER KKK,IERR
C
      S(1,1)=0.D0
      S(1,2)=0.D0
      S(2,2)=0.D0
      R(1,1)=0.D0
      R(1,2)=0.D0
      R(2,2)=0.D0
      DO 100 KKK=1,NNN
      S1=PSI0(KKK)
      S2=RES0(KKK)
      R1=-0.5D0*PSI2(KKK)+(V(KKK)-E)*PSI0(KKK)+XPSI(KKK)
      R2=-0.5D0*RES2(KKK)+(V(KKK)-E)*RES0(KKK)+XRES(KKK)
      WINT=WINTS(KKK)
      S(1,1)=S(1,1)+WINT*S1*S1
      S(1,2)=S(1,2)+WINT*S1*S2
      S(2,2)=S(2,2)+WINT*S2*S2
      R(1,1)=R(1,1)+WINT*R1*R1
      R(1,2)=R(1,2)+WINT*R1*R2
100   R(2,2)=R(2,2)+WINT*R2*R2
C
C     TODO revisit dimensions of WORK -- see 'man dsygv'
      CALL DSYGV(1,'V','U',2,R,2,S,2,W,WORK,5,IERR)
      IF(IERR.NE.0)THEN
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** MIXRES: ERROR IN ROUTINE DSYGV ***'
      ERRFLG = .TRUE.
      RETURN
      END IF
C
      IF (R(1,1).GT.0.D0) THEN
        CPSI=R(1,1)
        CRES=R(2,1)
      ELSE
        CPSI=-R(1,1)
        CRES=-R(2,1)
      END IF
      RATIO=CRES/CPSI
C
      DO 200 KKK=1,NNN
      PSI0(KKK)=CPSI*PSI0(KKK)+CRES*RES0(KKK)
      PSI2(KKK)=CPSI*PSI2(KKK)+CRES*RES2(KKK)
      PSIDC(KKK)=CPSI*PSIDC(KKK)+CRES*RESDC(KKK)
200   XPSI(KKK)=CPSI*XPSI(KKK)+CRES*XRES(KKK)
C
      RETURN
      END
C
C
C
      SUBROUTINE MOLIST(VECS)
C
      USE AtomicData
      USE CORES
      USE IOUNITS
      USE NuclearFramework
      USE SIZES
      IMPLICIT NONE
      ! Argument
      REAL(KIND(1D0)) VECS(NBASIS,NBASIS)
      ! Locals
      CHARACTER*12 OSYMB(0:3,7)
      DATA OSYMB(0,1)/'S           '/
      DATA OSYMB(1,1)/'P  (Z)      '/
      DATA OSYMB(1,2)/'P  (X)      '/
      DATA OSYMB(1,3)/'P  (Y)      '/
      DATA OSYMB(2,1)/'D  (Z2)     '/
      DATA OSYMB(2,2)/'D  (XZ)     '/
      DATA OSYMB(2,3)/'D  (YZ)     '/
      DATA OSYMB(2,4)/'D  (X2-Y2)  '/
      DATA OSYMB(2,5)/'D  (XY)     '/
      DATA OSYMB(3,1)/'F  (Z3)     '/
      DATA OSYMB(3,2)/'F  (XZ2)    '/
      DATA OSYMB(3,3)/'F  (YZ2)    '/
      DATA OSYMB(3,4)/'F  Z(X2-Y2) '/
      DATA OSYMB(3,5)/'F  Z(XY)    '/
      DATA OSYMB(3,6)/'F  X(X2-3Y2)'/
      DATA OSYMB(3,7)/'F  Y(3X2-Y2)'/
      INTEGER IMORB,IAORB,IFIRST,IATM,L,LDEG,IS,M
C
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'LCAO COEFFICIENTS:'
      DO 100 IMORB=1,MIN(NMORBS+3,NBASIS)
      IF(IMORB.EQ.NMORBS+1)THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'LUMOS:'
      END IF
      WRITE(IOUT,*)' '
      WRITE(IOUT,81)IMORB
  81  FORMAT(1X,'MOLECULAR ORBITAL',I4)
      IAORB=0
      IFIRST=1
      DO 100 IATM=1,NUCS
      DO 110 L=0,LMAX0(IATM)
      LDEG=2*L+1
      DO 110 IS=1,NORBS0(IATM,L)
      IF(ICORE(IATM,L,IS).EQ.1)GO TO 110
      DO 120 M=1,LDEG
      IAORB=IAORB+1
      IF(IFIRST.EQ.1)THEN
      WRITE(IOUT,82)IATM,ZSYMB(IATM),L+IS,OSYMB(L,M),VECS(IAORB,IMORB)
  82  FORMAT(1X,'ATOM',I6,':',4X,A2,I4,A12,F16.6)
      ELSE
      WRITE(IOUT,83)     ZSYMB(IATM),L+IS,OSYMB(L,M),VECS(IAORB,IMORB)
  83  FORMAT(1X,'    ',6X,' ',4X,A2,I4,A12,F16.6)
      END IF
      IFIRST=0
120   CONTINUE
110   CONTINUE
      IFIRST=1
100   CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE FCLEAN(FORCE)
C
C     ADJUST FORCES TO SATISFY TRANSLATIONAL AND ROTATIONAL
C     INVARIANCE CONSTRAINTS.
C
      USE ERRORHANDLING
      USE IOUNITS
      USE NuclearFramework
      IMPLICIT NONE
C     Argument
      REAL(KIND(1D0)) FORCE(NUCS,3)
C     Locals
      INTEGER XZERO,YZERO,ZZERO
      REAL(KIND(1D0)) ERRBAR(NUCS),WI(NUCS*3),C(6,NUCS*3)
      REAL(KIND(1D0)) A(6,6),B(6),X(NUCS*3)
      REAL(KIND(1D0)) X2SUM,Y2SUM,Z2SUM,SUM
      INTEGER I,N3,NZEROS,NC,J0,INUC,J,K,IXYZ,IERR
C
      DO 10 I=1,NUCS
10    ERRBAR(I)=DFLOAT(NUCZ(I))
C
      N3=3*NUCS
      DO 20 I=1,6
      B(I)=0.D0
      DO 20 J=1,N3
20    C(I,J)=0.D0
C
C     MOLECULE LINEAR OR PLANAR ?:
C
      XZERO=0
      YZERO=0
      ZZERO=0
      X2SUM=0.D0
      Y2SUM=0.D0
      Z2SUM=0.D0
      DO 30 I=1,NUCS
      X2SUM=X2SUM+XNUC(I)**2
      Y2SUM=Y2SUM+YNUC(I)**2
30    Z2SUM=Z2SUM+ZNUC(I)**2
      IF(X2SUM.LT.1.D-06)XZERO=1
      IF(Y2SUM.LT.1.D-06)YZERO=1
      IF(Z2SUM.LT.1.D-06)ZZERO=1
      NZEROS = XZERO+YZERO+ZZERO
C
      IF(NZEROS.EQ.3)THEN
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** FCLEAN: NZEROS=3! CHECK INPUT DATA! ***'
      ERRFLG = .TRUE.
      RETURN
C
      ELSE IF(NZEROS.EQ.2)THEN
      NC=3
      DO 100 INUC=1,NUCS
      J0=3*(INUC-1)
      C(1,J0+1)=1.D0
      C(2,J0+2)=1.D0
      C(3,J0+3)=1.D0
      B(1)=B(1)+FORCE(INUC,1)
      B(2)=B(2)+FORCE(INUC,2)
      B(3)=B(3)+FORCE(INUC,3)
      WI(J0+1)=0.5D0*(ERRBAR(INUC))**2
      WI(J0+2)=0.5D0*(ERRBAR(INUC))**2
100   WI(J0+3)=0.5D0*(ERRBAR(INUC))**2
C
      ELSE IF(NZEROS.EQ.1.AND.XZERO.EQ.1)THEN
      NC=4
      DO 110 INUC=1,NUCS
      J0=3*(INUC-1)
      C(1,J0+1)=1.D0
      C(2,J0+2)=1.D0
      C(3,J0+3)=1.D0
      C(4,J0+3)=YNUC(INUC)
      C(4,J0+2)=-ZNUC(INUC)
      WI(J0+1)=0.5D0*(ERRBAR(INUC))**2
      WI(J0+2)=0.5D0*(ERRBAR(INUC))**2
      WI(J0+3)=0.5D0*(ERRBAR(INUC))**2
      B(1)=B(1)+FORCE(INUC,1)
      B(2)=B(2)+FORCE(INUC,2)
      B(3)=B(3)+FORCE(INUC,3)
110   B(4)=B(4)+YNUC(INUC)*FORCE(INUC,3)-ZNUC(INUC)*FORCE(INUC,2)
C
      ELSE IF(NZEROS.EQ.1.AND.YZERO.EQ.1)THEN
      NC=4
      DO 120 INUC=1,NUCS
      J0=3*(INUC-1)
      C(1,J0+1)=1.D0
      C(2,J0+2)=1.D0
      C(3,J0+3)=1.D0
      C(4,J0+1)=ZNUC(INUC)
      C(4,J0+3)=-XNUC(INUC)
      WI(J0+1)=0.5D0*(ERRBAR(INUC))**2
      WI(J0+2)=0.5D0*(ERRBAR(INUC))**2
      WI(J0+3)=0.5D0*(ERRBAR(INUC))**2
      B(1)=B(1)+FORCE(INUC,1)
      B(2)=B(2)+FORCE(INUC,2)
      B(3)=B(3)+FORCE(INUC,3)
120   B(4)=B(4)+ZNUC(INUC)*FORCE(INUC,1)-XNUC(INUC)*FORCE(INUC,3)
C
      ELSE IF(NZEROS.EQ.1.AND.ZZERO.EQ.1)THEN
      NC=4
      DO 130 INUC=1,NUCS
      J0=3*(INUC-1)
      C(1,J0+1)=1.D0
      C(2,J0+2)=1.D0
      C(3,J0+3)=1.D0
      C(4,J0+2)=XNUC(INUC)
      C(4,J0+1)=-YNUC(INUC)
      WI(J0+1)=0.5D0*(ERRBAR(INUC))**2
      WI(J0+2)=0.5D0*(ERRBAR(INUC))**2
      WI(J0+3)=0.5D0*(ERRBAR(INUC))**2
      B(1)=B(1)+FORCE(INUC,1)
      B(2)=B(2)+FORCE(INUC,2)
      B(3)=B(3)+FORCE(INUC,3)
130   B(4)=B(4)+XNUC(INUC)*FORCE(INUC,2)-YNUC(INUC)*FORCE(INUC,1)
C
      ELSE
      NC=6
      DO 140 INUC=1,NUCS
      J0=3*(INUC-1)
      C(1,J0+1)=1.D0
      C(2,J0+2)=1.D0
      C(3,J0+3)=1.D0
      C(4,J0+3)=YNUC(INUC)
      C(4,J0+2)=-ZNUC(INUC)
      C(5,J0+1)=ZNUC(INUC)
      C(5,J0+3)=-XNUC(INUC)
      C(6,J0+2)=XNUC(INUC)
      C(6,J0+1)=-YNUC(INUC)
      WI(J0+1)=0.5D0*(ERRBAR(INUC))**2
      WI(J0+2)=0.5D0*(ERRBAR(INUC))**2
      WI(J0+3)=0.5D0*(ERRBAR(INUC))**2
      B(1)=B(1)+FORCE(INUC,1)
      B(2)=B(2)+FORCE(INUC,2)
      B(3)=B(3)+FORCE(INUC,3)
      B(4)=B(4)+YNUC(INUC)*FORCE(INUC,3)-ZNUC(INUC)*FORCE(INUC,2)
      B(5)=B(5)+ZNUC(INUC)*FORCE(INUC,1)-XNUC(INUC)*FORCE(INUC,3)
140   B(6)=B(6)+XNUC(INUC)*FORCE(INUC,2)-YNUC(INUC)*FORCE(INUC,1)
C
      END IF
C
      DO 200 I=1,NC
      DO 200 J=1,NC
      A(I,J)=0.D0
      DO 200 K=1,N3
200   A(I,J)=A(I,J)+C(I,K)*C(J,K)*WI(K)
C
      CALL DPOSV('U',NC,1,A,6,B,6,IERR)
      IF(IERR.NE.0)THEN
      WRITE(IOUT,'(/''*** FCLEAN: ERROR IN ROUTINE DPOSV ***'')')
      RETURN
      END IF

C
      DO 300 I=1,N3
      SUM=0.D0
      DO 310 K=1,NC
310   SUM=SUM+C(K,I)*B(K)
300   X(I)=-WI(I)*SUM
C
      I=0
      DO 400 INUC=1,NUCS
      DO 400 IXYZ=1,3
      I=I+1
400   FORCE(INUC,IXYZ)=FORCE(INUC,IXYZ)+X(I)
C
      RETURN
      END
C
C
C
      SUBROUTINE HSSIAN
C
C     TFD MODEL HESSIAN (CRUDE BUT FAST!).
C
      USE CONSTANTS
      USE ERRORHANDLING
      USE GRID
      USE HessianMatrix
      USE IOUNITS
      USE NuclearFramework
      USE SIZES
      USE SPLINES
      IMPLICIT NONE
      REAL(KIND(1D0)) XYZNN(3),XYZ1(3),XYZ2(3)
      REAL(KIND(1D0)) F1(0:MDATA),A1(0:MDATA),B1(0:MDATA),C1(0:MDATA)
      REAL(KIND(1D0)) F2(0:MDATA),A2(0:MDATA),B2(0:MDATA),C2(0:MDATA)
      REAL(KIND(1D0)) RHOM(NNN)
      REAL(KIND(1D0)) DARHO1(NNN),DAVEL1(NNN),DARHO2(NNN),DAVEL2(NNN)
      REAL(KIND(1D0)) CX,CX2,CK,CK2,H,CHRG,RMIDZ,ZRMD,X,Y,Z,R,Q,DRQ
      REAL(KIND(1D0)) DQ,ARHO,AVEL,DARHO,DAVEL,Z1,Z2,RNN,HESSNN,R1,R2
      REAL(KIND(1D0)) DVNUC1,DVNUC2,DRHO1,DRHO2,DVEL1,DVEL2
      INTEGER I,J,N3,NDATA,NDATA1,KKK,INTQ,INUC,JNUC,II,JJ,IXYZ,JXYZ
      INTEGER IATM,IBLK,JBLK,I0,J0
C
      F1(0)=0.D0; F2(0)=0.D0
      CX=0.75D0*(3.D0/PI)**THRD
      CX2=4.D0/9.D0*CX
      CK=0.3D0*(3.D0*PI*PI)**THRD2
      CK2=10.D0/9.D0*CK
C
      DO 10 I=1,NNN
10    RHOM(I)=0.D0
      N3=3*NUCS
      DO 20 I=1,N3
      DO 20 J=1,N3
20    HESS(I,J)=0.D0
C
      DO 100 IATM=1,NUCS
C
      NDATA=NSPL(IATM)
      NDATA1=NDATA+1
      H=1.D0/NDATA1
      CHRG=DFLOAT(NUCZ(IATM))
      RMIDZ=1.D0/CHRG**THRD
      DO 110 I=1,NDATA
      F1(I)=RARHO(I,IATM)
110   F2(I)=RAVEL(I,IATM)
      ZRMD=2.D0*CHRG*RMIDZ
      CALL SPLINE(H,F1,A1,B1,C1,NDATA,0,ZRMD,0.D0)
      CALL SPLINE(H,F2,A2,B2,C2,NDATA,2,0.D0,CHRG)
      if (ERRFLG) return
C
      DO 120 KKK=1,NNN
      X=XMESH(KKK)-XNUC(IATM)
      Y=YMESH(KKK)-YNUC(IATM)
      Z=ZMESH(KKK)-ZNUC(IATM)
      R=DSQRT(X**2+Y**2+Z**2)
      Q=R/(R+RMIDZ)
      DRQ=RMIDZ/(1.D0-Q)**2
      INTQ=IDINT(NDATA1*Q)
      DQ=Q-INTQ*H
      ARHO=(F1(INTQ)+DQ*(A1(INTQ)+DQ*(B1(INTQ)+DQ*C1(INTQ))))/R
      AVEL=(F2(INTQ)+DQ*(A2(INTQ)+DQ*(B2(INTQ)+DQ*C2(INTQ))))/R
      DARHO=(A1(INTQ)+2.D0*B1(INTQ)*DQ+3.D0*C1(INTQ)*DQ*DQ)/DRQ
      DAVEL=(A2(INTQ)+2.D0*B2(INTQ)*DQ+3.D0*C2(INTQ)*DQ*DQ)/DRQ
      DARHO1(KKK)=(DARHO-ARHO)/R
      DAVEL1(KKK)=(DAVEL-AVEL)/R
120   RHOM(KKK)=RHOM(KKK)+DABS(ARHO)
      WRITE(IBAS0,REC=IATM,ERR=9001)(DARHO1(I),I=1,NNN)
      WRITE(IBAS2,REC=IATM,ERR=9002)(DAVEL1(I),I=1,NNN)
C
100   CONTINUE
C
      DO 200 INUC=2,NUCS
      READ(IBAS0,REC=INUC)(DARHO1(I),I=1,NNN)
      READ(IBAS2,REC=INUC)(DAVEL1(I),I=1,NNN)
      DO 200 JNUC=1,INUC-1
      READ(IBAS0,REC=JNUC)(DARHO2(I),I=1,NNN)
      READ(IBAS2,REC=JNUC)(DAVEL2(I),I=1,NNN)
C
      Z1=DFLOAT(NUCZ(INUC))
      Z2=DFLOAT(NUCZ(JNUC))
      XYZNN(1)=XNUC(JNUC)-XNUC(INUC)
      XYZNN(2)=YNUC(JNUC)-YNUC(INUC)
      XYZNN(3)=ZNUC(JNUC)-ZNUC(INUC)
      RNN=DSQRT(XYZNN(1)**2+XYZNN(2)**2+XYZNN(3)**2)
C
      DO 200 IXYZ=1,3
      DO 200 JXYZ=1,3
      II=3*(INUC-1)+IXYZ
      JJ=3*(JNUC-1)+JXYZ
      HESSNN=-3.D0*XYZNN(IXYZ)*XYZNN(JXYZ)
      IF(IXYZ.EQ.JXYZ)HESSNN=HESSNN+RNN*RNN
      HESS(II,JJ) = Z1*Z2 * HESSNN / RNN**5
C
      DO 210 KKK=1,NNN
      XYZ1(1)=XMESH(KKK)-XNUC(INUC)
      XYZ1(2)=YMESH(KKK)-YNUC(INUC)
      XYZ1(3)=ZMESH(KKK)-ZNUC(INUC)
      XYZ2(1)=XMESH(KKK)-XNUC(JNUC)
      XYZ2(2)=YMESH(KKK)-YNUC(JNUC)
      XYZ2(3)=ZMESH(KKK)-ZNUC(JNUC)
      R1=DSQRT(XYZ1(1)**2+XYZ1(2)**2+XYZ1(3)**2)
      R2=DSQRT(XYZ2(1)**2+XYZ2(2)**2+XYZ2(3)**2)
      DVNUC1=-Z1*XYZ1(IXYZ)/R1**3
      DVNUC2=-Z2*XYZ2(JXYZ)/R2**3
      DRHO1=-DARHO1(KKK)*XYZ1(IXYZ)/R1
      DRHO2=-DARHO2(KKK)*XYZ2(JXYZ)/R2
      DVEL1=-DAVEL1(KKK)*XYZ1(IXYZ)/R1
      DVEL2=-DAVEL2(KKK)*XYZ2(JXYZ)/R2
210   HESS(II,JJ)=HESS(II,JJ)+WINTS(KKK)*
     +       (DRHO1*DRHO2*(CK2/RHOM(KKK)**THRD-CX2/RHOM(KKK)**THRD2)
     +       +DRHO1*DVNUC2+DRHO2*DVNUC1+0.5D0*(DRHO1*DVEL2+DRHO2*DVEL1))
C
200   CONTINUE
C
C     FILL UPPER TRIANGLE:
      DO 300 I=2,N3
      DO 300 J=1,I-1
300   HESS(J,I)=HESS(I,J)
C
C     FILL DIAGONAL BLOCKS
C     USING TRANSLATIONAL INVARIANCE:
      DO 400 IBLK=1,NUCS
      I0=3*(IBLK-1)
      DO 410 JBLK=1,NUCS
      IF(JBLK.EQ.IBLK)GO TO 410
      J0=3*(JBLK-1)
      HESS(I0+1,I0+1)=HESS(I0+1,I0+1)-HESS(I0+1,J0+1)
      HESS(I0+1,I0+2)=HESS(I0+1,I0+2)-HESS(I0+1,J0+2)
      HESS(I0+1,I0+3)=HESS(I0+1,I0+3)-HESS(I0+1,J0+3)
      HESS(I0+2,I0+1)=HESS(I0+2,I0+1)-HESS(I0+2,J0+1)
      HESS(I0+2,I0+2)=HESS(I0+2,I0+2)-HESS(I0+2,J0+2)
      HESS(I0+2,I0+3)=HESS(I0+2,I0+3)-HESS(I0+2,J0+3)
      HESS(I0+3,I0+1)=HESS(I0+3,I0+1)-HESS(I0+3,J0+1)
      HESS(I0+3,I0+2)=HESS(I0+3,I0+2)-HESS(I0+3,J0+2)
      HESS(I0+3,I0+3)=HESS(I0+3,I0+3)-HESS(I0+3,J0+3)
410   CONTINUE
400   CONTINUE
C
      RETURN
 9001 WRITE(IOUT,*) '*** Error writing IBAS0 in HSSIAN'
      GOTO 9999
 9002 WRITE(IOUT,*) '*** Error writing IBAS2 in HSSIAN'
 9999 ERRFLG = .TRUE.
      RETURN
      END
C
C
C
      SUBROUTINE NEWGEO(ITGEOM,DONE)
C
C     CALCULATION OF REVISED HESSIAN AND MOLECULAR GEOMETRY.
C
C    ********************************************************
C     NOTE:
C     ASSUMES THAT ATOM AT CARTESIAN ORIGIN IS FIXED.
C    ********************************************************
C
C     Modifies nuclear coordinates in NuclearFramework
C
C     TODO Review the variables in NEWGEO to see if any others
C     TODO   need to be saved from one invocation to the next.
C     TODO Refactor NEWGEO into separate routines dealing with
C     TODO   minimization and with saddle-point optimization.
C
      USE Constants
      USE ERRORHANDLING
      USE GeometryOptimization
      USE IOUNITS
      USE HessianMatrix
      USE KEYWORDS
      USE MPI, only : myRank
      USE NuclearFramework
      IMPLICIT NONE
      ! Arguments
      INTEGER,INTENT(IN) :: ITGEOM
      LOGICAL,INTENT(OUT):: DONE
      ! Locals
      INTEGER VAR(NUCS*3),I,J,N,INUC,IXYZ,MODE,MDPOS,IERR,K
      REAL(KIND(1D0)) COMPS(3),VECNEG(NUCS*3)
      REAL(KIND(1D0)) WORK(NUCS*9),WK1(NUCS*3),WK2(NUCS*3)
      REAL(KIND(1D0)) FVEC(NUCS*3),F0(NUCS*3),DELFRC(NUCS*3)
      REAL(KIND(1D0)) X,Y,Z,SHIFT,AVGPOS,AVGNEG,DOT,DOTMAX,RMS
      REAL(KIND(1D0)),PARAMETER :: FTOL = 1.D-03
C
      if (.not.allocated(HESS)) then
         ! first time through
         allocate(HESS(NUCS*3,NUCS*3))
         allocate(HEIGS(NUCS*3))
         allocate(HVECS(NUCS*3,NUCS*3))
      end if
C
      DONE=.TRUE.
C
      I=0
      N=0
      DO 100 INUC=1,NUCS
      DO 100 IXYZ=1,3
      I=I+1
      FVEC(I)=FORCE(INUC,IXYZ)
      IF(MOVE(INUC,IXYZ).EQ.1)THEN
      IF(DABS(FVEC(I)).GT.FTOL)DONE=.FALSE.
      N=N+1
      VAR(N)=I
      END IF
100   CONTINUE
C
      IF (myRank.EQ.0) THEN
      IF(DONE.AND.IOPTG.NE.-1)THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'=================================================='
      WRITE(IOUT,81)ITGEOM
  81  FORMAT(1X,'GEOMETRY CONVERGED AFTER',I3,' ITERATIONS.')
      WRITE(IOUT,*)'=================================================='
      RETURN
      END IF ! done and ioptg<>-1
C
      IF(ITGEOM.EQ.MGEOMS)THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'=================================================='
      WRITE(IOUT,82)ITGEOM
  82  FORMAT(1X,'GEOMETRY FAILS TO CONVERGE AFTER',I3,' ITERATIONS.')
      WRITE(IOUT,*)'=================================================='
      ERRFLG=.TRUE.
      RETURN
      END IF ! itgeom==mgeoms
      END IF ! myrank==0
C
      IF(ITGEOM.EQ.1)THEN
        CALL HSSIAN
        DO 210 I=1,N
        F0(I)=FVEC(VAR(I))
        DO 210 J=1,N
210     HVECS(I,J)=HESS(VAR(I),VAR(J))
        DO 211 I=1,N
        DO 211 J=1,N
211     HESS(I,J)=HVECS(I,J)
      ELSE
        DO 220 I=1,N
        F0(I)=FVEC(VAR(I))
220     DELFRC(I)=F0(I)-F1(I)
        IF(IOPTG.EQ.1) CALL UPDAT2(N,DELX,DELFRC)
        IF(IOPTG.EQ.-2)CALL UPDAT1(N,DELX,DELFRC)
        DO 221 I=1,N
        DO 221 J=1,N
221     HVECS(I,J)=HESS(I,J)
      END IF
C
      CALL DSYEV('V','U',N,HVECS,NUCS*3,HEIGS,WORK,NUCS*9,IERR)
      IF(IERR.NE.0)THEN
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** NEWGEO: ERROR IN ROUTINE DSYEV ***'
      ERRFLG = .TRUE.
      RETURN
      END IF
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'EIGENVALUES OF REDUCED HESSIAN:'
      DO 300 I=1,N
300   WRITE(IOUT,*)I,HEIGS(I)
      END IF ! myrank==0
C
      IF(IOPTG.EQ.-1)THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'EIGENVECTORS:'
C
      DO 310 MODE=1,N
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'MODE =',MODE
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'INUC   SYMB      X       Y       Z          EIGENVEC
     + COMPONENTS'
      END IF ! myrank==0
      I=0
      DO 320 INUC=1,NUCS
      DO 330 IXYZ=1,3
      COMPS(IXYZ)=0.D0
      IF(MOVE(INUC,IXYZ).EQ.1)THEN
      I=I+1
      COMPS(IXYZ)=HVECS(I,MODE)
      END IF ! MOVE()==1
330   CONTINUE
      X=XNUC(INUC)/DistanceScale
      Y=YNUC(INUC)/DistanceScale
      Z=ZNUC(INUC)/DistanceScale
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,'(1X,I4,4X,A2,3X,3F8.3,3X,3F8.4)')
     +       INUC,ZSYMB(INUC),X,Y,Z,(COMPS(K),K=1,3)
      END IF ! myrank==0
320   CONTINUE
310   CONTINUE
      RETURN
      END IF ! IOPTG==-1
C
      IF(ITGEOM.EQ.1.AND.IOPTG.EQ.1)THEN
      SHIFT=0.D0
      IF(HEIGS(1).LT.HMIN)SHIFT=HMIN-HEIGS(1)
      DO 400 I=1,N
      HEIGS(I)=HEIGS(I)+SHIFT
400   HESS(I,I)=HESS(I,I)+SHIFT
      AVGPOS=HEIGS(1)
      END IF ! ITGEOM==1 and IOPTG==1
C
      IF(ITGEOM.EQ.1.AND.IOPTG.EQ.-2)THEN
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'*** INITIAL TRANSITION VECTOR =',MDNEG
      END IF ! myrank==0
      MDPOS=2
      IF(MDNEG.NE.1)MDPOS=1
      SHIFT=0.D0
      IF(HEIGS(MDPOS).LT.HMIN)SHIFT=HMIN-HEIGS(MDPOS)
      DO 410 I=1,N
      HEIGS(I)=HEIGS(I)+SHIFT
410   HESS(I,I)=HESS(I,I)+SHIFT
      DO 420 I=1,N
      VECNEG(I)=HVECS(I,MDNEG)
      DO 420 J=1,N
420   HESS(I,J)=HESS(I,J)
     +         +(HTRANS-HEIGS(MDNEG))*HVECS(I,MDNEG)*HVECS(J,MDNEG)
      HEIGS(MDNEG)=HTRANS
      AVGNEG = HEIGS(MDNEG)
      AVGPOS = HEIGS(MDPOS)
      END IF ! itgeom==1 and ioptg==-2
C
      IF(ITGEOM.EQ.1)GO TO 1000
C
      IF(IOPTG.EQ.1)THEN
      IF(HEIGS(1).GT.0.D0)THEN
      AVGPOS=((ITGEOM-1)*AVGPOS+HEIGS(1))/DFLOAT(ITGEOM)
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'LOWEST MODE EIGENVALUE AVERAGE =',AVGPOS
      END IF ! myrank==0
      SHIFT=0.D0
      IF(HEIGS(1).LT.AVGPOS)SHIFT=AVGPOS-HEIGS(1)
      DO 500 I=1,N
500   HEIGS(I)=HEIGS(I)+SHIFT
      ELSE   ! HEIGS(1) <= 0
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'=================================================='
      WRITE(IOUT,*)'NON-POSITIVE-DEFINITE HESSIAN: SEARCH TERMINATED. '
      WRITE(IOUT,*)'=================================================='
      ERRFLG = .TRUE.
      RETURN
      END IF ! HEIGS...
      END IF ! IOPTG==1
C
      IF(IOPTG.EQ.-2)THEN
      MDNEG=0
      DOTMAX=0.D0
      DO 600 MODE=1,N
      DOT=0.D0
      DO 610 I=1,N
610   DOT=DOT+HVECS(I,MODE)*VECNEG(I)
      IF(DABS(DOT).GT.DOTMAX)THEN
      MDNEG=MODE
      DOTMAX=DABS(DOT)
      END IF  ! |DOT| > DOTMAX
600   CONTINUE
      IF(MDNEG.EQ.1.AND.HEIGS(1).LT.0.D0.AND.HEIGS(2).GT.0.D0)THEN
      DO 620 I=1,N
620   VECNEG(I)=HVECS(I,1)
      AVGNEG=((ITGEOM-1)*AVGNEG+HEIGS(1))/DFLOAT(ITGEOM)
      AVGPOS=((ITGEOM-1)*AVGPOS+HEIGS(2))/DFLOAT(ITGEOM)
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'TRANSITION VECTOR EIGENVALUE AVERAGE =',AVGNEG
      WRITE(IOUT,*)'LOWEST POSITIVE MODE EIGENVALUE AVERAGE =',AVGPOS
      END IF  ! myRank==1
      IF(HEIGS(1).GT.AVGNEG)HEIGS(1)=AVGNEG
      SHIFT=0.D0
      IF(HEIGS(2).LT.AVGPOS)SHIFT=AVGPOS-HEIGS(2)
      DO 630 I=2,N
630   HEIGS(I)=HEIGS(I)+SHIFT
      ELSE  ! MDNEG==1 etc
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'=================================================='
      WRITE(IOUT,*)'ILL-BEHAVED HESSIAN: SADDLE SEARCH TERMINATED.    '
      WRITE(IOUT,*)'=================================================='
      ERRFLG = .TRUE.
      RETURN
      END IF ! MDNEG==1 etc
      END IF ! IOPTG==-2
C
1000  CONTINUE
C
      DO 700 I=1,N
      WK1(I)=0.D0
      DO 710 J=1,N
C     EIGEN-COMPONENTS OF FORCE
710   WK1(I)=WK1(I)+HVECS(J,I)*F0(J)
C     EIGEN-COMPONENTS OF DELTA-X
700   WK2(I)=WK1(I)/HEIGS(I)
C
      RMS=0.D0
      DO 720 I=1,N
      DELX(I)=0.D0
      DO 730 J=1,N
730   DELX(I)=DELX(I)+HVECS(I,J)*WK2(J)
      F1(I)=F0(I)
720   RMS=RMS+DELX(I)**2
      RMS=DSQRT(RMS/N)
C
C     TODO find a more elegant way to code this
      I=0
      DO 800 INUC=1,NUCS
      IF(MOVE(INUC,1).EQ.1)THEN
        I=I+1
        XNUC(INUC)=XNUC(INUC)+DELX(I)
      END IF
      IF(MOVE(INUC,2).EQ.1)THEN
        I=I+1
        YNUC(INUC)=YNUC(INUC)+DELX(I)
      END IF
      IF(MOVE(INUC,3).EQ.1)THEN
        I=I+1
        ZNUC(INUC)=ZNUC(INUC)+DELX(I)
      END IF
800   CONTINUE
C
      IF (myRank.EQ.0) THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'=================================================='
      WRITE(IOUT,83)ITGEOM+1
  83  FORMAT(1X,'GEOMETRY OPTIMIZER:  GEOMETRY',I3)
      WRITE(IOUT,*)'=================================================='
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'INUC   SYMB      X        Y        Z'
      DO 810 INUC=1,NUCS
      X=XNUC(INUC)/DistanceScale
      Y=YNUC(INUC)/DistanceScale
      Z=ZNUC(INUC)/DistanceScale
810   WRITE(IOUT,'(1X,I4,4X,A2,3X,3F9.4)')
     +              INUC,ZSYMB(INUC),X,Y,Z
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'R.M.S. CHANGE =',RMS/DistanceScale
      END IF ! myrank==0
C
      RETURN
      END
C
C
C
      SUBROUTINE UPDAT1(N,DELX,DELFRC)
C
C     RANK ONE HESSIAN UPDATE (MURTAGH-SARGENT).
C
      USE HessianMatrix
      IMPLICIT NONE
      ! Arguments
      INTEGER,        INTENT(IN) :: N
      REAL(KIND(1D0)),INTENT(IN) :: DELX(N),DELFRC(N)
      ! Locals
      INTEGER I,J
      REAL(KIND(1D0)) U(N),CON,SUM
C
      CON=0.D0
      DO 100 I=1,N
      SUM=0.D0
      DO 110 J=1,N
110   SUM=SUM+HESS(I,J)*DELX(J)
      U(I)=DELFRC(I)+SUM
100   CON=CON+U(I)*DELX(I)
      CON=-1.D0/CON
C
      DO 200 I=1,N
      DO 200 J=1,N
200   HESS(I,J)=HESS(I,J)+CON*U(I)*U(J)
C
      RETURN
      END
C
C
C
      SUBROUTINE UPDAT2(N,DELX,DELFRC)
C
C     RANK TWO HESSIAN UPDATE (BROYDEN-FLETCHER-GOLDFARB-SHANNO)
C
      USE HessianMatrix
      IMPLICIT NONE
      ! Arguments
      INTEGER,        INTENT(IN) :: N
      REAL(KIND(1D0)),INTENT(IN) :: DELX(N),DELFRC(N)
      ! Locals
      INTEGER I,J
      REAL(KIND(1D0)) U(N),A,B
C
      A=0.D0
      B=0.D0
      DO 100 I=1,N
      U(I)=0.D0
      DO 110 J=1,N
110   U(I)=U(I)+HESS(I,J)*DELX(J)
      A=A+U(I)*DELX(I)
100   B=B+DELFRC(I)*DELX(I)
      A=-1.D0/A
      B=-1.D0/B
C
      DO 200 I=1,N
      DO 200 J=1,N
200   HESS(I,J)=HESS(I,J)+A*U(I)*U(J)+B*DELFRC(I)*DELFRC(J)
C
      RETURN
      END
C
C
C
      SUBROUTINE VIBRAN(ZPE)
C
C     NORMAL-MODE VIBRATIONAL ANALYSIS.
C     CURRENTLY USED ONLY AS TEST OF MODEL HESSIAN.
C
      USE ERRORHANDLING
      USE HessianMatrix
      USE IOUNITS
      USE NuclearFramework
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WORK(NUCS*3)
      REAL*8 MASS(100)
      DATA (MASS(I),I=1,100)/
     +  1.008,4.003,6.941,9.012,10.81,12.01,14.01,16.00,19.00,20.18,
     +  22.99,24.31,26.98,28.09,30.97,32.06,35.45,39.95,39.10,40.08,
     +  44.96,47.88,50.94,52.00,54.94,55.85,58.93,58.70,63.55,65.38,
     +  69.72,72.59,74.92,78.96,79.90,83.80,85.47,87.62,88.91,91.22,
     +  92.91,95.94,98.00,101.1,102.9,106.4,107.9,112.4,114.8,118.7,
     +  121.8,127.6,126.9,131.3,132.9,137.3,138.9,140.1,140.9,144.2,
     +  145.0,150.4,152.0,157.3,158.9,162.5,164.9,167.3,168.9,173.0,
     +  175.0,178.5,180.9,183.9,186.2,190.2,192.2,195.1,197.0,200.6,
     +  204.4,207.2,209.0,209.0,210.0,222.0,223.0,226.0,227.0,232.0,
     +  231.0,238.0,237.0,244.0,243.0,247.0,247.0,251.0,252.0,257.0/
C
C     S.I. UNITS!
C
C     ATOMIC MASSES TO KG:
      DO 10 I=1,100
10    MASS(I)=1.6605D0*MASS(I)
C
C     S.I. CONVERSION FACTOR FOR HESSIAN:
      CON=4.3597D0/5.2918D0**2
C
      CALL HSSIAN
C
      I=0
      DO 100 INUC=1,NUCS
      DO 100 IXYZ=1,3
      I=I+1
      J=0
      DO 200 JNUC=1,NUCS
      DO 200 JXYZ=1,3
      J=J+1
200   HESS(I,J)=CON*HESS(I,J)
     +         /DSQRT(MASS(NUCZ(INUC))*MASS(NUCZ(JNUC)))
100   CONTINUE
C
      N3=3*NUCS
      HVECS(:,:) = HESS(:,:)
      CALL DSYEV('V','U',N3,HVECS,N3,HEIGS,WORK,NUCS*9,IERR)
      IF(IERR.NE.0)THEN
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** VIBRAN: ERROR IN ROUTINE DSYEV ***'
      ERRFLG = .TRUE.
      RETURN
      END IF
C     WRITE(IOUT,*)' '
C     WRITE(IOUT,*)' '
C     WRITE(IOUT,*)'VIBRATIONAL FREQUENCIES (CM-1):'
C
      ZPE=0.D0
      IBEG=4
      IF(NUCS.EQ.2)IBEG=3
      DO 300 I=IBEG,N3
      HEIGS(I)=DABS(HEIGS(I))
      HFREQ=1.0546D-34*DSQRT(1.D31*HEIGS(I))
C    "HFREQ" IS IN JOULES!
      ZPE=ZPE+HFREQ
C     WRITE(IOUT,*)I,5.0341D22*HFREQ
300   CONTINUE
C
      ZPE=0.5D0*2.2937D17*ZPE
      WRITE(IOUT,*)' '
      WRITE(IOUT,81)ZPE
  81  FORMAT(1X,'VIBRATIONAL ZPE (IN ATOMIC UNITS)  =  ',F11.4)
C
      RETURN
      END
C
C
C
C$PRAGMA SUN OPT=5
      SUBROUTINE EXACTX(EXX,EXDENS)
C
C     EXACT EXCHANGE ENERGY AND EXCHANGE ENERGY DENSITY
C     TODO deallocate V(:,:) before coming in here?  Cf. XCMODS
C
      USE CONSTANTS
      USE CORES
      USE ERRORHANDLING
      USE GRID
      USE IOUNITS
      USE KEYWORDS
      USE MPI
      USE OccupationNumbers
      USE SIZES
      IMPLICIT NONE
      ! Arguments
      REAL(KIND(1D0)),INTENT(OUT):: EXX,EXDENS(NNN,2)
      ! Locals
      REAL(KIND(1D0)),PARAMETER :: SML = 1.D-03
      REAL(KIND(1D0)) PSII(NNN),PSIJ(NNN),PROD(NNN),COUL(NNN)
      INTEGER KKK,IMORB,JMORB,I,ISPIN,IREC,JREC,IO
      REAL(KIND(1D0)) SUM,FACTOR,EXXTOT
C
      EXX=0.D0
      DO 10 KKK=1,NNN
      EXDENS(KKK,1)=0.D0
10    EXDENS(KKK,2)=0.D0
C
C     CORE-CORE EXCHANGE:
      IF(NCORBS.NE.0.and.myRank.eq.0)THEN
      REWIND(ICOR0)
      DO 100 IMORB=1,NCORBS
      READ(ICOR0)(PSII(I),I=1,NNN)
      REWIND(ICOR0)
      DO 100 JMORB=1,IMORB
      READ(ICOR0)(PSIJ(I),I=1,NNN)
      IF(ICNTR(JMORB).NE.ICNTR(IMORB)) GO TO 100
      DO 110 KKK=1,NNN
110   PROD(KKK)=PSII(KKK)*PSIJ(KKK)
      CALL POISS1(PROD,COUL,ICNTR(IMORB))
      if (ERRFLG) return
      SUM=0.D0
      FACTOR=1.D0
      IF(JMORB.EQ.IMORB)FACTOR=0.5D0
      DO 120 KKK=1,NNN
      SUM=SUM+WINTS(KKK)*PROD(KKK)*COUL(KKK)
      EXDENS(KKK,1)=EXDENS(KKK,1)-FACTOR*PROD(KKK)*COUL(KKK)
120   EXDENS(KKK,2)=EXDENS(KKK,2)-FACTOR*PROD(KKK)*COUL(KKK)
      EXX=EXX-FACTOR*SUM
100   CONTINUE ! IMORB,JMORB
      IF(NSPINS.EQ.2)EXX=2.D0*EXX
      END IF ! NCORBS>0
C
C     VALENCE-VALENCE AND CORE-VALENCE EXCHANGE:
      DO 1000 ISPIN=1,NSPINS
      IF(ISPIN.EQ.1)IO=MOSA0
      IF(ISPIN.EQ.2)IO=MOSB0
      DO 200 IMORB=1,NMORBS
      IF(OCCS(IMORB,ISPIN).LT.SML)GO TO 200
      IF(myRank.EQ.IOWNER(IMORB))THEN
        IREC=MOSEQ(IMORB)
        READ(IO,REC=IREC)(PSII(I),I=1,NNN)
      END IF  
      CALL MPI_BCAST(PSII,NNN,MPI_DOUBLE_PRECISION,IOWNER(IMORB),
     +               MPI_COMM_WORLD,iError)
      IF(NCORBS.NE.0.and.myRank.EQ.IOWNER(IMORB))THEN 
C       This depends on every node having their own copy of ICOR0
        REWIND(ICOR0)
        DO 300 JMORB=1,NCORBS
        READ(ICOR0)(PSIJ(I),I=1,NNN)
        DO 310 KKK=1,NNN
310     PROD(KKK)=PSII(KKK)*PSIJ(KKK)
        CALL POISS1(PROD,COUL,ICNTR(JMORB))
        if (ERRFLG) return
        SUM=0.D0
        FACTOR=OCCS(IMORB,ISPIN)
        DO 320 KKK=1,NNN
        SUM=SUM+WINTS(KKK)*PROD(KKK)*COUL(KKK)
320     EXDENS(KKK,ISPIN)=EXDENS(KKK,ISPIN)-FACTOR*PROD(KKK)*COUL(KKK)
        EXX=EXX-FACTOR*SUM
300     CONTINUE ! JMORB
      END IF ! NCORBS>0
      DO 400 JMORB=1,IMORB
      IF(OCCS(JMORB,ISPIN).LT.SML)GO TO 400
      IF(myRank.EQ.IOWNER(JMORB))THEN
        JREC=MOSEQ(JMORB)
        READ(IO,REC=JREC)(PSIJ(I),I=1,NNN)
        DO 410 KKK=1,NNN
410     PROD(KKK)=PSII(KKK)*PSIJ(KKK)
        CALL POISS(PROD,COUL,0)
        if (ERRFLG) return
        SUM=0.D0
        FACTOR=OCCS(IMORB,ISPIN)*OCCS(JMORB,ISPIN)
        IF(JMORB.EQ.IMORB)FACTOR=0.5D0*FACTOR
        DO 420 KKK=1,NNN
        SUM=SUM+WINTS(KKK)*PROD(KKK)*COUL(KKK)
420     EXDENS(KKK,ISPIN)=EXDENS(KKK,ISPIN)-FACTOR*PROD(KKK)*COUL(KKK)
        EXX=EXX-FACTOR*SUM
      END IF ! myRank==IOWNER
400   CONTINUE ! JMORB
200   CONTINUE ! IMORB
C
C     Sum up EXDENS(:,ISPIN) using PROD as a temporary
      CALL MPI_REDUCE(EXDENS(1,ISPIN),PROD,NNN,MPI_DOUBLE_PRECISION,
     +                MPI_SUM,0,MPI_COMM_WORLD,iError)
      IF(myRank.EQ.0)THEN
        DO 20 KKK=1,NNN
20        EXDENS(KKK,ISPIN)=MIN(-EPS,PROD(KKK))
      END IF
C
1000  CONTINUE ! ISPIN
C
C     Sum up EXX
      CALL MPI_REDUCE(EXX,EXXTOT,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                0,MPI_COMM_WORLD,iError)
      EXX=EXXTOT 
C ... but no need to globalize it because root does all the remaining work.
c     CALL MPI_BCAST(EXX,1,MPI_DOUBLE_PRECISION,0,
c    +               MPI_COMM_WORLD,iError)
      IF(NSPINS.EQ.1)THEN
      EXX=2.D0*EXX
      DO 500 KKK=1,NNN
500   EXDENS(KKK,2)=EXDENS(KKK,1)
      END IF
C
      RETURN
      END
C
C
C
C$PRAGMA SUN OPT=5
      SUBROUTINE POISS1(RHO,VEL,ICNTR)
C
C     SOLUTION OF POISSON'S EQUATION FOR THE COULOMB POTENTIAL
C     OF A SINGLE-CENTRE SOURCE DENSITY.
C
      USE CONSTANTS
      USE ERRORHANDLING
      USE IOUNITS
      USE NuclearFramework
      USE NuclearMeshes
      USE SecondDerivativeMatrices
      USE SIZES
      USE SphericalHarmonics
      IMPLICIT NONE
      REAL(KIND(1D0)),INTENT(IN) :: RHO(NNN)
      REAL(KIND(1D0)),INTENT(OUT):: VEL(NNN)
      INTEGER,        INTENT(IN) :: ICNTR
      REAL(KIND(1D0)) ULM(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) A(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) B(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) C(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) BND(10,-2:MR+3)
      REAL(KIND(1D0)) RHOLM(MR),XX(MR)
      REAL(KIND(1D0)) RS(ML),XUNIT(ML),YUNIT(ML),ZUNIT(ML),DQS(ML)
      REAL(KIND(1D0)) H,CHARGE,SUM,R,BC,X,Y,Z,Q,DQ,RFAC
      INTEGER IPVT(MR),INTQS(ML),I,ICEL,NRADS,NANGS,NWAVES,KKK,IY,JY
      INTEGER IR,IL,L,L2,K,IERR,INUC,INTQ
C
      DO 10 I=1,NNN
10    VEL(I) = 0.D0
C
      ICEL=ICNTR
C
      NRADS=NR(ICEL)
      NANGS=NL(ICEL)
      NWAVES=LMAX(ICEL)
      H=1.D0/(NRADS+1)
C
      KKK=KKKTAB(ICEL)
      CHARGE=0.D0
      DO 101 IR=1,NRADS
      DO 101 IL=1,NANGS
      KKK=KKK+1
101   CHARGE=CHARGE+WRADS(IR,ICEL)*WANGS(IL,ICEL)*RHO(KKK)
C
C     SPHERICAL EXPANSION OF THE DENSITY:
C
      CALL YCALC(XANGS(:,ICEL),YANGS(:,ICEL),ZANGS(:,ICEL),NANGS,NWAVES)
C
      DO 200 IY=0,NWAVES
      DO 200 JY=0,NWAVES
C
      KKK=KKKTAB(ICEL)
      DO 210 IR=1,NRADS
      SUM=0.D0
      DO 211 IL=1,NANGS
211   SUM=SUM+WANGS(IL,ICEL)*YLM(IL,IY,JY)*RHO(KKK+IL)
      RHOLM(IR)=SUM
210   KKK=KKK+NANGS
C
C     SOLVE RADIAL LM EQUATION:
C
      L=MAX0(IY,JY)
      L2=L*(L+1)
C
C     L.H.S. (LAPACK BAND MATRIX):
C
      DO 220 K=1,7
      DO 220 IR=1,NRADS
220   BND(11-K,IR+K-4)=D2R(IR,K,ICEL)
      DO 221 IR=1,NRADS
      XX(IR)=-FOURPI*RADS(IR,ICEL)*RHOLM(IR)
221   BND(7,IR)=BND(7,IR)-L2/RADS(IR,ICEL)**2
C
C     INFINITE-R BOUNDARY CONDITION:
C
      IF(L.EQ.0)THEN
      XX(NRADS)=XX(NRADS)-CHARGE*D2R(NRADS,5,ICEL)
      XX(NRADS-1)=XX(NRADS-1)-CHARGE*D2R(NRADS-1,6,ICEL)
      XX(NRADS-2)=XX(NRADS-2)-CHARGE*D2R(NRADS-2,7,ICEL)
      END IF
C
      CALL DGBSV(NRADS,3,3,1,BND(1:10,1:),10,IPVT,XX,MR,IERR)
      IF(IERR.NE.0)THEN
      WRITE(IOUT,'(/''*** POISS1: ERROR IN ROUTINE DGBSV ***'')')
      ERRFLG = .TRUE.
      RETURN
      END IF
C
      ULM(0,IY,JY)=0.D0
      DO 230 IR=1,NRADS
230   ULM(IR,IY,JY)=XX(IR)
C
      KKK=KKKTAB(ICEL)
      DO 240 IR=1,NRADS
      R=RADS(IR,ICEL)
      DO 241 IL=1,NANGS
241   VEL(KKK+IL)=VEL(KKK+IL)+ULM(IR,IY,JY)*YLM(IL,IY,JY)/R
240   KKK=KKK+NANGS
C
      BC=0.D0
      IF(L.EQ.0)THEN
      CALL SPLINE(H,ULM(0,IY,JY),A(0,IY,JY),B(0,IY,JY),C(0,IY,JY),
     +                                        NRADS,2,0.D0,CHARGE)
      ELSE
      CALL SPLINE(H,ULM(0,IY,JY),A(0,IY,JY),B(0,IY,JY),C(0,IY,JY),
     +                                          NRADS,1,0.D0,0.D0)
      END IF
200   CONTINUE
C
      DO 300 INUC=1,NUCS
      IF(INUC.EQ.ICEL)GO TO 300
C
      KKK=KKKTAB(INUC)
      DO 310 IR=1,NR(INUC)
      DO 320 IL=1,NL(INUC)
C     CALC COORDS ON UNIT SPHERE ABOUT ICEL:
      X=XNUC(INUC)+XANGS(IL,INUC)*RADS(IR,INUC)-XNUC(ICEL)
      Y=YNUC(INUC)+YANGS(IL,INUC)*RADS(IR,INUC)-YNUC(ICEL)
      Z=ZNUC(INUC)+ZANGS(IL,INUC)*RADS(IR,INUC)-ZNUC(ICEL)
      RS(IL)=DSQRT(X**2+Y**2+Z**2)
      XUNIT(IL)=X/RS(IL)
      YUNIT(IL)=Y/RS(IL)
      ZUNIT(IL)=Z/RS(IL)
      Q=RS(IL)/(RS(IL)+RMID(ICEL))
      INTQS(IL)=IDINT((NRADS+1)*Q)
      DQS(IL)=Q-INTQS(IL)*H
320   CONTINUE
      CALL YCALC(XUNIT,YUNIT,ZUNIT,NL(INUC),NWAVES)
      DO 330 IY=0,NWAVES
      DO 330 JY=0,NWAVES
      DO 330 IL=1,NL(INUC)
      INTQ=INTQS(IL)
      DQ=DQS(IL)
      Y=ULM(INTQ,IY,JY)
      RFAC=Y+DQ*(A(INTQ,IY,JY)+DQ*(B(INTQ,IY,JY)+DQ*C(INTQ,IY,JY)))
330   VEL(KKK+IL)=VEL(KKK+IL)+RFAC*YLM(IL,IY,JY)/RS(IL)
310   KKK=KKK+NL(INUC)
C
300   CONTINUE
C
      RETURN
      END
C
C
C
C$PRAGMA SUN OPT=5
      SUBROUTINE FDPSI(PSIDC,IORDER,PSIXYZ,ICORB)
C
C     FIRST AND SECOND PARTIAL DERIVATIVES (CARTESIAN)
C     OF "DECOMPOSED" PSI BY FINITE DIFFERENCES.
C
C     ON OUTPUT, PSIXYZ(I,***) CONTAINS
C     FOR I=0      :  THE RECONSTRUCTED PSI
C     FOR I=1,2,3  :  1ST DERIVS  X,Y,Z     (IF IORDER.GE.1)
C     FOR I=4,5,6  :  2ND DERIVS  XX,YY,ZZ  (IF IORDER.GE.2)
C     FOR I=7,8,9  :  2ND DERIVS  YZ,XZ,XY  (IF IORDER.EQ.3)
C
      USE CORES
      USE NuclearFramework
      USE NuclearMeshes
      USE SIZES
      USE SphericalHarmonics
      IMPLICIT NONE
C     Arguments
      REAL(KIND(1D0)),INTENT(IN) :: PSIDC(NNN)
      INTEGER,        INTENT(IN) :: IORDER
      REAL(KIND(1D0)),INTENT(OUT):: PSIXYZ(0:3,NNN)
      INTEGER,        INTENT(IN) :: ICORB
C     Locals
      REAL(KIND(1D0)) PSILM(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) A(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) B(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) C(0:MR,0:MYLM,0:MYLM)
      REAL(KIND(1D0)) RS(ML),XUNIT(ML),YUNIT(ML),ZUNIT(ML),DQS(ML)
      REAL(KIND(1D0)) H,SUM,ZRMD,X,Y,Z,Q,DQ,RFAC,PRODRL
      INTEGER INTQS(ML),NXYZ,NDEL,I,K,I0,I1,ICEL,NRADS,NANGS,NWAVES
      INTEGER IY,JY,L,KKK,IR,IL,IDEL,INUC,INTQ,IXYZ
      REAL(KIND(1D0)) DELX(0:12),DELY(0:12),DELZ(0:12),COEF(0:9,0:12)
      REAL(KIND(1D0)),PARAMETER :: FDH=1.D-04
      INTEGER,PARAMETER :: 
     &         IDELX(0:12) = (/ 0,+1,-1, 0, 0, 0, 0, 0, 0,+1,-1,+1,-1/),
     &         IDELY(0:12) = (/ 0, 0, 0,+1,-1, 0, 0,+1,-1, 0, 0,+1,-1/),
     &         IDELZ(0:12) = (/ 0, 0, 0, 0, 0,+1,-1,+1,-1,+1,-1, 0, 0/)
      INTEGER ICOEF(1:9,0:12)
      DATA (ICOEF(1,I),I=0,12)/ 0,+1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(2,I),I=0,12)/ 0, 0, 0,+1,-1, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(3,I),I=0,12)/ 0, 0, 0, 0, 0,+1,-1, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(4,I),I=0,12)/-2,+1,+1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(5,I),I=0,12)/-2, 0, 0,+1,+1, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(6,I),I=0,12)/-2, 0, 0, 0, 0,+1,+1, 0, 0, 0, 0, 0, 0/
      DATA (ICOEF(7,I),I=0,12)/-2, 0, 0,+1,+1,+1,+1,-1,-1, 0, 0, 0, 0/
      DATA (ICOEF(8,I),I=0,12)/-2,+1,+1, 0, 0,+1,+1, 0, 0,-1,-1, 0, 0/
      DATA (ICOEF(9,I),I=0,12)/-2,+1,+1,+1,+1, 0, 0, 0, 0, 0, 0,-1,-1/
      DATA (COEF(0,I),I=0,12)/1.D0,12*0.D0/
C
      IF(IORDER.EQ.1)THEN
      NXYZ=3
      NDEL=6
      ELSE IF(IORDER.EQ.2)THEN
      NXYZ=6
      NDEL=6
      ELSE IF(IORDER.EQ.3)THEN
      NXYZ=9
      NDEL=12
      END IF
      DO 10 I=0,NDEL
      DELX(I)=FDH*IDELX(I)
      DELY(I)=FDH*IDELY(I)
      DELZ(I)=FDH*IDELZ(I)
      COEF(1,I)=0.5D0*ICOEF(1,I)/FDH
      COEF(2,I)=0.5D0*ICOEF(2,I)/FDH
      COEF(3,I)=0.5D0*ICOEF(3,I)/FDH
      COEF(4,I)=ICOEF(4,I)/FDH**2
      COEF(5,I)=ICOEF(5,I)/FDH**2
      COEF(6,I)=ICOEF(6,I)/FDH**2
      COEF(7,I)=-0.5D0*ICOEF(7,I)/FDH**2
      COEF(8,I)=-0.5D0*ICOEF(8,I)/FDH**2
10    COEF(9,I)=-0.5D0*ICOEF(9,I)/FDH**2
      DO 20 K=1,NNN
      DO 20 I=0,NXYZ
20    PSIXYZ(I,K)=0.D0
C
      I0=1
      I1=NUCS
      IF(ICORB.NE.0)THEN
      I0=ICNTR(ICORB)
      I1=ICNTR(ICORB)
      END IF
      DO 100 ICEL=I0,I1
C
      NRADS=NR(ICEL)
      NANGS=NL(ICEL)
      NWAVES=NY(ICEL)
      IF(ICORB.NE.0)NWAVES=LCOR(ICORB)
      H=1.D0/(NRADS+1)
C
C     SPHERICAL EXPANSION OF PSI:
C
      CALL YCALC(XANGS(:,ICEL),YANGS(:,ICEL),ZANGS(:,ICEL),
     +               NANGS,NWAVES)
C
      DO 200 IY=0,NWAVES
      DO 200 JY=0,NWAVES
      L=MAX0(IY,JY)
C
      PSILM(0,IY,JY)=0.D0
      KKK=KKKTAB(ICEL)
      DO 210 IR=1,NRADS
      SUM=0.D0
      DO 211 IL=1,NANGS
211   SUM=SUM+WANGS(IL,ICEL)*YLM(IL,IY,JY)*PSIDC(KKK+IL)
      PSILM(IR,IY,JY)=RADS(IR,ICEL)*SUM
210   KKK=KKK+NANGS
C
      IF(L.EQ.0)THEN
      ZRMD=DFLOAT(NUCZ(ICEL))*RMID(ICEL)
      CALL SPLINE(H,PSILM(0,IY,JY),A(0,IY,JY),B(0,IY,JY),C(0,IY,JY),
     +                                            NRADS,0,ZRMD,0.D0)
      ELSE
      CALL SPLINE(H,PSILM(0,IY,JY),A(0,IY,JY),B(0,IY,JY),C(0,IY,JY),
     +                                            NRADS,1,0.D0,0.D0)
      END IF
200   CONTINUE
C
      DO 12 IDEL=0,NDEL
C
      KKK=0
      DO 300 INUC=1,NUCS
C
      DO 310 IR=1,NR(INUC)
      DO 320 IL=1,NL(INUC)
C     CALC COORDS ON UNIT SPHERE ABOUT ICEL:
      X=XNUC(INUC)+XANGS(IL,INUC)*RADS(IR,INUC)-XNUC(ICEL)+DELX(IDEL)
      Y=YNUC(INUC)+YANGS(IL,INUC)*RADS(IR,INUC)-YNUC(ICEL)+DELY(IDEL)
      Z=ZNUC(INUC)+ZANGS(IL,INUC)*RADS(IR,INUC)-ZNUC(ICEL)+DELZ(IDEL)
      RS(IL)=DSQRT(X**2+Y**2+Z**2)
      XUNIT(IL)=X/RS(IL)
      YUNIT(IL)=Y/RS(IL)
      ZUNIT(IL)=Z/RS(IL)
      Q=RS(IL)/(RS(IL)+RMID(ICEL))
      INTQS(IL)=IDINT((NRADS+1)*Q)
      DQS(IL)=Q-INTQS(IL)*H
320   CONTINUE
      CALL YCALC(XUNIT,YUNIT,ZUNIT,NL(INUC),NWAVES)
      DO 330 IY=0,NWAVES
      DO 330 JY=0,NWAVES
      DO 330 IL=1,NL(INUC)
      INTQ=INTQS(IL)
      DQ=DQS(IL)
      Y=PSILM(INTQ,IY,JY)
      RFAC=Y+DQ*(A(INTQ,IY,JY)+DQ*(B(INTQ,IY,JY)+DQ*C(INTQ,IY,JY)))
      PRODRL=RFAC*YLM(IL,IY,JY)/RS(IL)
      DO 330 IXYZ=0,NXYZ
      PSIXYZ(IXYZ,KKK+IL)=PSIXYZ(IXYZ,KKK+IL)+COEF(IXYZ,IDEL)*PRODRL
330   CONTINUE
      KKK=KKK+NL(INUC)
310   CONTINUE
C
300   CONTINUE
C
12    CONTINUE
C
100   CONTINUE
C
      RETURN
      END
