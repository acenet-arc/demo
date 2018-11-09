C $Header: /home/ross/cvsroot/parnum/scfnumx.f,v 1.6 2005/04/08 15:09:57 ross Exp $
      SUBROUTINE SCFNUMX(IT,EKIN,RHO,RHOS,TAUS)
C
C     NUMERICAL SCF ITERATION FOR EXACT EXCHANGE POTENTIAL
C
      USE CORES
      USE ERRORHANDLING
      USE GRID
      USE IOUNITS
      USE KEYWORDS
      USE LCAOMatrices
      USE MPI
      USE OccupationNumbers
      USE POTENTIALS, VX=>VXC
      USE SCFConvergence
      USE SIZES
      IMPLICIT NONE
      ! Arguments
      INTEGER,         INTENT(IN) :: IT
      REAL(KIND(1D0)), INTENT(OUT):: EKIN
      REAL(KIND(1D0)), INTENT(OUT):: RHO(NNN),RHOS(NNN,2),TAUS(NNN,2)
      ! Locals
      REAL(KIND(1D0)) ERRS(NBASIS),RATS(NBASIS),BUF(NBASIS)
      REAL(KIND(1D0)) PSI0(NNN),PSI2(NNN),PSIDC(NNN)
      REAL(KIND(1D0)) RES0(NNN),RES2(NNN),RESDC(NNN),STOR(NNN)
      REAL(KIND(1D0)) XPSI(NNN),XRES(NNN)
      REAL(KIND(1D0)) DELEIG,EIG,PROJ,CNORM,SUMK,FACTOR,PSISQ,TTERM
      INTEGER I,KKK,ISPIN,IO0,IO2,IODC,IOX,IMORB,IREC,JMORB,JREC,IPROC
C
      IF(IT.EQ.NITS0+1)THEN
      DO 1 I=1,NMORBS
      RATSUM(I,1)=0.D0
1     RATSUM(I,2)=0.D0
      END IF
C
      EKIN=0.D0
      IF(MYRANK.EQ.0)THEN
      DO 10 KKK=1,NNN
        RHOS(KKK,1)=0.5D0*RHOCOR(KKK)
        RHOS(KKK,2)=0.5D0*RHOCOR(KKK)
        TAUS(KKK,1)=0.5D0*TAUCOR(KKK)
        TAUS(KKK,2)=0.5D0*TAUCOR(KKK)
10    CONTINUE
      ELSE
      DO 20 KKK=1,NNN
        RHOS(KKK,1)=0.D0
        RHOS(KKK,2)=0.D0
        TAUS(KKK,1)=0.D0
        TAUS(KKK,2)=0.D0
20    CONTINUE
      END IF
C
      DO 1000 ISPIN=1,NSPINS
      IF(ISPIN.EQ.1)THEN
      IO0=MOSA0
      IO2=MOSA2
      IODC=MODCA
      IOX=41
      ELSE
      IO0=MOSB0
      IO2=MOSB2
      IODC=MODCB
      IOX=42
      END IF
C
      DO 100 IMORB=1,NMORBS
      if (myrank.eq.IOwner(IMORB)) then

      irec=MOSeq(IMORB)
      READ(IO0 ,REC=irec,ERR=9001)(PSI0(I), I=1,NNN)
      READ(IO2 ,REC=irec,ERR=9001)(PSI2(I), I=1,NNN)
      READ(IODC,REC=irec,ERR=9001)(PSIDC(I),I=1,NNN)
C
      CALL XOP(ISPIN,IOX,PSI0,XPSI)
      if (ERRFLG) return
C
      DELEIG=0.D0
      DO 110 KKK=1,NNN
      STOR(KKK)=-0.5D0*PSI2(KKK)+XPSI(KKK)
     +          +(V(KKK,ISPIN)-EIGS(IMORB,ISPIN))*PSI0(KKK)
110   DELEIG=DELEIG+WINTS(KKK)*PSI0(KKK)*STOR(KKK)
      EIG=EIGS(IMORB,ISPIN)+DELEIG
      DO 120 KKK=1,NNN
120   STOR(KKK)=DELEIG*PSI0(KKK)-STOR(KKK)
      CALL SCHROD(V(1,ISPIN),VX(1,ISPIN),
     +                     EIG,STOR,RES0,RES2,RESDC,ERRS(IMORB))
      if (ERRFLG) return
      CALL XOP(ISPIN,IOX,RES0,XRES)
      if (ERRFLG) return
      CALL MIXRES(V(1,ISPIN),EIG,PSI0,PSI2,PSIDC,RES0,RES2,RESDC,
     +            RATS(IMORB),XPSI,XRES)
      if (ERRFLG) return
      WRITE(IBAS0,REC=irec,ERR=9001)(PSI0(I),I=1,NNN)
      WRITE(IBAS2,REC=irec,ERR=9001)(PSI2(I),I=1,NNN)
      WRITE(IBSDC,REC=irec,ERR=9001)(PSIDC(I),I=1,NNN)
      WRITE(31,   REC=irec,ERR=9001)(XPSI(I),I=1,NNN)
C
      end if ! myrank==IOwner
100   CONTINUE
C
      DO 200 IMORB=1,NMORBS
      IF(MYRANK.EQ.IOWNER(IMORB))THEN
        IREC=MOSEQ(IMORB)
        READ(IBAS0,REC=IREC,ERR=9001)(PSI0(I),I=1,NNN)
        READ(IBAS2,REC=IREC,ERR=9001)(PSI2(I),I=1,NNN)
        READ(IBSDC,REC=IREC,ERR=9001)(PSIDC(I),I=1,NNN)
        READ(31,   REC=IREC,ERR=9001)(XPSI(I),I=1,NNN)
      END IF ! myrank==IOwner
C
C     ORTHOGONALIZE:
C     We'll keep the highest already-orthogonalized MO in PSI*,
C     and we'll project it out of all *higher* MOs, loading them
C     successively into RES*.
C
      IF(IORTH.EQ.1)THEN
      CALL MPI_BCAST(PSI0,NNN,MPI_DOUBLE_PRECISION,IOWNER(IMORB),
     +               MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(PSI2,NNN,MPI_DOUBLE_PRECISION,IOWNER(IMORB),
     +               MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(PSIDC,NNN,MPI_DOUBLE_PRECISION,IOWNER(IMORB),
     +               MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(XPSI,NNN,MPI_DOUBLE_PRECISION,IOWNER(IMORB),
     +               MPI_COMM_WORLD,IERROR)
      DO 300 JMORB=IMORB+1,NMORBS
        IF(MYRANK.EQ.IOWNER(JMORB))THEN
          JREC=MOSEQ(JMORB)
          READ(IBAS0,REC=JREC,ERR=9001)(RES0(I),I=1,NNN)
          READ(IBAS2,REC=JREC,ERR=9001)(RES2(I),I=1,NNN)
          READ(IBSDC,REC=JREC,ERR=9001)(RESDC(I),I=1,NNN)
          READ(31,   REC=JREC,ERR=9001)(XRES(I),I=1,NNN)
          PROJ=0.D0
          DO 310 KKK=1,NNN
            PROJ=PROJ+WINTS(KKK)*PSI0(KKK)*RES0(KKK)
310       CONTINUE
          DO 320 KKK=1,NNN
            RES0(KKK) =RES0(KKK) -PROJ*PSI0(KKK)
            RES2(KKK) =RES2(KKK) -PROJ*PSI2(KKK)
            RESDC(KKK)=RESDC(KKK)-PROJ*PSIDC(KKK)
            XRES(KKK) =XRES(KKK) -PROJ*XPSI(KKK)
320       CONTINUE
          WRITE(IBAS0,REC=JREC,ERR=9001)(RES0(I),I=1,NNN)
          WRITE(IBAS2,REC=JREC,ERR=9001)(RES2(I),I=1,NNN)
          WRITE(IBSDC,REC=JREC,ERR=9001)(RESDC(I),I=1,NNN)
          WRITE(31,   REC=JREC,ERR=9001)(XRES(I),I=1,NNN)
        END IF ! MYRANK==JOWNER
300   CONTINUE ! JMORB
      END IF ! IORTH
C
C     NORMALIZE
C
      IF(MYRANK.EQ.IOWNER(IMORB))THEN
      CNORM=SUM(WINTS(:)*PSI0(:)**2)
      FACTOR=1.D0/DSQRT(CNORM)
      EIG=0.D0
      SUMK=0.D0
      PSI0(:) =FACTOR*PSI0(:)
      PSI2(:) =FACTOR*PSI2(:)
      PSIDC(:)=FACTOR*PSIDC(:)
      XPSI(:) =FACTOR*XPSI(:)
      DO 230 KKK=1,NNN
        PSISQ=PSI0(KKK)**2
        TTERM=PSI0(KKK)*PSI2(KKK)
        SUMK=SUMK-0.5D0*WINTS(KKK)*TTERM
        EIG=EIG+WINTS(KKK)*(-0.5D0*TTERM+V(KKK,ISPIN)*PSISQ
     +                      +PSI0(KKK)*XPSI(KKK))
        RHOS(KKK,ISPIN)=RHOS(KKK,ISPIN)+OCCS(IMORB,ISPIN)*PSISQ
        TAUS(KKK,ISPIN)=TAUS(KKK,ISPIN)+OCCS(IMORB,ISPIN)*TTERM
230   CONTINUE
      EIGS(IMORB,ISPIN)=EIG
      EKIN=EKIN+OCCS(IMORB,ISPIN)*SUMK
      WRITE(IO0, REC=IREC,ERR=9001)(PSI0(I),I=1,NNN)
      WRITE(IO2, REC=IREC,ERR=9001)(PSI2(I),I=1,NNN)
      WRITE(IODC,REC=IREC,ERR=9001)(PSIDC(I),I=1,NNN)
      WRITE(31,  REC=IREC,ERR=9001)(XPSI(I),I=1,NNN)
C
      READ(IOX,REC=IMORB,ERR=9001)(STOR(I),I=1,NNN)
      STOR(:)=FMIX*PSI0(:)+(1.D0-FMIX)*STOR(:)
      WRITE(IOX,REC=IMORB,ERR=9001)(STOR(I),I=1,NNN)
C
      END IF ! myrank==iowner
200   CONTINUE
C
C     Get EIGS from branches to root here.
C     Orbital errors and mix ratios too, although they are only
C     needed for output.
C     Each process sends everything in one call and then root picks
C     out that process' orbital energies from it.
C     This is one thing that could be done simpler if we distributed
C     the orbitals blockwise instead of round-robin.
C     TODO  Look into mask operations in F95 for this.
C     Alternatively we could loop over orbitals as we do elsewhere.
      DO 235 IPROC=1,NPROCS-1
        IF (MYRANK.EQ.IPROC) THEN

          call MPI_SEND(EIGS(1,ISPIN),NMORBS,MPI_DOUBLE_PRECISION,0,0,
     +                  MPI_COMM_WORLD,IERROR)
          call MPI_SEND(ERRS,NMORBS,MPI_DOUBLE_PRECISION,0,0,
     +                  MPI_COMM_WORLD,IERROR)
          call MPI_SEND(RATS,NMORBS,MPI_DOUBLE_PRECISION,0,0,
     +                  MPI_COMM_WORLD,IERROR)

        ELSE IF (MYRANK.EQ.0) THEN

          call MPI_RECV(BUF,NMORBS,MPI_DOUBLE_PRECISION,IPROC,0,
     +                  MPI_COMM_WORLD,ISTAT,IERROR)
          do 231 I=1,NMORBS
            IF(IPROC.EQ.IOWNER(I)) EIGS(I,ISPIN)=BUF(I)
231       continue

          call MPI_RECV(BUF,NMORBS,MPI_DOUBLE_PRECISION,IPROC,0,
     +                  MPI_COMM_WORLD,ISTAT,IERROR)
          do 232 I=1,NMORBS
            IF(IPROC.EQ.IOWNER(I)) ERRS(I)=BUF(I)
232       continue

          call MPI_RECV(BUF,NMORBS,MPI_DOUBLE_PRECISION,IPROC,0,
     +                  MPI_COMM_WORLD,ISTAT,IERROR)
          do 233 I=1,NMORBS
            IF(IPROC.EQ.IOWNER(I)) RATS(I)=BUF(I)
233       continue

        END IF ! myrank
235   CONTINUE ! iproc
C
C     Assemble the pieces of the spin- and k.e.-densities on root.
C     We re-use RES* as accumulator arrays here...
      CALL MPI_REDUCE(RHOS(1,ISPIN), RES0, NNN,
     +  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
      CALL MPI_REDUCE(TAUS(1,ISPIN), RES2, NNN,
     +  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
C
C     Root puts them together and re-distributes the result.
      IF (MYRANK.EQ.0) THEN
        RHOS(:,ISPIN)=RES0(:)
        TAUS(:,ISPIN)=RES2(:)
      END IF ! myrank==0
      CALL MPI_BCAST(RHOS(1,ISPIN),NNN,MPI_DOUBLE_PRECISION,0,
     +               MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(TAUS(1,ISPIN),NNN,MPI_DOUBLE_PRECISION,0,
     +               MPI_COMM_WORLD,IERROR)
C
C     Globalize the contents of units 41/42/IOX
C     That is:  Every process will have a complete copy of the PSIs
C     for the use of subroutine XOP.
      DO 250 IMORB=1,NMORBS
        IPROC=IOWNER(IMORB) ! owning process for this orbital
        ! get local XPSI...
        IF(IPROC.EQ.MYRANK)
     &    READ(IOX,REC=IMORB,ERR=9001)(XPSI(I),I=1,NNN)
        call MPI_BCAST(XPSI,NNN,MPI_DOUBLE_PRECISION,IPROC,
     &                 MPI_COMM_WORLD,IERROR)
        ! stash XPSI on all other processes...
        IF(IPROC.NE.MYRANK)
     &    WRITE(IOX,REC=IMORB,ERR=9001)(XPSI(I),I=1,NNN)
250   CONTINUE
C
      IF(MYRANK.EQ.0)THEN
      IF(MUTE.EQ.0.OR.IT.EQ.NITS)THEN
      IF(ISPIN.EQ.1)THEN
        WRITE(IOUT,81)IT
  81    FORMAT(//' ITERATION',I3,':  EIGS     ORB ERRORS     ',
     +           'MIX RATIOS----ACCUMULATED'/' ',64('-'))
      ELSE   ! ispin==1
        WRITE(IOUT,'('' SPIN DOWN'')')
      END IF ! ispin==2
      DO 400 I=1,NMORBS
      RATSUM(I,ISPIN)=RATSUM(I,ISPIN)+RATS(I)
      WRITE(IOUT,'(1X,I4,4F15.6)')I,EIGS(I,ISPIN),ERRS(I),RATS(I),
     +                                             RATSUM(I,ISPIN)
400   CONTINUE
      END IF ! mute==0 or it==nits
      END IF ! myrank==0
C
1000  CONTINUE ! ispin
C
      CALL MPI_REDUCE(EKIN, BUF(1), 1,
     +  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IERROR)
      EKIN=BUF(1)
C
      IF(NSPINS.EQ.1)THEN
        EKIN=2.D0*EKIN
        RHOS(:,2)=RHOS(:,1)
        TAUS(:,2)=TAUS(:,1)
      END IF ! nspins==1
      RHO(:)=RHOS(:,1)+RHOS(:,2)
      CALL MPI_BCAST(EKIN, 1, MPI_DOUBLE_PRECISION, 0,
     +  MPI_COMM_WORLD, IERROR)
C
      RETURN
C
 9001 WRITE(IOUT,*) '*** I/O error in SCFNUMX, node=',MYRANK
9999  ERRFLG=.TRUE.
      RETURN
      END
C
C
C
      SUBROUTINE LCAO1X(IT,EKIN,RHO,RHOS,TAUS)
C
C     LCAO INITIALIZATION OF MOLECULAR ORBITALS.
C
      USE CORES
      USE ERRORHANDLING
      USE GRID
      USE IOUNITS
      USE KEYWORDS
      USE LCAOMatrices
      USE MPI
      USE OccupationNumbers
      USE POTENTIALS, VX=>VXC
      USE SIZES
      IMPLICIT NONE
      ! Arguments
      INTEGER,        INTENT(IN) :: IT
      REAL(KIND(1D0)),INTENT(OUT):: EKIN
      REAL(KIND(1D0)),INTENT(OUT):: RHO(NNN),RHOS(NNN,2),TAUS(NNN,2)
      ! Locals
      REAL(KIND(1D0)) XPSI(NNN),BUF(NBASIS)
      REAL(KIND(1D0)) WORK(3*NBASIS-1)
      REAL(KIND(1D0)) HMAT(NBASIS,NBASIS),VECS(NBASIS,NBASIS)
      REAL(KIND(1D0)) BAS0(NNN),BAS2(NNN),PSI0(NNN),PSI2(NNN)
      REAL(KIND(1D0)) BASDC(NNN),PSIDC(NNN)
      REAL(KIND(1D0)) VSUM,PSISQ,TTERM,HFEIG
      INTEGER II,JJ,I,J,K,KKK,IERR,IMORB,IBAS,IREC,IPROC
C
      EKIN=0.D0
      DO 10 KKK=1,NNN
      RHOS(KKK,1)=0.5D0*RHOCOR(KKK)
      RHOS(KKK,2)=0.5D0*RHOCOR(KKK)
      TAUS(KKK,1)=0.5D0*TAUCOR(KKK)
10    TAUS(KKK,2)=0.5D0*TAUCOR(KKK)
C
      DO 100 II=1,NBASIS
      READ(IBAS0,REC=II)(BAS0(I),I=1,NNN)
      DO 100 JJ=1,II
      READ(IBAS0,REC=JJ)(BAS2(I),I=1,NNN)
      VSUM=0.D0
      DO 110 KKK=1,NNN
110   VSUM=VSUM+WINTS(KKK)*BAS0(KKK)*(V(KKK,1)+VX(KKK,1))*BAS2(KKK)
100   HMAT(JJ,II)=TMAT(JJ,II)+VSUM
C
      VECS(:,:) = HMAT(:,:)
c     TODO optimized dimension of WORK -- see 'man dsygv'
      CALL DSYGV(1,'V','U',NBASIS,VECS,NBASIS,SMAT,NBASIS,EIGS,
     &           WORK,3*NBASIS-1,IERR)
      IF(IERR.NE.0)THEN
      WRITE(IOUT,'(/''*** LCAO1: ERROR IN ROUTINE DSYGV ***'')')
      ERRFLG = .TRUE.
      RETURN
      END IF
C
      IF (MYRANK.EQ.0) THEN
      WRITE(IOUT,81)IT
  81  FORMAT(//'ITERATION',I3,':  LCAO EIGENVALUES'/' ',31('-'))
      END IF ! myrank==0
C
      DO 200 I=1,NMORBS
      EIGS(I,2)=EIGS(I,1)
  200 CONTINUE
C
      IF (MYRANK.EQ.0) THEN
      DO 201 I=1,NMORBS
      WRITE(IOUT,'(1X,I4,5X,2F10.5,F16.6)')
     +          I,OCCS(I,1),OCCS(I,2),EIGS(I,1)
  201 CONTINUE
      WRITE(IOUT,*)'LUMOS:'
      DO 210 I=NMORBS+1,NBASIS
      WRITE(IOUT,'(1X,I4,25X,F16.6)')I,EIGS(I,1)
  210 CONTINUE
      IF(MOPRN.EQ.1.AND.IT.EQ.NITS0)CALL MOLIST(VECS)
      END IF ! myrank==0
C
      DO 300 K=1,NMORBS
      DO 300 I=1,NBASIS
      DO 300 J=1,NBASIS
300   EKIN=EKIN+(OCCS(K,1)+OCCS(K,2))*VECS(I,K)*VECS(J,K)*TMAT(I,J)
C
      IF(IT.NE.NITS0)THEN
C
      DO 1400 IMORB=1,NMORBS
      PSI0(:)=0.D0
      DO 1410 IBAS=1,NBASIS
      READ(IBAS0,REC=IBAS)(BAS0(I),I=1,NNN)
      PSI0(:)=PSI0(:)+VECS(IBAS,IMORB)*BAS0(:)
1410  CONTINUE
      DO 1420 KKK=1,NNN
      PSISQ=PSI0(KKK)**2
      RHOS(KKK,1)=RHOS(KKK,1)+OCCS(IMORB,1)*PSISQ
1420  RHOS(KKK,2)=RHOS(KKK,2)+OCCS(IMORB,2)*PSISQ
1400  CONTINUE
C
      ELSE
C
      DO 2400 IMORB=1,NMORBS
      DO 2401 KKK=1,NNN
      PSI0(KKK)=0.D0
      PSI2(KKK)=0.D0
2401  PSIDC(KKK)=0.D0
      DO 2410 IBAS=1,NBASIS
      READ(IBAS0,REC=IBAS)(BAS0(I),I=1,NNN)
      READ(IBAS2,REC=IBAS)(BAS2(I),I=1,NNN)
      READ(IBSDC,REC=IBAS)(BASDC(I),I=1,NNN)
      DO 2411 KKK=1,NNN
      PSI0(KKK)=PSI0(KKK)+VECS(IBAS,IMORB)*BAS0(KKK)
      PSI2(KKK)=PSI2(KKK)+VECS(IBAS,IMORB)*BAS2(KKK)
2411  PSIDC(KKK)=PSIDC(KKK)+VECS(IBAS,IMORB)*BASDC(KKK)
2410  CONTINUE
      DO 2420 KKK=1,NNN
      PSISQ=PSI0(KKK)**2
      TTERM=PSI0(KKK)*PSI2(KKK)
      RHOS(KKK,1)=RHOS(KKK,1)+OCCS(IMORB,1)*PSISQ
      RHOS(KKK,2)=RHOS(KKK,2)+OCCS(IMORB,2)*PSISQ
      TAUS(KKK,1)=TAUS(KKK,1)+OCCS(IMORB,1)*TTERM
2420  TAUS(KKK,2)=TAUS(KKK,2)+OCCS(IMORB,2)*TTERM
      IF(MYRANK.EQ.IOWNER(IMORB))THEN
      IREC=MOSEQ(IMORB)
      WRITE(MOSA0,REC=IREC,ERR=9001)(PSI0(I),I=1,NNN)
      WRITE(MOSA2,REC=IREC,ERR=9001)(PSI2(I),I=1,NNN)
      WRITE(MODCA,REC=IREC,ERR=9001)(PSIDC(I),I=1,NNN)
      IF(NSPINS.EQ.2)THEN
      WRITE(MOSB0,REC=IREC,ERR=9001)(PSI0(I),I=1,NNN)
      WRITE(MOSB2,REC=IREC,ERR=9001)(PSI2(I),I=1,NNN)
      WRITE(MODCB,REC=IREC,ERR=9001)(PSIDC(I),I=1,NNN)
      END IF ! NSPINS==2
      END IF ! MYRANK==IOWNER
C     Units 41 and 42 are replicated everywhere for the benefit of XOP
      WRITE(41,REC=IMORB,ERR=9001)(PSI0(I),I=1,NNN)
      IF(NSPINS.EQ.2)WRITE(42,REC=IMORB,ERR=9001)(PSI0(I),I=1,NNN)
2400  CONTINUE
C
      DO 800 IMORB=1,NMORBS
      IF (MYRANK.EQ.IOWNER(IMORB)) THEN
      IREC=MOSEQ(IMORB)
      READ(MOSA0,REC=IREC,ERR=9001)(PSI0(I),I=1,NNN)
      READ(MOSA2,REC=IREC,ERR=9001)(PSI2(I),I=1,NNN)
      CALL XOP(1,41,PSI0,XPSI)
      HFEIG=0.D0
      DO 888 KKK=1,NNN
      HFEIG=HFEIG+WINTS(KKK)*PSI0(KKK)*(-0.5D0*PSI2(KKK)
     +                       +V(KKK,1)*PSI0(KKK)+XPSI(KKK))
888   CONTINUE
      EIGS(IMORB,1)=HFEIG
      EIGS(IMORB,2)=HFEIG
      END IF ! myrank==IOwner
800   CONTINUE
C
C     Globalize the eigenvalues just computed
C     TODO cf. loop 235 in subroutine SCFNUM 
      DO 235 IPROC=1,NPROCS-1
        IF (MYRANK.EQ.IPROC) THEN
          call MPI_SEND(EIGS(1,1),NMORBS,MPI_DOUBLE_PRECISION,0,0,
     +                  MPI_COMM_WORLD,IERROR)
        ELSE IF (MYRANK.EQ.0) THEN
          call MPI_RECV(BUF,NMORBS,MPI_DOUBLE_PRECISION,IPROC,0,
     +                  MPI_COMM_WORLD,ISTAT,IERROR)
          do 231 I=1,NMORBS
            IF(IPROC.EQ.IOWNER(I)) EIGS(I,1)=BUF(I)
231       continue
        END IF ! myrank
235   CONTINUE ! iproc
C
      END IF ! IT=NITS0
C
      DO 500 KKK=1,NNN
500   RHO(KKK)=RHOS(KKK,1)+RHOS(KKK,2)
C
      RETURN
C
 9001 WRITE(IOUT,*) '*** I/O error in LCAO1X, node=',MYRANK
 9999 ERRFLG = .TRUE.
      RETURN
      END
C
C
C
      SUBROUTINE LCAO2X(IT,EKIN,RHO,RHOS,TAUS)
C
C     LCAO INITIALIZATION OF MOLECULAR SPIN-ORBITALS.
C
      USE CORES
      USE ERRORHANDLING
      USE GRID
      USE IOUNITS
      USE KEYWORDS
      USE LCAOMatrices
      USE MPI
      USE OccupationNumbers
      USE POTENTIALS, VX=>VXC
      USE SIZES
      IMPLICIT REAL*8(A-H,O-Z)
      ! Arguments
      INTEGER,         INTENT(IN) :: IT
      REAL(KIND(1D0)), INTENT(OUT):: EKIN
      REAL(KIND(1D0)), INTENT(OUT):: RHO(NNN),RHOS(NNN,2),TAUS(NNN,2)
      ! Locals
      REAL(KIND(1D0)) XPSIA(NNN),XPSIB(NNN)
      REAL(KIND(1D0)) BUF(NBASIS),WORK(3*NBASIS-1),B(NBASIS,NBASIS)
      REAL(KIND(1D0)) HMAT(NBASIS,NBASIS,2),VECS(NBASIS,NBASIS,2)
      REAL(KIND(1D0)) BAS0(NNN),BAS2(NNN)
      REAL(KIND(1D0)) PSI0A(NNN),PSI2A(NNN),PSI0B(NNN),PSI2B(NNN)
      REAL(KIND(1D0)) BASDC(NNN),PSIDCA(NNN),PSIDCB(NNN)
      REAL(KIND(1D0)) VSUM1,VSUM2,HFEIG1,HFEIG2
      INTEGER II,JJ,I,J,K,KKK,IERR,IMORB,IBAS,IREC,IPROC,ISPIN
C
      EKIN=0.D0
      DO 10 KKK=1,NNN
      RHOS(KKK,1)=0.5D0*RHOCOR(KKK)
      RHOS(KKK,2)=0.5D0*RHOCOR(KKK)
      TAUS(KKK,1)=0.5D0*TAUCOR(KKK)
10    TAUS(KKK,2)=0.5D0*TAUCOR(KKK)
C
      DO 100 II=1,NBASIS
      READ(IBAS0,REC=II)(BAS0(I),I=1,NNN)
      DO 100 JJ=1,II
      READ(IBAS0,REC=JJ)(BAS2(I),I=1,NNN)
      VSUM1=0.D0
      VSUM2=0.D0
      DO 110 KKK=1,NNN
      VSUM1=VSUM1+WINTS(KKK)*BAS0(KKK)*(V(KKK,1)+VX(KKK,1))*BAS2(KKK)
110   VSUM2=VSUM2+WINTS(KKK)*BAS0(KKK)*(V(KKK,2)+VX(KKK,2))*BAS2(KKK)
      HMAT(JJ,II,1)=TMAT(JJ,II)+VSUM1
100   HMAT(JJ,II,2)=TMAT(JJ,II)+VSUM2
C
      DO 1000 ISPIN=1,2
C
      VECS(:,:,ISPIN) = HMAT(:,:,ISPIN)
      B(:,:) = SMAT(:,:)  ! protect SMAT from being overwritten by DSYGV
c     TODO optimized dimension of WORK -- see 'man dsygv'
      CALL DSYGV(1,'V','U',NBASIS,VECS(:,:,ISPIN),NBASIS,
     &           B,NBASIS,EIGS(:,ISPIN),WORK,3*NBASIS-1,IERR)
      IF(IERR.NE.0)THEN
      WRITE(IOUT,'(/''*** LCAO2X: ERROR IN ROUTINE DSYGV ***'')')
      ERRFLG = .TRUE.
      RETURN
      END IF
C
      IF (MYRANK.EQ.0) THEN
      IF(ISPIN.EQ.1)THEN
      WRITE(IOUT,81)IT
  81  FORMAT(//' ITERATION',I3,':  LCAO EIGENVALUES'/' ',32('-'))
      ELSE
      WRITE(IOUT,83)
  83  FORMAT(/' SPIN DOWN'/' ---------')
      END IF
      DO 200 I=1,NMORBS
      WRITE(IOUT,'(1X,I4,10X,F10.5,5X,F16.6)')
     +                 I,OCCS(I,ISPIN),EIGS(I,ISPIN)
  200 CONTINUE
      WRITE(IOUT,*)'LUMOS:'
      DO 210 I=NMORBS+1,NBASIS
      WRITE(IOUT,'(1X,I4,25X,F16.6)')I,EIGS(I,ISPIN)
  210 CONTINUE
      IF(MOPRN.EQ.1.AND.IT.EQ.NITS0)CALL MOLIST(VECS(1,1,ISPIN))
      END IF ! myrank==0
C
      DO 300 K=1,NMORBS
      DO 300 I=1,NBASIS
      DO 300 J=1,NBASIS
300   EKIN=EKIN+OCCS(K,ISPIN)*VECS(I,K,ISPIN)*VECS(J,K,ISPIN)*TMAT(I,J)
C
1000  CONTINUE
C
      IF(IT.NE.NITS0)THEN
C
      DO 1400 IMORB=1,NMORBS
      PSI0A(:)=0.D0
      PSI0B(:)=0.D0
      DO 1410 IBAS=1,NBASIS
      READ(IBAS0,REC=IBAS)(BAS0(I),I=1,NNN)
      DO 1411 KKK=1,NNN
      PSI0A(KKK)=PSI0A(KKK)+VECS(IBAS,IMORB,1)*BAS0(KKK)
1411  PSI0B(KKK)=PSI0B(KKK)+VECS(IBAS,IMORB,2)*BAS0(KKK)
1410  CONTINUE
      DO 1420 KKK=1,NNN
      RHOS(KKK,1)=RHOS(KKK,1)+OCCS(IMORB,1)*PSI0A(KKK)**2
1420  RHOS(KKK,2)=RHOS(KKK,2)+OCCS(IMORB,2)*PSI0B(KKK)**2
1400  CONTINUE
C
      ELSE
C
      DO 2400 IMORB=1,NMORBS
      IF(MYRANK.EQ.IOWNER(IMORB))THEN
      DO 2401 KKK=1,NNN
      PSI0A(KKK)=0.D0
      PSI2A(KKK)=0.D0
      PSI0B(KKK)=0.D0
      PSI2B(KKK)=0.D0
      PSIDCA(KKK)=0.D0
2401  PSIDCB(KKK)=0.D0
      DO 2410 IBAS=1,NBASIS
      READ(IBAS0,REC=IBAS)(BAS0(I),I=1,NNN)
      READ(IBAS2,REC=IBAS)(BAS2(I),I=1,NNN)
      READ(IBSDC,REC=IBAS)(BASDC(I),I=1,NNN)
      DO 2411 KKK=1,NNN
      PSI0A(KKK)=PSI0A(KKK)+VECS(IBAS,IMORB,1)*BAS0(KKK)
      PSI2A(KKK)=PSI2A(KKK)+VECS(IBAS,IMORB,1)*BAS2(KKK)
      PSI0B(KKK)=PSI0B(KKK)+VECS(IBAS,IMORB,2)*BAS0(KKK)
      PSI2B(KKK)=PSI2B(KKK)+VECS(IBAS,IMORB,2)*BAS2(KKK)
      PSIDCA(KKK)=PSIDCA(KKK)+VECS(IBAS,IMORB,1)*BASDC(KKK)
2411  PSIDCB(KKK)=PSIDCB(KKK)+VECS(IBAS,IMORB,2)*BASDC(KKK)
2410  CONTINUE
      DO 2420 KKK=1,NNN
      RHOS(KKK,1)=RHOS(KKK,1)+OCCS(IMORB,1)*PSI0A(KKK)**2
      RHOS(KKK,2)=RHOS(KKK,2)+OCCS(IMORB,2)*PSI0B(KKK)**2
      TAUS(KKK,1)=TAUS(KKK,1)+OCCS(IMORB,1)*PSI0A(KKK)*PSI2A(KKK)
2420  TAUS(KKK,2)=TAUS(KKK,2)+OCCS(IMORB,2)*PSI0B(KKK)*PSI2B(KKK)
      IREC=MOSEQ(IMORB)
      WRITE(MOSA0,REC=IREC,ERR=9001)(PSI0A(I),I=1,NNN)
      WRITE(MOSA2,REC=IREC,ERR=9001)(PSI2A(I),I=1,NNN)
      WRITE(MOSB0,REC=IREC,ERR=9001)(PSI0B(I),I=1,NNN)
      WRITE(MOSB2,REC=IREC,ERR=9001)(PSI2B(I),I=1,NNN)
      WRITE(MODCA,REC=IREC,ERR=9001)(PSIDCA(I),I=1,NNN)
      WRITE(MODCB,REC=IREC,ERR=9001)(PSIDCB(I),I=1,NNN)
      END IF ! MYRANK==IOWNER
C     Units 41 and 42 are replicated everywhere for the benefit of XOP
      WRITE(41,  REC=IMORB,ERR=9001)(PSI0A(I),I=1,NNN)
      WRITE(42,  REC=IMORB,ERR=9001)(PSI0B(I),I=1,NNN)
2400  CONTINUE
C
      DO 800 IMORB=1,NMORBS
      IF (MYRANK.EQ.IOWNER(IMORB)) THEN
      IREC=MOSEQ(IMORB)
      READ(MOSA0,REC=IREC,ERR=9001)(PSI0A(I),I=1,NNN)
      READ(MOSA2,REC=IREC,ERR=9001)(PSI2A(I),I=1,NNN)
      READ(MOSB0,REC=IREC,ERR=9001)(PSI0B(I),I=1,NNN)
      READ(MOSB2,REC=IREC,ERR=9001)(PSI2B(I),I=1,NNN)
      CALL XOP(1,41,PSI0A,XPSIA)
      CALL XOP(2,42,PSI0B,XPSIB)
      HFEIG1=0.D0
      HFEIG2=0.D0
      DO 888 KKK=1,NNN
      HFEIG1=HFEIG1+WINTS(KKK)*PSI0A(KKK)*(-0.5D0*PSI2A(KKK)
     +                       +V(KKK,1)*PSI0A(KKK)+XPSIA(KKK))
      HFEIG2=HFEIG2+WINTS(KKK)*PSI0B(KKK)*(-0.5D0*PSI2B(KKK)
     +                       +V(KKK,2)*PSI0B(KKK)+XPSIB(KKK))
888   CONTINUE
      EIGS(IMORB,1)=HFEIG1
      EIGS(IMORB,2)=HFEIG2
      END IF ! myrank==IOwner
800   CONTINUE
C
C     Globalize the eigenvalues just computed
C     TODO cf. loop 235 in subroutine SCFNUM -- is there a better way?
      DO 235 IPROC=1,NPROCS-1
        IF (MYRANK.EQ.IPROC) THEN
          call MPI_SEND(EIGS(1,1),NMORBS,MPI_DOUBLE_PRECISION,0,0,
     +                  MPI_COMM_WORLD,IERROR)
          call MPI_SEND(EIGS(1,2),NMORBS,MPI_DOUBLE_PRECISION,0,0,
     +                  MPI_COMM_WORLD,IERROR)
        ELSE IF (MYRANK.EQ.0) THEN
          call MPI_RECV(BUF,NMORBS,MPI_DOUBLE_PRECISION,IPROC,0,
     +                  MPI_COMM_WORLD,ISTAT,IERROR)
          do 231 I=1,NMORBS
            IF(IPROC.EQ.IOWNER(I)) EIGS(I,1)=BUF(I)
231       continue
          call MPI_RECV(BUF,NMORBS,MPI_DOUBLE_PRECISION,IPROC,0,
     +                  MPI_COMM_WORLD,ISTAT,IERROR)
          do 232 I=1,NMORBS
            IF(IPROC.EQ.IOWNER(I)) EIGS(I,2)=BUF(I)
232       continue
        END IF ! myrank
235   CONTINUE ! iproc
C
      END IF ! IT=NITS0
C
      DO 500 KKK=1,NNN
500   RHO(KKK)=RHOS(KKK,1)+RHOS(KKK,2)
C
      RETURN
C
 9001 WRITE(IOUT,*) '*** I/O error in LCAO2X, node=',MYRANK
 9999 ERRFLG = .TRUE.
      RETURN
      END
C
C
C
      SUBROUTINE XOP(ISPIN,IOX,FUN,XFUN)
C
C     HARTREE-FOCK EXCHANGE OPERATOR
C
C     Applies the exchange operator for the ISPIN orbitals stored in
C     file IOX to the function FUN, leaving the result in XFUN.
C
C     Each process has its own complete copy of IOX.  This implementation
C     choice ought to be tested:  Could we devise a scheme whereby PSIs
C     are passed around by the owning nodes when they're needed?  And
C     how much speed do we lose for that gain in space?
C
      USE ERRORHANDLING
      USE SIZES
      USE OccupationNumbers
      IMPLICIT NONE
      INTEGER ISPIN,IOX
      REAL(KIND(1D0)) FUN(NNN),XFUN(NNN)
      REAL(KIND(1D0)) PSI(NNN),PROD(NNN),COUL(NNN)
      INTEGER IMORB,I
C
      XFUN(:)=0.D0
      DO IMORB=1,NMORBS
        READ(IOX,REC=IMORB,ERR=999)(PSI(I),I=1,NNN)
        PROD(:)=PSI(:)*FUN(:)
        CALL POISS(PROD,COUL,0)
        XFUN(:)=XFUN(:)-OCCS(IMORB,ISPIN)*COUL(:)*PSI(:)
      END DO
      RETURN
C
999   CONTINUE
      ERRFLG=.TRUE.
      RETURN
      END
