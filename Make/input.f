C $Header: /home/ross/cvsroot/parnum/input.f,v 1.13 2005/04/26 14:44:54 ross Exp  $
      SUBROUTINE INPUT(NELS,NV0S)
C
C  This routine is in charge of reading and interpreting the user input
C  into Numol.  It will be called by the root process in a parallel run so
C  that only one process tries to read the input.  It will
C  - read and print the title, 
C  - read the keywords and set the appropriate flags in module KeyWords,
C  - read the molecular geometry and fill module NuclearFramework,
C  - optionally read the mesh sizes, and in any case fill NuclearMeshes,
C  - fill module SpinPolarization,
C  - fill module GeometryOptimization,
C  - sets the nucleus count NUCS and the total grid point count NNN.
C
C  TO DO:
C  - Do we really need to export NV0S?  Or could we reduce complexity  TODO
C    and just export NSPINS, simplifying the LCAO calls later?         TODO
C  - Can we encapsulate NELS somewhere?                                TODO
C 
      USE Constants
      USE ERRORHANDLING
      USE GeometryOptimization
      USE IOUNITS
      USE KEYWORDS
      USE NondefaultOccupancies
      USE NUCLEARFRAMEWORK
      USE NUCLEARMESHES
      USE SIZES
      USE SpinPolarization
      IMPLICIT NONE
      ! Arguments
      INTEGER NELS, NV0S
      ! Locals
      INTEGER I, INUC, NSPOLS, nWords, NL_AVAIL, L_AVAIL
      REAL(KIND(1D0)) SPOL1, SPOL2
      CHARACTER*128 :: InputLine
      CHARACTER*128 :: GeomBuffer(100)
      CHARACTER*1 FirstChar
      CHARACTER*2, PARAMETER :: PeriodicTable(100)= (/
     +  'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     +  'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     +  'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     +  'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     +  'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     +  'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     +  'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     +  'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     +  'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     +  'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'/)
      ! Functions
      INTEGER WordCount

      CALL TITLE
      if (ERRFLG) goto 9999
      CALL KWORDS
      if (ERRFLG) goto 9999
 
C     Count the nuclei so we can allocate the right space
      NUCS=0
      DO I=1,100 ! until we find a section divider or the end-of-file
         READ(IIN,'(A128)',END=9998) GeomBuffer(I)
         IF (GeomBuffer(I)(1:1).EQ.'/') GOTO 50  ! normal exit
         NUCS = NUCS + 1 
      END DO
      WRITE(IOUT,*)'Out of buffer space for geometry input.'
      ERRFLG = .TRUE.
      GOTO 9999
 
   50 CALL ALLOCATE_FRAMEWORK()
 
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'INPUT GEOMETRY AND ATOMIC MESH PARAMETERS:'
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'INUC   SYMB      X        Y        Z          NR'
     +,'    NL    NY   LMAX'
 
      NNN=0
      NELS=-NETCHG
      IF(IBOHR.EQ.1)THEN
         DistanceScale=1.D0
      ELSE
         DistanceScale=1.D0/AngstromsPerBohr
      END IF

      DO 100 I=1,NUCS
 
      InputLine = GeomBuffer(I)
      nWords = WordCount(InputLine)
      SELECT CASE (nWords)
      CASE (4)
        ! Read atomic number and X, Y, and Z coordinates
        READ(InputLine,*)NUCZ(I),XNUC(I),YNUC(I),ZNUC(I)
        NR(I)  = -1  ! Placeholders until default values supplied, below
        NL(I)  = -1
        NY(I)  = -1
        LMAX(I)= -1
      CASE (6)
        ! Read at.no., X, Y, Z coords, radial and angular mesh points
        READ(InputLine,*)NUCZ(I),XNUC(I),YNUC(I),ZNUC(I),NR(I),NL(I)
        NY(I)  = -1
        LMAX(I)= -1
      CASE (8)
        ! Read at.no., X, Y, Z coords, meshes and partial wave orders
        READ(InputLine,*)
     &    NUCZ(I),XNUC(I),YNUC(I),ZNUC(I),NR(I),NL(I),NY(I),LMAX(I)
      CASE DEFAULT
        WRITE(IOUT,*) '*** Wrong number of input items for atom',I
        WRITE(IOUT,*) '*** The offending line is:'
        WRITE(IOUT,*) InputLine
        GOTO 9999
      END SELECT
      ZSYMB(I)=PeriodicTable(NUCZ(I))

      IF (NR(I).LE.0) THEN
      ! Compute default radial mesh size from atomic number
      SELECT CASE (NUCZ(I))
        CASE (1:2);    NR(I)=40
        CASE (3:10);   NR(I)=60
        CASE (11:18);  NR(I)=80
        CASE (19:36);  NR(I)=100
        CASE (37:54);  NR(I)=120
        CASE (55:86);  NR(I)=140
        CASE (87:100); NR(I)=160 
        CASE DEFAULT
          WRITE(IOUT,*) '*** Bad atomic number for atom',I,':',NUCZ(I)
          GOTO 9999
      END SELECT
      END IF

      IF (NL(I).LE.0) THEN ! Supply default angular mesh size
        NL(I)=194
      END IF       
      ! Get closest smaller angular mesh size, and LMAX while we're at it
      CALL LEBEDEV_SIZE(NL(I),NL_AVAIL,L_AVAIL)
      NL(I)=NL_AVAIL

      IF (NY(I).LT.0) THEN
      ! Compute default Schrodinger partial wave order from angular mesh size
      ! TODO revise these default NYs downwards
      SELECT CASE (NL(I))
        CASE (  0:109); NY(I)=3
        CASE (110:193); NY(I)=5
        CASE (194:301); NY(I)=7
        CASE DEFAULT;   NY(I)=9
      END SELECT
      END IF

      IF (LMAX(I).LT.0) THEN
        ! Compute default Poisson-solver partial wave order ditto
        ! TODO revise these default LMAXs downwards for speed improvements
        LMAX(I)=L_AVAIL/2
      ELSE IF (LMAX(I).GT.L_AVAIL/2) THEN
        WRITE(IOUT,*) '*** WARNING! LMAX too large for integration grid'
        WRITE(IOUT,*) 'LMAX read:       ',LMAX(I)
        WRITE(IOUT,*) 'LMAX recommended:',L_AVAIL/2
      END IF

      WRITE(IOUT,'(1X,I4,4X,A2,3X,3F9.4,2X,4I6)')
     +      I,ZSYMB(I),XNUC(I),YNUC(I),ZNUC(I),NR(I),NL(I),NY(I),LMAX(I)
 
      NELS=NELS+NUCZ(I)
      NNN=NNN+NR(I)*NL(I)
      XNUC(I)=DistanceScale*XNUC(I)
      YNUC(I)=DistanceScale*YNUC(I)
      ZNUC(I)=DistanceScale*ZNUC(I)
      MOVE(I,1)=1
      MOVE(I,2)=1
      MOVE(I,3)=1
      IF(DABS(XNUC(I)).LT.1.D-06.AND.
     +   DABS(YNUC(I)).LT.1.D-06.AND.
     +   DABS(ZNUC(I)).LT.1.D-06)THEN
      MOVE(I,1)=0
      MOVE(I,2)=0
      MOVE(I,3)=0
      END IF
 
      SPOL(I,1)=0.D0
      SPOL(I,2)=0.D0
 
100   CONTINUE  ! I=1,NUCS
 
      MR = MAXVAL(NR)
      ML = MAXVAL(NL)

      NV0S=1
      IF(IVSTRT.EQ.1)THEN
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'INITIAL'
      WRITE(IOUT,*)'SPIN POLARIZATIONS:'
      READ(IIN,*)NSPOLS
      DO 110 I=1,NSPOLS
      READ(IIN,*)INUC,SPOL1,SPOL2
      SPOL(INUC,1)=SPOL1
      SPOL(INUC,2)=SPOL2
      IF(SPOL1.NE.SPOL2)NV0S=2
110   WRITE(IOUT,'(1X,I4,2F8.2)')INUC,SPOL1,SPOL2
      END IF
C
      IF(IOCCS.EQ.1)THEN
         READ(IIN,*)NCHNGS
         IF(NCHNGS.GT.MCHNGS)THEN
            WRITE(IOUT,*)'*** Too many orbital occupancy changes'
            WRITE(IOUT,*)'*** Redimension MCHNGS in OrbitalOccupancies'
            GOTO 9999
         END IF
         DO 130 I=1,NCHNGS
c           READ(IIN,*)IMORB,OCC1,OCC2
            READ(IIN,*)ICHANGE(I),CHNGOCC(I,1),CHNGOCC(I,2)
  130    CONTINUE
      END IF
      RETURN
C
 9998 CONTINUE
      WRITE(IOUT,*) '*** Unexpected end of input file...          ***'
      WRITE(IOUT,*) '*** Did you forget the section dividers? (/) ***'
 9999 CONTINUE
      ERRFLG = .TRUE.
      RETURN
      END
C
C
C
      SUBROUTINE TITLE
C
C     READS JOB DESCRIPTION (UP TO 10 LINES)
C     AND PLACES AT HEAD OF STANDARD OUTPUT FILE.
C
      USE IOUNITS
      USE ERRORHANDLING
      IMPLICIT NONE
      INTEGER LINES
      CHARACTER*70 STRING
C
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'=================================================='
      WRITE(IOUT,*)'NUMOL: PRE-RELEASE VERSION                        '
      WRITE(IOUT,*)'COPYRIGHT, 1988-2004, BY A.D. BECKE               '
      WRITE(IOUT,*)'PARALLELIZED BY R.M. DICKSON                      '
      WRITE(IOUT,*)'=================================================='
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'JOB:'
C
      LINES=0
   10 CONTINUE
      READ(IIN,'(A70)',END=99)STRING
      IF(STRING.EQ.'/')RETURN
      LINES=LINES+1
      IF(LINES.EQ.11.AND.STRING.NE.'/')THEN
        WRITE(IOUT,'(/)')
        WRITE(IOUT,*)'*** TITLE TRUNCATED TO TEN LINES ***'
      END IF
      IF(LINES.LE.10) WRITE(IOUT,'(1X,A70)')STRING
      GOTO 10
C  
   99 CONTINUE
      WRITE(IOUT,*) '*** Did you forget the section dividers? (/) ***'
      ERRFLG = .TRUE.
      RETURN
      END
C
C
C
      SUBROUTINE KWORDS
C
C     READS INPUT KEYWORDS AND PARAMETERS.
C
      USE GeometryOptimization
      USE IOUNITS
      USE KEYWORDS
      USE ERRORHANDLING
      USE NUCLEARFRAMEWORK
      USE SCFConvergence
      IMPLICIT NONE
      CHARACTER*40 STRING,KEY(11)
      INTEGER I
C
C     DEFAULT CONDITIONS:
C
      NETCHG=0
      KEY(1)='VERBOSE                                 '
      MUTE=0
      KEY(2)='ANGSTROMS                               '
      IBOHR=0
      KEY(3)='NO MOLIST                               '
      MOPRN=0
      KEY(4)='DEFAULT OCCS                            '
      IOCCS=0
      KEY(5)='DEFAULT MESH                            '
      IMESH=0
      KEY(6)='DEFAULT VSTART                          '
      IVSTRT=0
      KEY(7)='NO GEOMETRY SEARCH                      '
      IOPTG=0
      KEY(8)='FROZEN CORE ORBITALS                    '
      IFRZC=1
      KEY(9)='REORTHOGONALIZE ORBITALS                '
      IORTH=1
      KEY(10)='LDA X-C POTENTIAL                       '
      IPOT = LDA_POT
      KEY(11)='NO PRINT MESH                           '
      I_PRINT_MESH = 0
C
      NITS=30
      NITS0=1
      MGEOMS=30
C
      FMIX=0.3D0
      SING=0.25D0
      HMIN=0.05D0
      HTRANS=-0.05D0
C
      DO 100 I=1,111
      READ(IIN,'(A40)')STRING
C
      IF(STRING.EQ.'CHARGE                                  ')THEN
      READ(IIN,*)NETCHG
C
      ELSE IF(STRING.EQ.'VERBOSE                                 ')THEN
      KEY(1)=STRING
      MUTE=0
      ELSE IF(STRING.EQ.'MUTE                                    ')THEN
      KEY(1)=STRING
      MUTE=1
C
      ELSE IF(STRING.EQ.'ANGSTROMS                               ')THEN
      KEY(2)=STRING
      IBOHR=0
      ELSE IF(STRING.EQ.'BOHR                                    ')THEN
      KEY(2)=STRING
      IBOHR=1
C
      ELSE IF(STRING.EQ.'NO MOLIST                               ')THEN
      KEY(3)=STRING
      MOPRN=0
      ELSE IF(STRING.EQ.'MOLIST                                  ')THEN
      KEY(3)=STRING
      MOPRN=1
C
      ELSE IF(STRING.EQ.'DEFAULT OCCS                            ')THEN
      KEY(4)=STRING
      IOCCS=0
      ELSE IF(STRING.EQ.'NONDEFAULT OCCS                         ')THEN
      KEY(4)=STRING
      IOCCS=1
C
C TODO eliminate NONDEFAULT MESH keyword now that input ignores it
      ELSE IF(STRING.EQ.'DEFAULT MESH                            ')THEN
      KEY(5)=STRING
      IMESH=0
      ELSE IF(STRING.EQ.'NONDEFAULT MESH                         ')THEN
      KEY(5)=STRING
      IMESH=1
C
      ELSE IF(STRING.EQ.'DEFAULT VSTART                          ')THEN
      KEY(6)=STRING
      IVSTRT=0
      ELSE IF(STRING.EQ.'NONDEFAULT VSTART                       ')THEN
      KEY(6)=STRING
      IVSTRT=1
C
      ELSE IF(STRING.EQ.'NO GEOMETRY SEARCH                      ')THEN
      KEY(7)=STRING
      IOPTG=0
      ELSE IF(STRING.EQ.'MINIMUM SEARCH                          ')THEN
      KEY(7)=STRING
      IOPTG=1
      ELSE IF(STRING.EQ.'SADDLE SEARCH                           ')THEN
      KEY(7)=STRING
      IOPTG=-2
      READ(IIN,*)MDNEG
      IF(MDNEG.EQ.0)IOPTG=-1
      IF(MDNEG.EQ.0)KEY(7)='SADDLE SEARCH: INITIAL MODES            '
C
      ELSE IF(STRING.EQ.'FROZEN CORES                            ')THEN
      KEY(8)=STRING
      IFRZC=1
      ELSE IF(STRING.EQ.'UNFREEZE CORES                          ')THEN
      KEY(8)='CORES UNFROZEN                          '
      IFRZC=0
C
      ELSE IF(STRING.EQ.'ORTHOGONALIZE                           ')THEN
      KEY(9)=STRING
      IORTH=1
      ELSE IF(STRING.EQ.'NORTH                                   ')THEN
      KEY(9)='NO ORTHOGONALIZATION                    '
      IORTH=0
C
      ELSE IF(STRING.EQ.'LDA                                     ')THEN
      KEY(10)='LDA EXCHANGE-CORRELATION POTENTIAL      '
      IPOT = LDA_POT
      ELSE IF(STRING.EQ.'GGA                                     ')THEN
      KEY(10)='GGA X-C POTENTIAL (LDA CORES)           '
      IPOT = GGA_POT
      ELSE IF(STRING.EQ.'EXACT EXCHANGE                          ')THEN
      KEY(10)='EXACT EXCHANGE POTENTIAL (ALL-ELECTRON) '
      IPOT = XX_POT
C
      ELSE IF(STRING.EQ.'PRINT MESH                              ')THEN
      KEY(11)='PRINT MESH                              '
      I_PRINT_MESH = 1
C
      ELSE IF(STRING.EQ.'NITS                                    ')THEN
      READ(IIN,*)NITS
      ELSE IF(STRING.EQ.'NITS0                                   ')THEN
      READ(IIN,*)NITS0
      ELSE IF(STRING.EQ.'MGEOMS                                  ')THEN
      READ(IIN,*)MGEOMS
C
      ELSE IF(STRING.EQ.'FMIX                                    ')THEN
      READ(IIN,*)FMIX
      ELSE IF(STRING.EQ.'SING                                    ')THEN
      READ(IIN,*)SING
      ELSE IF(STRING.EQ.'HMIN                                    ')THEN
      READ(IIN,*)HMIN
      ELSE IF(STRING.EQ.'HTRANS                                  ')THEN
      READ(IIN,*)HTRANS
C
      ELSE IF(STRING.EQ.'/                                       ')THEN
      GO TO 101
C
      ELSE
      WRITE(IOUT,'(/)')
      WRITE(IOUT,*)'*** UNRECOGNIZED KEYWORD OR PHRASE ***'
      ERRFLG = .TRUE.
      RETURN
C
      END IF
100   CONTINUE
C
101   WRITE(IOUT,*)' '
C
      IF(IPOT.EQ.XX_POT.AND.IFRZC.NE.0)THEN
      IFRZC=0
      KEY(8)='CORES UNFROZEN (REQ''D BY EXACT EXCHANGE)'
      END IF
C
      IF(IOPTG.EQ.1)NITS=30
      IF(IOPTG.EQ.-2)NITS=30
C
      WRITE(IOUT,81)NETCHG
  81  FORMAT(1X,'CHARGE',I3)
C
      WRITE(IOUT,*)' '
      WRITE(IOUT,*)'RUN CONDITIONS:'
      DO 200 I=1,SIZE(KEY)
200   WRITE(IOUT,'(1X,A40)')KEY(I)
      IF(IOPTG.NE.-1)THEN
      WRITE(IOUT,*)'SCF ITERATIONS (TOT) =',NITS
      WRITE(IOUT,*)'MAX GEOM ITERATIONS (IF CALLED FOR) =',MGEOMS
      WRITE(IOUT,*)' '
      WRITE(IOUT,83)FMIX
  83  FORMAT(1X,'FMIX     = ',F6.2)
      WRITE(IOUT,84)SING
  84  FORMAT(1X,'SING     = ',F6.2)
      WRITE(IOUT,85)NITS0
  85  FORMAT(1X,'NITS0    = ',  I6)
      IF(IOPTG.NE.0)WRITE(IOUT,86)HMIN
  86  FORMAT(1X,'HMIN     = ',F6.2)
      IF(IOPTG.EQ.-2)WRITE(IOUT,87)HTRANS
  87  FORMAT(1X,'HTRANS   = ',F6.2)
      ! TODO collapse these (A8,' = ',F6.2) formats into one
      END IF
C
      RETURN
      END
C
C
C
      SUBROUTINE BCAST_INPUT_SCALARS(NELS,NV0S)
C
C  Broadcast scalar input data from root process to all others.
C
      USE ERRORHANDLING
      USE KEYWORDS
      USE MPI
      USE NondefaultOccupancies
      USE NUCLEARFRAMEWORK
      USE NuclearMeshes
      USE SCFConvergence
      USE SIZES
      IMPLICIT NONE
      INTEGER NELS, NV0S
 
C     TODO:  Figure out which keys we don't need on the subordinate nodes
C     and don't bother sending them.
 
C     KEYWORDS
      CALL MPI_BCAST(MUTE,  1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(IBOHR, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(MOPRN, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(IOCCS, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(IMESH, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(IVSTRT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(IPOT,  1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(IOPTG, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(IFRZC, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(IORTH, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(NITS0, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(NITS,  1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(MDNEG, 1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(MGEOMS,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
      CALL MPI_BCAST(NSPINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
 
C     SCFConvergence
      CALL MPI_BCAST( FMIX,  1, MPI_DOUBLE_PRECISION,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( SING,  1, MPI_DOUBLE_PRECISION,
     +                0, MPI_COMM_WORLD, IERROR )
 
C     NuclearFramework
      CALL MPI_BCAST( NUCS,  1, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( NETCHG,1, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )

C     NuclearMeshes
      CALL MPI_BCAST( MR,    1, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( ML,    1, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
 
C     NondefaultOccupancies
      CALL MPI_BCAST( NCHNGS,1, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )

C     other stuff
      CALL MPI_BCAST( NELS,  1, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( NV0S,  1, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( NNN,   1, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
      RETURN
      END
C
C
C
      SUBROUTINE BCAST_INPUT_ARRAYS()
C
C     Broadcast array data after all nodes have allocated space for them
C
      USE ERRORHANDLING
      USE GeometryOptimization
      USE MPI
      USE NondefaultOccupancies
      USE NUCLEARFRAMEWORK
      USE NUCLEARMESHES
      USE SpinPolarization
      IMPLICIT NONE
 
C     NuclearFramework
      CALL MPI_BCAST( NUCZ, NUCS, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( ZSYMB, NUCS*2, MPI_CHARACTER,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( XNUC, NUCS, MPI_DOUBLE_PRECISION,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( YNUC, NUCS, MPI_DOUBLE_PRECISION,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( ZNUC, NUCS, MPI_DOUBLE_PRECISION,
     +                0, MPI_COMM_WORLD, IERROR )
 
C     NuclearMeshes
      CALL MPI_BCAST( NR, NUCS, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( NL, NUCS, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( NY, NUCS, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( LMAX, NUCS, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )

C     GeometryOptimization
      CALL MPI_BCAST( HMIN,      1, MPI_DOUBLE_PRECISION,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( HTRANS,    1, MPI_DOUBLE_PRECISION,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( MOVE, NUCS*3, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
C
C     NondefaultOccupancies
      CALL MPI_BCAST( ICHANGE,   MCHNGS, MPI_INTEGER,
     +                0, MPI_COMM_WORLD, IERROR )
      CALL MPI_BCAST( CHNGOCC, MCHNGS*2, MPI_DOUBLE_PRECISION,
     +                0, MPI_COMM_WORLD, IERROR )
C
C     SpinPolarization
      CALL MPI_BCAST( SPOL, NUCS*2, MPI_DOUBLE_PRECISION,
     +                0, MPI_COMM_WORLD, IERROR )

      RETURN
      END
C
C
C
      SUBROUTINE ALLOCATE_FRAMEWORK()
      USE GeometryOptimization
      USE NuclearFramework
      USE NuclearMeshes
      USE SpinPolarization
      implicit none
      allocate(NUCZ(NUCS))
      allocate(ZSYMB(NUCS))
      allocate(XNUC(NUCS))
      allocate(YNUC(NUCS))
      allocate(ZNUC(NUCS))
      allocate(NR(NUCS))
      allocate(NY(NUCS))
      allocate(NL(NUCS))
      allocate(LMAX(NUCS))
      allocate(SPOL(NUCS,2))
      allocate(MOVE(NUCS,3))
      allocate(FORCE(NUCS,3))
      allocate(FORCE0(NUCS,3))
      allocate(FSTOR4(NUCS,3,4))
      allocate(DELX(NUCS*3))
      allocate(F1(NUCS*3))
      RETURN
      END
C
C
C
      integer function WordCount(String)
      ! Returns the number of words in the string, a word being defined
      ! crudely as a substring of non-blank characters separated by blanks.
      implicit none
      character*(*) String
      integer i
      logical inWord
      WordCount = 0
      inWord = .false.
      do i=1,len_trim(String)
        if (.not.inWord .and. String(i:i).ne.' ') then
          inWord = .true.
          WordCount = WordCount + 1
        else if (inWord .and. String(i:i).eq.' ') then
          inWord = .false.
        end if
      end do
      return
      end function
