! $Header: /home/ross/cvsroot/parnum/headers.f,v 1.7 2005/03/16 14:45:21 ross Exp $
! Modules in alphabetical order
! ----------------------------------------------------------------------
      module AtomicData

      ! supersedes commons ATOMS1 and ATOMS2
      integer, allocatable :: LMAX0(:)
      integer, allocatable :: NORBS0(:,:)
      real(kind(1d0)), allocatable :: OCC1(:,:,:)
      real(kind(1d0)), allocatable :: OCC2(:,:,:)
      real(kind(1d0)), allocatable :: OCC12(:,:,:)

      end module

! ----------------------------------------------------------------------

      module Constants

      real(kind(1d0)) PI
      real(kind(1d0)) FOURPI
      real(kind(1d0)),parameter :: THRD  = 1.d0/3.d0
      real(kind(1d0)),parameter :: THRD2 = 2.d0/3.d0
      real(kind(1d0)),parameter :: THRD4 = 4.d0/3.d0
      real(kind(1d0)) CVX
      real(kind(1d0)) CEX
      real(kind(1d0)) CTAU
      real(kind(1d0)) FAC(0:40)  ! table of factorials
      real(kind(1d0)) DFAC(0:20) ! double factorials
      real(kind(1d0)),parameter :: SMALL = 1.d-11
      real(kind(1d0)) EPS
      real(kind(1d0)),parameter :: AngstromsPerBohr = 0.52918D0
      real(kind(1d0)) DistanceScale

      contains

      subroutine InitializeConstants

      PI = 4.d0*datan(1.d0)
      FOURPI = 4.d0*PI
      CVX = 2.0d0*(0.75d0/PI)**THRD
      CEX = 1.5d0*(0.75d0/PI)**THRD
      CTAU = 0.6D0*(6.D0*PI*PI)**(2.D0/3.D0)

      ! TODO Raise the N=40 limit on factorials?  We'll run into trouble
      ! TODO with extra-large angular grids and NY/LMAX eventually...
      FAC(0)=1.d0
      do i=1,40
      FAC(i)=i*FAC(i-1)
      end do

      DFAC(0)=1.d0
      do i=1,20
      DFAC(i)=(2*i-1)*DFAC(i-1)
      end do

      EPS=1.d0
      do while (1.d0+EPS.gt.1.d0)
        EPS=EPS/2.d0
      end do
      EPS=3.d0*EPS

      return
      end subroutine

      end module

! ----------------------------------------------------------------------
      module Cores

      integer MCORBS ! max possible core orbitals, for allocation purposes
      real(kind(1d0)), allocatable :: RHOCOR(:)
      real(kind(1d0)), allocatable :: TAUCOR(:)
      real(kind(1d0)), allocatable :: DCOR(:,:)
      integer,         allocatable :: ICORE(:,:,:)
      integer,         allocatable :: ICNTR(:)
      integer,         allocatable :: LCOR(:)

      end module

! ----------------------------------------------------------------------
      module ErrorHandling

      logical ERRFLG 

      ! TODO:  Should provide an error exit routine here
      ! that, for instance, shuts down MPI cleanly.

      end module

! ----------------------------------------------------------------------
      module GeometryOptimization

      ! supersedes XYZFIX and NEWG
      ! and takes over some locals from NEWGEO which need to be saved

      real(kind(1d0)) HMIN
      real(kind(1d0)) HTRANS
      integer, allocatable :: MOVE(:,:) ! (MN,3) Is the nucleus movable?
      real(kind(1d0)),allocatable :: FORCE(:,:)    ! forces on nuclei
      real(kind(1d0)),allocatable :: FORCE0(:,:)   ! core contribution
      real(kind(1d0)),allocatable :: FSTOR4(:,:,:) ! stored forces from
                                     ! previous SCF cycles for averaging
      real(kind(1d0)),allocatable :: DELX(:)
      real(kind(1d0)),allocatable :: F1(:)

      end module

! ----------------------------------------------------------------------
      module Grid

      ! mesh point positions
      real(kind(1d0)), allocatable :: XMESH(:)
      real(kind(1d0)), allocatable :: YMESH(:)
      real(kind(1d0)), allocatable :: ZMESH(:)

      ! integration weights
      real(kind(1d0)), allocatable :: WINTS(:)
      real(kind(1d0)), allocatable :: WNUC(:)

      end module

! ----------------------------------------------------------------------
      module HessianMatrix

      ! supersedes common HESS

      real(kind(1d0)), allocatable :: HESS(:,:)  ! (MN*3,MN*3)
      real(kind(1d0)), allocatable :: HEIGS(:)
      real(kind(1d0)), allocatable :: HVECS(:,:)

      end module

! ----------------------------------------------------------------------
      module IOunits

C     Input/output units

      integer, parameter :: IIN   = 5  ! standard input
      integer, parameter :: IOUT  = 6  ! standard output
      integer, parameter :: IATMS = 8  ! atomic.data
      integer, parameter :: INRGS = 9  ! energy data for G1 tests
   
C     Temporary files

      integer, parameter :: IHRSH = 10 ! HIRSHFELD ATOMIC WEIGHT FUNCTIONS
      integer, parameter :: ICOR0 = 11 ! CORE ORBITALS
      integer, parameter :: ICOR2 = 12 ! LAPLACIANS OF CORE ORBITALS
      integer, parameter :: IBAS0 = 13 ! BASIS FUNCTIONS
      integer, parameter :: IBAS2 = 14 ! LAPLACIANS OF BASIS FUNCTIONS
      integer, parameter :: MOSA0 = 15 ! MOLECULAR ORBITALS (UP)
      integer, parameter :: MOSA2 = 16 ! LAPLACIANS OF MOLECULAR ORBITALS (UP)
      integer, parameter :: MOSB0 = 17 ! MOLECULAR ORBITALS (DOWN)
      integer, parameter :: MOSB2 = 18 ! LAPLACIANS OF MOLECULAR ORBITALS (DOWN)
      integer, parameter :: ICRDC = 19 ! "DECOMPOSED" CORE ORBITALS
      integer, parameter :: IBSDC = 20 ! "DECOMPOSED" BASIS FUNCTIONS
      integer, parameter :: MODCA = 21 ! "DECOMPOSED" MOLECULAR ORBITALS (UP)
      integer, parameter :: MODCB = 22 ! "DECOMPOSED" MOLECULAR ORBITALS (DOWN)

C     TODO add exact-exchange files to this list - 30,31,40,41,42?

      end module

! ----------------------------------------------------------------------
      module KeyWords
      
      integer MUTE  
      integer IBOHR
      integer MOPRN
      integer IOCCS
      integer IMESH
      integer IVSTRT
      integer IPOT
      integer, parameter :: LDA_POT = 0
      integer, parameter :: GGA_POT = 1
      integer, parameter :: XX_POT  = 2
      integer IOPTG
      integer IFRZC
      integer IORTH
      integer I_PRINT_MESH
      integer NITS0
      integer NITS
      integer MDNEG
      integer MGEOMS
      integer NSPINS

      ! TODO Add named constants for IOPTG={-2|-1|0|1}

      end module

! ----------------------------------------------------------------------
      module LCAOMatrices

      ! supersedes EIGS and MATS,
      ! but not LCAO which has been converted to automatic local arrays

      real(kind(1d0)), allocatable :: EIGS(:,:)
      real(kind(1d0)), allocatable :: RATSUM(:,:)
      real(kind(1d0)), allocatable :: SMAT(:,:)
      real(kind(1d0)), allocatable :: TMAT(:,:)
      ! TODO can any of these be deallocated at end of LCAO1/2/1X/2X ?

      end module

! ----------------------------------------------------------------------
      module MPI

      include 'mpif.h'
      integer iError
      integer myRank
      integer nProcs
      integer iStat(MPI_STATUS_SIZE)
      integer lenHost
      character(MPI_MAX_PROCESSOR_NAME) :: myHost

      contains

      function IOwner(iMOrb)
      ! Which process owns orbital iMOrb?
      ! Notice that since processes start at 0 and orbitals at 1, the root
      ! process (0) doesn't own the first orbital (1), it owns the Nth.
      implicit none
      integer IOwner, iMOrb
      IOwner = mod(iMOrb,nProcs)
      return
      end function
 
      function MOSeq(iMOrb)
      ! Returns the sequence number (e.g. direct access record number)
      ! of orbital iMOrb IF AND ONLY IF myrank==IOwner(iMOrb).
      implicit none
      integer MOSeq, iMOrb
      MOSeq = 1 + (iMOrb-1)/nProcs
      return
      end function

      end module

! ----------------------------------------------------------------------
      module NondefaultOccupancies 

      ! supersedes common ALTOCC

      integer, parameter :: MCHNGS = 5
      integer         NCHNGS
      real(kind(1d0)) CHNGOCC(MCHNGS,2)
      integer         ICHANGE(MCHNGS)

      end module

! ----------------------------------------------------------------------
      module NuclearFramework

      ! supersedes common blocks MOLEC1, MOLEC2 and ZSYMBS

      integer  NUCS   ! number of nuclei
      integer  NETCHG ! net charge on system TODO a better place for this?
      integer,        allocatable :: NUCZ(:)  ! nuclear charges
      character(2),   allocatable :: ZSYMB(:) ! atomic symbols
      real(kind(1d0)),allocatable :: XNUC(:)
      real(kind(1d0)),allocatable :: YNUC(:)  ! nuclear coordinates
      real(kind(1d0)),allocatable :: ZNUC(:)

      end module

! ----------------------------------------------------------------------
      module NuclearMeshes

      ! supersedes common blocks POINTS, MESHES, TABS

      integer MR    ! maximum no. of radial points per nucleus
      integer ML    ! maximum no. of angular points per nucleus

      integer, allocatable :: NR(:)    ! no. of radii around each nucleus
      integer, allocatable :: NL(:)    ! no. of angles "
      integer, allocatable :: NY(:)    ! no. of spherical harmonics used
      integer, allocatable :: LMAX(:)  ! max. no. of sph.harmonics usable
      ! LMAX is half the algebraic order of the Lebedev angular integration
      ! grid for each atom, and thus the maximum order of the expansion for
      ! rho used in solving Poisson's equation.  See p.44 of Dickson's
      ! thesis for the argument.
      ! NY is the maximum order of spherical harmonic used in our approximate
      ! solution to Schrodinger's equation.  NY is never larger than LMAX
      ! in a given atomic cell, and by default smaller. It can be thought
      ! of as analogous to the maximum polarization function in a basis set
      ! calculation --- NY=2 being equivalent to D functions, NY=3 to F, etc.

      real(kind(1d0)), allocatable :: RMID(:)    ! radial midpoints (1:nucs)
      real(kind(1d0)), allocatable :: RADS(:,:)  ! radii (1:nr,1:nucs)
      real(kind(1d0)), allocatable :: RQ(:,:)    ! dr/dq
      real(kind(1d0)), allocatable :: RQQ(:,:)   ! d2r/dq2
      real(kind(1d0)), allocatable :: WRADS(:,:) ! radial weights

      ! "Angular points" are represented by X-,Y-, and Z-coordinates
      ! on the unit sphere for computational convenience.
      real(kind(1d0)), allocatable :: XANGS(:,:) ! X-coords of angular pts
      real(kind(1d0)), allocatable :: YANGS(:,:) ! Y-coords " (1:nl,1:nucs)
      real(kind(1d0)), allocatable :: ZANGS(:,:) ! Z-coords "   "
      real(kind(1d0)), allocatable :: WANGS(:,:) ! angular weights

      integer, allocatable :: KKKTAB(:)  ! pointers to "atomic blocks"
                                         ! in global (NNN) arrays

      end module

! ----------------------------------------------------------------------
      module OccupationNumbers

      ! supersedes MOLOCC
      real(kind(1d0)), allocatable :: OCCS(:,:)

      end module

! ----------------------------------------------------------------------
      module Potentials

      real(kind(1d0)), allocatable :: V(:,:)
      real(kind(1d0)), allocatable :: VNUC(:)
      real(kind(1d0)), allocatable :: VEL(:)
      real(kind(1d0)), allocatable :: VXC(:,:)

      end module
  
! ----------------------------------------------------------------------
      module Promolecule

      real(kind(1d0)), allocatable :: RHO0(:)
      real(kind(1d0)), allocatable :: VEL0(:)

      end module

! ----------------------------------------------------------------------
      module SCFConvergence
 
      ! supersedes common TRIX
      real(kind(1d0)) FMIX
      real(kind(1d0)) SING

      end module

! ----------------------------------------------------------------------
      module SecondDerivativeMatrices

      ! supersedes common D2MATS
      ! TODO integrate SecondDerivativeMatrices into module NuclearMeshes

      real(kind(1d0)), allocatable :: D2R(:,:,:)    ! D2R(MR,7,MN)
      real(kind(1d0)), allocatable :: TRID2R(:,:,:) ! TRID2R(MR,3,MN)

      end module

! ----------------------------------------------------------------------
      module Sizes

      integer NNN       ! total number of mesh points
      integer NBASIS    ! number of LCAO basis fcns for SCF initialization
      integer NMORBS    ! number of active molecular orbitals
      integer NCORBS    ! number of frozen core orbitals

      end module

! ----------------------------------------------------------------------
      module SphericalHarmonics

      ! supersedes common YLMS
      ! associated common INTRPS has been replaced with automatic local arrays
      integer MYLM  ! max. order of spherical harmonic expansions
      real(kind(1d0)), allocatable :: YLM(:,:,:)  ! (ML,  0:MYLM,0:MYLM)

      end module

! ----------------------------------------------------------------------
      module SpinPolarization

      ! supersedes SPNPOL
      ! TODO integrate this with NuclearFramework maybe?
      ! TODO or make SPOL local to ATOMIN?
      real(kind(1d0)), allocatable :: SPOL(:,:)   ! (MN,2)

      end module
  
! ----------------------------------------------------------------------
      module Splines

      ! supersedes commons SPLIN1, SPLIN2

      integer, parameter :: MDATA = 1400
      integer,         allocatable :: NSPL(:)
      real(kind(1d0)), allocatable :: RARHO(:,:)
      real(kind(1d0)), allocatable :: RAVEL(:,:)

      end module
