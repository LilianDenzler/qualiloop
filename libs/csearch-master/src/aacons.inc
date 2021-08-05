**!!AACONS.INS
C
C
C -- Revised:-
C      20-Sep-1990      JPH    Declare all parameters as integers.
C       9-Oct-1990      JPH    Renamed from 'PARAM.INS' to 'AACONS.INS'.
C      17-Oct-1990      JPH    Define frequently used constants
C                                ( ZERO etc ).
C      22-Oct-1990      JPH    Add potential parameter pointers.
C      24-Oct-1990      JPH    Change KCTOR from 4 to 8 for new
C                                potential file.
C      24-OCT-1990      JPH    Define number of energies and times
C                                stored.
C      25-Oct-1990      JPH    Define some general physical constants.
C       5-Nov-1990      JPH    Define NBITWD.
C       7-Nov-1990      JPH    Define 'NXYZ' as 3 - no. of coordinates !
C      15-NOV-1990      JPH    Define 'MXSYM' - max no. of symmetry ops.
C
C -- Define parameters for maximum atoms, resiudes, groups etc
C
C -- Maximum number of atoms
      INTEGER MXTATM , MXATM
C -- Maximum number of groups
      INTEGER MXTGRP , MXGRP , MXNMGP
C -- Maximum number of residues
      INTEGER MXRES , MXNMRS
C -- Maximum number of molecules
      INTEGER MXMOL
C -- Maximum number of bonds
      INTEGER MXBND
C -- Maximum number of angles
      INTEGER MXTHET
C -- Maximum number of torsions
      INTEGER MXPHI
C -- Maximum number of planes
      INTEGER MXOPLN
C -- Maximum number of angle-angle interactions
      INTEGER MXTHTH
C -- Maximum number of excluded non-bonded interactions
      INTEGER MXECL
C -- Maximum number of neighbouring atoms
      INTEGER MXNEIB
C -- Maximum number of hydrogen bonds
      INTEGER MXHYDB
C -- Maximum number of ????
      INTEGER MXTABL
C -- Maximum number of symmetry operators.
      INTEGER MXSYM
      PARAMETER ( MXTATM =   6500 , MXATM  =   6500 , MXTGRP =   8000 ,
     +            MXGRP  =   2000 , MXNMGP =    100 , MXRES  =    900 ,
     +            MXNMRS =    100 , MXMOL  =    900 , MXBND  =   6560 ,
     +            MXTHET =  12000 , MXPHI  =  14000 , MXOPLN =   1200 ,
     +            MXTHTH =  16000 , MXECL  =  37120 , MXNEIB = 135460 ,
     +            MXHYDB =    100 , MXTABL =    900 , MXSYM  =      9 )
C
C -- Define number of coordinates
C
      INTEGER NXYZ
      PARAMETER ( NXYZ = 3 )
C
C -- Maximum number of ????
      INTEGER MXATM3
C -- Maximum number of ????
      INTEGER MXAT34
C -- Maximum number of ????
      INTEGER MXINT
C -- Maximum number of ????
      INTEGER MXIAT3
      PARAMETER ( MXATM3 = MXATM * NXYZ ,
     +            MXAT34 = MXATM3 * 4 ,
     +            MXINT  = MXBND + MXTHET + MXPHI + MXOPLN ,
     +            MXIAT3 = MXINT + MXATM3 )
C
C -- Maximum number of forcing potentials ( torsions )
      INTEGER MXFORC
C -- Maximum number of ????
      INTEGER MXVSET
C -- Maximum number of ????
      INTEGER MXLINT
C -- Maximum number of atom types
      INTEGER MXATYP
C -- Maximum number of ????
      INTEGER MXRSTP
C -- Maximum number of ????
      INTEGER MXPRTP
C -- Maximum number of charge parameters
      INTEGER MXCRGP
C -- Maximum number of bond parameters
      INTEGER MXBNDP
C -- Maximum number of angle parameters
      INTEGER MXTHTP
C -- Maximum number of torsion parameters
      INTEGER MXTORP
C -- Maximum number of out-of-plane parameters
      INTEGER MXOPLP
C -- Maximum number of angle-angle parameters
      INTEGER MXTTP
C -- Maximum number of ????
      INTEGER MXNBDP
C -- Maximum number of ????
      INTEGER MXPAR
C -- Maximum number of ????
      INTEGER MXTRPR
C -- Maximum number of ????
      INTEGER MXTORV
C -- Maximum number of ????
      INTEGER MXATRS
C -- Maximum number of iterations
      INTEGER MXITR
C -- Maximum number of ????
      INTEGER MXCBDM
C -- Maximum number of ????
      INTEGER MXCBLS
C -- Maximum number of ????
      INTEGER MXDD
      PARAMETER ( MXFORC = 100 , MXVSET =   30 , MXLINT =     20 ,
     +            MXATYP = 100 , MXRSTP =  100 , MXPRTP =   1000 ,
     +            MXCRGP =  40 , MXBNDP =  300 , MXTHTP =    300 ,
     +            MXTORP = 500 , MXOPLP =   50 , MXTTP  =    300 ,
     +            MXNBDP = 300 , MXPAR  = 3000 , MXTRPR =   2200 ,
     +            MXTORV = 500 , MXATRS =  100 , MXITR  =   1200 ,
     +            MXCBDM =  16 , MXCBLS = 6000 , MXDD   = 180000 )
C
C -- Maximum number of ????
      INTEGER     MXDSX
      PARAMETER ( MXDSX = MXDD / 12 )
C
C -- Maximum number of ????
      INTEGER MAXOBS
C -- Maximum number of ????
      INTEGER MAXPP
C -- Maximum number of ????
      INTEGER MAXZZ
C -- Maximum number of ????
      INTEGER MAXLSQ
C -- Maximum number of ????
      INTEGER MAXAT3
C -- Maximum number of ????
      INTEGER MXCEL3
C -- Maximum number of ????
      INTEGER MAXMAP
      PARAMETER ( MAXOBS = 400 , MAXPP  =  50 , MAXZZ  = 2500 ,
     +            MAXLSQ = 100 , MAXAT3 = 300 , MXCEL3 = 8000 ,
     +            MAXMAP = 100 )
C
C -- Maximum number of ????
      INTEGER     MXDUMY
      PARAMETER ( MXDUMY = 9 * MXTATM + MXATM3 )
C
C -- LSDPND is the max no of dependent parameters a single parameter
C    may have in the lsq procedure
C
      INTEGER     LSDPND
      PARAMETER ( LSDPND = 1 )
C
C -- Number of energies and times stored
C
      INTEGER     MXENGS      , MXTIMS
      PARAMETER ( MXENGS = 17 , MXTIMS = 16 )
C
C -- Define max no of bits that can be accessed within a word
C     - use 31 not 32 to avoid confusion with sign bit.
C
      INTEGER NBITWD
      PARAMETER ( NBITWD = 31 )
C
C -- Define:-
C      1. indices for each type of potential parameter term,
C      2. the number of terms associated with each,
C      3. the number of atom types for each type of parameter,
C      4. and the number of atom types for interactions between
C         parameters ( angle and torsion cross terms only ).
C
C -- 1 - Charge parameters
      INTEGER     KCHRG      , KCCHRG     , KACHRG
      PARAMETER ( KCHRG  = 1 , KCCHRG = 1 , KACHRG = 1 )
C -- 2 - Bond stretching parameters
      INTEGER     KBOND      , KCBOND     , KABOND
      PARAMETER ( KBOND  = 2 , KCBOND = 3 , KABOND = 2 )
C -- 3 - Angle/theta bending parameters
      INTEGER     KANGL      , KCANGL     , KAANGL     , KTANGL
      PARAMETER ( KANGL  = 3 , KCANGL = 5 , KAANGL = 3 , KTANGL = 2 )
C -- 4 - Torsion parameters
      INTEGER     KTOR       , KCTOR      , KATOR      , KTTOR
      PARAMETER ( KTOR   = 4 , KCTOR  = 8 , KATOR  = 4 , KTTOR = 2 )
C -- Parameters for wild-card and normal torsions
      INTEGER     KCTORW     , KCTORN
      PARAMETER ( KCTORW = 8 , KCTORN = 5 )
C -- 5 - Out of plane parameters
      INTEGER     KOPLN      , KCOPLN     , KAOPLN
      PARAMETER ( KOPLN  = 5 , KCOPLN = 4 , KAOPLN = 4 )
C -- number of elements in ICOP array
      INTEGER     KNICOP
      PARAMETER ( KNICOP = 2 )
C -- 6 - Theta*theta parameters
      INTEGER     KTHTHT     , KCTHTH     , KATHTH     , KTTHTH
      PARAMETER ( KTHTHT = 6 , KCTHTH = 1 , KATHTH = 4 , KTTHTH = 2 )
C -- 7 - Non-bond parameters
      INTEGER     KNBND      , KCNBND     , KANBND
      PARAMETER ( KNBND  = 7 , KCNBND = 3 , KANBND = 2 )
C -- Max no parameter types
      INTEGER KPARMX
      PARAMETER ( KPARMX = KNBND )
C
C -- Excluded non-bonded interactions:-
C
C -- Hydrogen-bond type
      INTEGER KHB
C -- No parameter
      INTEGER KZEROP
C -- ????
      INTEGER KZEROD
C -- Watts hydrogen bond potentials
      INTEGER KWHB
      PARAMETER ( KHB = -2 , KZEROP = -3 , KZEROD = -4 , KWHB = -5 )
C
C -- Define pointers to molecular dynamics information -
C      used to get information out of 'DXNDT'.
C
      INTEGER KCOOR , KVEL , KACLRN
      PARAMETER ( KCOOR  = 1 , KVEL   = 2 , KACLRN = 3 )
C
C -- Define some frequently used constants.
C
      DOUBLE PRECISION BIG2  , BIG3  , BIG4  , BIG5  , BIG6  , BIG7  ,
     +                 BIG8  , BIG9  , BIG10 , BIG15 ,
     +                 BIG30 , BIG32 , BIG33
      PARAMETER ( BIG2  = 1.0D2   , BIG3  = 1.0D3   , BIG4  = 1.0D4   ,
     +            BIG5  = 1.0D5   , BIG6  = 1.0D6   , BIG7  = 1.0D7   ,
     +            BIG8  = 1.0D8   , BIG9  = 1.0D9   , BIG10 = 1.0D10  ,
     +            BIG15 = 1.0D15  ,
     +            BIG30 = 1.0D30  , BIG32 = 1.0D32  , BIG33 = 1.0D33  )
C
      DOUBLE PRECISION SML1 , SML2 , SML3 , SML4 , SML5 , SML6 , SML7 ,
     +                 SML8 , SML9 , SML10 , SML15 , SML20
      PARAMETER ( SML1  = 1.0D-01 , SML2  = 1.0D-02 , SML3  = 1.0D-03 ,
     +            SML4  = 1.0D-04 , SML5  = 1.0D-05 , SML6  = 1.0D-06 ,
     +            SML7  = 1.0D-07 , SML8  = 1.0D-08 , SML9  = 1.0D-09 ,
     +            SML10 = 1.0D-10 , SML15 = 1.0D-15 , SML20 = 1.0D-20 )
C
      DOUBLE PRECISION ZERO   , ONE    , TWO    , THREE  , FOUR   ,
     +                 FIVE   , SIX    , SEVEN  , EIGHT  , NINE   ,
     +                 TEN    , ELEVEN , TWELV  , THIRTN , FOURTN ,
     +                 FIFTN  , SXTEEN , SEVNTN , EIGHTN , NINETN ,
     +                 TWENTY , TWNTY1 , TWNTY2 , TWNTY3 , TWNTY4 ,
     +                 TWNTY5 , TWNTY6 , TWNTY7 , TWNTY8 , TWNTY9 ,
     +                 THIRTY , THRTY2 , FORTY2 , FIFTY  , SIXTY  ,
     +                 EIGHT8 , NINTY6 , ONE2T  , ONE8T  , TWO51  ,
     +                 TWO98  , THR6T
      PARAMETER ( ZERO   =  0.0D0 , ONE    =  1.0D0 , TWO    =  2.0D0 ,
     +            THREE  =  3.0D0 , FOUR   =  4.0D0 , FIVE   =  5.0D0 ,
     +            SIX    =  6.0D0 , SEVEN  =  7.0D0 , EIGHT  =  8.0D0 ,
     +            NINE   =  9.0D0 , TEN    = 10.0D0 , ELEVEN = 11.0D0 ,
     +            TWELV  = 12.0D0 , THIRTN = 13.0D0 , FOURTN = 14.0D0 ,
     +            FIFTN  = 15.0D0 , SXTEEN = 16.0D0 , SEVNTN = 17.0D0 ,
     +            EIGHTN = 18.0D0 , NINETN = 19.0D0 , TWENTY = 20.0D0 ,
     +            TWNTY1 = 21.0D0 , TWNTY2 = 22.0D0 , TWNTY3 = 23.0D0 ,
     +            TWNTY4 = 24.0D0 , TWNTY5 = 25.0D0 , TWNTY6 = 26.0D0 ,
     +            TWNTY7 = 27.0D0 , TWNTY8 = 28.0D0 , TWNTY9 = 29.0D0 ,
     +            THIRTY = 30.0D0 , THRTY2 = 32.0D0 , FORTY2 = 42.0D0 ,
     +            FIFTY  = 50.0D0 , SIXTY  = 60.0D0 , EIGHT8 = 88.0D0 ,
     +            NINTY6 = 96.0D0 , ONE2T = 120.0D0 , ONE8T = 180.0D0 ,
     +            TWO51 = 251.0D0 , TWO98 = 298.0D0 , THR6T = 360.0D0 )
C
      DOUBLE PRECISION PT2 , HALF , PT9 , PT998 , THIRD ,
     +                 ONEPT2 , ONEPT5 ,
     +                 TWOPT5
      PARAMETER ( PT2 = 0.2D0 , HALF   =  0.5D0 , PT9 = 0.9D0 ,
     +            PT998 = 0.998D0 , THIRD  =  ONE/THREE ,
     +            ONEPT2 =  1.2D0 , ONEPT5 =  1.5D0 ,
     +            TWOPT5 =  2.5D0 )
C
      DOUBLE PRECISION PI , PI2 , DG2RAD , RAD2DG , R2PI
      PARAMETER ( PI = 3.14159250259399D0 , PI2 = 6.28318500518798D0 ,
     +            DG2RAD = 1.74532925199433D-2 , RAD2DG = ONE/DG2RAD ,
     +            R2PI = 0.1591549431D0 )
C
      DOUBLE PRECISION CONV , CUTINF , CONV1 , CONV2
      PARAMETER ( CONV  = 332.054D0 , CUTINF = 500.0D0 ,
     +            CONV1 = 0.208D0   , CONV2  = 18.225D0 )
C
C -- Various physical constants
C      Avogrado's number
C      Boltzmann constant
C      Gas constant
C
      DOUBLE PRECISION AVONO , BOLTZ , GASCON , GASKCN
      PARAMETER ( AVONO  = 6.02252D23 , BOLTZ = 3.2996D-27 ,
     +            GASCON = 1.9872D0   , GASKCN = GASCON*SML3 )
C
C -- Conversion factors
C      Mass conversion - atomic units to kg.
C      And vice versa.
C
      DOUBLE PRECISION AMASCV , RAMSCV
      PARAMETER ( AMASCV = 1.6737D-27 , RAMSCV = ONE/AMASCV )
C
C -- Energy
C      Conversion to joules from kcals.
C      Kinetic energy conversion to kcal mol-1
C          from gm. at. wt. angs**2 femtosecs**-2 mol-1
C          Also conversion to kcal per molecule from gm. at. wt. angs**2
C          femtosecs**-2 per molecule - ekoncv
C          Note. This is inverse of fconv
C -- Forces
C      Conversion from VFF units to newtons (fnwcnv)
C      Convert internal units (gm.at.wt. angs. femtosecs**-2) (fconv)
C
      DOUBLE PRECISION ECONV , EKCONV , FCONV , FNWCNV
      PARAMETER ( ECONV = 4184.0D0    , EKCONV = BIG7/ECONV ,
     +            FCONV = ECONV*SML7  , FNWCNV = ECONV*BIG10 )
C
C -- Temperature
C      Assuming kinetic energy in kcal, tmpacv converts temp.s for atoms
C
      DOUBLE PRECISION TMPACV
      PARAMETER ( TMPACV = TWO/(GASKCN*THREE) )
C
C -- Pressure
C      Convert pressure computed in mdprsr in kcal/mol/angstrom**3 to
C      bars (10**5 Newton mt**-2).
C      Note that the difference in kinetic energy and the virial
C      (which is an energy unit) is in kcal/mol but the volume is
C      the volume of 1 unit cell - we need the volume of a mole of unit
C      cells (to be consistent with the energy) - hence the use of
C      avogadro's number.
C
      DOUBLE PRECISION PRSCNV
      PARAMETER ( PRSCNV = (ECONV/AVONO)*BIG30*SML5 )
C
C -- Speed
C      acceleration conversion from angs femtosec**-2 to mt sec**-2
C      velocity conversion from ang femtosec**-1 to mt/sec
C
      DOUBLE PRECISION ACLCNV , VELCNV
      PARAMETER ( ACLCNV = 1.0D20 , VELCNV = 1.0D5 )
C
C
