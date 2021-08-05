**!!ENMOL.INS
C
C -- All these energies are intra-molecular.
C    See the 'E..MLR' arrays for inter-molecular energies.
C
C -- Total energy per molecule
      DOUBLE PRECISION ETOTML
C -- Bond strain energy per molecule
      DOUBLE PRECISION EBMOL
C -- Angle strain energy per molecule
      DOUBLE PRECISION ETMOL
C -- Torsion strain energy per molecule
      DOUBLE PRECISION EPMOL
C -- Out-of-plane energy per molecule
      DOUBLE PRECISION EOPMOL
C -- Bond-bond cross-term energy per molecule
      DOUBLE PRECISION EBBMOL
C -- Bond-angle cross-term energy per molecule
      DOUBLE PRECISION EBTMOL
C -- Angle-angle cross-term energy per molecule
      DOUBLE PRECISION ETTMOL
C -- Angle-angle-torsion cross-term energy per molecule
      DOUBLE PRECISION ETTPML
C -- Bond-bond-angle cross-term energy per molecule
      DOUBLE PRECISION EBPMOL
C -- Out-of-plane cross-term energy per molecule
      DOUBLE PRECISION EOOMOL
C -- Non-bond energy per molecule
      DOUBLE PRECISION ELJMOL
C -- Repulsion energy per molecule
      DOUBLE PRECISION ERPMOL
C -- Dispersion energy per molecule
      DOUBLE PRECISION EDPMOL
C -- Electrostatic energy per molecule
      DOUBLE PRECISION ESTMOL
C -- Urey-Bradley energy per molecule ( no longer used )
      DOUBLE PRECISION EUBMOL
C -- Hydrogen-bond energy per molecule
      DOUBLE PRECISION EHBMOL
C -- Buckingham energy per molecule
      DOUBLE PRECISION EBHMML
C
      COMMON /ENMOL/  ETOTML( MXMOL ), EBMOL ( MXMOL ), ETMOL ( MXMOL ),
     +                EPMOL ( MXMOL ), EOPMOL( MXMOL ), EBBMOL( MXMOL ),
     +                EBTMOL( MXMOL ), ETTMOL( MXMOL ), ETTPML( MXMOL ),
     +                EBPMOL( MXMOL ), EOOMOL( MXMOL ), ELJMOL( MXMOL ),
     +                ERPMOL( MXMOL ), EDPMOL( MXMOL ), ESTMOL( MXMOL ),
     +                EUBMOL( MXMOL ), EHBMOL( MXMOL ), EBHMML( MXMOL )
C
C -- Array to access all energies
C
      DOUBLE PRECISION ENMOLS ( MXMOL , MXENGS )
      EQUIVALENCE ( EBMOL(1) , ENMOLS(1,1) )
C
      SAVE /ENMOL/
C
