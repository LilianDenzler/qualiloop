**!!SYMDAT.INS
C
C -- Unit cell vectors - 3 unit cell lengths and 3 unit cell angles
      DOUBLE PRECISION UCVEC
C -- Unit cell vectors - 3 cartesian vectors describing the unit cell
      DOUBLE PRECISION UCVCAR
C -- Matrix to transform crystal to cartesian coordinates
      DOUBLE PRECISION CRYCAR
C -- Matrix to transform cartesian to crystal coordinates
      DOUBLE PRECISION CARCRY
C -- Array of symmetry operator matrices
      DOUBLE PRECISION SYMOPS
C -- Array of translations to go with each operator
      DOUBLE PRECISION TRANOP
C -- Array of matrices to generate symmetry generates atoms in cartesian
C      coordinates ( ie = CARCRY(n)*SYMOPS(n)*CRYCAR(n) )
      DOUBLE PRECISION SYMCR
C -- Coordinates to define basic unit of simulation if less than a unit
C      cell is being used as the "basic" cell.
      DOUBLE PRECISION ABOXA , ABOXB , ABOXC
C -- Cutoff for symmetry generated residues
      DOUBLE PRECISION BOXLIM
C -- Unit cell volume
      DOUBLE PRECISION VOLUM
C -- Minimum and maximum coordinates of periodic boundary box.
      DOUBLE PRECISION BOXMNX, BOXMXX, BOXMNY, BOXMXY, BOXMNZ, BOXMXZ
C
      COMMON /SYMDAT/  UCVEC ( 6 )          , UCVCAR ( 3 , 3 ) ,
     +                 CRYCAR ( 3 , 3 )     , CARCRY ( 3 , 3 ) ,
     +                 SYMOPS ( 3 , 3 , 9 ) , TRANOP ( 3 , 9 ) ,
     +                 SYMCR  ( 3 , 3 , 9 ) , ABOXA , ABOXB , ABOXC ,
     +                 BOXLIM , VOLUM , BOXMNX , BOXMXX , BOXMNY ,
     +                 BOXMXY , BOXMNZ , BOXMXZ
C
C -- Print flag ????
      INTEGER ISYMPR
C -- Number of symmetry operators, excluding identity
      INTEGER NSYMOP
C -- Total unit cells to be included along x-,y- and z-axes
      INTEGER NCELLX , NCELLY , NCELLZ
C -- No. unit cells to be included in 1 direction along x-,y- and z-axes
      INTEGER NCELX2 , NCELY2 , NCELZ2
C -- ????
      INTEGER NSYMCD
C -- ????
      INTEGER NSYMM2
C -- ????
      INTEGER NCLSZY
C -- ????
      INTEGER NCELSZ
      COMMON /SYMDAT/ ISYMPR ( MXCEL3 ) , NSYMOP , NCELLX , NCELLY ,
     +                NCELLZ , NCELX2 , NCELY2 , NCELZ2 , NSYMCD ,
     +                NSYMM2 , NCLSZY , NCELSZ
C
      SAVE /SYMDAT/
C
