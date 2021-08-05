**!!HBNDAT.INS
C
C -- Number of hydrogen types
      INTEGER NHTYP
C -- Number of acceptor types
      INTEGER NACTYP
      COMMON /HBNDAT/ NHTYP , NACTYP
C
C -- List of hydrogen types
      INTEGER*2 HYDTYP
C -- List of acceptor types
      INTEGER*2 ACPTYP
      COMMON /HBNDAT/ HYDTYP ( MXHYDB ) , ACPTYP ( MXHYDB )
C
C -- Maximum allowed distance for hydrogen bonding
      DOUBLE PRECISION HBDIST
      COMMON /HBNDAT/  HBDIST
C
      SAVE /HBNDAT/
C
