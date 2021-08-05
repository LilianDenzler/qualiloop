**!!FLXINF.INS
C
C -- ????
      INTEGER NMPPHI
C -- ????
      INTEGER NMPPSI
C -- ????
      INTEGER NMPCH1
C -- ????
      INTEGER NMPCH2
C -- ????
      INTEGER NMPCH3
C -- ????
      INTEGER NMPCH4
      COMMON /FLXINF/ NMPPHI , NMPPSI ,
     +                NMPCH1 , NMPCH2 , NMPCH3 , NMPCH4
C
C -- ????
      INTEGER*2 MAPPHI
C -- ????
      INTEGER*2 MAPPSI
C -- ????
      INTEGER*2 MAPCH1
C -- ????
      INTEGER*2 IC1ATF
C -- ????
      INTEGER*2 IC1ATL
C -- ????
      INTEGER*2 MAPCH2
C -- ????
      INTEGER*2 IC2ATF
C -- ????
      INTEGER*2 IC2ATL
C -- ????
      INTEGER*2 MAPCH3
C -- ????
      INTEGER*2 IC3ATF
C -- ????
      INTEGER*2 IC3ATL
C -- ????
      INTEGER*2 MAPCH4
C -- ????
      INTEGER*2 IC4ATF
C -- ????
      INTEGER*2 IC4ATL
      COMMON /FLXINF/ MAPPHI (MAXMAP), MAPPSI (MAXMAP),
     +                MAPCH1 (MAXMAP), IC1ATF (MAXMAP), IC1ATL (MAXMAP),
     +                MAPCH2 (MAXMAP), IC2ATF (MAXMAP), IC2ATL (MAXMAP),
     +                MAPCH3 (MAXMAP), IC3ATF (MAXMAP), IC3ATL (MAXMAP),
     +                MAPCH4 (MAXMAP), IC4ATF (MAXMAP), IC4ATL (MAXMAP)
C
      SAVE /FLXINF/
C
