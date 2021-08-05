**!!SIDINF.INS
C
C -- ????
      INTEGER NMOVSD
C -- ????
      INTEGER NAMVSD
C -- ????
      INTEGER NMV3SD
C -- ????
      INTEGER N32MSD
C -- ????
      INTEGER NBMVSD
C -- ????
      INTEGER NTMVSD
C -- ????
      INTEGER NPMVSD
C -- ????
      INTEGER NOMVSD
C -- ????
      INTEGER NTTMSD
C -- ????
      INTEGER NRMVSD
      COMMON /SIDINF/ NMOVSD , NAMVSD , NMV3SD , N32MSD , NBMVSD ,
     +                NTMVSD , NPMVSD , NOMVSD , NTTMSD , NRMVSD
C
C -- ????
      INTEGER*2 LMOVSD
C -- ????
      INTEGER*2 MVRSSD
      COMMON /SIDINF/ LMOVSD ( MXRES ) , MVRSSD ( MXRES )
C
      SAVE /SIDINF/
C
