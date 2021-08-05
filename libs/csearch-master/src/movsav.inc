**!!MOVSAV.INS
C
C -- ????
      INTEGER NMOVS
C -- ????
      INTEGER NATMVS
C -- ????
      INTEGER NAMV3S
C -- ????
      INTEGER NA32MS
C -- ????
      INTEGER NBMVS
C -- ????
      INTEGER NTHMVS
C -- ????
      INTEGER NPHMVS
C -- ????
      INTEGER NOPMVS
C -- ????
      INTEGER NTTMVS
C -- ????
      INTEGER NRSMVS
      COMMON /MOVSAV/ NMOVS  , NATMVS , NAMV3S , NA32MS , NBMVS ,
     +                NTHMVS , NPHMVS , NOPMVS , NTTMVS , NRSMVS
C
C -- ????
      INTEGER*2 LMOVSV
C -- ????
      INTEGER*2 MVRSSV
      COMMON /MOVSAV/ LMOVSV ( MXRES ) , MVRSSV ( MXRES )
C
      SAVE /MOVSAV/
C
