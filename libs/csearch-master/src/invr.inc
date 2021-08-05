**!!INVR.INS
C
C -- ????
      DOUBLE PRECISION HINC
C -- ????
      DOUBLE PRECISION BEGR2
      COMMON /INVR/    HINC , BEGR2
C
C -- ????
      REAL ATABLE
C -- ????
      REAL BTABLE
C -- ????
      REAL CTABLE
C -- ????
      REAL HINV
      COMMON /INVR/ ATABLE ( MXTABL ) , BTABLE ( MXTABL ) ,
     +              CTABLE ( MXTABL ) , HINV
C
C -- ????
      INTEGER       NI
      COMMON /INVR/ NI
C
      SAVE /INVR/
C
