**!!TORPRT.INS
C
C -- ????
      INTEGER NTORPR
      COMMON /TORPRT/ NTORPR
C
C -- ????
      INTEGER*2 KNDTOR
C -- ????
      INTEGER*2 KPHPTR
C -- ????
      INTEGER*2 KTORF
C -- ????
      INTEGER*2 KTORL
      COMMON /TORPRT/ KNDTOR ( MXTRPR ) , KPHPTR ( MXTRPR ) ,
     +                KTORF  ( MXRES  ) , KTORL  ( MXRES )
C
      SAVE /TORPRT/
C
