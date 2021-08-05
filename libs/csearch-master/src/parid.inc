**!!PARID.INS
C
C -- ????
      INTEGER*2 CGID
C -- ????
      INTEGER*2 BID
C -- ????
      INTEGER*2 TID
C -- ????
      INTEGER*2 PID
C -- ????
      INTEGER*2 OPID
C -- ????
      INTEGER*2 TTID
C -- ????
      INTEGER*2 PNBID
      COMMON /PARID/ CGID ( MXCRGP ) ,          BID ( KABOND, MXBNDP ) ,
     +               TID  ( KAANGL , MXTHTP ) , PID ( KATOR , MXTORP ) ,
     +               OPID ( KAOPLN , MXOPLP ) , TTID ( KATHTH, MXTTP ) ,
     +               PNBID ( KANBND , MXNBDP )
C
C -- ????
      INTEGER NPRM
      COMMON /PARID/ NPRM ( KPARMX )
C
C -- ????
      INTEGER NCHRGP
C -- ????
      INTEGER NBONDP
C -- ????
      INTEGER NTHETP
C -- ????
      INTEGER NPHIP
C -- ????
      INTEGER NOPLNP
C -- ????
      INTEGER NTTP
C -- ????
      INTEGER NNONBP
      EQUIVALENCE ( NPRM(KCHRG) , NCHRGP ) , (NPRM(KBOND)  , NBONDP ) ,
     +            ( NPRM(KANGL) , NTHETP ) , (NPRM(KTOR)   , NPHIP  ) ,
     +            ( NPRM(KOPLN) , NOPLNP ) , (NPRM(KTHTHT) , NTTP   ) ,
     +            ( NPRM(KNBND) , NNONBP )
C
C -- ????
      INTEGER NPHIPW
C -- ????
      INTEGER NPAR
C -- ????
      INTEGER NBC
      COMMON /PARID/ NPHIPW , NPAR , NBC
C
      SAVE /PARID/
C
