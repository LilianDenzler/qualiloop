**!!EURTMP.INC
C 
C INCLUDE file created by RAS on 20/11/91 to tidy up use of TEMPRY
C common block within EUREKA fortran code. After an idea by J.Holland
C
      DOUBLE PRECISION XDUMMY
      COMMON /TEMPRY/ XDUMMY(MXDUMY)
C
C
C
      REAL XCOOR (3 , MXATM ) , VEL ( 3, MXATM )
      EQUIVALENCE ( XCOOR , XDUMMY ) ,
     +            ( XDUMMY ( 3*MXATM + 1 ) , VEL )
C
      INTEGER*2 NNEIBP (MXTGRP) , NNBRSP (MXTGRP)
      LOGICAL INNEIB (MXTGRP)
      EQUIVALENCE ( NNEIBP , XDUMMY ) , 
     +            ( XDUMMY (MXTGRP + 1) , NNBRSP ) ,
     +            ( XDUMMY (2*MXTGRP + 1), INNEIB )
      DOUBLE PRECISION ERESB(MXRES) , EREST(MXRES) , ERESP(MXRES) ,
     +                 ERESOP(MXRES) , ERESNB(MXRES) , ERESLJ(MXRES) ,
     +                 ERESQ(MXRES) , ERESTO(MXRES)
      EQUIVALENCE ( ERESB , XDUMMY ) , 
     +            ( XDUMMY ( MXRES + 1 ) , EREST ) ,
     +            ( XDUMMY ( 2*MXRES + 1 ) , ERESP ) ,
     +            ( XDUMMY ( 3*MXRES + 1 ) , ERESOP ) ,
     +            ( XDUMMY ( 4*MXRES + 1 ) , ERESNB ) ,
     +            ( XDUMMY ( 5*MXRES + 1 ) , ERESLJ ) ,
     +            ( XDUMMY ( 6*MXRES + 1 ) , ERESQ ) ,
     +            ( XDUMMY ( 7*MXRES + 1 ) , ERESTO )
C
      INTEGER*2 IDUMY1(MXRES), IDUMY2(MXRES) ,
     +          IRESBF(MXRES) , JRESBF(MXRES)
      EQUIVALENCE ( IDUMY1 , XDUMMY ) ,
     +            ( XDUMMY ( MXRES + 1 ) , IDUMY2 ) ,
     +            ( XDUMMY ( 2*MXRES + 1 ) , IRESBF ) ,
     +            ( XDUMMY ( 3*MXRES + 1 ) , JRESBF )
C
      INTEGER*2 LMOVSV(MXRES) , MVRSSV(MXRES)
      EQUIVALENCE ( LMOVSV , XDUMMY ) ,
     +            ( XDUMMY ( MXRES + 1 ) , MVRSSV )
C
      DOUBLE PRECISION XTMP(3 , MXTATM) , DTMPP(3 , MXTATM) ,
     +                 DTMPM(3 , MXTATM)
      INTEGER*2 LSTAR(MXATM3)
      EQUIVALENCE ( XTMP , XDUMMY ) ,
     +            ( XDUMMY ( 3*MXTATM + 1 ) , DTMPP ) ,
     +            ( XDUMMY ( 6*MXTATM + 1 ) , DTMPM ) ,
     +            ( XDUMMY ( 9*MXTATM + 1 ) , LSTAR )
C
      INTEGER ATTORS(KATOR , MXFORC )
      EQUIVALENCE ( ATTORS , XDUMMY ) 
C
      DOUBLE PRECISION YA(MXATM3) , P(MXATM3)
      EQUIVALENCE ( YA , XDUMMY ) ,
     +            ( XDUMMY ( MXATM3 + 1 ) , P )
      SAVE /TEMPRY/
C Hope this does some good !!
C RAS 20/11/91
C
