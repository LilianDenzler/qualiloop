**!!NEIBRS.INS
C
C -- ????
      INTEGER NNEIBS
C -- ????
      INTEGER NNEIBR
C -- ????
      INTEGER NNEIB
C -- ????
      INTEGER INEIBF
C -- ????
      INTEGER INEIBL
C -- ????
      INTEGER NEIBOR
C -- ????
      INTEGER IBIT
      COMMON /NEIBRS/ NNEIBS , NNEIBR   , NNEIB ( MXGRP ) ,
     +                INEIBF ( MXGRP )  , INEIBL ( MXGRP ) ,
     +                NEIBOR ( MXNEIB ) , IBIT ( NBITWD )
C
      SAVE /NEIBRS/
C
