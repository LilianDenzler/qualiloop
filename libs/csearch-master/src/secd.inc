**!!SECD.INS
C
C -- ????
      DOUBLE PRECISION DD
      COMMON /SECD/    DD ( MXDD )
C
C -- ????
      DOUBLE PRECISION DSX ( NXYZ , 4 , MXDSX )
      EQUIVALENCE ( DD (1) , DSX (1,1,1) )
C
      SAVE /SECD/
C
