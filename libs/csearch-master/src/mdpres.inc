**!!MDPRES.INS
C
C -- ????
      DOUBLE PRECISION PRESUR
C -- ????
      DOUBLE PRECISION ETKV
C -- ????
      DOUBLE PRECISION VIRIAL
C -- ????
      DOUBLE PRECISION PRESST
C -- ????
      DOUBLE PRECISION PRESRM
C -- ????
      DOUBLE PRECISION ETKM
C -- ????
      DOUBLE PRECISION VIRILM
C -- ????
      DOUBLE PRECISION PRESTM
      COMMON /MDPRES/  PRESUR , ETKV ( 3 , 3 ) , VIRIAL ( 3 , 3 ) ,
     +                 PRESST ( 3 , 3 ) , PRESRM , ETKM ( 3 , 3 ) ,
     +                 VIRILM ( 3 , 3 ) , PRESTM ( 3 , 3 )
C
      SAVE /MDPRES/
C
