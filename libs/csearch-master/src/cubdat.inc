**!!CUBDAT.INS
C
C -- ????
      INTEGER*2 MCUBE
C -- ????
      INTEGER*2 LCUBGP
C -- ????
      INTEGER*2 ICUBGP
      COMMON /CUBDAT/ MCUBE  ( MXCBDM , MXCBDM , MXCBDM ) ,
     +                LCUBGP ( MXCBLS ) , ICUBGP ( MXCBLS )
C
      SAVE /CUBDAT/
C
