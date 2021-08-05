**!!INTXDT.INS
C
C -- Revised:-
C      10-Oct-1990      JPH    Use parameters for dimensions not numbers
C      25-Oct-1990      JPH    Added IBPW
C
C -- ????
      INTEGER*2 ITBW
C -- ????
      INTEGER*2 IPTW
C -- ????
      INTEGER*2 ITTTW
C -- ????
      INTEGER*2 IPBW
      COMMON /INTXDT/ ITBW  ( KTANGL, MXTHET ) , IPTW ( KTTOR, MXPHI ) ,
     +                ITTTW ( KTTHTH, MXTHTH ) , IPBW ( KTTOR, MXPHI )
C
      SAVE /INTXDT/
C
