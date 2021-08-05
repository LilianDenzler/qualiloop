**!!INTRA.INS
C
C -- List of values of bond lengths
      DOUBLE PRECISION B
C -- List of values of bond angles
      DOUBLE PRECISION TH
C -- List of values of torsion angles
      DOUBLE PRECISION PH
C -- List of values of out-of-plane "torsion" angles
      DOUBLE PRECISION OPLN
      COMMON /INTRA/   B    ( MXBND  ) , TH ( MXTHET ) , PH ( MXPHI ) ,
     +                 OPLN ( MXOPLN )
C
      SAVE /INTRA/
C
