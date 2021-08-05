**!!PARPTR.INS
C
C -- Function:-
C      These arrays provide pointers relating each internal variable
C      to the parameter to be used when calculating its energy.
C
C      These second element in 'ICOP' is used to indicates whether the
C      variable is involved in a cross-term or not. 1 if it is, 0 if not.
C
C -- Pointers to bond parameters
      INTEGER*2 ICB
C -- Pointers to angle parameters
      INTEGER*2 ICT
C -- Pointers to torsion parameters
      INTEGER*2 ICP
C -- Pointers to out-of-plane parameters
      INTEGER*2 ICOP
C -- Pointers to angle-angle parameters
      INTEGER*2 ICTT
      COMMON /PARPTR/ ICB  ( MXBND ) , ICT ( MXTHET ) , ICP ( MXPHI ) ,
     +                ICOP ( KNICOP , MXOPLN ) , ICTT ( MXTHTH )
C
      SAVE /PARPTR/
C
