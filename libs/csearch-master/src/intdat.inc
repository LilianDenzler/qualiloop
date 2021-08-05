**!!INTDAT.INS
C
C -- List of first atoms in each bondd
      INTEGER IB
C -- List of second atoms in each bondd
      INTEGER JB
      COMMON /INTDAT/ IB ( MXINT ) , JB ( MXINT )
C
C -- List of atoms in each bond
      INTEGER*2 IBW
C -- List of atoms in each angle
      INTEGER*2 ITW
C -- List of atoms in each torsion
      INTEGER*2 IPW
C -- List of atoms in each out-of-plane
      INTEGER*2 IOPW
C -- List of atoms in angle-angle interaction
      INTEGER*2 ITTW
      COMMON /INTDAT/ IBW  ( KABOND , MXBND  ), ITW  ( KAANGL, MXTHET ),
     +                IPW  ( KATOR  , MXPHI  ), IOPW ( KAOPLN, MXOPLN ),
     +                ITTW ( KATHTH , MXTHTH )
C
      SAVE /INTDAT/
C
