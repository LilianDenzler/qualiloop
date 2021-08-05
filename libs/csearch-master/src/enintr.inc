**!!ENINTR.INS
C -- Energy of each bond
      DOUBLE PRECISION EBOND
C -- Energy of each angle
      DOUBLE PRECISION ETHETA
C -- Energy of each torsion
      DOUBLE PRECISION EPHI
C -- Energy of each out-of-plane
      DOUBLE PRECISION EOPLN
      COMMON /ENINTR/  EBOND ( MXBND ) , ETHETA ( MXTHET ) ,
     +                 EPHI  ( MXPHI ) , EOPLN ( MXOPLN )
C
      SAVE /ENINTR/
C
