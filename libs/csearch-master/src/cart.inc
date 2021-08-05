**!!CART.INS
C
C -- Atomic coordinates
      DOUBLE PRECISION X
C -- Derivatives
      DOUBLE PRECISION D
      COMMON /CART/    X ( NXYZ , MXTATM ) , D ( NXYZ , MXTATM )
C
      SAVE /CART/
C
