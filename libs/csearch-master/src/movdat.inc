**!!MOVDAT.INS
C
C -- ????
      DOUBLE PRECISION XMOVDM
      COMMON /MOVDAT/  XMOVDM ( NXYZ , MXATM , 7 )
C
C
C -- ????
      DOUBLE PRECISION DXNDT  ( NXYZ , MXATM , 6 )
      EQUIVALENCE ( XMOVDM (1,1,1) , DXNDT(1,1,1) )
C
C -- ????
      DOUBLE PRECISION XMOVE ( NXYZ , MXATM ) , VELCTY( NXYZ , MXATM ) ,
     +                 ACLRN ( NXYZ , MXATM ) , XPREV ( NXYZ , MXATM ) ,
     +                 XNXT  ( NXYZ , MXATM ) , XTEMP ( NXYZ , MXATM ) ,
     +                 DMOVE ( NXYZ , MXATM )
      EQUIVALENCE ( XMOVE  (1,1) , XMOVDM(1,1,1) ) ,
     +            ( VELCTY (1,1) , XMOVDM(1,1,2) ) ,
     +            ( ACLRN  (1,1) , XMOVDM(1,1,3) ) ,
     +            ( XPREV  (1,1) , XMOVDM(1,1,4) ) ,
     +            ( XNXT   (1,1) , XMOVDM(1,1,5) ) ,
     +            ( XTEMP  (1,1) , XMOVDM(1,1,6) ) ,
     +            ( DMOVE  (1,1) , XMOVDM(1,1,7) )
C
C -- ????
      DOUBLE PRECISION ACLPRV ( NXYZ , MXATM ) , ACPRV1 ( NXYZ , MXATM )
      EQUIVALENCE ( ACLPRV (1,1) , XMOVDM(1,1,4) ) ,
     +            ( ACPRV1 (1,1) , XMOVDM(1,1,5) )
C
C -- ????
      DOUBLE PRECISION XOLD ( NXYZ , MXATM ) , VELOLD ( NXYZ , MXATM )
      EQUIVALENCE ( XOLD   (1,1) , XMOVDM(1,1,4) ) ,
     +            ( VELOLD (1,1) , XMOVDM(1,1,5) )
C
      SAVE /MOVDAT/
C
