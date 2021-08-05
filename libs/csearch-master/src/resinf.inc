**!!RESINF.INS
C
C -- ????
      INTEGER*2 LISRES
C -- ????
      INTEGER*2 IFIRST
C -- ????
      INTEGER*2 LAST
C -- ????
      INTEGER*2 MRESML
C -- ????
      INTEGER*2 IRSGPF
C -- ????
      INTEGER*2 IRSGPL
      COMMON /RESINF/ LISRES( MXRES ), IFIRST( MXRES ), LAST  ( MXRES ),
     +                MRESML( MXRES ), IRSGPF( MXRES ), IRSGPL( MXRES )
C
      SAVE /RESINF/
C
