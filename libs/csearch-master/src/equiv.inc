**!!EQUIV.INS
C
C -- ????
      INTEGER*2 ICRG
C -- Pointer to charge
      INTEGER*2 ICG
C -- Number and sign of charges
      INTEGER*2 IWG
C -- Pointer to bond parameter
      INTEGER*2 IBE
C -- Pointer to angle parameter
      INTEGER*2 ITE
C -- Pointer to torsion parameter
      INTEGER*2 IPE
C -- Pointer to out-of-plane parameter
      INTEGER*2 IOE
C -- Pointer to non-bond parameter
      INTEGER*2 INB
C -- Out-of-plane flag
      INTEGER*2 IOP
C -- Pointer to aatomic mass
      INTEGER*2 ISQM
      COMMON /EQUIV/ ICRG ( MXATYP ) , ICG ( MXATYP ) , IWG ( MXATYP ) ,
     +               IBE  ( MXATYP ) , ITE ( MXATYP ) , IPE ( MXATYP ) ,
     +               IOE  ( MXATYP ) , INB ( MXATYP ) , IOP ( MXATYP ) ,
     +               ISQM ( MXATYP )
C
C -- Is atom type hydrogen ?
      LOGICAL ISHYD
      COMMON /EQUIV/ ISHYD ( MXATYP )
      SAVE /EQUIV/
C
C
