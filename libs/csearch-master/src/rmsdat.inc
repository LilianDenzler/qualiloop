**!!RMSDAT.INS
C
C -- ????
      DOUBLE PRECISION XRMS
C -- ????
      DOUBLE PRECISION DXRMS
      COMMON /RMSDAT/  XRMS ( NXYZ , MXATM ) , DXRMS ( NXYZ , MXATM )
C
C -- ????
      INTEGER KNDRMS
C -- ????
      INTEGER NATRMS
C -- ????
      INTEGER NRSRMS
C -- ????
      INTEGER NMLRMS
C -- ????
      INTEGER NATFIT
      COMMON /RMSDAT/  KNDRMS ( MXATM ) , NATRMS , NRSRMS , NMLRMS ,
     +                 NATFIT
C
C -- ????
      INTEGER*2 MRMSRS
C -- ????
      INTEGER*2 MRMSX
C -- ????
      INTEGER*2 NRMSAT
C -- ????
      INTEGER*2 NRMSRS
C -- ????
      INTEGER*2 IFSRMS
C -- ????
      INTEGER*2 LSTRMS
C -- ????
      INTEGER*2 LISRMS
C -- ????
      INTEGER*2 LFIT
      COMMON /RMSDAT/  MRMSRS ( MXMOL ) , MRMSX  ( MXMOL ) ,
     +                 NRMSAT ( MXMOL ) , NRMSRS ( MXMOL ) ,
     +                 IFSRMS ( MXRES ) , LSTRMS ( MXRES ) ,
     +                 LISRMS ( MXRES ) , LFIT ( 2 , MXATM )
C
      SAVE /RMSDAT/
C
