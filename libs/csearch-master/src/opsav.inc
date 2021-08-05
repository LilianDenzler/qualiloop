**!!OPSAV.INS
C
C -- ????
      INTEGER IOPATS
      COMMON /OPSAV/ IOPATS ( 8 )
C
C -- ????
      DOUBLE PRECISION DCODXS
C -- ????
      DOUBLE PRECISION D2ODXS
C -- ????
      DOUBLE PRECISION OPNANS
C -- ????
      DOUBLE PRECISION COSOPS
C -- ????
      COMMON /OPSAV/   DCODXS ( NXYZ , KAOPLN ) ,
     +                 D2ODXS ( NXYZ*KAOPLN , NXYZ*KAOPLN ) ,
     +                 OPNANS , COSOPS
C
      SAVE /OPSAV/
C
