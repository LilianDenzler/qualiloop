**!!ATDAT.INS
C
C -- Names of the atom types
*-IF FPS
*-      WORD AT
*-ELSE
      INTEGER AT
*-ENDIF
C -- Number of atom types defined
      INTEGER NAT
      COMMON /ATDAT/ AT ( MXATYP ) , NAT
C -- Atomic massed of each atom type
      DOUBLE PRECISION SQM
      COMMON /ATDAT/   SQM ( MXATYP )
C -- Copy of 'AT' in character variable
      CHARACTER * 2  CATTYP
      COMMON /ATTYP/ CATTYP ( MXATYP )
C
      SAVE /ATTYP/ , /ATDAT/
C
