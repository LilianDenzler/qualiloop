**!!RESDAT.INS
C
C -- Number of names of residues.
      INTEGER NNMRES
C -- List of names of residues.
      INTEGER NAMRES
      COMMON /RESDAT/ NNMRES , NAMRES ( MXNMRS )
C
C -- Character equivalent of 'NAMRES'..
      CHARACTER * 4 CRESNM
      COMMON /RESNAM/ CRESNM ( MXNMRS )
C
      SAVE /RESDAT/ , /RESNAM/
C
