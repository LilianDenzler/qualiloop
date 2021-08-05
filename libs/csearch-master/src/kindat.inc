**!!KINDAT.INS
C
C -- Type of each atom
      INTEGER KND
C -- Group containing each atom
      INTEGER IGRPAT
C -- Residue containing each atom
      INTEGER IRESAT
      COMMON /KINDAT/ KND( MXTATM ), IGRPAT( MXTATM ), IRESAT( MXTATM )
C
      SAVE /KINDAT/
C
