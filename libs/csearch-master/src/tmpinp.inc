**!!TMPINP.INS
C
C
C -- This common block stores the number of residues, atoms etc
C    currently read into the last molecule during the input stage.
C
C -- Number of molecules currently read
      INTEGER  MMOL
C -- Number of molecules in current molecule
      INTEGER  MRES
C -- Number of groups in current molecule
      INTEGER MGROUP
C -- Number of atoms in current molecule
      INTEGER MATOM
C -- Number of bonds in current molecule
      INTEGER MBOND
C -- Number of angles in current molecule
      INTEGER MTHETA
C -- Number of torsions in current molecule
      INTEGER MPHI
C -- Number of out-of-planes in current molecule
      INTEGER MOPLN
C -- Number of angle-angle xross-terms in current molecule
      INTEGER MTHTH
C -- Number of non-bonded interactions in current molecule
      INTEGER MBNDED
C -- Nomove indicator
C      This is set to 1 initially to indicate moving molecules/residues
C      If 'STATIONARY' is specified, the value is set to 0 to indicate
C      all remaining molecules/residues top be read in will be fixed.
      INTEGER NOMOVE
C -- Moving molecule indicator
      INTEGER NOMOVL
C
      COMMON /TMPINP/ MMOL , MRES , MGROUP , MATOM , MBOND , MTHETA ,
     +                MPHI , MOPLN , MTHTH , MBNDED , NOMOVE , NOMOVL
C
      SAVE /TMPINP/
C
