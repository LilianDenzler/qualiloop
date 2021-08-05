**!!GRPDAT.INS
C
C -- Number of group names
      INTEGER NNMGRP
C -- Names of groups
      INTEGER NAMGRP
      COMMON /GRPDAT/ NNMGRP , NAMGRP ( MXNMGP )
C
      SAVE /GRPDAT/
C
