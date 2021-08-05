**!!GRPINF.INS
C
C -- ????
      INTEGER*2 LISGRP
C -- First atom in each group
      INTEGER*2 IGRPF
C -- Last atom in each group
      INTEGER*2 IGRPL
C -- ????
      INTEGER*2 MGRPRS
      COMMON /GRPINF/ LISGRP ( MXTGRP ) , IGRPF  ( MXTGRP ) ,
     +                IGRPL  ( MXTGRP ) , MGRPRS ( MXGRP )
C
C -- Switching atom in each group
      INTEGER IGRPSW
      COMMON /GRPINF/ IGRPSW ( MXTGRP )
C
      SAVE /GRPINF/
C
