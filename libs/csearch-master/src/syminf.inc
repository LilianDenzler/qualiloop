**!!SYMINF.INS
C
C -- ????
      INTEGER NBRES
C -- ????
      INTEGER NBGRP
C -- ????
      INTEGER NBATOM
C -- ????
      INTEGER NTGRP
C -- ????
      INTEGER NTATOM
C -- ????
      INTEGER NGGRP
C -- ????
      INTEGER NGATOM
C -- ????
      INTEGER IGATST
C -- ????
      INTEGER IGGPST
C -- ????
      INTEGER ISYMGP
      COMMON /SYMINF/ NBRES  , NBGRP  , NBATOM , NTGRP  , NTATOM ,
     +                NGGRP  , NGATOM , IGATST , IGGPST ,
     +                ISYMGP  ( MXTGRP )
C
C -- ????
      INTEGER*2 IBASIC
C -- ????
      INTEGER*2 IGGRPF
C -- ????
      INTEGER*2 IGGRPL
      COMMON /SYMINF/ IBASIC ( MXTGRP ) , IGGRPF ( MXGRP ) ,
     +                IGGRPL ( MXGRP )
C
      SAVE /SYMINF/
C
