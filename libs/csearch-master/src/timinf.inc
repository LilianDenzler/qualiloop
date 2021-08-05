**!!TIMINF.INS
C
C -- ????
      DOUBLE PRECISION TIMBND
C -- ????
      DOUBLE PRECISION TIMTHT
C -- ????
      DOUBLE PRECISION TIMPHI
C -- ????
      DOUBLE PRECISION TIMOP
C -- ????
      DOUBLE PRECISION TIMBB
C -- ????
      DOUBLE PRECISION TIMBT
C -- ????
      DOUBLE PRECISION TIMTT
C -- ????
      DOUBLE PRECISION TIMUB
C -- ????
      DOUBLE PRECISION TIMNBD
C -- ????
      DOUBLE PRECISION TIMNBS
C -- ????
      DOUBLE PRECISION TIMNBL
C -- ????
      DOUBLE PRECISION TIMNBA
C -- ????
      DOUBLE PRECISION TIMHBD
C -- ????
      DOUBLE PRECISION TIMGBX
C -- ????
      DOUBLE PRECISION TIMCBX
C -- ????
      DOUBLE PRECISION TIMNBG
      COMMON /TIMINF/  TIMBND , TIMTHT , TIMPHI , TIMOP  , TIMBB  ,
     +                 TIMBT  , TIMTT  , TIMUB  , TIMNBD , TIMNBS ,
     +                 TIMNBL , TIMNBA , TIMHBD , TIMGBX , TIMCBX ,
     +                 TIMNBG
C
C -- ????
      DOUBLE PRECISION TIMROU ( MXTIMS )
      EQUIVALENCE ( TIMROU ( 1 ) , TIMBND )
C
C -- ????
      INTEGER NBCALL
C -- ????
      INTEGER NTCALL
C -- ????
      INTEGER NPCALL
C -- ????
      INTEGER NOCALL
C -- ????
      INTEGER NBBCAL
C -- ????
      INTEGER NBTCAL
C -- ????
      INTEGER NTTCAL
C -- ????
      INTEGER NUBCAL
C -- ????
      INTEGER NNBCAL
C -- ????
      INTEGER NNBSCL
C -- ????
      INTEGER NNBLCL
C -- ????
      INTEGER NNBACL
C -- ????
      INTEGER NHBCAL
C -- ????
      INTEGER NCLGBX
C -- ????
      INTEGER NCLCBX
C -- ????
      INTEGER NCLNBG
       COMMON /TIMINF/ NBCALL , NTCALL , NPCALL , NOCALL , NBBCAL ,
     +                 NBTCAL , NTTCAL , NUBCAL , NNBCAL , NNBSCL ,
     +                 NNBLCL , NNBACL , NHBCAL , NCLGBX , NCLCBX ,
     +                 NCLNBG
C
C -- ????
      INTEGER NCALLS ( MXTIMS )
      EQUIVALENCE ( NCALLS ( 1 ) , NBCALL )
C
C -- ????
      INTEGER NROUT
      COMMON /TIMINF/ NROUT
C
      SAVE /TIMINF/
C
