**!!PRTOPT.INS
C
C -- ????
      INTEGER INTPR
C -- ????
      INTEGER KENPR
C -- ????
      INTEGER KDERPR
C -- ????
      INTEGER KCORPR
C -- ????
      INTEGER KPRMIN
C -- ????
      INTEGER ITRPR
C -- ????
      INTEGER IMODPR
C -- ????
      INTEGER KTIMPR
C -- ????
      INTEGER KPARTL
C -- ????
      INTEGER KNEIBP
C -- ????
      INTEGER KTORPR
C -- ????
      INTEGER KPRHB
C -- ????
      INTEGER KPRNMR
C -- ????
      INTEGER KRMSPR
C -- ????
      INTEGER KERESP
C -- ????
      INTEGER KRFPR
C -- ????
      INTEGER KCUBPR
C -- ????
      INTEGER KMVRSP
C -- ????
      INTEGER KALLPR
C -- ????
      INTEGER KNBGNP
C -- ????
      INTEGER KRSPRP
C -- ????
      INTEGER KGRPRP
C -- ????
      INTEGER KEIJPR
C -- ????
      INTEGER KBOXPR
      COMMON /PRTOPT/ INTPR  , KENPR  , KDERPR , KCORPR , KPRMIN ,
     +                ITRPR  , IMODPR , KTIMPR , KPARTL , KNEIBP ,
     +                KTORPR , KPRHB  , KPRNMR , KRMSPR , KERESP ,
     +                KRFPR  , KCUBPR , KMVRSP , KALLPR , KNBGNP ,
     +                KRSPRP , KGRPRP , KEIJPR , KBOXPR
C
C -- ????
      DOUBLE PRECISION EPRINT
C -- ????
      DOUBLE PRECISION DERPRT
      COMMON /PRTOPT/  EPRINT , DERPRT
C
      SAVE /PRTOPT/
C
