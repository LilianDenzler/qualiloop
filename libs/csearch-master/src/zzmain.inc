**!!ZZMAIN.INS
C
C -- Channel numbers
C -- ????
      INTEGER LTEMP
C -- ????
      INTEGER LOG1
C -- ????
      INTEGER LOG2
C -- ????
      INTEGER LENG
C -- ????
      INTEGER LINT1
C -- ????
      INTEGER LVIB
C -- ????
      INTEGER LDIFJA
C -- ????
      INTEGER LAB
C -- ????
      INTEGER LPOT
C -- ????
      INTEGER LOLDFR
C -- ????
      INTEGER LNEWFR
      COMMON /ZZMAIN/  LTEMP , LOG1 , LOG2 , LENG , LINT1 , LVIB ,
     +                 LDIFJA , LAB , LPOT , LOLDFR , LNEWFR
C -- CHANNEL NUMBERS
C -- ????
      INTEGER L1
C -- ????
      INTEGER L2
C -- ????
      INTEGER L3
C -- ????
      INTEGER L4
C -- ????
      INTEGER L9
C -- ????
      INTEGER L18
C -- ????
      INTEGER L19
C -- ????
      INTEGER L20
C -- ????
      INTEGER L21
C -- ????
      INTEGER L22
C -- ????
      INTEGER L30
      COMMON /ZZMAIN/  L1 , L2 , L3 , L4 , L9 , L18 , L19 , L20 ,
     +                 L21 , L22 , L30
C -- CALCULATION OPTIONS
C -- ????
      INTEGER ITER
C -- ????
      INTEGER ITIMP
C -- ????
      INTEGER KAP
C -- ????
      INTEGER KDBGPR
C -- ????
      INTEGER KGRANY
C -- ????
      INTEGER KLSQ
C -- ????
      INTEGER KMAP
C -- ????
      INTEGER KMIN
C -- ????
      INTEGER KMOLD
C -- ????
      INTEGER KMSF
C -- ????
      INTEGER KORTEP
C -- ????
      INTEGER KSETUP
C -- ????
      INTEGER KENERG
      COMMON /ZZMAIN/  ITER , ITIMP , KAP , KDBGPR , KGRANY , KLSQ ,
     +                 KMAP , KMIN , KMOLD , KMSF , KORTEP , KSETUP ,
     +                 KENERG
C -- ????
      INTEGER LAP
C -- ????
      INTEGER LFNLCR
C -- ????
      INTEGER LFREQS
C -- ????
      INTEGER LGRANY
C -- ????
      INTEGER LHFILE
C -- ????
      INTEGER LLIB
C -- ????
      INTEGER LLSQPR
C -- ????
      INTEGER LMAPCH
C -- ????
      INTEGER LMAPDM
C -- ????
      INTEGER LMAPFL
C -- ????
      INTEGER LMAPPR
C -- ????
      INTEGER LMINCR
C -- ????
      INTEGER LMSF
C -- ????
      INTEGER LNMODE
C -- ????
      INTEGER LOBS
C -- ????
      INTEGER LORTEP
C -- ????
      INTEGER LPFILE
C -- ????
      INTEGER NLSQ
C -- ????
      INTEGER NSETS
      COMMON /ZZMAIN/  LAP , LFNLCR , LFREQS , LGRANY , LHFILE , LLIB ,
     +                 LLSQPR , LMAPCH , LMAPDM , LMAPFL , LMAPPR ,
     +                 LMINCR , LMSF , LNMODE , LOBS , LORTEP , LPFILE ,
     +                 NLSQ , NSETS
C -- ????
      DOUBLE PRECISION EPSIL
C -- ????
      DOUBLE PRECISION EPSP
      DOUBLE PRECISION ZDERMX
      INTEGER IMINIM
      COMMON /ZZMAIN/  EPSIL , EPSP , ZDERMX , IMINIM
C
      SAVE /ZZMAIN/
C
