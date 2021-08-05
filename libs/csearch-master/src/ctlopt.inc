**!!CTLOPT.INS
C
C -- ????
      INTEGER IFREQ
C -- ????
      INTEGER IMODE
C -- ????
      INTEGER IASFR
C -- ????
      INTEGER ITS
C -- ????
      INTEGER IH
C -- ????
      INTEGER KSAVDD
C -- Regenerate symmetry-generated residues every 'NRBOX' iterations
      INTEGER NRBOX
C -- Regenerate neighbouring residue list every 'NRNEIB' iterations
      INTEGER NRNEIB
C -- ????
      INTEGER MXNINT
C -- ????
      INTEGER KRMS
C -- ????
      INTEGER KSQRT
C -- ????
      INTEGER KCHGRS
C -- ????
      INTEGER KSYM
C -- ????
      INTEGER KCRYST
C -- ????
      INTEGER KASYM
      COMMON /CTLOPT/ IFREQ  , IMODE  , IASFR  , ITS  , IH    , KSAVDD ,
     +                NRBOX  , NRNEIB , MXNINT , KRMS , KSQRT , KCHGRS ,
     +                KSYM   , KCRYST , KASYM
C -- ????
      DOUBLE PRECISION FRMS
      COMMON /CTLOPT/  FRMS
C
      SAVE /CTLOPT/
C
