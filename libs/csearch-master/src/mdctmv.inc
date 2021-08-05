**!!MDCTMV.INS
C
C -- ????
      DOUBLE PRECISION CTRMAS
C -- ????
      DOUBLE PRECISION CTRVEL
C -- ????
      DOUBLE PRECISION AMASML
      COMMON /MDCTMV/  CTRMAS ( NXYZ , MXMOL ) ,
     +                 CTRVEL ( NXYZ , MXMOL ) ,
     +                 AMASML ( MXMOL )
C
      SAVE /MDCTMV/
C
