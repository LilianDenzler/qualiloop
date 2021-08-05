**!!MDENG.INS
C
C -- Number of energies calculated
      INTEGER NENRG
      COMMON /MDENG/ NENRG
C
C -- Temperature
      DOUBLE PRECISION ATEMP
C -- Total energy
      DOUBLE PRECISION ETOT
C -- Potential energy
      DOUBLE PRECISION ETP
C -- Kinetic energy
      DOUBLE PRECISION ETK
      COMMON /MDENG/   ATEMP , ETOT , ETP , ETK
C
      SAVE /MDENG/
C
