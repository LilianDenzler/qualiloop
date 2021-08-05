**!!EPARTL.INS
C
C -- Total energy per atom
      DOUBLE PRECISION ENATEN
C -- Repulsion energy per atom
      DOUBLE PRECISION ENATRP
C -- Dispersion energy per atom
      DOUBLE PRECISION ENATDP
C -- Electrostatic energy per atom
      DOUBLE PRECISION ENATQ
      COMMON /EPARTL/  ENATEN ( MXATM ) , ENATRP ( MXATM ) ,
     +                 ENATDP ( MXATM ) , ENATQ  ( MXATM )
C
      SAVE /EPARTL/
C
