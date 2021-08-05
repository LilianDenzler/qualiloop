**!!ENERG.INS
C
C -- Total energy
      DOUBLE PRECISION ETOTL
C -- Bond strain energy
      DOUBLE PRECISION EB
C -- Angle strain energy
      DOUBLE PRECISION ET
C -- Torsion strain energy
      DOUBLE PRECISION EP
C -- Out-of-plane energy
      DOUBLE PRECISION EOP
C -- Bond-bond cross-term energy
      DOUBLE PRECISION EBB
C -- Bond-angle cross-term energy
      DOUBLE PRECISION EBT
C -- Angle-angle cross-term energy
      DOUBLE PRECISION ETT
C -- Angle-angle-torsion cross-term energy
      DOUBLE PRECISION ETTP
C -- Bond-bond-angle cross-term energy
      DOUBLE PRECISION EBP
C -- Out-of-plane cross-term energy
      DOUBLE PRECISION EOO
C -- Non-bond energy
      DOUBLE PRECISION ELJ
C -- Repulsion energy
      DOUBLE PRECISION ERP
C -- Dispersion energy
      DOUBLE PRECISION EDP
C -- Electrostatic energy
      DOUBLE PRECISION EST
C -- Urey-Bradley energy ( no longer used )
      DOUBLE PRECISION EUB
C -- Hydrogen-bond energy
      DOUBLE PRECISION EHB
C -- Buckingham energy
      DOUBLE PRECISION EBHM
      COMMON /ENERG/   ETOTL , EB  , ET  , EP  , EOP , EBB , EBT , ETT ,
     +                 ETTP  , EBP , EOO , ELJ , ERP , EDP , EST , EUB ,
     +                 EHB   , EBHM
C -- Array of energies
      DOUBLE PRECISION ENRGYS ( MXENGS )
      EQUIVALENCE ( EB , ENRGYS ( 1 ) )
C
      SAVE /ENERG/
C
