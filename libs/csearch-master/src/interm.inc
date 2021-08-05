**!!INTERM.INS
C
C -- Intermolecular repulsion energy for each molecule
      DOUBLE PRECISION ERPMLR
C -- Intermolecular dispersion energy for each molecule
      DOUBLE PRECISION EDPMLR
C -- Intermolecular non-bonded energy for each molecule
      DOUBLE PRECISION ELJMLR
C -- Intermolecular electrostatic energy for each molecule
      DOUBLE PRECISION ESTMLR
C -- Intermolecular hydrogen-bond energy for each molecule
      DOUBLE PRECISION EHBMLR
      COMMON /INTERM/  ERPMLR ( MXMOL ) , EDPMLR ( MXMOL ) ,
     +                 ELJMLR ( MXMOL ) , ESTMLR ( MXMOL ) ,
     +                 EHBMLR ( MXMOL )
C
      SAVE /INTERM/
C
