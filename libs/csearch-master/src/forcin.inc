**!!FORCIN.INS
C
C -- Forcing potentials included ?
C      = 0  No
C        1  Yes
      INTEGER KFORC
C -- Number of forcing potentials
      INTEGER NFOR
C -- Increase forcing torsion every 'KFINC' iterations
      INTEGER KFINC
      COMMON /FORCIN/ KFORC , NFOR , KFINC
C
C -- Pointer to torsion angles forced.
      INTEGER*2 IPPF
C -- Pointer to molecule containing forcing torsions.
      INTEGER*2 JFMOL
      COMMON /FORCIN/ IPPF ( MXFORC ) , JFMOL ( MXFORC )
C
C -- Value required for torsion ( in radians 0 to pi * 2 ).
      DOUBLE PRECISION FCNTOR
C -- Amount to increase forcing torsion each 'KFINC' steps.
      DOUBLE PRECISION AFINC
C -- Energy due to forcing torsions.
      DOUBLE PRECISION EFRCNG
C -- Gradients from forcing torsions.
      DOUBLE PRECISION DFRCNG
C -- Spring constant for forcing potential.
      DOUBLE PRECISION SPRK
C -- Energy per molecule due to forcing potentials.
      DOUBLE PRECISION EFRMOL
      COMMON /FORCIN/  FCNTOR ( MXFORC ) , AFINC , EFRCNG ,
     +                 DFRCNG ( NXYZ , KATOR , MXFORC ) ,
     +                 SPRK ( MXFORC ) , EFRMOL ( MXMOL )
C
      SAVE /FORCIN/
C
