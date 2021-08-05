**!!MOVINF.INS
C
C
C -- ????
      INTEGER NMOVE
C -- Number of moving atoms.
      INTEGER NATMOV
C -- Number of moving groups.
      INTEGER NGRPMV
C -- Number of moving coordinates.
      INTEGER NATMV3
C -- ????
      INTEGER NAT32M
C -- Number of moving bonds.
      INTEGER NBMV
C -- Number of moving angles.
      INTEGER NTHMV
C -- Number of moving torsions..
      INTEGER NPHMV
C -- Number of moving planes.
      INTEGER NOPMV
C -- Number of moving angle-angles.
      INTEGER NTTHMV
C -- Number of moving residues.
      INTEGER NRESMV
      COMMON /MOVINF/ NMOVE , NATMOV , NGRPMV , NATMV3 , NAT32M , NBMV ,
     +                NTHMV , NPHMV  , NOPMV  , NTTHMV , NRESMV
C
C -- Array of pointers to moving residues.
      INTEGER*2 LMOVE
C -- ????
      INTEGER*2 MOVRES
      COMMON /MOVINF/ LMOVE ( MXRES ) , MOVRES ( MXRES )
C
      SAVE /MOVINF/
C
