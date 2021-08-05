**!!NUMINF.INS
C
C
C -- This common block stores the number of molecules, residues, atoms
C    etc used for the calculations.
C
C -- Number of molecules.
      INTEGER NMOL
C -- Number of residues.
      INTEGER NRES
C -- Number of groups.
      INTEGER NGROUP
C -- Number of atoms.
      INTEGER NATOM
C -- Number of bonds.
      INTEGER NBOND
C -- Number of angles.
      INTEGER NTHETA
C -- Number of torsions.
      INTEGER NPHI
C -- Number of planes.
      INTEGER NOPLN
C -- Number of angle-anngles.
      INTEGER NTHTH
C -- ????
      INTEGER NINT
C -- ????
      INTEGER NBNDED
C -- ????
      INTEGER NIBNDD
C -- ????
      INTEGER NAT3
C -- ????
      INTEGER NAT32
 
      COMMON /NUMINF/ NMOL, NRES, NGROUP, NATOM, NBOND, NTHETA, NPHI ,
     +                NOPLN, NTHTH, NINT, NBNDED, NIBNDD, NAT3, NAT32
C
      SAVE /NUMINF/
C
